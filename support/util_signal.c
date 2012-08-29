/**
 * @file
 * 
 * @brief Signal handling with core dump creation.
 *
 * @author Luis Kornblueh, MPIM
 * @date August 2012
 *
 * @note This function has been developed out of an ECLIB function of ECMWF.
 */

#ifdef _SX

int signal_trap(int *create_dump, int *signals)
{
  /* 
   * Signal handling on SX not POSIX conform, so ignore its handling. 
   * Do not write a message, as this would polute the output when 
   * using MPI. Remove as soon SX is gone.
   */
  
  return 0;
}

#else

#if !(defined (__GNUC__)) && !(defined (__clang__))
#pragma STDC FENV_ACCESS ON
#endif

/** This code is POSIX: set independent of compiler options. */

#define _POSIX_C_SOURCE 200112L

#ifdef __linux__
#include <features.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <signal.h>
#include <unistd.h>

#if (defined __linux__ ) || (defined (__APPLE__) && defined (__MACH__))
#include <execinfo.h>
#include <sys/resource.h>
#ifdef __GLIBC__
#define __USE_GNU
#endif
#include <fenv.h>
#endif

#ifdef _AIX
#include <fptrap.h>
#endif

/* mimic gnu extension of fenv */
 
#if (defined (__APPLE__) && defined (__MACH__))

static int feenableexcept (int excepts)
{
  unsigned short fpu_mask, previous_fpu_mask;
  unsigned int simd_mask;

  excepts &= FE_ALL_EXCEPT;
  __asm__ ("fstcw %0" : "=m" (fpu_mask));
  previous_fpu_mask = (~fpu_mask) & FE_ALL_EXCEPT;
  fpu_mask &= ~excepts;
  __asm__ ("fldcw %0" : : "m" (fpu_mask));
  __asm__ ("stmxcsr %0" : "=m" (simd_mask));
  simd_mask &= ~(excepts << 7);
  __asm__ ("ldmxcsr %0" : : "m" (simd_mask));

  return previous_fpu_mask;
}

static int fedisableexcept (int excepts)
{
  unsigned short int fpu_mask, previous_fpu_mask;
  unsigned int simd_mask;

  excepts &= FE_ALL_EXCEPT;
  __asm__ ("fstcw %0" : "=m" (fpu_mask));
  previous_fpu_mask = (~fpu_mask) & FE_ALL_EXCEPT;
  fpu_mask |= excepts;
  __asm__ ("fldcw %0" : : "m" (fpu_mask));
  __asm__ ("stmxcsr %0" : "=m" (simd_mask));
  simd_mask |= excepts << 7;
  __asm__ ("ldmxcsr %0" : : "m" (simd_mask));

  return previous_fpu_mask;
}

#endif

/* generates poor man's backtrace */ 

#if (defined __linux__ ) || (defined (__APPLE__) && defined (__MACH__))

static void show_backtrace(void)
{
  void *callstack[32];
  char **symbols;
  int frames;
  int i;
  
  frames = backtrace(callstack, 32);
  symbols = backtrace_symbols(callstack, frames);
  
  for (i = 0; i < frames; i++) {
    fprintf(stderr,"%s\n", symbols[i]);
  }
  
  fprintf(stderr, "\nUse addr2line for address to line number conversion.\n\n");
  
  free(symbols);
  
  return;
}

#endif

/**
 * signal_trap_dump is a variable that is a copy of the value of the 
 * create_dump parameter passed into the signal_trap function.
 */

int signal_trap_dump = 0;

/**
 * signal_trap_incr is a special variable that is used to ensure that
 * if more than one thread of the process gets a signal more or less
 * concurrently, only one of them produces a traceback. This makes it
 * easier to see what happened. This variable is of type
 * static volatile sig_atomic_t (defined in <sys/signal.h>.
 *
 * @warning 
 * Some compilers give a warning message that says something like:
 * <tt> Duplicate type qualifier "volatile" ignored. </tt>
 * This warning message should be ignored. Do not be tempted to remove
 * the volatile from the definition.
 * }
 */

static volatile sig_atomic_t signal_trap_incr = 0;

extern void signal_trace(int signo, siginfo_t *sigcode, void *sigcontextptr);

/**
 * @brief Signal handling with core dump.
 *
 * The signal_trap subroutine installs a signal handler for various
 * signals that could cause a program to abort. The handler prints
 * out a traceback and optionally creates a core dump.
 *
 * sigaction is used rather than signal, since with signal
 * there is a timing window that allows other threads within the
 * process to run with the SIG_DFL handler and this could cause the
 * process to abort before the traceback has completed.
 *
 * The static volatile sig_atomic_t variable is used to ensure that
 * only 1 thread is updating it at any one time.
 *
 * Since the CPU does not normaly create a SIGFPE signal on a
 * floating-point error, the program initialises floating-point error
 * trapping (for TRP_INVALID, TRP_DIV_BY_ZERO and TRP_OVERFLOW) if
 * trapping for SIGFPE is requested. See man pages for fp_trap and
 * fp_enable.
 *
 * This funtion works for serial code and also for multi-threaded
 * codes, such as those produced using OpenMP.
 *
 * Using the environment variable USE_SIGNAL_HANDLING by setting to yes 
 * enables this singla handling. If not set or set to a value different 
 * than yes, signal handling will not be switched on. 
 *
 * @param[in] create_dump
 *            A pointer to an "int". If the integer has the value:
 *            0         no core dump is produced;
 *            non-zero  a core dump is produced.
 *
 * @param[in] signals
 *            An array of signals to trap. This list is terminated
 *            by an array element containing the value 0.
 *            If the first element of the array is zero, then this
 *            is taken to mean that the default set of signals,
 *            stored in 'signals_to_trap' below, is to be used.
 *
 * @return error code
 *             -2 will not be used as given by environment variable
 *             -1 an error occurred;
 *              0 floating-point error trapping was not set;
 *            > 0 the mode of floating-point trapping set.
 *                the mode returned is usually 1 (FP_TRAP_SYNC).
 *
 */

int signal_trap(int *create_dump, int *signals)
{
  /**
   * signals_to_trap is an array containing the signals that are to be
   * caught and passed to the signal handler by default. You can add more 
   * signals, just make sure that the list is terminated by 0.
   */
  
  int signals_to_trap[] = {
    SIGFPE,
    SIGILL,
    SIGBUS,
    SIGSEGV,
    SIGXCPU,
    /* <<<<< insert additional default signals between here >>>>>*/
    /* <<<<<                   and here                     >>>>>*/
    0
  };
  
  struct sigaction sa;
  int *sig, *sigs;
  int ret = 0;
  char errmsg[1024];

#if !((defined __linux__ ) || (defined _AIX) || (defined (__APPLE__) && defined (__MACH__)))
  fexcept_t fp_traps;
#endif

  const char *envname = "USE_SIGNAL_HANDLING";
  char *envval = NULL;

  envval = getenv(envname);
  if (envval == NULL)
    {
      return -2;
    }
  if (strncasecmp(envval, "y", (size_t) 1) != 0)
    {
      return -2;
    }
  
  signal_trap_dump = *create_dump;
  sigs = sig = (*signals != 0 ? signals : signals_to_trap);

  memset(&sa, 0, sizeof(struct sigaction));
  sa.sa_flags      = SA_SIGINFO;
  sa.sa_sigaction  = signal_trace;
  sigemptyset(&sa.sa_mask);
  
  /*
   * Add each signal to the signal mask, as we do not want multiple signals
   * interrupting the signal handler.
   */
  
  while (*sig != 0)
    {
      sigaddset(&sa.sa_mask, *sig++);
    }

  sig = sigs;
  while (*sig != 0)
    {
      if (sigaction(*sig, &sa, NULL) == -1)
	{
	  sprintf(errmsg, 
		  "Calling sigaction for signal %d failed\n (file: %s line: %d)\n",
		  (int) *sig, __FILE__, __LINE__);
	  perror(errmsg);
	  exit(1);
	}
      
      if (*sig++ == SIGFPE)
	{

#ifdef _AIX

	  /*
	   * If we are trapping Floating-Point Error, then set the processor in fast trap
	   * mode and enable TRP_INVALID, TRP_DIV_BY_ZERO and TRP_OVERFLOW.
	   */
	  if (((ret = fp_trap(FP_TRAP_FASTMODE)) == FP_TRAP_UNIMPL) || (ret == FP_TRAP_ERROR))
	    {
	      sprintf(errmsg, 
		      "Calling fp_trap failed (returned = %d)\n (file: %s line: %d)\n",
		      ret, __FILE__, __LINE__);
	      perror(errmsg);
	      exit(1);
	    }
	  fp_enable(TRP_INVALID | TRP_DIV_BY_ZERO | TRP_OVERFLOW);

#elif (defined (__GLIBC__)  || (defined (__APPLE__) && defined (__MACH__)))

          /*
	   * This overwrites the FP trap settings of gfortran and leads to a 
           *  working version on Linux systems
	   */
	  if (FE_ALL_EXCEPT != 0)
	    {
	      fedisableexcept(FE_ALL_EXCEPT);
	    }
	  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW); 	  

#else
	  /*
	   * It may work ... or not ...
	   */
	  fegetexceptflag(&fp_traps, FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW); 

#endif

	  ret = 1;
	}
    }
	  
  return ret;
}

/*
 * This is the actual signal-handler that is called by the Kernel when one
 * of the trapped SIGNALs occurs. All but the first thread that cause it to
 * be called are put to sleep for about 5 mins, giving the first thread time
 * to print out the traceback and optionally dump the memory before exiting.
 */

void signal_trace(int signo, siginfo_t *sigcode, void *sigcontextptr)
{
  struct rlimit core_dump_size;

  if (signal_trap_incr++)
    {
      sleep(300);
    }
  else
    {
      (void) setvbuf(stderr, NULL, _IOLBF, BUFSIZ);
      (void) setvbuf(stdout, NULL, _IOLBF, BUFSIZ);
#ifdef _AIX
      if (signal_trap_dump)
	{
	  xl__trcedump(signo, sigcode, sigcontextptr);
	}
      else
	{
	  xl__trce(signo, sigcode, sigcontextptr);
	}
#elif (defined (__GLIBC__) || (defined (__APPLE__) && defined (__MACH__)))
      show_backtrace();
      if (signal_trap_dump)
	{
	  core_dump_size.rlim_cur = RLIM_INFINITY;
	  core_dump_size.rlim_max = RLIM_INFINITY;
          setrlimit(RLIMIT_CORE, &core_dump_size);
	}
#else
      if (signal_trap_dump)
	{
	  core_dump_size.rlim_cur = RLIM_INFINITY;
	  core_dump_size.rlim_max = RLIM_INFINITY;
          setrlimit(RLIMIT_CORE, &core_dump_size);
	}
#endif
      raise(SIGQUIT);
    }
}

#endif
