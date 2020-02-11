/*
 *
 * Generate backtrace, need to get more working solutions. The current once
 * working are:
 *
 *     Linux/glibc  gcc/gfortran
 *     Linux/glibc  gcc/NAG f95
 *
 * Luis Kornblueh, MPIM, January 2008
 *
 */

// Fortan interface to the following functions is
// implemented in ../src/shared/mo_util_backtrace.f90

#include "config.h"

#if defined (HAVE_EXECINFO_H)

#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

#if defined (__GNUC__)

void util_backtrace(void)
{
  void *array[32];
  size_t size;
  char **strings;
  size_t i;

  fflush(stderr);
  fflush(stdout);
  raise(SIGSEGV);

  size = backtrace (array, 32);
  strings = backtrace_symbols (array, size);

  for (i = 0; i < size; i++) {
    fprintf (stderr,"%s\n", strings[i]);
  }

  fprintf (stderr, "Use addr2line for addresses to line number conversion.\n");

  free (strings);

  return;
}

#else

void util_backtrace(void)
{
  fprintf (stderr, "Stack traceback not available.\n");

  return;
}

#endif

#elif  defined (HAVE_UCONTEXT_H)

#include <ucontext.h>

#if (defined(sun) || defined(__sun)) && (defined(__SVR4) || defined(__svr4__))

void util_backtrace(void)
{
  (void) printstack(2); /* fd 2 is always stderr */

  return;
}

#else

void util_backtrace(void)
{
  fprintf (stderr, "Stack traceback not available.\n");

  return;
}

#endif
 
#else

void util_backtrace(void)
{
  fprintf (stderr, "Stack traceback not available.\n");

  return;
}

#endif
