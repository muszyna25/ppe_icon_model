#ifdef __linux
#include <fpu_control.h>

/*
 * disables exceptions due to underflow
 *
 * this seem to be necessary for the INTEL compiler
 * when using -fpe0 and MPI (for strange reasons)
 */

void disable_ufl_exception_(void)
{
  fpu_control_t cw;

  _FPU_GETCW (cw);

  cw |= _FPU_MASK_UM;
  _FPU_SETCW (cw);

  _FPU_GETCW (cw);
}
#else

/* dummy procedure */

void disable_ufl_exception_(void)
{
}

#endif
