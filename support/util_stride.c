#include <stddef.h>

#include "cfortran.h"

static void
cf_util_stride_1d(int *out, int elemsize, void *p1, void *p2)
{
  ptrdiff_t d = (unsigned char *)p2 - (unsigned char *)p1;
  d /= elemsize;
  *out = (int)d;
}

FCALLSCSUB4(cf_util_stride_1d, UTIL_STRIDE_1D, util_stride_1d,
            PINT, INT, PVOID, PVOID)

static void
cf_util_stride_2d(int *out, int elemsize,
                  const void *p1, const void *p2, const void *p3)
{
  ptrdiff_t d = (unsigned char *)p2 - (unsigned char *)p1;
  out[0] = (int)(d / elemsize);
  d = (unsigned char *)p3 - (unsigned char *)p1;
  out[1] = (int)(d / elemsize);
}

FCALLSCSUB5(cf_util_stride_2d, UTIL_STRIDE_2D, util_stride_2d,
            INTV, INT, PVOID, PVOID, PVOID)
