// Fortan interface to the following functions is
// implemented in ../src/shared/mo_util_stride.f90

#include <stddef.h>

void util_c_loc(const void *ptr_in, const void **ptr_out)
{
  *ptr_out = ptr_in;
}

void util_stride_1d(int *out, int elemsize, void *p1, void *p2)
{
  ptrdiff_t d = (unsigned char *)p2 - (unsigned char *)p1;
  d /= elemsize;
  *out = (int)d;
}

void util_stride_2d(int *out, int elemsize,
                    const void *p1, const void *p2, const void *p3)
{
  ptrdiff_t d = (unsigned char *)p2 - (unsigned char *)p1;
  out[0] = (int)(d / elemsize);
  d = (unsigned char *)p3 - (unsigned char *)p1;
  out[1] = (int)(d / elemsize);
}

