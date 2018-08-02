#include <string.h>

#include "cfortran.h"

static void
cf_util_memset(void *s, int c, size_t *n)
{
  memset(s, c, *n);
}

static int
cf_util_memcmp(void *s1, void *s2, size_t *n)
{
  return memcmp(s1, s2, *n) != 0;
}

/****************************************************************************/
/* cfortran prototypes:                                                     */

FCALLSCSUB3(cf_util_memset, UTIL_MEMSET, util_memset,
            PVOID, INT, PVOID);

FCALLSCFUN3(LOGICAL, cf_util_memcmp, UTIL_MEMCMP, util_memcmp,
            PVOID, PVOID, PVOID)
