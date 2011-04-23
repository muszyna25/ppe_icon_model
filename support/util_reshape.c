#include <string.h>

#include "cfortran.h"

/****************************************************************************/
/* cfortran prototypes:                                                     */

/*
  PDOUBLE on old CRAY PVPs has been 128bit, but that is not true anymore
*/
/* #if  defined  (CRAY) */
/* #  define  PDOUBLE  PFLOAT */
/* #endif */

void cf_util_reshape(double *dest, double *source, int amount);
FCALLSCSUB3(cf_util_reshape, UTIL_RESHAPE, util_reshape,
	    PDOUBLE, PDOUBLE, INT)

void cf_util_reshape2(double *dest, double *source, int ddim, int sdim);
FCALLSCSUB4(cf_util_reshape2, UTIL_RESHAPE2, util_reshape2,
	    PDOUBLE, PDOUBLE, INT, INT)

/****************************************************************************/

void cf_util_reshape(double *dest, double *source, int amount)
{
  size_t iamount;

  iamount = (size_t) amount * sizeof(double);

  memcpy (dest, source, iamount);

  return;
}

void cf_util_reshape2(double *dest, double *source, int ddim, int sdim)
{
  size_t iamount, izero;

  if ( ddim > sdim )
    {
      iamount =  sdim;    
      izero   =  ddim - sdim;
    }
  else
    {
      iamount =  ddim;    
      izero   =  0;
    }

  memcpy (dest, source, iamount*sizeof(double));

  if ( izero > 0 )
    memset (dest+iamount, 0, izero*sizeof(double));

  return;
}
