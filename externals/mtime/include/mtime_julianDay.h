/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_julianDay.h
 *
 * @brief Julian Day Calendar and some operations supported on julian dates.
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology.
 * @author  Rahul Sinha, Max Planck Institute for Meteorology.
 * 
 * @date March 2013
 */

#ifndef _MTIME_JULIANDAY_H
#define _MTIME_JULIANDAY_H

#include <stdint.h>
#include <stdbool.h>

/// provides a string length for toString.
#define MAX_JULIANDAY_STR_LEN 32


/**
 * @struct _julianday
 *
 * @brief Struct _julianday containing julian day parameters.
 */
struct _julianday
{
  int64_t day; 	///< the actual Julian day.
  int64_t ms;	///< the milisecond on that particular day.
};

/*! \cond PRIVATE */
/* Notice that Julian delta, as such, is not defined by any standard. In this lib,
 it serves as a proxy for equivalent deltas on julian time scale. For example, a 
 TD of P02DT12H is equivalend to (2 days and 43200000 ms) Julian delta and -P02DT12H
 is equivalent to (-2 days and -43200000 ms) */
struct _juliandelta
{
  char sign;
  int64_t day;
  int64_t ms;
};
/*! \endcond */

struct _julianday*
newJulianDay(int64_t _day, int64_t _ms);

void
deallocateJulianDay(struct _julianday* jd);

/*! \cond PRIVATE */
struct _juliandelta*
newJulianDelta(char _sign, int64_t _day, int64_t _ms);

void
deallocateJulianDelta(struct _juliandelta* jd);

struct _julianday*
addJulianDelta(struct _julianday* jd1, struct _juliandelta* jd2, struct _julianday* jd);

struct _julianday*
substractJulianDelta(struct _julianday* jd1, struct _juliandelta* jd2, struct _julianday* jd);

struct _juliandelta*
substractJulianDay(struct _julianday* jd1, struct _julianday* jd2, struct _juliandelta* jd);
/*! \endcond */

char*
juliandayToString(struct _julianday* jd, char* toStr);

/**
 * @}
 */

#endif
