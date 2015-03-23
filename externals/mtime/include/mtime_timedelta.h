/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_timedelta.h
 * 
 * @brief TimeDelta and some operations supported on TimeDelta.
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology.
 * @author  Rahul Sinha, Max Planck Institute for Meteorology. 
 * 
 * @date March 2013
 * 
 */

#ifndef _MTIME_TIMEDELTA_H
#define _MTIME_TIMEDELTA_H

#include <stdint.h>
#include <stdbool.h>

///provides a string length for toString
#define MAX_TIMEDELTA_STR_LEN 32

struct _timedelta;
struct _datetime;
struct _date;
struct _julianday;
struct _juliandelta;

/**
 * @struct _timedelta
 *
 * @brief Struct _timedelta containing timedelta and sign of year parameter. 
 */
struct _timedelta
{
  char sign;	///< sign of time delta. Sign can be '+' or '-'.

  int64_t year;	///< Year part of timedelta.
  int month;	///< Month part of timedelta.
  int day;	///< Day part of timedelta.
  int hour;	///< Hour part of timedelta.
  int minute;	///< Minute part of timedelta.
  int second;	///< Second part of timedelta.
  int ms;	///< Milli-Second part of timedelta.
};

struct _timedelta*
newTimeDelta(const char* timedelta_string);

struct _timedelta*
newRawTimeDelta(char _sign, int64_t _year, int _month, int _day, int _hour, int _minute, int _second, int _ms);

struct _timedelta*
constructAndCopyTimeDelta(struct _timedelta* td);

void
deallocateTimeDelta(struct _timedelta* td);

/*! \cond PRIVATE */
struct _juliandelta*
timeDeltaToJulianDelta(struct _timedelta *td, struct _datetime *dt, struct _juliandelta *jd);

struct _timedelta*
julianDeltaToTimeDelta(struct _juliandelta* jd, struct _datetime *dt, struct _timedelta* td_return);
/*! \endcond */

struct _timedelta*
getTimeDeltaFromDate(struct _date*, struct _date*, struct _timedelta*);

struct _timedelta*
getTimeDeltaFromDateTime(struct _datetime* dt1, struct _datetime* dt2, struct _timedelta* td_return);

int64_t
getTotalMilliSecondsTimeDelta(struct _timedelta *td, struct _datetime* dt);

int64_t
getTotalSecondsTimeDelta(struct _timedelta *td, struct _datetime* dt);

char*
timedeltaToString(struct _timedelta * td, char *toString);

struct _date*
addTimeDeltaToDate(struct _date* dt, struct _timedelta* td, struct _date* dt_return);

struct _datetime*
addTimeDeltaToDateTime(struct _datetime* dt, struct _timedelta* td, struct _datetime* dt_return);

struct _timedelta*
moduloTimeDeltaFromDateTime(struct _datetime* start_dt, struct _timedelta* timestep, struct _datetime* current_dt, struct _timedelta* modulo_td);

struct _timedelta*
elementwiseScalarMultiplyTimeDelta(struct _timedelta* base_td, int64_t lambda, struct _timedelta* scaled_td);

struct _timedelta*
elementwiseAddTimeDeltatoTimeDelta(struct _timedelta* td1, struct _timedelta* td2,  struct _timedelta* td_return);

char*
getPTStringFromMS(int64_t _ms, char* PTstr);

char*
getPTStringFromSeconds(int64_t _s, char* PTstr);

char*
getPTStringFromMinutes(int64_t _m, char* PTstr);

char*
getPTStringFromHours(int64_t _h, char* PTstr);

/*! \cond PRIVATE */
bool
testYearIsLeapYear(int64_t year);
/*! \endcond */

/**
 * @}
 */
#endif
