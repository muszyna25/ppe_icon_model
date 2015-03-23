/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_datetime.h
 * 
 * @brief DateTime and some operations supported on DateTime.
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology.
 * @author  Rahul Sinha, Max Planck Institute for Meteorology.
 * 
 * @date March 2013
 *
 */


#ifndef _MTIME_DATETIME_H
#define _MTIME_DATETIME_H

#include <stdint.h>

#include "mtime_date.h"
#include "mtime_date.h"
#include "mtime_time.h"
#include "mtime_calendar.h"

///provides a string length for toString.
#define MAX_DATETIME_STR_LEN 32


/**
 * @struct _datetime
 *
 * @brief struct _datetime contains a struct _date and a struct _time element. 
 */
struct _datetime
{
  struct _date date;	///< Date elements.
  struct _time time;	///< Time elements.
};

struct _datetime*
newDateTime(const char* datetime_string);

struct _datetime*
newRawDateTime(int64_t _year, int _month, int _day, int _hour, int _minute, int _second, int _ms);

struct _datetime*
constructAndCopyDateTime(struct _datetime* dt);

void
deallocateDateTime(struct _datetime* dt);

compare_return_val
compareDatetime(struct _datetime* dt1, struct _datetime* dt2);

struct _datetime*
replaceDatetime(struct _datetime* dtsrc, struct _datetime* dtdest);

char*
datetimeToString(struct _datetime* dt, char* toStr);

char *
datetimeToPosixString(struct _datetime* dt, char* toStr);

int
getNoOfDaysInMonthDateTime(struct _datetime* dt);

int
getNoOfDaysInYearDateTime(struct _datetime* dt);

int64_t
getNoOfSecondsElapsedInMonthDateTime(struct _datetime* dt);

int64_t
getNoOfSecondsElapsedInDayDateTime(struct _datetime* dt);

int
getDayOfYearFromDateTime(struct _datetime* currentdt);

struct _julianday*
getJulianDayFromDateTime(struct _datetime* dt, struct _julianday* jd);

double
datetimeDivideBySeconds(struct  _datetime* refdt, struct  _datetime* dt, double intvlsec);

struct _datetime*
datetimeAddSeconds(struct  _datetime* refdt, double intvlsec, struct _datetime* dt_return);

/**
 * @}
 */

#endif
