/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_date.h
 * 
 * @brief Date and some operations supported on Date.
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology.
 * @author  Rahul Sinha, Max Planck Institute for Meteorology.
 *
 * @date March 2013
 *
 */


#ifndef _MTIME_DATE_H
#define _MTIME_DATE_H

#include "mtime_calendar.h"

///provides a string length for toString.
#define MAX_DATE_STR_LEN 32

struct _timedelta;

/**
 * @struct _date
 *
 * @brief struct _date containing usual date parameters. 
 */
struct _date
{
  int64_t year; ///< Year of date. Can be both positive and negative.
  int month;	///< Month of date.	
  int day;	///< day of date.
};

struct _date*
newDate(const char* ds);

struct _date*
newRawDate(int64_t _year, int _month, int _day);

struct _date*
constructAndCopyDate(struct _date* d);

void
deallocateDate(struct _date* d);

/*! \cond PRIVATE */
struct _datetime *
convertDateToDateTime(struct _date* d, struct _datetime* dt_return);

struct _date *
convertDateTimeToDate(struct _datetime* dt, struct _date* d_return);

compare_return_val
compareDate(struct _date*, struct _date*);

struct _date*
replaceDate(struct _date*, struct _date*);
/*! \endcond  */

char*
dateToString(struct _date*, char* ds);

char *
dateToPosixString(struct _date* d, char* toStr);

/**
 * @}
 */ 
#endif
