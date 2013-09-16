/**
 * \addtogroup CBindings libmtime C language bindings
 * 
 * @file mtime_date.c
 * 
 * @brief Date and some operations supported on Date.
 *
 * @author Luis Kornblueh, Rahul Sinha. MPIM.
 * @date March 2013
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>

#include "mtime_date.h"

#include "mtime_julianDay.h"
#include "mtime_time.h"
#include "mtime_calendar.h"
#include "mtime_timedelta.h"
#include "mtime_datetime.h"
#include "mtime_iso8601.h"

/* MIN allowed year : -2147483648
   MAX allowed year :  2147483647 */


/**
 * @brief Construct new Date using an ISO 8601 conforming string. 
 *
 * MAX allowed year = 2147483647 and MIN allowed year = -2147483648.
 *
 * @param  ds
 *         A pointer to char. The string should contain parameters with which Date object is to be created.
 *	   A string literal can be accepted.
 * @return d
 *         A pointer to an initialized date object. 
 */

struct _date *
newDate(const char* ds)
{
  if ((ds != NULL) && (getCalendarType()))
    {
      /* Intialize a dummy ISO8601 object. */
      struct iso8601_datetime* isoDt = new_iso8601_datetime('+', 0, 0, 0, 0, 0, 0, 0, 'Z');
      if (isoDt == NULL)
        return NULL ;

      /* Check ISO8601 compliance. */
      if (verify_string_datetime(ds, isoDt) != DATETIME_MATCH)
        {
          deallocate_iso8601_datetime(isoDt);
          return NULL ;
        }

      struct _date *d = (struct _date*) calloc(1, sizeof(struct _date));
      if (d == NULL)
        {
          deallocate_iso8601_datetime(isoDt);
          return NULL ;
        }

      d->year = isoDt->year;
      d->month = isoDt->month;
      d->day = isoDt->day;

      /* Cleanup. */
      deallocate_iso8601_datetime(isoDt);

      return d;
    }
  else
    return NULL ;
}


/**
 * @brief Construct new Date using 'raw' numerical values.
 *
 * MAX allowed year = 2147483647 and MIN allowed year = -2147483648.
 *
 * @param  _year
 *         An "int64_t" value denoting the year part of Date. 
 * @param  _month
 *         An "int" value denoting the month part of Date.
 * @param  _day
 *         An "int" value denoting the day part of Date.
 *
 * @return d
 *         A pointer to an initialized Date object. 
 */

struct _date*
newRawDate(int64_t _year, int _month, int _day)
{
  char* ds = (char*) calloc(MAX_DATE_STR_LEN, sizeof(char));
  if (ds == NULL )
    {
      return NULL ;
    }

  snprintf(ds, MAX_DATE_STR_LEN, "%" PRIi64 "-%02d-%02d", _year, _month, _day);

  /* Resuse string interface and create object. */
  struct _date* d = newDate(ds);

  /* Cleanup. */
  free(ds);
  ds = NULL;

  return d;
}


/**
 * @brief Copy the values and construct a new date.
 *
 * @param  d
 *         A pointer to struct _date. Values of d are used to initialize the new date object being created. 
 *
 * @return _d
 *         A pointer to an initialized Date object. 
 */

struct _date*
constructAndCopyDate(struct _date* d)
{
  if( d != NULL )
    return newRawDate(d->year, d->month, d->day);
  else
    return NULL;
}


/**
 * @brief Destructor of Date. Free the Date object.
 *
 * @param  d
 *         A pointer to struct _date. d is deallocated.
 */

void
deallocateDate(struct _date* d)
{
  if (d != NULL )
    {
      free(d);
      d = NULL;
    }
}


/*
 * NOTE: Internal and not doxyfied.
 *
 * @brief Copy the contents of Date into a DateTime object.
 * 
 * Routine convertDateToDateTime copies the contents of date object into a datetime object 
 * and sets the time parameters to 0. For e.g A date 1999-12-01 can be copied into a DateTime
 * object resulting in 1999-12-01T00:00:00.000. This function can be used by Date objects
 * to get an equivalent DateTime object and then call routines supporting DateTime.
 *
 * @param  d
 *         A pointer to struct _date. Copy "FROM" date object.
 *
 * @param  dt_return
 *         A pointer to struct _datetime. Copy "TO" datetime object.
 *
 * @return dt_return
 *         A pointer to the copied datetime object.
 */

struct _datetime *
convertDateToDateTime(struct _date* d, struct _datetime* dt_return)
{
  if ((d != NULL) && (dt_return != NULL)){

    dt_return->date.year 	= d->year;
    dt_return->date.month 	= d->month;
    dt_return->date.day 	= d->day;

    /* Add dummy time. */
    dt_return->time.hour 	= 0;
    dt_return->time.minute 	= 0;
    dt_return->time.second 	= 0;
    dt_return->time.ms 		= 0;

    return dt_return;
  }
  else
    return NULL;
}

/*
 * Note: Internal and not doxyfied. 
 *
 * @brief Copy the contents of DateTime into a Date object.
 *
 * Routine convertDateTimeToDate copies the date elements of DateTime object(Year,Month,Day)
 * into a Date object and discards the time elements (Hour, minute, second, ms). For e.g 
 * A DateTime 1999-12-01T00:00:00.000 can be converted into a Date Object resulting in 
 * 1999-12-01. (This function is used to 'go-back-to' Date format! Date objects, which convert
 * to DateTime using the routine convertDateToDateTime can return to the original Date
 * format after relevant processing on the meta-DateTime object.)
 *
 * @param  dt
 *         A pointer to struct _datetime. Copy "FROM" datetime object.
 *
 * @param  d_return
 *         A pointer to struct _date. Copy "TO" date object.
 *
 * @return d_return
 *         A pointer to the copied date object.
 */

struct _date *
convertDateTimeToDate(struct _datetime* dt, struct _date* d_return)
{
  if ((dt != NULL) && (d_return != NULL)){

    /* Copy Date elements. Discard Time elements. */
    d_return->year = dt->date.year;
    d_return->month = dt->date.month;
    d_return->day = dt->date.day;

    return d_return;
  }
  else
    return NULL;
}


/**
 * @brief Compare two dates and return (d1 > d2) OR (d1 = d2) OR (d1 < d2).
 *
 * @param  d1
 *         A pointer to struct _date.
 *
 * @param  d2
 *         A pointer to struct _date.
 *
 * @return boolean
 *         if (d1 > d2), return 'greater_than'. If (d1 = d2), return 'equal_to'. If (d1 < d2), return less_than. Return compare_error indicating error.
 */

compare_return_val
compareDate(struct _date* d1, struct _date* d2)
{
  if ((d1 != NULL )&& (d2 != NULL) ){

    /* Initialize */
    compare_return_val boolean = compare_error;

    if (d1->year == d2->year)
      {
        if (d1->month == d2->month)
	  {
	    if (d1->day == d2->day)
	      {
	        boolean = equal_to;
	      }
	    else if (d1->day > d2->day)
	      {
	        boolean = greater_than;
	      }
	    else
	      {
	        boolean = less_than;
	      }
	  }
        else if (d1->month > d2->month)
	  {
	    boolean = greater_than;
	  }
        else
	  {
	    boolean = less_than;
  	  }

      }
    else if (d1->year > d2->year)
      {
        boolean = greater_than;
      }
    else
      {
        boolean = less_than;
      }

    return boolean;
  }
  else
    return compare_error;
}


/**
 * @brief COPY a Date object.
 *
 * Routine replaceDate copies the contents of source Date into a Destination Date object.
 *
 * @param  dsrc
 *         A pointer to struct _date. Copy "FROM" Date object.
 *
 * @param  ddest
 *         A pointer to struct _date. Copy "TO" Date object.
 *
 * @return ddest
 *         A pointer to 'copied' Date Object.
 */

struct _date*
replaceDate(struct _date* dsrc, struct _date* ddest)
{
  if ((ddest != NULL )&& ( dsrc!=NULL ) ){

    /* Copy. */
    ddest->year = dsrc->year;
    ddest->month = dsrc->month;
    ddest->day = dsrc->day;

    return ddest;
  }
  else
    return NULL;
}


/**
 * @brief Get Date as a string.
 *
 * DateToString returns a string in IS08601 compliant (and extended) format.
 *
 * @param  d
 *         A pointer to struct _date. The date to be converted to string.
 *
 * @param  toStr
 *         A pointer to char. String where date is to be written.
 *
 * @return toStr
 *         A pointer to the string containing date.
 */

char*
dateToString(struct _date* d, char* toStr)
{
  if ((d != NULL )&& ( toStr!=NULL ) ){

    memset(toStr,'\0',MAX_DATE_STR_LEN);
    snprintf(toStr,MAX_DATE_STR_LEN,"%" PRIi64 "-%02d-%02d", d->year, d->month, d->day);

    return toStr;
  }
  else
    return NULL;
}


/**
 * @brief Get Date in 'struct tm' format and return as a string.
 *
 * Only dates between and including 1582-10-15 TO 9999-12-31 supported.
 *
 * @param  d
 *         A pointer to struct _date. The date to be converted to string.
 *
 * @param  toStr
 *         A pointer to char. String where date is to be written.
 *
 * @return toStr
 *         A pointer to the string containing date.
 */

char *
dateToPosixString(struct _date* d, char* toStr)
{
  if ((d != NULL )&& ( toStr!=NULL ) ){

    struct tm tm_info;

    /* Range check. Return with NULL indicating Error */
    if ( d->year < POSIXSTRING_YEAR_LOWER_BOUND )
      {
        return NULL;
      }
    else if ( 
		((d->year == POSIXSTRING_YEAR_LOWER_BOUND) && (d->month < POSIXSTRING_MONTH_LOWER_BOUND)) 
		||  
		((d->year == POSIXSTRING_YEAR_LOWER_BOUND) && (d->month == POSIXSTRING_MONTH_LOWER_BOUND) && (d->day < POSIXSTRING_DAY_LOWER_BOUND) ) 
	    )
      { 
        return NULL;
      }
    else if ( d->year > POSIXSTRING_YEAR_UPPER_BOUND )
      {
        return NULL;
      }


    tm_info.tm_mday	= d->day;
    /* Range of month is from 0 to 11. */
    tm_info.tm_mon	= d->month - 1;
    /* tm's year is w.r.t the year 1900. */
    tm_info.tm_year	= d->year - 1900;

    memset(toStr,'\0',MAX_DATE_STR_LEN);
    strftime(toStr, MAX_DATE_STR_LEN , "%Y:%m:%d", &tm_info);

    return toStr;
  }
  else
    return NULL;
}
