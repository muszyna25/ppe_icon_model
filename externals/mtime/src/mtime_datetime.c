/**
 * @file mtime_datetime.c
 * 
 * @brief DateTime and some operations supported on DateTime.
 *
 * @author Luis Kornblueh, Rahul Sinha. MPIM.
 * @date March 2013
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>

#include "mtime_datetime.h"

#include "mtime_calendar.h"
#include "mtime_julianDay.h"
#include "mtime_timedelta.h"
#include "mtime_date.h"
#include "mtime_time.h"
#include "mtime_iso8601.h"


/* MIN allowed year : -2147483648
   MAX allowed year :  2147483647
   Upto MilliSecond resolution supported. */


/**
 * @brief Construct new DateTime using a string.
 *
 * @param  dts
 *         A pointer to char. The string should contain parameters with which DateTime is created.
 *	   A string literal can be accepted.
 *
 * @return dt
 *         A pointer to a filled datetime. 
 */

struct _datetime*
newDateTime(const char* dts)
{
  if ((dts != NULL) && (getCalendarType()))
    {
      /* Intialize a dummy ISO8601 object. */
      struct iso8601_datetime* isoDt = new_iso8601_datetime('+', 0, 0, 0, 0, 0, 0, 0, 'Z');
      if (isoDt == NULL)
        return NULL ;

      /* Check ISO8601 compliance. */
      if (verify_string_datetime(dts, isoDt) != DATETIME_MATCH)
        {
          deallocate_iso8601_datetime(isoDt);
          return NULL ;
        }

      struct _datetime *dt = (struct _datetime*)calloc(1,sizeof(struct _datetime));
      if (dt == NULL )
        {
          deallocate_iso8601_datetime(isoDt);
          return NULL ;
        }

      dt->date.year 	= isoDt->year;
      dt->date.month 	= isoDt->month;
      dt->date.day 	= isoDt->day;
      dt->time.hour 	= isoDt->hour;
      dt->time.minute	= isoDt->minute;
      dt->time.second	= isoDt->second;
      dt->time.ms	= isoDt->ms;

      //Cleanup.	
      deallocate_iso8601_datetime(isoDt);
    
      return dt;
    }
  else
    return NULL ;
}


/**
 * @brief Construct new DateTime using 'raw' numerical values.
 *
 * @param  _year
 *         An "int64_t" value denoting the year part of DateTime. 
 * @param  _month
 *         An "int" value denoting the month part of DateTime.
 * @param  _day
 *         An "int" value denoting the day part of DateTime.
 * @param  _hour
 *         An "int" value denoting the hour part of DateTime.
 * @param  _minute
 *         An "int" value denoting the minute part of DateTime.
 * @param  _second
 *         An "int" value denoting the second part of DateTime.
 * @param  _ms
 *         An "int" value denoting the milli-second part of DateTime.
 *
 * @return dt
 *         A pointer to a filled DateTime. 
 */

struct _datetime*
newRawDateTime(int64_t _year, int _month, int _day, int _hour, int _minute, int _second, int _ms)
{
  char* dts = (char*)calloc(MAX_DATETIME_STR_LEN, sizeof(char));
  if (dts == NULL )
    {
      return NULL ;
    }

  snprintf(dts, MAX_DATETIME_STR_LEN, "%" PRIi64 "-%02d-%02dT%02d:%02d:%02d.%03d", _year, _month, _day, _hour, _minute, _second,
      _ms);

  /* Resuse string interface and create object. */
  struct _datetime* dt = newDateTime(dts);

  free(dts);
  dts = NULL;

  return dt;
}


/**
 * @brief Copy the values and construct a new datetime.
 *
 * @param  dt
 *         A pointer to struct _datetime. Values of dt are used to initialize the new datetime being created. 
 *
 * @return _dt
 *         A pointer to an initialized DateTime. 
 */

struct _datetime*
constructAndCopyDateTime(struct _datetime* dt)
{
  if ( dt != NULL )
    return newRawDateTime(dt->date.year, dt->date.month, dt->date.day, dt->time.hour, dt->time.minute, dt->time.second, dt->time.ms);
  else
    return NULL;
}


/**
 * @brief Destructor of DateTime.
 *
 * @param  dt
 *         A pointer to struct _datetime. dt is deallocated.
 */

void
deallocateDateTime(struct _datetime* dt)
{
  if (dt != NULL )
    {
      free(dt);
      dt = NULL;
    }
}


/**
 * @brief Get the TimeDelta between two DateTimes dt1 and dt2 as (dt1-dt2).
 *
 * Routine getTimeDeltaFromDateTime 'substracts' two DateTime's and returns the TimeDelta between
 * them. Each datetime is converted to an equivalent Julian Date. Substraction is then performed
 * on Julian axis. The "Julian delta" is finally converted back to normal calendar delta. 
 * 
 * This routine handles all supported Calendar types; i.e. the translation from Calendar date 
 * to Julian date and conversion from Julian Delta to normal TimeDetla is Calendar-type dependent. 
 * For eg. for Calendar type Gregorian, the TimeDelta between 2001-02-01T00:00:00.000 and 
 * 2001-01-01T00:00:00.000 will be 1 month. Similarly, for Calendar of type 360-Day-Calendar, 
 * the TimeDelta will be 1 month. It must be noted however, that the two dates differ by 31 and 
 * 30 days respectively. 
 *
 * @param  dt1
 *         A pointer to struct _datetime. 
 *
 * @param  dt2
 *         A pointer to struct _datetime.
 *
 * @param  td_return
 *         A pointer to struct _timedelta. Copy the result of (dt1 - dt2) in td_return.
 *
 * @return td_return
 *	  A pointer to TimeDelta containing the result of substraction.
 */

struct _timedelta*
getTimeDeltaFromDateTime(struct _datetime* dt1, struct _datetime* dt2, struct _timedelta* td_return)
{
  if ((dt1 != NULL )&& (dt2 != NULL) && (td_return != NULL) ){
 
  /* Convert dt1 to Julian. */
  struct _julianday* jd1 = newJulianDay(0, 0);
  if (jd1 == NULL)
    return NULL;
  jd1 = date2julian(dt1, jd1);

  /* Convert dt2 to Julian. */
  struct _julianday* jd2 = newJulianDay(0, 0);
  if (jd2 == NULL)
    {
      deallocateJulianDay(jd1);
      return NULL;
    }
  jd2 = date2julian(dt2, jd2);

  /* Calculat Delta on Julian axis. "RULE: A - B = Delta". If A > B, Delta is positive. */
  struct _juliandelta* jd = newJulianDelta('+', 0, 0);
  if (jd == NULL)
    {
      deallocateJulianDay(jd1);
      deallocateJulianDay(jd2);
      return NULL;
    }

  /* Substract the 2 dates on julian axis. */
  jd = substractJulianDay(jd1, jd2, jd);

  /* Convert Julian-delta to TimeDelta. */
  td_return = julianDeltaToTimeDelta(jd, dt1, td_return);

  /* Cleanup. */
  deallocateJulianDay(jd1);
  deallocateJulianDay(jd2);
  deallocateJulianDelta(jd);

  return td_return;
}
else
  return NULL;
}


/**
 * @brief Compare two datetimes and return (dt1 > dt2) OR (dt1 = dt2) OR (dt1 < dt2).
 *
 * @param  dt1
 *         A pointer to struct _datetime.
 *
 * @param  dt2
 *         A pointer to struct _datetime.
 *
 * @return boolean
 *         if (dt1 > dt2), return greater_than. If (dt1 == dt2), return equal_to. If (dt1 < dt2), return less_than. 
 *	   Return compare_error indicating error.
 */

compare_return_val
compareDatetime(struct _datetime* dt1, struct _datetime* dt2)
{
  if ((dt1 != NULL) && (dt2 != NULL)){

  /* Initialize. */
  compare_return_val boolean = compare_error;

  if (dt1->date.year == dt2->date.year)
    {
      if (dt1->date.month == dt2->date.month)
        {
          if (dt1->date.day == dt2->date.day)
            {
              if (dt1->time.hour == dt2->time.hour)
                {
                  if (dt1->time.minute == dt2->time.minute)
                    {
                      if (dt1->time.second == dt2->time.second)
                        {
                          if (dt1->time.ms == dt2->time.ms)
                            {
                              boolean = equal_to;
                            }
                          else if (dt1->time.ms > dt2->time.ms)
                            {
                              boolean = greater_than;
                            }
                          else
                            {
                              boolean = less_than;
                            }

                        }
                      else if (dt1->time.second > dt2->time.second)
                        {
                          boolean = greater_than;
                        }
                      else
                        {
                          boolean = less_than;
                        }
                    }
                  else if (dt1->time.minute > dt2->time.minute)
                    {
                      boolean = greater_than;
                    }
                  else
                    {
                      boolean = less_than;
                    }
                }
              else if (dt1->time.hour > dt2->time.hour)
                {
                  boolean = greater_than;
                }
              else
                {
                  boolean = less_than;
                }

            }
          else if (dt1->date.day > dt2->date.day)
            {
              boolean = greater_than;
            }
          else
            {
              boolean = less_than;
            }
        }
      else if (dt1->date.month > dt2->date.month)
        {
          boolean = greater_than;
        }
      else
        {
          boolean = less_than;
        }

    }
  else if (dt1->date.year > dt2->date.year)
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
 * @brief COPY a DateTime object.
 *
 * Routine replaceDateTime copies the contents of source DateTime into a Destination DateTime object.
 *
 * @param  dtsrc         
 *         A pointer to struct _datetime. Copy "FROM" DateTime object.
 *
 * @param  dtdest
 *         A pointer to struct _datetime. Copy "TO" DateTime object.
 *
 * @return dtdest
 *         A pointer to 'copied' DateTime Object.
 */

struct _datetime*
replaceDatetime(struct _datetime* dtsrc, struct _datetime* dtdest)
{
  if ((dtdest != NULL) && (dtsrc != NULL)){

  dtdest->date.year 	= dtsrc->date.year;
  dtdest->date.month 	= dtsrc->date.month;
  dtdest->date.day 	= dtsrc->date.day;
  dtdest->time.hour 	= dtsrc->time.hour;
  dtdest->time.minute 	= dtsrc->time.minute;
  dtdest->time.second 	= dtsrc->time.second;
  dtdest->time.ms 	= dtsrc->time.ms;

  return dtdest;
}
else
  return NULL;
}


/**
 * @brief Get the 'day-of-year' value of a DateTime.
 *
 * Routine getDayOfYearFromDateTime returns Day of Year for the DateTime. This routine supports
 * all Calendar types.
 *
 * For eg. the day of year value for 2001-10-15T00:00:00.000 will be 288 for Gregorian Calendar. 
 * Similarly, this value will be 285 for Calendar of type 360 day-Calendar.
 *
 * @param  dt
 *         A pointer to struct _datetime. Retrieve the 'day-of-year' from this DT object.
 *
 * @return doy
 *         Integer value of doy. The value depends on the calendar type. Zero indicates error.
 */

int
getDayOfYearFromDateTime(struct _datetime* dt)
{
  if (dt != NULL )
    {
      int doy = 0;

      /* Get current Calendar type. */
      calendarType ct = getCalendarType();

      if (ct == PROLEPTIC_GREGORIAN)
        {
          if (testYearIsLeapYear(dt->date.year))
            {
	      /* Leap year. */
              doy = nofDaysAfterARGMonthsInLeapYear[dt->date.month - 1] + dt->date.day;
            }
          else
            {
	      /* Non-leap year. */
              doy = nofDaysAfterARGMonthsInNonLeapYear[dt->date.month - 1] + dt->date.day;
            }

        }
      else if (ct == YEAR_OF_365_DAYS)
        {
	  /* Non-leap year characteristics. */
          doy = nofDaysAfterARGMonthsInNonLeapYear[dt->date.month - 1] + dt->date.day;
        }
      else if (ct == YEAR_OF_360_DAYS)
        {
	  /* Each month has 30 days. */
          doy = (dt->date.month - 1) * NO_OF_DAYS_IN_A_MONTH_FOR_CAL_TYPE360 + dt->date.day;
        }

      return doy;
    }
  else
    return 0;
}


/**
 * @brief Get nod (number of Days) in the month of DateTime.
 *
 * Routine getNoOfDaysInMonthDateTime returns number of days in the month of DateTime. This routine  
 * supports all calendar types.
 *
 * For eg. the number of days for 2001-10-15T00:00:00.000 will be 31 for Gregorian Calendar. 
 * Similarly, this value will be 30 for Calendar of type 360 day-Calendar.
 *
 * @param  dt
 *         A pointer to struct _datetime.
 *
 * @return nod
 *         Integer value of nod. The value depends on the month and the calendar type. Zero indicates error.
 */
//TODO on Luis: Is this function doing the right thing?

int
getNoOfDaysInMonthDateTime(struct _datetime* dt)
{
  if (dt != NULL )
    {
      int nod = 0;

      /* Get current Calendar type. */
      calendarType ct = getCalendarType();

      if (ct == PROLEPTIC_GREGORIAN)
        {
          if (testYearIsLeapYear(dt->date.year))
            {
	      /* Leap year. */
              nod = nofDaysAfterARGMonthsInLeapYear[dt->date.month] - nofDaysAfterARGMonthsInLeapYear[dt->date.month - 1];
            }
          else
            {
	      /* Non leap year. */
              nod = nofDaysAfterARGMonthsInNonLeapYear[dt->date.month] - nofDaysAfterARGMonthsInNonLeapYear[dt->date.month - 1];
            }

        }
      else if (ct == YEAR_OF_365_DAYS)
        {
	  /* Non-leap year. */
          nod = nofDaysAfterARGMonthsInNonLeapYear[dt->date.month] - nofDaysAfterARGMonthsInNonLeapYear[dt->date.month - 1];
        }
      else if (ct == YEAR_OF_360_DAYS)
        {
	  /* Each month has 30 days. */
          nod = NO_OF_DAYS_IN_A_MONTH_FOR_CAL_TYPE360;
        }

      return nod;
    }
  else
    return 0;
}


/**
 * @brief Get number of days in the Year of DateTime.
 *
 * Routine getNoOfDaysInYearDateTime returns number of days in the Year of DateTime. This routine  
 * supports all calendar types.
 *
 * Number of days returned will depend on the calendar type and if applicable, leap v/s non leap year.
 *
 * @param  dt
 *         A pointer to struct _datetime.
 *
 * @return nod
 *         Integer value of nod. The value depends on the year and the calendar type. Zero indicates error.
 */
//TODO on Luis: Is this function doing the right thing?
int
getNoOfDaysInYearDateTime(struct _datetime* dt)
{
  if (dt != NULL )
    {
      int nod = 0;

      /* Get current Calendar type. */
      calendarType ct = getCalendarType();

      if (ct == PROLEPTIC_GREGORIAN)
        {
          if (testYearIsLeapYear(dt->date.year))
            {
	      /* Leap year. */
              nod = NO_OF_DAYS_IN_A_LEAP_YEAR;
            }
          else
            {
	      /* Non leap year. */
              nod = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
            }
        }
      else if (ct == YEAR_OF_365_DAYS)
        {
	  /* Non Leap year. */
          nod = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
        }
      else if (ct == YEAR_OF_360_DAYS)
        {
	  /* Year has 360 days. */
          nod = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE360;
        }

      return nod;
    }
  else
    return 0;
}


/**
 * @brief Get number of seconds elapsed in the month of DateTime.
 *
 * @param  dt
 *         A pointer to struct _datetime.
 *
 * @return no_of_seconds
 *         int64_t value of no_of_seconds. -1 indicates error.
 */

int64_t
getNoOfSecondsElapsedInMonthDateTime(struct _datetime* dt)
{
  if ( dt != NULL )
    return (	(dt->date.day-1) 	* (NO_OF_MS_IN_A_DAY / NO_OF_MS_IN_A_SECOND) 
		+ dt->time.hour 	* (NO_OF_MS_IN_A_HOUR / NO_OF_MS_IN_A_SECOND) 
		+ dt->time.minute 	* (NO_OF_MS_IN_A_MINUTE / NO_OF_MS_IN_A_SECOND) 
		+ dt->time.second
	   );
  else
    return -1;
}


/**
 * @brief Get number of seconds elapsed in the day of DateTime.
 *
 * @param  dt
 *         A pointer to struct _datetime.
 *
 * @return no_of_seconds
 *         int64_t value of no_of_seconds. -1 indicates error.
 */

int64_t
getNoOfSecondsElapsedInDayDateTime(struct _datetime* dt)
{
  if ( dt != NULL )
    return  (	  dt->time.hour 	* (NO_OF_MS_IN_A_HOUR / NO_OF_MS_IN_A_SECOND) 
		+ dt->time.minute 	* (NO_OF_MS_IN_A_MINUTE / NO_OF_MS_IN_A_SECOND) 
		+ dt->time.second
	    );
  else
    return -1;
}


/**
 * @brief Get the Julian Day from DateTime.
 *
 * The routine getJulianDayFromDateTime returns the equivalent Julian date to DateTime. Internally
 * it calls translation routines based on Calndar type.
 *
 * @param  dt
 *         A pointer to struct _datetime. The DT's value is converted to julian day value.
 *
 ** @param  jd
 *         A pointer to struct _julianday. JD where the converted value is stored.
 *
 * @return jd
 *         A pointer to struct _julianday containing a copy of the JD value corresponding to the DT.
 */

struct _julianday*
getJulianDayFromDateTime(struct _datetime* dt, struct _julianday* jd)
{
  /* Function pointer date2julian points to the correct translation routine. */
  return date2julian(dt,jd);
}


/**
 * @brief Get DateTime as a string.
 *
 * datetimeToString returns a string in IS08601 compliant (and extended) format.
 *
 * @param  dt
 *         A pointer to struct _datetime. The datetime to be converted to string.
 *
 * @param  toStr
 *         A pointer to char. String where datetime is to be written.
 *
 * @return toStr
 *         A pointer to the string containing datetime.
 */

char*
datetimeToString(struct _datetime* dt, char* toStr)
{
  if ((dt != NULL )&& (toStr != NULL)){

    memset(toStr,'\0',MAX_DATETIME_STR_LEN);

    snprintf(toStr,MAX_DATETIME_STR_LEN,"%" PRIi64 "-%02d-%02dT%02d:%02d:%02d.%03dZ",
      dt->date.year, dt->date.month, dt->date.day,
      dt->time.hour, dt->time.minute, dt->time.second, dt->time.ms);

  return toStr;
}
else
  return NULL;
}


/**
 * @brief Get DateTime in 'struct tm' format and return as a string.
 *
 * @param  dt
 *         A pointer to struct _datetime. The datetime to be converted to string.
 *
 * @param  toStr
 *         A pointer to char. String where datetime is to be written.
 *
 * @return toStr
 *         A pointer to the string containing datetime.
 */

char *
datetimeToPosixString(struct _datetime* dt, char* toStr)
{
  if ((dt != NULL ) && (toStr != NULL))
  {
    struct tm tm_info;

    /*Range check. Return with NULL indicating Error */
    if ( dt->date.year < POSIXSTRING_YEAR_LOWER_BOUND)
      {
        return NULL;
      }
    else if ( 
		(	(dt->date.year == POSIXSTRING_YEAR_LOWER_BOUND) 
			&& 
			(dt->date.month < POSIXSTRING_MONTH_LOWER_BOUND)
		) 
		||  
		(	(dt->date.year == POSIXSTRING_YEAR_LOWER_BOUND) 
			&& 
			(dt->date.month == POSIXSTRING_MONTH_LOWER_BOUND) 
			&& 	
			(dt->date.day <POSIXSTRING_DAY_LOWER_BOUND) 
		) 
	    )
      {
        return NULL;
      }
    else if ( dt->date.year > POSIXSTRING_YEAR_UPPER_BOUND )
      {
        return NULL;
      }

  
    tm_info.tm_sec       = dt->time.second;
    tm_info.tm_min       = dt->time.minute;
    tm_info.tm_hour       = dt->time.hour;

    tm_info.tm_mday       = dt->date.day;
    /* Range of month is from 0 to 11. */
    tm_info.tm_mon        = dt->date.month - 1;
    /* tm's year is w.r.t the year 1900. */
    tm_info.tm_year       = dt->date.year - 1900;

    memset(toStr,'\0',MAX_DATETIME_STR_LEN);
    strftime(toStr, MAX_DATETIME_STR_LEN, "%Y:%m:%d %H:%M:%S", &tm_info);
  
    return toStr;
  }
  else
    return NULL;
}
