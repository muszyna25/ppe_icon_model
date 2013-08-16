/**
 * @file mtime_timedelta.c
 * 
 * @brief TimeDelta and some operations supported on TimeDelta.
 * 
 * @author Luis Kornblueh, Rahul Sinha. MPIM.
 * @date March 2013
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>

#include "mtime_timedelta.h"

#include "mtime_calendar.h"
#include "mtime_julianDay.h"
#include "mtime_datetime.h"
#include "mtime_date.h"
#include "mtime_time.h"
#include "mtime_iso8601.h"



#define NO_OF_MS_IN_A_DAY 86400000
#define NO_OF_MS_IN_HALF_DAY 43200000
#define NO_OF_MS_IN_A_HOUR 3600000
#define NO_OF_MS_IN_A_MINUTE 60000
#define NO_OF_MS_IN_A_SECOND 1000
#define NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE360 360
#define NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365 365
#define NO_OF_DAYS_IN_A_LEAP_YEAR 366
#define NO_OF_DAYS_IN_A_MONTH_FOR_CAL_TYPE360 30
#define NO_OF_MONTHS_IN_A_YEAR 12
#define NO_OF_HOURS_IN_A_DAY 24

// MIN allowed year : -2147483648
// MAX allowed year :  2147483647
// Upto MilliSecond resolution supported.


/**
 * @brief Construct new TimeDelta using an ISO 8601 conforming string.
 *
 * @param  tds
 *         A pointer to char. The string should contain parameters with which TimeDelta is to be created.
 *
 * @return td
 *         A pointer to a filled TimeDelta. 
 */
struct _timedelta*
newTimeDelta(const char* tds)
{
  if (tds != NULL )
    {
      /* Verify string is ISO8601 compliant. */
      struct iso8601_duration* isoDuration = new_iso8601_duration('+', 0, 0, 0, 0, 0, 0, 0);
      if (isoDuration == NULL )
        return NULL ;

      if (verify_string_duration(tds, isoDuration) != DURATION_MATCH)
        {
          deallocate_iso8601_duration(isoDuration);
          return NULL ;
        }

      /* Create TimeDelta object. */
      struct _timedelta* td = (struct _timedelta *)calloc(1,sizeof(struct _timedelta));
      if (td == NULL )
        {
          deallocate_iso8601_duration(isoDuration);
          return NULL ;
        }

      /* IMPORTANT: Negative/Positive time delta is indicated using sign (-/+). year,month,day..etc are always positive integers or 0. */
      td->sign 	= isoDuration->sign;
      td->year 	= isoDuration->year;
      td->month = isoDuration->month;
      td->day 	= isoDuration->day;
      td->hour 	= isoDuration->hour;
      td->minute = isoDuration->minute;
      td->second = isoDuration->second;
      td->ms 	= isoDuration->ms;

      //Cleanup.
      deallocate_iso8601_duration(isoDuration);

      return td;
    }
  else
    return NULL ;
}

/**
 * @brief Construct new TimeDelta using 'raw' numerical values.
 *
 * @param  _sign
 *         An char value denoting positive('+') or negative('-') TimeDelta. 
 * @param  _year
 *         An "int64_t" value denoting the year part of TimeDelta.
 * @param  _month
 *         An "int" value denoting the month part of TimeDelta.
 * @param  _day
 *         An "int" value denoting the day part of TimeDelta.
 * @param  _hour
 *         An "int" value denoting the hour part of TimeDelta.
 * @param  _minute
 *         An "int" value denoting the minute part of TimeDelta.
 * @param  _second
 *         An "int" value denoting the second part of TimeDelta.
 * @param  _ms
 *         An "int" value denoting the milli-second part of TimeDelta.
 *
 * @return td
 *         A pointer to a filled TimeDelta. 
 */

struct _timedelta*
newRawTimeDelta(char _sign, int64_t _year, int _month, int _day, int _hour, int _minute, int _second, int _ms)
{
  char* tds = (char*)calloc(MAX_TIMEDELTA_STR_LEN,sizeof(char));
  if (tds == NULL )
    {
      return NULL ;
    }

  snprintf(tds, MAX_TIMEDELTA_STR_LEN, "%cP%" PRIi64 "Y%02dM%02dDT%02dH%02dM%02d.%03dS", _sign, _year, _month, _day, _hour, _minute,
      _second, _ms);

  /* Use existing string interface. */
  struct _timedelta* td = newTimeDelta(tds);

  free(tds);
  tds = NULL;

  return td;
}



/**
 * @brief Copy the values and construct a new TimeDelta.
 *
 * @param  td
 *         A pointer to struct _timedelta. Values of td are used to initialize the new timedelta being created. 
 *
 * @return _td
 *         A pointer to an initialized TimeDelta object. 
 */

struct _timedelta*
constructAndCopyTimeDelta(struct _timedelta* td)
{
  if ( td != NULL )
    return newRawTimeDelta(td->sign, td->year, td->month, td->day, td->hour, td->minute, td->second, td->ms);
  else
    return NULL;
}


/**
 * @brief Destructor of TimeDelta.
 *
 * @param  td
 *         A pointer to struct _timedelta. td is deallocated.
 */

void
deallocateTimeDelta(struct _timedelta* td)
{
  if (td != NULL )
    {
      free(td);
      td = NULL;
    }
}


/*! \cond PRIVATE */
/* Internal function. Test is year is a leap year. */
bool
testYearIsLeapYear(int64_t year)
{
  bool flag = false;

  if (!(year % 400))
    {
      flag = true;
    }
  else if (!(year % 100))
    {
      flag = false;
    }
  else if (!(year % 4))
    {
      flag = true;
    }
  else
    {
      flag = false;
    }

  return flag;
}


/* Internal function. */
/* Converts TimeDelta value to a lib defined juliandelta value. Juliadelta depends on the calendar type. 
 Notice that TimeDelta is not uniquely defined but depends on the definition of corresponding DateTime.

 The library assumes the following definition: Let A denote an anchor date and P a timedelta. For a 
 positive P, A + P = B where date B > A. Consequently, for a negative P, A + P = B where A > B. 
 Also, when P is positive, a delta of 1 month has as many days as in the month of anchor DateTime;
 a delta of 2 months corresponds to the number of days in the anchor date month and the next month and 
 so on. When P is negative, a delta of 1 month corresponds to as many days as in the month before the
 anchor date month; a delta of 2 month corresponds to as many days as in the month before the
 anchor date month and the month before that and so on. 
 
 For eg, TimeDelta of P01M, when the 'Anchor' DateTime is 2001-02-01T00:00:00.000 is equivalent to a
 Julian Delta of positive 28 days ( #days in February ) and 0 ms. Also, TimeDelta of P02M, when the 
 'Anchor' DateTime is 2001-02-01T00:00:00.000 is equivalent to a Julian Delta of positive 28+31 days
 ( #days in February PLUS #days in March) and 0 ms. Similarly, a TimeDelta of -P01M with 'Anchor' 
 DateTime of 2001-02-01T00:00:00.000 is equivalent to a Julian Delta of negative 31 days 
 ( #days in January ) and 0 ms. Likewise, TimeDelta of -P02M with 'Anchor' DateTime of 
 2001-02-01T00:00:00.000 is equivalent to a Julian Delta of negative 31+31 days ( #days in January
 PLUS #days in December) and 0 ms.
 */

struct _juliandelta*
timeDeltaToJulianDelta(struct _timedelta* td, struct _datetime* base_dt, struct _juliandelta* jd_return)
{
  if ((td != NULL) && (base_dt != NULL) && (jd_return != NULL)){

  /* Reset initial delta to 0.*/
  jd_return->day = 0;
  jd_return->ms = 0;

  /* No of days in a year. */
  int ndiny;

  /* Pointer to an array of month_specific_delta_in_months. */
  const int (*msdinm)[NO_OF_MONTHS_IN_A_YEAR+1];

  /* Set parameters according to calendar type. */
  switch (getCalendarType())
    {
      case YEAR_OF_365_DAYS:

      msdinm = month_specific_delta_in_months_365;
      ndiny = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
      break;

      case YEAR_OF_360_DAYS:

      msdinm = month_specific_delta_in_months_360;
      ndiny = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE360;
      break;

      case PROLEPTIC_GREGORIAN:
      /* Handle Gregorian here. */

      /* Gregorian will have 366 days and 365 days depending on Leap year. */

      if ( td->sign == '+' )
        {
          jd_return->sign = '+';

          /* In final year. The crucial point is the month of february. */
          if (
              ( (testYearIsLeapYear(base_dt->date.year + td->year + 1) ) && (base_dt->date.month >= 3) )
              ||
              ( (testYearIsLeapYear(base_dt->date.year + td->year) ) && (base_dt->date.month < 3) )
          )
            {
              /* If the base year + delta year is a year before a leap year and base month is >= 3 
               OR 
               base year + delta year is a leap year and month is < 3
               => An addition of leap-year specific delta for each month.
               */
              msdinm = month_specific_delta_in_months_leapyear;
              ndiny = NO_OF_DAYS_IN_A_LEAP_YEAR;
            }
          else
            {
              /*Otherwise addition of non-leap year specific month.*/
              msdinm = month_specific_delta_in_months_365;
              ndiny = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
            }

          int i = 0;
          /* The year from base date to base_date + delta years -1 */
          for (i = base_dt->date.year; i < base_dt->date.year + td->year; i++)
            {
              if (
                  ((testYearIsLeapYear(i + 1)) && (base_dt->date.month >= 3))
                  ||
                  ((testYearIsLeapYear(i)) && (base_dt->date.month < 3))
              )
                {
                  /* If the next year is a leap year and month is >= 3 OR 
                   this year is a leap year and month is less than 3 
                   => delta of 1 year corresponds to 366 day julian delta. 
                   */
                  jd_return->day = jd_return->day + NO_OF_DAYS_IN_A_LEAP_YEAR;
                }
              else
                {
                  /* Otherwise. */
                  jd_return->day = jd_return->day + NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
                }
            }

          jd_return->day = jd_return->day + (msdinm[base_dt->date.month - 1][td->month]);
          jd_return->day = jd_return->day + td->day;
          jd_return->ms = td->hour * NO_OF_MS_IN_A_HOUR + td->minute * NO_OF_MS_IN_A_MINUTE + td->second * NO_OF_MS_IN_A_SECOND + td->ms;
        }
      else if ( td->sign == '-' )
        {
          jd_return->sign = '-';

          /* In final year. The crucial point is the month of february. */
          if (
              ((testYearIsLeapYear(base_dt->date.year - td->year - 1)) && (base_dt->date.month < 3))
              ||
              ((testYearIsLeapYear(base_dt->date.year - td->year)) && (base_dt->date.month >= 3)))
            {
              /* If the base year - delta year is a year after leap year and base month is < 3 
               OR 
               base year - delta year is a leap year and month is >= 3
               => A substraction of leap-year specific delta for each month.
               */
              msdinm = month_specific_delta_in_months_leapyear;
              ndiny = NO_OF_DAYS_IN_A_LEAP_YEAR;
            }
          else
            {
              /* Otherwise. */
              msdinm = month_specific_delta_in_months_365;
              ndiny = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
            }

          int i = 0;
          /* The year from base date to base_date - delta years + 1 */
          for (i = base_dt->date.year; i > base_dt->date.year - td->year; i--)
            {
              if (
                  ( (testYearIsLeapYear(i - 1)) && (base_dt->date.month < 3) )
                  ||
                  ( (testYearIsLeapYear(i)) && (base_dt->date.month >= 3) )
              )
                {
                  /* If the previous year is a leap year and month is < 3 OR 
                   this year is a leap year and month is >= 3 
                   => delta of 1 year corresponds to 366 day julian delta. 
                   */
                  jd_return->day = jd_return->day - NO_OF_DAYS_IN_A_LEAP_YEAR;
                }
              else
                {
                  /* Otherwise. */
                  jd_return->day = jd_return->day - NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
                }
            }
          jd_return->day = jd_return->day - (ndiny - msdinm[base_dt->date.month - 1][NO_OF_MONTHS_IN_A_YEAR - td->month]);
          jd_return->day = jd_return->day - td->day;
          jd_return->ms = (-1)*td->hour * NO_OF_MS_IN_A_HOUR - td->minute * NO_OF_MS_IN_A_MINUTE - td->second * NO_OF_MS_IN_A_SECOND - td->ms;

        }
      else
        return NULL; /* ERROR: Sign not set. */
      //break;
      return jd_return;

      default:
      return NULL /* Calendar type not set. */;
    }

  /* Handle 360 and 365 calendar type here as they are similar. */

  /* Calender type is 360 or 365.*/
  if ( td->sign == '-' )
    {
      /* Negative TimeDelta. */

      jd_return->sign = '-';
      /* Year. */
      jd_return->day = (-1) * td->year * ndiny;
      /* Month. */
      jd_return->day = jd_return->day - (ndiny - msdinm[base_dt->date.month - 1][NO_OF_MONTHS_IN_A_YEAR - td->month]);
      /* Day. */
      jd_return->day = jd_return->day + (-1) * td->day;
      /* Rest. */
      jd_return->ms = (-1) * td->hour * NO_OF_MS_IN_A_HOUR + (-1) * td->minute * NO_OF_MS_IN_A_MINUTE + (-1) * td->second * NO_OF_MS_IN_A_SECOND + (-1) * td->ms;
    }
  else if ( td->sign == '+' )
    {
      /* Positive TimeDelta. */

      jd_return->sign = '+';
      /* Year.  */
      jd_return->day = td->year * ndiny;
      /* Month. No of days in a TimeDelta depends on Base date.*/
      jd_return->day = jd_return->day + (msdinm[base_dt->date.month - 1][ td->month ]);
      /* Day. */
      jd_return->day = jd_return->day + td->day;
      /* Rest. */
      jd_return->ms = td->hour * NO_OF_MS_IN_A_HOUR + td->minute * NO_OF_MS_IN_A_MINUTE + td->second * NO_OF_MS_IN_A_SECOND + td->ms;
    }
  else
  return NULL; /* ERROR: TD sign not defined. */

  return jd_return;
}
else
return NULL;
}

  /* Internal function. */
  /* Converts a lib defined Julian delta value to a TimeDelta value. TimeDelta depends on the calendar type. 
   Notice that TimeDelta is not uniquely defined but depends on the definition of corresponding DateTime.

   The library assumes the following definition: Let A denote an anchor date and P a timedelta. For a 
   positive P, A + P = B where date B > A. Consequently, for a negative P, A + P = B where A > B. 
   Also, when P is positive, a delta of 1 month has as many days as in the month of anchor DateTime;
   a delta of 2 months corresponds to the number of days in the anchor date month and the next month and 
   so on. When P is negative, a delta of 1 month corresponds to as many days as in the month before the
   anchor date month; a delta of 2 month corresponds to as many days as in the month before the
   anchor date month and the month before that and so on. 
   
   For eg, TimeDelta of P01M, when the 'Anchor' DateTime is 2001-02-01T00:00:00.000 is equivalent to a
   Julian Delta of positive 28 days ( #days in February ) and 0 ms. Also, TimeDelta of P02M, when the 
   'Anchor' DateTime is 2001-02-01T00:00:00.000 is equivalent to a Julian Delta of positive 28+31 days
   ( #days in February PLUS #days in March) and 0 ms. Similarly, a TimeDelta of -P01M with 'Anchor' 
   DateTime of 2001-02-01T00:00:00.000 is equivalent to a Julian Delta of negative 31 days 
   ( #days in January ) and 0 ms. Likewise, TimeDelta of -P02M with 'Anchor' DateTime of 
   2001-02-01T00:00:00.000 is equivalent to a Julian Delta of negative 31+31 days ( #days in January
   PLUS #days in December) and 0 ms.
   */

struct _timedelta*
julianDeltaToTimeDelta(struct _juliandelta* jd, struct _datetime* base_dt, struct _timedelta* td_return)
{
  if ((jd != NULL) && (base_dt != NULL) && (td_return != NULL)){

  int i = 0;
  /* No of days in a year. */
  int ndiny;

  /* Pointer to an array of no of days in month. */
  const int* month_days;
  /* Pointer to an array of month_specific_delta_in_months. */
  const int (*msdinm)[NO_OF_MONTHS_IN_A_YEAR+1];

  /* Set parameter according to calender. */
  switch (getCalendarType())
    {
      case YEAR_OF_365_DAYS:

      month_days = month_days_365;
      msdinm = month_specific_delta_in_months_365;
      ndiny = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
      break;

      case YEAR_OF_360_DAYS:

      month_days = month_days_360;
      msdinm = month_specific_delta_in_months_360;
      ndiny = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE360;
      break;

      case PROLEPTIC_GREGORIAN:
      /* Handle Gregorian here. */
      /* Gregorian will have 366 days and 365 days depending on Leap year. */

      if ( jd->sign == '-' )
        {
          td_return->sign = '-';

          /* No of days in the final year */
          int delta_final_year;
          int64_t days = (-1)*jd->day;
          /* Set counter to base year and then jump forward to get to the final year.
           For each loop forward, increment year by 1.
           */
          int j = base_dt->date.year;
          /* Reset. */
          td_return->year = 0;
          while (days >= 0)
            {

              /* Loop over and get to the final year by substracting 366/365 days depending
               on leap/non-leap year. For each substraction, increment year by 1.  
               */

              /* The crucial point is month of february. */
              delta_final_year = days;
              if (
                  ( (testYearIsLeapYear(j + 1)) && (base_dt->date.month >= 3) )
                  ||
                  ( (testYearIsLeapYear(j)) && (base_dt->date.month < 3) )
              )
                {
                  /* If next year is leap year and base month is >= 3 
                   OR 
                   this year is a leap year and month is < 3
                   => delta of 1 year corresponds to 366 day julian delta.  
                   */
                  days = days - NO_OF_DAYS_IN_A_LEAP_YEAR;
                }
              else
                {
                  /* Otherwise. */
                  days = days - NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
                }
              td_return->year++;
              j++;
            }
          /* The loop ran one time too much. */
          if (days < 0)
            {
              td_return->year--;
              j--;
            }

          /* In final year. The crucial point is the month of february. */
          if (
              ( (testYearIsLeapYear(j + 1)) && (base_dt->date.month >= 3) )
              ||
              ( (testYearIsLeapYear(j)) && (base_dt->date.month < 3) )
          )
            {
              /* If final year's next year is a leap year and base month is >= 3 
               OR 
               final year is a leap year and month is < 3
               => An addition of leap-year specific delta for each month.
               */
              month_days = month_days_leapyear;
              msdinm = month_specific_delta_in_months_leapyear;
              ndiny = NO_OF_DAYS_IN_A_LEAP_YEAR;
            }
          else
            {
              /* Otherwise. */
              month_days = month_days_365;
              msdinm = month_specific_delta_in_months_365;
              ndiny = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
            }

          for (i = 1; i <= NO_OF_MONTHS_IN_A_YEAR; i++)
            {
              if (delta_final_year < msdinm[base_dt->date.month - 1][i])
                {
                  /* Month  */
                  td_return->month = i - 1;
                  /* Day. */
                  td_return->day = delta_final_year - msdinm[base_dt->date.month - 1][i - 1];
                  break;
                }
            }

          /* Time. */
          td_return->hour = ( (-1)* jd->ms ) / NO_OF_MS_IN_A_HOUR;
          td_return->minute = ( (-1)*jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR) / NO_OF_MS_IN_A_MINUTE;
          td_return->second = ( (-1)*jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR - td_return->minute * NO_OF_MS_IN_A_MINUTE) / NO_OF_MS_IN_A_SECOND;
          td_return->ms = (-1)*jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR - td_return->minute * NO_OF_MS_IN_A_MINUTE - td_return->second * NO_OF_MS_IN_A_SECOND;
        }
      else if (jd->sign == '+')
        {
          td_return->sign = '+';

          /* No days in the final year. */
          int delta_final_year;
          int64_t days = jd->day;
          /* Set counter to base year and then loop back to get to the final year.
           For each loop back, increment year by 1.
           */
          int j = base_dt->date.year;
          /* Reset */
          td_return->year = 0;
          while (days >= 0)
            {
              /* Loop over and get the year by substracting 366/365 days depending
               on leap/non-leap year. For each substraction, increment year by 1.  
               */

              /* The crucial point is month of february. */
              delta_final_year = days;
              if (
                  ( (testYearIsLeapYear(j - 1)) && (base_dt->date.month < 3) )
                  ||
                  ( (testYearIsLeapYear(j)) && (base_dt->date.month >= 3) )
              )
                {
                  /* If previous year is leap year and base month is < 3 
                   OR 
                   this year is a leap year and month is >= 3
                   => delta of 1 year corresponds to 366 day julian delta.  
                   */
                  days = days - NO_OF_DAYS_IN_A_LEAP_YEAR;
                }
              else
                {
                  /* Otherwise. */
                  days = days - NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
                }

              td_return->year++;
              j--;
            }
          /* The loop ran one time too much. */
          if (days < 0)
            {
              td_return->year--;
              j++;
            }

          /* In final year. The crucial point is the month of february. */
          if (
              ( (testYearIsLeapYear(j - 1)) && (base_dt->date.month < 3) )
              ||
              ( (testYearIsLeapYear(j)) && (base_dt->date.month >= 3) )
          )
            {
              /* If final year is a leap year and base month is >= 3 
               OR 
               final year's previous year is a leap year and month is < 3
               => An addition of leap-year specific delta for each month.
               */
              month_days = month_days_leapyear;
              msdinm = month_specific_delta_in_months_leapyear;
              ndiny = NO_OF_DAYS_IN_A_LEAP_YEAR;
            }
          else
            {
              /* Otherwise. */
              month_days = month_days_365;
              msdinm = month_specific_delta_in_months_365;
              ndiny = NO_OF_DAYS_IN_A_YEAR_FOR_CAL_TYPE365;
            }

          for (i = NO_OF_MONTHS_IN_A_YEAR; i > 0; i--)
            {
              if (delta_final_year < (ndiny - msdinm[base_dt->date.month - 1][i - 1]))
                {
                  /* Month */
                  td_return->month = NO_OF_MONTHS_IN_A_YEAR - i;
                  /* Day.  */
                  td_return->day = delta_final_year - (ndiny - msdinm[base_dt->date.month - 1][i]);
                  break;
                }
            }

          /* Time */
          td_return->hour = jd->ms / NO_OF_MS_IN_A_HOUR;
          td_return->minute = (jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR) / NO_OF_MS_IN_A_MINUTE;
          td_return->second = (jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR - td_return->minute * NO_OF_MS_IN_A_MINUTE) / NO_OF_MS_IN_A_SECOND;
          td_return->ms = jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR - td_return->minute * NO_OF_MS_IN_A_MINUTE - td_return->second * NO_OF_MS_IN_A_SECOND;
        }

      return td_return;

      default:
      return NULL; /* Calender type not defined. */
    }

  /* Handle 360 and 365 calendar type here as they are similar. */

  /* 360 and 365 day calendars */
  if (jd->sign == '+')
    {
      /* Positive delta.*/
      td_return->sign = '+';

      /* Year. */
      td_return->year = jd->day / ndiny;

      int delta_final_year = jd->day % ndiny;

      for (i = NO_OF_MONTHS_IN_A_YEAR; i > 0; i--)
        {
          if (delta_final_year < (ndiny - msdinm[base_dt->date.month - 1][i - 1]))
            {
              /* Month. */
              td_return->month = NO_OF_MONTHS_IN_A_YEAR - i;
              /* Day. */
              td_return->day = delta_final_year - (ndiny - msdinm[base_dt->date.month - 1][i]);
              break;
            }
        }

      /* Time. */
      td_return->hour = jd->ms / NO_OF_MS_IN_A_HOUR;
      td_return->minute = (jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR) / NO_OF_MS_IN_A_MINUTE;
      td_return->second = (jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR - td_return->minute * NO_OF_MS_IN_A_MINUTE) / NO_OF_MS_IN_A_SECOND;
      td_return->ms = jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR - td_return->minute * NO_OF_MS_IN_A_MINUTE - td_return->second * NO_OF_MS_IN_A_SECOND;
    }
  else if (jd->sign == '-')
    {
      /* Negative delta. */
      td_return->sign = '-';

      /* Year. */
      td_return->year = ((-1)*jd->day) / ndiny;

      int delta_final_year = ((-1)*jd->day) % ndiny;

      for (i = 1; i <= NO_OF_MONTHS_IN_A_YEAR; i++)
        {
          if (delta_final_year < msdinm[base_dt->date.month - 1][i])
            {
              /* Month. */
              td_return->month = i - 1;
              /* Day.  */
              td_return->day = delta_final_year - msdinm[base_dt->date.month - 1][i - 1];
              break;
            }
        }

      /* Time. */
      td_return->hour = (-1)*jd->ms / NO_OF_MS_IN_A_HOUR;
      td_return->minute = ((-1)*jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR) / NO_OF_MS_IN_A_MINUTE;
      td_return->second = ((-1)*jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR - td_return->minute * NO_OF_MS_IN_A_MINUTE) / NO_OF_MS_IN_A_SECOND;
      td_return->ms = (-1)*jd->ms - td_return->hour * NO_OF_MS_IN_A_HOUR - td_return->minute * NO_OF_MS_IN_A_MINUTE - td_return->second * NO_OF_MS_IN_A_SECOND;
    }
  else
  return NULL; /* ERROR: Sign of julian delta not defined. */

  return td_return;
}
else
return NULL;
}
/*! \endcond */

/**
 * @brief Get the TimeDelta between two Dates d1 and d2 as (d1-d2).
 *
 * Routine getTimeDeltaFromDate 'substracts' two Dates and returns the TimeDelta between
 * them. Internally, Dates are converted to DateTimes and then delta is calculated using 
 * getTimeDeltaFromDateTime().
 * 
 * This routine  handles all supported Calendar types; i.e. the translation from Calendar date 
 * to Julian date and conversion from Julian Delta to normal TimeDetla is Calendar-type dependent. 
 * For eg. for Calendar type Gregorian, the TimeDelta between 2001-02-01 and 2001-01-01 will be 1 month. 
 * Similarly, for Calendar of type 360-Day-Calendar, the TimeDelta will be 1 month. It must be noted 
 * however, that the two dates differ by 31 and 30 days respectively.
 *
 * @param  d1
 *         A pointer to struct _date. 
 *
 * @param  d2
 *         A pointer to struct _date.
 *
 * @param  td_return
 *         A pointer to struct _timedelta. Copy the result of (d1 - d2) in td_return.
 *
 * @return td_return
 *         A pointer to TimeDelta containing the result of substraction.
 */

struct _timedelta*
getTimeDeltaFromDate(struct _date* d1, struct _date* d2, struct _timedelta* td_return)
{
  if ((d1 != NULL )&& (d2 != NULL) && (td_return != NULL) ){
  /* Convert Date to datetime and resuse the DateTime interface to calculate time delta. */
  struct _datetime *dt1 = newDateTime("0-01-01T00:00:00.000");
  if (dt1 == NULL)
    return NULL;
  dt1 = convertDateToDateTime(d1, dt1);

  struct _datetime *dt2 = newDateTime("0-01-01T00:00:00.000");
  if ( dt2 == NULL )
    {
      deallocateDateTime(dt1);
      return NULL;
    }
  dt2 = convertDateToDateTime(d2, dt2);

  /* Call the Datetime function to get TD. dt1 - dt2 */
  td_return = getTimeDeltaFromDateTime(dt1,dt2,td_return);

  /* Cleanup. */
  deallocateDateTime(dt1);
  deallocateDateTime(dt2);

  return td_return;
}
else
return NULL;
}

/**
 * @brief Get total number of milliseconds in timedelta (Absolute value).
 *
 * Routine getTotalMilliSecondsTimeDelta returns the total number of milliseconds in TimeDelta. 
 * Notice that TimeDelta is not uniquely defined but depends on the definition of corresponding 
 * DateTime. TimeDelta is first converted to corresponding delta on the Julian axis. Julian delta 
 * is finally converted to the correct millisecond value.
 *
 * @param  td
 *         A pointer to struct _timedelta. Retrieve the number of milliseconds in this TD object.
 *
 * @param  base_dt
 *         A pointer to struct _datetime. Reference Datetime for the TD.
 *
 * @return totalmilliSeconds
 *         Integer value of totalmilliSeconds. 0 indicates error. TODO on Luis: Is this ok?
 */
int64_t
getTotalMilliSecondsTimeDelta(struct _timedelta* td, struct _datetime* base_dt)
{
  if ((td != NULL )&& (base_dt != NULL) ){
  int64_t totalmilliSeconds = 0;

  struct _juliandelta* jd = newJulianDelta('+', 0, 0);
  if ( jd == NULL )
    return 0;
  jd = timeDeltaToJulianDelta(td, base_dt, jd);
  if ( jd == NULL )
    {
      deallocateJulianDelta(jd);
      return 0;
    }
  totalmilliSeconds = jd->day * NO_OF_MS_IN_A_DAY + jd->ms;

  deallocateJulianDelta(jd);

  return totalmilliSeconds;
}
else
return 0;
}

  /**
   * @brief Get total number of seconds in timedelta (Absolute value).
   *
   * Routine getTotalSecondsTimeDelta returns the total number of seconds in TimeDelta. Notice that TimeDelta 
   * is not uniquely defined but depends on the definition of corresponding DateTime. Internally, number of seconds
   * is calculated by calling the routine function_totalmilliSeconds and then converting the millisecond value 
   * to seconds by dividing it by 1000.
   *
   * @param  td
   *         A pointer to struct _timedelta. Retrieve the number of seconds in this TD object.
   *
   * @param  base_dt
   *         A pointer to struct _datetime. Reference Datetime for the TD.
   *
   * @return totalSeconds
   *         Integer value of totalSeconds. 0 indicates error.
   */

int64_t
getTotalSecondsTimeDelta(struct _timedelta* td, struct _datetime* base_dt)
{
  if ((td != NULL )&& (base_dt != NULL) ){
  return getTotalMilliSecondsTimeDelta(td, base_dt)/NO_OF_MS_IN_A_SECOND;
}
else
return 0;
}

/**
* @brief Get TimeDelta as a string.
*
* timedeltaToString returns a string in IS08601 compliant (and extended) format.
*
* @param  td
*         A pointer to struct _timedelta. The timedelta to be converted to string.
*
* @param  toStr
*         A pointer to char. String where timedelta is to be written.
*
* @return toStr
*         A pointer to the string containing timedelta.
*/
char *
timedeltaToString(struct _timedelta* td, char* toStr)
{
if ((td != NULL )&& (toStr != NULL) ){
/* Return only non-negative values. */

memset(toStr,'\0',MAX_TIMEDELTA_STR_LEN);

if (td->sign == '-')
sprintf(toStr,"%cP",td->sign);
else if (td->sign == '+')
strcpy (toStr,"P");
else
return NULL; /*ERROR: TD sign not set. */

if (td->year != 0)
{
sprintf(&(toStr[strlen(toStr)]),"%" PRIi64 "Y",td->year);
}
if (td->month != 0)
{
sprintf(&(toStr[strlen(toStr)]),"%02dM",td->month);
}
if (td->day != 0)
{
sprintf(&(toStr[strlen(toStr)]),"%02dD",td->day);
}

sprintf(&(toStr[strlen(toStr)]),"T");

if (td->hour != 0)
{
sprintf(&(toStr[strlen(toStr)]),"%02dH",td->hour);
}
if (td->minute != 0)
{
sprintf(&(toStr[strlen(toStr)]),"%02dM",td->minute);
}
if ((td->second != 0) || (td->ms != 0))
{
sprintf(&(toStr[strlen(toStr)]),"%02d.%03dS",td->second,td->ms);
}

//Discard T if all time values are 0.
if(toStr[strlen(toStr)-1] == 'T')
{
toStr[strlen(toStr)-1] = '\0';
}
//Return P00.000S if all delta values are 0.
if(toStr[strlen(toStr)-1] == 'P')
{
strcat(toStr,"00.000S");
}

return toStr;
}
else
return NULL;
}

/**
* @brief Add timedelta to Date.
*
* Routine addTimeDeltaToDate adds a timedelta to a Date and returns the new Date. Both Date 
* and TimeDetla are first converted to corresponding values on the Julian axis. Addition is performed on
* the julian axis and the resulting Julian Date is converted back to the corrsponding Date. 
*
* The library assumes the following definition: Let A denote an anchor date and P a timedelta. For a 
* positive P, A + P = B where date B > A. Consequently, for a negative P, A + P = B where A > B. 
* Also, when P is positive, a delta of 1 month has as many days as in the month of anchor DateTime;
* a delta of 2 months corresponds to the number of days in the anchor date month and the next month and 
* so on. When P is negative, a delta of 1 month corresponds to as many days as in the month before the
* anchor date month; a delta of 2 month corresponds to as many days as in the month before the
* anchor date month and the month before that and so on.
*
* @param  d
*         A pointer to struct _date. The base date.
*
* @param  td
*         A pointer to struct _timedelta. The time delta to be added to d.
*
* @param  d_return
*         A pointer to struct _date. The result of addition is copied here.
*
* @return d_return
*         A pointer to the struct _date contianing the result of addition.
*/

struct _date *
addTimeDeltaToDate(struct _date* d, struct _timedelta* td, struct _date* d_return)
{
if ((d != NULL )&& (td != NULL) && (d_return != NULL) ){
/* Convert Date to Datetime and reuse the DateTime interface for Calculating the sum.*/
struct _datetime *dt = newDateTime("0-01-01T00:00:00.000");
if ( dt == NULL )
return NULL;

dt = convertDateToDateTime(d, dt);

struct _datetime *dt_return = newDateTime("0-01-01T00:00:00.000");
if ( dt_return == NULL )
return NULL;

/* Call the DateTime interface to calculate the new Datetime. */
dt_return = addTimeDeltaToDateTime(dt,td, dt_return);

/* Get Date from Datetime. */
d_return = convertDateTimeToDate(dt_return, d_return);

deallocateDateTime(dt);
deallocateDateTime(dt_return);

return d_return;
}
else
return NULL;
}

/**
* @brief Add timedelta to DateTime.
*
* Routine addTimeDeltaToDateTime adds a timedelta to a DateTime and returns the new DateTime. Both DateTime 
* and TimeDetla are first converted to corresponding values on the Julian axis. Addition is performed on
* the julian axis and the resulting Julian Date is converted back to the corrsponding DateTime.
*
* The library assumes the following definition: Let A denote an anchor date and P a timedelta. For a 
* positive P, A + P = B where date B > A. Consequently, for a negative P, A + P = B where A > B. 
* Also, when P is positive, a delta of 1 month has as many days as in the month of anchor DateTime;
* a delta of 2 months corresponds to the number of days in the anchor date month and the next month and 
* so on. When P is negative, a delta of 1 month corresponds to as many days as in the month before the
* anchor date month; a delta of 2 month corresponds to as many days as in the month before the
* anchor date month and the month before that and so on.
*
* @param  dt
*         A pointer to struct _datetime. The base datetime.
*
* @param  td
*         A pointer to struct _timedelta. The time delta to be added to dt.
*
* @param  dt_return
*         A pointer to struct _datetime. The result of addition is copied here.
*
* @return dt_return
*         A pointer to the struct _datetime contianing the result of addition.
*/

struct _datetime*
addTimeDeltaToDateTime(struct _datetime* dt, struct _timedelta* td, struct _datetime* dt_return)
{
  if ((dt != NULL )&& (td != NULL) && (dt_return != NULL) ){
    /* Convert base datetime to Julian. */
    struct _julianday* jd1 = newJulianDay(0, 0);
    if ( jd1 == NULL)
      return NULL;
    jd1 = date2julian(dt, jd1);

    /* Get julian delta. */
    struct _juliandelta* jd2 = newJulianDelta('+', 0, 0);
    if ( jd2 == NULL )
      {
        deallocateJulianDay(jd1);
        return NULL;
      }

    jd2 = timeDeltaToJulianDelta(td, dt, jd2);

    struct _julianday* jd = newJulianDay(0, 0);
    if ( jd == NULL )
      {
        deallocateJulianDay(jd1);
        deallocateJulianDelta(jd2);
        return NULL;
      }

    if ( td->sign == '+' )
      {
        jd = addJulianDelta(jd1, jd2, jd);
      }
    else if ( td->sign == '-' )
      {
        jd = substractJulianDelta(jd1, jd2, jd);
      }
    else
      return NULL; /* ERROR: Sign of timedelta is not defined. */

    /* Get the Datetime */
    dt_return = julian2date(jd, dt_return);

    deallocateJulianDay(jd1);
    deallocateJulianDelta(jd2);
    deallocateJulianDay(jd);

    return dt_return;
  }
  else
    return NULL;
}

/**
* @brief Get the timedelta between current_dt and start_dt plus next integral-multiple-of-timestep (timedelta).
*
* Routine moduloTimeDeltaFromDateTime returns the timedelta between the current DateTime (current_dt) and the event's next-trigger time.
* The next trigger time is defined as the the Anchor DateTime (start_dt) + N * TimeDelta(timestep) 
* where N is the minimum positive integer for which this sum is >= Current DateTime. In case 
* Anchor DateTime > Current DateTime, TimeDelta is calculated as start_dt - current_dt.
* 
* Notice that this TimeDelta will always be positive.
*
* @param  start_dt
*         A pointer to struct _datetime. The base datetime.
*
* @param  timestep
*         A pointer to struct _timedelta. delta between two consecutive triggers.
*
* @param  current_dt
*         A pointer to struct _datetime. The Current Date time.
*
* @param  modulo_td
*         A pointer to struct _timedelta. The timedelta between 'current datetime' and 'Start Datetime plus next integral-multiple-of-timestep' is copied here.
*
* @return modulo_td
*         A pointer to the struct _timedelta contianing the modulo timedelta. If Start time is in the future, returns on start_time - current_time.
*/

struct _timedelta*
moduloTimeDeltaFromDateTime(struct _datetime* start_dt, struct _timedelta* timestep, struct _datetime* current_dt, struct _timedelta* modulo_td)
{
  if ((start_dt != NULL )&& (timestep != NULL) && (current_dt != NULL) && (modulo_td != NULL) ){
  struct _datetime* dt_tmp = newDateTime("0-01-01T00:00:00.000");
  if ( dt_tmp == NULL )
  return NULL;

  if (compareDatetime(start_dt,current_dt)==(less_than))
    {
      /* Loop over */
      replaceDatetime(start_dt,dt_tmp);
      while(compareDatetime((dt_tmp = addTimeDeltaToDateTime(dt_tmp,timestep,dt_tmp)),current_dt) == less_than);

      /* Return n*dt_tmp - current_dt */
      modulo_td = getTimeDeltaFromDateTime(dt_tmp,current_dt,modulo_td);
    }
  else
    {
      /* Start time is in the future, return start_time - current_time. */
      modulo_td = getTimeDeltaFromDateTime(start_dt,current_dt,modulo_td);
    }

  deallocateDateTime(dt_tmp);

  return modulo_td;
  }
  else
    return NULL;
}


/**
* @brief Return the element-wise product of a scalar and a timedelta.
*
* elementwiseScalarMultiplyTimeDelta multiplies scalar lambda with each element of timedelta and returns the result in scaled_td.
* Scalar can be both positive and negative and so can the timedelta. The timedelta can not have days,months or years however: Only
* Timedeltas upto hours should call this routine. Also scaled_td->hour >= 24 will lead to an error.
*
*
* @param  base_td
*         A pointer to struct _timedelta. The base timedelta.
*
* @param  lambda
*         A scalar to be multiplied.
*
* @param  scaled_td
*         A pointer to struct _timedelta. The element-wise product of lambda and base_td is stored here.
*
* @return scaled_td
*	  A pointer to struct _timedelta. The filled structure containing the scaled timedelta values.	  
*/

struct _timedelta*
elementwiseScalarMultiplyTimeDelta(struct _timedelta* base_td, int64_t lambda, struct _timedelta* scaled_td)
{
  if ((base_td != NULL) && (scaled_td != NULL) )
    {
      /* Scalar multiplication not supported for TimeDeltas consisting of day/month/year. */
      if ( base_td->day > 0 || base_td->month > 0 || base_td->year > 0 )
        return NULL;

      /* Reset scaled_td to 0. */
      memset(scaled_td,0,sizeof(struct _timedelta));

      /* Scalar can be positive or negative. */
      if ( (lambda < 0) && (base_td->sign == '+') || (lambda > 0) && (base_td->sign == '-'))
        scaled_td->sign = '-';
      else
        scaled_td->sign = '+';

      /* Sign already handled above. Make lambda positive. */
      if (lambda < 0)
        lambda *= -1;


      /* Multiply each element by scalar. */

      scaled_td->ms += lambda*base_td->ms;
      while ( scaled_td->ms >= NO_OF_MS_IN_A_SECOND )
        {
	  scaled_td->ms -= NO_OF_MS_IN_A_SECOND;
	  scaled_td->second += 1;
	}
 
      scaled_td->second += lambda*base_td->second;
      while ( scaled_td->second >= 60 )
        {       
          scaled_td->second -= 60;
          scaled_td->minute += 1;
        }
   
      scaled_td->minute += lambda*base_td->minute;
      while ( scaled_td->minute >= 60 )
        {
          scaled_td->minute -= 60;
          scaled_td->hour += 1;
        }

     scaled_td->hour += lambda*base_td->hour;

     /* Scalar multiplication can not give a value in excess of 24 hours. */
     if ( scaled_td->hour >= NO_OF_HOURS_IN_A_DAY )
       {
	 /* ERROR: Return on NULL. */
         return NULL;
       }

      scaled_td->day 	= 0;
      scaled_td->month 	= 0;
      scaled_td->year 	= 0;

      return scaled_td;
	
    }
  else
    return NULL;
}



/**
* @brief Return the element-wise sum of two timedeltas.
*
* elementwiseAddTimeDeltatoTimeDelta adds two timedeltas elementwise and returns the result.
* Timedeltas beind added must be of the same sign; Substraction is not supported. 
* The timedelta can not have days,months or years however: Only Timedeltas upto hours should call this routine. 
* Also td_return->hour >= 24 will lead to an error.
*
*
* @param  td1
*         A pointer to struct _timedelta.
*
* @param  td2
*         A pointer to struct _timedelta.
*
* @param  td_return
*         A pointer to struct _timedelta. The element-wise sum of td1 and td2 is stored here.
*
* @return td_return
*         A pointer to struct _timedelta. The filled structure containing the added timedelta values.   
*/

struct _timedelta*
elementwiseAddTimeDeltatoTimeDelta(struct _timedelta* td1, struct _timedelta* td2,  struct _timedelta* td_return)
{
  if ( (td1 != NULL) && (td2 != NULL) && (td_return != NULL) )
    {

      /* TD addition not supported for TimeDeltas consisting of day/month/year. */
      if ( td1->day > 0 || td1->month > 0 || td1->year > 0 || td2->day > 0 || td2->month > 0 || td2->year > 0)
        return NULL;

      /*Reset td_return to 0.*/
      memset(td_return,0,sizeof(struct _timedelta));
      
      if(td1->sign == td2->sign)
        {
          /* If signs match, do add. */
          td_return->sign = td1->sign;

          td_return->ms += (td1->ms + td2->ms);  
          if ( td_return->ms >= NO_OF_MS_IN_A_SECOND )
            {
              td_return->ms -= NO_OF_MS_IN_A_SECOND;
              td_return->second += 1;
            }

          td_return->second += (td1->second + td2->second);
          if ( td_return->second >= 60 )
            {
              td_return->second -= 60;
              td_return->minute += 1;
            }

          td_return->minute += (td1->minute + td2->minute);
          if ( td_return->minute >= 60 )
            {
              td_return->minute -= 60;
              td_return->hour += 1;
            }

          td_return->hour += (td1->hour + td2->hour);
          if ( td_return->hour >= NO_OF_HOURS_IN_A_DAY )
            {
	      /* Hour-sum can not be allowed to exceed 24 hours. */
	      return NULL;
            }

	  td_return->day 	= 0;
	  td_return->month 	= 0;
          td_return->year 	= 0;
        } 
      else
        {
	  /* Substraction not supported. */
          return NULL; 
        }

      return td_return;
    } 
  else
    return NULL;

}

