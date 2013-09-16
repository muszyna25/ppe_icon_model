/**
 * @addtogroup CBindings libmtime C language bindings
 *
 * @file mtime_time.c
 * 
 * @brief Time and some operations supported on Time.
 *
 * @author Luis Kornblueh, Rahul Sinha. MPIM.
 * @date March 2013
 *
 */


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"mtime_time.h"

#include"mtime_iso8601.h"
#include"mtime_datetime.h"

/* Upto MilliSecond resolution supported. */


/**
 * @brief Construct new Time using an ISO 8601 conforming string.
 *
 * @param  ts
 *         A pointer to char. The string contains parameters with which Time is created.
 *
 * @return t
 *         A pointer to a filled Time. 
 */

struct _time*
newTime(const char* ts)
{
  if ((ts != NULL) && (getCalendarType()))
    {
      /* Convert ts to dts by appending dummy Date 0-01-01 for testing with verify_string_datetime(). */
      char* dts = (char*)calloc(MAX_DATETIME_STR_LEN,sizeof(char));
      if (dts == NULL )
        return NULL;

      strcpy(dts, "0-01-01");
      strncat(dts, ts, MAX_DATETIME_STR_LEN - 8);

      struct iso8601_datetime* isoDt = new_iso8601_datetime('+', 0, 0, 0, 0, 0, 0, 0, 'Z');
      if (isoDt == NULL )
        {
          free(dts);
          dts = NULL;
          return NULL ;
        }


      /* Verify ISO 8601 compliance. */
      if (verify_string_datetime(dts, isoDt) != DATETIME_MATCH)
        {
          free(dts);
          dts = NULL;
          deallocate_iso8601_datetime(isoDt);
	  isoDt = NULL;
          return NULL ;
        }

      struct _time* t = (struct _time*)calloc(1,sizeof(struct _time));
      if (t == NULL )
        {
          free(dts);
          dts = NULL;
          deallocate_iso8601_datetime(isoDt);
	  isoDt = NULL;
          return NULL ;
        }

      t->hour = isoDt->hour;
      t->minute = isoDt->minute;
      t->second = isoDt->second;
      t->ms = isoDt->ms;

      //Cleanup.
      free(dts);
      dts = NULL;
      deallocate_iso8601_datetime(isoDt);

      return t;
    }
  else
    return NULL ;
}


/**
 * @brief Construct new Time using 'raw' numerical values.
 *
 * @param  _hour
 *	  An "int" value denoting the hour part of time. 
 * @param  _minute
 *	  An "int" value denoting the minute part of time.
 * @param  _second
 *	  An "int" value denoting the second part of time.
 * @param  _ms
 *	  An "int" value denoting the milli-second part of time.
 *
 * @return t
 *         A pointer to a filled Time. 
 */

struct _time*
newRawTime(int _hour, int _minute, int _second, int _ms)
{
  char* ts = (char*)calloc(MAX_TIME_STR_LEN,sizeof(char));
  if (ts == NULL )
    return NULL ;

  snprintf(ts, MAX_TIME_STR_LEN, "%02d:%02d:%02d.%03d", _hour, _minute, _second, _ms);

  /* Reuse the string interface. */
  struct _time* t = newTime(ts);

  free(ts);
  ts = NULL;

  return t;
}


/**
 * @brief Copy the values and construct a new Time.
 *
 * @param  t
 *         A pointer to struct _time. Values of t are used to initialize the new time being created. 
 *
 * @return _t
 *         A pointer to an initialized Time object. 
 */

struct _time*
constructAndCopyTime(struct _time* t)
{
  if ( t != NULL )
    return newRawTime( t->hour, t->minute, t->second, t->ms );
  else
    return NULL;
}


/**
 * @brief Destructor of Time.
 *
 * @param  t
 *         A pointer to struct _time. t is deallocated.
 */

void
deallocateTime(struct _time* t)
{
  if ( t != NULL )
    {
      free(t);
      t = NULL;
    }
}


/**
 * @brief COPY a time object.
 *
 * Routine replaceTime copies the contents of source Time into a Destination Time object.
 *
 * @param  tsrc
 *         A pointer to struct _time. Copy "FROM" time object.
 *
 * @param  tdest
 *	  A pointer to struct _time. Copy "TO" time object.
 *
 * @return tdest
 *         A pointer to 'copied' time Object.
 */

struct _time*
replaceTime(struct _time* tsrc, struct _time* tdest)
{
  if ((tdest != NULL )&& (tsrc != NULL) ){

  tdest->hour = tsrc->hour;
  tdest->minute = tsrc->minute;
  tdest->second = tsrc->second;
  tdest->ms = tsrc->ms;

  return tdest;
  }
  else
    return NULL;
}


/**
 * @brief Get time as a string.
 *
 * timetoString returns a string in IS08601 compliant (and extended) format.
 *
 * @param  t
 *         A pointer to struct _time. The time to be converted to string.
 *
 * @param  toStr
 *         A pointer to char. String where time is to be written.
 *
 * @return toStr
 *         A pointer to the string containing time.
 */

char*
timeToString(struct _time* t, char* toStr)
{
  if ((t != NULL )&& ( toStr != NULL ) )
  {
    memset(toStr,'\0',MAX_TIME_STR_LEN);

    snprintf(toStr,MAX_TIME_STR_LEN,"%02d:%02d:%02d.%03dZ", t->hour, t->minute, t->second, t->ms);

    return toStr;
  }
  else 
    return NULL;
}


/**
 * @brief Get Time in 'struct tm' format and return as a string.
 *
 * @param  t
 *         A pointer to struct _time. The time to be converted to string.
 *
 * @param  toStr
 *         A pointer to char. String where time is to be written.
 *
 * @return toStr
 *         A pointer to the string containing time.
 */

char *
timeToPosixString(struct _time* t, char* toStr)
{
  if ((t != NULL )&& ( toStr != NULL ) )
  {
    struct tm tm_info;

    tm_info.tm_sec        = t->second;
    tm_info.tm_min        = t->minute;
    tm_info.tm_hour       = t->hour;
    
    memset(toStr,'\0',MAX_TIME_STR_LEN);
    strftime(toStr, MAX_TIME_STR_LEN , "%H:%M:%S", &tm_info);

    return toStr;
  }
  else
    return NULL;
}
