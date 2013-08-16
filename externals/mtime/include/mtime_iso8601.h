/*! \cond PRIVATE */
/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_iso8601.h
 * 
 * @brief Routines and Data-Structures for checking ISO860 compliannce.
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology.
 * @author  Rahul Sinha, Max Planck Institute for Meteorology.
 *
 * @date March 2013
 *
 */


#ifndef _MTIME_ISO8601_H
#define _MTIME_ISO8601_H

#include <stdint.h>

typedef enum
{
  FAILURE = 0, DATETIME_MATCH = 1, DURATION_MATCH = 2
} ISO8601_STATUS;

struct iso8601_duration
{
  char sign;
  int64_t year;
  int month;
  int day;
  int hour;
  int minute;
  int second;
  int ms;
};

struct iso8601_datetime
{
  /* Year can be positive or negative. */
  char sign_of_year;
  int64_t year;
  int month;
  int day;
  int hour;
  int minute;
  int second;
  int ms;
  char timeZoneOffSet;
};

struct iso8601_datetime*
new_iso8601_datetime(char _sign_of_year, int64_t _year, int _month, int _day, int _hour, int _minute, int _second, int _ms,
    char _timeZoneOffSet);

void
deallocate_iso8601_datetime(struct iso8601_datetime* iso8601_datetimeObj);

struct iso8601_duration*
new_iso8601_duration(char _sign, int64_t _year, int _month, int _day, int _hour, int _minute, int _second, int _ms);

void
deallocate_iso8601_duration(struct iso8601_duration* iso8601_durationObj);

ISO8601_STATUS
verify_string_datetime(const char* test_string, struct iso8601_datetime* dummy_isoDtObj);

ISO8601_STATUS
verify_string_duration(const char* test_string, struct iso8601_duration* dummy_isoDObj);

/**
 * @}
 */
#endif
