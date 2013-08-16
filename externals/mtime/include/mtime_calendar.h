/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_calendar.h
 * 
 * @brief Singleton Calendar connecting all supported calendar types.
 *
 * @details
 *
 * Three calendar types are provided:
 *
 *   - a proleptic Gregorian calendar
 *   - a calendar with 365 days per year without leap years
 *   - a calendar with 360 days per year and each month having 30 days
 *
 * To use this library, a call to initCalendar() with the respective selector (enum calendarType) must be
 * done first. The implementation is based on a singleton concept meaning that only one calendar
 * can be active at a time. To release a calendar a call to freeCalendar() has to be done.
 *
 * 
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology
 * @author  Rahul Sinha, Max Planck Institute for Meteorology
 *
 * @date March 2013
 *
 * @note Calendar type, once initialized, should not be changed.
 */


#ifndef _MTIME_CALENDAR_H
#define _MTIME_CALENDAR_H

///provides a string length for toString.
#define MAX_CALENDAR_STR_LEN 32

#define YEAR_UPPER_BOUND 2147483647
#define YEAR_LOWER_BOUND -2147483648

struct _datetime;
struct _julianday;

extern const int nDaysNonLeap[13];
extern const int nDaysLeap[13]; 
extern const int month_days_leapyear[];
extern const int month_days_365[];
extern const int month_days_360[];
extern const int month_specific_delta_in_months_leapyear[12][13];
extern const int month_specific_delta_in_months_365[12][13];
extern const int month_specific_delta_in_months_360[12][13];

typedef enum
{
  equal_to      =  0,
  greater_than  =  1,
  less_than     = -1,
  compare_error =  -128
} compare_return_val;

/**
 * @enum calendarType
 *
 * @brief enum calendarType lists the calendarTypes supported. The values are used for selecting calendars.
 *
 */
typedef enum
  {
    CALENDAR_NOT_SET = 0,  	///< calendar is not defined yet.
    PROLEPTIC_GREGORIAN = 1,  	///< proleptic Gregorian calendar.
    YEAR_OF_365_DAYS = 2,  	///< 365 day year without leap years.
    YEAR_OF_360_DAYS = 3  	///< 360 day year with 30 day months.
  }calendarType;


///@brief Function pointer connecting the Calendar to Julian routine. The pointed-to Function depends on the selected Calendar type.
extern struct _julianday *
(*date2julian)(struct _datetime *date, struct _julianday *julian);

///@brief Function pointer connecting Julian to Calendar routine. The pointed to Function depends on the selected Calendar type.
extern struct _datetime *
(*julian2date)(struct _julianday *julian, struct _datetime *date);


void
initCalendar(calendarType ct);
void
freeCalendar(void);

calendarType
getCalendarType(void);

char*
calendarToString(char *calendar);


/**
 * @}
 */
#endif
