/**
 * @addtogroup CBindings libmtime C language bindings
 * @{
 *
 * @file mtime_eventHandling.h
 *
 * @brief  Definition of the basic event type and its methods. 
 *	   Definition of Even-Groups which contains a list of events.
 *
 * @author  Luis Kornblueh, Max Planck Institute for Meteorology.
 * @author  Rahul Sinha, Max Planck Institute for Meteorology.
 * 
 * @date March 2013
 */


#ifndef _MTIME_EVENTHANDLING_H
#define _MTIME_EVENTHANDLING_H

#include <stdint.h>
#include <stdbool.h>


///provides a Maximum string length for Group Names.
#define MAX_EVENT_STR_LEN 132 // TODO- On Luis. Need to agree on how a event_string should look.
///provides a Maximum string length for Group Names.
#define MAX_GROUPNAME_STR_LEN 132
///provides a Maximum string length for Event names.
#define MAX_EVENTNAME_STR_LEN 132

struct _datetime;
struct _timedelta;
struct _event;


/**
 * @struct _eventGroup
 *
 * @brief struct _eventGroup defines an Event-group. Each event group has an associated list of events.
 * 	  Event group is a place holder to 'group' events based on some user defined charteristics.
 * 	  	
 */
struct _eventGroup
{
  int64_t eventGroupId;		///< Event Group's ID.
  char* eventGroupName;		///< Event Group's name.

  struct _event* rootEvent;	///< Pointer to a List of events in this group.
};

struct _eventGroup*
newEventGroup(char* _eventGroupName);

void
deallocateEventGroup(struct _eventGroup* eg);

bool
addNewEventToEventGroup(struct _event* e, struct _eventGroup* eg);

bool
removeEventFromEventGroup(char* eventName, struct _eventGroup* eg);

int64_t
getEventGroupId(struct _eventGroup* eg);

char*
getEventGroupName(struct _eventGroup* eg, char* gname);

struct _event*
getEventGroupRootEvent(struct _eventGroup* eg);

/**
 * @struct _event
 *
 * @brief struct _event defines events. Events are set to be triggered at pre-specified intervals. 
 *              
 */
struct _event
{
  struct _event* nextEventInGroup;		///< Pointer to the next event in a given Event Group.

  int64_t eventId;				///< Auto generated. (For future use. Not in use yet).
  char* eventName;				///< Event's name.

  struct _datetime* eventReferenceDateTime;	///< Anchor datetime.

  struct _datetime* eventFirstDateTime;		///< Start datetime of the event. Can be empty string ? TODO: On Luis: we need to resolve this.
  struct _datetime* eventLastDateTime;		///< Last datetime of the event. Can be empty string  ? TODO: On Luis: we need to resolve this.

  struct _timedelta* eventInterval;		///< TimeDelta between succesive triggers.

  bool triggerCurrentEvent;			///< is this event active?


  //Properties of the event being triggered.

  bool nextEventIsFirst;			///< Is the next scheduled event the first to be triggered?

  bool eventisFirstInDay;	///< Is the event being triggered (i.e event when triggerCurrentEvent is true) First-In-Day?
  bool eventisFirstInMonth; 	///< Is the event being triggered (i.e event when triggerCurrentEvent is true) First-In-Month?
  bool eventisFirstInYear; 	///< Is the event being triggered (i.e event when triggerCurrentEvent is true) First-In-Year?
  bool eventisLastInDay; 	///< Is the event being triggered (i.e event when triggerCurrentEvent is true) Last-In-Day?
  bool eventisLastInMonth; 	///< Is the event being triggered (i.e event when triggerCurrentEvent is true) Last-In-Month?
  bool eventisLastInYear;	///< Is the event being triggered (i.e event when triggerCurrentEvent is true) Last-In-Year?

  struct _datetime* triggerNextEventDateTime;		///< Trigger the event next at this DT.
  struct _datetime* triggeredPreviousEventDateTime;	///< Last Trigger event happened at this DT.

};

struct _event*
newEvent(char* _eventName, char* _eventReferenceDate, char* _eventFirstDate, char* _eventLastDate, char* _eventInterval);

void
deallocateEvent(struct _event* e);

bool
isCurrentEventActive(struct _event* event, 
		     struct _datetime* current_dt, 
		     struct _timedelta* allowed_slack_add,
		     struct _timedelta* allowed_slack_subtract);

bool
iseventNextInNextDays(struct _event* e);

bool
iseventNextInNextMonths(struct _event* e);

bool
iseventNextInNextYears(struct _event* e);

char*
eventToString(struct _event* e, char* string);

struct _datetime*
getTriggerNextEventAtDateTime(struct _event* e, struct _datetime* dt_return);

struct _datetime*
getTriggeredPreviousEventAtDateTime(struct _event* e, struct _datetime* dt_return);

struct _event*
getNextEventInGroup(struct _event* e);

int64_t
getEventId(struct _event* e);

char*
getEventName(struct _event* e, char* ename);

struct _datetime*
getEventReferenceDateTime(struct _event* e);

struct _datetime*
getEventFirstDateTime(struct _event* e);

struct _datetime*
getEventLastDateTime(struct _event* e);

struct _timedelta*
getEventInterval(struct _event* e);

bool
getNextEventIsFirst(struct _event* e);

bool
getEventisFirstInDay(struct _event* e);

bool
getEventisFirstInMonth(struct _event* e);

bool
getEventisFirstInYear(struct _event* e);

bool
getEventisLastInDay(struct _event* e);

bool
getEventisLastInMonth(struct _event* e);

bool
getEventisLastInYear(struct _event* e);

/**
 * @}
 */

#endif
