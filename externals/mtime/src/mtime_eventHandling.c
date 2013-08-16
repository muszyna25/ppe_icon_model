/**
 * @file mtime_eventHandling.c
 *
 * @brief Event-groups which contains a list of events.
 *
 * @author Luis Kornblueh, Rahul Sinha. MPIM.
 * @date March 2013
 *
 * @note
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "mtime_eventHandling.h"

#include "mtime_date.h"
#include "mtime_time.h"
#include "mtime_datetime.h"
#include "mtime_timedelta.h"
#include "mtime_eventList.h"
#include "mtime_iso8601.h"

//TODO: Need to handle cases where some parameters are left out in the init. Discuss each value independently.

// The IDs are unique only when they are live. Once a group/event has been deleted, the IDs will be reused.
// Currently, the IDs are generated but not used.
static int64_t EVENTGROUPID = 0;
static int64_t EVENTID = 0;


/**
 * @brief Construct new event-Group using a string.
 *
 * @param  egn
 *         A pointer to char. This string contains the name of event group.
 *
 * @return eg
 *         A pointer to a initialized event-Group. 
 *
 */

struct _eventGroup*
newEventGroup(char* egn)
{
if ( egn != NULL )
{
  struct _eventGroup* eg = (struct _eventGroup*)calloc(1,sizeof(struct _eventGroup));
  if (eg == NULL )
    return NULL ;

  /* Copy Group name. */
  eg->eventGroupName = (char*)calloc(MAX_GROUPNAME_STR_LEN,sizeof(char));
  if (eg->eventGroupName == NULL )
    {
      free(eg);
      eg = NULL;
      return NULL ;
    }
  strncpy(eg->eventGroupName, egn, MAX_GROUPNAME_STR_LEN-1);
  eg->eventGroupName[MAX_GROUPNAME_STR_LEN - 1] = '\0';

  /* Generate eventGroupId. */
  EVENTGROUPID = EVENTGROUPID + 1;
  eg->eventGroupId = EVENTGROUPID;

  /* Initialize a NULL pointer. Stores a List of events associated with 'this' Event Group. */
  eg->rootEvent = NULL;

  return eg;
}
else
  return NULL;
}


/**
 * @brief Destructor of EventGroup.
 *
 * @param  eg
 *         A pointer to struct _eventGroup. eg is deallocated.
 */

void
deallocateEventGroup(struct _eventGroup* eg)
{
  if (eg != NULL )
    {
      EVENTGROUPID = EVENTGROUPID - 1;

      free(eg->eventGroupName);
      eg->eventGroupName = NULL;

      /* Deallocate all events in the list. */
      deallocateEventsInGroup(eg);

      free(eg);
      eg = NULL;
    }
}


/**
 * @brief Add new event to an eventgroup.
 *
 * @param  e
 *         A pointer to struct _event. The event to be added.
 *
 * @param  eg
 *         A pointer to struct _eventGroup. The eventGroup where the event is added.
 *
 * @return bool
 *         true/false indicating success or failure of addition.
 */

bool
addNewEventToEventGroup(struct _event* e, struct _eventGroup* eg)
{
  if ((e != NULL )&& (eg != NULL) )
    {
      if ( addNewEventToGroup(e,eg) == true )
        return true;
      else
        return false;
    }
  else
    return false;
}


/**
 * @brief Remove event from eventgroup. Also, deallocate the event.
 *
 * @param  en
 *         A pointer to char. The name of event to be removed.
 *
 * @param  eg
 *         A pointer to struct _eventGroup. The eventGroup to which this event belongs.
 *
 * @return bool
 *         true/false indicating success or failure of removal.
 */

bool
removeEventFromEventGroup(char* en, struct _eventGroup* eg)
{
  if ((en != NULL )&& (eg != NULL) )
    return removeEventWithNameFromGroup(en,eg);
  else
    return false;
}

int64_t
getEventGroupId(struct _eventGroup* eg)
{
  if (eg != NULL)
    return eg->eventGroupId;
  else
    return 0;
}

char*
getEventGroupName(struct _eventGroup* eg, char* gname)
{
  if ((eg != NULL) && (gname != NULL))
    {
      strncpy(gname,eg->eventGroupName,MAX_GROUPNAME_STR_LEN);
      return gname;
    }
  else
    return NULL;
}

struct _event*
getEventGroupRootEvent(struct _eventGroup* eg)
{
  if (eg != NULL)
    return eg->rootEvent;
  else
    return NULL;
}


  /**
   * @brief Construct new event using strings.
   *
   * @param  _en
   *         A pointer to char. This string contains the name of event.
   * @param  _eventReferenceDT
   *         A pointer to char. This string contains the First Reference date.
   * @param  _eventFirstDT
   *         A pointer to char. This string contains the Starting time.
   * @param  _eventLastDT
   *         A pointer to char. This string contains the Ending time.
   * @param  _eventInterval
   *         A pointer to char. This string contains the timestep.
   * @return e
   *         A pointer to an initialized event. 
   *
   */

struct _event*
newEvent(char* _en, char* _eventReferenceDT, char* _eventFirstDT, char* _eventLastDT, char* _eventInterval)
{
  struct _event* e = (struct _event*)calloc(1,sizeof(struct _event));
  if (e == NULL )
    return NULL ;

  /* Initialize a null pointer. Connectes to the 'next' event in the list of Events in an Event Group.  */
  e->nextEventInGroup = NULL;

  /* Copy event name. */
  e->eventName = (char*)calloc(MAX_EVENTNAME_STR_LEN,sizeof(char));
  if (e->eventName == NULL )
    {
      free(e);
      e = NULL;
      return NULL ;
    }
  strncpy(e->eventName, _en, MAX_EVENTNAME_STR_LEN-1);
  e->eventName[MAX_EVENTNAME_STR_LEN - 1] = '\0';

  //Generate eventId.
  EVENTID = EVENTID + 1;
  e->eventId = EVENTID;

  if (_eventReferenceDT)
    {
      e->eventReferenceDateTime = newDateTime(_eventReferenceDT);
      if (e->eventReferenceDateTime == NULL )
        {
          free(e->eventName);
          e->eventName = NULL;
          free(e);
          e = NULL;
          return NULL ;
        }
    }
  else
    {
      //TODO:
    }

  if (_eventFirstDT)
    {
      e->eventFirstDateTime = newDateTime(_eventFirstDT);
      if (e->eventFirstDateTime == NULL )
        {
          free(e->eventName);
          e->eventName = NULL;
          deallocateDateTime(e->eventReferenceDateTime);
          free(e);
          e = NULL;
          return NULL ;
        }
    }
  else
    {
      //TODO:
    }

  if (_eventLastDT)
    {
      e->eventLastDateTime = newDateTime(_eventLastDT);
      if (e->eventLastDateTime == NULL )
        {
          free(e->eventName);
          e->eventName = NULL;
          deallocateDateTime(e->eventReferenceDateTime);
          deallocateDateTime(e->eventFirstDateTime);
          free(e);
          e = NULL;
          return NULL ;
        }
    }
  else
    {
      //TODO:
    }

  if (_eventInterval)
    {
      e->eventInterval = newTimeDelta(_eventInterval);
      if (e->eventInterval == NULL )
        {
          free(e->eventName);
          e->eventName = NULL;
          deallocateDateTime(e->eventReferenceDateTime);
          deallocateDateTime(e->eventFirstDateTime);
          deallocateDateTime(e->eventLastDateTime);
          free(e);
          e = NULL;
          return NULL ;
        }
    }
  else
    {
      //TODO:
    }

  /* Initialize with 'some' value. */
  e->triggerNextEventDateTime = newDateTime("0-01-01T00:00:00.000");
  e->triggeredPreviousEventDateTime = newDateTime("0-01-01T00:00:00.000");

  /* Intialize the next trigger to be the first trigger. */
  e->triggerNextEventDateTime = replaceDatetime(e->eventFirstDateTime, e->triggerNextEventDateTime);

  e->triggerCurrentEvent = false;
  e->nextEventIsFirst = true;

  e->eventisFirstInDay = false;
  e->eventisFirstInMonth = false;
  e->eventisFirstInYear = false;
  e->eventisLastInDay = false;
  e->eventisLastInMonth  = false;
  e->eventisLastInYear = false;

  return e;
}


/**
 * @brief Destructor of Event.
 *
 * @param  e
 *         A pointer to struct _event. e is deallocated.
 */

void
deallocateEvent(struct _event* e)
{
  if (e != NULL )
    {
      EVENTID = EVENTID - 1;

      free(e->eventName);
      e->eventName = NULL;

      deallocateDateTime(e->eventReferenceDateTime);
      deallocateDateTime(e->eventFirstDateTime);
      deallocateDateTime(e->eventLastDateTime);
      deallocateDateTime(e->triggerNextEventDateTime);
      deallocateDateTime(e->triggeredPreviousEventDateTime);

      deallocateTimeDelta(e->eventInterval);

      e->nextEventInGroup = NULL;

      free(e);
      e = NULL;
    }
}


/*! \cond PRIVATE */
/*  Internal event. Should not be exposed to the user.*/
static
void
setEvent(struct _event* e)
{
  if (e != NULL )
    e->triggerCurrentEvent = true;
}
/*! \endcond */

/**
 * @brief Check if this event is active by comparing event's trigger time with current_dt.
 *
 * @param  event
 *         A pointer to struct _event. This is the event being tested.
 *
 * @param  current_dt
 *         A pointer to struct _datetime. This is the 'current' datetime of the system.
 *
 * @return bool
 *         true/false indicating if the event is active.
 */

bool
isCurrentEventActive(struct _event* event, struct _datetime* current_dt)
{
  if ((event != NULL )&& (current_dt != NULL) ){
  if(compareDatetime(current_dt,event->triggerNextEventDateTime) >= equal_to)
    {
      /* If current Datetime is ahead of next trigger datetime, Event is active. */


      /* If the event being triggred is the first event to be triggered, it is all of FirstIn*  */
      /* If previous triggered event had a different day/month/year, this event must be FirstIn* */
      if( (event->nextEventIsFirst == true) || (iseventNextInNextDays(event) == true))
	{
	  event->eventisFirstInDay = true;
	}
      if( (event->nextEventIsFirst == true)  || (iseventNextInNextMonths(event) == true))
	{
	  event->eventisFirstInMonth = true;
	}
      if( (event->nextEventIsFirst == true)  || (iseventNextInNextYears(event) == true))
	{
	  event->eventisFirstInYear = true;
	}



      /* Set the new next-trigger datetime and the new previous-triggered-datetime. */
      replaceDatetime(event->triggerNextEventDateTime, event->triggeredPreviousEventDateTime);
      event->triggerNextEventDateTime = addTimeDeltaToDateTime(event->triggerNextEventDateTime,event->eventInterval,event->triggerNextEventDateTime);

      /* Set event. */
      setEvent(event);



      /* If the future event (not the current event being triggered) has a different day/month/year, current event must be LastIn* */
      if(iseventNextInNextDays(event) == true)
	{
	  event->eventisLastInDay = true;
	}
      
      if(iseventNextInNextMonths(event) == true)
	{
	  event->eventisLastInMonth = true;
	}
      if(iseventNextInNextYears(event) == true)
	{
	  event->eventisLastInYear = true;
	}

      return true;
    }
  else if ((event->triggerCurrentEvent == true) && (compareDatetime(current_dt,event->triggerNextEventDateTime) == less_than))
    {
      /* This is the next iteration after event was set in last call. Reset event. Return false. */
      event->triggerCurrentEvent = false;

      event->eventisFirstInDay = false;
      event->eventisFirstInMonth = false;
      event->eventisFirstInYear = false;
      event->eventisLastInDay = false;
      event->eventisLastInMonth	 = false;
      event->eventisLastInYear = false;
      
      event->nextEventIsFirst = false;
    }

  return false;
}
else
return false;
}


/* Is next event in next day(s) of previous event. */
bool
iseventNextInNextDays(struct _event* e)
{
  if( e != NULL )
    {
      bool ret = false;

      struct _date* d_next = newRawDate(0,1,1);
      struct _date* d_prev = newRawDate(0,1,1);
      d_next = convertDateTimeToDate(e->triggerNextEventDateTime,d_next);
      d_prev = convertDateTimeToDate(e->triggeredPreviousEventDateTime,d_prev);

      if(compareDate(d_next,d_prev) == greater_than)
	{
	  ret = true;
	}
      else
	{
	  ret = false;
	}

      deallocateDate(d_next);
      deallocateDate(d_prev);

      return ret;
    }
  else
    return false;
}


/* Is next event in next month(s) of previous event. */
bool
iseventNextInNextMonths(struct _event* e)
{
  if( e != NULL )
    {
      bool ret = false;

      struct _date* d_next = newRawDate(0,1,1);
      struct _date* d_prev = newRawDate(0,1,1);
      d_next = convertDateTimeToDate(e->triggerNextEventDateTime,d_next);
      d_prev = convertDateTimeToDate(e->triggeredPreviousEventDateTime,d_prev);

      d_next->day = 1;
      d_prev->day = 1;

      if(compareDate(d_next,d_prev) == greater_than)
	{
	  ret = true;
	}
      else
	{
	  ret = false;
	}

      deallocateDate(d_next);
      deallocateDate(d_prev);

      return ret;
    }
  else
    return false;
}

/* Is next event in next year(s) of previous event. */
bool
iseventNextInNextYears(struct _event* e)
{
  if( e != NULL )
    {
      bool ret = false;

      struct _date* d_next = newRawDate(0,1,1);
      struct _date* d_prev = newRawDate(0,1,1);
      d_next = convertDateTimeToDate(e->triggerNextEventDateTime,d_next);
      d_prev = convertDateTimeToDate(e->triggeredPreviousEventDateTime,d_prev);

      d_next->day = 1;
      d_next->month = 1;
      
      d_prev->day = 1;
      d_prev->month = 1;


      if(compareDate(d_next,d_prev) == greater_than)
	{
	  ret = true;
	}
      else
      {
	ret = false;
      }

      deallocateDate(d_next);
      deallocateDate(d_prev);

      return ret;
    }
  else
    return false;
}


  /**
   * @brief Get Event as a string.
   *
   * @param  e
   *         A pointer to struct _event. The event to be converted to string.
   *
   * @param  string
   *         A pointer to char. String where event is to be written.
   *
   * @return string
   *         A pointer to the string containing event description.
   */
//TODO on Luis. Exact format.
char*
eventToString(struct _event* e, char* string)
{
  if ((e != NULL )&& (string != NULL) )
    {
      memset(string,'\0',MAX_EVENT_STR_LEN);
      if (e->eventName != NULL)
	{
      	  sprintf(string, "%s", e->eventName);
          return string;
	}
      else
	{  
	   //TODO.
	   return NULL;
	}
    }
  else
    return NULL;
}

  /**
   * @brief Get the Datetime when 'this' event will be triggered next.
   *
   * @param  e
   *         A pointer to struct _event. This is the event being queried.
   *
   * @param  dt_return
   *         A pointer to struct _datetime. The next trigger datetime is copied here.
   *
   * @return dt_return
   *         A pointer to DateTime with next-trigger datetime.
   */

struct _datetime*
getTriggerNextEventAtDateTime(struct _event* e, struct _datetime* dt_return)
{
  if ((e != NULL )&& (dt_return != NULL) ){
  replaceDatetime(e->triggerNextEventDateTime, dt_return);
  return dt_return;
}
else
return NULL;
}

  /**
   * @brief Get the Datetime when 'this' event was triggered last.
   *
   * @param  e
   *         A pointer to struct _event. This is the event being queried.
   *
   * @param  dt_return
   *         A pointer to struct _datetime. The last trigger datetime is copied here.
   *
   * @return dt_return
   *         A pointer to DT with previous-trigger datetime.
   */
  /* If the event was never tiggered, 0-01-01T00:00:00.000 is returned. */
struct _datetime*
getTriggeredPreviousEventAtDateTime(struct _event* e, struct _datetime* dt_return)
{
  if ((e != NULL )&& (dt_return != NULL) ){
  replaceDatetime(e->triggeredPreviousEventDateTime, dt_return);
  return dt_return;
}
else
return NULL;
}


struct _event*
getNextEventInGroup(struct _event* e)
{
  if (e != NULL)
    return e->nextEventInGroup;   
  else
    return NULL;
}

int64_t
getEventId(struct _event* e)
{
  if (e != NULL)
    return e->eventId;
  else
    return NULL;
}

char*
getEventName(struct _event* e, char* ename)
{
  if ((e != NULL) && (ename != NULL))
    {
      strncpy(ename, e->eventName, MAX_EVENTNAME_STR_LEN);    
      return ename;
    }
  else
    return NULL;
}

struct _datetime*
getEventReferenceDateTime(struct _event* e)
{
  if (e != NULL)
    return e->eventReferenceDateTime;
  else
    return NULL;
}

struct _datetime*
getEventFirstDateTime(struct _event* e)
{
  if (e != NULL)
    return e->eventFirstDateTime;
  else
    return NULL;
}

struct _datetime*
getEventLastDateTime(struct _event* e)
{
  if (e != NULL)
    return e->eventLastDateTime;
  else
    return NULL;
}

struct _timedelta*
getEventInterval(struct _event* e)
{
  if (e != NULL)
    return e->eventInterval;
  else
    return NULL;
}

bool
getNextEventIsFirst(struct _event* e)
{
  if (e != NULL)
    return e->nextEventIsFirst;
  else
    return NULL;
}

bool
getEventisFirstInDay(struct _event* e)
{
  if (e != NULL)
    return e->eventisFirstInDay;
  else
    return NULL;
}

bool
getEventisFirstInMonth(struct _event* e)
{
  if (e != NULL)
    return e->eventisFirstInMonth;
  else
    return NULL;
}

bool
getEventisFirstInYear(struct _event* e)
{
  if (e != NULL)
    return e->eventisFirstInYear;
  else
    return NULL;
}
bool
getEventisLastInDay(struct _event* e)
{
  if (e != NULL)
    return e->eventisLastInDay;
  else
    return NULL;
}

bool
getEventisLastInMonth(struct _event* e)
{
  if (e != NULL)
    return e->eventisLastInMonth;
  else
    return NULL;
}

bool
getEventisLastInYear(struct _event* e)
{
  if (e != NULL)
    return e->eventisLastInYear;
  else
    return NULL;
}
