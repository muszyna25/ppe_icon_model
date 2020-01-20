'''\
Object wrapper for libmtime Python bindings

This should eventually replace the original, very low-level bindings
'''

import ctypes
import os
import sys
import libmtime

CALENDAR_TYPE = libmtime.CALENDAR_TYPE

def setCalendar(calendar_type):
    '''(Re-)initialize calendar used for date and time calculations'''
    if libmtime.getCalendarType().cType != CALENDAR_TYPE.calendar_not_set:
        libmtime.resetCalendar()
    libmtime.setCalendar(calendar_type)

calendarToString = libmtime.calendarToString

class MTime (object):
    '''\
    Base class for all other data classes
    
    If the calendar has not been set before object creation,
    it defaults to Proleptic Gregorian
    '''
    def __init__(self):
        if '_lib' not in self.__dict__:
            self._lib = libmtime
            if libmtime.getCalendarType().cType == CALENDAR_TYPE.calendar_not_set:
                setCalendar(CALENDAR_TYPE.proleptic_gregorian)

class Date (MTime):
    def __init__(self, spec):
        super(Date, self).__init__()
        self._lib = libmtime
        self._my = libmtime.newDate(str(spec))
    def __del__(self):
        self._lib.deallocateDate(self._my)
    def __str__(self):
        return libmtime.dateToString(self._my)

class Time (MTime):
    def __init__(self, spec):
        super(Time, self).__init__()
        self._lib = libmtime
        self._my = libmtime.newTime(str(spec))
    def __del__(self):
        self._lib.deallocateTime(self._my)
    def __str__(self):
        return libmtime.timeToString(self._my)

class DateTime (MTime):
    def __init__(self, spec):
        super(DateTime, self).__init__()
        self._my = libmtime.newDateTime(str(spec))
    def __del__(self):
        self._lib.deallocateDateTime(self._my)
    def __str__(self):
        return libmtime.datetimeToString(self._my)
    def __cmp__(self, other):
        return libmtime.compareDatetime(self._my, other._my)
    def __add__(self, delta):
        result = DateTime(self)
        result._my = libmtime.addTimeDeltaToDateTime(
            result._my, delta._my, result._my)
        return result
    @property
    def date(self):
        return Date(libmtime.dateToString(ctypes.byref(self._my.contents.date)))
    @property
    def time(self):
        return Time(libmtime.timeToString(ctypes.byref(self._my.contents.time)))

class TimeDelta (MTime):
    def __init__(self, spec):
        super(TimeDelta, self).__init__()
        self._lib = libmtime
        self._my = libmtime.newTimeDelta(str(spec))
    def __del__(self):
        self._lib.deallocateTimeDelta(self._my)
    def __str__(self):
        return libmtime.timedeltaToString(self._my)
    def items(self):
        for (name, ctype) in self._my.contents._fields_:
            yield (name, self._my.contents.__getattribute__(name))

# Testing section

if __name__ == '__main__':
    import mtime

    # Test initialization and time delta arithmetics
    print('-'*80)
    initial_date = mtime.DateTime('1850-01-01')
    start_date = initial_date + mtime.TimeDelta('P1Y')
    print(initial_date)    
    print(start_date)

    # Test before/after tests on datetime objects
    print('-'*80)
    print(str(start_date) + ' == ' + str(start_date) + ': ' + str(start_date == start_date))
    print(str(start_date) + ' > ' + str(initial_date) + ': ' + str(start_date > initial_date))
    print(str(start_date) + ' == ' + str(initial_date) + ': ' + str(start_date == initial_date))
    print(str(start_date) + ' < ' + str(initial_date) + ': ' + str(start_date < initial_date))

    # Test re-initializing of calendar, showing Feb 28th and successing days

    # Use default calendar (Proleptic Gregorian)
    print('-'*80)
    start_day = mtime.DateTime('1852-02-28')
    one_day = mtime.TimeDelta('P1D')
    print(mtime.calendarToString())
    print(start_day)
    next_day = start_day + one_day
    print(next_day)
    next_day += one_day
    print(next_day)
    next_day += mtime.TimeDelta('P1D')
    print(next_day)

    # Use 365 day calendar (no Feb 29th)
    print('-'*80)
    mtime.setCalendar(mtime.CALENDAR_TYPE.year_of_365_days)
    print(mtime.calendarToString())
    print(start_day)
    next_day = start_day + mtime.TimeDelta('P1D')
    print(next_day)
    next_day += mtime.TimeDelta('P1D')
    print(next_day)
    next_day += mtime.TimeDelta('P1D')
    print(next_day)

    # Use 360 day calendar (twelve 30-day months)
    print('-'*80)
    mtime.setCalendar(mtime.CALENDAR_TYPE.year_of_360_days)
    print(mtime.calendarToString())
    print(start_day)
    next_day = start_day + mtime.TimeDelta('P1D')
    print(next_day)
    next_day += mtime.TimeDelta('P1D')
    print(next_day)
    next_day += mtime.TimeDelta('P1D')
    print(next_day)

    # Test selection and string conversion for date and time only parts
    print('-'*80)
    print(start_date)
    print(start_date.date)
    print(start_date.time)

