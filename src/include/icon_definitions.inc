
!--------------------------------------------------
! timers definition
!needs:
!   USE mo_timer, ONLY: timer_start, timer_stop, timers_level, <timers_names>...
!
#if defined(__TIMERS__) || defined(__ALL_TIMERS__)
#define start_timer(timer,level) IF (timers_level >= level) CALL timer_start(timer)
#define stop_timer(timer,level) IF (timers_level >= level) CALL timer_stop(timer)
#else
#define start_timer(timer,level) ! no timer
#define stop_timer(timer,level)  ! no timer
#endif

#if defined(__ALL_TIMERS__)
#define start_detail_timer(timer,level) IF (timers_level >= level) CALL timer_start(timer)
#define stop_detail_timer(timer,level) IF (timers_level >= level) CALL timer_stop(timer)
#else
#define start_detail_timer(timer,level) ! no timer
#define stop_detail_timer(timer,level)  ! no timer
#endif


