MODULE mo_bcs_time_interpolation

  USE mo_kind, ONLY: wp
  USE mtime
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: t_time_interpolation_weights
  PUBLIC :: calculate_time_interpolation_weights
  PUBLIC :: store_time_interpolation_weight
  PUBLIC :: current_time_interpolation_weights
  
  TYPE t_time_interpolation_weights
    TYPE(datetime) :: reference_date
    REAL(wp)       :: weight1, weight2
    INTEGER        :: month1, month2
    INTEGER        :: year1, year2
    LOGICAL        :: initialized = .FALSE.
  END TYPE t_time_interpolation_weights

  TYPE(t_time_interpolation_weights), SAVE, PROTECTED :: current_time_interpolation_weights
  
CONTAINS

  SUBROUTINE store_time_interpolation_weight(time_interpolation_weight)
    TYPE(t_time_interpolation_weights), INTENT(in) :: time_interpolation_weight
    current_time_interpolation_weights = time_interpolation_weight
  END SUBROUTINE store_time_interpolation_weight

  FUNCTION calculate_time_interpolation_weights(current_date) RESULT(time_interpolation_weight)
    TYPE(t_time_interpolation_weights) :: time_interpolation_weight
    TYPE(datetime), POINTER, INTENT(in) :: current_date

    TYPE(datetime), POINTER :: next_month => NULL()
    TYPE(datetime), POINTER :: previous_month => NULL()
    TYPE(timedelta), POINTER :: one_month => NULL()

    INTEGER :: seconds_in_month
    INTEGER :: seconds_in_middle_of_previous_month, seconds_in_middle_of_month, seconds_in_middle_of_next_month
    INTEGER :: days_in_previous_month, days_in_month, days_in_next_month
    
    time_interpolation_weight%reference_date = current_date 
    
    days_in_month = getNoOfDaysInMonthDateTime(current_date) 
    seconds_in_middle_of_month = 43200 * days_in_month          ! = 86400 * my_month_len / 2
    
    seconds_in_month = getNoOfSecondsElapsedInMonthDateTime(current_date)

    IF (seconds_in_month <= seconds_in_middle_of_month) THEN

      ! first half of month 

      one_month => newTimedelta('-P1M')
      previous_month => newDatetime(current_date)
      previous_month = current_date + one_month
      days_in_previous_month = getNoOfDaysInMonthDateTime(previous_month)
      seconds_in_middle_of_previous_month = 43200 * days_in_previous_month          ! = 86400 * my_month_len / 2

      ! simple linear interpolation

      time_interpolation_weight%weight1 = REAL(seconds_in_middle_of_month - seconds_in_month,wp) &
           &                             /REAL(seconds_in_middle_of_month + seconds_in_middle_of_previous_month,wp)
      time_interpolation_weight%weight2 = 1.0_wp - time_interpolation_weight%weight1
      time_interpolation_weight%month1 = current_date%date%month - 1 ! does indexing only, so do not use previous_month%date%month
      time_interpolation_weight%month2 = current_date%date%month
      time_interpolation_weight%year1 = previous_month%date%year
      time_interpolation_weight%year2 = current_date%date%year

    ELSE
      
      ! second half of month

      one_month => newTimedelta('P1M')
      next_month => newDatetime(current_date)
      next_month = current_date + one_month      
      days_in_next_month = getNoOfDaysInMonthDateTime(next_month)
      seconds_in_middle_of_next_month = 43200 * days_in_next_month         ! = 86400 * my_month_len / 2

      ! simple linear interpolation

      time_interpolation_weight%weight2 = REAL(seconds_in_month - seconds_in_middle_of_month,wp) &
           &                             /REAL(seconds_in_middle_of_month + seconds_in_middle_of_next_month,wp)
      time_interpolation_weight%weight1 = 1.0_wp - time_interpolation_weight%weight2
      time_interpolation_weight%month1 = current_date%date%month
      time_interpolation_weight%month2 = current_date%date%month + 1 ! does indexing only, so do not use next_month%date%month
      time_interpolation_weight%year1 = current_date%date%year
      time_interpolation_weight%year2 = next_month%date%year

    ENDIF

    CALL deallocateTimedelta(one_month)
    IF (ASSOCIATED(previous_month)) CALL deallocateDatetime(previous_month)
    IF (ASSOCIATED(next_month)) CALL deallocateDatetime(next_month)

    time_interpolation_weight%initialized = .TRUE.
    
  END FUNCTION calculate_time_interpolation_weights
    
END MODULE mo_bcs_time_interpolation
