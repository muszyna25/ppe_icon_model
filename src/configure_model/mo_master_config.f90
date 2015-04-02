module mo_master_config

  use mtime, only: datetime, timedelta, newDatetime, newTimedelta

  implicit none
  
  private
  
  type(datetime), pointer, protected :: tc_exp_refdate => null()
                                                        
  type(datetime), pointer, protected :: tc_exp_startdate => null()
  type(datetime), pointer, protected :: tc_exp_stopdate => null()
                                                        
  type(datetime), pointer, protected :: tc_startdate => null()
  type(datetime), pointer, protected :: tc_stopdate => null()

  type(timedelta), pointer, protected :: tc_dt_checkpoint => null()
  type(timedelta), pointer, protected :: tc_dt_restart => null()

  logical :: lrestart
  
contains

  subroutine setRestart(lr)
    logical, intent(in) :: lr
    lrestart = lr
  end subroutine setRestart

  logical function isRestart()
    isRestart = lrestart
  end function isRestart

  subroutine setExpRefdate(experimentReferenceDate)   
    character(len=*), intent(in) :: experimentReferenceDate   
    tc_exp_refdate => newDatetime(experimentReferenceDate)
  end subroutine setExpRefdate

  subroutine setExpStartdate(experimentStartDate)   
    character(len=*), intent(in) :: experimentStartDate   
    tc_exp_startdate => newDatetime(experimentStartDate)   
  end subroutine setExpStartdate

  subroutine setExpStopdate(experimentStopDate)
    character(len=*), intent(in) :: experimentStopDate
    tc_exp_stopdate => newDatetime(experimentStopDate)
  end subroutine setExpStopdate

  subroutine setStartdate(startdate)
    character(len=*), intent(in) :: startdate
    tc_startdate => newDateime(startdate)
  end subroutine setStartdate

  subroutine setStopdate(stopdate)
    character(len=*), intent(in) :: stopdate
    tc_stopdate => newDateime(stopdate)
  end subroutine setStopdate

  subroutine setCheckpointTimeInterval(checkpointTimeIntval)
    character(len=*), intent(in) :: checkpointTimeIntval
    tc_dt_checkpoint => newTimedelta(checkpointTimeIntval)
  end subroutine setCheckpointTimeInterval
  
  subroutine setRestartTimeInterval(restartTimeIntval)
    character(len=*), intent(in) :: restartTimeIntval   
    tc_dt_restart => newTimedelta(restartTimeIntval)
  end subroutine setRestartTimeInterval

end module mo_master_config
