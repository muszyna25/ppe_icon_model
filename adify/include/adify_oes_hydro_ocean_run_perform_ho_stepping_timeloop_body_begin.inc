#if defined(__COMPAD_CPTEST__) || defined(__COMPAD_ADJ_IO_TEST__)
  CALL STORE_CHECKPOINTS(jstep+1,jg,ocean_state,sea_ice,p_oce_sfc, p_atm_f,&
       & current_time,time_config,output_file)
  WRITE(istr1,'(i010)') jstep+1
  WRITE(*,*) 'CPTEST_BEGIN OCE-TimeStep-',TRIM(ADJUSTL(istr1)),'-1'
#else         
#  ifdef __DCO_A1S__
     WRITE(istr1,'(i010)') jstep+1
     WRITE(*,*) 'ICON_AD ADJ_TIMELOOP_BEGIN OCE-TimeStep-',TRIM(ADJUSTL(istr1)),'-1'
#    ifndef __DCO_A1S_NOCP__
       CALL STORE_CHECKPOINTS(jstep+1,jg,ocean_state,sea_ice,p_oce_sfc, p_atm_f,&
         & current_time,time_config,output_file)
#    endif /* NOT __DCO_A1S_NOCP__ */
#  endif /* __DCO_A1S__ */
#endif /*  __COMPAD_CPTEST__ */
