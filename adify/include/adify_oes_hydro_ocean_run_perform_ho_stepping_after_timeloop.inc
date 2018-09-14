#ifdef __COMPAD_CPTEST__
      !! store number of performed forward steps
      nsteps = jstep

      PRINT *,'  NSTEPS :: ', nsteps
      !! store additional cp for verification
      call store_checkpoints(jstep0+nsteps+1,jg,ocean_state,sea_ice,p_oce_sfc,p_atm_f, &
          & current_time,time_config,output_file)

      !!  prepare time step counter
      ALLOCATE(tscnt(nsteps))
      tscnt = 1
      
      !!  prepare pseudo random start CP
      ALLOCATE(iperm(nsteps))
      iperm = (/ (i,i=1,nsteps) /)
      DO i = jstep0+1, nsteps
         CALL RANDOM_NUMBER( r )
         iswap(1) = 1 + INT(REAL(nsteps-1) * r ) 
         CALL RANDOM_NUMBER( r )
         iswap(2) = 1 + INT(REAL(nsteps-1) * r ) 
         iswap(3)         = iperm(iswap(1))
         iperm(iswap(1))  = iperm(iswap(2))
         iperm(iswap(2))  = iswap(3)
      END DO
  !!$    DO i = jstep0+1, nsteps
  !!$       iperm(i)  = i
  !!$    END DO
      WRITE(0,*)' CPTEST IPERM ', iperm
      !! do checkpoint test run  
      DO i = jstep0+1, nsteps
         j = iperm(i)
         WRITE(0,*) ' CPTEST  start from CP ', j
         call restore_checkpoint(j,jg,ocean_state,sea_ice, p_oce_sfc, p_atm_f, current_time, time_config, output_file)
         !      j =7
         jstep = j-1
         PRINT *,' DO ',jstep,' = ',jstep0+j-1, (jstep0+nsteps)-1, '     !!!!!!!!!!!!!!!  CHECK!!!'
         DO j = jstep0+j-1, (jstep0+nsteps)-1
            PRINT *, '            CPTEST tloop to ',jstep+1,' to  ',jstep0+nsteps
            WRITE(istr1,'(i010)') jstep+1
            tscnt(jstep+1) = tscnt(jstep+1)+1
            WRITE(istr2,'(i010)') tscnt(jstep+1)
            WRITE(*,*) 'CPTEST_BEGIN OCE-TimeStep-',TRIM(ADJUSTL(istr1)),'-',TRIM(ADJUSTL(istr2))
            call ocean_time_step()
            WRITE(*,*) 'CPTEST_END OCE-TimeStep-',TRIM(ADJUSTL(istr1)),'-',TRIM(ADJUSTL(istr2))
         ENDDO
         call store_checkpoints(jstep0+nsteps+2,jg,ocean_state,sea_ice,p_oce_sfc,p_atm_f, &
             & current_time,time_config,output_file)
         CALL compare_checkpoints( jg, jstep0+nsteps+1, jstep0+nsteps+2, j, ocean_state )
      END DO
#endif /* __COMPAD_CPTEST__ */


#ifdef __COMPAD_ADJ_IO_TEST__
      !! store number of performed forward steps
      nsteps = jstep
      
  !!$    CALL DCO_A1S_GET( ADJINP_h(:,:)         , ocean_state(jg)%p_adjoint%h(:,:)        , -1 )
  !!$    CALL DCO_A1S_GET( ADJINP_vn(:,:,:)      , ocean_state(jg)%p_adjoint%vn(:,:,:)      , -1 )
  !!$    CALL DCO_A1S_GET( ADJINP_tr(:,:,:,1)    , ocean_state(jg)%p_adjoint%t(:,:,:)    , -1 )
  !!$    CALL DCO_A1S_GET( ADJINP_tr(:,:,:,2)    , ocean_state(jg)%p_adjoint%s(:,:,:)    , -1 )
      DO i = nsteps-1, jstep0+1,-1
         ocean_state(jg)%p_adjoint%h(:,:)      =  -1 * i
         ocean_state(jg)%p_adjoint%vn(:,:,:)   =   2 * i
         ocean_state(jg)%p_adjoint%t(:,:,:)    =  -3 * i
         ocean_state(jg)%p_adjoint%s(:,:,:)    =   4 * i
      
         call restore_checkpoint(i,jg,ocean_state,sea_ice, p_oce_sfc, p_atm_f, current_time, time_config, output_file)
         ! write netcdf output
         CALL restartDescriptor%updatePatch(patch_2d, &
                                           &opt_nice_class=1, &
                                           &opt_ocean_zlevels=n_zlev, &
                                           &opt_ocean_zheight_cellmiddle = patch_3d%p_patch_1d(1)%zlev_m(:), &
                                           &opt_ocean_zheight_cellinterfaces = patch_3d%p_patch_1d(1)%zlev_i(:))
         CALL restartDescriptor%writeRestart(current_time, jstep)
         !CALL create_restart_file( patch = patch_2d,       &
         !     & current_time= i * ocean_time_step ,      &
         !     & jstep=i,            &
         !     & model_type="oce",       &
         !     & opt_sim_time=time_config%sim_time(1),&
         !     & opt_nice_class=1,       &
         !     & ocean_zlevels=n_zlev,                                         &
         !     & ocean_zheight_cellmiddle = patch_3d%p_patch_1d(1)%zlev_m(:),  &
         !     & ocean_zheight_cellinterfaces = patch_3d%p_patch_1d(1)%zlev_i(:))
      END DO
#endif /* __COMPAD_ADJ_IO_TEST__ */


