!>
!! This module is the interface between nwp_nh_interface to the radiation schemes
!! (RRTM or Ritter-Geleyn).
!!
!! @author Thorsten Reinhardt, AGeoBw, Offenbach
!!
!! @par Revision History
!! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_nwp_mpiomp_rrtm_interface

  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_exception,            ONLY: message,  finish !message_tex
  USE mo_ext_data,             ONLY: t_external_data
  USE mo_parallel_config,      ONLY: nproma, p_test_run, test_parallel_radiation
  USE mo_run_config,           ONLY: msg_level, iqv, iqc, iqi, &
    &                                io3, ntracer, ntracer_static
  USE mo_grf_interpolation,    ONLY: t_gridref_state
  USE mo_impl_constants,       ONLY: min_rlcell_int, icc!, min_rlcell 
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c, grf_ovlparea_start_c
  USE mo_interpolation,        ONLY: t_int_state
  USE mo_kind,                 ONLY: wp
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_nwp_lnd_state,        ONLY: t_lnd_prog, t_lnd_diag
  USE mo_model_domain,         ONLY: t_patch, p_patch_local_parent
  USE mo_mpi,                  ONLY: my_process_is_mpi_seq
  USE mo_phyparam_soil,        ONLY: csalb, csalb_snow_min, csalb_snow_max, &
    &                                csalb_snow_fe, csalb_snow_fd, csalb_p, cf_snow
  USE mo_phys_nest_utilities,  ONLY: upscale_rad_input, downscale_rad_output, &
    &                                upscale_rad_input_rg, downscale_rad_output_rg
  USE mo_nonhydro_state,       ONLY: t_nh_prog, t_nh_diag
  USE mo_nwp_phy_state,        ONLY: t_nwp_phy_diag !,prm_diag
  USE mo_radiation,            ONLY: radiation, pre_radiation_nwp_steps
  USE mo_radiation_config,     ONLY: irad_o3, irad_aero, vmr_co2, rad_csalbw
  USE mo_satad,                ONLY: qsat_rho
!   USE mo_sync,                 ONLY: SYNC_C, sync_patch_array_mult

  USE mo_nwp_rrtm_interface,   ONLY:  nwp_rrtm_radiation, nwp_rrtm_radiation_reduced

  USE mo_parallel_config,      ONLY: radiation_ompthreads
  USE mo_timer,                ONLY: timer_start, timer_stop, print_timer, &
    & new_timer, delete_timer, timer_omp_radiation

  USE mo_ompthreads

  
  IMPLICIT NONE


  PRIVATE
  
  PUBLIC :: nwp_start_radiation_ompthread, model_end_ompthread
  PUBLIC :: nwp_omp_rrtm_interface
  PUBLIC :: init_ompthread_radiation

  
  INTEGER :: radiation_ompthread_status, model_ompthread_status
  
  LOGICAL:: continue_radiation_ompthread = .true.
!  LOGICAL:: radiation_barier = .false.
  INTEGER:: no_radiation_barriers = 0
  INTEGER:: no_model_barriers = 0

  
  TYPE(t_ompthread) :: radiation_ompthread
  TYPE(t_ompthread) :: model_ompthread
  
  TYPE t_nwp_omp_radiation_data

    ! input
    TYPE(t_patch),        POINTER :: pt_patch     !<grid/patch info.
    TYPE(t_external_data),POINTER :: ext_data
    TYPE(t_lnd_diag),     POINTER :: lnd_diag      !< diag vars for sfc
    TYPE(t_nh_prog),      POINTER :: pt_prog_rcf !<the prognostic variables (with
    TYPE(t_nh_diag),      POINTER :: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag), POINTER :: prm_diag
    TYPE(t_lnd_prog),     POINTER :: lnd_prog_now
    
    REAL(wp)              :: p_sim_time
    
    REAL(wp), ALLOCATABLE::  fr_land_smt(:,:)   !< fraction land in a grid element        [ ]
                                                 !  = smoothed fr_land
    REAL(wp), ALLOCATABLE::  fr_glac_smt(:,:)   !< fraction land glacier in a grid element [ ]
                                                 ! = smoothed fr_glac
    REAL(wp), ALLOCATABLE::  cosmu0(:,:)        ! cosine of solar zenith angle
    REAL(wp), ALLOCATABLE::  emis_rad(:,:)      ! lw sfc emissivity
    REAL(wp), ALLOCATABLE::  tsfctrad(:,:)      ! surface temperature at trad [K]
    REAL(wp), ALLOCATABLE::  pres_ifc(:,:,:)    ! pressure at interfaces (nproma,nlevp1,nblks_c)  [Pa]
    REAL(wp), ALLOCATABLE::  pres(:,:,:)        ! pressure (nproma,nlev,nblks_c)                  [Pa]
    REAL(wp), ALLOCATABLE::  temp(:,:,:)        ! temperature (nproma,nlev,nblks_c)                 [K]
    REAL(wp), ALLOCATABLE::  tot_cld(:,:,:,:)   ! total cloud variables (cc,qv,qc,qi)
    REAL(wp), ALLOCATABLE::  qm_o3(:,:,:)       ! in o3 mass mixing ratio
    REAL(wp), ALLOCATABLE::  acdnc(:,:,:)       ! cloud droplet number concentration [1/m**3]

    ! output
    REAL(wp), ALLOCATABLE::  lwflxclr(:,:,:)    ! longwave clear-sky net flux [W/m2]
    REAL(wp), ALLOCATABLE::  trsolclr(:,:,:)    ! shortwave clear-sky net tranmissivity []
    REAL(wp), ALLOCATABLE::  lwflxall(:,:,:)    ! longwave net flux           [W/m2]
    REAL(wp), ALLOCATABLE::  trsolall(:,:,:)    ! shortwave net tranmissivity []
    
  END TYPE t_nwp_omp_radiation_data

  TYPE(t_nwp_omp_radiation_data):: omp_radiation_data
    

  CHARACTER(len=*), PARAMETER:: version = '$Id$'

  REAL(wp), PARAMETER::  &
    & zaeops = 0.05_wp,   &
    & zaeopl = 0.2_wp,    &
    & zaeopu = 0.1_wp,    &
    & zaeopd = 1.9_wp,    &
    & ztrpt  = 30.0_wp,   &
    & ztrbga = 0.03_wp  / (101325.0_wp - 19330.0_wp), &
    & zvobga = 0.007_wp /  19330.0_wp , &
    & zstbga = 0.045_wp /  19330.0_wp!, &
!      & zaeadk(1:3) = (/0.3876E-03_wp,0.6693E-02_wp,0.8563E-03_wp/)

CONTAINS

  !-----------------------------------------
  !>
  SUBROUTINE init_ompthread_radiation()

    INTEGER init_result
    
    init_result = init_ompthreads()
    radiation_ompthread = new_ompthread()
    model_ompthread = new_ompthread()
    
  END SUBROUTINE init_ompthread_radiation
  !-----------------------------------------


  !-----------------------------------------
  !>
  SUBROUTINE allocate_omp_radiation_data()

    INTEGER:: nblks_c, nlev, nlevp1
    INTEGER:: ist, total_status

    ! wait until omp_radiation_data%pt_patch is initialized
    ist=radiation_barrier()
    
    nblks_c = omp_radiation_data%pt_patch%nblks_c
!     write(0,*) "nproma,nblks_c=",nproma,nblks_c
    ! number of vertical levels
    nlev   = omp_radiation_data%pt_patch%nlev
    nlevp1 = omp_radiation_data%pt_patch%nlevp1
!     write(0,*) "nlev,nlevp1=",nlev,nlevp1

    total_status=0
    ! input    
    ALLOCATE(omp_radiation_data%fr_land_smt(nproma,nblks_c),STAT=ist)
    total_status=total_status+ist
    ALLOCATE(omp_radiation_data%fr_glac_smt(nproma,nblks_c),STAT=ist)
    total_status=total_status+ist
    ALLOCATE(omp_radiation_data%emis_rad(nproma,nblks_c),STAT=ist)
    total_status=total_status+ist
    ALLOCATE(omp_radiation_data%cosmu0(nproma,nblks_c), STAT=ist ) 
    total_status=total_status+ist
    ALLOCATE(omp_radiation_data%tsfctrad(nproma,nblks_c), STAT=ist ) 
    total_status=total_status+ist
    ALLOCATE(omp_radiation_data%pres_ifc(nproma,nlevp1,nblks_c), STAT=ist )
    total_status=total_status+ist
    ALLOCATE(omp_radiation_data%pres(nproma,nlev,nblks_c), STAT=ist )
    total_status=total_status+ist
    ALLOCATE(omp_radiation_data%temp(nproma,nlev,nblks_c), STAT=ist )
    total_status=total_status+ist
    ALLOCATE(omp_radiation_data%tot_cld(nproma,nlev,nblks_c,4), STAT=ist )
    total_status=total_status+ist
    ALLOCATE(omp_radiation_data%qm_o3(nproma,nlev,nblks_c), STAT=ist)
    total_status=total_status+ist
    ALLOCATE(omp_radiation_data%acdnc(nproma,nlev,nblks_c), STAT=ist ) 
    total_status=total_status+ist

    ! output
    ALLOCATE(omp_radiation_data%lwflxclr(nproma,nlevp1,nblks_c), STAT=ist ) 
    total_status=total_status+ist
    ALLOCATE(omp_radiation_data%trsolclr(nproma,nlevp1,nblks_c), STAT=ist ) 
    total_status=total_status+ist
    ALLOCATE(omp_radiation_data%lwflxall(nproma,nlevp1,nblks_c), STAT=ist ) 
    total_status=total_status+ist
    ALLOCATE(omp_radiation_data%trsolall(nproma,nlevp1,nblks_c), STAT=ist ) 
    total_status=total_status+ist
    
!     write(0,*) "allocate omp_radiation_data is done!"
  
    IF (total_status /= 0 ) THEN
      CALL finish ("allocate_omp_radiation_data",'allocating omp_radiation_data failed')
    ENDIF

!     write(0,*) "zero in..."
    omp_radiation_data%p_sim_time = 0.0_wp
    omp_radiation_data%fr_land_smt(:,:) = 0.0_wp   !< fraction land in a grid element        [ ]
                                                 !  = smoothed fr_land
    omp_radiation_data%fr_glac_smt(:,:) = 0.0_wp   !< fraction land glacier in a grid element [ ]
                                                 ! = smoothed fr_glac
    omp_radiation_data%cosmu0(:,:) = 0.0_wp        ! cosine of solar zenith angle
    omp_radiation_data%emis_rad(:,:) = 0.0_wp      ! lw sfc emissivity
    omp_radiation_data%tsfctrad(:,:) = 0.0_wp      ! surface temperature at trad [K]
    omp_radiation_data%pres_ifc(:,:,:) = 0.0_wp    ! pressure at interfaces (nproma,nlevp1,nblks_c)  [Pa]
    omp_radiation_data%pres(:,:,:)  = 0.0_wp      ! pressure (nproma,nlev,nblks_c)                  [Pa]
    omp_radiation_data%temp(:,:,:)  = 0.0_wp       ! temperature (nproma,nlev,nblks_c)                 [K]
    omp_radiation_data%tot_cld(:,:,:,:) = 0.0_wp   ! total cloud variables (cc,qv,qc,qi)
    omp_radiation_data%qm_o3(:,:,:)  = 0.0_wp      ! in o3 mass mixing ratio
    omp_radiation_data%acdnc(:,:,:) = 0.0_wp       ! cloud droplet number concentration [1/m**3]

!     write(0,*) "zero out..."
    ! output
    omp_radiation_data%lwflxclr(:,:,:) = 0.0_wp   ! longwave clear-sky net flux [W/m2]
    omp_radiation_data%trsolclr(:,:,:) = 0.0_wp    ! shortwave clear-sky net tranmissivity []
    omp_radiation_data%lwflxall(:,:,:) = 0.0_wp   ! longwave net flux           [W/m2]
    omp_radiation_data%trsolall(:,:,:) = 0.0_wp   ! shortwave net tranmissivity []
    
    ! notify that omp_radiation_data is allocated
    ist=radiation_barrier()
    
!     write(0,*) "allocate_omp_radiation_data exits"
  END SUBROUTINE allocate_omp_radiation_data
  !-----------------------------------------

    
  !-----------------------------------------
  !>
  ! This is the radiation thread barrier
  INTEGER FUNCTION radiation_barrier()
    
    write(0,*) "radiation_barrier starts..."
    radiation_barrier = radiation_ompthread.syncto.model_ompthread
    
    IF (ompthread_has_endrequest(radiation_ompthread)) THEN
      write(0,*) "Reached stop_radiation_ompthread."
      CALL end_ompthread(radiation_ompthread)
    ENDIF
    
    no_radiation_barriers = no_radiation_barriers + 1
    write(0,*) "no_radiation_barriers=", no_radiation_barriers
         
    RETURN
    
  END FUNCTION radiation_barrier
  !-----------------------------------------


  !-----------------------------------------
  !>
  ! This is the model thread barrier
  INTEGER FUNCTION model_barrier()

    write(0,*) "model_barrier starts..."
    model_barrier = model_ompthread.syncto.radiation_ompthread
    
    no_model_barriers = no_model_barriers + 1
    write(0,*) "no_model_barriers=", no_model_barriers
      
  END FUNCTION model_barrier
  !-----------------------------------------

  !-----------------------------------------
  !>
  SUBROUTINE model_end_ompthread()

    INTEGER :: sync_status
    
    CALL request_end_ompthread(model_ompthread)
    write(0,*) 'model_ompthread requests end'
    
    sync_status = model_ompthread.syncto.radiation_ompthread
    
    CALL end_ompthread(model_ompthread)
    write(0,*) 'model_ompthread ends'

  END SUBROUTINE model_end_ompthread
  !-----------------------------------------

  !-----------------------------------------
  !>
  SUBROUTINE receive_in_omp_radiation_data()

    INTEGER :: wait_cnt
!     write(0,*) 'receive_in_omp_radiation_data starts'
    ! start receving input
    wait_cnt=radiation_barrier()
    ! end of receive input
    wait_cnt=radiation_barrier()
!     write(0,*) 'receive_in_omp_radiation_data exits'
  
  END SUBROUTINE receive_in_omp_radiation_data
  !-----------------------------------------


  !-----------------------------------------
  !>
  SUBROUTINE send_out_omp_radiation_data()
  
    INTEGER :: wait_cnt
    ! write(0,*) 'send_out_omp_radiation_data starts'
    ! start receving input
    wait_cnt=radiation_barrier()
    ! end of receive input
    wait_cnt=radiation_barrier()
!     write(0,*) 'send_out_omp_radiation_data exits'
  
  END SUBROUTINE send_out_omp_radiation_data
  !-----------------------------------------

  !-----------------------------------------
  !>
  ! Send the input omp_radiation_data
  SUBROUTINE send_in_omp_radiation_data(p_sim_time )
    
    REAL(wp), INTENT(in)  :: p_sim_time
    INTEGER :: wait_cnt
    
!     write(0,*) 'send_in_omp_radiation_data starts'
    wait_cnt=model_barrier()
!$OMP PARALLEL WORKSHARE
    omp_radiation_data%p_sim_time       = p_sim_time
    omp_radiation_data%fr_land_smt(:,:) = omp_radiation_data%ext_data%atm%fr_land_smt(:,:)
    omp_radiation_data%fr_glac_smt(:,:) = omp_radiation_data%ext_data%atm%fr_glac_smt(:,:)
    omp_radiation_data%emis_rad(:,:)    = omp_radiation_data%ext_data%atm%emis_rad(:,:)
    omp_radiation_data%tsfctrad(:,:)    = omp_radiation_data%prm_diag%tsfctrad(:,:)
    omp_radiation_data%tsfctrad(:,:)    = omp_radiation_data%lnd_prog_now%t_g(:,:)
    ! fill also the prm_diag%tsfctrad(:,:) !
    omp_radiation_data%prm_diag%tsfctrad(:,:) = omp_radiation_data%lnd_prog_now%t_g(:,:)
    omp_radiation_data%pres_ifc(:,:,:)  = omp_radiation_data%pt_diag%pres_ifc(:,:,:)
    omp_radiation_data%pres(:,:,:)      = omp_radiation_data%pt_diag%pres(:,:,:)
    omp_radiation_data%temp(:,:,:)      = omp_radiation_data%pt_diag%temp(:,:,:)
    omp_radiation_data%tot_cld(:,:,:,:) = omp_radiation_data%prm_diag%tot_cld(:,:,:,:)
    omp_radiation_data%qm_o3(:,:,:)     = omp_radiation_data%pt_prog_rcf%tracer(:,:,:,io3)
    omp_radiation_data%acdnc(:,:,:)     = omp_radiation_data%prm_diag%acdnc(:,:,:)
!$OMP END PARALLEL WORKSHARE
    wait_cnt=model_barrier()
!     write(0,*) 'send_in_omp_radiation_data ends'

  END SUBROUTINE send_in_omp_radiation_data
  !-----------------------------------------
  
  !-----------------------------------------
  !>
  ! Receive the output omp_radiation_data
  SUBROUTINE receive_out_omp_radiation_data()
        
    INTEGER :: wait_cnt
!     write(0,*) 'receive_out_omp_radiation_data starts'
    wait_cnt=model_barrier()

    
    IF (test_parallel_radiation) THEN
      ! compare to the sequential version
!!$      CALL nwp_rrtm_radiation ( omp_radiation_data%p_sim_time,omp_radiation_data%pt_patch, &
!!$        & omp_radiation_data%ext_data,omp_radiation_data%lnd_diag,omp_radiation_data%pt_prog_rcf,&
!!$        & omp_radiation_data%pt_diag,omp_radiation_data%prm_diag, omp_radiation_data%lnd_prog_now )

      IF (MAXVAL(ABS(omp_radiation_data%prm_diag%lwflxclr(:,:,:) &
                  - omp_radiation_data%lwflxclr(:,:,:))) /= 0.0_wp) THEN
        CALL finish("receive_out_omp_radiation_data","lwflxclr differs")
      ENDIF

      IF (MAXVAL(ABS(omp_radiation_data%prm_diag%trsolclr(:,:,:) &
                  - omp_radiation_data%trsolclr(:,:,:))) /= 0.0_wp) THEN
        CALL finish("receive_out_omp_radiation_data","trsolclr differs")
      ENDIF

      IF (MAXVAL(ABS(omp_radiation_data%prm_diag%lwflxall(:,:,:) &
                  - omp_radiation_data%lwflxall(:,:,:))) /= 0.0_wp) THEN
        CALL finish("receive_out_omp_radiation_data","lwflxall differs")
      ENDIF

      IF (MAXVAL(ABS(omp_radiation_data%prm_diag%trsolall(:,:,:) &
                  - omp_radiation_data%trsolall(:,:,:))) /= 0.0_wp) THEN
        CALL finish("receive_out_omp_radiation_data","trsolall differs")
      ENDIF
    ENDIF
    
!$OMP PARALLEL WORKSHARE
    omp_radiation_data%prm_diag%lwflxclr(:,:,:) =omp_radiation_data%lwflxclr(:,:,:)
    omp_radiation_data%prm_diag%trsolclr(:,:,:) =omp_radiation_data%trsolclr(:,:,:)
    omp_radiation_data%prm_diag%lwflxall(:,:,:) =omp_radiation_data%lwflxall(:,:,:)
    omp_radiation_data%prm_diag%trsolall(:,:,:) =omp_radiation_data%trsolall(:,:,:)
!$OMP END PARALLEL WORKSHARE
    wait_cnt=model_barrier()
!     write(0,*) 'receive_out_omp_radiation_data ends'

  END SUBROUTINE receive_out_omp_radiation_data
  !-----------------------------------------
  
  !-----------------------------------------
  !>
  SUBROUTINE nwp_start_radiation_ompthread()
!     INTEGER :: timer_omp_radiation!, timer_omp_model, timer_omp_model_waits

!$  INTEGER omp_get_num_threads

!$  CALL omp_set_num_threads(radiation_ompthreads)

!$     write(0,*) 'Entering nwp_start_radiation_ompthread, threads=',&
!$       omp_get_num_threads()
!$OMP PARALLEL
!$    write(0,*) 'nwp_start_radiation_ompthread parallel, threads=',&
!$      omp_get_num_threads()
!$OMP END PARALLEL

    !--------------------------------------
    ! first call of the radiation ompthread
    ! we need to allocate the omp_radiation_data

!     timer_omp_radiation = new_timer("omp_radiation")
!     timer_omp_model = new_timer("omp_model")
!     timer_omp_model_waits = new_timer("omp_model_waits")
    
    CALL allocate_omp_radiation_data()
    
    ! receive input data
    CALL receive_in_omp_radiation_data()

    ! calculate radiation
    CALL nwp_rrtm_radiation_ompthread()

    ! sent the result
    CALL send_out_omp_radiation_data()
    
    !---------------------------
    !  radiation loop
    !  loop until the stepping is done
    
    DO WHILE (ompthread_is_alive(radiation_ompthread))
      ! calculate radiation
      ! the first itearation will calculate with the same input
      ! but extended radiation timestep 
      CALL timer_start(timer_omp_radiation)
      CALL nwp_rrtm_radiation_ompthread()
      CALL timer_stop(timer_omp_radiation)

      
! #ifndef __TEST_OMP_RADIATION__
!       IF (test_parallel_radiation) THEN
!       CALL timer_start(timer_omp_radiation)
!       CALL nwp_rrtm_radiation_ompthread()
!       CALL timer_stop(timer_omp_radiation)
! #endif

      CALL receive_in_omp_radiation_data()

      IF (ompthread_is_dead(radiation_ompthread)) EXIT
            
      CALL send_out_omp_radiation_data()

    ENDDO
    !---------------------------

    write(0,*) 'Leaving nwp_start_radiation_ompthread'
!     CALL print_timer(timer_omp_radiation)
!     CALL delete_timer(timer_omp_radiation)
    
  END SUBROUTINE nwp_start_radiation_ompthread
  


  !---------------------------------------------------------------------------------------
  !>
  SUBROUTINE nwp_omp_rrtm_interface ( p_sim_time,pt_patch, &
    & ext_data,lnd_diag,pt_prog_rcf,pt_diag,prm_diag, &
    & lnd_prog_now )

!    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
!      &  routine = 'mo_nwp_rad_interface:'

    REAL(wp),INTENT(in)         :: p_sim_time

    TYPE(t_patch),        TARGET, INTENT(in)    :: pt_patch     !<grid/patch info.
    TYPE(t_external_data),TARGET, INTENT(in)    :: ext_data
    TYPE(t_lnd_diag),     TARGET, INTENT(in)    :: lnd_diag    !< diag vars for sfc
    TYPE(t_nh_prog),      TARGET, INTENT(inout) :: pt_prog_rcf !<the prognostic variables (with
    TYPE(t_nh_diag),      TARGET, INTENT(in)    :: pt_diag     !<the diagnostic variables
    TYPE(t_nwp_phy_diag), TARGET, INTENT(inout) :: prm_diag
    TYPE(t_lnd_prog),     TARGET, INTENT(inout) :: lnd_prog_now
    
    INTEGER :: wait_cnt

    write(0,*) "Entering nwp_omp_rrtm_interface..."
    IF (no_model_barriers == 0) THEN
      ! this is the first call to the nwp_omp_rrtm_interface
      ! first pass the patch info and wait to allocate the buffers
      omp_radiation_data%pt_patch => pt_patch
      omp_radiation_data%ext_data =>  ext_data
      omp_radiation_data%lnd_diag =>  lnd_diag
      omp_radiation_data%pt_prog_rcf  =>  pt_prog_rcf
      omp_radiation_data%pt_diag =>  pt_diag
      omp_radiation_data%prm_diag =>  prm_diag
      omp_radiation_data%lnd_prog_now =>  lnd_prog_now
      
      wait_cnt=model_barrier()
      ! wait until the omp_radiation_data is allocated      
      wait_cnt=model_barrier()
      
!     ELSE
!       CALL timer_stop(timer_omp_model)    
    ENDIF

!     CALL timer_start(timer_omp_model_waits)
    ! now sent the input data
    CALL send_in_omp_radiation_data(p_sim_time)

    ! and receive the output data
    CALL receive_out_omp_radiation_data()    
!     CALL timer_stop(timer_omp_model_waits)
!     CALL timer_start(timer_omp_model)

    write(0,*) "Leaving nwp_omp_rrtm_interface..."
    
  END SUBROUTINE nwp_omp_rrtm_interface
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  !>
  SUBROUTINE nwp_rrtm_radiation_ompthread ( )

    REAL(wp), PARAMETER::  &
      & cosmu0_dark =  -1.e-9_wp  ! minimum cosmu0, for smaller values no shortwave calculations
    
    REAL(wp):: albvisdir     (nproma,omp_radiation_data%pt_patch%nblks_c) !<
    REAL(wp):: albnirdir     (nproma,omp_radiation_data%pt_patch%nblks_c) !<
    REAL(wp):: albvisdif     (nproma,omp_radiation_data%pt_patch%nblks_c) !<
    REAL(wp):: albnirdif     (nproma,omp_radiation_data%pt_patch%nblks_c) !<
    REAL(wp):: aclcov        (nproma,omp_radiation_data%pt_patch%nblks_c) !<

    REAL(wp):: zaeq1(nproma,omp_radiation_data%pt_patch%nlev,omp_radiation_data%pt_patch%nblks_c)
    REAL(wp):: zaeq2(nproma,omp_radiation_data%pt_patch%nlev,omp_radiation_data%pt_patch%nblks_c)
    REAL(wp):: zaeq3(nproma,omp_radiation_data%pt_patch%nlev,omp_radiation_data%pt_patch%nblks_c)
    REAL(wp):: zaeq4(nproma,omp_radiation_data%pt_patch%nlev,omp_radiation_data%pt_patch%nblks_c)
    REAL(wp):: zaeq5(nproma,omp_radiation_data%pt_patch%nlev,omp_radiation_data%pt_patch%nblks_c)

    INTEGER :: itype(nproma)   !< type of convection

    ! Local scalars:
    REAL(wp):: zsct        ! solar constant (at time of year)
    INTEGER:: jc,jk,jb
    INTEGER:: jg                !domain id
    INTEGER:: nlev, nlevp1      !< number of full and half levels

    INTEGER:: rl_start, rl_end
    INTEGER:: i_startblk, i_endblk    !> blocks
    INTEGER:: i_startidx, i_endidx    !< slices
    INTEGER:: i_nchdom                !< domain index
    INTEGER:: i_chidx
!     LOGICAL:: l_parallel

    i_nchdom  = MAX(1,omp_radiation_data%pt_patch%n_childdom)
    jg        = omp_radiation_data%pt_patch%id

    ! number of vertical levels
    nlev   = omp_radiation_data%pt_patch%nlev
    nlevp1 = omp_radiation_data%pt_patch%nlevp1

    !-------------------------------------------------------------------------
    !> Radiation setup
    !-------------------------------------------------------------------------

    ! Calculation of zenith angle optimal during dt_rad.
    ! (For radheat, actual zenith angle is calculated separately.)
    CALL pre_radiation_nwp_steps (                        &
      & kbdim        = nproma,                            &
      & cosmu0_dark  = cosmu0_dark,                       &
      & p_inc_rad    = atm_phy_nwp_config(jg)%dt_rad,     &
      & p_inc_radheat= atm_phy_nwp_config(jg)%dt_radheat, &
      & p_sim_time   = omp_radiation_data%p_sim_time,     &
      & pt_patch     = omp_radiation_data%pt_patch,       &
      & zsmu0        = omp_radiation_data%cosmu0(:,:),    &
      & zsct         = zsct )


    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------


    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = omp_radiation_data%pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = omp_radiation_data%pt_patch%cells%end_blk(rl_end,i_nchdom)

    IF (msg_level >= 12) &
      &           CALL message('mo_nwp_rad_interface', 'RRTM radiation on full grid')

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,itype),SCHEDULE(guided)
    !
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(omp_radiation_data%pt_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)


      ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM
      itype(1:i_endidx) = 0 !INT(field%rtype(1:i_endidx,jb))

      albvisdir(1:i_endidx,jb) = 0.07_wp ! ~ albedo of water
      albnirdir(1:i_endidx,jb) = 0.07_wp ! ~ albedo of water
      albvisdif(1:i_endidx,jb) = 0.07_wp ! ~ albedo of water
      albnirdif(1:i_endidx,jb) = 0.07_wp ! ~ albedo of water
      zaeq1(1:i_endidx,:,jb)   = 0.0_wp
      zaeq2(1:i_endidx,:,jb)   = 0.0_wp
      zaeq3(1:i_endidx,:,jb)   = 0.0_wp
      zaeq4(1:i_endidx,:,jb)   = 0.0_wp
      zaeq5(1:i_endidx,:,jb)   = 0.0_wp

      CALL radiation(               &
                              !
                              ! input
                              ! -----
                              !
                              ! indices and dimensions
        & jce        =i_endidx             ,&!< in  end   index for loop over block
        & kbdim      =nproma               ,&!< in  dimension of block over cells
        & klev       =nlev                 ,&!< in  number of full levels = number of layers
        & klevp1     =nlevp1               ,&!< in  number of half levels = number of layer ifcs
                              !
        & ktype      =itype                ,&!< in     type of convection
                              !
                              ! surface: albedo + temperature
        & zland      =omp_radiation_data%fr_land_smt(:,jb)   ,&!< in     land fraction
        & zglac      =omp_radiation_data%fr_glac_smt(:,jb)   ,&!< in     land glacier fraction
                              !
        & cos_mu0    =omp_radiation_data%cosmu0  (:,jb) ,&!< in  cos of zenith angle mu0
        & alb_vis_dir=albvisdir        (:,jb) ,&!< in surface albedo for visible range, direct
        & alb_nir_dir=albnirdir        (:,jb) ,&!< in surface albedo for near IR range, direct
        & alb_vis_dif=albvisdif        (:,jb) ,&!< in surface albedo for visible range, diffuse
        & alb_nir_dif=albnirdif        (:,jb) ,&!< in surface albedo for near IR range, diffuse
        & emis_rad   =omp_radiation_data%emis_rad(:,jb), & !< lw sfc emissivity
        & tk_sfc     =omp_radiation_data%tsfctrad(:,jb) ,&!< in surface temperature
                              !
                              ! atmosphere: pressure, tracer mixing ratios and temperature
        & pp_hl  =omp_radiation_data%pres_ifc(:,:,jb)    ,&!< in pres at half levels at t-dt [Pa]
        & pp_fl  =omp_radiation_data%pres    (:,:,jb)    ,&!< in pres at full levels at t-dt [Pa]
        & tk_fl  =omp_radiation_data%temp    (:,:,jb)    ,&!< in temperature at full level at t-dt
        & qm_vap =omp_radiation_data%tot_cld (:,:,jb,iqv),&!< in water vapor mass mix ratio at t-dt
        & qm_liq =omp_radiation_data%tot_cld (:,:,jb,iqc),&!< in cloud water mass mix ratio at t-dt
        & qm_ice =omp_radiation_data%tot_cld (:,:,jb,iqi),&!< in cloud ice mass mixing ratio at t-dt
        & qm_o3  =omp_radiation_data%qm_o3   (:,:,jb) ,   &!< in o3 mass mixing ratio at t-dt
        & cdnc   =omp_radiation_data%acdnc   (:,:,jb)    ,&!< in cloud droplet numb conc. [1/m**3]
        & cld_frc=omp_radiation_data%tot_cld (:,:,jb,icc),&!< in cloud fraction [m2/m2]
        & zaeq1   = zaeq1(:,:,jb)                        ,&!< in aerosol continental
        & zaeq2   = zaeq2(:,:,jb)                        ,&!< in aerosol maritime
        & zaeq3   = zaeq3(:,:,jb)                        ,&!< in aerosol urban
        & zaeq4   = zaeq4(:,:,jb)                        ,&!< in aerosol volcano ashes
        & zaeq5   = zaeq5(:,:,jb)                        ,&!< in aerosol stratospheric background
                              !
                              ! output
                              ! ------
                              !
        & cld_cvr    =aclcov             (:,jb),&!< out cloud cover in a column [m2/m2]
        & emter_clr  =omp_radiation_data%lwflxclr(:,:,jb),&!< out terrestrial flux, clear sky, net down
        & trsol_clr  =omp_radiation_data%trsolclr(:,:,jb),&!< out sol. transmissivity, clear sky, net down
        & emter_all  =omp_radiation_data%lwflxall(:,:,jb),&!< out terrestrial flux, all sky, net down
        & trsol_all  =omp_radiation_data%trsolall(:,:,jb),&!< out solar transmissivity, all sky, net down
        & opt_halo_cosmu0 = .FALSE. )
      ENDDO ! blocks

!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE nwp_rrtm_radiation_ompthread
  !---------------------------------------------------------------------------------------

  

END MODULE mo_nwp_mpiomp_rrtm_interface
! 
!  ----------------------------------------------------------------------------------------------------
!  Timer report:
!  calls  t_min       t_average   t_max       t_total
!  ----------------------------------------------------------------------------------------------------
!  0: total                              1      02m09s      02m09s      02m09s      02m09s   129.831
!  0: solve_nh                          34     .40207s      1.848s      4.250s      01m02s    62.842
!  0: physics                            9      2.074s      4.575s     19.706s     41.172s    41.172
!  0: nwp_radiation                      5     .02977s      3.777s     18.763s     18.885s    18.885
!  0: physic_acc                         9     .14167s     .17915s     .26425s      1.612s     1.612
!  0: physic_acc_2                       9     .08113s     .09902s     .13016s     .89114s     0.891
!  0: nwp_turbulence                     8     .39173s     .50518s     .69060s      4.041s     4.041
!  0: nwp_microphysics                   8     .18578s     .24390s     .39189s      1.951s     1.951
!  0: phys_sync_patch                    9     .00000s     .00000s     .00000s     .00001s     0.000
!  0: fast_phys                          8     .12697s     .17952s     .25356s      1.436s     1.436
!  0: pre_radiation_nwp                  9     .00030s     .00032s     .00036s     .00288s     0.003
!  0: intp                               1     .01288s     .01288s     .01288s     .01288s     0.013
!  0: transport                          8      2.167s      2.955s      4.312s     23.644s    23.644
!  ----------------------------------------------------------------------------------------------------
 
