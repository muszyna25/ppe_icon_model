!>
!!  Initialize the physical schemes at start time
!!
!! @author Kristina  Froehlich, DWD
!!
!! @par Revision History
!! First implementations by Kristina Froehlich, DWD, 2010-070-20
!! Include initialitation of SST for APE experiments by P. Ripodas, DWD,2010-11
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_nwp_phy_init

  USE mo_kind,               ONLY: wp
  USE mo_math_constants,     ONLY: pi
  USE mo_physical_constants, ONLY: re
  USE mo_nwp_phy_state,      ONLY: t_nwp_phy_diag,t_nwp_phy_tend
  USE mo_nwp_lnd_state,      ONLY: t_lnd_prog,t_lnd_diag
  USE mo_ext_data,           ONLY: t_external_data
  USE mo_nonhydro_state,     ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_exception,          ONLY: message, finish,message_text
  USE mo_vertical_coord_table,ONLY: vct_a, vct
  USE mo_model_domain,       ONLY: t_patch
  USE mo_model_domain_import,ONLY: nroot 
  USE mo_impl_constants,     ONLY: min_rlcell
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_run_nml,            ONLY: ltestcase, iqv, iqc,       &
       &                           iqt, ntracer, nproma, &
       &                           msg_level
  USE mo_atm_phy_nwp_nml,    ONLY: inwp_gscp, inwp_convection,&
       &                           inwp_radiation, inwp_turb,inwp_surface
  !radiation
  USE mo_newcld_optics,        ONLY: setup_newcld_optics
  USE mo_lrtm_setup,           ONLY: lrtm_setup
  USE mo_radiation_nml,        ONLY: ssi, tsi, irad_aero
  USE mo_srtm_config,          ONLY: setup_srtm, ssi_amip
  USE mo_radiation_rg_par,     ONLY: rad_aibi, init_aerosol, zaef
  ! microphysics
  USE mo_gscp_cosmo,           ONLY: hydci_pp_init
  ! convection
  USE mo_cuparameters,         ONLY: sucst,  sucumf,    &
    &                                su_yoethf,         &
    &                                sucldp, suphli,    &
    &                                suvdf , suvdfs
  USE mo_convect_tables,       ONLY: init_convect_tables
  !turbulence
!  USE mo_turbdiff_ras,         ONLY: init_canopy, organize_turbdiff
  USE src_turbdiff,            ONLY: init_canopy, organize_turbdiff
! for APE_nh experiments

  ! air-sea-land interface
  USE mo_icoham_sfc_indices,   ONLY: nsfc_type, iwtr, iice, ilnd !, &
 !   &                                init_sfc_indices
 ! vertical diffusion
  USE mo_echam_vdiff_params,   ONLY: init_vdiff_params, z0m_min, &
    &                                z0m_oce , tke_min 
  USE mo_vdiff_solver,         ONLY: init_vdiff_solver


  USE mo_satad,                ONLY: sat_pres_water, &  !! saturation vapor pressure w.r.t. water
    &                                spec_humi,qsat_rho !! Specific humidity

  USE mo_nh_testcases,         ONLY: nh_test_name, ape_sst_case
  USE mo_ape_params,           ONLY: ape_sst

  IMPLICIT NONE

  PUBLIC  :: init_nwp_phy

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS

SUBROUTINE init_nwp_phy ( pdtime                         , &
                       &  p_patch, p_metrics,              &
                       &  p_prog_now,  p_prog,  p_diag,    &
                       &  prm_diag,prm_nwp_tend,           &
                       &  p_prog_lnd_now, p_prog_lnd_new,  &
                       &  p_diag_lnd,                      &
                       &  ext_data, mean_charlen )

  TYPE(t_patch),        TARGET,INTENT(in)    :: p_patch
  TYPE(t_nh_metrics),          INTENT(in)    :: p_metrics
  TYPE(t_nh_prog),      TARGET,INTENT(inout) :: p_prog_now !!the prognostic variables
  TYPE(t_nh_prog),      TARGET,INTENT(inout) :: p_prog  !!the prognostic variables
  TYPE(t_nh_diag),      TARGET,INTENT(inout) :: p_diag  !!the diagostic variables
  TYPE(t_external_data),       INTENT(in)    :: ext_data
  TYPE(t_nwp_phy_diag),        INTENT(inout) :: prm_diag
  TYPE(t_nwp_phy_tend), TARGET,INTENT(inout) :: prm_nwp_tend
  TYPE(t_lnd_prog),            INTENT(inout) :: p_prog_lnd_now, p_prog_lnd_new
  TYPE(t_lnd_diag),            INTENT(inout) :: p_diag_lnd

  REAL(wp),INTENT(OUT)::  mean_charlen
  INTEGER             :: jk , nsmax
  REAL(wp)            :: pdtime
  REAL(wp)            :: pref(p_patch%nlevp1)
  REAL(wp), PARAMETER :: h_scal = 10000._wp    ! [m]      scale height
  REAL(wp), PARAMETER :: p0sl   = 101325._wp   ! [Pa]     sea level pressure
  REAL(wp)            :: zlat, zprat, zn1, zn2, zcdnc
  
  LOGICAL  :: lland, lglac
  
  INTEGER :: jb,jc
  INTEGER :: nlev, nlevp1            !< number of full and half levels
!  INTEGER :: jg
  INTEGER :: rl_start, rl_end
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !! slices
  INTEGER :: i_nchdom                !! domain index
  INTEGER  :: nexp

  INTEGER :: ierrstat=0
  CHARACTER (LEN=25) :: eroutine=''
  CHARACTER (LEN=80) :: errormsg=''

  INTEGER :: khydromet, ktrac
!>JH  
!  REAL(wp), DIMENSION(9):: czmls=(/ 0.,0.005,0.02,0.06,0.18,0.54,1.62,4.86,14.58 /)
!<JH
    i_nchdom  = MAX(1,p_patch%n_childdom)

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1


    IF ( ltestcase )THEN 

    rl_start = 1 ! Initialization should be done for all points
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)
    
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
             &  i_startidx, i_endidx, rl_start, rl_end)

        IF( nh_test_name == 'APE_nh') THEN

          ! t_g = ape_sst1
          
          DO jc = i_startidx, i_endidx
            zlat = p_patch%cells%center(jc,jb)%lat
            p_prog_lnd_now%t_g (jc,jb) = ape_sst(ape_sst_case,zlat) ! set SST
            p_prog_lnd_new%t_g (jc,jb) = ape_sst(ape_sst_case,zlat) ! set SST
            ! Humidity at water surface = humidity at saturation
            p_diag_lnd%qv_s    (jc,jb) = &
  !        & qsat_rho(p_prog_lnd_now%t_g (jc,jb),p_prog%rho(jc,nlev,jb))
          & spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))
          END DO

        ELSE ! any other testcase

          ! t_g  =  t(nlev)
          ! qv_ s= qv(nlev)
          ! KF increase the surface values to obtain fluxes          

          DO jc = i_startidx, i_endidx
            p_prog_lnd_now%t_g (jc,jb) = p_diag%temp  (jc,nlev,jb)!+0.2_wp
            p_prog_lnd_new%t_g (jc,jb) = p_diag%temp  (jc,nlev,jb)!+0.2_wp
            ! KF NOTE: as long as we have only water as lower boundary
            ! this is the same setting as for APE
           p_diag_lnd%qv_s    (jc,jb) = &
!                & qsat_rho(p_prog_lnd_now%t_g (jc,jb),p_prog%rho(jc,nlev,jb))
          & spec_humi(sat_pres_water(p_prog_lnd_now%t_g (jc,jb)),p_diag%pres_sfc(jc,jb))

          END DO
        ENDIF
        
      END DO
        CALL message('mo_nwp_phy_init:', 'initialized surface temp and humidity')
    END IF

    !--------------------------------------------------------------
    !< characteristic gridlength needed by convection and turbulence
    !--------------------------------------------------------------
      CALL mean_domain_values (p_patch,mean_charlen)

    !------------------------------------------
    !< call for cloud microphysics
    !------------------------------------------
  IF ( inwp_gscp == 1 )  THEN

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init microphysics')
    CALL hydci_pp_init

  ENDIF
    !------------------------------------------
    !< radiation
    !------------------------------------------
  IF ( inwp_radiation == 1 ) THEN

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init RRTM')

    ! Note (GZ): irad_aero=2 does no action but is the default in radiation_nml
    ! and therefore should not cause the model to stop
    IF ( irad_aero /= 0 .AND. irad_aero /= 2 ) THEN
      CALL finish('mo_nwp_phy_init: init_nwp_phy',  &
        &         'Wrong irad_aero. For RRTM, currently only irad_aero=0 is implemented.')
    ENDIF
    
!    prm_diag%lfglac (:,:) = ext_data%atm%soiltyp(:,:) == 1  !soiltyp=ice

    ssi(:) = ssi_amip(:)
    tsi    = SUM(ssi(:))

    IF ( nh_test_name == 'APE_nh' ) THEN
      ssi(:) = ssi(:)*1365._wp/tsi
      tsi = 1365._wp
    ENDIF
    
    CALL setup_srtm

    CALL lrtm_setup

    CALL setup_newcld_optics
    
    rl_start = 1  ! Initialization should be done for all points
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)
    
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx, i_endidx,jk)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
           &  i_startidx, i_endidx, rl_start, rl_end)

      ! Initialize cloud droplet number concentration (acdnc)
      ! like in mo_echam_phy_init.f90
      nexp=2
      DO jk = 1,nlev
        DO jc = i_startidx, i_endidx
          zprat=(MIN(8._wp,80000._wp/p_diag%pres(jc,jk,jb)))**nexp
          lland = ext_data%atm%lsm_atm_c(jc,jb) > 0 !
          lglac = ext_data%atm%soiltyp(jc,jb) == 1
          IF (lland.AND.(.NOT.lglac)) THEN
            zn1= 50._wp
            zn2=220._wp
          ELSE
            zn1= 50._wp
            zn2= 80._wp
          ENDIF
          IF (p_diag%pres(jc,jk,jb) < 80000._wp) THEN
            zcdnc=1.e6_wp*(zn1+(zn2-zn1)*(EXP(1._wp-zprat)))
          ELSE
            zcdnc=zn2*1.e6_wp
          ENDIF
          prm_diag%acdnc(jc,jk,jb) = zcdnc
        END DO !jc
      END DO   !jk
    ENDDO      !jb
!$OMP END DO
!$OMP END PARALLEL    
    
  ELSEIF ( inwp_radiation == 2 ) THEN

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init Ritter Geleyn')

    ! Note (GZ): irad_aero=2 does no action but is the default in radiation_nml
    ! and therefore should not cause the model to stop
    SELECT CASE ( irad_aero )
    CASE (0,2,5)
      !ok
    CASE DEFAULT
      CALL finish('mo_nwp_phy_init: init_nwp_phy',  &
        &      'Wrong irad_aero. For Ritter-Geleyn radiation, this irad_aero is not implemented.')
    END SELECT

    
    ssi(:) = ssi_amip(:)
    tsi    = SUM(ssi(:))

    IF ( nh_test_name == 'APE_nh' ) THEN
      ssi(:) = ssi(:)*1365._wp/tsi
      tsi = 1365._wp
    ENDIF
    
    CALL rad_aibi

    zaef(:,:)= 0.0_wp
    
    IF ( irad_aero == 5 ) THEN

      CALL init_aerosol (             &
        & kbdim    = nproma,          & !in
        & pt_patch = p_patch,         & !in
        & aersea   = prm_diag%aersea, & !out
        & aerlan   = prm_diag%aerlan, & !out
        & aerurb   = prm_diag%aerurb, & !out
        & aerdes   = prm_diag%aerdes )
      
    ENDIF
    
  ENDIF
  !------------------------------------------
  !< call for convection
  !------------------------------------------

  IF ( inwp_convection == 1 ) THEN

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init convection')
! Please take care for scale-dependent initializations!

!>reference pressure
    DO jk = nlevp1, 1, -1
        pref(jk)= p0sl * EXP( -vct_a (jk)/h_scal)
    ENDDO

      nsmax = INT(2._wp*pi*re/mean_charlen)

!        WRITE(message_text,'(i3,i10,f20.10)') jg,nsmax,mean_charlen
!       CALL message('nwp_phy_init, nsmax=', TRIM(message_text))
        CALL sucst(54,20020211,0,0)
        CALL su_yoethf
        CALL sucumf(nsmax,nlevp1,pref)
        CALL suphli
        CALL suvdf
        CALL suvdfs
        CALL sucldp
        CALL message('mo_nwp_phy_init:', 'convection initialized')
  ENDIF

  !------------------------------------------
  !< setup for turbulence
  !------------------------------------------

  ! nsfc_type is used for dimensioning local variables in the NWP interface;
  ! thus, it must be set even if no turbulence scheme called
  nsfc_type = 1

  IF ( inwp_turb == 1 ) THEN

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init COSMO turbulence')
    
      rl_start = 1 ! Initialization should be done for all points
      rl_end   = min_rlcell

      i_startblk = p_patch%cells%start_blk(rl_start,1)
      i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

 
!$OMP PARALLEL

!$OMP DO PRIVATE(jb,i_startidx,i_endidx)

  DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)

        CALL init_canopy( ie=nproma, je=1, ke=nlev, ke1=nlevp1, kcm=nlevp1, &
!
         &  istartpar=i_startidx, iendpar=i_endidx, jstartpar=1, jendpar=1, &
!
!         &  hhl=p_metrics%z_ifc(:,:,jb), &
         &  fr_land=ext_data%atm%fr_land(:,jb), plcov=ext_data%atm%plcov_mx(:,jb), & 
         &  sai=prm_diag%sai(:,jb), lai=ext_data%atm%lai_mx(:,jb), &
         &  tai=prm_diag%tai(:,jb), eai=prm_diag%eai(:,jb) )

        CALL organize_turbdiff(action='tran_diff', iini=1, lstfnct=.TRUE., &
!
         &  dt_var=pdtime, dt_tke=pdtime, nprv=1, ntur=1, ntim=1, &
!
         &  ie=nproma, je=1, ke=nlev, ke1=nlevp1, kcm=nlevp1, vst=0, &
!
         &  istart   =i_startidx, iend   =i_endidx, istartu=i_startidx, iendu=i_endidx, &
         &  istartpar=i_startidx, iendpar=i_endidx, istartv=i_startidx, iendv=i_endidx, &
!           
         &  jstart   =1,          jend   =1       , jstartu=1         , jendu=1       , &
         &  jstartpar=1         , jendpar=1       , jstartv=1         , jendv=1       , &
!
         &  l_hori=mean_charlen, hhl=p_metrics%z_ifc(:,:,jb), dp0=p_diag%dpres_mc(:,:,jb), &   
!
         &  fr_land=ext_data%atm%fr_land(:,jb), depth_lk=ext_data%atm%depth_lk(:,jb), &
         &  sai=prm_diag%sai(:,jb), h_ice=prm_diag%h_ice (:,jb), &
!
         &  ps=p_diag%pres_sfc(:,jb), t_g=p_prog_lnd_now%t_g(:,jb), qv_s=p_diag_lnd%qv_s(:,jb), &
!
         &  u=p_diag%u(:,:,jb), v=p_diag%v(:,:,jb), w=p_prog%w(:,:,jb), T=p_diag%temp(:,:,jb), &
         &  qv=p_prog%tracer(:,:,jb,iqv), qc=p_prog%tracer(:,:,jb,iqc), & 
!
         &  prs=p_diag%pres(:,:,jb), rho=p_prog%rho(:,:,jb), epr=p_prog%exner(:,:,jb), &
!
         &  gz0=prm_diag%gz0(:,jb), tcm=prm_diag%tcm(:,jb), tch=prm_diag%tch(:,jb), &
         &  tfm=prm_diag%tfm(:,jb), tfh=prm_diag%tfh(:,jb), tfv=prm_diag%tfv(:,jb), &
!
         &  tke=p_prog_now%tke(:,:,jb), & !  edr=prm_diag%edr(:,:,jb), &
         &  tkvm=prm_diag%tkvm(:,:,jb), tkvh=prm_diag%tkvh (:,:,jb), rcld=prm_diag%rcld(:,:,jb), &
!
         &  u_tens=prm_nwp_tend%ddt_u_turb(:,:,jb), v_tens=prm_nwp_tend%ddt_v_turb(:,:,jb), & 
         &  tketens=prm_nwp_tend%ddt_tke(:,:,jb), &
         &  ut_sso=prm_nwp_tend%ddt_u_sso(:,:,jb), vt_sso=prm_nwp_tend%ddt_v_sso(:,:,jb) ,& 
!
         &  t_2m=prm_diag%t_2m(:,jb), qv_2m=prm_diag%qv_2m(:,jb), td_2m=prm_diag%td_2m (:,jb), &
         &  rh_2m=prm_diag%rh_2m(:,jb), u_10m=prm_diag%u_10m(:,jb), v_10m=prm_diag%v_10m (:,jb), &
         &  shfl_s=prm_diag%shfl_s(:,jb), lhfl_s=prm_diag%lhfl_s(:,jb), &     
!
         &  ierrstat=ierrstat, errormsg=errormsg, eroutine=eroutine )
  ENDDO
!$OMP END DO

!$OMP PARALLEL WORKSHARE
        p_prog %tke (:,:,:) =  p_prog_now%tke (:,:,:)
!$OMP END PARALLEL WORKSHARE

!$OMP END PARALLEL

        CALL message('mo_nwp_phy_init:', 'cosmo turbulence initialized')

    ELSE IF ( inwp_turb == 2) THEN  !ECHAM vdiff

    IF (msg_level >= 12)  CALL message('mo_nwp_phy_init:', 'init ECHAM turbulence')
      ! Currently the tracer indices are sorted such that we count
      ! the water substances first, and then other species like 
      ! aerosols and their precursors. "ntracer" is the total number 
      ! of tracers (including water substances) handled in the model;
      ! "iqt" is the starting index for non-water species.
      ! Before more sophisticated meta-data structure becomes available, 
      ! it is assumed here that all tracers are subject to turbulent mixing.

    ! For surface processes: 
    ! nsfc_type, iwtr, etc. are set in this subroutine. 
    ! See mo_icoham_sfc_indicies.f90 for further details.

      !<KF temporarly set in, has to moved to general place
      CALL init_convect_tables

!      CALL init_sfc_indices( ltestcase, 'APE' ) !call of a hydrostatic testcase
                                                ! to obtain the demanded parameters

      khydromet = 2 !iqt - 1        ! # of hydrometeors
      ktrac = 1   !ntracer - iqt + 1  ! # of non-water species 

     !IF (p_patch%id == 1) CALL init_vdiff_solver( khydromet, ktrac, nproma, nlev, nsfc_type )
      IF (p_patch%id == 1) CALL init_vdiff_solver( khydromet, ktrac, nlev )

      CALL init_vdiff_params( nlev, nlevp1, nlevp1, vct )

      !KF special setting for ICONAM
       tke_min = 1.e-4_wp
        
!$OMP PARALLEL
!$OMP PARALLEL WORKSHARE
        prm_diag% ustar (:,:)   = 1._wp
        prm_diag% kedisp(:,:)   = 0._wp
        prm_diag% thvvar(:,:,:) = 1.e-4_wp
!$OMP END PARALLEL WORKSHARE
!$OMP END PARALLEL

        IF (iwtr<=nsfc_type) prm_diag%z0m_tile(:,iwtr,:) = 1.e-3_wp !see init_surf in echam (or z0m_oce?)
        IF (iice<=nsfc_type) prm_diag%z0m_tile(:,iice,:) = 1.e-3_wp !see init_surf in echam (or z0m_ice?)
        IF (ilnd<=nsfc_type) prm_diag%z0m_tile(:,ilnd,:) = z0m_min ! or maybe a larger value?

!    ENDIF
        
        WRITE(message_text,'(a,3I4)') 'init sfc inidces = ',iwtr,iice,nsfc_type
            CALL message('', TRIM(message_text))

        CALL message('mo_nwp_phy_init:', 'echam turbulence initialized')

! ELSE IF ( inwp_turb == 3) THEN  !DUALM
  ENDIF

  IF ( inwp_surface == 1 ) THEN  ! TERRA



  END IF





END SUBROUTINE init_nwp_phy

!-------------------------------------------------------------------------
!
!
!>
!! Calculates the domain mean of the characteristical
!! length scale of the area. Needed for physics.
!!
!! @par Revision History
!! Implemented by Kristina Froehlich, DWD (2020-10-29).
!!
!!
SUBROUTINE mean_domain_values( p_patch, mean_charlen ) ! output
!
TYPE(t_patch),       INTENT(IN) :: p_patch ! patch on specific level
REAL(wp),            INTENT(OUT)   :: mean_charlen


!mean_charlen (jg) = SQRT (pi*re**2 /REAL(20*nroot**2*4**(p_patch%level),wp))
 mean_charlen      = SQRT (pi*re**2 /REAL(20*nroot**2*4**(p_patch%level),wp))




END SUBROUTINE mean_domain_values



END MODULE mo_nwp_phy_init

