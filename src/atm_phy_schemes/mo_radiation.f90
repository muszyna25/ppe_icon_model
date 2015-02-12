!>
!! @brief Module to provide interface to radiation routines.
!!
!! @remarks
!!   This module contains routines that provide the interface between ECHAM
!!   and the radiation code.  Mostly it organizes and calculates the
!!   information necessary to call the radiative transfer solvers.
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2009-09-19):
!!
!!         Hauke Schmidt, MPI-M, Hamburg (2009-12-18): Few modifications to
!!              allow specific solar irradiance for AMIP-type and preindustrial
!!              simulations.
!!         Luis Kornblueh, MPI-M, Hamburg (2010-04-06): Never ever use write
!!              directly
!!         Martin Schultz, FZJ, Juelich (2010-04-13):
!!              Extracted public parameters into new module mo_radiation_parameters
!!              to avoid circular dependencies in submodels
!!                                      (2010-06-03):
!!              Added submodel calls, decl_sun_cur
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Major segments of this code combines and rewrites (for the ICON standard)
!!   code previously contained in the ECHAM5 routines rad_int.f90,
!!   radiation.f90 and prerad.f90.  Modifications were also made to provide
!!   a cleaner interface to the aerosol and cloud properties. Contributors to
!!   the code from which the present routines were derived include:  M. Jarraud,
!!   ECMWF (1983-06); M.A. Giorgetta, MPI-M (2002-05); U. Schulzweida,  MPI-M
!!   (2002-05); P. Stier MPI-M \& Caltech (2004-04, 2006-07), M. Thomas MPI-M
!!   (2007-06); U. Schlese, MPI-M (2007-06); M. Esch, MPI-M (2007-06); S.J.
!!   Lorenz, MPI-M (2007-11); T. Raddatz, MPI-M (2006-05); I. Kirchner.
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!
MODULE mo_radiation

  USE mo_aerosol_util,         ONLY: zaea_rrtm,zaes_rrtm,zaeg_rrtm
  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish
  USE mo_run_config,           ONLY: ltimer

  USE mo_model_domain,         ONLY: t_patch

  USE mo_math_constants,       ONLY: pi
  USE mo_physical_constants,   ONLY: grav,  rd,    avo,   amd,  amw,  &
    &                                amco2, amch4, amn2o, amo3, amo2, &
    &                                stbo,  vpp_ch4, vpp_n2o

  USE mo_datetime,             ONLY: rdaylen

  USE mo_radiation_config,     ONLY: tsi_radt,   ssi_radt,    &
    &                                irad_co2,   mmr_co2,     &
    &                                irad_ch4,   mmr_ch4,     &
    &                                irad_n2o,   mmr_n2o,     &
    &                                irad_o2,    mmr_o2,      &
    &                                irad_cfc11, vmr_cfc11,   &
    &                                irad_cfc12, vmr_cfc12,   &
    &                                irad_aero,               &
    &                                lrad_aero_diag,          &
    &                                izenith
  USE mo_lnd_nwp_config,       ONLY: isub_seaice, isub_lake

  USE mo_newcld_optics,        ONLY: newcld_optics
  USE mo_bc_aeropt_kinne,      ONLY: set_bc_aeropt_kinne
  USE mo_bc_aeropt_stenchikov, ONLY: add_bc_aeropt_stenchikov

  USE mo_lrtm_par,             ONLY: jpband => nbndlw, jpxsec => maxxsec
  USE mo_lrtm,                 ONLY: lrtm

  USE mo_srtm_config,          ONLY: jpsw, jpinpx
  USE mo_srtm,                 ONLY: srtm_srtm_224gp
  USE mo_get_utc_date_tr,      ONLY: get_utc_date_tr
  USE mo_timer,                ONLY: timer_radiation, timer_start, timer_stop

  USE mo_echam_phy_memory,     ONLY: prm_field
  USE mo_nh_testcases_nml,     ONLY: zenithang

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pre_radiation_nwp, radiation, radiation_nwp, radheat, pre_radiation_nwp_steps


  ! --- radiative transfer parameters
  !
  REAL(wp), PARAMETER :: diff   = 1.66_wp   !< LW Diffusivity Factor

CONTAINS

  SUBROUTINE pre_radiation_nwp_steps( &
    & kbdim,cosmu0_dark,p_inc_rad,p_inc_radheat,p_sim_time,pt_patch,zsmu0,zsct)

    INTEGER, INTENT(IN)   :: &
      & kbdim

    REAL(wp), INTENT(IN)   :: &
      & cosmu0_dark, &
      & p_inc_rad, & !radiation time step in seconds
      & p_inc_radheat, &
      & p_sim_time

    TYPE(t_patch),      INTENT(IN)    :: pt_patch    ! Patch

    REAL(wp), INTENT(OUT), OPTIONAL   :: zsct                  ! solar constant (at time of year)
    REAL(wp), INTENT(OUT)             :: zsmu0(kbdim,pt_patch%nblks_c)   ! Cosine of zenith angle

    REAL(wp) ::                    &
      & p_sim_time_rad,            &
      & zstunde,                   & ! output from routine get_utc_date_tr
      & ztwo, ztho  ,              & 
      & zdek,                      & 
      & zsocof, zeit0,             &
      & zsct_h

    REAL(wp), SAVE ::              &
      & zsct_save, zdtzgl,         &
      & zdeksin,zdekcos


    REAL(wp) ::                          &
      & zsinphi(kbdim,pt_patch%nblks_c) ,&
      & zcosphi(kbdim,pt_patch%nblks_c) ,&
      & zeitrad(kbdim,pt_patch%nblks_c) ,&
      & z_cosmu0(kbdim,pt_patch%nblks_c)

    INTEGER :: jj, itaja   ! output from routine get_utc_date_tr

    INTEGER :: ie,jb,jc,jmu0,n_zsct,nsteps

    INTEGER :: n_cosmu0pos(kbdim,pt_patch%nblks_c)

    INTEGER , SAVE :: itaja_zsct_previous = 0

    ie = kbdim

    IF (izenith == 0) THEN
    ! local insolation = constant = global mean insolation (ca. 340 W/m2)
    ! zenith angle = 0,
      DO jb = 1, pt_patch%nblks_c
        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c
        zsmu0(1:ie,jb) = 1._wp ! sun in zenith everywhere
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt/4._wp ! scale ztsi to get the correct global mean insolation
    ELSEIF(izenith == 1) THEN
    ! circular non-seasonal orbit,
    ! perpetual equinox,
    ! no diurnal cycle,
    ! local time always 12:00
    ! --> sin(time of day)=1 ) and zenith angle depends on latitude only
      DO jb = 1, pt_patch%nblks_c
        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c      
        zsmu0(1:ie,jb) = COS( pt_patch%cells%center(1:ie,jb)%lat )
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt/pi ! because sun is always in local noon, the TSI needs to be
      ! scaled by 1/pi to get the correct global mean insolation
    ELSEIF (izenith == 2) THEN
    ! circular non-seasonal orbit,
    ! perpetual equinox,
    ! no diurnal cycle,
    ! local time always  07:14:15 or 16:45:45
    ! --> sin(time of day)=1/pi and zenith angle depends on latitude only
      DO jb = 1, pt_patch%nblks_c
        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c      
        zsmu0(1:ie,jb) = COS( pt_patch%cells%center(1:ie,jb)%lat )/pi
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt
    ELSEIF (izenith == 3) THEN  !Second: case izenith==3 (time (but no date) needed)
    ! circular non-seasonal orbit,
    ! perpetual equinox,
    ! with diurnal cycle,

      zsmu0(:,:)=0.0_wp
      n_cosmu0pos(:,:) = 0
      nsteps = NINT(p_inc_rad/p_inc_radheat)

      DO jmu0=1,nsteps

        p_sim_time_rad = p_sim_time + (REAL(jmu0,wp)-0.5_wp)*p_inc_radheat

        CALL get_utc_date_tr (                 &
          &   p_sim_time     = p_sim_time_rad, &
          &   itype_calendar = 0,              &
          &   iyear          = jj,             &
          &   nactday        = itaja,          &
          &   acthour        = zstunde )

        ie = kbdim
        DO jb = 1, pt_patch%nblks_c

          IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c

          z_cosmu0(1:ie,jb) = -COS( pt_patch%cells%center(1:ie,jb)%lat ) &
            & *COS( pt_patch%cells%center(1:ie,jb)%lon                &
            &      +zstunde/24._wp* 2._wp*pi )

          DO jc = 1,ie
            IF ( z_cosmu0(jc,jb) > -1.e-5_wp ) THEN
              zsmu0(jc,jb) = zsmu0(jc,jb) + MAX(1.e-3_wp,z_cosmu0(jc,jb))**2
              n_cosmu0pos(jc,jb) = n_cosmu0pos(jc,jb) + 1
            ENDIF
          ENDDO

        ENDDO !jb

      ENDDO!jmu0

      ie = kbdim

      DO jb = 1, pt_patch%nblks_c

        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c

        DO jc = 1,ie
          IF ( n_cosmu0pos(jc,jb) > 0 ) THEN
            zsmu0(jc,jb) = SQRT(zsmu0(jc,jb)/REAL(n_cosmu0pos(jc,jb),wp))
          ELSE
            zsmu0(jc,jb) = cosmu0_dark
          ENDIF
        ENDDO

      ENDDO !jb

      IF (PRESENT(zsct)) zsct = tsi_radt

    ELSEIF (izenith == 4) THEN
    ! elliptical seasonal orbit,
    !  with diurnal cycle

      zsct_h = 0.0_wp
      zsmu0(:,:)=0.0_wp
      n_cosmu0pos(:,:) = 0
      n_zsct = 0

      nsteps = NINT(p_inc_rad/p_inc_radheat)

      DO jmu0=1,nsteps

        p_sim_time_rad = p_sim_time + (REAL(jmu0,wp)-0.5_wp)*p_inc_radheat

        CALL get_utc_date_tr (                 &
          &   p_sim_time     = p_sim_time_rad, &
          &   itype_calendar = 0,              &
          &   iyear          = jj,             &
          &   nactday        = itaja,          &
          &   acthour        = zstunde )

        IF ( itaja /= itaja_zsct_previous ) THEN
          itaja_zsct_previous = itaja

          ztwo    = 0.681_wp + 0.2422_wp*REAL(jj-1949,wp)-REAL((jj-1949)/4,wp)
          ztho    = 2._wp*pi*( REAL(itaja, wp) -1.0_wp + ztwo )/365.2422_wp
          zdtzgl  = 0.000075_wp + 0.001868_wp*COS(      ztho) - 0.032077_wp*SIN(      ztho) &
            - 0.014615_wp*COS(2._wp*ztho) - 0.040849_wp*SIN(2._wp*ztho)
          zdek    = 0.006918_wp - 0.399912_wp*COS(      ztho) + 0.070257_wp*SIN(      ztho) &
            - 0.006758_wp*COS(2._wp*ztho) + 0.000907_wp*SIN(2._wp*ztho) &
            - 0.002697_wp*COS(3._wp*ztho) + 0.001480_wp*SIN(3._wp*ztho)

          zdeksin = SIN (zdek)
          zdekcos = COS (zdek)

          IF ( PRESENT(zsct) ) THEN

            zsocof  = 1.000110_wp + 0.034221_wp*COS(   ztho) + 0.001280_wp*SIN(   ztho) &
              + 0.000719_wp*COS(2._wp*ztho) + 0.000077_wp*SIN(2._wp*ztho)
            zsct_save = zsocof*tsi_radt
            zsct_h = zsct_h + zsct_save
            n_zsct = n_zsct + 1

          ENDIF

        ENDIF

        zeit0   = pi*(zstunde-12._wp)/12._wp + zdtzgl

        ie = kbdim
        DO jb = 1, pt_patch%nblks_c
          IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c

          zsinphi(1:ie,jb)      = SIN (pt_patch%cells%center(1:ie,jb)%lat)
          zcosphi(1:ie,jb)      = SQRT(1.0_wp - zsinphi(1:ie,jb)**2)
          zeitrad(1:ie,jb)      = zeit0 + pt_patch%cells%center(1:ie,jb)%lon
          z_cosmu0(1:ie,jb)     = zdeksin * zsinphi(1:ie,jb) + zdekcos * zcosphi(1:ie,jb) * &
            COS(zeitrad(1:ie,jb))

          DO jc = 1,ie
            IF ( z_cosmu0(jc,jb) > -1.e-5_wp ) THEN
              zsmu0(jc,jb) = zsmu0(jc,jb) + MAX(1.e-3_wp,z_cosmu0(jc,jb))**2
              n_cosmu0pos(jc,jb) = n_cosmu0pos(jc,jb) + 1
            ENDIF
          ENDDO

        ENDDO

      ENDDO !jmu0

      ie = kbdim
      DO jb = 1, pt_patch%nblks_c

        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c

        DO jc = 1,ie
          IF ( n_cosmu0pos(jc,jb) > 0 ) THEN
            ! The averaged cosine of zenith angle is limited to 0.05 in order to avoid 
            ! numerical trouble along the day-night boundary near the model top
            zsmu0(jc,jb) = MAX(0.05_wp,SQRT(zsmu0(jc,jb)/REAL(n_cosmu0pos(jc,jb),wp)))
          ELSE
            zsmu0(jc,jb) = cosmu0_dark
          ENDIF
        ENDDO

      ENDDO !jb

      IF ( PRESENT(zsct) ) THEN
        IF ( n_zsct > 0 ) THEN
          zsct = zsct_h/REAL(n_zsct,wp)
        ELSE
          zsct = zsct_save
        ENDIF
      ENDIF

    ELSEIF (izenith == 5) THEN
     ! Radiative convective equilibrium
     ! circular non-seasonal orbit,
     ! perpetual equinox,
     ! no diurnal cycle,
     ! the product tsi*cos(zenith angle) should equal 340 W/m2
     ! see Popke et al. 2013 and Cronin 2013
      DO jb = 1, pt_patch%nblks_c
        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c
        zsmu0(1:ie,jb) = COS(zenithang*pi/180._wp)
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt ! no rescale tsi was adjstd in atm_phy_nwp w ssi_rce

    ENDIF

  END SUBROUTINE pre_radiation_nwp_steps

  SUBROUTINE pre_radiation_nwp(kbdim,p_inc_rad,p_sim_time,pt_patch,zsmu0,zsct)

    INTEGER, INTENT(IN)   :: &
      & kbdim

    REAL(wp), INTENT(IN)   :: &
      & p_inc_rad, & !radiation time step in seconds
      & p_sim_time

    TYPE(t_patch),      INTENT(IN)    :: pt_patch    ! Patch

    REAL(wp), INTENT(OUT), OPTIONAL   :: zsct                  ! solar constant (at time of year)
    REAL(wp), INTENT(OUT)             :: zsmu0(kbdim,pt_patch%nblks_c)   ! Cosine of zenith angle

    REAL(wp) ::                     &
      & p_sim_time_rad,  &
      & zstunde,                   & ! output from routine get_utc_date_tr
      & ztwo  , ztho  ,            & !
      & zdtzgl, zdek  ,            & !
      & zsocof, zeit0 ,            & !
      & zdeksin,zdekcos!,           & !

    REAL(wp) ::                     &
      & zsinphi(kbdim,pt_patch%nblks_c) ,&
      & zcosphi(kbdim,pt_patch%nblks_c) ,&
      & zeitrad(kbdim,pt_patch%nblks_c)! ,&

    INTEGER :: &
      & jj, itaja, jb, ie !& ! output from routine get_utc_date_tr

    INTEGER , SAVE :: itaja_zsct_previous = 0
    REAL(wp), SAVE :: zsct_save

    ie = kbdim

    !First: cases izenith==0 to izenith==2 (no date and time needed)
    IF (izenith == 0) THEN
     ! for testing: provisional setting of cos(zenith angle) and TSI
     ! The global mean insolation is TSI/4 (ca. 340 W/m2)
      DO jb = 1, pt_patch%nblks_c
        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c
        zsmu0(1:ie,jb) = 1._wp ! sun in zenith everywhere
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt/4._wp ! scale ztsi to get the correct global mean insolation
      RETURN
    ELSEIF(izenith == 1) THEN
      ! circular non-seasonal orbit, zenith angle dependent on latitude only,
      ! no diurnal cycle (always at 12:00 local time --> sin(time of day)=1 )
      DO jb = 1, pt_patch%nblks_c
        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c      
        zsmu0(1:ie,jb) = COS( pt_patch%cells%center(1:ie,jb)%lat )
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt/pi ! because sun is always in local noon, the TSI needs to be
                                       ! scaled by 1/pi to get the correct global mean insolation
      RETURN
    ELSEIF (izenith == 2) THEN
      ! circular non-seasonal orbit, no diurnal cycle
      ! at 07:14:15 or 16:45:45 local time (--> sin(time of day)=1/pi )
      DO jb = 1, pt_patch%nblks_c
        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c      
        zsmu0(1:ie,jb) = COS( pt_patch%cells%center(1:ie,jb)%lat )/pi
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt
      RETURN
    ENDIF

    p_sim_time_rad = p_sim_time + 0.5_wp*p_inc_rad

    CALL get_utc_date_tr (                 &
      &   p_sim_time     = p_sim_time_rad, &
      &   itype_calendar = 0,              &
      &   iyear          = jj,             &
      &   nactday        = itaja,          &
      &   acthour        = zstunde )

    !Second case izenith==3 (time (but no date) needed)
    IF (izenith == 3) THEN
      
      DO jb = 1, pt_patch%nblks_c
        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c
        zsmu0(1:ie,jb) = -COS( pt_patch%cells%center(1:ie,jb)%lat ) &
          & *COS( pt_patch%cells%center(1:ie,jb)%lon                &
          &      +zstunde/24._wp* 2._wp*pi )
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt

    !Third: case izenith=4 (time and date needed)
    ELSEIF (izenith == 4) THEN
    
      ztwo    = 0.681_wp + 0.2422_wp*REAL(jj-1949,wp)-REAL((jj-1949)/4,wp)
      ztho    = 2._wp*pi*( REAL(itaja, wp) -1.0_wp + ztwo )/365.2422_wp
      zdtzgl  = 0.000075_wp + 0.001868_wp*COS(      ztho) - 0.032077_wp*SIN(      ztho) &
        - 0.014615_wp*COS(2._wp*ztho) - 0.040849_wp*SIN(2._wp*ztho)
      zdek    = 0.006918_wp - 0.399912_wp*COS(      ztho) + 0.070257_wp*SIN(      ztho) &
        - 0.006758_wp*COS(2._wp*ztho) + 0.000907_wp*SIN(2._wp*ztho) &
        - 0.002697_wp*COS(3._wp*ztho) + 0.001480_wp*SIN(3._wp*ztho)
      zeit0   = pi*(zstunde-12._wp)/12._wp + zdtzgl
      zdeksin = SIN (zdek)
      zdekcos = COS (zdek)

      IF ( PRESENT(zsct) ) THEN
        !decide whether new zsct calculation is necessary
        IF ( itaja /= itaja_zsct_previous ) THEN
          itaja_zsct_previous = itaja
          zsocof  = 1.000110_wp + 0.034221_wp*COS(   ztho) + 0.001280_wp*SIN(   ztho) &
            + 0.000719_wp*COS(2._wp*ztho) + 0.000077_wp*SIN(2._wp*ztho)
          zsct_save = zsocof*tsi_radt
        ENDIF
        zsct = zsct_save
      ENDIF

      DO jb = 1, pt_patch%nblks_c
        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c
        zsinphi(1:ie,jb)      = SIN (pt_patch%cells%center(1:ie,jb)%lat)
        zcosphi(1:ie,jb)      = SQRT(1.0_wp - zsinphi(1:ie,jb)**2)
        zeitrad(1:ie,jb)      = zeit0 + pt_patch%cells%center(1:ie,jb)%lon
        zsmu0(1:ie,jb)        = zdeksin * zsinphi(1:ie,jb) + zdekcos * zcosphi(1:ie,jb) * &
          COS(zeitrad(1:ie,jb))

      ENDDO
      
    ELSEIF (izenith == 5) THEN
     ! Radiative convective equilibrium
     ! circular non-seasonal orbit,
     ! perpetual equinox,
     ! no diurnal cycle,
     ! the product tsi*cos(zenith angle) should equal 340 W/m2
     ! see Popke et al. 2013 and Cronin 2013
      DO jb = 1, pt_patch%nblks_c
        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c
        zsmu0(1:ie,jb) = COS(zenithang*pi/180._wp)
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi_radt ! no rescale tsi was adjstd in atm_phy_nwp w ssi_rce

    ENDIF

  END SUBROUTINE pre_radiation_nwp

  !-----------------------------------------------------------------------------
  !>
  !! @brief Organizes the calls to the ratiation solver
  !!
  !! @remarks This routine organises the input/output for the radiation
  !! computation.  The state of radiatively active constituents is set as the
  !! input. Output are flux transmissivities (ratio solar flux/solar input) 
  !! and thermal fluxes at all the half levels of the grid. This output will be 
  !! used in radheat at all time steps until the next full radiation time step.
  !
  SUBROUTINE radiation(                                                    &
    ! input
    &  jg, jb                                                              &
    & ,jce               ,kbdim           ,klev             ,klevp1        &
    & ,ktype             ,zland           ,zglac            ,cos_mu0       &
    & ,alb_vis_dir       ,alb_nir_dir     ,alb_vis_dif      ,alb_nir_dif   &
    & ,emis_rad                                                            &
    & ,tk_sfc            ,z_mc            ,pp_hl            ,pp_fl         &
    & ,tk_fl             ,qm_vap          ,qm_liq           ,qm_ice        &
    & ,qm_o3                                                               &
    & ,cdnc              ,cld_frc                                          &
    & ,zaeq1, zaeq2, zaeq3, zaeq4, zaeq5 , dt_rad                          &
    ! output
    & ,cld_cvr                                                             &
    & ,flx_lw_net_clr    ,trm_sw_net_clr  ,flx_lw_net       ,trm_sw_net    &
    & ,opt_halo_cosmu0  )

    ! input
    ! -----
    !
    INTEGER, INTENT(in)   :: &
      &  jg,                 & !< domain index
      &  jb,                 & !< block index
      &  jce,                & !< end   index for loop over block
      &  kbdim,              & !< dimension of block over cells
      &  klev,               & !< number of full levels = number of layers
      &  klevp1                !< number of half levels = number of layer interfaces
    INTEGER, INTENT(in)   :: &
      &  ktype(kbdim)          !< type of convection

    REAL(wp), INTENT(in)  :: &
      &  zland(kbdim),       & !< land-sea mask. (1. = land, 0. = sea/lakes)
      &  zglac(kbdim),       & !< fraction of land covered by glaciers
      &  cos_mu0(kbdim),     & !< cos of zenith angle
      &  alb_vis_dir(kbdim), & !< surface albedo for visible range and direct light
      &  alb_nir_dir(kbdim), & !< surface albedo for NIR range and direct light
      &  alb_vis_dif(kbdim), & !< surface albedo for visible range and diffuse light
      &  alb_nir_dif(kbdim), & !< surface albedo for NIR range and diffuse light
      &  emis_rad(kbdim),    & !< longwave surface emissivity
      &  tk_sfc(kbdim),      & !< Surface temperature
      &  z_mc(kbdim,klev),   & !< height at full levels [m]
      &  pp_hl(kbdim,klevp1),& !< pressure at half levels [Pa]
      &  pp_fl(kbdim,klev),  & !< Pressure at full levels [Pa]
      &  tk_fl(kbdim,klev),  & !< Temperature on full levels [K]
      &  qm_vap(kbdim,klev), & !< Water vapor mixing ratio
      &  qm_liq(kbdim,klev), & !< Liquid water mixing ratio
      &  qm_ice(kbdim,klev), & !< Ice water mixing ratio
      &  cdnc(kbdim,klev),   & !< Cloud drop number concentration
      &  cld_frc(kbdim,klev),& !< Cloud fraction
      &  zaeq1(kbdim,klev) , & !< aerosol continental
      &  zaeq2(kbdim,klev) , & !< aerosol maritime
      &  zaeq3(kbdim,klev) , & !< aerosol urban
      &  zaeq4(kbdim,klev) , & !< aerosol volcano ashes
      &  zaeq5(kbdim,klev) , & !< aerosol stratospheric background
      &  dt_rad                !< radiation time step
    
    LOGICAL, INTENT(in), OPTIONAL :: opt_halo_cosmu0    

    ! output
    ! ------
    !
    REAL(wp), INTENT(out) ::           &
      &  cld_cvr(kbdim),               & !< Cloud cover in a column
      &  flx_lw_net_clr(kbdim,klevp1), & !< Net clear-sky longwave radiative flux [W/m**2] (positive downward)
      &  trm_sw_net_clr(kbdim,klevp1), & !< Net clear-sky solar transmissivity (= net clear-sky shortwave radiative flux normalized by irradiance)
      &  flx_lw_net(kbdim,klevp1),     & !< Net longwave radiative flux [W/m**2] (positive downward)
      &  trm_sw_net(kbdim,klevp1)        !< Net solar transmissivity (= net shortwave radiative flux normalized by irradiance)


    INTEGER  :: jk, jl

    REAL(wp) ::                     &
    &    cos_mu0_halo,              & !< cos(zenith angle) value delimiting the halo
    &    cos_mu0_mod(kbdim)           !< modified cos(zenith angle)

    REAL(wp) ::                       &
      &  pp_sfc(kbdim),               & !< surface pressure [Pa}
      &  tk_hl(kbdim,klevp1),         & !< Tempeature at half levels [Pa]
      &  xq_vap(kbdim,klev),          & !< Water vapor mixing ratio
      &  xq_liq(kbdim,klev),          & !< Liquid water mixing ratio
      &  xq_ice(kbdim,klev),          & !< Ice mixing ratio
      &  cld_frc_sec(kbdim,klev),     & !< secure cloud fraction [m2/m2]
      &  xm_co2(kbdim,klev),          & !< CO2 mixing ratio
      &  qm_o3(kbdim,klev),           & !< Ozone mixing ratio
      &  xm_o2(kbdim,klev),           & !< O2 mixing ratio
      &  xm_ch4(kbdim,klev),          & !< Methane mixing ratio
      &  xm_n2o(kbdim,klev),          & !< Nitrous Oxide mixing ratio
      &  xm_cfc11(kbdim,klev),        & !< CFC 11 mixing ratio
      &  xm_cfc12(kbdim,klev),        & !< CFC 12 mixing ratio
      &  flx_uplw_sfc(kbdim),         & !< Srfc upward lw flux  [Wm2]
      &  flx_upsw_sfc(kbdim),         & !< Srfc upward sw flux  [Wm2]
      &  flx_uplw_clr_sfc(kbdim),     & !< Srfc upward lw flux (clear sky) [Wm2]
      &  flx_upsw_clr_sfc(kbdim),     & !< Srfc upward sw flux (clear sky) [Wm2]
      &  flx_sw_net(kbdim,klevp1),    & !< Net SW flux [Wm2] (positive down)
      &  flx_sw_net_clr(kbdim,klevp1)   !< Net SW flux (clear sky) [Wm2] (positive down)

    LOGICAL :: l_halo_cosmu0

    IF (ltimer) CALL timer_start(timer_radiation)

    ! check for optional arguments
    IF ( PRESENT(opt_halo_cosmu0) ) THEN
      l_halo_cosmu0 = opt_halo_cosmu0
    ELSE
      l_halo_cosmu0 = .TRUE.
    ENDIF

    IF (l_halo_cosmu0) THEN
      !
      ! 1.0 Add halo to sun-lit area
      ! ----------------------------
      !
      ! --- Add a halo to the sun-lit hemisphere in order to include all points,
      !     which are sun-lit at any time step, at which the heating rate is
      !     calculated using this radiative transfer calculation.
      !
      !     The width of the halo is set to the change in cos(zenith angle) over
      !     half of the radiation time step dt_rad at the equator at equinox.
      !
      cos_mu0_halo = -SIN(pi*dt_rad/rdaylen)
      !
      WHERE (cos_mu0(1:jce) > cos_mu0_halo)
        ! Within the sun-lit hemisphere and the halo: use a minimum value of
        ! 0.1 for the SW computations.
        cos_mu0_mod(1:jce) = MAX(cos_mu0(1:jce),0.1_wp)
      ELSEWHERE
        ! Elsewhere keep the negative cos(zenith angle). No SW computations
        ! will be made in this area.
        cos_mu0_mod(1:jce) = cos_mu0(1:jce)
      END WHERE
    ELSE
      cos_mu0_mod(1:jce) = cos_mu0(1:jce)
    ENDIF !l_halo_cosmu0
    
    !
    ! 1.1 p, T, q(vap,liq,ice) and clouds
    ! -----------------------------------
    !
    ! --- Pressure (surface and distance between half levels)
    !
    pp_sfc(1:jce)        = pp_hl(1:jce,klevp1)
    !
    ! --- temperature at half levels
    !
    DO jk=2,klev
      DO jl = 1, jce
        tk_hl(jl,jk) = (tk_fl(jl,jk-1)*pp_fl(jl,jk-1)*( pp_fl(jl,jk)          &
          &    - pp_hl(jl,jk) ) + tk_fl(jl,jk)*pp_fl(jl,jk)*( pp_hl(jl,jk)    &
          &    - pp_fl(jl,jk-1))) /(pp_hl(jl,jk)*(pp_fl(jl,jk) -pp_fl(jl,jk-1)))
      END DO
    END DO
    DO jl = 1, jce
      tk_hl(jl,klevp1) = tk_sfc(jl)
      tk_hl(jl,1)      = tk_fl(jl,1)-pp_fl(jl,1)*(tk_fl(jl,1) - tk_hl(jl,2))  &
        &                / (pp_fl(jl,1)-pp_hl(jl,2))
    END DO
    !
    ! --- phases of water substance
    !
    xq_vap(1:jce,:) = MAX(qm_vap(1:jce,:),0.0_wp)
    xq_liq(1:jce,:) = MAX(qm_liq(1:jce,:),0.0_wp)       ! cloud liquid
    xq_ice(1:jce,:) = MAX(qm_ice(1:jce,:),0.0_wp)       ! cloud ice
    !
    ! --- cloud cover
    !
    DO jk = 1, klev
      DO jl = 1, jce

        ! no cloud water -> no cloud fraction
        IF (xq_liq(jl,jk) > 0.0_wp .OR. xq_ice(jl,jk) > 0.0_wp) THEN
          cld_frc_sec(jl,jk) = cld_frc(jl,jk)
        ELSE
          cld_frc_sec(jl,jk) = 0.0_wp
        END IF

        ! cloud fraction <= 100% !
        cld_frc_sec(jl,jk) = MIN( cld_frc_sec(jl,jk), 1.0_wp )

      END DO
    END DO

    cld_cvr(:) = 0._wp
    cld_cvr(1:jce) = 1.0_wp - cld_frc_sec(1:jce,1)
    DO jk = 2, klev
      cld_cvr(1:jce) = cld_cvr(1:jce)                                                &
        &              *(1.0_wp-MAX(cld_frc_sec(1:jce,jk),cld_frc_sec(1:jce,jk-1)))  &
        &              /(1.0_wp-MIN(cld_frc_sec(1:jce,jk-1),1.0_wp-EPSILON(1.0_wp)))
    END DO
    cld_cvr(1:jce) = 1.0_wp-cld_cvr(1:jce)
    !
    ! 1.2 Non-water tracers
    ! ---------------------
    !
    ! --- gas profiles in [ppm]
    !
    xm_co2   (1:jce,:) = gas_profile(jce, klev, irad_co2,    &
      &                              mmr_gas = mmr_co2   )
    xm_ch4   (1:jce,:) = gas_profile(jce, klev, irad_ch4,    &
      &                              mmr_gas = mmr_ch4,      &
      &                              pressure = pp_fl,       &
      &                              xp = vpp_ch4        )
    xm_n2o   (1:jce,:) = gas_profile(jce, klev, irad_n2o,    &
      &                              mmr_gas = mmr_n2o,      &
      &                              pressure = pp_fl,       &
      &                              xp = vpp_n2o        )
    xm_o2    (1:jce,:) = gas_profile(jce, klev, irad_o2,     &
      &                              mmr_gas = mmr_o2    )
#ifdef __SX__
    xm_cfc11 (1:jce,:) = gas_profile(jce, klev, irad_cfc11,  &
      &                              mmr_gas = REAL(vmr_cfc11,wp) )
    xm_cfc12 (1:jce,:) = gas_profile(jce, klev, irad_cfc12,  &
      &                              mmr_gas = REAL(vmr_cfc12,wp) )
#else
    xm_cfc11 (1:jce,:) = gas_profile(jce, klev, irad_cfc11,  &
      &                              mmr_gas = vmr_cfc11 )
    xm_cfc12 (1:jce,:) = gas_profile(jce, klev, irad_cfc12,  &
      &                              mmr_gas = vmr_cfc12 )
#endif
    !

    ! 2.0 Call interface to radiation solver
    ! --------------------------------------
    !
    CALL rrtm_interface(                                                    &
      ! input
      & jg              ,jb                                                ,&
      & jce             ,kbdim           ,klev                             ,&
      & ktype           ,zland           ,zglac                            ,&
      & cos_mu0_mod                                                        ,&
      & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
      & emis_rad        ,z_mc                                              ,&
      & pp_fl           ,pp_hl           ,pp_sfc          ,tk_fl           ,&
      & tk_hl           ,tk_sfc          ,xq_vap                           ,&
      & xq_liq          ,xq_ice                                            ,&
      & cdnc                                                               ,&
      & cld_frc_sec                                                        ,&
      & qm_o3           ,xm_co2          ,xm_ch4                           ,&
      & xm_n2o          ,xm_cfc11        ,xm_cfc12        ,xm_o2           ,&
      & zaeq1,zaeq2,zaeq3,zaeq4,zaeq5                                      ,&
      ! output
      & flx_lw_net      ,flx_sw_net      ,flx_lw_net_clr  ,flx_sw_net_clr  ,&
      & flx_uplw_sfc    ,flx_upsw_sfc    ,flx_uplw_clr_sfc,flx_upsw_clr_sfc )

    !
    ! 3.0 Additional diagnostics
    ! ---------------
    !
    ! --- Net solar transmissivity
    trm_sw_net (:,:) = 0._wp
    trm_sw_net (1:jce,1:klevp1) = flx_sw_net (1:jce,1:klevp1)                       &
      &                          / SPREAD(cos_mu0_mod(1:jce)*tsi_radt,2,klevp1)
    !
    ! --- Net clear-sky solar transmissivity
    trm_sw_net_clr(:,:) = 0._wp
    trm_sw_net_clr(1:jce,1:klevp1) = flx_sw_net_clr(1:jce,1:klevp1)                 &
      &                         / SPREAD(cos_mu0_mod(1:jce)*tsi_radt,2,klevp1)

    IF (ltimer) CALL timer_stop(timer_radiation)

  END SUBROUTINE radiation

  !-----------------------------------------------------------------------------
  !>
  !! @brief Organizes the calls to the ratiation solver
  !!
  !! @remarks This routine organises the input/output for the radiation
  !! computation.  The state of radiatively active constituents is set as the
  !! input. Output are flux transmissivities and emissivities at all the half
  !! levels of the grid (respectively ratio solar flux/solar input and ratio
  !! thermal flux/local black-body flux). This output will be used in radheat
  !! at all time steps until the next full radiation time step.
  !
  SUBROUTINE radiation_nwp(                                                &
    ! input
    &  jg, jb                                                              &
    & ,jce               ,kbdim           ,klev             ,klevp1        &
    & ,ktype             ,zland           ,zglac            ,cos_mu0       &
    & ,alb_vis_dir       ,alb_nir_dir     ,alb_vis_dif      ,alb_nir_dif   &
    & ,emis_rad                                                            &
    & ,tk_sfc            ,z_mc            ,pp_hl            ,pp_fl         &
    & ,tk_fl             ,qm_vap          ,qm_liq           ,qm_ice        &
    & ,qm_o3                                                               &
    & ,cdnc              ,cld_frc                                          &
    & ,zaeq1, zaeq2, zaeq3, zaeq4, zaeq5 , dt_rad                          &
    ! output
    & ,cld_cvr, flx_lw_net, flx_uplw_sfc, trsol_net, trsol_up_toa,         &
    &  trsol_up_sfc, trsol_dn_sfc_diffus, trsol_clr_sfc, trsol_par_sfc     )

    ! input
    ! -----
    !
    INTEGER, INTENT(in)   :: &
      &  jg,                 & !< domain index
      &  jb,                 & !< block index
      &  jce,                & !< end   index for loop over block
      &  kbdim,              & !< dimension of block over cells
      &  klev,               & !< number of full levels = number of layers
      &  klevp1                !< number of half levels = number of layer interfaces
    INTEGER, INTENT(in)   :: &
      &  ktype(kbdim)          !< type of convection

    REAL(wp), INTENT(in)  :: &
      &  zland(kbdim),       & !< land-sea mask. (1. = land, 0. = sea/lakes)
      &  zglac(kbdim),       & !< fraction of land covered by glaciers
      &  cos_mu0(kbdim),     & !< cos of zenith angle
      &  alb_vis_dir(kbdim), & !< surface albedo for visible range and direct light
      &  alb_nir_dir(kbdim), & !< surface albedo for NIR range and direct light
      &  alb_vis_dif(kbdim), & !< surface albedo for visible range and diffuse light
      &  alb_nir_dif(kbdim), & !< surface albedo for NIR range and diffuse light
      &  emis_rad(kbdim),    & !< longwave surface emissivity
      &  tk_sfc(kbdim),      & !< Surface temperature
      &  z_mc(kbdim,klev),   & !< height at full levels [m]
      &  pp_hl(kbdim,klevp1),& !< pressure at half levels [Pa]
      &  pp_fl(kbdim,klev),  & !< Pressure at full levels [Pa]
      &  tk_fl(kbdim,klev),  & !< Temperature on full levels [K]
      &  qm_vap(kbdim,klev), & !< Water vapor mixing ratio
      &  qm_liq(kbdim,klev), & !< Liquid water mixing ratio
      &  qm_ice(kbdim,klev), & !< Ice water mixing ratio
      &  cdnc(kbdim,klev),   & !< Cloud drop number concentration
      &  cld_frc(kbdim,klev),& !< Cloud fraction
      &  zaeq1(kbdim,klev) , & !< aerosol continental
      &  zaeq2(kbdim,klev) , & !< aerosol maritime
      &  zaeq3(kbdim,klev) , & !< aerosol urban
      &  zaeq4(kbdim,klev) , & !< aerosol volcano ashes
      &  zaeq5(kbdim,klev) , & !< aerosol stratospheric background
      &  dt_rad                !< radiation time step


    ! output
    ! ------
    !
    REAL(wp), INTENT(inout) ::      &
      &  cld_cvr(kbdim),            & !< Cloud cover in a column
      &  flx_lw_net(kbdim,klevp1),  & !< Net longwave radiative flux [W/m**2]
      &  flx_uplw_sfc(kbdim),       & !< Srfc upward lw flux  [Wm2]
      &  trsol_net (kbdim,klevp1),  & !< Net solar transmissivity (=net radiative flux normalized by irradiance)
      &  trsol_up_toa(kbdim),       & !< TOA upward shortwave transmissivity (normalized upward flux)
      &  trsol_up_sfc(kbdim),       & !< Surface upward shortwave transmissivity (normalized upward flux)
      &  trsol_dn_sfc_diffus(kbdim),& !< Surface downward diffuse shortwave transmissivity (normalized diffuse downward flux)
      &  trsol_clr_sfc(kbdim),      & !< Surface net clear-sky solar transmissivity
      &  trsol_par_sfc(kbdim)         !< Surface transmissivity for photosynthetically active part of solar radiation

    INTEGER  :: jk, jl

    REAL(wp) ::                     &
      &  pp_sfc(kbdim),             & !< surface pressure [Pa}
      &  tk_hl(kbdim,klevp1),       & !< Tempeature at half levels [Pa]
      &  xq_vap(kbdim,klev),        & !< Water vapor mixing ratio
      &  xq_liq(kbdim,klev),        & !< Liquid water mixing ratio
      &  xq_ice(kbdim,klev),        & !< Ice mixing ratio
      &  cld_frc_sec(kbdim,klev),   & !< secure cloud fraction [m2/m2]
      &  xm_co2(kbdim,klev),        & !< CO2 mixing ratio
      &  qm_o3(kbdim,klev),         & !< Ozone mixing ratio
      &  xm_o2(kbdim,klev),         & !< O2 mixing ratio
      &  xm_ch4(kbdim,klev),        & !< Methane mixing ratio
      &  xm_n2o(kbdim,klev),        & !< Nitrous Oxide mixing ratio
      &  xm_cfc11(kbdim,klev),      & !< CFC 11 mixing ratio
      &  xm_cfc12(kbdim,klev),      & !< CFC 12 mixing ratio
      &  flx_uplw_sfc_clr(kbdim),   & !< Srfc upward lw flux (clear sky) [Wm2]
      &  flx_upsw_sfc_clr(kbdim),   & !< Srfc upward sw flux (clear sky) [Wm2]
      &  flx_upsw_sfc(kbdim),       & !< Srfc upward sw flux (total) [Wm2]
      &  flx_upsw_toa(kbdim),       & !< TOA  upward sw flux (total) [Wm2]
      &  flx_dnsw_diff_sfc(kbdim),  & !< Srfc diffuse downward sw flux (total) [Wm2]
      &  flx_par_sfc(kbdim),        & !< Srfc downward PAR flux [Wm2]
      &  flx_sw_net(kbdim,klevp1),  & !< Net dwnwrd SW flux [Wm2]
      &  flx_lw_net_clr(kbdim,klevp1),& !< Net dn LW flux (clear sky) [Wm2]
      &  flx_sw_net_clr(kbdim,klevp1)   !< Net dn SW flux (clear sky) [Wm2]


    IF (ltimer) CALL timer_start(timer_radiation)
    
    !
    ! 1.1 p, T, q(vap,liq,ice) and clouds
    ! -----------------------------------
    !
    ! --- Pressure (surface and distance between half levels)
    !
    pp_sfc(1:jce)        = pp_hl(1:jce,klevp1)
    !
    ! --- temperature at half levels
    !
    DO jk=2,klev
      DO jl = 1, jce
        tk_hl(jl,jk) = (tk_fl(jl,jk-1)*pp_fl(jl,jk-1)*( pp_fl(jl,jk)          &
          &    - pp_hl(jl,jk) ) + tk_fl(jl,jk)*pp_fl(jl,jk)*( pp_hl(jl,jk)    &
          &    - pp_fl(jl,jk-1))) /(pp_hl(jl,jk)*(pp_fl(jl,jk) -pp_fl(jl,jk-1)))
      END DO
    END DO
    DO jl = 1, jce
      tk_hl(jl,klevp1) = tk_sfc(jl)
      tk_hl(jl,1)      = tk_fl(jl,1)-pp_fl(jl,1)*(tk_fl(jl,1) - tk_hl(jl,2))  &
        &                / (pp_fl(jl,1)-pp_hl(jl,2))
    END DO
    !
    ! --- phases of water substance
    !
    xq_vap(1:jce,:) = MAX(qm_vap(1:jce,:),0.0_wp)
    xq_liq(1:jce,:) = MAX(qm_liq(1:jce,:),0.0_wp)       ! cloud liquid
    xq_ice(1:jce,:) = MAX(qm_ice(1:jce,:),0.0_wp)       ! cloud ice
    !
    ! --- cloud cover
    !
    DO jk = 1, klev
      DO jl = 1, jce

        ! no cloud water -> no cloud fraction
        IF (xq_liq(jl,jk) > 0.0_wp .OR. xq_ice(jl,jk) > 0.0_wp) THEN
          cld_frc_sec(jl,jk) = cld_frc(jl,jk)
        ELSE
          cld_frc_sec(jl,jk) = 0.0_wp
        END IF

        ! cloud fraction <= 100% !
        cld_frc_sec(jl,jk) = MIN( cld_frc_sec(jl,jk), 1.0_wp )

      END DO
    END DO

    cld_cvr(1:jce) = 1.0_wp - cld_frc_sec(1:jce,1)
    DO jk = 2, klev
      cld_cvr(1:jce) = cld_cvr(1:jce)                                                &
        &              *(1.0_wp-MAX(cld_frc_sec(1:jce,jk),cld_frc_sec(1:jce,jk-1)))  &
        &              /(1.0_wp-MIN(cld_frc_sec(1:jce,jk-1),1.0_wp-EPSILON(1.0_wp)))
    END DO
    cld_cvr(1:jce) = 1.0_wp-cld_cvr(1:jce)
    !
    ! 1.2 Non-water tracers
    ! ---------------------
    !
    ! --- gas profiles in [ppm]
    !
    xm_co2   (1:jce,:) = gas_profile(jce, klev, irad_co2,    &
      &                              mmr_gas = mmr_co2   )
    xm_ch4   (1:jce,:) = gas_profile(jce, klev, irad_ch4,    &
      &                              mmr_gas = mmr_ch4,      &
      &                              pressure = pp_fl,       &
      &                              xp = vpp_ch4        )
    xm_n2o   (1:jce,:) = gas_profile(jce, klev, irad_n2o,    &
      &                              mmr_gas = mmr_n2o,      &
      &                              pressure = pp_fl,       &
      &                              xp = vpp_n2o        )
    xm_o2    (1:jce,:) = gas_profile(jce, klev, irad_o2,     &
      &                              mmr_gas = mmr_o2    )
#ifdef __SX__
    xm_cfc11 (1:jce,:) = gas_profile(jce, klev, irad_cfc11,  &
      &                              mmr_gas = REAL(vmr_cfc11,wp) )
    xm_cfc12 (1:jce,:) = gas_profile(jce, klev, irad_cfc12,  &
      &                              mmr_gas = REAL(vmr_cfc12,wp) )
#else
    xm_cfc11 (1:jce,:) = gas_profile(jce, klev, irad_cfc11,  &
      &                              mmr_gas = vmr_cfc11 )
    xm_cfc12 (1:jce,:) = gas_profile(jce, klev, irad_cfc12,  &
      &                              mmr_gas = vmr_cfc12 )
#endif
    !

    ! 2.0 Call interface to radiation solver
    ! --------------------------------------
    !
    CALL rrtm_interface(                                                    &
      ! input
      & jg              ,jb                                                ,&
      & jce             ,kbdim           ,klev                             ,&
      & ktype           ,zland           ,zglac                            ,&
      & cos_mu0                                                            ,&
      & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
      & emis_rad        ,z_mc                                              ,&
      & pp_fl           ,pp_hl           ,pp_sfc          ,tk_fl           ,&
      & tk_hl           ,tk_sfc          ,xq_vap                           ,&
      & xq_liq          ,xq_ice                                            ,&
      & cdnc                                                               ,&
      & cld_frc_sec                                                        ,&
      & qm_o3           ,xm_co2          ,xm_ch4                           ,&
      & xm_n2o          ,xm_cfc11        ,xm_cfc12        ,xm_o2           ,&
      & zaeq1,zaeq2,zaeq3,zaeq4,zaeq5                                      ,&
      ! output
      & flx_lw_net      ,flx_sw_net      ,flx_lw_net_clr  ,flx_sw_net_clr  ,&
      & flx_uplw_sfc    ,flx_upsw_sfc    ,flx_uplw_sfc_clr,flx_upsw_sfc_clr,&
      & flx_dnsw_diff_sfc, flx_upsw_toa  ,flx_par_sfc                       )


    !
    ! 3.0 Diagnostics
    ! ---------------
    !
    ! --- Diagnose transmissivities from solar fluxes

    trsol_net (1:jce,1:klevp1) = flx_sw_net (1:jce,1:klevp1)/SPREAD(cos_mu0(1:jce)*tsi_radt,2,klevp1)
    trsol_clr_sfc(1:jce)       = flx_sw_net_clr(1:jce,klevp1)/(cos_mu0(1:jce)*tsi_radt)
    trsol_up_sfc(1:jce)        = flx_upsw_sfc  (1:jce)/(cos_mu0(1:jce)*tsi_radt)
    trsol_up_toa(1:jce)        = flx_upsw_toa  (1:jce)/(cos_mu0(1:jce)*tsi_radt)
    trsol_dn_sfc_diffus(1:jce) = flx_dnsw_diff_sfc(1:jce)/(cos_mu0(1:jce)*tsi_radt)
    trsol_par_sfc(1:jce)       = flx_par_sfc(1:jce)/(cos_mu0(1:jce)*tsi_radt)

    IF (ltimer) CALL timer_stop(timer_radiation)

  END SUBROUTINE radiation_nwp

  !---------------------------------------------------------------------------
  !>
  !! GAS_PROFILE:  Determines Gas distributions based on case specification
  !!
  !! @par Revsision History
  !! B. Stevens (2009-08).
  !!
  !! Description: This routine calculates the gas distributions for one of
  !! five cases:  (0) no gas present; (1) prognostic gas; (2) specified
  !! mixing ratio; (3) mixing ratio decaying with height given profile;
  !! (4) scenario run with different mixing ratio.
  !
  FUNCTION gas_profile (jce, klev, igas, mmr_gas, gas_scenario, mmr_gas_v, &
    &                   gas_scenario_v, gas_val, xp, pressure)

    INTEGER,  INTENT (in)           :: jce, klev, igas
    REAL(wp), INTENT (in), OPTIONAL :: mmr_gas, gas_scenario
    REAL(wp), INTENT (in), OPTIONAL :: pressure(:,:), xp(3)
    REAL(wp), INTENT (in), OPTIONAL :: mmr_gas_v(:,:)
    REAL(wp), INTENT (in), OPTIONAL :: gas_scenario_v(:,:)
    REAL(wp), INTENT (in), OPTIONAL :: gas_val(:,:)

    REAL(wp) :: gas_profile(jce,klev)
    REAL(wp) :: zx_d, zx_m
    LOGICAL  :: gas_initialized

    gas_initialized = .FALSE.
    SELECT CASE (igas)
    CASE (0)
      gas_profile(1:jce,:) = 0.0_wp
      gas_initialized = .TRUE.
    CASE (1)
      IF (PRESENT(gas_val)) THEN
        gas_profile(1:jce,:) = MAX(gas_val(1:jce,:), EPSILON(1.0_wp))
        gas_initialized = .TRUE.
      END IF
    CASE (2)
      IF (PRESENT(mmr_gas)) THEN
        gas_profile(1:jce,:) = mmr_gas
        gas_initialized = .TRUE.
      ELSE IF (PRESENT(mmr_gas_v)) THEN
        gas_profile(1:jce,:) = mmr_gas_v(1:jce,:)
        gas_initialized = .TRUE.
      END IF
    CASE (3)
      IF (PRESENT(mmr_gas) .AND. PRESENT(xp) .AND. PRESENT(pressure)) THEN
        zx_m = (mmr_gas+xp(1)*mmr_gas)*0.5_wp
        zx_d = (mmr_gas-xp(1)*mmr_gas)*0.5_wp
        gas_profile(1:jce,:)=(1._wp-(zx_d/zx_m)*TANH(LOG(pressure(1:jce,:)   &
          &                     /xp(2)) /xp(3))) * zx_m
        gas_initialized = .TRUE.
      END IF
    CASE (4)
      IF (PRESENT(gas_scenario)) THEN
        ! TODO: (Hauke Schmidt)
        ! If the respective parameters are present, a vertical
        ! profile is calculated as in option (3). This allows a seamless
        ! continuation of preindustrial control with scenarios. The treatment here is
        ! inconsistent with having two different options for the constant
        ! concentration cases (2 without and 3 with profile). However, instead
        ! of adding a fifth option, it seems more advisable to clean up the
        ! complete handling of radiation switches (including ighg), later.
        IF (PRESENT(xp) .AND. PRESENT(pressure)) THEN
          zx_m = (gas_scenario+xp(1)*gas_scenario)*0.5_wp
          zx_d = (gas_scenario-xp(1)*gas_scenario)*0.5_wp
          gas_profile(1:jce,:)=(1._wp-(zx_d/zx_m)*TANH(LOG(pressure(1:jce,:)   &
            &                     /xp(2)) /xp(3))) * zx_m
        ELSE
          gas_profile(1:jce,:) = gas_scenario
        END IF
        gas_initialized = .TRUE.
      ELSE IF (PRESENT(gas_scenario_v)) THEN
        gas_profile(1:jce,:) = gas_scenario_v(1:jce,:)
        gas_initialized = .TRUE.
      END IF
    END SELECT
    IF (.NOT. gas_initialized) &
      CALL finish('radiation','gas_profile options not supported')

  END FUNCTION gas_profile
  !-----------------------------------------------------------------------------
  !>
  !! @brief arranges input and calls rrtm sw and lw routines
  !!
  !! @par Revision History
  !! Original Source Rewritten and renamed by B. Stevens (2009-08)
  !!
  !! @remarks
  !!   Because the RRTM indexes vertical levels differently than ECHAM a chief
  !!   function of this routine is to reorder the input in the vertical.  In
  !!   addition some cloud physical properties are prescribed, which are
  !!   required to derive cloud optical properties
  !!
  !! @par The gases are passed into RRTM via two multi-constituent arrays:
  !!   zwkl and wx_r. zwkl has JPINPX species and  wx_r has JPXSEC species
  !!   The species are identifed as follows.
  !!     ZWKL [#/cm2]          WX_R [#/cm2]
  !!    index = 1 => H20     index = 1 => n/a
  !!    index = 2 => CO2     index = 2 => CFC11
  !!    index = 3 =>  O3     index = 3 => CFC12
  !!    index = 4 => N2O     index = 4 => n/a
  !!    index = 5 => n/a
  !!    index = 6 => CH4
  !!    index = 7 => O2
  !


  SUBROUTINE rrtm_interface(                                              &
    ! input
    & jg              ,jb                                                ,&
    & jce             ,kbdim           ,klev                             ,&
    & ktype           ,zland           ,zglac                            ,&
    & pmu0                                                               ,&
    & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
    & emis_rad        ,z_mc                                              ,&
    & pp_fl           ,pp_hl           ,pp_sfc          ,tk_fl           ,&
    & tk_hl           ,tk_sfc          ,xm_vap                           ,&
    & xm_liq          ,xm_ice                                            ,&
    & cdnc                                                               ,&
    & cld_frc                                                            ,&
    & xm_o3           ,xm_co2          ,xm_ch4                           ,&
    & xm_n2o          ,xm_cfc11        ,xm_cfc12        ,xm_o2           ,&
    & zaeq1, zaeq2, zaeq3, zaeq4, zaeq5                                  ,&
    ! output
    & flx_lw_net      ,flx_sw_net      ,flx_lw_net_clr  ,flx_sw_net_clr  ,&
    & flx_uplw_sfc    ,flx_upsw_sfc    ,flx_uplw_sfc_clr,flx_upsw_sfc_clr,&
    & flx_dnsw_diff_sfc, flx_upsw_toa  ,flx_par_sfc                       )

    INTEGER,INTENT(in)  ::                &
      &  jg,                              & !< domain index
      &  jb,                              & !< block index
      &  jce,                             & !< number of columns
      &  kbdim,                           & !< first dimension of 2-d arrays
      &  klev                               !< number of levels

    INTEGER,INTENT(in)  ::                &
      &  ktype(kbdim)                       !< type of convection

    REAL(wp),INTENT(in) ::                &
      &  zland(kbdim),                    & !< land-sea mask. (1. = land, 0. = sea/lakes)
      &  zglac(kbdim),                    & !< fraction of land covered by glaciers
      &  pmu0(kbdim),                     & !< mu0 for solar zenith angle
      &  alb_vis_dir(kbdim),              & !< surface albedo for vis range and dir light
      &  alb_nir_dir(kbdim),              & !< surface albedo for NIR range and dir light
      &  alb_vis_dif(kbdim),              & !< surface albedo for vis range and dif light
      &  alb_nir_dif(kbdim),              & !< surface albedo for NIR range and dif light
      &  emis_rad(kbdim),                 & !< longwave surface emissivity
      &  z_mc(kbdim,klev),                & !< height at full levels [m]
      &  pp_fl(kbdim,klev),               & !< full level pressure in Pa
      &  pp_hl(kbdim,klev+1),             & !< half level pressure in Pa
      &  pp_sfc(kbdim),                   & !< surface pressure in Pa
      &  tk_fl(kbdim,klev),               & !< full level temperature in K
      &  tk_hl(kbdim,klev+1),             & !< half level temperature in K
      &  tk_sfc(kbdim),                   & !< surface temperature in K
      &  xm_vap(kbdim,klev),              & !< specific humidity in g/g
      &  xm_liq(kbdim,klev),              & !< specific liquid water content
      &  xm_ice(kbdim,klev),              & !< specific ice content in g/g
      &  cdnc(kbdim,klev),                & !< cloud nuclei concentration
      &  cld_frc(kbdim,klev),             & !< fractional cloud cover
      &  xm_o3(kbdim,klev),               & !< o3 mass mixing ratio
      &  xm_co2(kbdim,klev),              & !< co2 mass mixing ratio
      &  xm_ch4(kbdim,klev),              & !< ch4 mass mixing ratio
      &  xm_n2o(kbdim,klev),              & !< n2o mass mixing ratio
      &  xm_cfc11(kbdim,klev),            & !< cfc 11 volume mixing ratio
      &  xm_cfc12(kbdim,klev),            & !< cfc 12 volume mixing ratio
      &  xm_o2(kbdim,klev),               & !< o2  mass mixing ratio
      &  zaeq1(kbdim,klev),               & !< aerosol continental
      &  zaeq2(kbdim,klev),               & !< aerosol maritime
      &  zaeq3(kbdim,klev),               & !< aerosol urban
      &  zaeq4(kbdim,klev),               & !< aerosol volcano ashes
      &  zaeq5(kbdim,klev)                  !< aerosol stratospheric background


    REAL(wp), INTENT(out) ::              &
      &  flx_lw_net(kbdim,klev+1),        & !< net downward LW flux profile,
      &  flx_sw_net(kbdim,klev+1),        & !< net downward SW flux profile,
      &  flx_lw_net_clr(kbdim,klev+1),    & !< clrsky downward LW flux profile,
      &  flx_sw_net_clr(kbdim,klev+1),    & !< clrsky downward SW flux profile,
      &  flx_uplw_sfc(kbdim),             & !< sfc LW upward flux,
      &  flx_upsw_sfc(kbdim),             & !< sfc SW upward flux,
      &  flx_uplw_sfc_clr(kbdim),         & !< clrsky sfc LW upward flux,
      &  flx_upsw_sfc_clr(kbdim)            !< clrsky sfc SW upward flux,

    REAL(wp), INTENT(out), OPTIONAL ::    &
      &  flx_dnsw_diff_sfc(kbdim),        & !< sfc SW diffuse downward flux,
      &  flx_upsw_toa(kbdim),             & !< TOA SW upward flux,
      &  flx_par_sfc(kbdim)                 !< PAR downward sfc flux

    INTEGER  :: jk, jl, jp, jkb, jspec,   & !< loop indicies
      &  icldlyr(kbdim,klev)                !< index for clear or cloudy

    REAL(wp) ::                           &
      &  zsemiss(kbdim,jpband),           & !< LW surface emissivity by band
      &  ppd_hl(kbdim,klev),              & !< pressure thickness in Pa
      &  pm_sfc(kbdim),                   & !< surface pressure in hPa
      &  amm,                             & !< molecular weight of moist air
      &  delta,                           & !< pressure thickness
      &  zscratch                           !< scratch array

    REAL(wp) :: z_sum_aea, z_sum_aes !help variables for aerosol
    
    !
    ! --- vertically reversed _vr variables
    !
    REAL(wp) ::                           &
      &  col_dry_vr(kbdim,klev),          & !< number of molecules/cm2 of
      &  pm_fl_vr(kbdim,klev),            & !< full level pressure [hPa]
      &  pm_hl_vr(kbdim,klev+1),          & !< half level pressure [hPa]
      &  tk_fl_vr(kbdim,klev),            & !< full level temperature [K]
      &  tk_hl_vr(kbdim,klev+1),          & !< half level temperature [K]
      &  cdnc_vr(kbdim,klev),             & !< cloud nuclei concentration
      &  cld_frc_vr(kbdim,klev),          & !< secure cloud fraction
      &  ziwgkg_vr(kbdim,klev),           & !< specific ice water content
      &  ziwc_vr(kbdim,klev),             & !< ice water content per volume
      &  ziwp_vr(kbdim,klev),             & !< ice water path in g/m2
      &  zlwgkg_vr(kbdim,klev),           & !< specific liquid water content
      &  zlwp_vr(kbdim,klev),             & !< liquid water path in g/m2
      &  zlwc_vr(kbdim,klev),             & !< liquid water content per
      &  wkl_vr(kbdim,jpinpx,klev),       & !< number of molecules/cm2 of
      &  wx_vr(kbdim,jpxsec,klev),        & !< number of molecules/cm2 of
      &  cld_tau_lw_vr(kbdim,klev,jpband),& !< LW optical thickness of clouds
      &  cld_tau_sw_vr(kbdim,jpsw,klev),  & !< extincion
      &  cld_cg_sw_vr(kbdim,jpsw,klev),   & !< asymmetry factor
      &  cld_piz_sw_vr(kbdim,jpsw,klev),  & !< single scattering albedo
      &  aer_tau_lw_vr(kbdim,klev,jpband),& !< LW optical thickness of aerosols
      &  aer_tau_sw_vr(kbdim,klev,jpsw),  & !< aerosol optical thickness
      &  aer_cg_sw_vr(kbdim,klev,jpsw),   & !< aerosol asymmetry factor
      &  aer_piz_sw_vr(kbdim,klev,jpsw),  & !< aerosol single scattering albedo
      &  flx_uplw_vr(kbdim,klev+1),       & !< upward flux, total sky
      &  flx_uplw_clr_vr(kbdim,klev+1),   & !< upward flux, clear sky
      &  flx_dnlw_vr(kbdim,klev+1),       & !< downward flux, total sky
      &  flx_dnlw_clr_vr(kbdim,klev+1),   & !< downward flux, clear sky
      &  flx_upsw(kbdim,klev+1),          & !< upward flux total sky
      &  flx_upsw_clr(kbdim,klev+1),      & !< upward flux clear sky
      &  flx_dnsw(kbdim,klev+1),          & !< downward flux total sky
      &  flx_dnsw_clr(kbdim,klev+1)         !< downward flux clear sky

    CHARACTER(LEN=3)     :: c_irad_aero

    ! Initialize output variables
    flx_lw_net(:,:)     = 0._wp
    flx_lw_net_clr(:,:) = 0._wp
    flx_sw_net(:,:)     = 0._wp
    flx_sw_net_clr(:,:) = 0._wp
    flx_uplw_sfc(:)     = 0._wp
    flx_uplw_sfc_clr(:) = 0._wp
    flx_upsw_sfc(:)     = 0._wp
    flx_upsw_sfc_clr(:) = 0._wp
    IF (PRESENT(flx_dnsw_diff_sfc)) flx_dnsw_diff_sfc(:) = 0._wp
    IF (PRESENT(flx_upsw_toa))      flx_upsw_toa(:)      = 0._wp
    IF (PRESENT(flx_par_sfc))       flx_par_sfc(:)       = 0._wp

    !
    ! 1.0 Constituent properties
    !--------------------------------

    !
    ! --- control for infintesimal cloud fractions
    !
    DO jk = 1, klev
      !
      jkb = klev+1-jk
      cld_frc_vr(1:jce,jk)  = cld_frc(1:jce,jkb)

! LL: Apparently xlf gets confused with this where, replace it by IFs
#ifdef __xlC__
      icldlyr  (:,jk) = 0
      ziwgkg_vr(:,jk) = 0.0_wp
      zlwgkg_vr(:,jk) = 0.0_wp

      DO jl=1,jce
        IF (cld_frc_vr(jl,jk) > 2.0_wp*EPSILON(1.0_wp)) THEN
          icldlyr  (jl,jk) = 1
          ziwgkg_vr(jl,jk) = xm_ice(jl,jkb)*1000.0_wp/cld_frc_vr(jl,jk)
          zlwgkg_vr(jl,jk) = xm_liq(jl,jkb)*1000.0_wp/cld_frc_vr(jl,jk)
        ENDIF
      END DO
#else
      !
      WHERE (cld_frc_vr(1:jce,jk) > 2.0_wp*EPSILON(1.0_wp))
        ! only clouds > 2 epsilon are made visible to radiation
        icldlyr  (1:jce,jk) = 1
        ziwgkg_vr(1:jce,jk) = xm_ice(1:jce,jkb)*1000.0_wp/cld_frc_vr(1:jce,jk)
        zlwgkg_vr(1:jce,jk) = xm_liq(1:jce,jkb)*1000.0_wp/cld_frc_vr(1:jce,jk)
      ELSEWHERE
        ! clouds <= 2 epsilon are ade invisble to radiation
        icldlyr  (1:jce,jk) = 0
        ziwgkg_vr(1:jce,jk) = 0.0_wp
        zlwgkg_vr(1:jce,jk) = 0.0_wp
      END WHERE
#endif
    END DO
    !
    ! --- main constituent reordering
    !
    DO jl = 1, jce
      pm_hl_vr(jl,klev+1) = 0.01_wp*pp_hl(jl,1)
      tk_hl_vr(jl,klev+1) = tk_hl(jl,1)
      pm_sfc(jl)          = 0.01_wp*pp_sfc(jl)
    END DO

    DO jk = 1, klev
      jkb = klev+1-jk
      ! initialization
      wkl_vr(:,:,jk) = 0.0_wp
      wx_vr (:,:,jk) = 0.0_wp
      DO jl = 1, jce
        delta = pp_hl(jl,jkb+1)-pp_hl(jl,jkb)
        !
        ! --- thermodynamic arrays
        !
        pm_hl_vr(jl,jk) = 0.01_wp*pp_hl(jl,jkb+1)
        pm_fl_vr(jl,jk) = 0.01_wp*pp_fl(jl,jkb)
        tk_hl_vr(jl,jk) = tk_hl(jl,jkb+1)
        tk_fl_vr(jl,jk) = tk_fl(jl,jkb)
        !
        ! --- cloud properties
        !
        zscratch       = pp_fl(jl,jkb)/tk_fl(jl,jkb)
        ziwc_vr(jl,jk) = ziwgkg_vr(jl,jk)*zscratch/rd
        ziwp_vr(jl,jk) = ziwgkg_vr(jl,jk)*delta/grav
        zlwc_vr(jl,jk) = zlwgkg_vr(jl,jk)*zscratch/rd
        zlwp_vr(jl,jk) = zlwgkg_vr(jl,jk)*delta/grav
        cdnc_vr(jl,jk) = cdnc(jl,jkb)*1.e-6_wp
        !
        ! --- radiatively active gases
        !
        wkl_vr(jl,1,jk)   = xm_vap(jl,jkb)*amd/amw
        wkl_vr(jl,2,jk)   = xm_co2(jl,jkb)*amd/amco2
        wkl_vr(jl,3,jk)   = xm_o3(jl,jkb) *amd/amo3
        wkl_vr(jl,4,jk)   = xm_n2o(jl,jkb)*amd/amn2o
        wkl_vr(jl,6,jk)   = xm_ch4(jl,jkb)*amd/amch4
        wkl_vr(jl,7,jk)   = xm_o2 (jl,jkb)*amd/amo2
        amm               = (1.0_wp-wkl_vr(jl,1,jk))*amd + wkl_vr(jl,1,jk)*amw
        col_dry_vr(jl,jk) = (0.01_wp*delta)*10.0_wp*avo/grav/amm / (1.0_wp+wkl_vr(jl,1,jk))
        !
        ! --- alternate treatment for cfcs
        !
        wx_vr(jl,2,jk) = col_dry_vr(jl,jk)*xm_cfc11(jl,jkb)*1.e-20_wp
        wx_vr(jl,3,jk) = col_dry_vr(jl,jk)*xm_cfc12(jl,jkb)*1.e-20_wp
      END DO
    END DO
    DO jp = 1, 7
      wkl_vr(1:jce,jp,:)=col_dry_vr(1:jce,:)*wkl_vr(1:jce,jp,:)
    END DO
    !
    ! 2.0 Surface Properties
    ! --------------------------------
    zsemiss(1:jce,:) = SPREAD(emis_rad(1:jce),2,jpband)
    !
    ! 3.0 Particulate Optical Properties
    ! --------------------------------
    ppd_hl(1:jce,:) = pp_hl(1:jce,2:klev+1)-pp_hl(1:jce,1:klev)

    SELECT CASE (irad_aero)
    CASE (0,2)
      aer_tau_lw_vr(:,:,:) = 0.0_wp
      aer_tau_sw_vr(:,:,:) = 0.0_wp
      aer_piz_sw_vr(:,:,:) = 1.0_wp
      aer_cg_sw_vr(:,:,:)  = 0.0_wp
    CASE (5,6)
      DO jspec=1,jpband
        DO jk=1,klev
          jkb = klev+1-jk
          DO jl = 1,jce
            ! LW opt thickness of aerosols
            aer_tau_lw_vr(jl,jk,jspec) =  zaeq1(jl,jkb) * zaea_rrtm(jspec,1) &
              &                         + zaeq2(jl,jkb) * zaea_rrtm(jspec,2) &
              &                         + zaeq3(jl,jkb) * zaea_rrtm(jspec,3) &
              &                         + zaeq4(jl,jkb) * zaea_rrtm(jspec,4) &
              &                         + zaeq5(jl,jkb) * zaea_rrtm(jspec,5)
          ENDDO
        ENDDO
      ENDDO
      DO jspec=1+jpband,jpband+jpsw
        DO jk=1,klev
          jkb = klev+1-jk
          DO jl = 1,jce

            z_sum_aea = zaeq1(jl,jkb) * zaea_rrtm(jspec,1) &
              &       + zaeq2(jl,jkb) * zaea_rrtm(jspec,2) &
              &       + zaeq3(jl,jkb) * zaea_rrtm(jspec,3) &
              &       + zaeq4(jl,jkb) * zaea_rrtm(jspec,4) &
              &       + zaeq5(jl,jkb) * zaea_rrtm(jspec,5)

            z_sum_aes = zaeq1(jl,jkb) * zaes_rrtm(jspec,1) &
              &       + zaeq2(jl,jkb) * zaes_rrtm(jspec,2) &
              &       + zaeq3(jl,jkb) * zaes_rrtm(jspec,3) &
              &       + zaeq4(jl,jkb) * zaes_rrtm(jspec,4) &
              &       + zaeq5(jl,jkb) * zaes_rrtm(jspec,5)

            ! sw aerosol optical thickness
            aer_tau_sw_vr(jl,jk,jspec-jpband) = z_sum_aea + z_sum_aes

            ! sw aerosol single scattering albedo
            aer_piz_sw_vr(jl,jk,jspec-jpband) = z_sum_aes / ( z_sum_aea + z_sum_aes ) 

            ! sw aerosol asymmetry factor
            aer_cg_sw_vr(jl,jk,jspec-jpband) =                                  &
              & (   zaeq1(jl,jkb) * zaes_rrtm(jspec,1) * zaeg_rrtm(jspec,1)   &
              &   + zaeq2(jl,jkb) * zaes_rrtm(jspec,2) * zaeg_rrtm(jspec,2)   &
              &   + zaeq3(jl,jkb) * zaes_rrtm(jspec,3) * zaeg_rrtm(jspec,3)   &
              &   + zaeq4(jl,jkb) * zaes_rrtm(jspec,4) * zaeg_rrtm(jspec,4)   &
              &   + zaeq5(jl,jkb) * zaes_rrtm(jspec,5) * zaeg_rrtm(jspec,5) ) / z_sum_aes
          ENDDO
        ENDDO
      ENDDO
    CASE (13)
      CALL set_bc_aeropt_kinne( jg,                                  &
        & jce              ,kbdim                 ,klev             ,&
        & jb               ,jpband                ,jpsw             ,&
        & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
        & aer_cg_sw_vr     ,ppd_hl                ,pp_fl            ,&
        & tk_fl                                                      )
    CASE (14)
      ! set zero aerosol before adding Stenchikov aerosols
      aer_tau_lw_vr(:,:,:) = 0.0_wp
      aer_tau_sw_vr(:,:,:) = 0.0_wp
      aer_piz_sw_vr(:,:,:) = 1.0_wp
      aer_cg_sw_vr(:,:,:)  = 0.0_wp
      CALL add_bc_aeropt_stenchikov( jg,                             &
        & jce              ,kbdim                 ,klev             ,&
        & jb               ,jpband                ,jpsw             ,&
        & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
        & aer_cg_sw_vr     ,ppd_hl                ,pp_fl            ,&
        & tk_fl                                                      )
    CASE (15)
      CALL set_bc_aeropt_kinne( jg,                                  &
        & jce              ,kbdim                 ,klev             ,&
        & jb               ,jpband                ,jpsw             ,&
        & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
        & aer_cg_sw_vr     ,ppd_hl                ,pp_fl            ,&
        & tk_fl                                                      )
      CALL add_bc_aeropt_stenchikov( jg,                             &
        & jce              ,kbdim                 ,klev             ,&
        & jb               ,jpband                ,jpsw             ,&
        & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
        & aer_cg_sw_vr     ,ppd_hl                ,pp_fl            ,&
        & tk_fl                                                      )
    CASE DEFAULT
      WRITE (c_irad_aero,'(i3)') irad_aero
      CALL finish ('rrtm_interface of mo_radition','irad_aero= '// &
                   TRIM(ADJUSTL(c_irad_aero))//' does not exist')
    END SELECT
    IF (lrad_aero_diag) THEN
      CALL rad_aero_diag (                                  &
      & jg              ,jb              ,jce             , &
      & kbdim           ,klev            ,jpband          , &
      & jpsw            ,aer_tau_lw_vr   ,aer_tau_sw_vr   , &
      & aer_piz_sw_vr   ,aer_cg_sw_vr                       )
    END IF


    CALL newcld_optics(                                                       &
      & jce          ,kbdim        ,klev         ,jpband       ,jpsw         ,&
      & zglac        ,zland        ,ktype        ,icldlyr      ,tk_fl_vr     ,&
      & zlwp_vr      ,ziwp_vr      ,zlwc_vr      ,ziwc_vr      ,cdnc_vr      ,&
      & cld_tau_lw_vr,cld_tau_sw_vr,cld_piz_sw_vr,cld_cg_sw_vr                )


    !
    ! 4.0 Radiative Transfer Routines
    ! --------------------------------
    CALL lrtm(                                                                &
      !    input
      &    jce             ,klev                                             ,&
      &    pm_fl_vr        ,pm_sfc          ,tk_fl_vr        ,tk_hl_vr       ,&
      &    tk_sfc          ,wkl_vr          ,wx_vr           ,col_dry_vr     ,&
      &    zsemiss         ,cld_frc_vr      ,cld_tau_lw_vr   ,aer_tau_lw_vr  ,&
      !    output
      &    flx_uplw_vr     ,flx_dnlw_vr     ,flx_uplw_clr_vr,flx_dnlw_clr_vr )


    CALL srtm_srtm_224gp(                                                     &
      !    input
      &    jce             ,kbdim           ,klev            ,jpsw           ,&
      &    alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif    ,&
      &    pm_fl_vr        ,tk_fl_vr        ,pmu0                            ,&
      &    col_dry_vr      ,wkl_vr                                           ,&
      &    cld_frc_vr      ,cld_tau_sw_vr   ,cld_cg_sw_vr    ,cld_piz_sw_vr  ,&
      &    aer_tau_sw_vr   ,aer_cg_sw_vr    ,aer_piz_sw_vr                   ,&
      &    ssi_radt        ,z_mc                                             ,&
      !    output
      &    flx_dnsw        ,flx_upsw        ,flx_dnsw_clr    ,flx_upsw_clr,   &
      &    flx_dnsw_diff_sfc ,flx_par_sfc                                     )


    ! 5.0 Post Processing
    ! --------------------------------
    DO jk = 1, klev+1
      jkb = klev+2-jk
      DO jl = 1, jce
        flx_lw_net(jl,jk)     = flx_dnlw_vr(jl,jkb)-flx_uplw_vr(jl,jkb)
        flx_lw_net_clr(jl,jk) = flx_dnlw_clr_vr(jl,jkb)-flx_uplw_clr_vr(jl,jkb)
        flx_sw_net(jl,jk)     = flx_dnsw(jl,jk) - flx_upsw(jl,jk)
        flx_sw_net_clr(jl,jk) = flx_dnsw_clr(jl,jk)-flx_upsw_clr(jl,jk)
      END DO
    END DO
    flx_uplw_sfc(1:jce)     = flx_uplw_vr(1:jce,1)
    flx_uplw_sfc_clr(1:jce) = flx_uplw_clr_vr(1:jce,1)
    flx_upsw_sfc(1:jce)     = flx_upsw(1:jce,klev+1)
    flx_upsw_sfc_clr(1:jce) = flx_upsw_clr(1:jce,klev+1)
    IF (PRESENT(flx_upsw_toa)) flx_upsw_toa(1:jce) = flx_upsw(1:jce,1)
!!$    sw_irr_toa(1:jce)       = flx_dnsw(1:jce,1)
    !

  END SUBROUTINE rrtm_interface

  !-----------------------------------------------------------------------------
  !>
  !! Compute shortwave and longwave heating rates
  !!
  !! The radheat subroutine computes the radiative heating rates resulting from
  !! the divergence of the vertical profiles of longwave and shortwave net fluxes.
  !!
  !! - Shortwave net flux profiles are computed from:
  !!   - the vertical profiles of net transmissivity
  !!   - the solar incoming flux at TOA
  !! - Longwave net flux profiles are given as input
  !! - Specific heat depends on the moisture in the air
  !!
  !! @author Marco Giorgetta, Max Planck Institute for Meteorology
  !!
  !!
  !! @par Revision History
  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
  !!

  SUBROUTINE radheat (jcs, jce, kbdim, &
    &                 klev  , klevp1,  &
    &                 ntiles        ,  &
    &                 ntiles_wtr    ,  &
    &                 pmair         ,  &
    &                 pqv           ,  &
    &                 pcd           ,  &
    &                 pcv           ,  &
    &                 pi0           ,  &
    &                 pemiss        ,  &
    &                 ptsfc         ,  &
    &                 ptsfc_t       ,  & ! optional: tile-specific ground temperature
    &                 ptsfctrad     ,  &
    &                 pqc           ,  & ! optional ! must be present if opt_nh_corr=.true.
    &                 pqi           ,  & ! optional ! must be present if opt_nh_corr=.true.
    &                 ppres_ifc     ,  & ! optional ! must be present if opt_nh_corr=.true.
    &                 albedo, albedo_t,& ! optional: albedo fields
    &                 lwflx_up_sfc_rs, & ! optional: longwave upward flux at surface
    &                 trsol_up_toa,    & ! optional: normalized shortwave upward flux at the top of the atmosphere
    &                 trsol_up_sfc,    & ! optional: normalized shortwave upward flux at the surface
    &                 trsol_par_sfc,   & ! optional: normalized photosynthetically active downward flux at the surface
    &                 trsol_dn_sfc_diff,&! optional: normalized shortwave diffuse downward radiative flux at the surface
    &                 trsol_clr_sfc,   & ! optional: normalized shortwave clear-sky net radiative flux at the surface
    &                 lp_count,        & ! optional: number of land points
    &                 gp_count_t,      & ! optional: number of land points per tile
    &                 spi_count,       & ! optional: number of seaice points
    &                 fp_count,        & ! optional: number of lake points
    &                 idx_lst_lp,      & ! optional: index list of land points
    &                 idx_lst_t,       & ! optional: index list of land points per tile
    &                 idx_lst_spi,     & ! optional: index list of seaice points
    &                 idx_lst_fp,      & ! optional: index list of (f)lake points
    &                 cosmu0,          & ! optional: cosine of zenith angle
    &                 opt_nh_corr   ,  & ! optional: switch for applying corrections for NH model
    &                 use_trsolclr_sfc,& ! optional: use clear-sky surface transmissivity passed on input
    &                 ptrmsw        ,  &
    &                 pflxlw        ,  &
    &                 pdtdtradsw    ,  &
    &                 pdtdtradlw    ,  &
    &                 pflxsfcsw     ,  &
    &                 pflxsfclw     ,  &
    &                 pflxsfcsw_t   ,  &
    &                 pflxsfclw_t   ,  &
    &                 pflxtoasw     ,  &
    &                 pflxtoalw     ,  &
    &                 lwflx_up_sfc  ,  &
    &                 swflx_up_toa  ,  &
    &                 swflx_up_sfc  ,  &
    &                 swflx_par_sfc ,  &
    &                 swflx_dn_sfc_diff)

    INTEGER,  INTENT(in)  ::    &
      &     jcs, jce, kbdim,    &
      &     klev,   klevp1, ntiles, ntiles_wtr

    REAL(wp), INTENT(in)  ::         &
      &     pmair      (kbdim,klev), & ! mass of air in layer                     [kg/m2]
      &     pqv        (kbdim,klev), & ! specific humidity at t-dt                [kg/kg]
      &     pcd,                     & ! specific heat of dry air                 [J/kg/K]
      &     pcv,                     & ! specific heat of vapor                   [J/kg/K]
      &     pi0        (kbdim),      & ! local solar incoming flux at TOA         [W/m2]
      &     pemiss     (kbdim),      & ! lw sfc emissivity
      &     ptsfc      (kbdim),      & ! surface temperature at t                 [K]
      &     ptsfctrad  (kbdim),      & ! surface temperature at trad              [K]
      &     ptrmsw     (kbdim,klevp1), & ! shortwave transmissivity at trad         []
      &     pflxlw     (kbdim,klevp1)    ! longwave net flux at trad                [W/m2]

    REAL(wp), INTENT(in), OPTIONAL  ::         &
      &     pqc       (kbdim,klev),  & ! specific cloud water               [kg/kg]
      &     pqi       (kbdim,klev),  & ! specific cloud ice                 [kg/kg]
      &     ppres_ifc (kbdim,klevp1),& ! pressure at interfaces             [Pa]
      &     ptsfc_t   (kbdim,ntiles+ntiles_wtr),& ! tile-specific surface temperature at t  [K]
      &     cosmu0    (kbdim),       & ! cosine of solar zenith angle
      &     albedo    (kbdim),       & ! grid-box average albedo
      &     albedo_t  (kbdim,ntiles+ntiles_wtr), &   ! tile-specific albedo
      &     lwflx_up_sfc_rs(kbdim),& ! longwave upward flux at surface calculated at radiation time steps
      &     trsol_up_toa(kbdim),   & ! normalized shortwave upward flux at the top of the atmosphere
      &     trsol_up_sfc(kbdim),   & ! normalized shortwave upward flux at the surface
      &     trsol_par_sfc(kbdim),  & ! normalized photosynthetically active downward flux at the surface
      &     trsol_clr_sfc(kbdim),  & ! normalized shortwave clear-sky net radiative flux at the surface
      &     trsol_dn_sfc_diff(kbdim) ! normalized shortwave diffuse downward radiative flux at the surface

    INTEGER, INTENT(in), OPTIONAL  ::     &
      &     lp_count, gp_count_t(ntiles), &  ! number of land points
      &     spi_count,                    &  ! number of seaice points
      &     fp_count,                     &  ! number of lake points
      &     idx_lst_lp(kbdim), idx_lst_t(kbdim,ntiles),& ! corresponding index lists
      &     idx_lst_spi(kbdim),           &  ! sea-ice point index list
      &     idx_lst_fp(kbdim)                ! (f)lake point index list

    LOGICAL, INTENT(in), OPTIONAL   ::  &
      &     opt_nh_corr, use_trsolclr_sfc

    REAL(wp), INTENT(inout) ::       &
      &     pdtdtradsw (kbdim,klev), & ! shortwave temperature tendency           [K/s]
      &     pdtdtradlw (kbdim,klev)    ! longwave temperature tendency            [K/s]

    REAL(wp), INTENT(inout), OPTIONAL :: &
      &     pflxsfcsw (kbdim), &       ! shortwave surface net flux [W/m2]
      &     pflxsfclw (kbdim), &       ! longwave  surface net flux [W/m2]
      &     pflxsfcsw_t(kbdim,ntiles+ntiles_wtr), & ! tile-specific shortwave 
                                                    ! surface net flux [W/m2]
      &     pflxsfclw_t(kbdim,ntiles+ntiles_wtr), & ! tile-specific longwave 
                                                    ! surface net flux [W/m2]
      &     pflxtoasw (kbdim), &       ! shortwave toa net flux [W/m2]
      &     pflxtoalw (kbdim), &       ! longwave  toa net flux [W/m2]
      &     lwflx_up_sfc(kbdim), &     ! longwave upward flux at surface [W/m2]
      &     swflx_up_toa(kbdim), &     ! shortwave upward flux at the top of the atmosphere [W/m2]
      &     swflx_up_sfc(kbdim), &     ! shortwave upward flux at the surface [W/m2]
      &     swflx_par_sfc(kbdim), &    ! photosynthetically active downward flux at the surface [W/m2]
      &     swflx_dn_sfc_diff(kbdim)   ! shortwave diffuse downward radiative flux at the surface [W/m2]

    ! Local arrays
    REAL(wp) ::                    &
      &     zflxsw (kbdim,klevp1), &
      &     zflxlw (kbdim,klevp1), &
      &     zconv  (kbdim,klev)  , &
      &     tqv    (kbdim)       , &
      &     dlwem_o_dtg(kbdim)   , &
      &     lwfac1 (kbdim)       , &
      &     lwfac2 (kbdim)       , &
      &     intclw (kbdim,klevp1), &
      &     intcli (kbdim,klevp1), &
      &     dlwflxall_o_dtg(kbdim,klevp1)

    REAL(wp) :: swfac1(kbdim), swfac2(kbdim), dflxsw_o_dalb(kbdim), trsolclr(kbdim), logtqv(kbdim)

    ! local scalars
    REAL(wp) :: dpresg, pfaclw, intqctot, dlwflxclr_o_dtg

    REAL(wp), PARAMETER  :: pscal = 1._wp/4000._wp ! pressure scale for longwave correction

    INTEGER :: jc,jk,jt,ic

    LOGICAL  :: l_nh_corr, lcalc_trsolclr

    IF ( PRESENT(opt_nh_corr) ) THEN
      l_nh_corr = opt_nh_corr
      IF (l_nh_corr .AND. .NOT.(PRESENT(pqc).AND.PRESENT(pqi).AND.PRESENT(ppres_ifc))) THEN
        CALL finish('radheat','error: if opt_nh_corr, pqc, pqi, and ppres_ifc must be present.')
      ENDIF
    ELSE
      l_nh_corr = .FALSE.
    ENDIF

    lcalc_trsolclr = .TRUE.
    IF (PRESENT(use_trsolclr_sfc) .AND. PRESENT(trsol_clr_sfc)) THEN
      IF (use_trsolclr_sfc) lcalc_trsolclr = .FALSE.
    ENDIF

    ! Conversion factor for heating rates
    zconv(jcs:jce,1:klev) = 1._wp/(pmair(jcs:jce,1:klev)*(pcd+(pcv-pcd)*pqv(jcs:jce,1:klev)))

    ! Shortwave fluxes = transmissivity * local solar incoming flux at TOA
    ! ----------------
    ! - TOA
    zflxsw(jcs:jce,1)      = ptrmsw(jcs:jce,1)      *        pi0(jcs:jce)
    ! - Atmosphere
    zflxsw(jcs:jce,2:klev) = ptrmsw(jcs:jce,2:klev) * SPREAD(pi0(jcs:jce),2,klev-1)
    ! - Surface
    zflxsw(jcs:jce,klevp1) = ptrmsw(jcs:jce,klevp1) *        pi0(jcs:jce)
    ! Longwave fluxes
    ! - TOA
!    zflxlw(jcs:jce,1)      = pflxlw(jcs:jce,1)
!    ! - Atmosphere
!    zflxlw(jcs:jce,2:klev) = pflxlw(jcs:jce,2:klev)

    IF (l_nh_corr) THEN !

      ! Additional shortwave fluxes for NWP requirements
      DO jc = jcs, jce
        swflx_up_toa(jc)      = pi0(jc)*trsol_up_toa(jc)
        swflx_up_sfc(jc)      = pi0(jc)*trsol_up_sfc(jc)
        swflx_dn_sfc_diff(jc) = pi0(jc)*trsol_dn_sfc_diff(jc)
        swflx_par_sfc(jc)     = pi0(jc)*trsol_par_sfc(jc)
      ENDDO

      ! Correction of longwave fluxes for changes in ground temperature
      tqv(:)            = 0._wp
      intclw(:,klevp1)  = 0._wp
      intcli(:,klevp1)  = 0._wp
      DO jk = klev,1,-1
        DO jc = jcs, jce
          dpresg        = (ppres_ifc(jc,jk+1) - ppres_ifc(jc,jk))/grav
          tqv(jc)       = tqv(jc)+pqv(jc,jk)*dpresg
          intclw(jc,jk) = intclw(jc,jk+1)+pqc(jc,jk)*dpresg
          intcli(jc,jk) = intcli(jc,jk+1)+pqi(jc,jk)*dpresg
        ENDDO
      ENDDO

      DO jc = jcs, jce
        logtqv(jc) = LOG(MAX(1._wp,tqv(jc)))
        dlwem_o_dtg(jc) = pemiss(jc)*4._wp*stbo*ptsfc(jc)**3
        lwflx_up_sfc(jc) = lwflx_up_sfc_rs(jc) + dlwem_o_dtg(jc)*(ptsfc(jc) - ptsfctrad(jc))
        lwfac2(jc) = 0.92_wp*EXP(-0.07_wp*logtqv(jc))
      ENDDO
      DO jc = jcs, jce
        IF (tqv(jc) > 15._wp) then
          lwfac1(jc) = 1.677_wp*EXP(-0.72_wp*logtqv(jc))
        ELSE
          lwfac1(jc) = 0.4388_wp*EXP(-0.225_wp*logtqv(jc))
        ENDIF
      ENDDO

      DO jk = 1,klevp1
        DO jc = jcs, jce
          pfaclw = lwfac1(jc)+(lwfac2(jc)-lwfac1(jc))*EXP(-SQRT((ppres_ifc(jc,klevp1)- &
            ppres_ifc(jc,jk))*pscal))
          intqctot = MIN(0.30119_wp,MAX(1.008e-3_wp,intclw(jc,jk)+0.2_wp*intcli(jc,jk)))

          dlwflxclr_o_dtg = -dlwem_o_dtg(jc)*pfaclw

          ! derivative of LW flux w.r.t. ground temperature
          dlwflxall_o_dtg(jc,jk) = dlwflxclr_o_dtg*(1._wp-(6.9_wp+LOG(intqctot))/5.7_wp)
          ! Now apply the correction
          zflxlw(jc,jk) =  pflxlw(jc,jk) + dlwflxall_o_dtg(jc,jk) * (ptsfc(jc) - ptsfctrad(jc))
        ENDDO
      ENDDO

      ! Disaggregation of longwave and shortwave fluxes for tile approach
      IF (ntiles > 1) THEN
        IF (lcalc_trsolclr) THEN ! (relevant for Ritter-Geleyn radiation scheme only)
          ! parameterization of clear-air solar transmissivity in order to use the same 
          ! formulation as in mo_phys_nest_utilities:downscale_rad_output
          DO jc = jcs, jce
            trsolclr(jc) = MAX(0.02_wp,0.8_wp*cosmu0(jc)/(0.25_wp*tqv(jc))**0.15)**0.333_wp*&
             (1._wp-albedo(jc))**(1._wp-0.2_wp*cosmu0(jc)+0.1_wp*MIN(10._wp,tqv(jc))**0.33_wp)
          ENDDO
        ELSE
          ! use clear-air solar transmissivity passed as argument
          trsolclr(jcs:jce) = trsol_clr_sfc(jcs:jce)
        ENDIF

        DO jc = jcs, jce
          swfac1(jc) = EXP(0.36_wp*LOG( MAX(1.e-3_wp,ptrmsw(jc,klevp1))/MAX(1.e-3_wp,trsolclr(jc)) ))
          swfac2(jc) = EXP(0.1_wp*LOG( MAX(0.25_wp,3._wp*cosmu0(jc)) ))

          ! derivative of SW surface flux w.r.t. albedo
          dflxsw_o_dalb(jc) = - zflxsw(jc,klevp1)*swfac1(jc)/((1._wp-albedo(jc))*swfac2(jc))
        ENDDO

        DO jt = 1,ntiles
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, gp_count_t(jt)
            jc = idx_lst_t(ic,jt)
            pflxsfcsw_t(jc,jt) = MAX(0.1_wp*zflxsw(jc,klevp1), zflxsw(jc,klevp1) + &
                                 dflxsw_o_dalb(jc)*(albedo_t(jc,jt)-albedo(jc)))
            pflxsfclw_t(jc,jt) = zflxlw(jc,klevp1) + dlwflxall_o_dtg(jc,klevp1)* &
                                 (ptsfc_t(jc,jt)-ptsfc(jc))
          ENDDO
        ENDDO

        ! seaice points
        !
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, spi_count
          jc = idx_lst_spi(ic)
          pflxsfcsw_t(jc,isub_seaice) = MAX(0.1_wp*zflxsw(jc,klevp1), zflxsw(jc,klevp1) &
            &                  + dflxsw_o_dalb(jc)*(albedo_t(jc,isub_seaice)-albedo(jc)))
          pflxsfclw_t(jc,isub_seaice) = zflxlw(jc,klevp1) + dlwflxall_o_dtg(jc,klevp1) &
            &                  * (ptsfc_t(jc,isub_seaice)-ptsfc(jc))
        ENDDO

        ! lake points
        !
!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, fp_count
          jc = idx_lst_fp(ic)
          pflxsfcsw_t(jc,isub_lake) = MAX(0.1_wp*zflxsw(jc,klevp1), zflxsw(jc,klevp1) &
            &                  + dflxsw_o_dalb(jc)*(albedo_t(jc,isub_lake)-albedo(jc)))
          pflxsfclw_t(jc,isub_lake) = zflxlw(jc,klevp1) + dlwflxall_o_dtg(jc,klevp1) &
            &                  * (ptsfc_t(jc,isub_lake)-ptsfc(jc))
        ENDDO

        ! (open) water points
        ! not needed, yet

      ELSE IF (PRESENT(pflxsfcsw_t) .AND. PRESENT(pflxsfclw_t)) THEN

!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, lp_count
          jc = idx_lst_lp(ic)
          pflxsfcsw_t(jc,1) = zflxsw(jc,klevp1)
          pflxsfclw_t(jc,1) = zflxlw(jc,klevp1) 
        ENDDO

!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, spi_count
          jc = idx_lst_spi(ic)
          pflxsfcsw_t(jc,1) = zflxsw(jc,klevp1)
          pflxsfclw_t(jc,1) = zflxlw(jc,klevp1)
        ENDDO

!CDIR NODEP,VOVERTAKE,VOB
        DO ic = 1, fp_count
          jc = idx_lst_fp(ic)
          pflxsfcsw_t(jc,1) = zflxsw(jc,klevp1)
          pflxsfclw_t(jc,1) = zflxlw(jc,klevp1)
        ENDDO

      ENDIF ! ntiles


    ELSE ! hydrostatic version

      ! For ECHAM physics: Returns flux divergences in [W/m2]
      ! instead of the heating rates in [K/s].
      zconv(:,:) = 1._wp

      ! Longwave fluxes: For now keep fluxes fixed at TOA and in atmosphere,
      ! but adjust flux from surface to the current surface temperature.
      ! - TOA
      zflxlw(jcs:jce,1)      = pflxlw(jcs:jce,1)
      ! - Atmosphere
      zflxlw(jcs:jce,2:klev) = pflxlw(jcs:jce,2:klev)
      
      ! - Surface
      !   Adjust for changed surface temperature (ptsfc) with respect to the
      !   surface temperature used for the longwave flux computation (ptsfctrad).
      !   --> modifies heating in lowermost layer only (is this smart?)
      zflxlw(jcs:jce,klevp1) = pflxlw(jcs:jce,klevp1)                      &
        &                   + pemiss(jcs:jce)*stbo * ptsfctrad(jcs:jce)**4 &
        &                   - pemiss(jcs:jce)*stbo * ptsfc    (jcs:jce)**4


    ENDIF
    
    !
    !
    !     4.2  Fluxes and heating rates except for lowest layer
    !
    pdtdtradsw(jcs:jce,1:klev) = (zflxsw(jcs:jce,1:klev)-zflxsw(jcs:jce,2:klev+1)) * &
      & zconv(jcs:jce,1:klev)
    pdtdtradlw(jcs:jce,1:klev) = (zflxlw(jcs:jce,1:klev)-zflxlw(jcs:jce,2:klev+1)) * &
      & zconv(jcs:jce,1:klev)

    !
    !     4.3 net fluxes at surface
    !
    IF ( PRESENT(pflxsfcsw) ) pflxsfcsw(jcs:jce) = zflxsw(jcs:jce,klevp1)
    IF ( PRESENT(pflxsfclw) ) pflxsfclw(jcs:jce) = zflxlw(jcs:jce,klevp1)

    !
    !     4.4 net sw flux at toa
    !
    IF ( PRESENT(pflxtoasw) ) pflxtoasw(jcs:jce) = zflxsw(jcs:jce,1)
    IF ( PRESENT(pflxtoalw) ) pflxtoalw(jcs:jce) = zflxlw(jcs:jce,1)

    
  END SUBROUTINE radheat

!>
!! SUBROUTINE rad_aero_diag writes actual aerosol optical properties to output stream
!!
!! @author J.S.Rast, MPI-Met Hamburg
!!
!! @par Revision History
!! Origianl Source by J.S.Rast, MPI-Met Hamburg, (2013-08-30)
!!
!!
SUBROUTINE rad_aero_diag (                                  &
      & kg              ,kb              ,kce             , &
      & kbdim           ,klev            ,kpband          , &
      & kpsw            ,paer_tau_lw_vr  ,paer_tau_sw_vr  , &
      & paer_piz_sw_vr  ,paer_cg_sw_vr                      )
      
      INTEGER, INTENT(in)    :: kg      ! domain index
      INTEGER, INTENT(in)    :: kb      ! block index
      INTEGER, INTENT(in)    :: kce     ! actual block length
      INTEGER, INTENT(in)    :: kbdim   ! declaration block length
      INTEGER, INTENT(in)    :: klev    ! levels
      INTEGER, INTENT(in)    :: kpband  ! number of lw bands
      INTEGER, INTENT(in)    :: kpsw    ! number of sw bands
      REAL(wp), INTENT(in)   :: paer_tau_lw_vr(kbdim,klev,kpband) ! aod thermal wavelengths
      REAL(wp), INTENT(in)   :: paer_tau_sw_vr(kbdim,klev,kpsw)   ! aod solar wavelengths
      REAL(wp), INTENT(in)   :: paer_piz_sw_vr(kbdim,klev,kpsw)   ! ssa solar wavelengths
      REAL(wp), INTENT(in)   :: paer_cg_sw_vr(kbdim,klev,kpsw)    ! asy solar wavelengths

      prm_field(kg)%aer_aod_9731(1:kce,1:klev,kb) = &
                   paer_tau_lw_vr(1:kce,klev:1:-1,7)
      prm_field(kg)%aer_aod_533 (1:kce,1:klev,kb) = &
                   paer_tau_sw_vr(1:kce,klev:1:-1,10)
      prm_field(kg)%aer_ssa_533 (1:kce,1:klev,kb) = &
                   paer_piz_sw_vr(1:kce,klev:1:-1,10)
      prm_field(kg)%aer_asy_533 (1:kce,1:klev,kb) = &
                   paer_cg_sw_vr(1:kce,klev:1:-1,10)
      prm_field(kg)%aer_aod_2325(1:kce,1:klev,kb) = &
                   paer_tau_sw_vr(1:kce,klev:1:-1,3)
      prm_field(kg)%aer_ssa_2325(1:kce,1:klev,kb) = &
                   paer_piz_sw_vr(1:kce,klev:1:-1,3)
      prm_field(kg)%aer_asy_2325(1:kce,1:klev,kb) = &
                   paer_cg_sw_vr(1:kce,klev:1:-1,3)
END SUBROUTINE rad_aero_diag
END MODULE mo_radiation
