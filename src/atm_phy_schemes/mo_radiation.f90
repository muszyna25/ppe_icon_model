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
!! @par Copyright
!!   2002-2009 by the Deutsche Wetterdienst (DWD) and the Max-Planck-Institut
!!   for Meteorology (MPI-M).  This software is provided for non-commerical
!!   use only.  See the LICENSE and the WARRANTY conditions
!!
!! @par License
!!   The use of ICON is hereby granted free of charge for an unlimited time,
!!   provided:
!!   <ol>
!!    <li> Its use is limited to own non-commercial and non-violent purposes;
!!    <li> The code is not re-distributed without the consent of DWD and MPI-M;
!!    <li> This header appears in all copies of the code;
!!    <li> You accept the warranty conditions (see WARRANTY).
!!   </ol>
!!   Commericial use of the code is allowed subject to a separate licensing
!!   agreement with the DWD and MPI-M
!!
!! @par Warranty
!!   This code is distributed in the home that it will be useful, but WITHOUT
!!   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!!   FITNESS FOR A PARTICULAR PURPOSE.
!
MODULE mo_radiation

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish
  USE mo_run_nml,              ONLY: ltimer

  USE mo_model_domain,         ONLY: t_patch

  USE mo_math_constants,       ONLY: pi
  USE mo_physical_constants,   ONLY: grav,  rd,    avo,   amd,  amw,  &
    &                                amco2, amch4, amn2o, amo3, amo2, &
    &                                stbo,  rcpd,  vtmpc2

  USE mo_datetime,             ONLY: rdaylen

  USE mo_radiation_nml,        ONLY: tsi,        ssi,                 &
    &                                irad_co2,   mmr_co2,             &
    &                                irad_ch4,   mmr_ch4, vpp_ch4,    &
    &                                irad_n2o,   mmr_n2o, vpp_n2o,    &
    &                                irad_o3,                         &
    &                                irad_o2,    mmr_o2,              &
    &                                irad_cfc11, vmr_cfc11,           &
    &                                irad_cfc12, vmr_cfc12,           &
!!$    &                                irad_aero,                       &
    &                                dt_rad,                          &
    &                                izenith

!!$  USE mo_greenhouse_gases,     ONLY: ghg_co2mmr, ghg_mmr_ch4, ghg_n2ommr, ghg_cfcvmr
!!$  USE mo_o3clim,               ONLY: o3clim

  USE mo_newcld_optics,        ONLY: newcld_optics
!!$  USE mo_aero_kinne,           ONLY: set_aop_kinne
!!$  USE mo_aero_volc,            ONLY: add_aop_volc

  USE mo_lrtm_par,             ONLY: jpband => nbndlw, jpxsec => maxxsec
  USE mo_lrtm,                 ONLY: lrtm

  USE mo_srtm_config,          ONLY: jpsw, jpinpx
  USE mo_srtm,                 ONLY: srtm_srtm_224gp
  USE mo_get_utc_date_tr,      ONLY: get_utc_date_tr
  USE mo_timer,                ONLY: timer_radiation, timer_start, timer_stop

!!$  USE mo_echam_phy_memory,     ONLY: prm_field

!!$  USE mo_radiation_forcing,    ONLY: prepare_forcing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pre_radiation_nwp, radiation, radheat, pre_radiation_nwp_steps

  ! --- radiative transfer parameters
  !
  REAL(wp), PARAMETER :: cemiss = 0.996_wp  !< LW Emissivity Factor
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

    !First: cases izenith==0 to izenith==2 (no date and time needed)
    IF (izenith == 0) THEN
      ! for testing: provisional setting of cos(zenith angle) and TSI
      ! The global mean insolation is TSI/4 (ca. 340 W/m2)
      DO jb = 1, pt_patch%nblks_c
        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c
        zsmu0(1:ie,jb) = 1._wp ! sun in zenith everywhere
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi/4._wp ! scale ztsi to get the correct global mean insolation
    ELSEIF(izenith == 1) THEN
      ! circular non-seasonal orbit, zenith angle dependent on latitude only,
      ! no diurnal cycle (always at 12:00 local time --> sin(time of day)=1 )
      DO jb = 1, pt_patch%nblks_c
        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c      
        zsmu0(1:ie,jb) = COS( pt_patch%cells%center(1:ie,jb)%lat )
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi/pi ! because sun is always in local noon, the TSI needs to be
      ! scaled by 1/pi to get the correct global mean insolation
    ELSEIF (izenith == 2) THEN
      ! circular non-seasonal orbit, no diurnal cycle
      ! at 07:14:15 or 16:45:45 local time (--> sin(time of day)=1/pi )
      DO jb = 1, pt_patch%nblks_c
        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c      
        zsmu0(1:ie,jb) = COS( pt_patch%cells%center(1:ie,jb)%lat )/pi
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi
    ELSEIF (izenith == 3) THEN  !Second: case izenith==3 (time (but no date) needed)

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
            IF ( z_cosmu0(jc,jb) > 1.e-6_wp ) THEN
              zsmu0(jc,jb) = zsmu0(jc,jb) + z_cosmu0(jc,jb)
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
            zsmu0(jc,jb) = zsmu0(jc,jb) / REAL(n_cosmu0pos(jc,jb),wp)
          ELSE
            zsmu0(jc,jb) = cosmu0_dark
          ENDIF
        ENDDO

      ENDDO !jb

      IF (PRESENT(zsct)) zsct = tsi

      !Third: case izenith=4 (time and date needed)
    ELSEIF (izenith == 4) THEN

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
            zsct_save = zsocof*tsi
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
            IF ( z_cosmu0(jc,jb) > 1.e-6_wp ) THEN
              zsmu0(jc,jb) = zsmu0(jc,jb) + z_cosmu0(jc,jb)
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
            zsmu0(jc,jb) = zsmu0(jc,jb) / REAL(n_cosmu0pos(jc,jb),wp)
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
      IF (PRESENT(zsct)) zsct = tsi/4._wp ! scale ztsi to get the correct global mean insolation
      RETURN
    ELSEIF(izenith == 1) THEN
      ! circular non-seasonal orbit, zenith angle dependent on latitude only,
      ! no diurnal cycle (always at 12:00 local time --> sin(time of day)=1 )
      DO jb = 1, pt_patch%nblks_c
        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c      
        zsmu0(1:ie,jb) = COS( pt_patch%cells%center(1:ie,jb)%lat )
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi/pi ! because sun is always in local noon, the TSI needs to be
                                       ! scaled by 1/pi to get the correct global mean insolation
      RETURN
    ELSEIF (izenith == 2) THEN
      ! circular non-seasonal orbit, no diurnal cycle
      ! at 07:14:15 or 16:45:45 local time (--> sin(time of day)=1/pi )
      DO jb = 1, pt_patch%nblks_c
        IF (jb == pt_patch%nblks_c) ie = pt_patch%npromz_c      
        zsmu0(1:ie,jb) = COS( pt_patch%cells%center(1:ie,jb)%lat )/pi
      ENDDO
      IF (PRESENT(zsct)) zsct = tsi
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
      IF (PRESENT(zsct)) zsct = tsi

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
          zsct_save = zsocof*tsi
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
      
    ENDIF

  END SUBROUTINE pre_radiation_nwp

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
  SUBROUTINE radiation(                                                    &
    ! input
!!$    &  jg, jb,                                                             &
    &  jce               ,kbdim           ,klev             ,klevp1        &
    & ,ktype             ,zland           ,zglac            ,cos_mu0       &
    & ,alb_vis_dir       ,alb_nir_dir     ,alb_vis_dif      ,alb_nir_dif   &
    & ,tk_sfc            ,pp_hl           ,pp_fl            ,tk_fl         &
    & ,qm_vap            ,qm_liq          ,qm_ice                          &
    & ,qm_o3                                                               &
!!$    & ,pgeom1                                                              &
    & ,cdnc              ,cld_frc                                          &
    ! output
    & ,cld_cvr                                                             &
!!$    & ,nir_sfc           ,nir_dff_sfc     ,vis_sfc          ,vis_dff_sfc   &
!!$    & ,dpar_sfc          ,par_dff_sfc                                      &
    & ,emter_clr         ,trsol_clr       ,emter_all        ,trsol_all     &
    & ,opt_halo_cosmu0  )

    ! input
    ! -----
    !
    INTEGER, INTENT(in)   :: &
!!$      &  jg,                 & !< domain index
!!$      &  jb,                 & !< block index
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
      &  tk_sfc(kbdim),      & !< Surface temperature
      &  pp_hl(kbdim,klevp1),& !< pressure at half levels [Pa]
      &  pp_fl(kbdim,klev),  & !< Pressure at full levels [Pa]
      &  tk_fl(kbdim,klev),  & !< Temperature on full levels [K]
      &  qm_vap(kbdim,klev), & !< Water vapor mixing ratio
      &  qm_liq(kbdim,klev), & !< Liquid water mixing ratio
      &  qm_ice(kbdim,klev), & !< Ice water mixing ratio
!!$      &  pgeom1(:,:),        & !< geopotential above ground
      &  cdnc(kbdim,klev),   & !< Cloud drop number concentration
      &  cld_frc(kbdim,klev)   !< Cloud fraction
    
    LOGICAL, INTENT(in), OPTIONAL :: opt_halo_cosmu0    

    ! output
    ! ------
    !
    REAL(wp), INTENT(out) ::      &
      &  cld_cvr(kbdim),          & !< Cloud cover in a column
!!$      &  nir_sfc(:),           & !< SW Near Infrared
!!$      &  nir_dff_sfc(:),       & !< Near-Infrared diffuse fraction of SW NIR
!!$      &  vis_sfc(:),           & !< SW Visible (250-680 nm)
!!$      &  vis_dff_sfc(:),       & !< Near-Infrared diffuse fraction of SW vis
!!$      &  dpar_sfc(:),          & !< surf. PAR downw.
!!$      &  par_dff_sfc(:),       & !< fraction of diffuse PAR
      &  emter_clr(kbdim,klevp1), & !< Surface terrestrial emissivity
      &  trsol_clr(kbdim,klevp1), & !< Surface solar transmissivity
      &  emter_all(kbdim,klevp1), & !< Terrestrial emissivity
      &  trsol_all(kbdim,klevp1)    !< Solar transmissivity


    INTEGER  :: jk, jl

    REAL(wp) ::                     &
    &    cos_mu0_halo,              & !< cos(zenith angle) value delimiting the halo
    &    cos_mu0_mod(kbdim)           !< modified cos(zenith angle)

    REAL(wp) ::                     &
      &  ppd_hl(kbdim,klev),        & !< pressure diff between half levels [Pa]
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
!!$      &  toa_solar_irr(kbdim),      & !< TOA solar irradiance
      &  flx_uplw_sfc(kbdim),       & !< Srfc upward lw flux  [Wm2]
      &  flx_upsw_sfc(kbdim),       & !< Srfc upward sw flux  [Wm2]
      &  flx_uplw_clr_sfc(kbdim),   & !< Srfc upward lw flux (clear sky) [Wm2]
      &  flx_upsw_clr_sfc(kbdim),   & !< Srfc upward sw flux (clear sky) [Wm2]
      &  flx_dnlw(kbdim,klevp1),    & !< Net dwnwrd LW flux [Wm2]
      &  flx_dnsw(kbdim,klevp1),    & !< Net dwnwrd SW flux [Wm2]
      &  flx_dnlw_clr(kbdim,klevp1),& !< Net dn LW flux (clear sky) [Wm2]
      &  flx_dnsw_clr(kbdim,klevp1)   !< Net dn SW flux (clear sky) [Wm2]
!!$      &  flx_dnsw_clr(kbdim,klevp1),& !< Net dn SW flux (clear sky) [Wm2]
!!$      &  dvis_sfc(kbdim)              !< surface downward visible flux

!!$    ! for debugging
!!$    REAL(wp) :: rloland(kbdim) ,rloglac(kbdim)

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
    ppd_hl(1:jce,1:klev) = pp_hl(1:jce,2:klev+1)-pp_hl(1:jce,1:klev)
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
    cld_frc_sec(1:jce,:) = MAX(cld_frc(1:jce,:),0.0_wp)
    !
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
    xm_cfc11 (1:jce,:) = gas_profile(jce, klev, irad_cfc11,  &
      &                              mmr_gas = vmr_cfc11 )
    xm_cfc12 (1:jce,:) = gas_profile(jce, klev, irad_cfc12,  &
      &                              mmr_gas = vmr_cfc12 )
    !
    ozon: SELECT CASE (irad_o3)
    CASE (0)
      qm_o3(1:jce,:) = 0.0_wp
    CASE (6)
!!$    CASE (4)
!!$      xm_o3(1:jce,:) = o3clim(krow,jce,kbdim,klev,pp_hl,pp_fl)
    CASE default
      CALL finish('radiation','o3: this "irad_o3" is not supported')
    END SELECT ozon

!!$    ! debug rrtmin1
!!$    ! -------------
!!$    DO jl = 1, jce
!!$      IF (loland(jl)) THEN
!!$        rloland(jl) = 1.0_wp
!!$      ELSE
!!$        rloland(jl) = 0.0_wp
!!$      END IF
!!$      IF (loglac(jl)) THEN
!!$        rloglac(jl) = 1.0_wp
!!$      ELSE
!!$        rloglac(jl) = 0.0_wp
!!$      END IF
!!$    END DO
!!$    prm_field(jg)% debug_2d_1(:,jb) = REAL(ktype(:),wp)
!!$    prm_field(jg)% debug_2d_2(:,jb) = rloland(:)
!!$    prm_field(jg)% debug_2d_3(:,jb) = rloglac(:)
!!$    prm_field(jg)% debug_2d_4(:,jb) = cos_mu0_mod(:)
!!$    prm_field(jg)% debug_2d_5(:,jb) = alb_vis_dir(:)
!!$    prm_field(jg)% debug_2d_6(:,jb) = alb_nir_dir(:)
!!$    prm_field(jg)% debug_2d_7(:,jb) = alb_vis_dif(:)
!!$    prm_field(jg)% debug_2d_8(:,jb) = alb_nir_dif(:)
!!$    !
!!$    prm_field(jg)% debug_3d_1(:,1:klev,jb) = pp_fl(:,1:klev)
!!$    prm_field(jg)% debug_3d_2(:,1:klev,jb) = pp_hl(:,1:klev)
!!$    prm_field(jg)% debug_3d_3(:,1:klev,jb) = tk_fl(:,1:klev)
!!$    prm_field(jg)% debug_3d_4(:,1:klev,jb) = tk_hl(:,1:klev)
!!$    prm_field(jg)% debug_3d_5(:,1:klev,jb) = xq_vap(:,1:klev)
!!$    prm_field(jg)% debug_3d_6(:,1:klev,jb) = xq_liq(:,1:klev)
!!$    prm_field(jg)% debug_3d_7(:,1:klev,jb) = xq_ice(:,1:klev)
!!$    prm_field(jg)% debug_3d_8(:,1:klev,jb) = cdnc(:,1:klev)
    !
!!$    ! debug rrtmin2
!!$    ! -------------
!!$    prm_field(jg)% debug_2d_1(:,jb) = pp_sfc(:)
!!$    prm_field(jg)% debug_2d_2(:,jb) = tk_sfc(:)
!!$    !
!!$    prm_field(jg)% debug_3d_1(:,1:klev,jb) = xm_o3(:,1:klev)
!!$    prm_field(jg)% debug_3d_2(:,1:klev,jb) = xm_co2(:,1:klev)
!!$    prm_field(jg)% debug_3d_3(:,1:klev,jb) = xm_ch4(:,1:klev)
!!$    prm_field(jg)% debug_3d_4(:,1:klev,jb) = xm_n2o(:,1:klev)
!!$    prm_field(jg)% debug_3d_5(:,1:klev,jb) = xm_cfc11(:,1:klev)
!!$    prm_field(jg)% debug_3d_6(:,1:klev,jb) = xm_cfc12(:,1:klev)
!!$    prm_field(jg)% debug_3d_7(:,1:klev,jb) = xm_o2(:,1:klev)
    !
    ! 2.0 Call interface to radiation solver
    ! --------------------------------------
    !
    CALL rrtm_interface(                                                    &
      ! input
!!$      & irad_aero                                                          ,&
!!$      & jg              ,jb                                                ,&
      & jce             ,kbdim           ,klev                             ,&
      & ktype           ,zland           ,zglac                            ,&
      & cos_mu0_mod                                                        ,&
!!$      & pgeom1                                                             ,&
      & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
      & pp_fl           ,pp_hl           ,pp_sfc          ,tk_fl           ,&
      & tk_hl           ,tk_sfc          ,xq_vap                           ,&
      & xq_liq          ,xq_ice                                            ,&
      & cdnc                                                               ,&
      & cld_frc_sec                                                        ,&
      & qm_o3           ,xm_co2          ,xm_ch4                           ,&
      & xm_n2o          ,xm_cfc11        ,xm_cfc12        ,xm_o2           ,&
      ! output
      & flx_dnlw        ,flx_dnsw        ,flx_dnlw_clr    ,flx_dnsw_clr    ,&
      & flx_uplw_sfc    ,flx_upsw_sfc    ,flx_uplw_clr_sfc,flx_upsw_clr_sfc)
!!$      & flx_uplw_sfc    ,flx_upsw_sfc    ,flx_uplw_clr_sfc,flx_upsw_clr_sfc,&
!!$      & toa_solar_irr   ,nir_sfc         ,nir_dff_sfc                      ,&
!!$      & vis_sfc         ,vis_dff_sfc     ,dvis_sfc        ,dpar_sfc        ,&
!!$      & par_dff_sfc                                                         )

!!$    ! debug rrtmout
!!$    ! -------------
!!$    prm_field(jg)% debug_2d_1(:,jb) = flx_uplw_sfc(:)
!!$    prm_field(jg)% debug_2d_2(:,jb) = flx_upsw_sfc(:)
!!$    prm_field(jg)% debug_2d_3(:,jb) = flx_uplw_clr_sfc(:)
!!$    prm_field(jg)% debug_2d_4(:,jb) = flx_upsw_clr_sfc(:)
!!$    !
!!$    prm_field(jg)% debug_3d_1(:,1:klev,jb) = flx_dnlw(:,1:klev)
!!$    prm_field(jg)% debug_3d_2(:,1:klev,jb) = flx_dnsw(:,1:klev)
!!$    prm_field(jg)% debug_3d_3(:,1:klev,jb) = flx_dnlw_clr(:,1:klev)
!!$    prm_field(jg)% debug_3d_4(:,1:klev,jb) = flx_dnsw_clr(:,1:klev)

    !
    ! 3.0 Diagnostics
    ! ---------------
    !
    ! --- Total fluxes
    emter_all (1:jce,1:klevp1) = flx_dnlw (1:jce,1:klevp1)
    trsol_all (1:jce,1:klevp1) = flx_dnsw (1:jce,1:klevp1)                           &
      &                          / SPREAD(cos_mu0_mod(1:jce)*tsi,2,klevp1)
    !
!!$    ! --- fluxes for JSBACH
!!$    nir_sfc  (1:jce) = nir_sfc  (1:jce) / cos_mu0_mod (1:jce)
!!$    vis_sfc  (1:jce) = vis_sfc  (1:jce) / cos_mu0_mod (1:jce)
!!$    dpar_sfc (1:jce) = dpar_sfc (1:jce) / cos_mu0_mod (1:jce)
    !
    ! --- Clear sky fluxes
    emter_clr(1:jce,1:klevp1) = flx_dnlw_clr(1:jce,1:klevp1)
    trsol_clr(1:jce,1:klevp1) = flx_dnsw_clr(1:jce,1:klevp1)                        &
      &                         / SPREAD(cos_mu0_mod(1:jce)*tsi,2,klevp1)

    IF (ltimer) CALL timer_stop(timer_radiation)

  END SUBROUTINE radiation
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
        gas_profile(1:jce,:) = MAX(gas_val(1:jce,:), 0.0_wp)
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
        gas_profile(1:jce,:) = gas_scenario
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
!!$    & irad_aero                                                          ,&
!!$    & jg              ,jb                                                ,&
    & jce             ,kbdim           ,klev                             ,&
    & ktype           ,zland           ,zglac                            ,&
    & pmu0                                                               ,&
!!$    & pgeom1                                                             ,&
    & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
    & pp_fl           ,pp_hl           ,pp_sfc          ,tk_fl           ,&
    & tk_hl           ,tk_sfc          ,xm_vap                           ,&
    & xm_liq          ,xm_ice                                            ,&
    & cdnc                                                               ,&
    & cld_frc                                                            ,&
    & xm_o3           ,xm_co2          ,xm_ch4                           ,&
    & xm_n2o          ,xm_cfc11        ,xm_cfc12        ,xm_o2           ,&
    ! output
    & flx_lw_net      ,flx_sw_net      ,flx_lw_net_clr  ,flx_sw_net_clr  ,&
    & flx_uplw_sfc    ,flx_upsw_sfc    ,flx_uplw_sfc_clr,flx_upsw_sfc_clr)
!!$    & flx_uplw_sfc    ,flx_upsw_sfc    ,flx_uplw_sfc_clr,flx_upsw_sfc_clr,&
!!$    & sw_irr_toa      ,nir_sfc         ,nir_dff_sfc                      ,&
!!$    & vis_sfc         ,vis_dff_sfc     ,dvis_sfc        ,dpar_sfc        ,&
!!$    & par_dff_sfc                                                         )

    INTEGER,INTENT(in)  ::                &
!!$      &  irad_aero,                       & !< aerosol control
!!$      &  jg,                              & !< domain index
!!$      &  jb,                              & !< block index
      &  jce,                             & !< number of columns
      &  kbdim,                           & !< first dimension of 2-d arrays
      &  klev                               !< number of levels

    INTEGER,INTENT(in)  ::                &
      &  ktype(kbdim)                       !< type of convection

    REAL(wp),INTENT(in) ::                &
      &  zland(kbdim),                    & !< land-sea mask. (1. = land, 0. = sea/lakes)
      &  zglac(kbdim),                    & !< fraction of land covered by glaciers
      &  pmu0(kbdim),                     & !< mu0 for solar zenith angle
!!$      &  pgeom1(kbdim,klev),              & !< geopotential above ground
      &  alb_vis_dir(kbdim),              & !< surface albedo for vis range and dir light
      &  alb_nir_dir(kbdim),              & !< surface albedo for NIR range and dir light
      &  alb_vis_dif(kbdim),              & !< surface albedo for vis range and dif light
      &  alb_nir_dif(kbdim),              & !< surface albedo for NIR range and dif light
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
      &  xm_o2(kbdim,klev)                  !< o2  mass mixing ratio

    REAL(wp), INTENT(out) ::              &
      &  flx_lw_net(kbdim,klev+1),        & !< net downward LW flux profile,
      &  flx_sw_net(kbdim,klev+1),        & !< net downward SW flux profile,
      &  flx_lw_net_clr(kbdim,klev+1),    & !< clrsky downward LW flux profile,
      &  flx_sw_net_clr(kbdim,klev+1),    & !< clrsky downward SW flux profile,
      &  flx_uplw_sfc(kbdim),             & !< sfc LW upward flux,
      &  flx_upsw_sfc(kbdim),             & !< sfc SW upward flux,
      &  flx_uplw_sfc_clr(kbdim),         & !< clrsky sfc LW upward flux,
      &  flx_upsw_sfc_clr(kbdim)            !< clrsky sfc SW upward flux,
!!$      &  flx_upsw_sfc_clr(kbdim),         & !< clrsky sfc SW upward flux,
!!$      &  sw_irr_toa(kbdim),               & !< top of atmosphere solar irradiation
!!$      &  nir_sfc(kbdim),                  & !< net surface NIR flux
!!$      &  nir_dff_sfc(kbdim),              & !< fraction of diffuse NIR
!!$      &  vis_sfc(kbdim),                  & !< net surface visible flux
!!$      &  vis_dff_sfc(kbdim),              & !< fraction of diffuse visible
!!$      &  dvis_sfc(kbdim),                 & !< surf. visible downward
!!$      &  dpar_sfc(kbdim),                 & !< surf. PAR downw.
!!$      &  par_dff_sfc(kbdim)                 !< fraction of diffuse PAR

    INTEGER  :: jk, jl, jp, jkb,          & !< loop indicies
      &  icldlyr(kbdim,klev)                !< index for clear or cloudy

    REAL(wp) ::                           &
      &  zsemiss(kbdim,jpband),           & !< LW surface emissivity by band
!!$      &  ppd_hl(kbdim,klev),              & !< pressure thickness in Pa
      &  pm_sfc(kbdim),                   & !< surface pressure in hPa
!!$      &  dnir_sfc(kbdim),                 & !< surf. NIR downw.
!!$      &  unir_toa(kbdim),                 & !< top NIR upw.
!!$      &  uvis_toa(kbdim),                 & !< top vis. upw.
!!$      &  dnir_sfc_cld(kbdim),             & !< surf. NIR downw. cloudy sky
!!$      &  dvis_sfc_cld(kbdim),             & !< surf. vis. downw. cloudy sky
!!$      &  unir_toa_cld(kbdim),             & !< top NIR upw. cloudy sky
!!$      &  uvis_toa_cld(kbdim),             & !< top vis. upw. cloudy sky
!!$      &  zsudu(kbdim),                    & !< direct beam SW
      &  amm,                             & !< molecular weight of moist air
      &  delta,                           & !< pressure thickness
      &  zscratch                           !< scratch array
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

!!$    ! for debugging
!!$    REAL(wp) :: rloland(kbdim) ,rloglac(kbdim)

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
    zsemiss(1:jce,:) = cemiss
    !
    ! 3.0 Particulate Optical Properties
    ! --------------------------------
!!$    ppd_hl(1:jce,:) = pp_hl(1:jce,2:klev+1)-pp_hl(1:jce,1:klev)
!!$
!!$    SELECT CASE (irad_aero)
!!$    CASE (0)
      aer_tau_lw_vr(:,:,:) = 0.0_wp
      aer_tau_sw_vr(:,:,:) = 0.0_wp
      aer_piz_sw_vr(:,:,:) = 1.0_wp
      aer_cg_sw_vr(:,:,:)  = 0.0_wp
!!$    CASE (3)
!!$      CALL set_aop_kinne( &
!!$        & jce              ,kbdim                 ,klev             ,&
!!$        & krow             ,jpband                ,jpsw             ,&
!!$        & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
!!$        & aer_cg_sw_vr     ,ppd_hl                ,pp_fl            ,&
!!$        & tk_fl            ,pgeom1 )
!!$    CASE (5)
!!$      CALL set_aop_kinne( &
!!$        & jce              ,kbdim                 ,klev             ,&
!!$        & krow             ,jpband                ,jpsw             ,&
!!$        & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
!!$        & aer_cg_sw_vr     ,ppd_hl                ,pp_fl            ,&
!!$        & tk_fl            ,pgeom1 )
!!$      CALL add_aop_volc( &
!!$        & jce              ,kbdim                 ,klev             ,&
!!$        & krow             ,jpband                ,jpsw             ,&
!!$        & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
!!$        & aer_cg_sw_vr     ,ppd_hl                ,pp_fl            ,&
!!$        & tk_fl )
!!$    CASE DEFAULT
!!$    END SELECT

!!$    ! debug newcldin
!!$    ! --------------
!!$    DO jl = 1, jce
!!$      IF (laland(jl)) THEN
!!$        rloland(jl) = 1.0_wp
!!$      ELSE
!!$        rloland(jl) = 0.0_wp
!!$      END IF
!!$      IF (laglac(jl)) THEN
!!$        rloglac(jl) = 1.0_wp
!!$      ELSE
!!$        rloglac(jl) = 0.0_wp
!!$      END IF
!!$    END DO
!!$    prm_field(jg)% debug_2d_1(:,jb) = rloglac(:)
!!$    prm_field(jg)% debug_2d_2(:,jb) = rloland(:)
!!$    prm_field(jg)% debug_2d_3(:,jb) = REAL(ktype(:),wp)
!!$    !
!!$    prm_field(jg)% debug_3d_1(:,1:klev,jb) = REAL(icldlyr(:,1:klev),wp)
!!$    prm_field(jg)% debug_3d_2(:,1:klev,jb) = zlwp_vr(:,1:klev)
!!$    prm_field(jg)% debug_3d_3(:,1:klev,jb) = ziwp_vr(:,1:klev)
!!$    prm_field(jg)% debug_3d_4(:,1:klev,jb) = zlwc_vr(:,1:klev)
!!$    prm_field(jg)% debug_3d_5(:,1:klev,jb) = ziwc_vr(:,1:klev)
!!$    prm_field(jg)% debug_3d_6(:,1:klev,jb) = cdnc_vr(:,1:klev)

      CALL newcld_optics(                                                       &
        & jce          ,kbdim        ,klev         ,jpband       ,jpsw         ,&
        & zglac        ,zland        ,ktype        ,icldlyr                    ,&
        & zlwp_vr      ,ziwp_vr      ,zlwc_vr      ,ziwc_vr      ,cdnc_vr      ,&
        & cld_tau_lw_vr,cld_tau_sw_vr,cld_piz_sw_vr,cld_cg_sw_vr                )

!!$    ! debug newcldout
!!$    ! ---------------
!!$    ! for band 5
!!$    prm_field(jg)% debug_3d_1(:,1:klev,jb) = cld_tau_lw_vr(:,1:klev,5)
!!$    prm_field(jg)% debug_3d_2(:,1:klev,jb) = cld_tau_sw_vr(:,5,1:klev)
!!$    prm_field(jg)% debug_3d_3(:,1:klev,jb) = cld_piz_sw_vr(:,5,1:klev)
!!$    prm_field(jg)% debug_3d_4(:,1:klev,jb) = cld_cg_sw_vr (:,5,1:klev)

!!$    ! debug lrtmin1
!!$    ! -------------
!!$    prm_field(jg)% debug_2d_1(:,jb) = pm_sfc(:)
!!$    prm_field(jg)% debug_2d_2(:,jb) = tk_sfc(:)
!!$    !
!!$    prm_field(jg)% debug_3d_1(:,1:klev,jb) = pm_fl_vr(:,1:klev)
!!$    prm_field(jg)% debug_3d_2(:,1:klev,jb) = tk_fl_vr(:,1:klev)
!!$    prm_field(jg)% debug_3d_3(:,1:klev,jb) = tk_hl_vr(:,1:klev)
!!$    prm_field(jg)% debug_3d_4(:,1:klev,jb) = wx_vr(:,2,1:klev)
!!$    prm_field(jg)% debug_3d_5(:,1:klev,jb) = wx_vr(:,3,1:klev)
!!$    prm_field(jg)% debug_3d_6(:,1:klev,jb) = cld_frc_vr(:,1:klev)
!!$    prm_field(jg)% debug_3d_7(:,1:klev,jb) = cld_tau_lw_vr(:,1:klev,5) ! band 5
!!$    prm_field(jg)% debug_3d_8(:,1:klev,jb) = aer_tau_lw_vr(:,1:klev,5) ! band 5

!!$    ! debug lrtmin2
!!$    ! -------------
!!$    prm_field(jg)% debug_3d_1(:,1:klev,jb) = wkl_vr(:,1,1:klev)   ! H2O
!!$    prm_field(jg)% debug_3d_2(:,1:klev,jb) = wkl_vr(:,2,1:klev)   ! CO2
!!$    prm_field(jg)% debug_3d_3(:,1:klev,jb) = wkl_vr(:,3,1:klev)   ! O3
!!$    prm_field(jg)% debug_3d_4(:,1:klev,jb) = wkl_vr(:,4,1:klev)   ! N2O
!!$    prm_field(jg)% debug_3d_5(:,1:klev,jb) = wkl_vr(:,5,1:klev)   ! 0.0
!!$    prm_field(jg)% debug_3d_6(:,1:klev,jb) = wkl_vr(:,6,1:klev)   ! CH4
!!$    prm_field(jg)% debug_3d_7(:,1:klev,jb) = wkl_vr(:,7,1:klev)   ! O2
!!$    prm_field(jg)% debug_3d_8(:,1:klev,jb) = col_dry_vr(:,1:klev) ! air

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

!!$    ! debug lrtmout
!!$    ! -------------
!!$    prm_field(jg)% debug_3d_1(:,1:klev,jb) = flx_uplw_vr(:,1:klev)
!!$    prm_field(jg)% debug_3d_2(:,1:klev,jb) = flx_dnlw_vr(:,1:klev)
!!$    prm_field(jg)% debug_3d_3(:,1:klev,jb) = flx_uplw_clr_vr(:,1:klev)
!!$    prm_field(jg)% debug_3d_4(:,1:klev,jb) = flx_dnlw_clr_vr(:,1:klev)

!!$    ! debug srtmin (without fields checked in lrtmin and lrtmin2)
!!$    ! ------------
!!$    prm_field(jg)% debug_2d_1(:,jb) = alb_vis_dir(:)
!!$    prm_field(jg)% debug_2d_2(:,jb) = alb_nir_dir(:)
!!$    prm_field(jg)% debug_2d_3(:,jb) = alb_vis_dif(:)
!!$    prm_field(jg)% debug_2d_4(:,jb) = alb_nir_dif(:)
!!$    prm_field(jg)% debug_2d_5(:,jb) = pmu0(:)
!!$    !
!!$    ! band 5
!!$    prm_field(jg)% debug_3d_2(:,1:klev,jb) = cld_tau_sw_vr(:,5,1:klev)
!!$    prm_field(jg)% debug_3d_3(:,1:klev,jb) = cld_cg_sw_vr (:,5,1:klev)
!!$    prm_field(jg)% debug_3d_4(:,1:klev,jb) = cld_piz_sw_vr(:,5,1:klev)
!!$    prm_field(jg)% debug_3d_5(:,1:klev,jb) = aer_tau_sw_vr(:,1:klev,5)
!!$    prm_field(jg)% debug_3d_6(:,1:klev,jb) = aer_cg_sw_vr (:,1:klev,5)
!!$    prm_field(jg)% debug_3d_7(:,1:klev,jb) = aer_piz_sw_vr(:,1:klev,5)
!!$    ! 14 reals on 14 levels
!!$    prm_field(jg)% debug_3d_8(:,1:14,jb) = SPREAD(ssi(1:14),1,kbdim)

    CALL srtm_srtm_224gp(                                                     &
      !    input
      &    jce             ,kbdim           ,klev            ,jpsw           ,&
      &    alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif    ,&
      &    pm_fl_vr        ,tk_fl_vr        ,pmu0                            ,&
      &    col_dry_vr      ,wkl_vr                                           ,&
      &    cld_frc_vr      ,cld_tau_sw_vr   ,cld_cg_sw_vr    ,cld_piz_sw_vr  ,&
      &    aer_tau_sw_vr   ,aer_cg_sw_vr    ,aer_piz_sw_vr                   ,&
      &    ssi                                                               ,&
      !    output
      &    flx_dnsw        ,flx_upsw        ,flx_dnsw_clr    ,flx_upsw_clr)
!!$      &    flx_dnsw        ,flx_upsw        ,flx_dnsw_clr    ,flx_upsw_clr   ,&
!!$      &    dnir_sfc        ,dvis_sfc        ,unir_toa        ,uvis_toa       ,&
!!$      &    dnir_sfc_cld    ,dvis_sfc_cld    ,unir_toa_cld    ,uvis_toa_cld   ,&
!!$      &    zsudu           ,nir_sfc         ,nir_dff_sfc     ,vis_sfc        ,&
!!$      &    vis_dff_sfc     ,dpar_sfc        ,par_dff_sfc)

!!$    ! debug srtmout
!!$    ! -------------
!!$    prm_field(jg)% debug_3d_1(:,1:klev,jb) = flx_dnsw(:,1:klev)
!!$    prm_field(jg)% debug_3d_2(:,1:klev,jb) = flx_upsw(:,1:klev)
!!$    prm_field(jg)% debug_3d_3(:,1:klev,jb) = flx_dnsw_clr(:,1:klev)
!!$    prm_field(jg)% debug_3d_4(:,1:klev,jb) = flx_upsw_clr(:,1:klev)

    !
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
    &                 pmair         ,  &
    &                 pq            ,  &
    &                 pi0           ,  &
    &                 ptsfc         ,  &
    &                 ptsfctrad     ,  &
    &                 ptrmsw        ,  &
    &                 pflxlw        ,  &
    &                 pdtdtradsw    ,  &
    &                 pdtdtradlw    ,  &
    &                 pflxsfcsw     ,  &
    &                 pflxsfclw     ,  &
    &                 pflxtoasw )

    INTEGER,  INTENT(in)  ::    &
      &     jcs, jce, kbdim,    &
      &     klev,   klevp1

    REAL(wp), INTENT(in)  ::         &
      &     pmair      (kbdim,klev), & ! mass of air in layer                     [kg/m2]
      &     pq         (kbdim,klev), & ! specific humidity at t-dt                [kg/kg]
      &     pi0        (kbdim),      & ! local solar incoming flux at TOA         [W/m2]
      &     ptsfc      (kbdim),      & ! surface temperature at t                 [K]
      &     ptsfctrad  (kbdim),      & ! surface temperature at trad              [K]
      &     ptrmsw     (kbdim,klevp1), & ! shortwave transmissivity at trad         []
      &     pflxlw     (kbdim,klevp1)    ! longwave net flux at trad                [W/m2]

    REAL(wp), INTENT(out) ::         &
      &     pdtdtradsw (kbdim,klev), & ! shortwave temperature tendency           [K/s]
      &     pdtdtradlw (kbdim,klev)    ! longwave temperature tendency            [K/s]

    REAL(wp), INTENT(inout), OPTIONAL :: &
      &     pflxsfcsw (kbdim), &       ! shortwave surface net flux [W/m2]
      &     pflxsfclw (kbdim), &       ! longwave surface net flux [W/m2]
      &     pflxtoasw (kbdim)          ! shortwave toa net flux [W/m2]


    
    ! Local arrays
    REAL(wp) ::                    &
      &     zflxsw (kbdim,klevp1), &
      &     zflxlw (kbdim,klevp1), &
      &     zconv  (kbdim,klev)

    ! Conversion factor for heating rates:
    zconv(jcs:jce,1:klev) = rcpd/(pmair(jcs:jce,1:klev)*(1._wp+vtmpc2*pq(jcs:jce,1:klev)))

    ! Shortwave fluxes = transmissivity * loacal solar incoming flux at TOA
    ! ----------------
    ! - TOA
    zflxsw(jcs:jce,1)      = ptrmsw(jcs:jce,1)      *        pi0(jcs:jce)
    ! - Atmosphere
    zflxsw(jcs:jce,2:klev) = ptrmsw(jcs:jce,2:klev) * SPREAD(pi0(jcs:jce),2,klev-1)
    ! - Surface
    zflxsw(jcs:jce,klevp1) = ptrmsw(jcs:jce,klevp1) *        pi0(jcs:jce)

    ! Longwave fluxes
    ! - TOA
    zflxlw(jcs:jce,1)      = pflxlw(jcs:jce,1)
    ! - Atmosphere
    zflxlw(jcs:jce,2:klev) = pflxlw(jcs:jce,2:klev)
    ! - Surface
    !   Adjust for changed surface temperature (ptsfc) with respect to the
    !   surface temperature used for the longwave flux computation (ptsfctrad).
    !   --> modifies heating in lowermost layer only (is this smart?)
    zflxlw(jcs:jce,klevp1) = pflxlw(jcs:jce,klevp1)               &
      &                   + cemiss*stbo * ptsfctrad(jcs:jce)**4 &
      &                   - cemiss*stbo * ptsfc    (jcs:jce)**4
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
    
  END SUBROUTINE radheat

END MODULE mo_radiation
