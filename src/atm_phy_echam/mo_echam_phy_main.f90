!>
!! @brief Subroutine echam_phy_main calls all the parameterization schemes
!!
!! @author Hui Wan, MPI-M
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version from ECHAM6 (revision 2028)
!!  Modified for ICOHAM by Hui Wan and Marco Giorgetta (2010)
!!  Modified for ICONAM by Marco Giorgetta (2014)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
@PROCESS SPILLSIZE(5000)
#endif
!OCL NOALIAS

MODULE mo_echam_phy_main

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_math_constants,      ONLY: pi
  USE mo_physical_constants,  ONLY: grav, cpd, cpv, cvd, cvv
  USE mo_impl_constants,      ONLY: inh_atmosphere
  USE mo_run_config,          ONLY: ntracer, nlev, nlevm1, nlevp1,    &
    &                               iqv, iqc, iqi, iqt
  USE mo_dynamics_config,     ONLY: iequations
  USE mo_ext_data_state,      ONLY: ext_data
  USE mo_echam_phy_config,    ONLY: phy_config => echam_phy_config
  USE mo_echam_conv_config,   ONLY: echam_conv_config
  USE mo_echam_cloud_config,  ONLY: echam_cloud_config
  USE mo_cumastr,             ONLY: cucall
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field, prm_field,     &
    &                               t_echam_phy_tend,  prm_tend
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop,                &
    &                               timer_cover, timer_radiation, timer_radheat,    &
    &                               timer_vdiff_down, timer_surface,timer_vdiff_up, &
    &                               timer_gw_hines, timer_ssodrag,                  &
    &                               timer_cucall, timer_cloud
  USE mo_datetime,            ONLY: t_datetime
  USE mo_ham_aerosol_params,  ONLY: ncdnc, nicnc
  USE mo_echam_sfc_indices,   ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_surface,             ONLY: update_surface
  USE mo_surface_diag,        ONLY: nsurf_diag
  USE mo_cloud,               ONLY: cloud
  USE mo_cover,               ONLY: cover
  USE mo_radheating,          ONLY: radheating
  USE mo_psrad_radiation,     ONLY: psrad_radiation
  USE mo_psrad_radiation_parameters, ONLY: psctm
  USE mo_radiation_config,    ONLY: tsi, izenith, irad_o3
  USE mo_vdiff_config,        ONLY: vdiff_config
  USE mo_vdiff_downward_sweep,ONLY: vdiff_down
  USE mo_vdiff_upward_sweep,  ONLY: vdiff_up
  USE mo_vdiff_solver,        ONLY: nvar_vdiff, nmatrix, imh, imqv,   &
                                  & ih_vdiff=>ih, iqv_vdiff=>iqv
  USE mo_gw_hines,            ONLY: gw_hines
  USE mo_ssortns,             ONLY: ssodrag
  ! provisional to get coordinates
  USE mo_model_domain,        ONLY: p_patch
  USE mo_util_dbg_prnt,      ONLY: dbg_print

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: echam_phy_main

CONTAINS
  !>
  !!
  SUBROUTINE echam_phy_main( jg,jb,jcs,jce,nbdim,      &
    &                        datetime,pdtime,psteplen, &
    &                        ltrig_rad                 )

    INTEGER         ,INTENT(IN) :: jg             !< grid level/domain index
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block

    TYPE(t_datetime),INTENT(IN) :: datetime       !< time step
    REAL(wp)        ,INTENT(IN) :: pdtime         !< time step
    REAL(wp)        ,INTENT(IN) :: psteplen       !< 2*pdtime in case of leapfrog

    LOGICAL         ,INTENT(IN) :: ltrig_rad      !< perform radiative transfer computation

    ! Local variables

    TYPE(t_echam_phy_field),   POINTER :: field
    TYPE(t_echam_phy_tend) ,   POINTER :: tend

    REAL(wp) :: zlat_deg(nbdim)           !< latitude in deg N

!!$    REAL(wp) :: zhmixtau   (nbdim,nlev)   !< timescale of mixing for horizontal eddies
    REAL(wp) :: zvmixtau   (nbdim,nlev)   !< timescale of mixing for vertical turbulence
    REAL(wp) :: zqtvar_prod(nbdim,nlev)   !< production rate of total water variance
                                          !< due to turbulence. Computed in "vdiff",
                                          !< used by "cloud"
    INTEGER  :: itype(nbdim)              !< type of convection
    INTEGER  :: ictop (nbdim)             !< from massflux

    REAL(wp) :: zfrl (nbdim)              !< fraction of land in the grid box
    REAL(wp) :: zfrw (nbdim)              !< fraction of water (without ice) in the grid point
    REAL(wp) :: zfri (nbdim)              !< fraction of ice in the grid box
    REAL(wp) :: zfrc (nbdim,nsfc_type)    !< zfrl, zfrw, zfrc combined
    REAL(wp) :: zri_tile(nbdim,nsfc_type) !< Richardson number

    INTEGER  :: ilab   (nbdim,nlev)
!    REAL(wp) :: zcvcbot(nbdim)
!    REAL(wp) :: zwcape (nbdim)

    REAL(wp) :: zqtec  (nbdim,nlev)       !< tracer tendency due to entrainment/detrainment

    REAL(wp) :: flux_factor(nbdim)
    REAL(wp) :: ztsi                      !< total solar irradiation at 1 AU   [W/m2]
    REAL(wp) :: zcd                       !< specific heat of dry air          [J/K/kg]
    REAL(wp) :: zcv                       !< specific heat of water vapor      [J/K/kg]
    REAL(wp) :: zcair  (nbdim,nlev)       !< specific heat of moist air        [J/K/kg]
    REAL(wp) :: zcpair (nbdim,nlev)       !< specific heat of moist air at const. pressure [J/K/kg]
    REAL(wp) :: zcvair (nbdim,nlev)       !< specific heat of moist air at const. volume   [J/K/kg]
    REAL(wp) :: zconv  (nbdim,nlev)       !< conversion factor q-->dT/dt       [(K/s)/(W/m2)]
    REAL(wp) :: zdelp  (nbdim,nlev)       !< layer thickness in pressure coordinate  [Pa]

    REAL(wp) :: zq_phy (nbdim,nlev)       !< heating by whole ECHAM physics    [W/m2]
    REAL(wp) :: zq_rsw (nbdim,nlev)       !< heating by short wave radiation   [W/m2]
    REAL(wp) :: zq_rlw (nbdim,nlev)       !< heating by long  wave radiation   [W/m2]
    REAL(wp) :: zq_rlw_impl (nbdim)       !< additional heating by LW rad. due to impl. coupling in surface energy balance [W/m2]
!!$    REAL(wp) :: zq_vdf (nbdim,nlev)       !< heating by vertical diffusion     [W/m2]
    REAL(wp) :: zq_sso (nbdim,nlev)       !< heating by subgrid scale orogr.   [W/m2]
    REAL(wp) :: zq_gwh (nbdim,nlev)       !< heating by atm. gravity waves     [W/m2]
!!$    REAL(wp) :: zq_cnv (nbdim,nlev)       !< heating by convection             [W/m2]
!!$    REAL(wp) :: zq_cld (nbdim,nlev)       !< heating by stratiform clouds      [W/m2]

    INTEGER  :: ihpbl  (nbdim)            !< location of PBL top given as vertical level index
    REAL(wp) :: zxt_emis(nbdim,ntracer-iqt+1)  !< tracer tendency due to surface emission
                                               !< and dry deposition. "zxtems" in ECHAM5
    INTEGER  :: jk
    INTEGER  :: jks   !< start index for vertical loops
    INTEGER  :: nc    !< number of cells/columns from (jce-jcs+1)
    INTEGER  :: jc
    INTEGER  :: ntrac !< # of tracers excluding water vapour and hydrometeors
                      !< (handled by sub-models, e.g., chemical species)

    ! Coefficient matrices and right-hand-side vectors for the turbulence solver
    ! _btm refers to the lowest model level (i.e., full level "klev", not the surface)

    REAL(wp) :: zaa    (nbdim,nlev,3,nmatrix)       !< coeff. matrices, all variables
    REAL(wp) :: zaa_btm(nbdim,3,nsfc_type,imh:imqv) !< last row of coeff. matrix of heat and moisture
    REAL(wp) :: zbb    (nbdim,nlev,nvar_vdiff)  !< r.h.s., all variables
    REAL(wp) :: zbb_btm(nbdim,nsfc_type,ih_vdiff:iqv_vdiff) !< last row of r.h.s. of heat and moisture

    ! Temporary arrays used by VDIFF, JSBACH

    REAL(wp) :: zfactor_sfc(nbdim)
    REAL(wp) :: zcpt_sfc_tile(nbdim,nsfc_type)  !< dry static energy at surface

    REAL(wp) :: zcptgz   (nbdim,nlev) !< dry static energy
    REAL(wp) :: zrhoh    (nbdim,nlev) !< air density at half levels
    REAL(wp) :: zqshear  (nbdim,nlev) !<
    REAL(wp) :: zthvvar  (nbdim,nlev) !< intermediate value of thvvar
    REAL(wp) :: ztkevn   (nbdim,nlev) !< intermediate value of tke
    REAL(wp) :: zch_tile (nbdim,nsfc_type)   !<  for "nsurf_diag"
!    REAL(wp) :: zchn_tile(nbdim,nsfc_type)   !<  for "nsurf_diag"
!    REAL(wp) :: zcdn_tile(nbdim,nsfc_type)   !<  for "nsurf_diag"
!    REAL(wp) :: zcfnc_tile(nbdim,nsfc_type)  !<  for "nsurf_diag"
    REAL(wp) :: zbn_tile (nbdim,nsfc_type)   !<  for "nsurf_diag"
    REAL(wp) :: zbhn_tile(nbdim,nsfc_type)   !<  for "nsurf_diag"
    REAL(wp) :: zbm_tile (nbdim,nsfc_type)   !<  for "nsurf_diag"
    REAL(wp) :: zbh_tile (nbdim,nsfc_type)   !<  for "nsurf_diag"

    REAL(wp) :: ztte_corr(nbdim)      !< tte correction for snow melt over land (JSBACH)

    ! Temporary array used by GW_HINES

    REAL(wp) :: zdis_gwh(nbdim,nlev)  !<  out, energy dissipation rate [J/s/kg]


    ! Temporary array used by SSODRAG

    REAL(wp) :: zdis_sso(nbdim,nlev)  !<  out, energy dissipation rate [J/s/kg]


    ! Temporary variables used for zenith angle

    REAL(wp) :: zleapfrac
    REAL(wp) :: zyearfrac
    REAL(wp) :: zdeclination_sun
    REAL(wp) :: ztime_dateline

    
    ! Temporary variables used for cloud droplet number concentration

    REAL(wp) :: zprat, zn1, zn2, zcdnc
    LOGICAL  :: lland(nbdim), lglac(nbdim)

    CHARACTER(len=12)  :: str_module = 'e_phy_main'        ! Output of module for 1 line debug
    INTEGER            :: idt_src               ! Determines level of detail for 1 line debug
    idt_src=4

    ! number of cells/columns from index jcs to jce
    nc = jce-jcs+1

    ! start index for vertical loops
    jks=1

    ! 1. Associate pointers

    field  => prm_field(jg)
    tend   => prm_tend (jg)

    ! provisionally copy the incoming tedencies

    tend% temp_phy (jcs:jce,:,jb)   = tend% temp (jcs:jce,:,jb)
    tend%    u_phy (jcs:jce,:,jb)   = tend%    u (jcs:jce,:,jb)
    tend%    v_phy (jcs:jce,:,jb)   = tend%    v (jcs:jce,:,jb)
    tend% qtrc_phy (jcs:jce,:,jb,:) = tend% qtrc (jcs:jce,:,jb,:)

    ! initialize physics heating
    zq_phy(:,:) = 0._wp

    ! 2. local switches and parameters

    ntrac = ntracer-iqt+1  !# of tracers excluding water vapour and hydrometeors

    !------------------------------------------------------------
    ! 3. COMPUTE SOME FIELDS NEEDED BY THE PHYSICAL ROUTINES.
    !------------------------------------------------------------

    ! Use constant volume or constant pressure specific heats for dry air and vapor
    ! for the computation of temperature tendencies.
!!$    IF ( iequations == inh_atmosphere ) THEN
!!$      zcd = cvd
!!$      zcv = cvv
!!$    ELSE
      zcd = cpd
      zcv = cpv
!!$    END IF

    DO jk = 1,nlev
      DO jc = jcs,jce
        !
        ! 3.2 Thickness of model layer in pressure coordinate
        !
        zdelp   (jc,jk) = field% presi_old (jc,jk+1,jb) - field% presi_old (jc,jk,jb)
        !
        ! 3.2b Specific heat of moist air
        !
        zcair   (jc,jk) = zcd+(zcv-zcd)*field%qtrc(jc,jk,jb,iqv)
        zconv   (jc,jk) = 1._wp/(field%mair(jc,jk,jb)*zcair(jc,jk))
        !
        zcpair  (jc,jk) = cpd+(cpv-cpd)*field%qtrc(jc,jk,jb,iqv)
        zcvair  (jc,jk) = cvd+(cvv-cvd)*field%qtrc(jc,jk,jb,iqv)
        !
      END DO
    END DO

    ! 3.3 Weighting factors for fractional surface coverage
    !     Accumulate ice portion for diagnostics

    DO jc=jcs,jce

      ! fraction of land in the grid box. lsmask: land-sea mask, 1.= land

      ! TBD: use fractional mask here
      zfrl(jc) = field% lsmask(jc,jb)

      ! fraction of sea/lake in the grid box
      ! * (1. - fraction of sea ice in the sea/lake part of the grid box)
      ! => fraction of open water in the grid box

      zfrw(jc) = (1._wp-zfrl(jc))*(1._wp-field%seaice(jc,jb))

      ! fraction of sea ice in the grid box
      zfri(jc) = 1._wp-zfrl(jc)-zfrw(jc)
    END DO

    ! 3.4 Merge three pieces of information into one array for vdiff
    IF (ilnd.LE.nsfc_type) zfrc(jcs:jce,ilnd) = zfrl(jcs:jce)
    IF (iwtr.LE.nsfc_type) zfrc(jcs:jce,iwtr) = zfrw(jcs:jce)
    IF (iice.LE.nsfc_type) zfrc(jcs:jce,iice) = zfri(jcs:jce)

    !---------------------------------------------------------------------
    ! 3.9 DETERMINE TROPOPAUSE HEIGHT AND MASS BUDGETS
    !     (Needed only for sub-models. Note: sequence of arguments
    !      different from the original ECHAM6)
    !---------------------------------------------------------------------
    !
    !CALL WMO_tropopause( jce, nbdim, nlev,         &! in
    !                   & ncctop, nccbot, lresum,   &! in
    !                   & field% temp(:,:,jb),      &! in
    !                   & field% presm_old(:,:,jb), &! in
    !                   & field% tropo(:,jb),       &! out for diagnostics
    !                   & itrpwmo, itrpwmop1        )! out for submodel

    !---------------------------------------------------------------------
    ! 3.12 INITIALISATION OF CLOUD DROPLET NUMBER CONCENTRATION 
    !      (1/M**3) USED IN RADLSW AND CLOUD
    !---------------------------------------------------------------------

    DO jc=jcs,jce
      lland(jc) = field%lfland(jc,jb)
      lglac(jc) = lland(jc).AND.field%glac(jc,jb).GT.0._wp
    END DO

    DO jk = 1,nlev
      DO jc = jcs,jce
        !
        zprat=(MIN(8._wp,80000._wp/field%presm_old(jc,jk,jb)))**2

        IF (lland(jc).AND.(.NOT.lglac(jc))) THEN
          zn1= echam_cloud_config% cn1lnd
          zn2= echam_cloud_config% cn2lnd
        ELSE
          zn1= echam_cloud_config% cn1sea
          zn2= echam_cloud_config% cn2sea
        ENDIF
        IF (field%presm_old(jc,jk,jb).LT.80000._wp) THEN
          zcdnc=1.e6_wp*(zn1+(zn2-zn1)*(EXP(1._wp-zprat)))
        ELSE
          zcdnc=zn2*1.e6_wp
        ENDIF
        field% acdnc(jc,jk,jb) = zcdnc
        !
      END DO
    END DO

    !-------------------------------------------------------------------
    ! 3.13 DIAGNOSE CURRENT CLOUD COVER
    !-------------------------------------------------------------------

    itype(jcs:jce) = NINT(field%rtype(jcs:jce,jb))

    IF (phy_config%lcond) THEN
      IF (ltimer) CALL timer_start(timer_cover)

      CALL cover( jce, nbdim, jks,          &! in
        &         nlev, nlevp1,             &! in
        &         itype,  zfrw, zfri,       &! in
        &         field% zf(:,:,jb),        &! in
        &         field% presi_old(:,:,jb), &! in
        &         field% presm_old(:,:,jb), &! in
        &         field%  temp(:,:,jb),     &! in    tm1
        &         field%  qtrc(:,:,jb,iqv), &! in    qm1
        &         field%  qtrc(:,:,jb,iqi), &! in    xim1
        &         field%  aclc(:,:,jb),     &! out   (for "radiation" and "vdiff_down")
        &         field% rintop(:,  jb)    ) ! out   (for output)

      IF (ltimer) CALL timer_stop(timer_cover)
    ENDIF ! lcond

    !-------------------------------------------------------------------
    ! 4. RADIATION PARAMETERISATION
    !-------------------------------------------------------------------
    IF (phy_config%lrad) THEN

       SELECT CASE(izenith)
       CASE(0)
       ! local insolation = constant = global mean insolation (ca. 340 W/m2)
       ! zenith angle = 0,
         ztsi = tsi/4._wp ! scale ztsi by 1/4 to get the correct global mean insolation
       CASE(1)
       ! circular non-seasonal orbit,
       ! perpetual equinox,
       ! no diurnal cycle,
       ! local time always 12:00
       ! --> sin(time of day)=1 ) and zenith angle depends on latitude only
         ztsi = tsi/pi ! because sun is always in local noon, the TSI needs to be
       !               ! scaled by 1/pi to get the correct global mean insolation
       CASE(2)
       ! circular non-seasonal orbit,
       ! perpetual equinox,
       ! no diurnal cycle,
       ! local time always  07:14:15 or 16:45:45
       ! --> sin(time of day)=1/pi and zenith angle depends on latitude only
         ztsi = tsi
       CASE(3)
       ! circular non-seasonal orbit,
       ! perpetual equinox,
       ! with diurnal cycle,
         ztsi = tsi
       CASE(4)
       ! elliptical seasonal orbit, with diurnal cycle
         zleapfrac = 0.681_wp + 0.2422_wp * REAL(datetime%year - 1949,wp) - &
                        REAL((datetime%year - 1949) / 4,wp)
         zyearfrac = 2._wp * pi * (REAL(datetime%yeaday,wp) - 1.0_wp + zleapfrac) / 365.2422_wp
         ztsi = (1.000110_wp + 0.034221_wp * COS(zyearfrac) + 0.001280_wp * SIN(zyearfrac) &
            + 0.000719_wp * COS(2._wp * zyearfrac) + 0.000077_wp * SIN(2._wp * zyearfrac)) * tsi
       CASE(5)
       ! Radiative convective equilibrium
       ! circular non-seasonal orbit,
       ! perpetual equinox,
       ! no diurnal cycle,
       ! the product tsi*cos(zenith angle) should equal 340 W/m2
         ztsi = tsi ! no rescale because tsi has been adjusted in echam_phy_init with ssi_rce
       END SELECT

       ! 4.1 RADIATIVE TRANSFER
       !-----------------------
       IF (ltrig_rad) THEN

          ! store tsfc_rad of this radiatiative transfer timestep in tsfc_radt,
          ! so that it can be reused in radheat in the other timesteps
          field%tsfc_radt(jcs:jce,jb) = field%tsfc_rad(jcs:jce,jb)

          ! to do (for implementing seasonal cycle):
          ! - compute orbit position at datetime_radtran

          SELECT CASE(izenith)

          CASE(0)
          ! local insolation = constant = global mean insolation (ca. 340 W/m2)
          ! zenith angle = 0,

!!$            field%cosmu0(jcs:jce,jb) = 1._wp ! sun in zenith everywhere

          CASE(1)
          ! circular non-seasonal orbit,
          ! perpetual equinox,
          ! no diurnal cycle,
          ! local time always 12:00
          ! --> sin(time of day)=1 ) and zenith angle depends on latitude only

!!$            field%cosmu0(jcs:jce,jb) = COS( p_patch(jg)%cells%center(jcs:jce,jb)%lat )

          CASE(2)
          ! circular non-seasonal orbit,
          ! perpetual equinox,
          ! no diurnal cycle,
          ! local time always  07:14:15 or 16:45:45
          ! --> sin(time of day)=1/pi and zenith angle depends on latitude only

!!$            field%cosmu0(jcs:jce,jb) = COS( p_patch(jg)%cells%center(jcs:jce,jb)%lat )/pi

          CASE(3)
          ! circular non-seasonal orbit,
          ! perpetual equinox,
          ! with diurnal cycle,

!!$            field%cosmu0(jcs:jce,jb) = -COS( p_patch(jg)%cells%center(jcs:jce,jb)%lat ) &
!!$                                     & *COS( p_patch(jg)%cells%center(jcs:jce,jb)%lon   &
!!$                                     &      +2._wp*pi*datetime_radtran%daytim )

          CASE(4)
          ! elliptical seasonal orbit,
          !  with diurnal cycle

            zleapfrac = 0.681_wp + 0.2422_wp * REAL(datetime%year - 1949,wp) - &
                        REAL((datetime%year - 1949) / 4,wp)
            zyearfrac = 2._wp * pi * (REAL(datetime%yeaday,wp) - 1.0_wp + zleapfrac) / 365.2422_wp
            zdeclination_sun = 0.006918_wp - 0.399912_wp * COS(zyearfrac) + &
                               0.070257_wp * SIN(zyearfrac) -               &
                               0.006758_wp * COS(2._wp * zyearfrac) +       &
                               0.000907_wp * SIN(2._wp * zyearfrac) -       &
                               0.002697_wp * COS(3._wp * zyearfrac) +       &
                               0.001480_wp * SIN(3._wp * zyearfrac)
            ztime_dateline = ((REAL(datetime%hour,wp) * 3600._wp + &
                              REAL(datetime%minute,wp) * 60._wp +  &
                              REAL(datetime%second,wp)) /          &
                              REAL(datetime%daylen,wp)) - 0.5_wp
            ztime_dateline = ztime_dateline * 2._wp * pi + 0.000075_wp +              &
               0.001868_wp * COS(zyearfrac) - 0.032077_wp * SIN(zyearfrac) -          &
               0.014615_wp * COS(2._wp * zyearfrac) - 0.040849_wp * SIN(2._wp * zyearfrac)

!!$            field%cosmu0(jcs:jce,jb) = SIN(zdeclination_sun) * SIN(p_patch(jg)%cells%center(jcs:jce,jb)%lat) + &
!!$                                       COS(zdeclination_sun) * COS(p_patch(jg)%cells%center(jcs:jce,jb)%lat) * &
!!$                                       COS(ztime_dateline + p_patch(jg)%cells%center(jcs:jce,jb)%lon)
          CASE(5)
          ! Radiative convective equilibrium
          ! circular non-seasonal orbit,
          ! perpetual equinox,
          ! no diurnal cycle,
          ! see Popke et al. 2013 and Cronin 2013
          !cosmu0 = pi/4._wp ! zenith = 45 deg
          !cosmu0 = 2._wp/3._wp ! Cronin: zenith = 48.19

!!$            field%cosmu0(jcs:jce,jb) = 0.7854_wp ! Popke: zenith = 38

          END SELECT

        IF (ltimer) CALL timer_start(timer_radiation)

        CALL psrad_radiation(      &
        & jg                      ,&!< in  domain index
        & jb                      ,&!< in  block index
        & kproma     = jce        ,&!< in  end index for loop over block
        & kbdim      = nbdim      ,&!< in  dimension of block over cells
        & klev       = nlev       ,&!< in  number of full levels = number of layers
        & klevp1     = nlevp1     ,&!< in  number of half levels = number of layer interfaces
        & ktype      = itype(:)   ,&!< in  type of convection
        & loland     = lland      ,&!< in  land-sea mask. (logical)
        & loglac     = lglac      ,&!< in  glacier mask (logical)
        & datetime   = datetime   ,&!< in  actual time step
        & pcos_mu0   = field%cosmu0_rad(:,jb)  ,&!< in  solar zenith angle
        & alb_vis_dir= field%albvisdir(:,jb)   ,&!< in  surface albedo for visible range, direct
        & alb_nir_dir= field%albnirdir(:,jb)   ,&!< in  surface albedo for near IR range, direct
        & alb_vis_dif= field%albvisdif(:,jb)   ,&!< in  surface albedo for visible range, diffuse
        & alb_nir_dif= field%albnirdif(:,jb)   ,&!< in  surface albedo for near IR range, diffuse
        & tk_sfc     = field%tsfc_radt(:,jb)   ,&!< in  grid box mean surface temperature
        & zf         = field%zf(:,:,jb)        ,&!< in  geometric height at full level      [m]
        & zh         = field%zh(:,:,jb)        ,&!< in  geometric height at half level      [m]
        & dz         = field%dz(:,:,jb)        ,&!< in  geometric height thickness of layer [m]
        & pp_hl      = field%presi_old(:,:,jb) ,&!< in  pressure at half levels at t-dt [Pa]
        & pp_fl      = field%presm_old(:,:,jb) ,&!< in  pressure at full levels at t-dt [Pa]
        & tk_fl      = field%temp(:,:,jb)      ,&!< in  tk_fl  = temperature at full level at t-dt
        & xm_dry     = field%mdry(:,:,jb)      ,&!< in  dry air mass in layer [kg/m2]
        & xm_trc     = field%mtrc(:,:,jb,:)    ,&!< in  tracer  mass in layer [kg/m2]
        & xm_ozn     = field%o3(:,:,jb)        ,&!< inout  ozone  mass mixing ratio [kg/kg]
        & cdnc       = field% acdnc(:,:,jb)    ,&!< in   cloud droplet number conc
        & cld_frc    = field% aclc(:,:,jb)     ,&!< in   cloud fraction [m2/m2]
        & cld_cvr    = field%aclcov(:,jb)      ,&!< out  total cloud cover
        & vis_frc_sfc= field%visfrcsfc(:,jb)   ,&!< out  visible (250-680nm) fraction of net surface radiation
        & par_dn_sfc = field%partrmdnsfc(:,jb) ,&!< out  downward photosynthetically active radiation (par) at surface
        & nir_dff_frc= field%nirdffsfc(:,jb)   ,&!< out  diffuse fraction of downward surface near-infrared radiation
        & vis_dff_frc= field%visdffsfc(:,jb)   ,&!< out  diffuse fraction of downward surface visible radiation
        & par_dff_frc= field%pardffsfc(:,jb)   ,&!< out  diffuse fraction of downward surface par
        & lw_dnw_clr = field%rldcs(:,:,jb)     ,&!< out  Clear-sky net longwave  at all levels
        & lw_upw_clr = field%rlucs(:,:,jb)     ,&!< out  Clear-sky net longwave  at all levels
        & sw_dnw_clr = field%rsdcs(:,:,jb)     ,&!< out  Clear-sky net shortwave at all levels
        & sw_upw_clr = field%rsucs(:,:,jb)     ,&!< out  Clear-sky net shortwave at all levels
        & lw_dnw     = field%rld  (:,:,jb)     ,&!< out  All-sky net longwave  at all levels
        & lw_upw     = field%rlu  (:,:,jb)     ,&!< out  All-sky net longwave  at all levels
        & sw_dnw     = field%rsd  (:,:,jb)     ,&!< out  All-sky net longwave  at all levels
        & sw_upw     = field%rsu  (:,:,jb)      &!< out  All-sky net longwave  at all levels
        &                           )

        flux_factor(1:jce) = 1._wp / (psctm*field%cosmu0_rad(1:jce,jb))
        field%rsn  (1:jce,1:nlevp1,jb) = (field%rsd  (1:jce,1:nlevp1,jb) - field%rsu  (1:jce,1:nlevp1,jb)) * &
             &                           SPREAD(flux_factor(1:jce),2,nlevp1)
        field%rsncs(1:jce,1:nlevp1,jb) = (field%rsdcs(1:jce,1:nlevp1,jb) - field%rsucs(1:jce,1:nlevp1,jb)) * &
             &                           SPREAD(flux_factor(1:jce),2,nlevp1)
        field%partrmdnsfc(1:jce,jb)    = field%partrmdnsfc(1:jce,jb) * flux_factor(1:jce)
        field%rln  (1:jce,1:nlevp1,jb) = (field%rld  (1:jce,1:nlevp1,jb) - field%rlu  (1:jce,1:nlevp1,jb))
        field%rlncs(1:jce,1:nlevp1,jb) = (field%rldcs(1:jce,1:nlevp1,jb) - field%rlucs(1:jce,1:nlevp1,jb))
        
        IF (ltimer) CALL timer_stop(timer_radiation)

      END IF ! ltrig_rad

      ! 4.2 RADIATIVE HEATING
      !----------------------

      ! to do:
      ! - compute orbit position at datetime

      ! - solar incoming flux at TOA
      field% rsdt(jcs:jce,jb) = MAX(0._wp,field%cosmu0(jcs:jce,jb)) * ztsi  ! instantaneous for radheat

      ! radheat first computes the shortwave and longwave radiation for the current time step from transmissivity and
      ! the longwave flux at the radiation time step and, from there, the radiative heating due to sw and lw radiation.
      ! If radiation is called every time step, the longwave flux is not changed.

      IF (ltimer) CALL timer_start(timer_radheat)

      CALL radheating (                                &
        !
        ! input
        ! -----
        !
        & jcs        = jcs                            ,&! in    loop start index
        & jce        = jce                            ,&! in    loop end index
        & kbdim      = nbdim                          ,&! in    dimension size
        & klev       = nlev                           ,&! in    vertical dimension size
        & klevp1     = nlevp1                         ,&! in    vertical dimension size
        !
        & cosmu0_rad = field%cosmu0_rad(:,jb)         ,&! in    solar zenith angle at radiation time
        & cosmu0     = field%cosmu0    (:,jb)         ,&! in    solar zenith angle at current   time
        !
        & prsdt      = field%rsdt               (:,jb),&! in    solar incoming flux at TOA [W/m2]
        & pemiss     = ext_data(jg)%atm%emis_rad(:,jb),&! in    lw sfc emissivity
        & ptsfc      = field%tsfc_rad (:,jb)          ,&! in    rad. surface temperature now         [K]
        & ptsfctrad  = field%tsfc_radt(:,jb)          ,&! in    rad. surface temp. at last rad. step [K]
        & lwflx_up_sfc_rs = field%rlu    (:,nlevp1,jb),&! in    surface longwave upward flux at last rad. step [W/m2]
        !
        & rsd        = field%rsd              (:,:,jb),&! in    all-sky   shortwave downward flux at last radiation step [W/m2]
        & rsu        = field%rsu              (:,:,jb),&! in    all-sky   shortwave upward   flux at last radiation step [W/m2]
        & rld        = field%rsd              (:,:,jb),&! in    all-sky   longwave  downward flux at last radiation step [W/m2]
        & rlu        = field%rsu              (:,:,jb),&! in    all-sky   longwave  upward   flux at last radiation step [W/m2]
        !
        & rsdcs      = field%rsdcs            (:,:,jb),&! in    clear-sky shortwave downward flux at last radiation step [W/m2]
        & rsucs      = field%rsucs            (:,:,jb),&! in    clear-sky shortwave upward   flux at last radiation step [W/m2]
        & rldcs      = field%rsdcs            (:,:,jb),&! in    clear-sky longwave  downward flux at last radiation step [W/m2]
        & rlucs      = field%rsucs            (:,:,jb),&! in    clear-sky longwave  upward   flux at last radiation step [W/m2]
        !
        & ptrmsw     = field%rsn              (:,:,jb),&! in    shortwave net transmissivity at last rad. step []
        & rln        = field%rln              (:,:,jb),&! in    longwave net flux at last rad. step [W/m2]
        & ptrmswclr  = field%rsncs            (:,:,jb),&! in    shortwave net transmissivity at last rad. step clear sky []
        & rlncs      = field%rlncs            (:,:,jb),&! in    longwave net flux at last rad. step clear sky [W/m2]
        !
        ! output
        ! ------
        !
        & pq_rsw     = zq_rsw                   (:,:) ,&! out   rad. heating by SW           [W/m2]
        & pq_rlw     = zq_rlw                   (:,:) ,&! out   rad. heating by LW           [W/m2]
        !
        & rsns       = field%rsns               (:,jb),&! out   shortwave surface net flux   [W/m2]
        & rlns       = field%rlns               (:,jb),&! out   longwave surface net flux    [W/m2]
        & rsnt       = field%rsnt               (:,jb),&! out   shortwave toa net flux       [W/m2]
        & rlnt       = field%rlnt               (:,jb),&! out   longwave toa net flux        [W/m2]
        & lwflx_up_sfc = field%rlus             (:,jb)) ! out   longwave surface upward flux [W/m2]

      IF (ltimer) CALL timer_stop(timer_radheat)

      ! heating accumulated
      zq_phy(jcs:jce,:) = zq_phy(jcs:jce,:) + zq_rsw(jcs:jce,:) + zq_rlw(jcs:jce,:)
      
      ! tendencies
      tend%temp_rsw(jcs:jce,:,jb) = zq_rsw(jcs:jce,:) * zconv(jcs:jce,:)
      tend%temp_rlw(jcs:jce,:,jb) = zq_rlw(jcs:jce,:) * zconv(jcs:jce,:)

      ! tendencies accumulated
      tend% temp(jcs:jce,:,jb) = tend% temp     (jcs:jce,:,jb) &
        &                      + tend% temp_rsw (jcs:jce,:,jb) &
        &                      + tend% temp_rlw (jcs:jce,:,jb)

    ELSE   ! If computation of radiative heating is by-passed

      tend%temp_rsw(jcs:jce,:,jb) = 0.0_wp
      tend%temp_rlw(jcs:jce,:,jb) = 0.0_wp

      field%rsdt(jcs:jce,jb)= 0.0_wp

    END IF ! lrad

    ! Compute VIS/NIR/PAR shortwave fluxes for surface processes

    field%vissfc   (jcs:jce,jb) = field%rsns(jcs:jce,jb) * field%visfrcsfc   (jcs:jce,jb)
    field%nirsfc   (jcs:jce,jb) = field%rsns(jcs:jce,jb) - field%vissfc      (jcs:jce,jb)
    field%parsfcdn (jcs:jce,jb) = field%rsdt(jcs:jce,jb) * field%partrmdnsfc (jcs:jce,jb)

    !-------------------------------------------------------------------
    ! 5. BOUNDARY LAYER AND SURFACE PROCESSES
    !-------------------------------------------------------------------
    ! Note: In ECHAM this part is located between "CALL radiation" and
    !       "CALL radheat".
    !
    ! 5.1 Emission of aerosols or other tracers. Not implemented yet.

      IF (ntrac>0) THEN
        !CALL tracer_emission()
        zxt_emis(jcs:jce,:) = 0._wp
      ENDIF
    !
    ! 5.2 Dry deposition of aerosols or other tracers. Not implemented yet.
    ! CALL dry_deposition()
    !

    ! 5.3 Turbulent mixing, part I:
    !     computation of exchange coefficients in the atmosphere and at the surface;
    !     build up the tridiagonal linear algebraic system;
    !     downward sweep (Gaussian elimination from top till level nlev-1)

    IF (phy_config%lvdiff) THEN
      IF (ltimer) CALL timer_start(timer_vdiff_down)


      CALL vdiff_down( vdiff_config%lsfc_mom_flux,      &! in
                     & vdiff_config%lsfc_heat_flux,     &! in
                     & jce, nbdim, nlev, nlevm1, nlevp1,&! in
                     & ntrac, nsfc_type,                &! in
                     & iwtr, iice, ilnd,                &! in, indices of different surface types
                     & psteplen,                        &! in, time step (2*dt if leapfrog)
                     & field%coriol(:,jb),              &! in, Coriolis parameter
                     & zfrc(:,:),                       &! in, area fraction of each sfc type
                     & field% tsfc_tile(:,jb,:),        &! in, surface temperature
                     & field% ocu (:,jb),               &! in, ocean sfc velocity, u-component
                     & field% ocv (:,jb),               &! in, ocean sfc velocity, v-component
                     & field% presi_old(:,nlevp1,jb),   &! in, sfc pressure
                     & field%    u(:,:,jb),             &! in, um1
                     & field%    v(:,:,jb),             &! in, vm1
                     & field% temp(:,:,jb),             &! in, tm1
                     & field% qtrc(:,:,jb,iqv),         &! in, qm1
                     & field% qtrc(:,:,jb,iqc),         &! in, xlm1
                     & field% qtrc(:,:,jb,iqi),         &! in, xim1
                     & field%   qx(:,:,jb),             &! in, xlm1 + xim1
                     & field% qtrc(:,:,jb,iqt:),        &! in, xtm1
                     & field% presi_old(:,:,jb),        &! in, aphm1
                     & field% presm_old(:,:,jb),        &! in, apm1
                     & zdelp(:,:),                      &! in, layer thickness [Pa]
                     & field% geom(:,:,jb),             &! in, pgeom1 = geopotential above ground
                     & field% geoi(:,:,jb),             &! in, pgeohm1 = half-level geopotential
                     & field%   tv(:,:,jb),             &! in, virtual temperaturea
                     & field% aclc(:,:,jb),             &! in, cloud fraction
                     & zxt_emis,                        &! in, zxtems
                     & field% thvvar(:,:,jb),           &! in, variance of theta_v at step t-dt
                     & field%   xvar(:,:,jb),           &! in
                     & field% z0m_tile(:,jb,:),         &! in
                     & field%  tkem1(:,:,jb),           &! in, TKE at step t-dt
                     & field%  ustar(:,  jb),           &! inout
                     & field%  wstar(:,  jb),           &! out, convective velocity scale
                     & field%  wstar_tile(:,jb,:),      &! inout, convective velocity scale (each sfc type)
                     & field% qs_sfc_tile(:,jb,:),      &! out, sfc specific humidity at saturation
                     & ihpbl(:),                        &! out, for "vdiff_up"
                     & field%    ghpbl(:,jb),           &! out, for output
                     & field%      ri (:,:,jb),         &! out, for output
                     & zri_tile (:,:),                  &! out, for nsurf_diag
                     & field%  mixlen (:,:,jb),         &! out, for output
                     & field% cfm     (:,:,jb),         &! out, for output
                     & field% cfm_tile(:,jb,:),         &! out, for output and "vdiff_up"
                     & field% cfh     (:,:,jb),         &! out, for output
                     & field% cfh_tile(:,jb,:),         &! out, for output and "vdiff_up"
                     & field% cfv     (:,:,jb),         &! out, for output
                     & field% cftke   (:,:,jb),         &! out, for output
                     & field% cfthv   (:,:,jb),         &! out, for output
                     & zaa, zaa_btm, zbb, zbb_btm,      &! out, for "vdiff_up"
                     & zfactor_sfc(:),                  &! out, for "vdiff_up"
                     & zcpt_sfc_tile(:,:),              &! out, for "vdiff_up"
                     & zcptgz(:,:), zrhoh(:,:),         &! out, for "vdiff_up"
                     & zqshear(:,:),                    &! out, for "vdiff_up"
                     & zthvvar(:,:),                    &! out, for "vdiff_up"
                     & field%   thvsig(:,  jb),         &! out, for "cucall"
                     & ztkevn (:,:),                    &! out, for "vdiff_up"
                     & zch_tile(:,:),                   &! out, for "nsurf_diag"
!                     & zchn_tile(:,:),                  &! out, for "nsurf_diag"
!                     & zcdn_tile(:,:),                  &! out, for "nsurf_diag"
!                     & zcfnc_tile(:,:),                 &! out, for "nsurf_diag"
                     & zbn_tile(:,:),                   &! out, for "nsurf_diag"
                     & zbhn_tile(:,:),                  &! out, for "nsurf_diag"
                     & zbm_tile(:,:),                   &! out, for "nsurf_diag"
                     & zbh_tile(:,:),                   &! out, for "nsurf_diag"
                     & pcsat = field% csat(:,jb),       &! in, optional, area fraction with wet land surface
                     & pcair = field% cair(:,jb),       &! in, optional, area fraction with wet land surface (air)
                     & paz0lh = field% z0h_lnd(:,jb))     ! in, optional, roughness length for heat over land

      IF (ltimer) CALL timer_stop(timer_vdiff_down)

    ! 5.4 Surface processes that provide time-dependent lower boundary
    !     condition for wind, temperature, tracer concentraion, etc.
    !KF To avoid
        field% lhflx_tile(:,jb,:) = 0._wp
        field% shflx_tile(:,jb,:) = 0._wp
        field% evap_tile (:,jb,:) = 0._wp

        IF (ltimer) CALL timer_start(timer_surface)

        CALL update_surface( vdiff_config%lsfc_heat_flux,  &! in
          & vdiff_config%lsfc_mom_flux,   &! in
          & pdtime, psteplen,             &! in, time steps
          & jg, jce, nbdim, field%kice,   &! in
          & nlev, nsfc_type,              &! in
          & iwtr, iice, ilnd,             &! in, indices of surface types
          & zfrc(:,:),                    &! in, area fraction
          & field% cfh_tile(:,jb,:),      &! in, from "vdiff_down"
          & field% cfm_tile(:,jb,:),      &! in, from "vdiff_down"
          & zfactor_sfc(:),               &! in, from "vdiff_down"
          & field% ocu (:,jb),            &! in, ocean sfc velocity, u-component
          & field% ocv (:,jb),            &! in, ocean sfc velocity, v-component
          & zaa, zaa_btm, zbb, zbb_btm,   &! inout
          & zcpt_sfc_tile(:,:),           &! inout, from "vdiff_down", for "vdiff_up"
          & field%qs_sfc_tile(:,jb,:),    &! inout, from "vdiff_down", for "vdiff_up"
          & field% tsfc_tile(:,jb,:),     &! inout
          & field%u_stress    (:,  jb),   &! out
          & field%v_stress    (:,  jb),   &! out
          & field% lhflx      (:,  jb),   &! out
          & field% shflx      (:,  jb),   &! out
          & field%  evap      (:,  jb),   &! out, for "cucall"
          & field%u_stress_tile  (:,jb,:),   &! out
          & field%v_stress_tile  (:,jb,:),   &! out
          & field% lhflx_tile    (:,jb,:),   &! out
          & field% shflx_tile    (:,jb,:),   &! out
          & field% dshflx_dT_tile(:,jb,:),   &! out for Sea ice
          & field%  evap_tile    (:,jb,:),   &! out
                                !! optional
          & nblock = jb,                  &! in
          & lsm = field%lsmask(:,jb), &!< in, land-sea mask
          & pu = field% u(:,nlev,jb),     &! in, um1
          & pv = field% v(:,nlev,jb),     &! in, vm1
          & ptemp = field% temp(:,nlev,jb), &! in, tm1
          & pq = field% qtrc(:,nlev,jb,iqv),  &! in, qm1
          & prsfl = field% rsfl(:,jb),    &! in, rain surface large scale (from cloud)
          & prsfc = field% rsfc(:,jb),    &! in, rain surface concective (from cucall)
          & pssfl = field% ssfl(:,jb),    &! in, snow surface large scale (from cloud)
          & pssfc = field% ssfc(:,jb),    &! in, snow surface concective (from cucall)
          & plw         = field% rlns (:,jb), &! inout, net surface longwave flux [W/m2]
          & plw_down    = field% rlns (:,jb) + field%rlus(:,jb),     &! in, downward surface longwave flux [W/m2]
          & psw         = field% rsns     (:,jb), &! inout, net surface shortwave flux [W/m2]
          & pswvis      = field% vissfc   (:,jb), &! in, net surface shortwave flux in visible range [W/m2]
          & pswnir      = field% nirsfc   (:,jb), &! in, net surface shortwave flux in NIR range [W/m2]
          & pswpar_down = field% parsfcdn (:,jb), &! in, downward surface shortwave flux in PAR range [W/m2]
          & pvisdff     = field% visdffsfc(:,jb), &! in, diffuse fraction in visible shortwave surface flux
          & pnirdff     = field% nirdffsfc(:,jb), &! in, diffuse fraction in NIR shortwave surface flux
          & ppardff     = field% pardffsfc(:,jb), &! in, diffuse fraction in PAR shortwave surface flux
          & presi_old = field% presi_old(:,:,jb),&! in, paphm1, half level pressure
          & pcosmu0 = field% cosmu0(:,jb),&! in, amu0_x, cos of zenith angle
          & pch_tile = zch_tile(:,:),     &! in, from "vdiff_down" for JSBACH
          & pcsat = field%csat(:,jb),      &! inout, area fraction with wet land surface
          & pcair = field%cair(:,jb),      &! inout, area fraction with wet land surface (air)
          & tte_corr = ztte_corr(:),       &! out, tte correction for snow melt over land
          & z0m_tile = field% z0m_tile(:,jb,:), &! inout, roughness length for momentum over tiles
          & z0h_lnd  = field% z0h_lnd (:,jb),   &! out, roughness length for heat over land
          & albvisdir      = field% albvisdir     (:,jb)  ,                    &! inout
          & albnirdir      = field% albnirdir     (:,jb)  ,                    &! inout
          & albvisdif      = field% albvisdif     (:,jb)  ,                    &! inout
          & albnirdif      = field% albnirdif     (:,jb)  ,                    &! inout
          & albvisdir_tile = field% albvisdir_tile(:,jb,:),                    &! inout
          & albnirdir_tile = field% albnirdir_tile(:,jb,:),                    &! inout
          & albvisdif_tile = field% albvisdif_tile(:,jb,:),                    &! inout
          & albnirdif_tile = field% albnirdif_tile(:,jb,:),                    &! inout
          & albedo         = field% albedo        (:,jb)  ,                    &! inout
          & albedo_tile    = field% albedo_tile(:,jb,:),                       &! inout
          & ptsfc     = field%tsfc    (:,jb),                      &! out
          & ptsfc_rad = field%tsfc_rad(:,jb),                      &! out
          & plwflx_tile = field%lwflxsfc_tile(:,jb,:),             &! out (for coupling)
          & pswflx_tile = field%swflxsfc_tile(:,jb,:),             &! out (for coupling)
          & Tsurf = field% Tsurf(:,:,jb),  &! inout, for sea ice
          & T1    = field% T1   (:,:,jb),  &! inout, for sea ice
          & T2    = field% T2   (:,:,jb),  &! inout, for sea ice
          & hi    = field% hi   (:,:,jb),  &! in, for sea ice
          & hs    = field% hs   (:,:,jb),  &! in, for sea ice
          & conc  = field% conc (:,:,jb),  &! in, for sea ice
          & Qtop  = field% Qtop (:,:,jb),  &! out, for sea ice
          & Qbot  = field% Qbot (:,:,jb),  &! out, for sea ice
          & albvisdir_ice = field% albvisdir_ice(:,:,jb), &! inout ice albedos
          & albnirdir_ice = field% albnirdir_ice(:,:,jb), &! inout
          & albvisdif_ice = field% albvisdif_ice(:,:,jb), &! inout
          & albnirdif_ice = field% albnirdif_ice(:,:,jb))  ! inout

        IF (ltimer) CALL timer_stop(timer_surface)

    ! 5.5 Turbulent mixing, part II:
    !     - Elimination for the lowest model level using boundary conditions
    !       provided by the surface model(s);
    !     - Back substitution to get solution of the tridiagonal system;
    !     - Compute tendencies and additional diagnostics.

      IF (ltimer) CALL timer_start(timer_vdiff_up)

      CALL vdiff_up( jce, nbdim, nlev, nlevm1, nlevp1,&! in
                   & ntrac, nsfc_type,                &! in
                   & iwtr,                            &! in, indices of different sfc types
                   & pdtime, psteplen,                &! in, time steps
                   & zfrc(:,:),                       &! in, area fraction of each sfc type
                   & field% cfm_tile(:,jb,:),         &! in
                   & zaa,                             &! in, from "vdiff_down"
                   &   ihpbl(:),                      &! in, from "vdiff_down"
                   &  zcptgz(:,:),                    &! in, from "vdiff_down"
                   &   zrhoh(:,:),                    &! in, from "vdiff_down"
                   & zqshear(:,:),                    &! in, from "vdiff_down"
                   & field%    u(:,:,jb),             &! in, um1
                   & field%    v(:,:,jb),             &! in, vm1
                   & field% temp(:,:,jb),             &! in, tm1
                   & field% qtrc(:,:,jb,iqv),         &! in, qm1
                   & field% qtrc(:,:,jb,iqc),         &! in, xlm1
                   & field% qtrc(:,:,jb,iqi),         &! in, xim1
                   & field% qtrc(:,:,jb,iqt:),        &! in, xtm1
                   & zcd,                             &! in, specific heat of dry air
                   & zcv,                             &! in, specific heat of water vapor
                   & zdelp(:,:),                      &! in, layer thickness [Pa]
                   & field% geom(:,:,jb),             &! in, pgeom1 = geopotential above ground
                   &      ztkevn(:,:),                &! in, tke at intermediate time step
                   & field%tkem1(:,:,jb),             &! in, TKE at step t-dt
                   & ztte_corr(:),                    &! in
                   & zbb,                             &! inout
                   & zthvvar(:,:),                    &! inout
                   & field%   xvar(:,:,jb),           &! inout
                   & field% z0m_tile(:,jb,:),         &! inout
                   & field% kedisp(:,  jb),           &! inout, "vdis" in ECHAM
                   &  tend%    u_vdf(:,:,jb),         &! out
                   &  tend%    v_vdf(:,:,jb),         &! out
                   &  tend% temp_vdf(:,:,jb),         &! out
                   &  tend% qtrc_vdf(:,:,jb,iqv),     &! out
                   &  tend% qtrc_vdf(:,:,jb,iqc),     &! out
                   &  tend% qtrc_vdf(:,:,jb,iqi),     &! out
                   &  tend% qtrc_vdf(:,:,jb,iqt:),    &! out
                   &  zqtvar_prod,                    &! out, for "cloud" ("zvdiffp" in echam)
                   &  zvmixtau,                       &! out, for "cloud"
                   & field%   z0m   (:,  jb),         &! out, for the next step
                   & field%   thvvar(:,:,jb),         &! out, for the next step
                   & field%   thvsig(:,  jb),         &! out, for "cucall"
                   & field%      tke(:,:,jb),         &! out
                   & field%   sh_vdiff(:,  jb),       &! out, for energy diagnostic
                   & field%   qv_vdiff(:,  jb)        )! out, for energy diagnostic

      IF (ltimer) CALL timer_stop(timer_vdiff_up)

      ! tendencies accumulated
      tend%    u(jcs:jce,:,jb)      = tend%    u(jcs:jce,:,jb)      + tend%    u_vdf(jcs:jce,:,jb)
      tend%    v(jcs:jce,:,jb)      = tend%    v(jcs:jce,:,jb)      + tend%    v_vdf(jcs:jce,:,jb)
      tend% temp(jcs:jce,:,jb)      = tend% temp(jcs:jce,:,jb)      + tend% temp_vdf(jcs:jce,:,jb)
      tend% qtrc(jcs:jce,:,jb,iqv)  = tend% qtrc(jcs:jce,:,jb,iqv)  + tend% qtrc_vdf(jcs:jce,:,jb,iqv)
      tend% qtrc(jcs:jce,:,jb,iqc)  = tend% qtrc(jcs:jce,:,jb,iqc)  + tend% qtrc_vdf(jcs:jce,:,jb,iqc)
      tend% qtrc(jcs:jce,:,jb,iqi)  = tend% qtrc(jcs:jce,:,jb,iqi)  + tend% qtrc_vdf(jcs:jce,:,jb,iqi)
      tend% qtrc(jcs:jce,:,jb,iqt:) = tend% qtrc(jcs:jce,:,jb,iqt:) + tend% qtrc_vdf(jcs:jce,:,jb,iqt:)

!    ! TIME FILTER FOR TURBULENT KINETIC ENERGY
!
!    IF(.NOT.lstart) THEN
!      zeps=eps
!    ELSE
!      zeps=0._wp
!    END IF
!    DO 397 jk=ktdia,klev
!      DO 396 jl=1,kproma
!        ptkem1(jl,jk)=ptkem(jl,jk)                                    &
!                  +zeps*(ptkem1(jl,jk)-2._wp*ptkem(jl,jk)+ptke(jl,jk))
!        ptkem(jl,jk)=ptke(jl,jk)
!396   END DO
!397 END DO

      IF (ABS(pdtime*2._wp-psteplen)<1e-6_wp) THEN
        ! Leapfrog scheme. No Asselin filter. Just swap time steps
        field% tkem1(jcs:jce,:,jb) = field% tkem0(jcs:jce,:,jb)
        field% tkem0(jcs:jce,:,jb) = field% tke  (jcs:jce,:,jb)
      ELSE
        ! 2-tl-scheme
        field% tkem1(jcs:jce,:,jb) = field% tke  (jcs:jce,:,jb)
      ENDIF

    ! 5.6 Turbulent mixing, part III:
    !     - Further diagnostics.

    CALL nsurf_diag( jce, nbdim, nsfc_type,           &! in
                   & ilnd,                            &! in
                   & zfrc(:,:),                       &! in
                   & field%  qtrc(:,nlev,jb,iqv),     &! in humidity qm1
                   & field%  temp(:,nlev,jb),         &! in tm1
                   & field% presm_old(:,nlev,jb),     &! in, apm1
                   & field% presi_old(:,nlevp1,jb),   &! in, aphm1
                   & field%   qx(:,nlev,jb),          &! in, xlm1 + xim1
                   & field%    u(:,nlev,jb),          &! in, um1
                   & field%    v(:,nlev,jb),          &! in, vm1
                   & field% ocu (:,jb),               &! in, ocean sfc velocity, u-component
                   & field% ocv (:,jb),               &! in, ocean sfc velocity, v-component
                   & field%  geom(:,nlev,jb),         &! in geopotential above surface
                   & zcptgz(:,nlev),                  &! in dry static energy
                   & zcpt_sfc_tile(:,:),              &! in dry static energy
                   & zbn_tile(:,:),                   &! in for diagnostic
                   & zbhn_tile(:,:),                  &! in for diagnostic
                   & zbh_tile(:,:),                   &! in for diagnostic
                   & zbm_tile(:,:),                   &! in for diagnostic
                   & zri_tile(:,:),                   &! in 
                   & field%sfcWind(:,  jb),           &! out 10m windspeed
                   & field%    tas(:,  jb),           &! out temperature in 2m
                   & field%   dew2(:,  jb),           &! out dew point temperature in 2m
                   & field%    uas(:,  jb),           &! out zonal wind in 10m
                   & field%    vas(:,  jb),           &! out meridional wind in 10m
                   & field%tasmax (:,  jb),           &! out max 2m temperature
                   & field%tasmin (:,  jb),           &! out min 2m temperature
                   & field%sfcWind_tile(:,jb,:),      &! out 10m windspeed on tiles
                   & field%    tas_tile(:,jb,:),      &! out temperature in 2m on tiles
                   & field%   dew2_tile(:,jb,:),      &! out dew point temperature in 2m on tiles
                   & field%    uas_tile(:,jb,:),      &! out zonal wind in 10m on tiles
                   & field%    vas_tile(:,jb,:)       )! out meridional wind in 10m on tiles

    ELSE
      zvmixtau   (jcs:jce,:) = 0._wp
      field% evap(jcs:jce,jb)= 0._wp
      zqtvar_prod(jcs:jce,:) = 0._wp

      tend%    u_vdf(jcs:jce,:,jb)      = 0._wp
      tend%    v_vdf(jcs:jce,:,jb)      = 0._wp
      tend% temp_vdf(jcs:jce,:,jb)      = 0._wp
      tend% qtrc_vdf(jcs:jce,:,jb,iqv)  = 0._wp
      tend% qtrc_vdf(jcs:jce,:,jb,iqc)  = 0._wp
      tend% qtrc_vdf(jcs:jce,:,jb,iqi)  = 0._wp
      tend% qtrc_vdf(jcs:jce,:,jb,iqt:) = 0._wp

    ENDIF !lvdiff

    IF (phy_config%lrad) THEN

      ! Heating due to the fact that surface model only used part of longwave radiation to compute new surface temperature
      zq_rlw_impl(jcs:jce) =                                            &
        & ( (field%rln(jcs:jce,nlev,jb) - field%rlns(jcs:jce,jb)) ) &  ! new heating from new rlns
        & - zq_rlw(jcs:jce,nlev)                                       ! old heating from radheat

      ! Heating accumulated
      zq_phy(jcs:jce,nlev) = zq_phy(jcs:jce,nlev) + zq_rlw_impl(jcs:jce)

      ! Tendency
      tend%temp_rlw_impl(jcs:jce,jb) = zq_rlw_impl(jcs:jce) * zconv(jcs:jce,nlev)

      ! Tendencies accumulated
      tend%temp(jcs:jce,nlev,jb) = tend%temp(jcs:jce,nlev,jb) + tend%temp_rlw_impl(jcs:jce,jb)

    ELSE

      tend%temp_rlw_impl(jcs:jce,jb) = 0._wp

    END IF

    !-------------------------------------------------------------------
    ! 6. ATMOSPHERIC GRAVITY WAVES
    !-------------------------------------------------------------------

    ! 6.1   CALL SUBROUTINE GW_HINES

    IF (phy_config%lgw_hines) THEN

      zlat_deg(jcs:jce) = p_patch(jg)%cells%center(jcs:jce,jb)%lat * 180._wp/pi

      IF (ltimer) call timer_start(timer_gw_hines)

      CALL gw_hines ( jg                       ,&
        &             nbdim                    ,&
        &             jcs                      ,&
        &             jce                      ,&
        &             nc                       ,&
        &             nlev                     ,&
        &             field% presi_old(:,:,jb) ,&
        &             field% presm_old(:,:,jb) ,&
        &             field% temp(:,:,jb)      ,&
        &             field%    u(:,:,jb)      ,&
        &             field%    v(:,:,jb)      ,&
        &             zlat_deg(:)              ,&
!!$        &             aprflux(:,krow)          ,&
        &             zdis_gwh(:,:)            ,&
        &             tend%    u_gwh(:,:,jb)   ,&
        &             tend%    v_gwh(:,:,jb) )

      IF (ltimer) call timer_stop(timer_gw_hines)

      ! heating
      zq_gwh(jcs:jce,:) = zdis_gwh(jcs:jce,:) * field%mair(jcs:jce,:,jb)

      ! heating accumulated
      zq_phy(jcs:jce,:) = zq_phy(jcs:jce,:) + zq_gwh(jcs:jce,:)

      ! tendency
      tend% temp_gwh(jcs:jce,:,jb) = zq_gwh(jcs:jce,:)*zconv(jcs:jce,:)

      ! tendencies accumulated
      tend% temp(jcs:jce,:,jb) = tend% temp(jcs:jce,:,jb) + tend% temp_gwh(jcs:jce,:,jb)
      tend%    u(jcs:jce,:,jb) = tend%    u(jcs:jce,:,jb) + tend%    u_gwh(jcs:jce,:,jb)
      tend%    v(jcs:jce,:,jb) = tend%    v(jcs:jce,:,jb) + tend%    v_gwh(jcs:jce,:,jb)

    ELSE ! NECESSARY COMPUTATIONS IF GW_HINES IS BY-PASSED

      tend% temp_gwh(jcs:jce,:,jb) = 0._wp
      tend%    u_gwh(jcs:jce,:,jb) = 0._wp
      tend%    v_gwh(jcs:jce,:,jb) = 0._wp

    END IF !lgw_hines


    ! 6.2   CALL SUBROUTINE SSODRAG

    IF (phy_config%lssodrag) THEN

      IF (ltimer) call timer_start(timer_ssodrag)

       CALL ssodrag( nc                                        ,& ! in,  number of cells/columns in loop (jce-jcs+1)
                     nbdim                                     ,& ! in,  dimension of block of cells/columns
                     nlev                                      ,& ! in,  number of levels
                     !
                     p_patch(jg)%cells%center(:,jb)%lat        ,& ! in,  Latitude in radians
                     psteplen                                  ,& ! in,  time step length, usually 2*delta_time
                     !
                     field% presi_old(:,:,jb)                  ,& ! in,  p at half levels
                     field% presm_old(:,:,jb)                  ,& ! in,  p at full levels
                     field% geom(:,:,jb)                       ,& ! in,  geopotential above surface (t-dt)
                     field% temp(:,:,jb)                       ,& ! in,  T
                     field%    u(:,:,jb)                       ,& ! in,  u
                     field%    v(:,:,jb)                       ,& ! in,  v
                     !
                     field% oromea(:,jb)                       ,& ! in,  Mean Orography (m)
                     field% orostd(:,jb)                       ,& ! in,  SSO standard deviation (m)
                     field% orosig(:,jb)                       ,& ! in,  SSO slope
                     field% orogam(:,jb)                       ,& ! in,  SSO Anisotropy
                     field% orothe(:,jb)                       ,& ! in,  SSO Angle
                     field% oropic(:,jb)                       ,& ! in,  SSO Peaks elevation (m)
                     field% oroval(:,jb)                       ,& ! in,  SSO Valleys elevation (m)
                     !
                     field% u_stress_sso(:,jb)                 ,& ! out, u-gravity wave stress
                     field% v_stress_sso(:,jb)                 ,& ! out, v-gravity wave stress
                     field% dissipation_sso(:,jb)              ,& ! out, dissipation by gravity wave drag
                     !
                     zdis_sso(:,:)                             ,& ! out, energy dissipation rate
                     tend%    u_sso(:,:,jb)                    ,& ! out, tendency of zonal wind
                     tend%    v_sso(:,:,jb)                     ) ! out, tendency of meridional wind

      IF (ltimer) call timer_stop(timer_ssodrag)

      ! heating
      zq_sso(jcs:jce,:) = zdis_sso(jcs:jce,:) * field%mair(jcs:jce,:,jb)

      ! heating accumulated
      zq_phy(jcs:jce,:) = zq_phy(jcs:jce,:) + zq_sso(jcs:jce,:)

      ! tendency
      tend% temp_sso(jcs:jce,:,jb) = zq_sso(jcs:jce,:)*zconv(jcs:jce,:)

      ! tendencies accumulated
      tend% temp(jcs:jce,:,jb) = tend% temp(jcs:jce,:,jb) + tend% temp_sso(jcs:jce,:,jb)
      tend%    u(jcs:jce,:,jb) = tend%    u(jcs:jce,:,jb) + tend%    u_sso(jcs:jce,:,jb)
      tend%    v(jcs:jce,:,jb) = tend%    v(jcs:jce,:,jb) + tend%    v_sso(jcs:jce,:,jb)

    ELSE ! NECESSARY COMPUTATIONS IF SSODRAG IS BY-PASSED

      tend% temp_sso(jcs:jce,:,jb) = 0._wp
      tend%    u_sso(jcs:jce,:,jb) = 0._wp
      tend%    v_sso(jcs:jce,:,jb) = 0._wp

    END IF ! SSODRAG

    !-------------------------------------------------------------------
    ! 7. CONVECTION PARAMETERISATION
    !-------------------------------------------------------------------
    itype(jcs:jce) = 0

    ! 7.1   INITIALIZE ARRAYS FOR CONVECTIVE PRECIPITATION
    !       AND COPY ARRAYS FOR CONVECTIVE CLOUD PARAMETERS

    tend% xl_dtr(jcs:jce,:,jb) = 0._wp
    tend% xi_dtr(jcs:jce,:,jb) = 0._wp
    zqtec  (jcs:jce,:) = 0._wp

    field% rsfc(:,jb) = 0._wp
    field% ssfc(:,jb) = 0._wp

    ! 7.2   CALL SUBROUTINE CUCALL FOR CUMULUS PARAMETERIZATION

    IF (phy_config%lconv) THEN

      IF (ltimer) call timer_start(timer_cucall)

      CALL cucall( jce, nbdim, nlev,          &! in
        &          nlevp1, nlevm1,            &! in
        &          ntrac,                     &! in     tracers
!        &          jb,                        &! in     row index
        &          psteplen,                  &! in
        &          field% lfland(:,jb),       &! in     loland
        &          field% temp(:,:,jb),       &! in     tm1
        &          field% u(:,:,jb),          &! in     um1
        &          field% v(:,:,jb),          &! in     vm1
        &          field% qtrc(:,:,jb,iqv),   &! in     qm1
        &          field% qtrc(:,:,jb,iqc),   &! in     xlm1
        &          field% qtrc(:,:,jb,iqi),   &! in     xim1
        &          field% qtrc(:,:,jb,iqt:),  &! in     xtm1
        &          tend% qtrc(:,:,jb,iqv),    &! in     qte  for internal updating
        &          tend% qtrc(:,:,jb,iqc),    &! in     xlte
        &          tend% qtrc(:,:,jb,iqi),    &! in     xite
        &          field% omega(:,:,jb),      &! in     vervel
        &          field% evap(:,jb),         &! in     qhfla (from "vdiff")
        &          field% geom(:,:,jb),       &! in     geom1
        &          field% presm_new(:,:,jb),  &! in     app1
        &          field% presi_new(:,:,jb),  &! in     aphp1
        &          field% thvsig(:,jb),       &! in           (from "vdiff")
        &          tend% temp(:,:,jb),        &! in     tte  for internal updating
        &          tend% u(:,:,jb),           &! in     vom  for internal updating
        &          tend% v(:,:,jb),           &! in     vol  for internal updating
        &          tend% qtrc(:,:,jb,iqt:),   &! in     xtte for internal updating
        &          zqtec,                     &! inout
        &          field% ch_concloud(:,jb),  &! inout condensational heat
        &          field% cw_concloud(:,jb),  &! inout condensational heat
        &          field% rsfc(:,jb),         &! out
        &          field% ssfc(:,jb),         &! out
        &          tend% xl_dtr(:,:,jb),      &! inout  xtecl
        &          tend% xi_dtr(:,:,jb),      &! inout  xteci
        &          itype,                     &! inout
        &          ictop,                     &! out
        &          ilab,                      &! out
        &          field% topmax(:,jb),       &! inout
        &          echam_conv_config%cevapcu, &! in
        &          zcd, zcv,                  &! in
        &          tend% qtrc_dyn(:,:,jb,iqv),&! in     qte by transport
        &          tend% qtrc_phy(:,:,jb,iqv),&! in     qte by physics
        &          field% con_dtrl(:,jb),     &! inout detrained liquid
        &          field% con_dtri(:,jb),     &! inout detrained ice
        &          field% con_iteqv(:,jb),    &! inout v. int. tend of water vapor within conv
        &          tend%temp_cnv(:,:,jb),     &! out
        &          tend%   u_cnv(:,:,jb),     &! out
        &          tend%   v_cnv(:,:,jb),     &! out
        &          tend%qtrc_cnv(:,:,jb,iqv), &! out
        &          tend%qtrc_cnv(:,:,jb,iqt:) )! out

      IF (ltimer) CALL timer_stop(timer_cucall)

      field% rtype(jcs:jce,jb) = REAL(itype(jcs:jce),wp)

      ! tendencies accumulated
      tend%    u(jcs:jce,:,jb)      = tend%    u(jcs:jce,:,jb)      + tend%    u_cnv(jcs:jce,:,jb)
      tend%    v(jcs:jce,:,jb)      = tend%    v(jcs:jce,:,jb)      + tend%    v_cnv(jcs:jce,:,jb)
      tend% temp(jcs:jce,:,jb)      = tend% temp(jcs:jce,:,jb)      + tend% temp_cnv(jcs:jce,:,jb)
      tend% qtrc(jcs:jce,:,jb,iqv)  = tend% qtrc(jcs:jce,:,jb,iqv)  + tend% qtrc_cnv(jcs:jce,:,jb,iqv)
      tend% qtrc(jcs:jce,:,jb,iqt:) = tend% qtrc(jcs:jce,:,jb,iqt:) + tend% qtrc_cnv(jcs:jce,:,jb,iqt:)


    ELSE ! NECESSARY COMPUTATIONS IF MASSFLUX IS BY-PASSED

      ilab(jcs:jce,1:nlev) = 0
      ictop(jcs:jce)       = nlev-1

      tend%    u_cnv(jcs:jce,:,jb)      = 0._wp
      tend%    v_cnv(jcs:jce,:,jb)      = 0._wp
      tend% temp_cnv(jcs:jce,:,jb)      = 0._wp
      tend% qtrc_cnv(jcs:jce,:,jb,iqv)  = 0._wp
      tend% qtrc_cnv(jcs:jce,:,jb,iqt:) = 0._wp

    ENDIF !lconv

    !-------------------------------------------------------------
    ! 7. LARGE SCALE CONDENSATION.
    !-------------------------------------------------------------
    IF(phy_config%lcond) THEN

      !IF (lcotra) CALL get_col_pol( tend%temp(:,:,jb),tend%qtrc(:,:,jb,iqv),jb )

      IF (ncdnc==0 .AND. nicnc==0) THEN

        field% rsfl(:,jb) = 0._wp
        field% ssfl(:,jb) = 0._wp

  CALL dbg_print('xl_dtr',tend%xl_dtr,str_module,idt_src,in_subset=p_patch(1)%cells%owned)
  CALL dbg_print('xi_dtr',tend%xi_dtr,str_module,idt_src,in_subset=p_patch(1)%cells%owned)

        IF (ltimer) CALL timer_start(timer_cloud)

        CALL cloud(jce, nbdim, jks, nlev, nlevp1, &! in
          &        psteplen,                  &! in
          &        ictop,                     &! in (from "cucall")
          &        field% presi_old(:,:,jb),  &! in
          &        field% omega (:,:,jb),     &! in. vervel
          &        field% presm_old(:,:,jb),  &! in
!          &        field% presm_new(:,:,jb), &! in
          &        field% acdnc (:,:,jb),     &! in. acdnc
          &        field% qtrc  (:,:,jb,iqv), &! in.  qm1
          &        field% temp  (:,:,jb),     &! in. tm1
          &        field%   tv  (:,:,jb),     &! in. ztvm1
          &        field% qtrc  (:,:,jb,iqc), &! in. xlm1
          &        field% qtrc  (:,:,jb,iqi), &! in. xim1
          &        zcair(:,:),                &! in
          &        field% geom  (:,:,jb),     &! in. geom1
          &        field% aclcov(:,  jb),     &! out
          &        field%  qvi  (:,  jb),     &! out
          &        field% xlvi  (:,  jb),     &! out
          &        field% xivi  (:,  jb),     &! out
          &        itype,                     &!
          &        field% ch_concloud(:,jb),  &! inout condens. heat
          &        field% cw_concloud(:,jb),  &! inout condens. heat
          &         tend% xl_dtr(:,:,jb),     &! inout  xtecl
          &         tend% xi_dtr(:,:,jb),     &! inout  xteci
          &        zqtec,                     &! inout (there is a clip inside)
          &         tend% qtrc  (:,:,jb,iqv), &! inout.  qte
          &         tend% temp  (:,:,jb),     &! inout.  tte
          &         tend% qtrc  (:,:,jb,iqc), &! inout. xlte
          &         tend% qtrc  (:,:,jb,iqi), &! inout. xite
          &        field% cld_dtrl(:,jb),     &! inout detrained liquid
          &        field% cld_dtri(:,jb),     &! inout detrained ice
          &        field% cld_iteq(:,jb),     &! inout v. int. tend of qv,qc, and qi within cloud
!          &         tend% x_dtr(:,:,jb),      &! inout (there is a clip inside)
          &        field% aclc  (:,:,jb),     &! inout
          &        field% ssfl  (:,  jb),     &! out
          &        field% rsfl  (:,  jb),     &! out
          &        field% relhum(:,:,jb),     &! out
          &        tend%temp_cld(:,:,jb),     &! out
          &        tend%qtrc_cld(:,:,jb,iqv), &! out
          &        tend%qtrc_cld(:,:,jb,iqc), &! out
          &        tend%qtrc_cld(:,:,jb,iqi)  )! out

        IF (ltimer) CALL timer_stop(timer_cloud)

      ELSE IF (ncdnc>0 .AND. nicnc>0) THEN
!0      CALL cloud_cdnc_icnc(...) !!skipped in ICON
      ELSE
        IF (my_process_is_stdio()) CALL finish('echam_phy_main', ' check setting of ncdnc and nicnc.')
      END IF

    ELSE ! NECESSARY COMPUTATIONS IF *CLOUD* IS BY-PASSED.

      field% rsfl (jcs:jce,  jb) = 0._wp
      field% ssfl (jcs:jce,  jb) = 0._wp
      field% aclc (jcs:jce,:,jb) = 0._wp

      tend% temp_cld(jcs:jce,:,jb)      = 0._wp
      tend% qtrc_cld(jcs:jce,:,jb,iqv)  = 0._wp
      tend% qtrc_cld(jcs:jce,:,jb,iqc)  = 0._wp
      tend% qtrc_cld(jcs:jce,:,jb,iqi)  = 0._wp
      tend% qtrc_cld(jcs:jce,:,jb,iqt:) = 0._wp

    ENDIF !lcond

    ! KF accumulate fields for diagnostics

    !  total precipitation flux
       field% totprec (jcs:jce,jb)     =  field% rsfl (jcs:jce,jb) & ! rain large scale
            &                            +field% ssfl (jcs:jce,jb) & ! snow large scale
            &                            +field% rsfc (jcs:jce,jb) & ! rain convection
            &                            +field% ssfc (jcs:jce,jb)   ! snow convection

    ! accumulated total precipitation flux => average when output
       field% totprec_avg (jcs:jce,jb) =  field% totprec_avg (jcs:jce,jb)          &
            &                            +field% totprec     (jcs:jce,jb) * pdtime

    ! Now compute tendencies from physics alone

    tend% temp_phy (jcs:jce,:,jb)   = tend% temp (jcs:jce,:,jb)   - tend% temp_phy (jcs:jce,:,jb)
    tend%    u_phy (jcs:jce,:,jb)   = tend%    u (jcs:jce,:,jb)   - tend%    u_phy (jcs:jce,:,jb)
    tend%    v_phy (jcs:jce,:,jb)   = tend%    v (jcs:jce,:,jb)   - tend%    v_phy (jcs:jce,:,jb)
    tend% qtrc_phy (jcs:jce,:,jb,:) = tend% qtrc (jcs:jce,:,jb,:) - tend% qtrc_phy (jcs:jce,:,jb,:)

    IF ( iequations == inh_atmosphere ) THEN
      tend% temp_phy (jcs:jce,:,jb) = tend% temp_phy(jcs:jce,:,jb)*zcpair(jcs:jce,:)/zcvair(jcs:jce,:)
    END IF

    ! Done. Disassociate pointers.
    NULLIFY(field,tend)

  END SUBROUTINE echam_phy_main
  !-------------

END MODULE mo_echam_phy_main
