!>
!! @brief Subroutine physc calls all the parameterization schemes
!!
!! @author Hui Wan, MPI-M
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version from ECHAM6 (revision 2028)
!!  Modified for ICOHAM by Hui Wan and Marco Giorgetta (2010)
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
#ifdef __xlC__
@PROCESS HOT
@PROCESS SPILLSIZE(5000)
#endif
!OCL NOALIAS

MODULE mo_echam_phy_main

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_ext_data,            ONLY: ext_data
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_math_constants,      ONLY: pi
  USE mo_physical_constants,  ONLY: grav
  USE mo_run_config,          ONLY: ntracer, nlev, nlevp1, ltestcase, &
    &                               iqv, iqc, iqi, io3, iqt, ltimer
  USE mo_ha_testcases,        ONLY: ctest_name
  USE mo_vertical_coord_table,ONLY: nlevm1
  USE mo_echam_phy_config,    ONLY: echam_phy_config
  USE mo_echam_conv_config,   ONLY: echam_conv_config
  USE mo_cucall,              ONLY: cucall
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field, prm_field,     &
    &                               t_echam_phy_tend,  prm_tend
  USE mo_timer,               ONLY: timer_start, timer_stop,          &
    &                               timer_cover, timer_cloud,         &
    &                               timer_radheat,                    &
    &                               timer_cucall, timer_vdiff,        &
    &                               timer_gw_hines
  USE mo_ham_aerosol_params,  ONLY: ncdnc, nicnc
  USE mo_icoham_sfc_indices,  ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_cloud,               ONLY: cloud
  USE mo_cover,               ONLY: cover
  USE mo_echam_cloud_params,  ONLY: ctaus, ctaul, ctauk !, ncctop, nccbot
  USE mo_radiation,           ONLY: radiation, radheat
  USE mo_radiation_config,    ONLY: tsi, izenith
  USE mo_vdiff_config,        ONLY: vdiff_config
  USE mo_vdiff_downward_sweep,ONLY: vdiff_down
  USE mo_vdiff_upward_sweep,  ONLY: vdiff_up
  USE mo_vdiff_solver,        ONLY: nvar_vdiff, nmatrix,        &
                                  & ih_vdiff=>ih, iqv_vdiff=>iqv
  USE mo_gw_hines,            ONLY: gw_hines
  ! provisional to get coordinates
  USE mo_atmo_control,        ONLY: p_patch

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: physc

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !>
  !!
  SUBROUTINE physc( jg,jb,jcs,jce,nbdim,pdtime,psteplen,  &
                  & ltrig_rad,ptime_radtran,ptime_radheat )

    INTEGER, INTENT(IN) :: jg             !< grid level/domain index 
    INTEGER, INTENT(IN) :: jb             !< block index
    INTEGER, INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER, INTENT(IN) :: nbdim          !< size of this block
    REAL(wp),INTENT(IN) :: pdtime         !< time step
    REAL(wp),INTENT(IN) :: psteplen       !< 2*pdtime in case of leapfrog

    LOGICAL, INTENT(IN) :: ltrig_rad      !< perform radiative transfer computation
    REAL(wp),INTENT(IN) :: ptime_radtran  !< time instance of the radiative transfer
                                          !< computation, scaled into radians
    REAL(wp),INTENT(IN) :: ptime_radheat  !< time instance of the radiative heating 
                                          !< computation, scaled into radians
    ! Local variables

    TYPE(t_echam_phy_field),POINTER :: field 
    TYPE(t_echam_phy_tend) ,POINTER :: tend

    REAL(wp) :: zlat_deg(nbdim)           !< latitude in deg N

    REAL(wp) :: zbetaa (nbdim,nlev)       !< qt distribution minimum in beta
    REAL(wp) :: zbetab (nbdim,nlev)       !< qt distribution maximum in beta
    REAL(wp) :: zbetass(nbdim,nlev)

    REAL(wp) :: zhmixtau   (nbdim,nlev)   !< timescale of mixing for horizontal eddies
    REAL(wp) :: zvmixtau   (nbdim,nlev)   !< timescale of mixing for vertical turbulence
    REAL(wp) :: zqtvar_prod(nbdim,nlev)   !< production rate of total water variance
                                          !< due to turbulence. Computed in "vdiff",
                                          !< used by "cloud"
    INTEGER  :: itype(nbdim)              !< type of convection
    INTEGER  :: invb (nbdim)

    REAL(wp) :: zfrl (nbdim)              !< fraction of land in the grid box
    REAL(wp) :: zfrw (nbdim)              !< fraction of water (without ice) in the grid point
    REAL(wp) :: zfri (nbdim)              !< fraction of ice in the grid box
    REAL(wp) :: zfrc (nbdim,nsfc_type)    !< zfrl, zfrw, zfrc combined

    REAL(wp) :: zqhflx (nbdim)
    INTEGER  :: ilab   (nbdim,nlev)
    REAL(wp) :: zcvcbot(nbdim)
    REAL(wp) :: zwcape (nbdim)

    REAL(wp) :: zqtec  (nbdim,nlev)       !< tracer tendency due to entrainment/detrainment
    REAL(wp) :: zxtecl (nbdim,nlev)       !< tracer tendency due to entrainment/detrainment
    REAL(wp) :: zxteci (nbdim,nlev)       !< tracer tendency due to entrainment/detrainment
    REAL(wp) :: zxtecnl(nbdim,nlev)       !< tracer tendency due to entrainment/detrainment
    REAL(wp) :: zxtecni(nbdim,nlev)       !< tracer tendency due to entrainment/detrainment

    REAL(wp) :: ztsi                      !< total solar irradiation at 1 AU   [W/m2]
    REAL(wp) :: zi0    (nbdim)            !< solar incoming radiation at TOA   [W/m2]
    REAL(wp) :: zmair  (nbdim,nlev)       !< mass of air                       [kg/m2]
    REAL(wp) :: zdelp  (nbdim,nlev)       !< layer thickness in pressure coordinate  [Pa]

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

    REAL(wp) :: zaa    (nbdim,nlev,3,nmatrix)   !< coeff. matrices, all variables
    REAL(wp) :: zaa_btm(nbdim,3,nsfc_type)      !< last row of coeff. matrix of heat and moisture
    REAL(wp) :: zbb    (nbdim,nlev,nvar_vdiff)  !< r.h.s., all variables
    REAL(wp) :: zbb_btm(nbdim,nsfc_type,ih_vdiff:iqv_vdiff) !< last row of r.h.s. of heat and moisture

    ! Temporary arrays used by VDIFF

    REAL(wp) :: zprfac_sfc(nbdim)
    REAL(wp) :: zcpt_sfc_tile(nbdim,nsfc_type)  !< dry static energy at surface

    REAL(wp) :: zcptgz  (nbdim,nlev) !< dry static energy
    REAL(wp) :: zrhoh   (nbdim,nlev) !< air density at half levels
    REAL(wp) :: zqshear (nbdim,nlev) !<
    REAL(wp) :: zthvvar (nbdim,nlev) !< intermediate value of thvvar
    REAL(wp) :: ztkevn  (nbdim,nlev) !< intermediate value of tke

!!$ REAL(wp) :: rlfland (nbdim), rlfglac (nbdim)

    ! number of cells/columns from index jcs to jce
    nc = jce-jcs+1

    ! start index for vertical loops
    jks=1

    ! 1. Associate pointers

    field => prm_field(jg)
    tend  => prm_tend (jg)

    ! 2. local switches and parameters

    ntrac = ntracer-iqt+1  !# of tracers excluding water vapour and hydrometeors

    !------------------------------------------------------------
    ! 3. COMPUTE SOME FIELDS NEEDED BY THE PHYSICAL ROUTINES.
    !------------------------------------------------------------
    DO jk = 1,nlev
      DO jc = jcs,jce
        !
        ! 3.1 Timescale of mixing for horizontal eddies
        !
        zhmixtau(jc,jk) = ctauk*ABS( field%vor(jc,jk,jb) )
        zhmixtau(jc,jk) = MIN( ctaus,MAX(ctaul,zhmixtau(jc,jk)) )
        !
        ! 3.2 Thickness of model layer in pressure coordinate; mass of air
        !
        zdelp   (jc,jk) = field% presi_old (jc,jk+1,jb) - field% presi_old (jc,jk,jb)
        zmair   (jc,jk) = zdelp(jc,jk)/grav
      END DO
    END DO

    ! 3.3 Weighting factors for fractional surface coverage
    !     Accumulate ice portion for diagnostics

    DO jc=jcs,jce

      ! fraction of land in the grid box. lsmask: land-sea mask, 1.= land

      zfrl(jc) = field% lsmask(jc,jb)

      ! fraction of sea/lake in the grid box
      ! * (1. - fraction of sea ice in the sea/lake part of the grid box)
      ! => fraction of open water in the grid box

      zfrw(jc) = (1._wp-zfrl(jc))*(1._wp-field%seaice(jc,jb))

      ! fraction of sea ice in the grid box
      zfri(jc) = 1._wp-zfrl(jc)-zfrw(jc)
    END DO

    ! 3.4 Merge three pieces of information into one array for vdiff
    IF (echam_phy_config%lvdiff) THEN
      IF (ilnd.LE.nsfc_type) zfrc(jcs:jce,ilnd) = zfrl(jcs:jce)
      IF (iwtr.LE.nsfc_type) zfrc(jcs:jce,iwtr) = zfrw(jcs:jce)
      IF (iice.LE.nsfc_type) zfrc(jcs:jce,iice) = zfri(jcs:jce)
    ENDIF

    !!3.9 DETERMINE TROPOPAUSE HEIGHT AND MASS BUDGETS--------------------
    !     (Needed only for sub-models. Note: sequence of arguments
    !      different from the original ECHAM6)
    !
    !CALL WMO_tropopause( jce, nbdim, nlev,         &! in
    !                   & ncctop, nccbot, lresum,   &! in
    !                   & field% temp(:,:,jb),      &! in
    !                   & field% presm_old(:,:,jb), &! in
    !                   & field% tropo(:,jb),       &! out for diagnostics
    !                   & itrpwmo, itrpwmop1        )! out for submodel
    !---------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! 3.13 DIAGNOSE CURRENT CLOUD COVER
    !-------------------------------------------------------------------
    itype(jcs:jce) = NINT(field%rtype(jcs:jce,jb))

    IF (echam_phy_config%lcond) THEN
      IF (ltimer) CALL timer_start(timer_cover)

      CALL cover( jce, nbdim, jks,          &! in
        &         nlev, nlevp1,             &! in
        &         echam_phy_config%lcover,  &! in
        &         itype,  zfrw,             &! in
        &         field% presi_old(:,:,jb), &! in
        &         field% presm_old(:,:,jb), &! in
        &         field% omega(:,:,jb),     &! in    vervel
        &         field%  temp(:,:,jb),     &! in    tm1
        &         field%     q(:,:,jb,iqv), &! in    qm1
        &         field%     q(:,:,jb,iqc), &! in    xlm1
        &         field%     q(:,:,jb,iqi), &! in    xim1
        &         field%  aclc(:,:,jb),     &! inout
        &         field%  xvar(:,:,jb),     &! inout (for "vdiff"+"cloud")
        &         field% xskew (:,:,jb),    &! inout (for "cloud")
        &         invb,                     &! out   (for "cloud")
        &         field% rintop(:,  jb),    &! out   (for output)
        &         zbetaa, zbetab, zbetass )  ! out   (for "cloud")

      IF (ltimer) CALL timer_stop(timer_cover)
    ENDIF ! lcond

    !-------------------------------------------------------------------
    ! 4. RADIATION PARAMETERISATION
    !-------------------------------------------------------------------
    IF (echam_phy_config%lrad) THEN

       SELECT CASE(izenith)  
       CASE(0)
       ! for testing: provisional setting of cos(zenith angle) and TSI
       ! The global mean insolation is TSI/4 (ca. 340 W/m2)
         ztsi = tsi/4._wp ! scale ztsi by 1/4 to get the correct global mean insolation
       CASE(1)
       ! circular non-seasonal orbit, zenith angle dependent on latitude only,
       ! no diurnal cycle (always at 12:00 local time --> sin(time of day)=1 )
         ztsi = tsi/pi ! because sun is always in local noon, the TSI needs to be
       !               ! scaled by 1/pi to get the correct global mean insolation
       CASE(2)
       ! circular non-seasonal orbit, no diurnal cycle
       ! at 07:14:15 or 16:45:45 local time (--> sin(time of day)=1/pi )
         ztsi = tsi
       CASE(3) 
       ! circular non-seasonal orbit, with diurnal cycle
         ztsi = tsi
       END SELECT

       ! 4.1 RADIATIVE TRANSFER
       !-----------------------
       IF (ltrig_rad) THEN

          ! to do (for implementing seasonal cycle):
          ! - compute orbit position at ptime_radtran

          SELECT CASE(izenith)  
          CASE(0)
          ! for testing: provisional setting of cos(zenith angle) and TSI
          ! The global mean insolation is TSI/4 (ca. 340 W/m2)

            field%cosmu0(jcs:jce,jb) = 1._wp ! sun in zenith everywhere

          CASE(1)
          ! circular non-seasonal orbit, zenith angle dependent on latitude only,
          ! no diurnal cycle (always at 12:00 local time --> sin(time of day)=1 )

            field%cosmu0(jcs:jce,jb) = COS( p_patch(jg)%cells%center(jcs:jce,jb)%lat )
           
          CASE(2)
          ! circular non-seasonal orbit, no diurnal cycle
          ! at 07:14:15 or 16:45:45 local time (--> sin(time of day)=1/pi )

            field%cosmu0(jcs:jce,jb) = COS( p_patch(jg)%cells%center(jcs:jce,jb)%lat )/pi

          CASE(3) 
          ! circular non-seasonal orbit, with diurnal cycle

            field%cosmu0(jcs:jce,jb) = -COS( p_patch(jg)%cells%center(jcs:jce,jb)%lat ) &
                                     & *COS( p_patch(jg)%cells%center(jcs:jce,jb)%lon   &
                                     &      +ptime_radtran )

          END SELECT

          field% albvisdir(jcs:jce,jb) = 0.07_wp ! ~ albedo of water
          field% albnirdir(jcs:jce,jb) = 0.07_wp ! ~ albedo of water
          field% albvisdif(jcs:jce,jb) = 0.07_wp ! ~ albedo of water
          field% albnirdif(jcs:jce,jb) = 0.07_wp ! ~ albedo of water

          
          !  as a preparation step:
          !  bring ozone climatology via time weigths onto a 3D field 
!
!        DO jl=1,nproma
!          z_o3_3d(jcs:jce,:,jb) = wgt1_m*ext_data(jg)%atm_td_O3(jcs:jce,:,nmw1_m)&
!                                 & +wgt2_m*ext_data(jg)%atm_td_O3ozone(jcs:jce,:nmw2_m)
!        END DO


!!$        ! debug fields "radin"
!!$        !
!!$        DO jc = jcs,jce
!!$          !
!!$          IF (field% lfland(jc,jb)) THEN
!!$            rlfland(jc) = 1._wp
!!$          ELSE
!!$            rlfland(jc) = 0._wp
!!$          END IF
!!$          !
!!$          IF (field% lfglac(jc,jb)) THEN
!!$            rlfglac(jc) = 1._wp
!!$          ELSE
!!$            rlfglac(jc) = 0._wp
!!$          END IF
!!$          !
!!$        END DO
!!$        !
!!$        field% debug_2d_1(:,jb) = field% albvisdir(:,jb)
!!$        field% debug_2d_2(:,jb) = field% albnirdir(:,jb)
!!$        field% debug_2d_3(:,jb) = field% albvisdif(:,jb)
!!$        field% debug_2d_4(:,jb) = field% albnirdif(:,jb)
!!$        field% debug_2d_5(:,jb) = REAL(itype(:),wp)
!!$        field% debug_2d_6(:,jb) = field% tsfc(:,jb)
!!$        field% debug_2d_7(:,jb) = rlfland(:)
!!$        field% debug_2d_8(:,jb) = rlfglac(:)
!!$        !
!!$        field% debug_3d_1(:,1:nlev,jb) = field% presi_old(:,1:nlev,jb)
!!$        field% debug_3d_2(:,1:nlev,jb) = field% presm_old(:,1:nlev,jb)
!!$        field% debug_3d_3(:,1:nlev,jb) = field% temp     (:,1:nlev,jb)
!!$        field% debug_3d_4(:,1:nlev,jb) = field% q        (:,1:nlev,jb,iqv)
!!$        field% debug_3d_5(:,1:nlev,jb) = field% q        (:,1:nlev,jb,iqc)
!!$        field% debug_3d_6(:,1:nlev,jb) = field% q        (:,1:nlev,jb,iqi)
!!$        field% debug_3d_7(:,1:nlev,jb) = field% acdnc    (:,1:nlev,jb)
!!$        field% debug_3d_8(:,1:nlev,jb) = field% aclc     (:,1:nlev,jb)

!!        IF (ltimer) CALL timer_start(timer_radiation)

        CALL radiation(               &
          !
          ! argument                   !  INTENT comment
          !
          ! input
          ! -----
          !
          ! indices and dimensions
!!$          & jg                       ,&!< in     domain index
!!$          & jb                       ,&!< in     block index
          & jce                      ,&!< in     end   index for loop over block
          & nbdim                    ,&!< in     dimension of block over cells
          & nlev                     ,&!< in     number of full levels = number of layers
          & nlevp1                   ,&!< in     number of half levels = number of layer interfaces
          !
          & itype(:)                 ,&!< in     type of convection
          !
          ! surface: albedo + temperature
          & field% lsmask(:,jb)      ,&!< in     land-sea mask. (1. = land, 0. = sea/lakes)
          & field% glac(:,jb)        ,&!< in     fraction of land covered by glaciers
          !
          & field% cosmu0(:,jb)      ,&!< in     cos of zenith angle mu0 for rad. transfer calc.
          & field% albvisdir(:,jb)   ,&!< in     surface albedo for visible range, direct
          & field% albnirdir(:,jb)   ,&!< in     surface albedo for near IR range, direct
          & field% albvisdif(:,jb)   ,&!< in     surface albedo for visible range, diffuse
          & field% albnirdif(:,jb)   ,&!< in     surface albedo for near IR range, diffuse
          & ext_data(jg)%atm%emis_rad(:,jb), & !< in longwave surface emissivity
          & field% tsfc(:,jb)        ,&!< in     grid box mean surface temperature
          !
          ! atmopshere: pressure, tracer mixing ratios and temperature
          & field% presi_old(:,:,jb) ,&!< in     pressure at half levels at t-dt [Pa]
          & field% presm_old(:,:,jb) ,&!< in     pressure at full levels at t-dt [Pa]
          & field% temp (:,:,jb)     ,&!< in     tk_fl  = temperature at full level at t-dt
          & field% q(:,:,jb,iqv)     ,&!< in     qm_vap = water vapor mass mixing ratio at t-dt
          & field% q(:,:,jb,iqc)     ,&!< in     qm_liq = cloud water mass mixing ratio at t-dt
          & field% q(:,:,jb,iqi)     ,&!< in     qm_ice = cloud ice mass mixing ratio at t-dt
          & field% q(:,:,jb,io3)     ,&!< in     qm_o3 = o3 mass mixing ratio at t-dt
!!$       & field% geom(:,:,jb)     ,&!< in     pgeom1 = geopotential above ground at t-dt [m2/s2]
          & field% acdnc(:,:,jb)     ,&!< in     cld_frac = cloud fraction [m2/m2]
          & field% aclc(:,:,jb)      ,&!< in     cld_frac = cloud fraction [m2/m2]
          !
          ! output
          ! ------
          !
          & field% aclcov(:,jb)      ,&!< out    cloud cover in a column [m2/m2]
          !
!!$          ! surface shortwave
!!$          ! - transmissivities in spectral range w.r.t. total solar irradiation
!!$          ! - diffuse fractions in a spectral range w.r.t. total flux in the same spectral range
!!$          & field% nirsfc(:,jb)      ,&!< out    solar transmissivity in NIR, net downward
!!$          & field% nirdffsfc(:,jb)   ,&!< out    diffuse fraction in NIR net downw. flux
!!$          & field% vissfc(:,jb)      ,&!< out    solar transmissivity in VIS, net downward
!!$          & field% visdffsfc(:,jb)   ,&!< out    diffuse fraction in VIS net downw. flux
!!$          & field% parsfc(:,jb)      ,&!< out    solar transmissivity in PAR, downward
!!$          & field% pardffsfc(:,jb)   ,&!< out    diffuse fraction in PAR net downw. flux
          !
          ! atmospheric profiles
          & field% emterclr(:,:,jb)  ,&!< out    terrestrial emissivity, clear sky, net downward
          & field% trsolclr(:,:,jb)  ,&!< out    solar transmissivity  , clear sky, net downward
          !
          & field% emterall(:,:,jb)  ,&!< out    terrestrial flux      , all   sky, net downward
          & field% trsolall(:,:,jb)   )!< out    solar transmissivity  , all   sky, net downward

!!        IF (ltimer) CALL timer_stop(timer_radiation)

!!$        ! debug fields "radout"
!!$        !
!!$        field% debug_2d_1(:,jb) = field% aclcov(:,jb)
!!$        field% debug_2d_2(:,jb) = field% nirsfc(:,jb)
!!$        field% debug_2d_3(:,jb) = field% nirdffsfc(:,jb)
!!$        field% debug_2d_4(:,jb) = field% vissfc(:,jb)
!!$        field% debug_2d_5(:,jb) = field% visdffsfc(:,jb)
!!$        field% debug_2d_6(:,jb) = field% parsfc(:,jb)
!!$        field% debug_2d_7(:,jb) = field% pardffsfc(:,jb)
!!$        !
!!$        field% debug_3d_1(:,1:nlev,jb) = field% emterclr(:,1:nlev,jb)
!!$        field% debug_3d_2(:,1:nlev,jb) = field% trsolclr(:,1:nlev,jb)
!!$        field% debug_3d_3(:,1:nlev,jb) = field% emterall(:,1:nlev,jb)
!!$        field% debug_3d_4(:,1:nlev,jb) = field% trsolall(:,1:nlev,jb)
!!$
         END IF ! ltrig_rad
!!$
      ! 4.2 RADIATIVE HEATING
      !----------------------

      ! to do:
      ! - compute orbit position at ptime_radheat

      ! - solar incoming flux at TOA
      field% cosmu0(jcs:jce,jb) = -COS( p_patch(jg)%cells%center(jcs:jce,jb)%lat ) &
                                & *COS( p_patch(jg)%cells%center(jcs:jce,jb)%lon   &
                                &      +ptime_radheat )


      zi0(jcs:jce) = MAX(0._wp,field%cosmu0(jcs:jce,jb)) * ztsi  ! instantaneous for radheat
      field% flxdwswtoa(jcs:jce,jb) = zi0 (jcs:jce)               ! (to be accumulated for output)

      IF (ltimer) CALL timer_start(timer_radheat)

      CALL radheat (                   &
        !
        ! input
        ! -----
        !
        & jcs, jce ,nbdim             ,&! in     loop start/end indices and dimension size
        & nlev  ,nlevp1               ,&! in     vertical dimension size
        & zmair           (:,:)       ,&! in     layer air mass                [kg/m2]
        & field% q(:,:,jb,iqv)        ,&! in     specific moisture             [kg/kg]
        & zi0             (:)         ,&! in     solar incoming flux at TOA    [W/m2]
        & ext_data(jg)%atm%emis_rad(:,jb), & !in lw sfc emissivity
        & field% tsfc     (:,jb)      ,&! in     surface temperature           [K]
        & field% tsfc     (:,jb)      ,&! in     sfc temp. used in "radiation" [K]
        & field% trsolall (:,:,jb)    ,&! in     shortwave net tranmissivity   []
        & field% emterall (:,:,jb)    ,&! in     longwave net flux             [W/m2]
        !
        ! output
        ! ------
        !
        & tend%temp_radsw (:,:,jb)    ,&! out    rad. heating by SW         [K/s]
        & tend%temp_radlw (:,:,jb)     )! out    rad. heating by LW         [K/s]

      IF (ltimer) CALL timer_stop(timer_radheat)

      ! Add shortwave and longwave heating rate to total heating rate
      tend% temp(jcs:jce,:,jb) = tend% temp       (jcs:jce,:,jb) &
        &                      + tend% temp_radsw (jcs:jce,:,jb) &
        &                      + tend% temp_radlw (jcs:jce,:,jb)

    ELSE   ! If computation of radiative heating is by-passed

      tend%temp_radsw(jcs:jce,:,jb) = 0.0_wp
      tend%temp_radlw(jcs:jce,:,jb) = 0.0_wp

    END IF ! lrad


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

    IF (echam_phy_config%lvdiff) THEN
      IF (ltimer) CALL timer_start(timer_vdiff)

      CALL vdiff_down( vdiff_config%lsfc_mom_flux,      &! in
                     & vdiff_config%lsfc_heat_flux,     &! in
                     & jce, nbdim, nlev, nlevm1, nlevp1,&! in
                     & ntrac, nsfc_type,                &! in
                     & iwtr, iice, ilnd,                &! in, indices of different surface types
                     & psteplen,                        &! in, time step (2*dt if leapfrog)
                     & field%coriol(:,jb),              &! in, Coriolis parameter
                     & zfrc(:,:),                       &! in, area fraction of each sfc type
                    !& field% tsfc_tile(:,:,jb),        &! in, surface temperature
                     & field% tsfc_tile(:,jb,:),        &! in, surface temperature
                     & field% ocu (:,jb),               &! in, ocean sfc velocity, u-component
                     & field% ocv (:,jb),               &! in, ocean sfc velocity, v-component
                     & field% presi_old(:,nlevp1,jb),   &! in, sfc pressure
                     & field%    u(:,:,jb),             &! in, um1
                     & field%    v(:,:,jb),             &! in, vm1
                     & field% temp(:,:,jb),             &! in, tm1
                     & field%    q(:,:,jb,iqv),         &! in, qm1
                     & field%    q(:,:,jb,iqc),         &! in, xlm1
                     & field%    q(:,:,jb,iqi),         &! in, xim1
                     & field%   qx(:,:,jb),             &! in, xlm1 + xim1
                     & field%    q(:,:,jb,iqt:),        &! in, xtm1
                     & field% presi_old(:,:,jb),        &! in, aphm1
                     & field% presm_old(:,:,jb),        &! in, apm1
                     & zdelp(:,:),                      &! in, layer thickness [Pa]
                     & field% geom(:,:,jb),             &! in, pgeom1 = geopotential above ground
                     & field%   tv(:,:,jb),             &! in, virtual temperaturea
                     & field% aclc(:,:,jb),             &! in, cloud fraction
                     & zxt_emis,                        &! in, zxtems
                     & field% thvvar(:,:,jb),           &! in, variance of theta_v at step t-dt
                     & field%   xvar(:,:,jb),           &! in
                    !& field% z0m_tile(:,:,jb),         &! in
                     & field% z0m_tile(:,jb,:),         &! in
                     & field%  tkem1(:,:,jb),           &! in, TKE at step t-dt
                     & field%  ustar(:,  jb),           &! inout
                    !& field% qs_sfc_tile(:,:,jb),      &! out, sfc specific humidity at saturation
                     & field% qs_sfc_tile(:,jb,:),      &! out, sfc specific humidity at saturation
                     & ihpbl(:),                        &! out, for "vdiff_up"
                     & field%    ghpbl(:,jb),           &! out, for output
                     & field%      ri (:,:,jb),         &! out, for output
                     & field%  mixlen (:,:,jb),         &! out, for output
                     & field% cfm     (:,:,jb),         &! out, for output
                    !& field% cfm_tile(:,:,jb),         &! out, for output and "vdiff_up"
                     & field% cfm_tile(:,jb,:),         &! out, for output and "vdiff_up"
                     & field% cfh     (:,:,jb),         &! out, for output
                    !& field% cfh_tile(:,:,jb),         &! out, for output and "vdiff_up"
                     & field% cfh_tile(:,jb,:),         &! out, for output and "vdiff_up"
                     & field% cfv     (:,:,jb),         &! out, for output
                     & field% cftke   (:,:,jb),         &! out, for output
                     & field% cfthv   (:,:,jb),         &! out, for output
                     & zaa, zaa_btm, zbb, zbb_btm,      &! out, for "vdiff_up"
                     & zprfac_sfc(:),                   &! out, for "vdiff_up"
                     & zcpt_sfc_tile(:,:),              &! out, for "vdiff_up"
                     & zcptgz(:,:), zrhoh(:,:),         &! out, for "vdiff_up"
                     & zqshear(:,:),                    &! out, for "vdiff_up"
                     & zthvvar(:,:),                    &! out, for "vdiff_up"
                     & ztkevn (:,:)                     )! out, for "vdiff_up"

      IF (ltimer) CALL timer_stop(timer_vdiff)
    END IF !lvdiff

    ! 5.4 Surface processes that provide time-dependent lower boundary
    !     condition for wind, temperature, tracer concentraion, etc.

    IF (ltestcase)THEN
      SELECT CASE (ctest_name)
      CASE('APE','JWw-Moist','LDF-Moist')
      ! Surface temperature is fixed; Surface wind is set to zero
      ! during model initialization; Nothing to be done here.
        CONTINUE
      CASE DEFAULT
        CALL finish('physc','Unknown choice of ctest_name')
      END SELECT
    ELSE
      ! Get new SST by interpolating the AMIP data or calling the
      ! ocean model; Call land and sea ice models.
      CALL finish('physc','surface processes not implemented')
    ENDIF

    ! 5.5 Turbulent mixing, part II:
    !     - Elimination for the lowest model level using boundary conditions
    !       provided by the surface model(s);
    !     - Back substitution to get solution of the tridiagonal system;
    !     - Compute tendencies and additional diagnostics.

    IF (echam_phy_config%lvdiff) THEN
      IF (ltimer) CALL timer_start(timer_vdiff)

      CALL vdiff_up( vdiff_config%lsfc_heat_flux,     &! in
                   & jce, nbdim, nlev, nlevm1, nlevp1,&! in
                   & ntrac, nsfc_type,                &! in
                   & iwtr, iice, ilnd,                &! in, indices of different sfc types
                   & pdtime, psteplen,                &! in, time steps
                   & zfrc(:,:),                       &! in, area fraction of each sfc type
                   & field% ocu (:,jb),               &! in, ocean sfc velocity, u-component
                   & field% ocv (:,jb),               &! in, ocean sfc velocity, v-component
                  !& field% cfm_tile(:,:,jb),         &! in
                  !& field% cfh_tile(:,:,jb),         &! in
                   & field% cfm_tile(:,jb,:),         &! in
                   & field% cfh_tile(:,jb,:),         &! in
                  !& field% qs_sfc_tile(:,:,jb),      &! in, sfc spec. humidity at saturation
                   & field% qs_sfc_tile(:,jb,:),      &! in, sfc spec. humidity at saturation
                   & zaa, zaa_btm, zprfac_sfc(:),     &! in, from "vdiff_down"
                   & zcpt_sfc_tile(:,:),              &! in, from "vdiff"down"
                   &   ihpbl(:),                      &! in, from "vdiff_down"
                   &  zcptgz(:,:),                    &! in, from "vdiff_down"
                   &   zrhoh(:,:),                    &! in, from "vdiff_down"
                   & zqshear(:,:),                    &! in, from "vdiff_down"
                   & field%    u(:,:,jb),             &! in, um1
                   & field%    v(:,:,jb),             &! in, vm1
                   & field% temp(:,:,jb),             &! in, tm1
                   & field%    q(:,:,jb,iqv),         &! in, qm1
                   & field%    q(:,:,jb,iqc),         &! in, xlm1
                   & field%    q(:,:,jb,iqi),         &! in, xim1
                   & field%    q(:,:,jb,iqt:),        &! in, xtm1
                   & zdelp(:,:),                      &! in, layer thickness [Pa]
                   & field% geom(:,:,jb),             &! in, pgeom1 = geopotential above ground
                   &      ztkevn(:,:),                &! in, tke at intermediate time step
                   & field%tkem1(:,:,jb),             &! in, TKE at step t-dt
                   & zbb, zbb_btm,                    &! inout
                   & zthvvar(:,:),                    &! inout
                   & field%   xvar(:,:,jb),           &! inout
                  !& field% z0m_tile(:,:,jb),         &! inout
                   & field% z0m_tile(:,jb,:),         &! inout
                   & field% kedisp(:,  jb),           &! inout, "vdis" in ECHAM
                   &  tend%    u(:,:,jb),             &! inout
                   &  tend%    v(:,:,jb),             &! inout
                   &  tend% temp(:,:,jb),             &! inout
                   &  tend%    q(:,:,jb,iqv),         &! inout
                   &  tend%    q(:,:,jb,iqc),         &! inout
                   &  tend%    q(:,:,jb,iqi),         &! inout
                   &  tend%    q(:,:,jb,iqt:),        &! inout
                   &  tend%    u_vdf(:,:,jb),         &! out
                   &  tend%    v_vdf(:,:,jb),         &! out
                   &  tend% temp_vdf(:,:,jb),         &! out
                   &  tend%    q_vdf(:,:,jb,iqv),     &! out
                   &  tend%    q_vdf(:,:,jb,iqc),     &! out
                   &  tend%    q_vdf(:,:,jb,iqi),     &! out
                   &  tend%    q_vdf(:,:,jb,iqt:),    &! out
                   &  zqtvar_prod,                    &! out, for "cloud" ("zvdiffp" in echam)
                   &  zvmixtau,                       &! out, for "cloud"
                   &  zqhflx,                         &! out, for "cloud"
                   & field% z0m   (:,  jb),           &! out, for the next step
                   & field% thvvar(:,:,jb),           &! out, for the next step
                   & field% thvsig(:,  jb),           &! out, for "cucall"
                   & field%    tke(:,:,jb)            )! out

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

      IF (ltimer) CALL timer_stop(timer_vdiff)
    ELSE
      zvmixtau   (jcs:jce,:) = 0._wp
      zqhflx     (jcs:jce)   = 0._wp
      zqtvar_prod(jcs:jce,:) = 0._wp

      tend% u_vdf(jcs:jce,:,jb) = 0._wp
      tend% v_vdf(jcs:jce,:,jb) = 0._wp
    ENDIF !lvdiff

    !-------------------------------------------------------------------
    ! 6. ATMOSPHERIC GRAVITY WAVES
    !-------------------------------------------------------------------

    ! 6.1   CALL SUBROUTINE GW_HINES

    IF (echam_phy_config%lgw_hines) THEN

      IF (ltimer) call timer_start(timer_gw_hines)

      zlat_deg(jcs:jce) = p_patch(jg)%cells%center(jcs:jce,jb)%lat * 180._wp/pi

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
        &             tend% temp_gwh(:,:,jb)   ,&
        &             tend%    u_gwh(:,:,jb)   ,&
        &             tend%    v_gwh(:,:,jb) )

      IF (ltimer) call timer_stop(timer_gw_hines)

      tend% temp(jcs:jce,:,jb) = tend% temp(jcs:jce,:,jb) + tend% temp_gwh(jcs:jce,:,jb)
      tend%    u(jcs:jce,:,jb) = tend%    u(jcs:jce,:,jb) + tend%    u_gwh(jcs:jce,:,jb)
      tend%    v(jcs:jce,:,jb) = tend%    v(jcs:jce,:,jb) + tend%    v_gwh(jcs:jce,:,jb)

    ELSE ! NECESSARY COMPUTATIONS IF GW_HINES IS BY-PASSED

      tend% temp_gwh(jcs:jce,:,jb) = 0._wp
      tend%    u_gwh(jcs:jce,:,jb) = 0._wp
      tend%    v_gwh(jcs:jce,:,jb) = 0._wp

    END IF !lgw_hines

    !-------------------------------------------------------------------
    ! 7. CONVECTION PARAMETERISATION
    !-------------------------------------------------------------------
    itype(jcs:jce) = 0

    ! 7.1   INITIALIZE ARRAYS FOR CONVECTIVE PRECIPITATION
    !       AND COPY ARRAYS FOR CONVECTIVE CLOUD PARAMETERS

    tend% x_dtr(jcs:jce,:,jb) = 0._wp
    zqtec  (jcs:jce,:) = 0._wp
    zxtecl (jcs:jce,:) = 0._wp
    zxteci (jcs:jce,:) = 0._wp
    zxtecnl(jcs:jce,:) = 0._wp
    zxtecni(jcs:jce,:) = 0._wp

    field% rsfc(jcs:jce,jb) = 0._wp
    field% ssfc(jcs:jce,jb) = 0._wp

    ! 7.2   CALL SUBROUTINE CUCALL FOR CUMULUS PARAMETERIZATION

    IF (echam_phy_config%lconv) THEN

      IF (ltimer) call timer_start(timer_cucall)

      CALL cucall( echam_conv_config%ncvmicro,&! in
        &          echam_conv_config%iconv,   &! in
        &          echam_conv_config%lmfdudv, &! in
        &          echam_conv_config%lmfdd,   &! in
        &          echam_conv_config%lmfmid,  &! in
        &          echam_conv_config%dlev,    &! in
        &          echam_conv_config%cmftau,  &! in
        &          echam_conv_config%cmfctop, &! in
        &          echam_conv_config%cprcon,  &! in
        &          echam_conv_config%cminbuoy,&! in
        &          echam_conv_config%entrpen, &! in
        &          echam_conv_config%nmctop,  &! in
        &          echam_conv_config%cevapcu, &! in
        &          jce, nbdim, nlev,          &! in
        &          nlevp1, nlevm1,            &! in
        &          ntrac,                     &! in     tracers
!0      &          jb,                        &! in     row index
        &          pdtime, psteplen,          &! in
        &          field% lfland(:,jb),       &! in     loland
        &          field% temp(:,:,jb),       &! in     tm1
        &          field% u(:,:,jb),          &! in     um1
        &          field% v(:,:,jb),          &! in     vm1
        &          field% omega(:,:,jb),      &! in     vervel
        &          field% geom(:,:,jb),       &! in     geom1
        &          field% q(:,:,jb,iqv),      &! in     qm1
        &          field% q(:,:,jb,iqc),      &! in     xlm1
        &          field% q(:,:,jb,iqi),      &! in     xim1
        &          field% q(:,:,jb,iqt:),     &! in     xtm1
        &           tend% q(:,:,jb,iqc),      &! in     xlm1
        &           tend% q(:,:,jb,iqi),      &! in     xim1
        &          field% presm_new(:,:,jb),  &! in     app1
        &          field% presi_new(:,:,jb),  &! in     aphp1
        &          zqhflx,                    &! in     qhfla (from "vdiff")
!0      &          field% tke,                &! in     tkem1 (from "vdiff")
        &          field% thvsig(:,jb),       &! in           (from "vdiff")
        &          tend% temp(:,:,jb),        &! inout  tte
        &          tend% u(:,:,jb),           &! inout  vom
        &          tend% v(:,:,jb),           &! inout  vol
        &          tend% q(:,:,jb,iqv),       &! inout  qte
        &          tend% q(:,:,jb,iqt:),      &! inout  xtte
        &          zqtec,                     &! inout
        &          tend% x_dtr(:,:,jb),       &! inout  xtec
        &          zxtecl,  zxteci,           &! inout
        &          zxtecnl, zxtecni,          &! inout
        &          field% rsfc(:,jb),         &! inout
        &          field% ssfc(:,jb),         &! inout
        &          field% aprc(:,jb),         &! inout
        &          field% aprs(:,jb),         &! inout
        &          field% topmax(:,jb),       &! inout
        &          itype,                     &! inout
        &          ilab,                      &! out
        &          zcvcbot,                   &! out
        &          zwcape,                    &! out
        &          tend%temp_cnv(:,:,jb),     &! out
        &          tend%   u_cnv(:,:,jb),     &! out
        &          tend%   v_cnv(:,:,jb),     &! out
        &          tend%   q_cnv(:,:,jb,iqv), &! out
        &          tend%   q_cnv(:,:,jb,iqt:) )! out

      field% rtype(jcs:jce,jb) = REAL(itype(jcs:jce),wp)

      IF (ltimer) call timer_stop(timer_cucall)

    ELSE ! NECESSARY COMPUTATIONS IF MASSFLUX IS BY-PASSED

      ilab(jcs:jce,1:nlev) = 0

      tend% u_cnv(jcs:jce,:,jb) = 0._wp
      tend% v_cnv(jcs:jce,:,jb) = 0._wp

    ENDIF !lconv

    !-------------------------------------------------------------
    ! 7. LARGE SCALE CONDENSATION.
    !-------------------------------------------------------------
    IF(echam_phy_config%lcond) THEN

      !IF (lcotra) CALL get_col_pol( tend%temp(:,:,jb),tend%q(:,:,jb,iqv),jb )
      IF (ltimer) CALL timer_start(timer_cloud)

      IF (ncdnc==0 .AND. nicnc==0) THEN

        CALL cloud(jce, nbdim, jks, nlev, nlevp1, ntrac,  &! in
!0        &        jb,                                    &! in
          &        pdtime, psteplen, echam_phy_config%lcover, &! in
          &        field% presi_old(:,:,jb), &! in
!0        &        field% presi_new(:,:,jb), &! in
          &        field% presm_old(:,:,jb), &! in
          &        field% presm_new(:,:,jb), &! in
          &        field% temp (:,:,jb),     &! in. tm1
          &        field%   tv (:,:,jb),     &! in. ztvm1
          &        field% geom (:,:,jb),     &! in. geom1
          &        field% omega(:,:,jb),     &! in. vervel
          &        field% acdnc(:,:,jb),     &! in. acdnc
          &        field% q(:,:,jb,iqv),     &! in.  qm1
          &        field% q(:,:,jb,iqc),     &! in. xlm1
          &        field% q(:,:,jb,iqi),     &! in. xim1
!0        &        field% q(:,:,jb,iqt:),    &! in. xtm1
          &        zqtvar_prod, zhmixtau, zvmixtau, &! in
          &        zbetaa, zbetab, zbetass, invb,   &! in (from "cover")
          &        zqtec,                     &! inout (there is a clip inside)
          &         tend% x_dtr(:,:,jb),      &! inout (there is a clip inside)
          &         tend% temp(:,:,jb),       &! inout.  tte
          &         tend% q(:,:,jb,iqv),      &! inout.  qte
          &         tend% q(:,:,jb,iqc),      &! inout. xlte
          &         tend% q(:,:,jb,iqi),      &! inout. xite
          &         tend% q(:,:,jb,iqt:),     &! inout. xtte
          &        field% xvar  (:,:,jb),     &! inout
          &        field% xskew (:,:,jb),     &! inout
          &        field% aclc  (:,:,jb),     &! inout
          &        field% aclcac(:,:,jb),     &! inout
          &        field% aclcov(:,  jb),     &! inout
          &        field%  qvi  (:,  jb),     &! inout
          &        field% xlvi  (:,  jb),     &! inout
          &        field% xivi  (:,  jb),     &! inout
          &        field% aprl  (:,  jb),     &! inout
          &        field% aprs  (:,  jb),     &! inout
          &        field% rsfl  (:,  jb),     &! out
          &        field% ssfl  (:,  jb),     &! out
          &        field% relhum(:,:,jb),     &! out
          &        tend%temp_cld(:,:,jb),     &! out
          &        tend%   q_cld(:,:,jb,iqv), &! out
          &        tend%   q_cld(:,:,jb,iqc), &! out
          &        tend%   q_cld(:,:,jb,iqi), &! out
          &        tend%   q_cld(:,:,jb,iqt:) )! out

      ELSE IF (ncdnc>0 .AND. nicnc>0) THEN
!0      CALL cloud_cdnc_icnc(...) !!skipped in ICON
      ELSE
        IF (my_process_is_stdio()) CALL finish('physc', ' check setting of ncdnc and nicnc.')
      END IF

      IF (ltimer) CALL timer_stop(timer_cloud)

    ELSE ! NECESSARY COMPUTATIONS IF *CLOUD* IS BY-PASSED.

      field% rsfl (jcs:jce,  jb) = 0._wp
      field% ssfl (jcs:jce,  jb) = 0._wp
      field% aclc (jcs:jce,:,jb) = 0._wp

    ENDIF !lcond

    ! Done. Disassociate pointers.
    NULLIFY(field,tend)

  END SUBROUTINE physc
  !-------------

END MODULE mo_echam_phy_main
