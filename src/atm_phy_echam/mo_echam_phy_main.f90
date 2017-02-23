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

!----------------------------
#include "omp_definitions.inc"
!----------------------------
#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
@PROCESS SPILLSIZE(5000)
#endif
!OCL NOALIAS

MODULE mo_echam_phy_main

  USE mo_kind,                ONLY: wp
  USE mo_physical_constants,  ONLY: cpd, cpv, cvd, cvv, Tf, tmelt
  USE mo_impl_constants      ,ONLY: inh_atmosphere
  USE mo_run_config,          ONLY: ntracer, nlev, nlevp1,    &
    &                               iqv, iqi, iqt
  USE mo_echam_phy_config,    ONLY: echam_phy_config
  USE mo_echam_cloud_config,  ONLY: echam_cloud_config
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field, prm_field,     &
    &                               t_echam_phy_tend,  prm_tend,      &
    &                               cdimissval
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop,                &
    &                               timer_cover
  USE mtime,                  ONLY: datetime
  USE mo_echam_sfc_indices,   ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_cover,               ONLY: cover

  USE mo_interface_echam_radiation,    ONLY: interface_echam_radiation
  USE mo_interface_echam_radheating,   ONLY: interface_echam_radheating
  USE mo_interface_echam_vdiff_surface,ONLY: interface_echam_vdiff_surface
  USE mo_interface_echam_o3_cariolle,  ONLY: interface_echam_o3_cariolle
  USE mo_interface_echam_gravity_waves,ONLY: interface_echam_gravity_waves
  USE mo_interface_echam_convection,   ONLY: interface_echam_convection
  USE mo_interface_echam_condensation, ONLY: interface_echam_condensation

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_model_domain        ,ONLY: t_patch
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: echam_phy_main

  ! variables that are time-independent, potentailly  shared by the routines of this module 
  ! they could have been calculated at initializatiop phase 
  !    these can be shared through physcics parameters
  INTEGER  :: ntrac !< # of tracers excluding water vapour and hydrometeors
                      !< (handled by sub-models, e.g., chemical species)

  INTEGER  :: jks   !< start index for vertical loops


CONTAINS

  !>
  !!
  SUBROUTINE echam_phy_main( patch,                         &
    &                        rl_start, rl_end,              &
    &                        this_datetime,pdtime,        &
    &                        ltrig_rad                      )

    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(IN)  :: rl_start, rl_end

    TYPE(datetime),  POINTER     :: this_datetime    !< date and time 
    REAL(wp)        ,INTENT(IN)  :: pdtime         !< time step

    LOGICAL         ,INTENT(IN)  :: ltrig_rad      !< perform radiative transfer computation

    ! Local variables
    TYPE(t_echam_phy_field),   POINTER :: field
    TYPE(t_echam_phy_tend) ,   POINTER :: tend

    REAL(wp) :: zconv  (nproma,nlev,patch%nblks_c)       !< conversion factor q-->dT/dt       [(K/s)/(W/m2)]
    REAL(wp) :: zcpair (nproma,nlev,patch%nblks_c)       !< specific heat of moist air at const. pressure [J/K/kg]
    REAL(wp) :: zcvair (nproma,nlev,patch%nblks_c)       !< specific heat of moist air at const. volume   [J/K/kg]

    REAL(wp) :: zq_phy (nproma,nlev,patch%nblks_c)       !< heating by whole ECHAM physics    [W/m2]
    REAL(wp) :: zq_rlw (nproma,nlev,patch%nblks_c)       !< heating by long  wave radiation   [W/m2]
    REAL(wp) :: zq_rlw_impl (nproma)       !< additional heating by LW rad. due to impl. coupling in surface energy balance [W/m2]

    INTEGER  :: ictop (nproma,patch%nblks_c)             !< from massflux

    INTEGER  :: jg             !< grid level/domain index
    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block

    jg         = patch%id
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    ! start index for vertical loops
    jks=1

    ntrac = ntracer-iqt+1  !# of tracers excluding water vapour and hydrometeors

    ! 1. Associate pointers
    field  => prm_field(jg)
    tend   => prm_tend (jg)

    ! provisionally copy the incoming tedencies
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      tend%   ta_phy (jcs:jce,:,jb)   = tend%   ta (jcs:jce,:,jb)
      tend%   ua_phy (jcs:jce,:,jb)   = tend%   ua (jcs:jce,:,jb)
      tend%   va_phy (jcs:jce,:,jb)   = tend%   va (jcs:jce,:,jb)
      tend% qtrc_phy (jcs:jce,:,jb,:) = tend% qtrc (jcs:jce,:,jb,:)
    ENDDO
!$OMP END PARALLEL DO 

    ! initialize physics accumulated heating
    zq_phy(:,:,:) = 0._wp

    !------------------------------------------------------------
    ! 3. COMPUTE SOME FIELDS NEEDED BY THE PHYSICAL ROUTINES.
    !------------------------------------------------------------

    ! 3.2b Specific heat of moist air
    !
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      CALL calculate_zcpair_zconv( jb,jcs,jce, nproma, field, zcpair(:,:,jb), zcvair(:,:,jb), zconv(:,:,jb))
    ENDDO
!$OMP END PARALLEL DO 

    !-------------------------------------------------------------------
    ! 3.13 DIAGNOSE CURRENT CLOUD COVER
    !
    !  and land fraction (combined here for convenience)
    !-------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      CALL landFraction_cloudCover( jb,jcs,jce, nproma, field)
    ENDDO
!$OMP END PARALLEL DO 

    !---------------------------------------------------------------------
    ! 3.9 DETERMINE TROPOPAUSE HEIGHT AND MASS BUDGETS
    !     (Needed only for sub-models. Note: sequence of arguments
    !      different from the original ECHAM6)
    !---------------------------------------------------------------------
    !
    !CALL WMO_tropopause( jce, nproma, nlev,         &! in
    !                   & ncctop, nccbot, lresum,   &! in
    !                   & field% ta(:,:,jb),        &! in
    !                   & field% presm_old(:,:,jb), &! in
    !                   & field% tropo(:,jb),       &! out for diagnostics
    !                   & itrpwmo, itrpwmop1        )! out for submodel
    !---------------------------------------------------------------------


    !---------------------------------------------------------------------
    ! 3.12 INITIALISATION OF CLOUD DROPLET NUMBER CONCENTRATION 
    !      (1/M**3) USED IN RADLSW AND CLOUD
    !---------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      CALL calculate_droplet_number( jb,jcs,jce, nproma, field)
    ENDDO
!$OMP END PARALLEL DO 


    !-------------------------------------------------------------------
    ! 4. RADIATION PARAMETERISATION
    !-------------------------------------------------------------------
    IF (echam_phy_config%lrad .AND. ltrig_rad) THEN

      ! 4.1 RADIATIVE TRANSFER
      !----------------------- 
      ! store ts_rad of this radiatiative transfer timestep in ts_rad_rt,
      ! so that it can be reused in radheat in the other timesteps
      field%ts_rad_rt(:,:) = field%ts_rad(:,:)

      CALL interface_echam_radiation( patch, field, this_datetime)

    END IF 

    !-------------------------------------------------------------------
    ! 4.2 RADIATIVE HEATING
    !----------------------
    CALL  interface_echam_radheating(patch, rl_start, rl_end, field, tend, zconv, zq_rlw, zq_phy)


    !-------------------------------------------------------------------
    ! 5. BOUNDARY LAYER AND SURFACE PROCESSES
    !-------------------------------------------------------------------
    ! Note: In ECHAM this part is located between "CALL radiation" and
    !       "CALL radheat".
    CALL interface_echam_vdiff_surface(patch, rl_start, rl_end, field, tend, zconv, zq_phy, pdtime)
    !-----------------------


    !-----------------------
    ! This can probably be moved to the radheating 
    IF (echam_phy_config%lrad) THEN

!$OMP PARALLEL DO PRIVATE(jcs,jce, zq_rlw_impl)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
        ! Heating due to the fact that surface model only used part of longwave radiation to compute new surface temperature
        zq_rlw_impl(jcs:jce) =                                                  &
          &   ( (field%rld_rt(jcs:jce,nlev,jb)-field%rlu_rt(jcs:jce,nlev,jb))   & ! (( rlns from "radiation"
          &    -(field%rlds  (jcs:jce,jb)     -field%rlus  (jcs:jce,jb)     ))  & !   -rlns from "radheating" and "update_surface")
          &  -zq_rlw(jcs:jce,nlev,jb)                                                ! old heating from radheat

        ! Heating accumulated
        zq_phy(jcs:jce,nlev,jb) = zq_phy(jcs:jce,nlev,jb) + zq_rlw_impl(jcs:jce)

        ! Tendency
        tend%ta_rlw_impl(jcs:jce,jb) = zq_rlw_impl(jcs:jce) * zconv(jcs:jce,nlev,jb)

        ! Tendencies accumulated
        tend%ta(jcs:jce,nlev,jb) = tend%ta(jcs:jce,nlev,jb) + tend%ta_rlw_impl(jcs:jce,jb)
      ENDDO
!$OMP END PARALLEL DO 

    END IF
    !---------------------

    !---------------------
    IF (echam_phy_config%lcariolle) THEN
      CALL  interface_echam_o3_cariolle(patch, rl_start, rl_end, field, tend, this_datetime) 
    END IF
    !---------------------

    !-------------------------------------------------------------------
    ! 6. ATMOSPHERIC GRAVITY WAVES
    !-------------------------------------------------------------------
    CALL  interface_echam_gravity_waves(patch, rl_start, rl_end, field, tend, zconv, zq_phy, pdtime)

    !-------------------------------------------------------------------
    ! 7. CONVECTION PARAMETERISATION
    !-------------------------------------------------------------------
    CALL interface_echam_convection(patch, rl_start, rl_end, field, tend, &
      & ictop, zconv, zq_phy, ntrac, pdtime)
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    CALL interface_echam_condensation(patch, rl_start, rl_end, jks, field, tend, &
      & zcpair, zconv, zq_phy, ictop, pdtime)
    !-------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

    ! KF accumulate fields for diagnostics

    !  total precipitation flux
      field% totprec (jcs:jce,jb)     =  field% rsfl (jcs:jce,jb) & ! rain large scale
            &                            +field% ssfl (jcs:jce,jb) & ! snow large scale
            &                            +field% rsfc (jcs:jce,jb) & ! rain convection
            &                            +field% ssfc (jcs:jce,jb)   ! snow convection

      ! Now compute tendencies from physics alone

      tend%   ta_phy (jcs:jce,:,jb)   = tend%   ta (jcs:jce,:,jb)   - tend%   ta_phy (jcs:jce,:,jb)
      tend%   ua_phy (jcs:jce,:,jb)   = tend%   ua (jcs:jce,:,jb)   - tend%   ua_phy (jcs:jce,:,jb)
      tend%   va_phy (jcs:jce,:,jb)   = tend%   va (jcs:jce,:,jb)   - tend%   va_phy (jcs:jce,:,jb)
      tend% qtrc_phy (jcs:jce,:,jb,:) = tend% qtrc (jcs:jce,:,jb,:) - tend% qtrc_phy (jcs:jce,:,jb,:)

      tend% ta_phy (jcs:jce,:,jb) = tend% ta_phy(jcs:jce,:,jb)*zcpair(jcs:jce,:,jb)/zcvair(jcs:jce,:,jb)
    ENDDO
!$OMP END PARALLEL DO 

    ! Done. Disassociate pointers.
    NULLIFY(field,tend)

  END SUBROUTINE echam_phy_main
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  SUBROUTINE landFraction_cloudCover( jb,jcs,jce, nbdim, field)

    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block
    TYPE(t_echam_phy_field),   POINTER :: field

    REAL(wp) :: zfrl (nbdim)              !< fraction of land in the grid box
    REAL(wp) :: zfrw (nbdim)              !< fraction of water (without ice) in the grid point
    REAL(wp) :: zfri (nbdim)              !< fraction of ice in the grid box
    INTEGER  :: itype(nbdim)              !< type of convection

    INTEGER  :: jc
 
 
     ! 3.3 Weighting factors for fractional surface coverage
    !     Accumulate ice portion for diagnostics

    DO jc=jcs,jce

      ! fraction of land in the grid box.
      ! lsmask: land-sea mask, depends on input data, either:
      ! fractional, including or excluding lakes in the land part or
      ! non-fractional, each grid cell is either land, sea, or sea-ice.
      ! See mo_echam_phy_init or input data set for details.

      zfrl(jc) = field% lsmask(jc,jb)

      ! fraction of sea/lake in the grid box
      ! * (1. - fraction of sea ice in the sea/lake part of the grid box)
      ! => fraction of open water in the grid box

      zfrw(jc) = (1._wp-zfrl(jc))*(1._wp-field%seaice(jc,jb))

      ! fraction of sea ice in the grid box
      zfri(jc) = 1._wp-zfrl(jc)-zfrw(jc)
      ! security for ice temperature with changing ice mask
      !
      IF(zfri(jc) > 0._wp .AND. field%ts_tile(jc,jb,iice) == cdimissval ) THEN
         field% ts_tile(jc,jb,iice)  = tmelt + Tf    ! = 271.35 K
      ENDIF
    END DO

    ! 3.4 Merge three pieces of information into one array for vdiff
    IF (ilnd.LE.nsfc_type) field%frac_tile(jcs:jce,jb,ilnd) = zfrl(jcs:jce)
    IF (iwtr.LE.nsfc_type) field%frac_tile(jcs:jce,jb,iwtr) = zfrw(jcs:jce)
    IF (iice.LE.nsfc_type) field%frac_tile(jcs:jce,jb,iice) = zfri(jcs:jce)


    !-------------------------------------------------------------------
    ! 3.13 DIAGNOSE CURRENT CLOUD COVER
    !-------------------------------------------------------------------

    itype(jcs:jce) = NINT(field%rtype(jcs:jce,jb))
    IF (echam_phy_config%lcond) THEN
      IF (ltimer) CALL timer_start(timer_cover)

      CALL cover( jce, nbdim, jks,          &! in
        &         nlev, nlevp1,             &! in
        &         itype,  zfrw, zfri,       &! in
        &         field% zf(:,:,jb),        &! in
        &         field% presi_old(:,:,jb), &! in
        &         field% presm_old(:,:,jb), &! in
        &         field%  ta(:,:,jb),       &! in    tm1
        &         field%  qtrc(:,:,jb,iqv), &! in    qm1
        &         field%  qtrc(:,:,jb,iqi), &! in    xim1
        &         field%  aclc(:,:,jb),     &! out   (for "radiation" and "vdiff_down")
        &         field% rintop(:,  jb)    ) ! out   (for output)

      IF (ltimer) CALL timer_stop(timer_cover)
    ENDIF ! lcond

  END SUBROUTINE landFraction_cloudCover
  !----------------------------------------------------------------

  !---------------------------------------------------------------------
  ! 3.12 INITIALISATION OF CLOUD DROPLET NUMBER CONCENTRATION 
  !      (1/M**3) USED IN RADLSW AND CLOUD
  !---------------------------------------------------------------------
  SUBROUTINE calculate_droplet_number( jb,jcs,jce, nbdim, field)
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block
    TYPE(t_echam_phy_field),   POINTER :: field

    INTEGER  :: jc, jk
    LOGICAL  :: lland(nbdim), lglac(nbdim)
    REAL(wp) :: zprat, zn1, zn2, zcdnc


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

  END SUBROUTINE calculate_droplet_number
  !----------------------------------------------------------------

 
  !---------------------------------------------------------------------
  SUBROUTINE calculate_zcpair_zconv( jb,jcs,jce, nbdim, field, zcpair, zcvair, zconv)
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block
    TYPE(t_echam_phy_field),   POINTER :: field    ! in
    REAL(wp)         ,INTENT(INOUT) :: zcpair  (nbdim,nlev)       !< specific heat of moist air        [J/K/kg]
    REAL(wp)         ,INTENT(INOUT) :: zcvair  (nbdim,nlev)     
    REAL(wp)         ,INTENT(INOUT) :: zconv   (nbdim,nlev)      

    INTEGER  :: jc, jk

    zcpair(:,:) = 0.0_wp
    zcvair(:,:) = 0.0_wp
    zconv(:,:) = 0.0_wp

    DO jk = 1,nlev
      DO jc = jcs,jce
        !
        ! 3.2b Specific heat of moist air
        !
        zcpair  (jc,jk) = cpd+(cpv-cpd)*field%qtrc(jc,jk,jb,iqv)
        zcvair  (jc,jk) = cvd+(cvv-cvd)*field%qtrc(jc,jk,jb,iqv)
        !
        zconv   (jc,jk) = 1._wp/(field%mair(jc,jk,jb)*zcpair(jc,jk))
     END DO
    END DO
  END SUBROUTINE calculate_zcpair_zconv
  !---------------------------------------------------------------------



END MODULE mo_echam_phy_main
