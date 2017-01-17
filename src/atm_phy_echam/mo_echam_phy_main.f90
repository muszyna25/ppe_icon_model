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
  USE mo_exception,           ONLY: finish
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_math_constants,      ONLY: pi
  USE mo_physical_constants,  ONLY: cpd, cpv, cvd, cvv
  USE mo_impl_constants      ,ONLY: inh_atmosphere
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
    &                               timer_cover, timer_radheat,                     &
    &                               timer_gw_hines, timer_ssodrag,                  &
    &                               timer_cucall, timer_cloud
  USE mtime,                  ONLY: datetime
  USE mo_ham_aerosol_params,  ONLY: ncdnc, nicnc
  USE mo_echam_sfc_indices,   ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_cloud,               ONLY: cloud
  USE mo_cover,               ONLY: cover
  USE mo_gw_hines,            ONLY: gw_hines
  USE mo_ssortns,             ONLY: ssodrag

  USE mo_interface_echam_radiation,    ONLY: interface_echam_radiation
  USE mo_interface_echam_radheating,   ONLY: interface_echam_radheating
  USE mo_interface_echam_vdiff_surface,ONLY: interface_echam_vdiff_surface
  USE mo_interface_echam_o3_cariolle,  ONLY: interface_echam_o3_cariolle

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
  REAL(wp) :: zcd                       !< specific heat of dry air          [J/K/kg]
  REAL(wp) :: zcv                       !< specific heat of water vapor      [J/K/kg]

  INTEGER  :: jks   !< start index for vertical loops

  ! the following would depend on the nesting, ie jg
  REAL(wp) :: pdtime         !< time step

CONTAINS

  !>
  !!
  SUBROUTINE echam_phy_main( patch,                         &
    &                        rl_start, rl_end,              &
    &                        this_datetime,in_pdtime,        &
    &                        ltrig_rad                      )

    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(IN)  :: rl_start, rl_end

    TYPE(datetime),  POINTER     :: this_datetime    !< date and time 
    REAL(wp)        ,INTENT(IN)  :: in_pdtime         !< time step

    LOGICAL         ,INTENT(IN)  :: ltrig_rad      !< perform radiative transfer computation

    ! Local variables
    TYPE(t_echam_phy_field),   POINTER :: field
    TYPE(t_echam_phy_tend) ,   POINTER :: tend

    REAL(wp) :: zfrc (nproma,nsfc_type,patch%nblks_c)    !< zfrl, zfrw, zfrc combined

    REAL(wp) :: zcair  (nproma,nlev,patch%nblks_c)       !< specific heat of moist air        [J/K/kg]
    REAL(wp) :: zconv  (nproma,nlev,patch%nblks_c)       !< conversion factor q-->dT/dt       [(K/s)/(W/m2)]
    REAL(wp) :: zcpair (nproma,nlev)       !< specific heat of moist air at const. pressure [J/K/kg]
    REAL(wp) :: zcvair (nproma,nlev)       !< specific heat of moist air at const. volume   [J/K/kg]

    REAL(wp) :: zq_phy (nproma,nlev,patch%nblks_c)       !< heating by whole ECHAM physics    [W/m2]
    REAL(wp) :: zq_rlw (nproma,nlev,patch%nblks_c)       !< heating by long  wave radiation   [W/m2]
    REAL(wp) :: zq_rlw_impl (nproma)       !< additional heating by LW rad. due to impl. coupling in surface energy balance [W/m2]

    INTEGER  :: jg             !< grid level/domain index
    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block

    jg         = patch%id
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    ! fill the module timestep lengths
    pdtime = in_pdtime
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

    ! Use constant volume or constant pressure specific heats for dry air and vapor
    ! for the computation of temperature tendencies.
!!$    IF ( iequations == inh_atmosphere ) THEN
!!$      zcd = cvd
!!$      zcv = cvv
!!$    ELSE
      zcd = cpd
      zcv = cpv
!!$    END IF

    !-------------------------------------------------------------------
    ! 3.13 DIAGNOSE CURRENT CLOUD COVER
    !
    !  and land fraction zfrc used in vdiff. (combined here for convinience)
    !-------------------------------------------------------------------
 
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      CALL landFraction_cloudCover( jb,jcs,jce, nproma, field, zfrc(:,:,jb))
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
    IF (phy_config%lrad .AND. ltrig_rad) THEN

      ! 4.1 RADIATIVE TRANSFER
      !----------------------- 
      ! store ts_rad of this radiatiative transfer timestep in ts_rad_rt,
      ! so that it can be reused in radheat in the other timesteps
      field%ts_rad_rt(:,:) = field%ts_rad(:,:)

      CALL interface_echam_radiation( patch, field, this_datetime)

    END IF 

!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      CALL calculate_zcair_zconv( jb,jcs,jce, nproma, field, zcair(:,:,jb), zconv(:,:,jb))
    ENDDO
!$OMP END PARALLEL DO 

    ! initialize physics accumulated heating
    zq_phy(:,:,:) = 0._wp

    !-------------------------------------------------------------------
    ! 4.2 RADIATIVE HEATING
    !----------------------
    CALL  interface_echam_radheating(patch, rl_start, rl_end, field, tend, zconv, zq_rlw, zq_phy)


    !-------------------------------------------------------------------
    ! 5. BOUNDARY LAYER AND SURFACE PROCESSES
    !-------------------------------------------------------------------
    ! Note: In ECHAM this part is located between "CALL radiation" and
    !       "CALL radheat".
    CALL interface_echam_vdiff_surface(patch, rl_start, rl_end, field, tend, zfrc, in_pdtime, zcd, zcd)
    !-----------------------


    !-----------------------
    ! This can probably be moved to the radheating 
    IF (phy_config%lrad) THEN

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

    ELSE

      tend%ta_rlw_impl(:,:) = 0._wp

    END IF
    !---------------------

    !---------------------
    IF (phy_config%lcariolle) THEN
      CALL  interface_echam_o3_cariolle(patch, rl_start, rl_end, field, tend, this_datetime) 
    END IF
    !---------------------

    !-------------------------------------------------------------------
    ! 6. ATMOSPHERIC GRAVITY WAVES
    !-------------------------------------------------------------------

    ! 6.1   CALL SUBROUTINE GW_HINES

    IF (phy_config%lgw_hines) THEN

      IF (ltimer) call timer_start(timer_gw_hines)
!$OMP PARALLEL DO PRIVATE(jcs,jce)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

        CALL echam_gw_hinesg(patch, jg, jb,jcs,jce, nproma, field, tend, zconv(:,:,jb), zq_phy(:,:,jb))
      ENDDO
!$OMP END PARALLEL DO 

      IF (ltimer) call timer_stop(timer_gw_hines)
      
    ELSE ! NECESSARY COMPUTATIONS IF GW_HINES IS BY-PASSED
      ! this should not be necessary
!$OMP PARALLEL DO PRIVATE(jcs,jce)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
        ! this should not be necessary, if are initialize to zero
        tend%   ta_gwh(jcs:jce,:,jb) = 0._wp
        tend%   ua_gwh(jcs:jce,:,jb) = 0._wp
        tend%   va_gwh(jcs:jce,:,jb) = 0._wp
      ENDDO
!$OMP END PARALLEL DO 

    END IF !lgw_hines


    ! 6.2   CALL SUBROUTINE SSODRAG

    IF (phy_config%lssodrag) THEN

      IF (ltimer) call timer_start(timer_ssodrag)

!$OMP PARALLEL DO PRIVATE(jcs,jce)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

        CALL echam_ssodrag(patch, jb,jcs,jce, nproma, field, tend, zconv(:,:,jb), zq_phy(:,:,jb))
      ENDDO
!$OMP END PARALLEL DO 

      IF (ltimer) call timer_stop(timer_ssodrag)

    ELSE ! NECESSARY COMPUTATIONS IF SSODRAG IS BY-PASSED

      ! not necessary, if initialized with zeroes
!$OMP PARALLEL DO PRIVATE(jcs,jce)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
        tend%   ta_sso(jcs:jce,:,jb) = 0._wp
        tend%   ua_sso(jcs:jce,:,jb) = 0._wp
        tend%   va_sso(jcs:jce,:,jb) = 0._wp
      ENDDO
!$OMP END PARALLEL DO 

    END IF ! SSODRAG

    !-------------------------------------------------------------------
    ! 7. CONVECTION PARAMETERISATION
    !-------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

      CALL echam_cumulus_condensation(jb,jcs,jce, nproma, field, tend, zcair(:,:,jb))

    ENDDO
!$OMP END PARALLEL DO 
    !-------------------------------------------------------------------


!$OMP PARALLEL DO PRIVATE(jcs,jce, zcpair, zcvair)
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

      IF ( iequations == inh_atmosphere ) THEN
        !
        zcpair  (jcs:jce,:) = cpd+(cpv-cpd)*field%qtrc(jcs:jce,:,jb,iqv)
        zcvair  (jcs:jce,:) = cvd+(cvv-cvd)*field%qtrc(jcs:jce,:,jb,iqv)
        !
        tend% ta_phy (jcs:jce,:,jb) = tend% ta_phy(jcs:jce,:,jb)*zcpair(jcs:jce,:)/zcvair(jcs:jce,:)
      END IF
    ENDDO
!$OMP END PARALLEL DO 

    ! Done. Disassociate pointers.
    NULLIFY(field,tend)

  END SUBROUTINE echam_phy_main
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  SUBROUTINE landFraction_cloudCover( jb,jcs,jce, nbdim, field, zfrc)

    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block
    TYPE(t_echam_phy_field),   POINTER :: field

    REAL(wp),        INTENT(inout) :: zfrc (nbdim,nsfc_type)    !< zfrl, zfrw, zfrc combined, output

    REAL(wp) :: zfrl (nbdim)              !< fraction of land in the grid box
    REAL(wp) :: zfrw (nbdim)              !< fraction of water (without ice) in the grid point
    REAL(wp) :: zfri (nbdim)              !< fraction of ice in the grid box
    INTEGER  :: itype(nbdim)              !< type of convection

    INTEGER  :: jc
 
 
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
  SUBROUTINE calculate_zcair_zconv( jb,jcs,jce, nbdim, field, zcair, zconv)
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block
    TYPE(t_echam_phy_field),   POINTER :: field    ! in
    REAL(wp)         ,INTENT(INOUT) :: zcair  (nbdim,nlev)       !< specific heat of moist air        [J/K/kg]
    REAL(wp)         ,INTENT(INOUT) :: zconv  (nbdim,nlev)       !< specific heat of moist air        [J/K/kg]

    INTEGER  :: jc, jk

    DO jk = 1,nlev
      DO jc = jcs,jce
        !
        ! 3.2b Specific heat of moist air
        !
        zcair   (jc,jk) = zcd+(zcv-zcd)*field%qtrc(jc,jk,jb,iqv)
        zconv   (jc,jk) = 1._wp/(field%mair(jc,jk,jb)*zcair(jc,jk))
      END DO
    END DO
  END SUBROUTINE calculate_zcair_zconv
  !---------------------------------------------------------------------



  !---------------------------------------------------------------------
  SUBROUTINE echam_gw_hinesg(patch, jg, jb,jcs,jce, nbdim, field, tend, zconv, zq_phy)
    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(IN) :: jg
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block  
    REAL(wp)        ,INTENT(IN) :: zconv  (nbdim,nlev)       !< specific heat of moist air        [J/K/kg]
    TYPE(t_echam_phy_field),   POINTER :: field
    TYPE(t_echam_phy_tend) ,   POINTER :: tend
    REAL(wp)        ,INTENT(INOUT) :: zq_phy (nbdim,nlev)       !< heating by whole ECHAM physics    [W/m2]

    ! Temporary array used by GW_HINES
    REAL(wp) :: zdis_gwh(nbdim,nlev)  !<  out, energy dissipation rate [J/s/kg]

    INTEGER  :: nc    !< number of cells/columns from (jce-jcs+1)
    REAL(wp) :: zlat_deg(nbdim)           !< latitude in deg N
    REAL(wp) :: zq_gwh (nbdim,nlev)       !< heating by atm. gravity waves     [W/m2]
 
    ! number of cells/columns from index jcs to jce
    nc = jce-jcs+1

    zlat_deg(jcs:jce) = patch%cells%center(jcs:jce,jb)%lat * 180._wp/pi

    CALL gw_hines ( jg                       ,&
      &             nbdim                    ,&
      &             jcs                      ,&
      &             jce                      ,&
      &             nc                       ,&
      &             nlev                     ,&
      &             field% presi_old(:,:,jb) ,&
      &             field% presm_old(:,:,jb) ,&
      &             field%   ta(:,:,jb)      ,&
      &             field%   ua(:,:,jb)      ,&
      &             field%   va(:,:,jb)      ,&
      &             zlat_deg(:)              ,&
!!$        &             aprflux(:,krow)          ,&
      &             zdis_gwh(:,:)            ,&
      &             tend%   ua_gwh(:,:,jb)   ,&
      &             tend%   va_gwh(:,:,jb) )


    ! heating
    zq_gwh(jcs:jce,:) = zdis_gwh(jcs:jce,:) * field%mair(jcs:jce,:,jb)

    ! heating accumulated
    zq_phy(jcs:jce,:) = zq_phy(jcs:jce,:) + zq_gwh(jcs:jce,:)

    ! tendency
    tend% ta_gwh(jcs:jce,:,jb) = zq_gwh(jcs:jce,:)*zconv(jcs:jce,:)

    ! tendencies accumulated
    tend%   ta(jcs:jce,:,jb) = tend%   ta(jcs:jce,:,jb) + tend%   ta_gwh(jcs:jce,:,jb)
    tend%   ua(jcs:jce,:,jb) = tend%   ua(jcs:jce,:,jb) + tend%   ua_gwh(jcs:jce,:,jb)
    tend%   va(jcs:jce,:,jb) = tend%   va(jcs:jce,:,jb) + tend%   va_gwh(jcs:jce,:,jb)

  END SUBROUTINE echam_gw_hinesg
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE echam_ssodrag(patch, jb,jcs,jce, nbdim, field, tend, zconv, zq_phy)
    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block 
    REAL(wp)        ,INTENT(IN) :: zconv  (nbdim,nlev)       !< specific heat of moist air        [J/K/kg]
    TYPE(t_echam_phy_field),   POINTER :: field
    TYPE(t_echam_phy_tend) ,   POINTER :: tend
    REAL(wp)        ,INTENT(INOUT) :: zq_phy (nbdim,nlev)       !< heating by whole ECHAM physics    [W/m2]

    ! Temporary array used by SSODRAG
    REAL(wp) :: zdis_sso(nbdim,nlev)  !<  out, energy dissipation rate [J/s/kg]
    REAL(wp) :: zq_sso (nbdim,nlev)       !< heating by subgrid scale orogr.   [W/m2]
    INTEGER :: nc

    ! number of cells/columns from index jcs to jce
    nc = jce-jcs+1

    CALL ssodrag( nc                                        ,& ! in,  number of cells/columns in loop (jce-jcs+1)
                  nbdim                                     ,& ! in,  dimension of block of cells/columns
                  nlev                                      ,& ! in,  number of levels
                  !
                  patch%cells%center(:,jb)%lat              ,& ! in,  Latitude in radians
                  pdtime                                    ,& ! in,  time step length
                  !
                  field% presi_old(:,:,jb)                  ,& ! in,  p at half levels
                  field% presm_old(:,:,jb)                  ,& ! in,  p at full levels
                  field% geom(:,:,jb)                       ,& ! in,  geopotential above surface (t-dt)
                  field%   ta(:,:,jb)                       ,& ! in,  T
                  field%   ua(:,:,jb)                       ,& ! in,  u
                  field%   va(:,:,jb)                       ,& ! in,  v
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
                  tend%   ua_sso(:,:,jb)                    ,& ! out, tendency of zonal wind
                  tend%   va_sso(:,:,jb)                     ) ! out, tendency of meridional wind


      ! heating
      zq_sso(jcs:jce,:) = zdis_sso(jcs:jce,:) * field%mair(jcs:jce,:,jb)

      ! heating accumulated
      zq_phy(jcs:jce,:) = zq_phy(jcs:jce,:) + zq_sso(jcs:jce,:)

      ! tendency
      tend% ta_sso(jcs:jce,:,jb) = zq_sso(jcs:jce,:)*zconv(jcs:jce,:)

      ! tendencies accumulated
      tend%   ta(jcs:jce,:,jb) = tend%   ta(jcs:jce,:,jb) + tend%   ta_sso(jcs:jce,:,jb)
      tend%   ua(jcs:jce,:,jb) = tend%   ua(jcs:jce,:,jb) + tend%   ua_sso(jcs:jce,:,jb)
      tend%   va(jcs:jce,:,jb) = tend%   va(jcs:jce,:,jb) + tend%   va_sso(jcs:jce,:,jb)

  END SUBROUTINE echam_ssodrag
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE echam_cumulus_condensation(jb,jcs,jce, nbdim, field, tend, zcair)
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block 
    REAL(wp)        ,INTENT(IN) :: zcair  (nbdim,nlev)       !< specific heat of moist air        [J/K/kg]
    TYPE(t_echam_phy_field),   POINTER :: field
    TYPE(t_echam_phy_tend) ,   POINTER :: tend

    ! local 
    REAL(wp) :: zqtec  (nbdim,nlev)       !< tracer tendency due to entrainment/detrainment
    INTEGER  :: itype(nbdim)              !< type of convection
    INTEGER  :: ictop (nbdim)             !< from massflux
    INTEGER  :: ilab   (nproma,nlev)

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
        &          pdtime,                    &! in
        &          field% lfland(:,jb),       &! in     loland
        &          field% ta(:,:,jb),         &! in     tm1
        &          field% ua(:,:,jb),         &! in     um1
        &          field% va(:,:,jb),         &! in     vm1
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
        &          tend% ta(:,:,jb),          &! in     tte  for internal updating
        &          tend% ua(:,:,jb),          &! in     vom  for internal updating
        &          tend% va(:,:,jb),          &! in     vol  for internal updating
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
        &          tend%  ta_cnv(:,:,jb),     &! out
        &          tend%  ua_cnv(:,:,jb),     &! out
        &          tend%  va_cnv(:,:,jb),     &! out
        &          tend%qtrc_cnv(:,:,jb,iqv), &! out
        &          tend%qtrc_cnv(:,:,jb,iqt:) )! out

      IF (ltimer) CALL timer_stop(timer_cucall)

      field% rtype(jcs:jce,jb) = REAL(itype(jcs:jce),wp)

!!$      ! heating accumulated
!!$      zq_phy(jcs:jce,:) = zq_phy(jcs:jce,:) + zq_cnv(jcs:jce,:)
!!$
!!$      ! tendency
!!$      tend% temp_cnv(jcs:jce,:,jb) = zq_cnv(jcs:jce,:)*zconv(jcs:jce,:)
!!$
      ! tendencies accumulated
      tend%   ua(jcs:jce,:,jb)      = tend%   ua(jcs:jce,:,jb)      + tend%   ua_cnv(jcs:jce,:,jb)
      tend%   va(jcs:jce,:,jb)      = tend%   va(jcs:jce,:,jb)      + tend%   va_cnv(jcs:jce,:,jb)
      tend%   ta(jcs:jce,:,jb)      = tend%   ta(jcs:jce,:,jb)      + tend%   ta_cnv(jcs:jce,:,jb)
      tend% qtrc(jcs:jce,:,jb,iqv)  = tend% qtrc(jcs:jce,:,jb,iqv)  + tend% qtrc_cnv(jcs:jce,:,jb,iqv)
      tend% qtrc(jcs:jce,:,jb,iqt:) = tend% qtrc(jcs:jce,:,jb,iqt:) + tend% qtrc_cnv(jcs:jce,:,jb,iqt:)


    ELSE ! NECESSARY COMPUTATIONS IF MASSFLUX IS BY-PASSED

      ilab(jcs:jce,1:nlev) = 0
      ictop(jcs:jce)       = nlev-1

      tend%   ua_cnv(jcs:jce,:,jb)      = 0._wp
      tend%   va_cnv(jcs:jce,:,jb)      = 0._wp
      tend%   ta_cnv(jcs:jce,:,jb)      = 0._wp
      tend% qtrc_cnv(jcs:jce,:,jb,iqv)  = 0._wp
      tend% qtrc_cnv(jcs:jce,:,jb,iqt:) = 0._wp

    ENDIF !lconv

    !-------------------------------------------------------------
    ! 7. LARGE SCALE CONDENSATION.
    !-------------------------------------------------------------
    IF(phy_config%lcond) THEN

      !IF (lcotra) CALL get_col_pol( tend%ta(:,:,jb),tend%qtrc(:,:,jb,iqv),jb )

      IF (ncdnc==0 .AND. nicnc==0) THEN

        field% rsfl(:,jb) = 0._wp
        field% ssfl(:,jb) = 0._wp

        IF (ltimer) CALL timer_start(timer_cloud)

        CALL cloud(jce, nproma, jks, nlev, nlevp1, &! in
          &        pdtime,                    &! in
          &        ictop,                     &! in (from "cucall")
          &        field% presi_old(:,:,jb),  &! in
          &        field% presm_old(:,:,jb),  &! in
!          &        field% presm_new(:,:,jb), &! in
          &        field% acdnc (:,:,jb),     &! in. acdnc
          &        field%   ta  (:,:,jb),     &! in. tm1
          &        field%   tv  (:,:,jb),     &! in. ztvm1
          &        field% qtrc  (:,:,jb,iqv), &! in.  qm1
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
          &         tend% ta  (:,:,jb),     &! inout.  tte
          &         tend% qtrc  (:,:,jb,iqv), &! inout.  qte
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
          &        tend%  ta_cld(:,:,jb),     &! out
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

      tend%   ta_cld(jcs:jce,:,jb)      = 0._wp
      tend% qtrc_cld(jcs:jce,:,jb,iqv)  = 0._wp
      tend% qtrc_cld(jcs:jce,:,jb,iqc)  = 0._wp
      tend% qtrc_cld(jcs:jce,:,jb,iqi)  = 0._wp
      tend% qtrc_cld(jcs:jce,:,jb,iqt:) = 0._wp

    ENDIF !lcond

  END SUBROUTINE echam_cumulus_condensation
  !---------------------------------------------------------------------

END MODULE mo_echam_phy_main
