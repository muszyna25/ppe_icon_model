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
  USE mo_physical_constants,  ONLY: cpd, cpv, cvd, cvv, &
    &                               amd, amo3
  USE mo_impl_constants      ,ONLY: inh_atmosphere, min_rlcell_int, grf_bdywidth_c
  USE mo_run_config,          ONLY: ntracer, nlev, nlevm1, nlevp1,    &
    &                               iqv, iqc, iqi, iqt, io3
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
  USE mtime,                  ONLY: datetime
  USE mo_ham_aerosol_params,  ONLY: ncdnc, nicnc
  USE mo_echam_sfc_indices,   ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_surface,             ONLY: update_surface
  USE mo_surface_diag,        ONLY: nsurf_diag
  USE mo_cloud,               ONLY: cloud
  USE mo_cover,               ONLY: cover
  USE mo_radheating,          ONLY: radheating
  USE mo_psrad_radiation,     ONLY: psrad_radiation
  USE mo_psrad_radiation_parameters, ONLY: psctm
  USE mo_vdiff_config,        ONLY: vdiff_config
  USE mo_vdiff_downward_sweep,ONLY: vdiff_down
  USE mo_vdiff_upward_sweep,  ONLY: vdiff_up
  USE mo_vdiff_solver,        ONLY: nvar_vdiff, nmatrix, imh, imqv,   &
    &                               ih_vdiff=>ih, iqv_vdiff=>iqv
  USE mo_gw_hines,            ONLY: gw_hines
  USE mo_ssortns,             ONLY: ssodrag
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights, &
    &                                  calculate_time_interpolation_weights
  USE mo_echam_radiation,     ONLY: echam_radiation

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_lcariolle_types,     ONLY: t_avi, t_time_interpolation
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: echam_phy_main

  ! variables that are time-independent, potentailly  shared by the routines of this module 
  ! they could have been calculated at initializatiop phase
  INTEGER  :: ntrac !< # of tracers excluding water vapour and hydrometeors
                      !< (handled by sub-models, e.g., chemical species)
  INTEGER  :: jks   !< start index for vertical loops
  REAL(wp) :: zcd                       !< specific heat of dry air          [J/K/kg]
  REAL(wp) :: zcv                       !< specific heat of water vapor      [J/K/kg]

  ! the following would depend on the nesting, ie jg
  REAL(wp) :: pdtime         !< time step

CONTAINS
  !>
  !!
  SUBROUTINE echam_phy_main( patch,                         &
    &                        rl_start, rl_end,              &
    &                        this_datetime,inpdtime,        &
    &                        ltrig_rad                      )

    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(IN)  :: rl_start, rl_end

    TYPE(datetime),  POINTER     :: this_datetime    !< date and time 
    REAL(wp)        ,INTENT(IN)  :: inpdtime         !< time step

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
    REAL(wp) :: zq_rsw (nproma,nlev)       !< heating by short wave radiation   [W/m2]
    REAL(wp) :: zq_rlw (nproma,nlev,patch%nblks_c)       !< heating by long  wave radiation   [W/m2]
    REAL(wp) :: zq_rlw_impl (nproma)       !< additional heating by LW rad. due to impl. coupling in surface energy balance [W/m2]

    REAL(wp) :: zxt_emis(nproma,ntracer-iqt+1)  !< tracer tendency due to surface emission
                                               !< and dry deposition. "zxtems" in ECHAM5
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
    pdtime = inpdtime
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
    IF (phy_config%lrad) THEN

       ! 4.1 RADIATIVE TRANSFER
       !-----------------------
       IF (ltrig_rad) THEN

          ! store ts_rad of this radiatiative transfer timestep in ts_rad_rt,
          ! so that it can be reused in radheat in the other timesteps
          field%ts_rad_rt(:,:) = field%ts_rad(:,:)

          CALL echam_radiation( patch, jg, field, this_datetime)

      END IF ! ltrig_rad
    END IF ! phy_config%lrad

!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      CALL calculate_zcair_zconv( jb,jcs,jce, nproma, field, zcair(:,:,jb), zconv(:,:,jb))
    ENDDO
!$OMP END PARALLEL DO 

    ! initialize physics heating
    zq_phy(:,:,:) = 0._wp

    IF (phy_config%lrad) THEN
      IF (ltimer) CALL timer_start(timer_radheat)
      ! 4.2 RADIATIVE HEATING
      !----------------------
      ! radheat first computes the shortwave and longwave radiation for the current time step from transmissivity and
      ! the longwave flux at the radiation time step and, from there, the radiative heating due to sw and lw radiation.
      ! If radiation is called every time step, the longwave flux is not changed.

!$OMP PARALLEL DO PRIVATE(jcs,jce, zq_rsw)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

        CALL echam_radheating( jg, jb,jcs,jce, nproma, field, zq_rsw, zq_rlw(:,:,jb))

        ! heating accumulated
        zq_phy(jcs:jce,:,jb) = zq_phy(jcs:jce,:,jb) + zq_rsw(jcs:jce,:) + zq_rlw(jcs:jce,:,jb)
        
        ! tendencies
        tend%ta_rsw(jcs:jce,:,jb) = zq_rsw(jcs:jce,:)    * zconv(jcs:jce,:,jb)
        tend%ta_rlw(jcs:jce,:,jb) = zq_rlw(jcs:jce,:,jb) * zconv(jcs:jce,:,jb)

        ! tendencies accumulated
        tend% ta(jcs:jce,:,jb) = tend% ta     (jcs:jce,:,jb) &
          &                    + tend% ta_rsw (jcs:jce,:,jb) &
          &                    + tend% ta_rlw (jcs:jce,:,jb)


      ENDDO
!$OMP END PARALLEL DO 
      IF (ltimer) CALL timer_stop(timer_radheat)

    ELSE   
      ! If computation of radiative heating is by-passed
      ! this shound not be needed, as these field are zeroed in the initialization
!$OMP PARALLEL DO PRIVATE(jcs,jce)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

        tend%ta_rsw(jcs:jce,:,jb) = 0.0_wp
        tend%ta_rlw(jcs:jce,:,jb) = 0.0_wp

        field%rsdt(jcs:jce,jb)= 0.0_wp
      ENDDO
!$OMP END PARALLEL DO 

    END IF ! lrad
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! 5. BOUNDARY LAYER AND SURFACE PROCESSES
    !-------------------------------------------------------------------
    ! Note: In ECHAM this part is located between "CALL radiation" and
    !       "CALL radheat".
    !

    IF (phy_config%lvdiff) THEN

!$OMP PARALLEL DO PRIVATE(jcs,jce, zxt_emis)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

        ! 5.1 Emission of aerosols or other tracers. Not implemented yet.
    ! 
    !       IF (ntrac>0) THEN
    !         !CALL tracer_emission()
            zxt_emis(jcs:jce,:) = 0._wp
    !       ENDIF
        !
        ! 5.2 Dry deposition of aerosols or other tracers. Not implemented yet.
        ! CALL dry_deposition()
        !

        CALL echam_vdiffDownUp_surf( jg, jb,jcs,jce, nproma, zfrc(:,:,jb), zxt_emis, field, tend)

      ENDDO
!$OMP END PARALLEL DO 

    ELSE

!$OMP PARALLEL DO PRIVATE(jcs,jce)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
        field% evap(jcs:jce,jb)= 0._wp

        tend%   ua_vdf(jcs:jce,:,jb)      = 0._wp
        tend%   va_vdf(jcs:jce,:,jb)      = 0._wp
        tend%   ta_vdf(jcs:jce,:,jb)      = 0._wp
        tend% qtrc_vdf(jcs:jce,:,jb,iqv)  = 0._wp
        tend% qtrc_vdf(jcs:jce,:,jb,iqc)  = 0._wp
        tend% qtrc_vdf(jcs:jce,:,jb,iqi)  = 0._wp
        tend% qtrc_vdf(jcs:jce,:,jb,iqt:) = 0._wp
      ENDDO
!$OMP END PARALLEL DO 

    ENDIF !lvdiff
    !-----------------------

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
      CALL  echam_lcariolle(patch, rl_start, rl_end, field, this_datetime, tend) 
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

        CALL echam_gw_hinesg(patch, jg, jb,jcs,jce, nproma, zconv(:,:,jb), field, tend, zq_phy(:,:,jb))
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

        CALL echam_ssodrag(patch, jb,jcs,jce, nproma, zconv(:,:,jb), field, tend, zq_phy(:,:,jb))
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

      CALL echam_cumulus_condensation(jb,jcs,jce, nproma, zcair(:,:,jb), field, tend)

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
  SUBROUTINE echam_radheating( jg, jb,jcs,jce, nbdim, field, zq_rsw, zq_rlw)
    INTEGER         ,INTENT(IN) :: jg
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block
    TYPE(t_echam_phy_field),   POINTER :: field    ! in

    REAL(wp)        ,INTENT(INOUT) :: zq_rsw (nbdim,nlev)       !< out, heating by short wave radiation   [W/m2]
    REAL(wp)        ,INTENT(INOUT) :: zq_rlw (nbdim,nlev)      !<  out, heating by long  wave radiation   [W/m2]

    CALL radheating (                                &
      !
      ! input
      ! -----
      !
      & jcs        = jcs                            ,&! loop start index
      & jce        = jce                            ,&! loop end index
      & kbdim      = nbdim                          ,&! dimension size
      & klev       = nlev                           ,&! vertical dimension size
      & klevp1     = nlevp1                         ,&! vertical dimension size
      !
      & rsdt0      = psctm                          ,&! toa incident shortwave radiation for sun in zenith
      & cosmu0     = field%cosmu0    (:,jb)         ,&! solar zenith angle at current time
      !
      & emiss      = ext_data(jg)%atm%emis_rad(:,jb),&! lw sfc emissivity
      & tsr        = field%ts_rad (:,jb)            ,&! radiative surface temperature at current   time [K]
      & tsr_rt     = field%ts_rad_rt(:,jb)          ,&! radiative surface temperature at radiation time [K]
      !
      & rsd_rt     = field%rsd_rt           (:,:,jb),&! all-sky   shortwave downward flux at radiation time [W/m2]
      & rsu_rt     = field%rsu_rt           (:,:,jb),&! all-sky   shortwave upward   flux at radiation time [W/m2]
      !
      & rsdcs_rt   = field%rsdcs_rt         (:,:,jb),&! clear-sky shortwave downward flux at radiation time [W/m2]
      & rsucs_rt   = field%rsucs_rt         (:,:,jb),&! clear-sky shortwave upward   flux at radiation time [W/m2]
      !
      & rld_rt     = field%rld_rt           (:,:,jb),&! all-sky   longwave  downward flux at radiation time [W/m2]
      & rlu_rt     = field%rlu_rt           (:,:,jb),&! all-sky   longwave  upward   flux at radiation time [W/m2]
      !
      & rldcs_rt   = field%rldcs_rt         (:,:,jb),&! clear-sky longwave  downward flux at radiation time [W/m2]
      & rlucs_rt   = field%rlucs_rt         (:,:,jb),&! clear-sky longwave  upward   flux at radiation time [W/m2]
      !
      & rvds_dir_rt= field%rvds_dir_rt        (:,jb),&!< out  all-sky downward direct visible radiation at surface
      & rpds_dir_rt= field%rpds_dir_rt        (:,jb),&!< out  all-sky downward direct PAR     radiation at surface
      & rnds_dir_rt= field%rnds_dir_rt        (:,jb),&!< out  all-sky downward direct near-IR radiation at surface
      & rvds_dif_rt= field%rvds_dif_rt        (:,jb),&!< out  all-sky downward diffuse visible radiation at surface
      & rpds_dif_rt= field%rpds_dif_rt        (:,jb),&!< out  all-sky downward diffuse PAR     radiation at surface
      & rnds_dif_rt= field%rnds_dif_rt        (:,jb),&!< out  all-sky downward diffuse near-IR radiation at surface
      & rvus_rt    = field%rvus_rt            (:,jb),&!< out  all-sky upward visible radiation at surface
      & rpus_rt    = field%rpus_rt            (:,jb),&!< out  all-sky upward PAR     radiation at surfac
      & rnus_rt    = field%rnus_rt            (:,jb),&!< out  all-sky upward near-IR radiation at surface
      !
      ! output
      ! ------
      !
      & rsdt       = field%rsdt               (:,jb),&! all-sky   shortwave downward flux at current   time [W/m2]
      & rsut       = field%rsut               (:,jb),&! all-sky   shortwave upward   flux at current   time [W/m2]
      & rsds       = field%rsds               (:,jb),&! all-sky   shortwave downward flux at current   time [W/m2]
      & rsus       = field%rsus               (:,jb),&! all-sky   shortwave upward   flux at current   time [W/m2]
      !
      & rsutcs     = field%rsutcs             (:,jb),&! clear-sky shortwave upward   flux at current   time [W/m2]
      & rsdscs     = field%rsdscs             (:,jb),&! clear-sky shortwave downward flux at current   time [W/m2]
      & rsuscs     = field%rsuscs             (:,jb),&! clear-sky shortwave upward   flux at current   time [W/m2]
      !
      & rvds_dir   = field%rvds_dir           (:,jb),&!< out  all-sky downward direct visible radiation at surface
      & rpds_dir   = field%rpds_dir           (:,jb),&!< out  all-sky downward direct PAR     radiation at surface
      & rnds_dir   = field%rnds_dir           (:,jb),&!< out  all-sky downward direct near-IR radiation at surface
      & rvds_dif   = field%rvds_dif           (:,jb),&!< out  all-sky downward diffuse visible radiation at surface
      & rpds_dif   = field%rpds_dif           (:,jb),&!< out  all-sky downward diffuse PAR     radiation at surface
      & rnds_dif   = field%rnds_dif           (:,jb),&!< out  all-sky downward diffuse near-IR radiation at surface
      & rvus       = field%rvus               (:,jb),&!< out  all-sky upward visible radiation at surface
      & rpus       = field%rpus               (:,jb),&!< out  all-sky upward PAR     radiation at surfac
      & rnus       = field%rnus               (:,jb),&!< out  all-sky upward near-IR radiation at surface
      !
      & rlut       = field%rlut               (:,jb),&! all-sky   longwave  upward   flux at current   time [W/m2]
      & rlds       = field%rlds               (:,jb),&! all-sky   longwave  downward flux at current   time [W/m2]
      & rlus       = field%rlus               (:,jb),&! all-sky   longwave  upward   flux at current   time [W/m2]
      !
      & rlutcs     = field%rlutcs             (:,jb),&! clear-sky longwave  upward   flux at current   time [W/m2]
      & rldscs     = field%rldscs             (:,jb),&! clear-sky longwave  downward flux at current   time [W/m2]
      !
      & q_rsw      = zq_rsw                   (:,:) ,&! rad. heating by SW           [W/m2]
      & q_rlw      = zq_rlw                   (:,:) ) ! rad. heating by LW           [W/m2]

  END SUBROUTINE echam_radheating
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! 5.3 Turbulent mixing, part I:
  !     computation of exchange coefficients in the atmosphere and at the surface;
  !     build up the tridiagonal linear algebraic system;
  !     downward sweep (Gaussian elimination from top till level nlev-1)
  SUBROUTINE echam_vdiffDownUp_surf(jg, jb,jcs,jce, nbdim, zfrc, zxt_emis, field, tend)
    INTEGER         ,INTENT(IN) :: jg             
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block
    REAL(wp),        INTENT(IN) :: zfrc (nbdim,nsfc_type)    !< zfrl, zfrw, zfrc combined, output
    REAL(wp)        ,INTENT(IN) :: zxt_emis(nbdim,ntracer-iqt+1)  !< tracer tendency due to surface emission
    TYPE(t_echam_phy_field),   POINTER :: field    ! inout
    TYPE(t_echam_phy_tend) ,   POINTER :: tend


    REAL(wp) :: zdelp  (nbdim,nlev)       !< layer thickness in pressure coordinate  [Pa]
    INTEGER  :: ihpbl  (nbdim)            !< location of PBL top given as vertical level index
    REAL(wp) :: zcpt_sfc_tile(nbdim,nsfc_type)  !< dry static energy at surface
   
    ! Coefficient matrices and right-hand-side vectors for the turbulence solver
    ! _btm refers to the lowest model level (i.e., full level "klev", not the surface)
    REAL(wp) :: zri_tile(nbdim,nsfc_type) !< Richardson number
    REAL(wp) :: zaa    (nbdim,nlev,3,nmatrix)       !< coeff. matrices, all variables
    REAL(wp) :: zaa_btm(nbdim,3,nsfc_type,imh:imqv) !< last row of coeff. matrix of heat and moisture
    REAL(wp) :: zbb    (nbdim,nlev,nvar_vdiff)  !< r.h.s., all variables
    REAL(wp) :: zbb_btm(nbdim,nsfc_type,ih_vdiff:iqv_vdiff) !< last row of r.h.s. of heat and moisture
    ! Temporary arrays used by VDIFF, JSBACH

    REAL(wp) :: zfactor_sfc(nbdim)

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

    REAL(wp) :: zvmixtau   (nproma,nlev)   !< timescale of mixing for vertical turbulence
    REAL(wp) :: zqtvar_prod(nproma,nlev)   !< production rate of total water variance
                                          !< due to turbulence. Computed in "vdiff",
                                          !< used by "cloud"

    INTEGER  :: jc, jk

      IF (ltimer) CALL timer_start(timer_vdiff_down)

      DO jk = 1,nlev
        DO jc = jcs,jce
          ! 3.2 Thickness of model layer in pressure coordinate
           zdelp   (jc,jk) = field% presi_old (jc,jk+1,jb) - field% presi_old (jc,jk,jb)
        END DO
      END DO

      CALL vdiff_down( vdiff_config%lsfc_mom_flux,      &! in
                     & vdiff_config%lsfc_heat_flux,     &! in
                     & jce, nbdim, nlev, nlevm1, nlevp1,&! in
                     & ntrac, nsfc_type,                &! in
                     & iwtr, iice, ilnd,                &! in, indices of different surface types
                     & pdtime,                          &! in, time step 
                     & field%coriol(:,jb),              &! in, Coriolis parameter
                     & zfrc(:,:),                    &! in, area fraction of each sfc type
                     & field% ts_tile(:,jb,:),          &! in, surface temperature
                     & field% ocu (:,jb),               &! in, ocean sfc velocity, u-component
                     & field% ocv (:,jb),               &! in, ocean sfc velocity, v-component
                     & field% presi_old(:,nlevp1,jb),   &! in, sfc pressure
                     & field%   ua(:,:,jb),             &! in, um1
                     & field%   va(:,:,jb),             &! in, vm1
                     & field%   ta(:,:,jb),             &! in, tm1
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
          & pdtime,                       &! in, time step
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
          & field% ts_tile(:,jb,:),       &! inout
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
          & pu    = field% ua(:,nlev,jb), &! in, um1
          & pv    = field% va(:,nlev,jb), &! in, vm1
          & ptemp = field% ta(:,nlev,jb), &! in, tm1
          & pq = field% qtrc(:,nlev,jb,iqv),  &! in, qm1
          & prsfl = field% rsfl(:,jb),    &! in, rain surface large scale (from cloud)
          & prsfc = field% rsfc(:,jb),    &! in, rain surface concective (from cucall)
          & pssfl = field% ssfl(:,jb),    &! in, snow surface large scale (from cloud)
          & pssfc = field% ssfc(:,jb),    &! in, snow surface concective (from cucall)
          & rlds        = field% rlds (:,jb), &! in,  downward surface  longwave flux [W/m2]
          & rlus        = field% rlus (:,jb), &! inout, upward surface  longwave flux [W/m2]
          & rsds        = field% rsds (:,jb), &! in,  downward surface shortwave flux [W/m2]
          & rsus        = field% rsus (:,jb), &! inout, upward surface shortwave flux [W/m2]
          !
          & rvds_dir   = field%rvds_dir   (:,jb), &! in, all-sky downward direct visible radiation at surface
          & rpds_dir   = field%rpds_dir   (:,jb), &! in, all-sky downward direct PAR     radiation at surface
          & rnds_dir   = field%rnds_dir   (:,jb), &! in, all-sky downward direct near-IR radiation at surface
          & rvds_dif   = field%rvds_dif   (:,jb), &! in, all-sky downward diffuse visible radiation at surface
          & rpds_dif   = field%rpds_dif   (:,jb), &! in, all-sky downward diffuse PAR     radiation at surface
          & rnds_dif   = field%rnds_dif   (:,jb), &! in, all-sky downward diffuse near-IR radiation at surface
          !
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
          & ptsfc     = field%ts    (:,jb),                        &! out
          & ptsfc_rad = field%ts_rad(:,jb),                        &! out
          & rlns_tile = field%lwflxsfc_tile(:,jb,:),               &! out (for coupling)
          & rsns_tile = field%swflxsfc_tile(:,jb,:),               &! out (for coupling)
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
                   & pdtime,                          &! in, time step
                   & zfrc(:,:),                       &! in, area fraction of each sfc type
                   & field% cfm_tile(:,jb,:),         &! in
                   & zaa,                             &! in, from "vdiff_down"
                   &   ihpbl(:),                      &! in, from "vdiff_down"
                   &  zcptgz(:,:),                    &! in, from "vdiff_down"
                   &   zrhoh(:,:),                    &! in, from "vdiff_down"
                   & zqshear(:,:),                    &! in, from "vdiff_down"
                   & field%   ua(:,:,jb),             &! in, um1
                   & field%   va(:,:,jb),             &! in, vm1
                   & field%   ta(:,:,jb),             &! in, tm1
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
                   &  tend%   ua_vdf(:,:,jb),         &! out
                   &  tend%   va_vdf(:,:,jb),         &! out
                   &  tend%   ta_vdf(:,:,jb),         &! out
!!$                   &          zq_vdf(:,:),            &! out   heating W/m2
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

!!$      ! heating accumulated
!!$      zq_phy(jcs:jce,:) = zq_phy(jcs:jce,:) + zq_vdf(jcs:jce,:)
!!$
!!$      ! tendency
!!$      tend% temp_vdf(jcs:jce,:,jb) = zq_vdf(jcs:jce,:)*zconv(jcs:jce,:)
!!$
      ! tendencies accumulated
      tend%   ua(jcs:jce,:,jb)      = tend%   ua(jcs:jce,:,jb)      + tend%   ua_vdf(jcs:jce,:,jb)
      tend%   va(jcs:jce,:,jb)      = tend%   va(jcs:jce,:,jb)      + tend%   va_vdf(jcs:jce,:,jb)
      tend%   ta(jcs:jce,:,jb)      = tend%   ta(jcs:jce,:,jb)      + tend%   ta_vdf(jcs:jce,:,jb)
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

    ! 2-tl-scheme
    field% tkem1(jcs:jce,:,jb) = field% tke  (jcs:jce,:,jb)

    ! 5.6 Turbulent mixing, part III:
    !     - Further diagnostics.

    CALL nsurf_diag( jce, nbdim, nsfc_type,           &! in
                   & ilnd,                            &! in
                   & zfrc(:,:),                       &! in
                   & field%  qtrc(:,nlev,jb,iqv),     &! in humidity qm1
                   & field%    ta(:,nlev,jb),         &! in tm1
                   & field% presm_old(:,nlev,jb),     &! in, apm1
                   & field% presi_old(:,nlevp1,jb),   &! in, aphm1
                   & field%   qx(:,nlev,jb),          &! in, xlm1 + xim1
                   & field%   ua(:,nlev,jb),          &! in, um1
                   & field%   va(:,nlev,jb),          &! in, vm1
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

  END SUBROUTINE echam_vdiffDownUp_surf
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE echam_lcariolle(patch, rl_start, rl_end, field, this_datetime, tend)
    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(in) :: rl_start, rl_end
    TYPE(t_echam_phy_field),   POINTER :: field
    TYPE(t_echam_phy_tend) ,   POINTER :: tend
    TYPE(datetime),  POINTER     :: this_datetime  !< time step

     ! Temporary variables for Cariolle scheme (ozone)
    REAL(wp)    :: do3dt(nproma,nlev)
    TYPE(t_time_interpolation) :: time_interpolation
    EXTERNAL       lcariolle_lat_intp_li, lcariolle_pres_intp_li
    TYPE(t_time_interpolation_weights) :: current_time_interpolation_weights
    TYPE(t_avi) :: avi

    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block
 
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    current_time_interpolation_weights = calculate_time_interpolation_weights(this_datetime)
    time_interpolation%imonth1=current_time_interpolation_weights%month1_index
    time_interpolation%imonth2=current_time_interpolation_weights%month2_index
    time_interpolation%weight1=current_time_interpolation_weights%weight1
    time_interpolation%weight2=current_time_interpolation_weights%weight2

!  NOTE: something is wrong with the avi dimensions; cannot be run in OpenMP
!$OMP PARALLEL PRIVATE(avi)
    ALLOCATE(avi%o3_vmr(nproma,nlev), avi%vmr2molm2(nproma,nlev), avi%cell_center_lat(nproma), &
        & avi%lday(nproma))
!$OMP DO PRIVATE(jcs,jce,do3dt)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
 
      avi%ldown=.TRUE.
      avi%o3_vmr(jcs:jce,:)        = field%qtrc(jcs:jce,:,jb,io3)*amd/amo3
      avi%tmprt                    => field%ta(:,:,jb)
      avi%vmr2molm2(jcs:jce,:)     = field%mdry(jcs:jce,:,jb) / amd * 1.e3_wp
      avi%pres                     => field%presm_old(jcs:jce,:,jb)
      avi%cell_center_lat(jcs:jce) = patch%cells%center(jcs:jce,jb)%lat
      avi%lday(jcs:jce)            = field%cosmu0(jcs:jce,jb) > 1.e-3_wp

      CALL lcariolle_do3dt(                                                    &
          & jcs,                    jce,                nproma,                  &
          & nlev,                   time_interpolation, lcariolle_lat_intp_li,  &
          & lcariolle_pres_intp_li, avi,                do3dt                   )
      tend% qtrc(jcs:jce,:,jb,io3) = tend% qtrc(jcs:jce,:,jb,io3) + do3dt(jcs:jce,:)*amo3/amd
    ENDDO
!$OMP END DO 
    DEALLOCATE(avi%o3_vmr, avi%vmr2molm2, avi%cell_center_lat, avi%lday)
!$OMP END PARALLEL

  END SUBROUTINE echam_lcariolle
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE echam_gw_hinesg(patch, jg, jb,jcs,jce, nbdim, zconv, field, tend, zq_phy)
    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(IN) :: jg
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block  !---------------------------------------------------------------------
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
  SUBROUTINE echam_ssodrag(patch, jb,jcs,jce, nbdim, zconv, field, tend, zq_phy)
    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block  !---------------------------------------------------------------------
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
  SUBROUTINE echam_cumulus_condensation(jb,jcs,jce, nbdim, zcair, field, tend)
    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block  !---------------------------------------------------------------------
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
