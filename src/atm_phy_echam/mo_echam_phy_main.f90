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
  USE mo_exception,           ONLY: message
  USE mo_physical_constants,  ONLY: cpd, cpv, cvd, cvv, Tf, tmelt
  USE mo_run_config,          ONLY: ntracer, nlev, nlevp1,    &
    &                               iqv, iqc, iqi, iqt
  USE mo_mpi_phy_config,      ONLY: mpi_phy_tc, dt_zero
  USE mo_echam_cloud_config,  ONLY: echam_cloud_config
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field, prm_field,     &
    &                               t_echam_phy_tend,  prm_tend,      &
    &                               cdimissval
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop,  &
    &                               timer_cover
  USE mtime,                  ONLY: datetime, isCurrentEventActive, &
    &                               OPERATOR(<=), OPERATOR(>)
  USE mo_echam_sfc_indices,   ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_cover,               ONLY: cover

  USE mo_interface_echam_radiation,    ONLY: interface_echam_radiation
  USE mo_interface_echam_radheating,   ONLY: interface_echam_radheating
  USE mo_interface_echam_vdiff_surface,ONLY: interface_echam_vdiff_surface
  USE mo_interface_echam_o3_cariolle,  ONLY: interface_echam_o3_cariolle
  USE mo_interface_echam_convection,   ONLY: interface_echam_convection
  USE mo_interface_echam_gwhines,      ONLY: interface_echam_gwhines
  USE mo_interface_echam_sso,          ONLY: interface_echam_sso
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
  SUBROUTINE echam_phy_main( patch,            &
    &                        rl_start, rl_end, &
    &                        datetime_old,     &
    &                        pdtime            )

    TYPE(t_patch)  ,TARGET ,INTENT(in) :: patch
    INTEGER                ,INTENT(in) :: rl_start, rl_end

    TYPE(datetime)         ,POINTER    :: datetime_old    !< old date and time (at start of this step)
    REAL(wp)               ,INTENT(in) :: pdtime          !< time step

    ! Local variables
    TYPE(t_echam_phy_field),   POINTER :: field
    TYPE(t_echam_phy_tend) ,   POINTER :: tend

    REAL(wp) :: zq_rlw_impl (nproma)                     !< additional heating by LW rad. due to implicit coupling
                                                         !  in surface energy balance [W/m2]

    INTEGER  :: ictop (nproma,patch%nblks_c)             !< from massflux

    INTEGER  :: jg                                       !< grid level/domain index
    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb                                       !< block index
    INTEGER  :: jcs, jce                                 !< start/end column index within this block

    LOGICAL  :: is_in_sd_ed_interval                     !< time is in process interval [sd,ed[
    LOGICAL  :: is_active                                !< process is active

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

    ictop(:,:)   = nlev-1
    ! provisionally copy the incoming tedencies
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      tend%   ta_phy (jcs:jce,:,jb)   = tend%   ta (jcs:jce,:,jb)
      tend%   ua_phy (jcs:jce,:,jb)   = tend%   ua (jcs:jce,:,jb)
      tend%   va_phy (jcs:jce,:,jb)   = tend%   va (jcs:jce,:,jb)
      tend% qtrc_phy (jcs:jce,:,jb,:) = tend% qtrc (jcs:jce,:,jb,:)
    END DO
!$OMP END PARALLEL DO 

    ! initialize physics accumulated heating
    field% q_phy(:,:,:) = 0._wp

    !------------------------------------------------------------
    ! 3. COMPUTE SOME FIELDS NEEDED BY THE PHYSICAL ROUTINES.
    !------------------------------------------------------------

    ! 3.2b Specific heat of moist air
    !
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      CALL calculate_cpair_cvair_qconv(jb,jcs,jce, field)
    END DO
!$OMP END PARALLEL DO 

    !-------------------------------------------------------------------
    ! Calculate surface fraction 
    !------------------------------------------------------------------- 
!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
      CALL surface_fractions( jb,jcs,jce, nproma, field)
    END DO
!$OMP END PARALLEL DO 
 
    !-------------------------------------------------------------------
    ! Cloud cover (diagnostic)
    !------------------------------------------------------------------- 
    IF ( mpi_phy_tc(jg)%dt_cld > dt_zero ) THEN
       !
       IF ( (mpi_phy_tc(jg)%sd_cld <= datetime_old) .AND. &
            & (mpi_phy_tc(jg)%ed_cld >  datetime_old) .AND. &
            & (mpi_phy_tc(jg)%dt_cld >  dt_zero     )       ) THEN
          !
!$OMP PARALLEL DO PRIVATE(jcs,jce)
          DO jb = i_startblk,i_endblk
             CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
             CALL cloud_cover(jb,jcs,jce, nproma, field)
          END DO
!$OMP END PARALLEL DO 
          !
       END IF
       !
    END IF

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
    END DO
!$OMP END PARALLEL DO 


    !-------------------------------------------------------------------
    ! Radiation (one interface for LW+SW)
    !-------------------------------------------------------------------
    !
    IF ( (mpi_phy_tc(jg)%dt_rad > dt_zero) ) THEN
       !
       is_in_sd_ed_interval =          (mpi_phy_tc(jg)%sd_rad <= datetime_old) .AND. &
            &                          (mpi_phy_tc(jg)%ed_rad >  datetime_old)
       is_active = isCurrentEventActive(mpi_phy_tc(jg)%ev_rad,   datetime_old)
       !
       CALL message_forcing_action('LW and SW radiation (rad)',  &
            &                      is_in_sd_ed_interval,is_active)
       !
       ! radiative fluxes
       CALL interface_echam_radiation(is_in_sd_ed_interval,     &
            &                         is_active,                &
            &                         patch,                    &
            &                         field,                    &
            &                         datetime_old              )
       !
       ! radiative heating
       CALL interface_echam_radheating(is_in_sd_ed_interval,    &
            &                          patch, rl_start, rl_end, &
            &                          field, tend              )
       !
    END IF


    !-------------------------------------------------------------------
    ! Vertical diffusion, boundary layer and surface
    !-------------------------------------------------------------------
    !
    IF ( mpi_phy_tc(jg)%dt_vdf >  dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (mpi_phy_tc(jg)%sd_vdf <= datetime_old) .AND. &
            &                          (mpi_phy_tc(jg)%ed_vdf >  datetime_old)
       is_active = isCurrentEventActive(mpi_phy_tc(jg)%ev_vdf,   datetime_old)
       !
       CALL message_forcing_action('vertical diffusion (vdf)',   &
            &                      is_in_sd_ed_interval,is_active)
       !
       CALL interface_echam_vdiff_surface(is_in_sd_ed_interval,     &
            &                             is_active,                &
            &                             patch,  rl_start, rl_end, &
            &                             field,  tend,             &
            &                             pdtime                    )
       !
    END IF
    !
    !-----------------------

    !-----------------------
    IF ( mpi_phy_tc(jg)%dt_rad > dt_zero ) THEN
!$OMP PARALLEL DO PRIVATE(jcs,jce, zq_rlw_impl)
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
        ! Heating due to the fact that surface model only used part of longwave radiation to compute new surface temperature
        zq_rlw_impl(jcs:jce) =                                                  &
          &   ( (field%rld_rt(jcs:jce,nlev,jb)-field%rlu_rt(jcs:jce,nlev,jb))   & ! (( rlns from "radiation"
          &    -(field%rlds  (jcs:jce,jb)     -field%rlus  (jcs:jce,jb)     ))  & !   -rlns from "radheating" and "update_surface")
          &  -field%q_rlw(jcs:jce,nlev,jb)                                                ! old heating from radheat

        ! Heating accumulated
        field% q_phy(jcs:jce,nlev,jb) = field% q_phy(jcs:jce,nlev,jb) + zq_rlw_impl(jcs:jce)

        ! Tendency
        tend%ta_rlw_impl(jcs:jce,jb) = zq_rlw_impl(jcs:jce) * field% qconv(jcs:jce,nlev,jb)

        ! Tendencies accumulated
        tend%ta(jcs:jce,nlev,jb) = tend%ta(jcs:jce,nlev,jb) + tend%ta_rlw_impl(jcs:jce,jb)
      END DO
!$OMP END PARALLEL DO 

    END IF
    !---------------------

    !-------------------------------------------------------------------
    ! Linearized ozone chemistry of Cariolle
    !-------------------------------------------------------------------
    !
    IF ( mpi_phy_tc(jg)%dt_car > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (mpi_phy_tc(jg)%sd_car <= datetime_old) .AND. &
            &                          (mpi_phy_tc(jg)%ed_car >  datetime_old)
       is_active = isCurrentEventActive(mpi_phy_tc(jg)%ev_car,   datetime_old)
       !
       CALL message_forcing_action('lin. Cariolle ozone chem. (car)', &
            &                      is_in_sd_ed_interval,is_active     )
       !
       CALL  interface_echam_o3_cariolle(is_in_sd_ed_interval,    &
            &                            is_active,               &
            &                            patch, rl_start, rl_end, &
            &                            field, tend,             &
            &                            datetime_old             )
       !
    END IF

    !-------------------------------------------------------------------
    ! Atmospheric gravity wave drag
    !-------------------------------------------------------------------
    !
    IF ( mpi_phy_tc(jg)%dt_gwd > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (mpi_phy_tc(jg)%sd_gwd <= datetime_old) .AND. &
            &                          (mpi_phy_tc(jg)%ed_gwd >  datetime_old)
       is_active = isCurrentEventActive(mpi_phy_tc(jg)%ev_gwd,   datetime_old)
       !
       CALL message_forcing_action('Hines gravity wave drag (gwd)', &
            &                      is_in_sd_ed_interval,is_active   )
       !
       CALL interface_echam_gwhines(is_in_sd_ed_interval,    &
            &                       is_active,               &
            &                       patch, rl_start, rl_end, &
            &                       field, tend              )
       !
    END IF

    !-------------------------------------------------------------------
    ! Sub grid scale orographic effects: blocking and orog. gravit waves
    !-------------------------------------------------------------------
    !
    IF ( mpi_phy_tc(jg)%dt_sso > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (mpi_phy_tc(jg)%sd_sso <= datetime_old) .AND. &
            &                          (mpi_phy_tc(jg)%ed_sso >  datetime_old)
       is_active = isCurrentEventActive(mpi_phy_tc(jg)%ev_sso,   datetime_old)
       !
       CALL message_forcing_action('sub grid scale orographic effects (oro)', &
            &                      is_in_sd_ed_interval,is_active             )
       !
       CALL  interface_echam_sso(is_in_sd_ed_interval,    &
            &                    is_active,               &
            &                    patch, rl_start, rl_end, &
            &                    field, tend,             &
            &                    pdtime                   )
       !
    END IF

    !-------------------------------------------------------------------
    ! Cumulus convection
    !-------------------------------------------------------------------
    !
    IF ( mpi_phy_tc(jg)%dt_cnv > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (mpi_phy_tc(jg)%sd_cnv <= datetime_old) .AND. &
            &                          (mpi_phy_tc(jg)%ed_cnv >  datetime_old)
       is_active = isCurrentEventActive(mpi_phy_tc(jg)%ev_cnv,   datetime_old)
       !
       CALL message_forcing_action('cumulus convection (cnv)',   &
            &                      is_in_sd_ed_interval,is_active)
       !
       CALL interface_echam_convection(is_in_sd_ed_interval,    &
            &                          is_active,               &
            &                          patch, rl_start, rl_end, &
            &                          field, tend,             &
            &                          ntrac, ictop,            &
            &                          pdtime                   )
       !
    ELSE
       field% rtype(:,:) = 0.0_wp ! this probably is needed as is io in the echam_condensation
    END IF

    !-------------------------------------------------------------
    ! Update provisional physics state
    field% qtrc(:,:,:,iqc) = field% qtrc(:,:,:,iqc) + tend% qtrc(:,:,:,iqc)*pdtime
    field% qtrc(:,:,:,iqi) = field% qtrc(:,:,:,iqi) + tend% qtrc(:,:,:,iqi)*pdtime

    !-------------------------------------------------------------------
    ! Cloud processes
    !-------------------------------------------------------------------
    !
    IF ( mpi_phy_tc(jg)%dt_cld > dt_zero ) THEN
       !
       is_in_sd_ed_interval =          (mpi_phy_tc(jg)%sd_cld <= datetime_old) .AND. &
            &                          (mpi_phy_tc(jg)%ed_cld >  datetime_old)
       is_active = isCurrentEventActive(mpi_phy_tc(jg)%ev_cld,   datetime_old)
       !
       CALL message_forcing_action('cloud microphysics (cld)',   &
            &                      is_in_sd_ed_interval,is_active)
       !
       CALL interface_echam_condensation(is_in_sd_ed_interval,    &
            &                            is_active,               &
            &                            patch, rl_start, rl_end, &
            &                            jks,                     &
            &                            field, tend,             &
            &                            ictop,                   &
            &                            pdtime                   )
       !
    END IF

    !-------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)

      !  total precipitation flux
      field% totprec (jcs:jce,jb)     =  field% rsfl (jcs:jce,jb) & ! rain large scale
            &                           +field% ssfl (jcs:jce,jb) & ! snow large scale
            &                           +field% rsfc (jcs:jce,jb) & ! rain convection
            &                           +field% ssfc (jcs:jce,jb)   ! snow convection

      ! vertical integral
      field% q_phy_vi(jcs:jce,jb) = SUM(field% q_phy(jcs:jce,:,jb),DIM=2)

      ! Now compute tendencies from physics alone
      tend%   ta_phy (jcs:jce,:,jb)   = tend%   ta (jcs:jce,:,jb)   - tend%   ta_phy (jcs:jce,:,jb)
      tend%   ua_phy (jcs:jce,:,jb)   = tend%   ua (jcs:jce,:,jb)   - tend%   ua_phy (jcs:jce,:,jb)
      tend%   va_phy (jcs:jce,:,jb)   = tend%   va (jcs:jce,:,jb)   - tend%   va_phy (jcs:jce,:,jb)
      tend% qtrc_phy (jcs:jce,:,jb,:) = tend% qtrc (jcs:jce,:,jb,:) - tend% qtrc_phy (jcs:jce,:,jb,:)

      tend% ta_phy (jcs:jce,:,jb) = tend% ta_phy(jcs:jce,:,jb)*field%cpair(jcs:jce,:,jb)/field%cvair(jcs:jce,:,jb)
    END DO
!$OMP END PARALLEL DO 

    ! Done. Disassociate pointers.
    NULLIFY(field,tend)

  END SUBROUTINE echam_phy_main
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  SUBROUTINE surface_fractions( jb,jcs,jce, nbdim, field)

    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block
    TYPE(t_echam_phy_field),   POINTER :: field

    REAL(wp) :: zfrw (nbdim)              !< fraction of water (without ice) in the grid point
    REAL(wp) :: zfri (nbdim)              !< fraction of ice in the grid box
    REAL(wp) :: zfrl (nbdim)              !< fraction of land in the grid box

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
      END IF
    END DO

    ! 3.4 Merge three pieces of information into one array for vdiff
    IF (ilnd.LE.nsfc_type) field%frac_tile(jcs:jce,jb,ilnd) = zfrl(jcs:jce)
    IF (iwtr.LE.nsfc_type) field%frac_tile(jcs:jce,jb,iwtr) = zfrw(jcs:jce)
    IF (iice.LE.nsfc_type) field%frac_tile(jcs:jce,jb,iice) = zfri(jcs:jce)

  END SUBROUTINE surface_fractions
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  SUBROUTINE cloud_cover(jb,jcs,jce, nbdim, field)

    INTEGER         ,INTENT(IN) :: jb             !< block index
    INTEGER         ,INTENT(IN) :: jcs, jce       !< start/end column index within this block
    INTEGER         ,INTENT(IN) :: nbdim          !< size of this block
    TYPE(t_echam_phy_field),   POINTER :: field

    INTEGER  :: itype(nbdim)              !< type of convection

    !-------------------------------------------------------------------
    ! 3.13 DIAGNOSE CURRENT CLOUD COVER
    !-------------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_cover)

    itype(jcs:jce) = NINT(field%rtype(jcs:jce,jb))
      
    CALL cover( jce, nbdim, jks,              &! in
         &         nlev, nlevp1,              &! in
         &         itype,                     &! zfrw, zfri,       &! in
         &         field%frac_tile(:,jb,iwtr),&
         &         field%frac_tile(:,jb,iice),&
         &         field% zf(:,:,jb),         &! in
         &         field% presi_old(:,:,jb),  &! in
         &         field% presm_old(:,:,jb),  &! in
         &         field%  ta(:,:,jb),        &! in    tm1
         &         field%  qtrc(:,:,jb,iqv),  &! in    qm1
         &         field%  qtrc(:,:,jb,iqi),  &! in    xim1
         &         field%  aclc(:,:,jb),      &! out   (for "radiation" and "vdiff_down")
         &         field% rintop(:,  jb)     ) ! out   (for output)

    IF (ltimer) CALL timer_stop(timer_cover)

  END SUBROUTINE cloud_cover
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
        END IF
        IF (field%presm_old(jc,jk,jb).LT.80000._wp) THEN
          zcdnc=1.e6_wp*(zn1+(zn2-zn1)*(EXP(1._wp-zprat)))
        ELSE
          zcdnc=zn2*1.e6_wp
        END IF
        field% acdnc(jc,jk,jb) = zcdnc
        !
      END DO
    END DO

  END SUBROUTINE calculate_droplet_number
  !----------------------------------------------------------------

 
  !---------------------------------------------------------------------
  SUBROUTINE calculate_cpair_cvair_qconv( jb,jcs,jce, field)
    INTEGER                 ,INTENT(in) :: jb
    INTEGER                 ,INTENT(in) :: jcs, jce
    TYPE(t_echam_phy_field) ,POINTER    :: field

    INTEGER  :: jc, jk

    DO jk = 1,nlev
       DO jc = jcs,jce
          !
          field%cpair(jc,jk,jb) = cpd+(cpv-cpd)*field%qtrc(jc,jk,jb,iqv)
          field%cvair(jc,jk,jb) = cvd+(cvv-cvd)*field%qtrc(jc,jk,jb,iqv)
          !
          field%qconv(jc,jk,jb) = 1._wp/(field%mair(jc,jk,jb)*field%cpair(jc,jk,jb))
          !
       END DO
    END DO
    
  END SUBROUTINE calculate_cpair_cvair_qconv
  !----------------------------------------------------------------

 
  !---------------------------------------------------------------------
  SUBROUTINE message_forcing_action(process,is_in_sd_ed_interval,is_active)
    CHARACTER(LEN=*) ,INTENT(in) :: process
    LOGICAL          ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL          ,INTENT(in) :: is_active

    IF (is_in_sd_ed_interval) THEN
       IF (is_active) THEN
          CALL message('echam_phy_main','compute forcing by '//process)
       ELSE
          CALL message('echam_phy_main','recycle forcing by '//process)
       END IF
    ELSE
       CALL    message('echam_phy_main','no      forcing by '//process)
    END IF

  END SUBROUTINE message_forcing_action
  !---------------------------------------------------------------------



END MODULE mo_echam_phy_main
