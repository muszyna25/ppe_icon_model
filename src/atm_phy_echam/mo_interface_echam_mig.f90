!>
!! @brief Subroutine interface_echam_mig calls NWP graupel scheme
!!
!! @author Monika Esch, MPI-M
!!
!! @par Revision History
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_interface_echam_mig

  USE mo_kind                ,ONLY: wp
  USE mtime                  ,ONLY: datetime

  USE mo_echam_phy_config    ,ONLY: echam_phy_config
  USE mo_echam_cld_config    ,ONLY: echam_cld_config
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
       &                            t_echam_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_mig

  USE mo_run_config          ,ONLY: iqv, iqc, iqi, iqr, iqs, iqg, msg_level
  USE gscp_graupel           ,ONLY: graupel
  USE gscp_data              ,ONLY: cloud_num
  USE mo_satad               ,ONLY: satad_v_3d        ! new saturation adjustment
  USE mo_echam_mig_config    ,ONLY: echam_mig_config
  !$ser verbatim USE mo_ser_echam_mig, ONLY: serialize_mig_input,&
  !$ser verbatim                             serialize_mig_output,&
  !$ser verbatim                             serialize_mig_after_satad1,&
  !$ser verbatim                             serialize_mig_before_satad2

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: interface_echam_mig

CONTAINS

  SUBROUTINE interface_echam_mig(jg,jb,jcs,jce        ,&
       &                         nproma,nlev,ntracer  ,& 
       &                         is_in_sd_ed_interval ,&
       &                         is_active            ,&
       &                         datetime_old         ,&
       &                         pdtime               )

    ! Arguments
    !
    INTEGER                 ,INTENT(in) :: jg,jb,jcs,jce
    INTEGER                 ,INTENT(in) :: nproma,nlev,ntracer
    LOGICAL                 ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in) :: is_active
    TYPE(datetime)          ,POINTER    :: datetime_old
    REAL(wp)                ,INTENT(in) :: pdtime

    ! Pointers
    !
    LOGICAL                 ,POINTER    :: lparamcpl
    INTEGER                 ,POINTER    :: fc_mig, jks
    !REAL                    ,POINTER    :: qi0_nwp, qc0_nwp
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    !
    REAL(wp) :: zprec_r(nproma), zprec_i(nproma), zprec_s(nproma)
    REAL(wp) :: zprec_g(nproma), zqrsflux(nproma,nlev)
    REAL(wp) :: zqnc(nproma)
    !
    REAL(wp)                            :: tend_ta_mig  (nproma,nlev)
    REAL(wp)                            :: tend_qtrc_mig(nproma,nlev,ntracer)
    !
    ! Local variables for security
    !
    REAL(wp) :: xlta(nproma,nlev)
    REAL(wp) :: xlqv(nproma,nlev), xlqc(nproma,nlev)
    REAL(wp) :: xlqi(nproma,nlev), xlqr(nproma,nlev)
    REAL(wp) :: xlqs(nproma,nlev), xlqg(nproma,nlev)

    REAL(wp) :: zdtr      ! reciprocal of timestep

    INTEGER :: jk

    IF (ltimer) call timer_start(timer_mig)

    zqnc     = cloud_num

    ! associate pointers
    lparamcpl => echam_phy_config(jg)%lparamcpl
    fc_mig    => echam_phy_config(jg)%fc_mig

    jks       => echam_cld_config(jg)% jks

    !qi0_nwp   => echam_mig_config(jg)%qi0_nwp
    !qc0_nwp   => echam_mig_config(jg)%qc0_nwp

    field     => prm_field(jg)
    tend      => prm_tend (jg)

    !$ser verbatim call serialize_mig_input(jg, jb, jcs, jce, nproma, nlev, field)

    ! preset local mig tendencies to avoid NaNs

     tend_ta_mig(:,:)     = 0.0_wp
     tend_qtrc_mig(:,:,:) = 0.0_wp

    ! reciprocal of timestep

     zdtr = 1.0_wp/pdtime

    IF ( is_in_sd_ed_interval ) THEN
       !
       IF ( is_active ) THEN

       ! Copy fields to local variables
       ! This prevents the fields to be updated directly in graupel scheme
       !  and is used for the calculation of tendencies at the end of the complete
       !  satad - graupel - satad -cycle
       !

         xlta     = field% ta       (:,:,jb)
         xlqv     = field% qtrc     (:,:,jb,iqv)
         xlqc     = field% qtrc     (:,:,jb,iqc)
         xlqi     = field% qtrc     (:,:,jb,iqi)
         xlqr     = field% qtrc     (:,:,jb,iqr)
         xlqs     = field% qtrc     (:,:,jb,iqs)
         xlqg     = field% qtrc     (:,:,jb,iqg)

      !!-------------------------------------------------------------------------
      !> Initial saturation adjustment (a second one follows at the end of the microphysics)
      !!-------------------------------------------------------------------------
        CALL satad_v_3D( &
               & maxiter  = 10                             ,& !> IN
               & tol      = 1.e-3_wp                       ,& !> IN
               & te       = xlta               (:,:)       ,& !> INOUT
               & qve      = xlqv               (:,:)       ,& !> INOUT
               & qce      = xlqc               (:,:)       ,& !> INOUT
               & rhotot   = field% rho         (:,:,jb)    ,& !> IN
               & idim     = nproma                         ,& !> IN
               & kdim     = nlev                           ,& !> IN
               & ilo      = jcs                            ,& !> IN
               & iup      = jce                            ,& !> IN
               & klo      = jks                            ,& !> IN
               & kup      = nlev                            & !> IN
               )
!
    !$ser verbatim call serialize_mig_after_satad1(jg, jb, jcs, jce, nproma, nlev, field,&
    !$ser verbatim        xlta, xlqv, xlqc, xlqi, xlqr, xlqs, xlqg, zqnc,&
    !$ser verbatim        zprec_r, zprec_s, zprec_g, zqrsflux)
          !
            CALL graupel (                            &
              & nvec   =nproma                      , & !> in:  actual array size
              & ke     =nlev                        , & !< in:  actual array size
              & ivstart=jcs                         , & !< in:  start index of calculation
              & ivend  =jce                         , & !< in:  end index of calculation
              & kstart =jks                         , & !< in:  vertical start index ! optional
              !& idbg   = 51                         , &
              & zdt    =pdtime                      , & !< in:  timestep
              & qi0    =echam_mig_config(jg)%qi0_nwp, & !< in: cloud ice threshold for autoconversion
              & qc0    =echam_mig_config(jg)%qc0_nwp, & !< in: cloud water threshold for autoconversion
              & dz     =field% dz       (:,:,jb)    , & !< in: vertical layer thickness
              & t      =xlta            (:,:)       , & !< inout: temp
              & p      =field% presm_old(:,:,jb)    , & !< in:  pressure
              & rho    =field% rho      (:,:,jb)    , & !< in:  density
              & qv     =xlqv            (:,:)       , & !< inout:sp humidity
              & qc     =xlqc            (:,:)       , & !< inout:cloud water
              & qi     =xlqi            (:,:)       , & !< inout: ice
              & qr     =xlqr            (:,:)       , & !< inout:rain
              & qs     =xlqs            (:,:)       , & !< inout: snow
              & qg     =xlqg            (:,:)       , & !< inout: graupel
              & qnc    =zqnc                        , & !< inout: cloud droplet number
              & prr_gsp=zprec_r                     , & !< inout precp rate rain
              & prs_gsp=zprec_s                     , & !< inout precp rate snow
              & prg_gsp=zprec_g                     , & !< inout precp rate graupel
              & qrsflux=zqrsflux                    , & !< out  precipitation flux
              & l_cv=.TRUE.                          )

    !$ser verbatim call serialize_mig_before_satad2(jg, jb, jcs, jce, nproma, nlev, field,&
    !$ser verbatim      & xlta, xlqv, xlqc, xlqi, xlqr, xlqs, xlqg, zqnc,&
    !$ser verbatim      & zprec_r, zprec_s, zprec_g, zqrsflux)

      !!-------------------------------------------------------------------------
      !> Final saturation adjustment (as this has been removed from the end of the microphysics)
      !!-------------------------------------------------------------------------
        CALL satad_v_3D( &
               & maxiter  = 10                             ,& !> IN
               & tol      = 1.e-3_wp                       ,& !> IN
               & te       = xlta (:,:)                     ,& !> INOUT
               & qve      = xlqv (:,:)                     ,& !> INOUT
               & qce      = xlqc (:,:)                     ,& !> INOUT
               & rhotot   = field% rho         (:,:,jb)    ,& !> IN
               & idim     = nproma                         ,& !> IN
               & kdim     = nlev                           ,& !> IN
               & ilo      = jcs                            ,& !> IN
               & iup      = jce                            ,& !> IN
               & klo      = jks                            ,& !> IN
               & kup      = nlev                            & !> IN
               )
      !
      ! Calculate rain and snow
      !
              field% rsfl   (jcs:jce,jb)      = zprec_r (jcs:jce)
              field% ssfl   (jcs:jce,jb)      = zprec_s (jcs:jce) + zprec_g (jcs:jce)
      !
      ! Calculate tendencies
      !
            DO jk = jks, nlev
             tend_ta_mig(jcs:jce,jk) = (xlta(jcs:jce,jk)-field% ta(jcs:jce,jk,jb))*zdtr
             tend_qtrc_mig(jcs:jce,jk,iqv) = MAX(-field% qtrc(jcs:jce,jk,jb,iqv)*zdtr,            &
                                         & (xlqv(jcs:jce,jk)-field% qtrc(jcs:jce,jk,jb,iqv))*zdtr)
             tend_qtrc_mig(jcs:jce,jk,iqc) = MAX(-field% qtrc(jcs:jce,jk,jb,iqc)*zdtr,            &
                                         & (xlqc(jcs:jce,jk)-field% qtrc(jcs:jce,jk,jb,iqc))*zdtr)
             tend_qtrc_mig(jcs:jce,jk,iqi) = MAX(-field% qtrc(jcs:jce,jk,jb,iqi)*zdtr,            &
                                         & (xlqi(jcs:jce,jk)-field% qtrc(jcs:jce,jk,jb,iqi))*zdtr)
             tend_qtrc_mig(jcs:jce,jk,iqr) = MAX(-field% qtrc(jcs:jce,jk,jb,iqr)*zdtr,            &
                                         & (xlqr(jcs:jce,jk)-field% qtrc(jcs:jce,jk,jb,iqr))*zdtr)
             tend_qtrc_mig(jcs:jce,jk,iqs) = MAX(-field% qtrc(jcs:jce,jk,jb,iqs)*zdtr,            &
                                         & (xlqs(jcs:jce,jk)-field% qtrc(jcs:jce,jk,jb,iqs))*zdtr)
             tend_qtrc_mig(jcs:jce,jk,iqg) = MAX(-field% qtrc(jcs:jce,jk,jb,iqg)*zdtr,            &
                                         & (xlqg(jcs:jce,jk)-field% qtrc(jcs:jce,jk,jb,iqg))*zdtr)
            ENDDO
          !
          ! store in memory for output or recycling
          !
          IF (ASSOCIATED(tend% ta_mig)) tend% ta_mig(jcs:jce,:,jb) = tend_ta_mig(jcs:jce,:)
          !
          IF (ASSOCIATED(tend% qtrc_mig )) THEN
             tend% qtrc_mig(jcs:jce,:,jb,iqv) = tend_qtrc_mig(jcs:jce,:,iqv)
             tend% qtrc_mig(jcs:jce,:,jb,iqc) = tend_qtrc_mig(jcs:jce,:,iqc)
             tend% qtrc_mig(jcs:jce,:,jb,iqi) = tend_qtrc_mig(jcs:jce,:,iqi)
             tend% qtrc_mig(jcs:jce,:,jb,iqr) = tend_qtrc_mig(jcs:jce,:,iqr)
             tend% qtrc_mig(jcs:jce,:,jb,iqs) = tend_qtrc_mig(jcs:jce,:,iqs)
             tend% qtrc_mig(jcs:jce,:,jb,iqg) = tend_qtrc_mig(jcs:jce,:,iqg)
          END IF
          !
       ELSE
          !
          ! retrieve from memory for recycling
          !
          IF (ASSOCIATED(tend% ta_mig)) tend_ta_mig(jcs:jce,:) = tend% ta_mig(jcs:jce,:,jb)
          !
          IF (ASSOCIATED(tend% qtrc_mig )) THEN
             tend_qtrc_mig(jcs:jce,:,iqv) = tend% qtrc_mig(jcs:jce,:,jb,iqv)
             tend_qtrc_mig(jcs:jce,:,iqc) = tend% qtrc_mig(jcs:jce,:,jb,iqc)
             tend_qtrc_mig(jcs:jce,:,iqi) = tend% qtrc_mig(jcs:jce,:,jb,iqi)
             tend_qtrc_mig(jcs:jce,:,iqr) = tend% qtrc_mig(jcs:jce,:,jb,iqr)
             tend_qtrc_mig(jcs:jce,:,iqs) = tend% qtrc_mig(jcs:jce,:,jb,iqs)
             tend_qtrc_mig(jcs:jce,:,iqg) = tend% qtrc_mig(jcs:jce,:,jb,iqg)
          END IF
          !

       END IF
       !

       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_mig)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          ! use tendency to update the model state
          IF (ASSOCIATED(tend% ta_mig))           &
          & tend%   ta_phy(jcs:jce,:,jb)      = tend%   ta_phy(jcs:jce,:,jb)     + tend_ta_mig  (jcs:jce,:)
          IF (ASSOCIATED(tend% qtrc_mig )) THEN
            tend% qtrc_phy(jcs:jce,:,jb,iqv)  = tend% qtrc_phy(jcs:jce,:,jb,iqv) + tend_qtrc_mig(jcs:jce,:,iqv)
            tend% qtrc_phy(jcs:jce,:,jb,iqc)  = tend% qtrc_phy(jcs:jce,:,jb,iqc) + tend_qtrc_mig(jcs:jce,:,iqc)
            tend% qtrc_phy(jcs:jce,:,jb,iqi)  = tend% qtrc_phy(jcs:jce,:,jb,iqi) + tend_qtrc_mig(jcs:jce,:,iqi)
            tend% qtrc_phy(jcs:jce,:,jb,iqr)  = tend% qtrc_phy(jcs:jce,:,jb,iqr) + tend_qtrc_mig(jcs:jce,:,iqr)
            tend% qtrc_phy(jcs:jce,:,jb,iqs)  = tend% qtrc_phy(jcs:jce,:,jb,iqs) + tend_qtrc_mig(jcs:jce,:,iqs)
            tend% qtrc_phy(jcs:jce,:,jb,iqg)  = tend% qtrc_phy(jcs:jce,:,jb,iqg) + tend_qtrc_mig(jcs:jce,:,iqg)
          END IF
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
       END SELECT
       !
       IF (lparamcpl) THEN
         SELECT CASE(fc_mig)
         CASE(0)
            ! diagnostic, do not use tendency
         CASE(1)
            ! use tendency to update physics state for input to the next physics process
            !
            IF (ASSOCIATED(tend% ta_mig)) &
            & field%   ta(jcs:jce,:,jb)      = field%   ta(jcs:jce,:,jb)      + tend%   ta_mig(jcs:jce,:,jb)    *pdtime
            IF (ASSOCIATED(tend% qtrc_mig )) THEN
              field% qtrc(jcs:jce,:,jb,iqv)  = field% qtrc(jcs:jce,:,jb,iqv)  + tend% qtrc_mig(jcs:jce,:,jb,iqv)*pdtime
              field% qtrc(jcs:jce,:,jb,iqc)  = field% qtrc(jcs:jce,:,jb,iqc)  + tend% qtrc_mig(jcs:jce,:,jb,iqc)*pdtime
              field% qtrc(jcs:jce,:,jb,iqi)  = field% qtrc(jcs:jce,:,jb,iqi)  + tend% qtrc_mig(jcs:jce,:,jb,iqi)*pdtime
              field% qtrc(jcs:jce,:,jb,iqr)  = field% qtrc(jcs:jce,:,jb,iqr)  + tend% qtrc_mig(jcs:jce,:,jb,iqr)*pdtime
              field% qtrc(jcs:jce,:,jb,iqs)  = field% qtrc(jcs:jce,:,jb,iqs)  + tend% qtrc_mig(jcs:jce,:,jb,iqs)*pdtime
              field% qtrc(jcs:jce,:,jb,iqg)  = field% qtrc(jcs:jce,:,jb,iqg)  + tend% qtrc_mig(jcs:jce,:,jb,iqg)*pdtime
            END IF
!!$         CASE(2)
!!$            ! use tendency as forcing in the dynamics
!!$            ...
         END SELECT
       END IF
       !
    ELSE
       !
       IF (ASSOCIATED(tend% ta_mig)) tend% ta_mig(jcs:jce,:,jb) = 0.0_wp
       !
       IF (ASSOCIATED(tend% qtrc_mig )) THEN
          tend% qtrc_mig(jcs:jce,:,jb,iqv) = 0.0_wp
          tend% qtrc_mig(jcs:jce,:,jb,iqc) = 0.0_wp
          tend% qtrc_mig(jcs:jce,:,jb,iqi) = 0.0_wp
          tend% qtrc_mig(jcs:jce,:,jb,iqr) = 0.0_wp
          tend% qtrc_mig(jcs:jce,:,jb,iqs) = 0.0_wp
          tend% qtrc_mig(jcs:jce,:,jb,iqg) = 0.0_wp
       END IF
       !
    END IF

    !$ser verbatim call serialize_mig_input(jg, jb, jcs, jce, nproma, nlev, field)
    ! disassociate pointers
    NULLIFY(lparamcpl)
    NULLIFY(fc_mig)
    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) call timer_stop(timer_mig)

  END SUBROUTINE interface_echam_mig

END MODULE mo_interface_echam_mig
