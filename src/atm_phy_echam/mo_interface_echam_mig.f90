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
  USE mo_echam_cov_config    ,ONLY: echam_cov_config
  USE mo_echam_phy_memory    ,ONLY: t_echam_phy_field, prm_field, &
       &                            t_echam_phy_tend,  prm_tend

  USE mo_timer               ,ONLY: ltimer, timer_start, timer_stop, timer_mig

  USE mo_physical_constants  ,ONLY: cvd
  USE mo_run_config          ,ONLY: iqv, iqc, iqi, iqr, iqs, iqg, msg_level
  USE gscp_graupel           ,ONLY: graupel
  USE gscp_data              ,ONLY: cloud_num
  USE mo_satad               ,ONLY: satad_v_3d, satad_v_3D_gpu        ! new saturation adjustment
  USE mo_echam_mig_config    ,ONLY: echam_mig_config
  USE mo_fortran_tools       ,ONLY: init
  USE mo_exception           ,ONLY: warning

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
    LOGICAL                             :: lparamcpl
    INTEGER                             :: fc_mig, jkscov
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    !
    REAL(wp) :: zqrsflux(nproma,nlev)
    REAL(wp) :: zqnc(nproma)
    !
    REAL(wp) :: tend_ta_mig, tend_qtrc_mig_iqv, tend_qtrc_mig_iqc, tend_qtrc_mig_iqi, &
                tend_qtrc_mig_iqr, tend_qtrc_mig_iqs, tend_qtrc_mig_iqg
    LOGICAL  :: associated_ta_mig, associated_qtrc_mig
    !
    ! Local variables for security
    !
    REAL(wp) :: xlta(nproma,nlev)
    REAL(wp) :: xlqv(nproma,nlev), xlqc(nproma,nlev)
    REAL(wp) :: xlqi(nproma,nlev), xlqr(nproma,nlev)
    REAL(wp) :: xlqs(nproma,nlev), xlqg(nproma,nlev)

    REAL(wp) :: zdtr      ! reciprocal of timestep

    INTEGER :: jk, jl

    IF (ltimer) call timer_start(timer_mig)

    ! associate pointers
    lparamcpl =  echam_phy_config(jg)%lparamcpl
    fc_mig    =  echam_phy_config(jg)%fc_mig

    jkscov    =  echam_cov_config(jg)% jkscov

    field     => prm_field(jg)
    tend      => prm_tend (jg)

    associated_ta_mig   = ASSOCIATED(tend% ta_mig)
    associated_qtrc_mig = ASSOCIATED(tend% qtrc_mig)

    !$ACC DATA PCREATE( zqrsflux, zqnc           , &
    !$ACC               xlta                     , &
    !$ACC               xlqv, xlqc               , &
    !$ACC               xlqi, xlqr               , &
    !$ACC               xlqs, xlqg               ) &
    !$ACC      PRESENT( field, tend              )

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(NONE) ASYNC(1)
    DO jl = 1, nproma
      zqnc(jl) = cloud_num
    END DO

    ! reciprocal of timestep
    zdtr = 1.0_wp/pdtime

    IF ( is_in_sd_ed_interval ) THEN
      !
      IF ( is_active ) THEN
        !
        ! Copy fields to local variables
        ! This prevents the fields to be updated directly in graupel scheme
        !  and is used for the calculation of tendencies at the end of the complete
        !  satad - graupel - satad -cycle
        !
        !$ACC PARALLEL LOOP GANG VECTOR TILE(128,1) DEFAULT(NONE) ASYNC(1)
        DO jk = 1,nlev
          DO jl = jcs,jce
            xlta(jl,jk) = field% ta   (jl,jk,jb)
            xlqv(jl,jk) = field% qtrc (jl,jk,jb,iqv)
            xlqc(jl,jk) = field% qtrc (jl,jk,jb,iqc)
            xlqi(jl,jk) = field% qtrc (jl,jk,jb,iqi)
            xlqr(jl,jk) = field% qtrc (jl,jk,jb,iqr)
            xlqs(jl,jk) = field% qtrc (jl,jk,jb,iqs)
            xlqg(jl,jk) = field% qtrc (jl,jk,jb,iqg)
          END DO
        END DO
        !$ACC END PARALLEL

      !!-------------------------------------------------------------------------
      !> Initial saturation adjustment (a second one follows at the end of the microphysics)
      !!-------------------------------------------------------------------------
#ifdef _OPENACC
        !$ACC WAIT
        CALL satad_v_3d_gpu(                               &
#else
        CALL satad_v_3d(                                   &
#endif
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
              & klo      = jkscov                         ,& !> IN
              & kup      = nlev                            & !> IN
              )
        !
        CALL graupel (                                &
              & nvec   =nproma                      , & !> in:  actual array size
              & ke     =nlev                        , & !< in:  actual array size
              & ivstart=jcs                         , & !< in:  start index of calculation
              & ivend  =jce                         , & !< in:  end index of calculation
              & kstart =jkscov                      , & !< in:  vertical start index ! optional
              & idbg   =msg_level                   , &
              & zdt    =pdtime                      , & !< in:  timestep
              & qi0    =echam_mig_config(jg)%qi0    , & !< in: cloud ice threshold for autoconversion
              & qc0    =echam_mig_config(jg)%qc0    , & !< in: cloud water threshold for autoconversion
              & dz     =field% dz       (:,:,jb)    , & !< in: vertical layer thickness
              & t      =xlta            (:,:)       , & !< inout: temp
              & p      =field% pfull    (:,:,jb)    , & !< in:  pressure
              & rho    =field% rho      (:,:,jb)    , & !< in:  density
              & qv     =xlqv            (:,:)       , & !< inout:sp humidity
              & qc     =xlqc            (:,:)       , & !< inout:cloud water
              & qi     =xlqi            (:,:)       , & !< inout: ice
              & qr     =xlqr            (:,:)       , & !< inout:rain
              & qs     =xlqs            (:,:)       , & !< inout: snow
              & qg     =xlqg            (:,:)       , & !< inout: graupel
              & qnc    =zqnc                        , & !< inout: cloud droplet number
              & prr_gsp=field%rain_gsp_rate (:, jb) , & !< inout precp rate rain
              & prs_gsp=field%snow_gsp_rate (:, jb) , & !< inout precp rate snow
              & prg_gsp=field%graupel_gsp_rate (:, jb) , & !< inout precp rate graupel
              & qrsflux=zqrsflux                    , & !< out  precipitation flux
              & l_cv=.TRUE.                         , &
              & ldiag_ttend = echam_mig_config(jg)%ldiag_ttend , & !< in:  if temp.  tendency shall be diagnosed
              & ldiag_qtend = echam_mig_config(jg)%ldiag_qtend )   !< in:  if moisture tendencies shall be diagnosed
                                                                   ! IF true tendencies have to be re-implemented

        !!-------------------------------------------------------------------------
        !> Final saturation adjustment (as this has been removed from the end of the microphysics)
        !!-------------------------------------------------------------------------
#ifdef _OPENACC
        CALL satad_v_3d_gpu(                                &
#else
        CALL satad_v_3d(                                    &
#endif
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
              & klo      = jkscov                         ,& !> IN
              & kup      = nlev                            & !> IN
              )

        !
        ! Calculate rain and snow
        !
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jl = jcs, jce
          field% rsfl(jl,jb) = field%rain_gsp_rate (jl,jb)
          field% ssfl(jl,jb) = field%snow_gsp_rate (jl,jb) + field%graupel_gsp_rate (jl,jb)
        END DO
        !$ACC END PARALLEL
        !
        !
        ! is_active true: calculate new tendencies, update memory if required, update model state
        ! cases by is_active are deliberately separated to avoid temporary memory
        !
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR TILE(128,1)
        DO jk = jkscov, nlev
          DO jl = jcs, jce
            !
            ! Calculate tendencies and convert temperature tendency, as computed for constant volume in satad/graupel
            !  to constant pressure
            !
            tend_ta_mig = (xlta(jl,jk)-field% ta(jl,jk,jb))*zdtr                &
                                     & * field% cvair(jl,jk,jb)/field% cpair(jl,jk,jb)
            tend_qtrc_mig_iqv = MAX(-field% qtrc(jl,jk,jb,iqv)*zdtr,            &
                                        & (xlqv(jl,jk)-field% qtrc(jl,jk,jb,iqv))*zdtr)
            tend_qtrc_mig_iqc = MAX(-field% qtrc(jl,jk,jb,iqc)*zdtr,            &
                                        & (xlqc(jl,jk)-field% qtrc(jl,jk,jb,iqc))*zdtr)
            tend_qtrc_mig_iqi = MAX(-field% qtrc(jl,jk,jb,iqi)*zdtr,            &
                                        & (xlqi(jl,jk)-field% qtrc(jl,jk,jb,iqi))*zdtr)
            tend_qtrc_mig_iqr = MAX(-field% qtrc(jl,jk,jb,iqr)*zdtr,            &
                                        & (xlqr(jl,jk)-field% qtrc(jl,jk,jb,iqr))*zdtr)
            tend_qtrc_mig_iqs = MAX(-field% qtrc(jl,jk,jb,iqs)*zdtr,            &
                                        & (xlqs(jl,jk)-field% qtrc(jl,jk,jb,iqs))*zdtr)
            tend_qtrc_mig_iqg = MAX(-field% qtrc(jl,jk,jb,iqg)*zdtr,            &
                                        & (xlqg(jl,jk)-field% qtrc(jl,jk,jb,iqg))*zdtr)
            !
            ! store in memory for output or recycling
            !
            IF (associated_ta_mig) tend% ta_mig(jl,jk,jb) = tend_ta_mig
            !
            IF (associated_qtrc_mig) THEN
              tend% qtrc_mig(jl,jk,jb,iqv) = tend_qtrc_mig_iqv
              tend% qtrc_mig(jl,jk,jb,iqc) = tend_qtrc_mig_iqc
              tend% qtrc_mig(jl,jk,jb,iqi) = tend_qtrc_mig_iqi
              tend% qtrc_mig(jl,jk,jb,iqr) = tend_qtrc_mig_iqr
              tend% qtrc_mig(jl,jk,jb,iqs) = tend_qtrc_mig_iqs
              tend% qtrc_mig(jl,jk,jb,iqg) = tend_qtrc_mig_iqg
            ENDIF
            !
            ! Accumulate tendencies for later updating the model state
            !
            SELECT CASE(fc_mig)
            CASE(0)
            ! diagnostic, do not use tendency
            CASE(1)
            ! use tendency to update the model state
              tend%   ta_phy(jl,jk,jb)      = tend%   ta_phy(jl,jk,jb)     + tend_ta_mig  
              tend% qtrc_phy(jl,jk,jb,iqv)  = tend% qtrc_phy(jl,jk,jb,iqv) + tend_qtrc_mig_iqv
              tend% qtrc_phy(jl,jk,jb,iqc)  = tend% qtrc_phy(jl,jk,jb,iqc) + tend_qtrc_mig_iqc
              tend% qtrc_phy(jl,jk,jb,iqi)  = tend% qtrc_phy(jl,jk,jb,iqi) + tend_qtrc_mig_iqi
              tend% qtrc_phy(jl,jk,jb,iqr)  = tend% qtrc_phy(jl,jk,jb,iqr) + tend_qtrc_mig_iqr
              tend% qtrc_phy(jl,jk,jb,iqs)  = tend% qtrc_phy(jl,jk,jb,iqs) + tend_qtrc_mig_iqs
              tend% qtrc_phy(jl,jk,jb,iqg)  = tend% qtrc_phy(jl,jk,jb,iqg) + tend_qtrc_mig_iqg
            END SELECT
            !
            ! update physics state for input to the next physics process
            IF (lparamcpl) THEN
              SELECT CASE(fc_mig)
              CASE(0)
                  ! diagnostic, do not use tendency
              CASE(1)
                field%   ta(jl,jk,jb)      = field%   ta(jl,jk,jb)      + tend_ta_mig      *pdtime
                field% qtrc(jl,jk,jb,iqv)  = field% qtrc(jl,jk,jb,iqv)  + tend_qtrc_mig_iqv*pdtime
                field% qtrc(jl,jk,jb,iqc)  = field% qtrc(jl,jk,jb,iqc)  + tend_qtrc_mig_iqc*pdtime
                field% qtrc(jl,jk,jb,iqi)  = field% qtrc(jl,jk,jb,iqi)  + tend_qtrc_mig_iqi*pdtime
                field% qtrc(jl,jk,jb,iqr)  = field% qtrc(jl,jk,jb,iqr)  + tend_qtrc_mig_iqr*pdtime
                field% qtrc(jl,jk,jb,iqs)  = field% qtrc(jl,jk,jb,iqs)  + tend_qtrc_mig_iqs*pdtime
                field% qtrc(jl,jk,jb,iqg)  = field% qtrc(jl,jk,jb,iqg)  + tend_qtrc_mig_iqg*pdtime
              END SELECT
            END IF
          ENDDO
        ENDDO
        !$ACC END PARALLEL
        !
      ELSE    ! is_active
        !
        ! is_active false: reuse tendencies from memory if available, update model state
        !
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
        !$ACC LOOP GANG VECTOR TILE(128,1)
        DO jk = jkscov, nlev
          DO jl = jcs, jce
            !
            ! retrieve from memory for recycling
            !
            tend_ta_mig       = 0.0_wp
            tend_qtrc_mig_iqv = 0.0_wp
            tend_qtrc_mig_iqc = 0.0_wp
            tend_qtrc_mig_iqi = 0.0_wp
            tend_qtrc_mig_iqr = 0.0_wp
            tend_qtrc_mig_iqs = 0.0_wp
            tend_qtrc_mig_iqg = 0.0_wp
            !
            IF (associated_ta_mig) tend_ta_mig = tend% ta_mig(jl,jk,jb)
            !
            IF (associated_qtrc_mig) THEN
              tend_qtrc_mig_iqv = tend% qtrc_mig(jl,jk,jb,iqv)
              tend_qtrc_mig_iqc = tend% qtrc_mig(jl,jk,jb,iqc)
              tend_qtrc_mig_iqi = tend% qtrc_mig(jl,jk,jb,iqi)
              tend_qtrc_mig_iqr = tend% qtrc_mig(jl,jk,jb,iqr)
              tend_qtrc_mig_iqs = tend% qtrc_mig(jl,jk,jb,iqs)
              tend_qtrc_mig_iqg = tend% qtrc_mig(jl,jk,jb,iqg)
            ENDIF
            !
            ! Accumulate tendencies for later updating the model state
            !
            SELECT CASE(fc_mig)
            CASE(0)
            ! diagnostic, do not use tendency
            CASE(1)
            ! use tendency to update the model state
              tend%   ta_phy(jl,jk,jb)      = tend%   ta_phy(jl,jk,jb)     + tend_ta_mig  
              tend% qtrc_phy(jl,jk,jb,iqv)  = tend% qtrc_phy(jl,jk,jb,iqv) + tend_qtrc_mig_iqv
              tend% qtrc_phy(jl,jk,jb,iqc)  = tend% qtrc_phy(jl,jk,jb,iqc) + tend_qtrc_mig_iqc
              tend% qtrc_phy(jl,jk,jb,iqi)  = tend% qtrc_phy(jl,jk,jb,iqi) + tend_qtrc_mig_iqi
              tend% qtrc_phy(jl,jk,jb,iqr)  = tend% qtrc_phy(jl,jk,jb,iqr) + tend_qtrc_mig_iqr
              tend% qtrc_phy(jl,jk,jb,iqs)  = tend% qtrc_phy(jl,jk,jb,iqs) + tend_qtrc_mig_iqs
              tend% qtrc_phy(jl,jk,jb,iqg)  = tend% qtrc_phy(jl,jk,jb,iqg) + tend_qtrc_mig_iqg
            END SELECT
            !
            ! update physics state for input to the next physics process
            IF (lparamcpl) THEN
              SELECT CASE(fc_mig)
              CASE(0)
                  ! diagnostic, do not use tendency
              CASE(1)
                field%   ta(jl,jk,jb)      = field%   ta(jl,jk,jb)      + tend_ta_mig      *pdtime
                field% qtrc(jl,jk,jb,iqv)  = field% qtrc(jl,jk,jb,iqv)  + tend_qtrc_mig_iqv*pdtime
                field% qtrc(jl,jk,jb,iqc)  = field% qtrc(jl,jk,jb,iqc)  + tend_qtrc_mig_iqc*pdtime
                field% qtrc(jl,jk,jb,iqi)  = field% qtrc(jl,jk,jb,iqi)  + tend_qtrc_mig_iqi*pdtime
                field% qtrc(jl,jk,jb,iqr)  = field% qtrc(jl,jk,jb,iqr)  + tend_qtrc_mig_iqr*pdtime
                field% qtrc(jl,jk,jb,iqs)  = field% qtrc(jl,jk,jb,iqs)  + tend_qtrc_mig_iqs*pdtime
                field% qtrc(jl,jk,jb,iqg)  = field% qtrc(jl,jk,jb,iqg)  + tend_qtrc_mig_iqg*pdtime
              END SELECT
            END IF
          ENDDO
        ENDDO
        !$ACC END PARALLEL
        !
      END IF
      !
    ELSE       ! is_in_sd_ed_interval
       !
       !$ACC PARALLEL DEFAULT(NONE) ASYNC(1)
       IF (associated_ta_mig) THEN 
         !$ACC LOOP GANG(static:1) VECTOR TILE(128,1)
         DO jk = 1, nlev
           DO jl = jcs, jce
             tend% ta_mig(jl,jk,jb) = 0.0_wp
           END DO
         END DO
       END IF
       !
       IF (associated_qtrc_mig) THEN
          !$ACC LOOP GANG(static:1) VECTOR TILE(128,1)
          DO jk = 1, nlev
            DO jl = jcs, jce
              tend% qtrc_mig(jl,jk,jb,iqv) = 0.0_wp
              tend% qtrc_mig(jl,jk,jb,iqc) = 0.0_wp
              tend% qtrc_mig(jl,jk,jb,iqi) = 0.0_wp
              tend% qtrc_mig(jl,jk,jb,iqr) = 0.0_wp
              tend% qtrc_mig(jl,jk,jb,iqs) = 0.0_wp
              tend% qtrc_mig(jl,jk,jb,iqg) = 0.0_wp
            END DO
          END DO
       END IF
       !$ACC END PARALLEL
       !
    END IF     ! is_in_sd_ed_interval

    !$ACC WAIT
    !$ACC END DATA

    ! disassociate pointers

    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) call timer_stop(timer_mig)

  END SUBROUTINE interface_echam_mig

END MODULE mo_interface_echam_mig
