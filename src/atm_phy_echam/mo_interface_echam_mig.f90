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
  USE mo_satad               ,ONLY: satad_v_3d        ! new saturation adjustment
  USE mo_echam_mig_config    ,ONLY: echam_mig_config
  USE mo_fortran_tools       ,ONLY: init
  USE mo_exception           ,ONLY: warning

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
    LOGICAL                             :: lparamcpl
    INTEGER                             :: fc_mig, jkscov
    TYPE(t_echam_phy_field) ,POINTER    :: field
    TYPE(t_echam_phy_tend)  ,POINTER    :: tend

    ! Local variables
    !
    REAL(wp) :: zqrsflux(nproma,nlev)
    REAL(wp) :: zqnc(nproma)
    !
    REAL(wp)                            :: tend_ta_mig  (nproma,nlev)
    REAL(wp)                            :: tend_qtrc_mig(nproma,nlev,ntracer)
    REAL(wp)                            :: ddt_tend_t   (nproma,nlev)
    REAL(wp)                            :: ddt_tend_qv  (nproma,nlev)
    REAL(wp)                            :: ddt_tend_qc  (nproma,nlev)
    REAL(wp)                            :: ddt_tend_qi  (nproma,nlev)
    REAL(wp)                            :: ddt_tend_qr  (nproma,nlev)
    REAL(wp)                            :: ddt_tend_qs  (nproma,nlev)
    REAL(wp)                            :: ddt_tend_qg  (nproma,nlev)
    !
    ! Local variables for security
    !
    REAL(wp) :: xlta(nproma,nlev)
    REAL(wp) :: xlqv(nproma,nlev), xlqc(nproma,nlev)
    REAL(wp) :: xlqi(nproma,nlev), xlqr(nproma,nlev)
    REAL(wp) :: xlqs(nproma,nlev), xlqg(nproma,nlev)

    REAL(wp) :: zdtr      ! reciprocal of timestep

    INTEGER :: jk, jl

    !IF (ltimer) call timer_start(timer_mig)

    ! associate pointers
    lparamcpl =  echam_phy_config(jg)%lparamcpl
    fc_mig    =  echam_phy_config(jg)%fc_mig

    jkscov    =  echam_cov_config(jg)% jkscov

    field     => prm_field(jg)
    tend      => prm_tend (jg)

    !$ser verbatim call serialize_mig_input(jg, jb, jcs, jce, nproma, nlev, field)

    !$ACC DATA PCREATE( zqrsflux, zqnc           , &
    !$ACC               tend_ta_mig              , &
    !$ACC               tend_qtrc_mig            , &
    !$ACC               xlta                     , &
    !$ACC               xlqv, xlqc               , &
    !$ACC               xlqi, xlqr               , &
    !$ACC               xlqs, xlqg               , &
    !$ACC               ddt_tend_t               , &
    !$ACC               ddt_tend_qv, ddt_tend_qc , &
    !$ACC               ddt_tend_qi, ddt_tend_qr , &
    !$ACC               ddt_tend_qs, ddt_tend_qg )

    !$ACC PARALLEL
    !$ACC LOOP GANG VECTOR
    DO jl = 1, nproma
      zqnc(jl) = cloud_num
    END DO
    !$ACC END PARALLEL

    ! preset local mig tendencies to avoid NaNs

    CALL init(tend_ta_mig)
    CALL init(tend_qtrc_mig)
    CALL init(ddt_tend_t )
    CALL init(ddt_tend_qv)
    CALL init(ddt_tend_qc)
    CALL init(ddt_tend_qi)
    CALL init(ddt_tend_qr)
    CALL init(ddt_tend_qs)
    CALL init(ddt_tend_qg)

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
       !$ACC DATA PRESENT(field%ta, field%qtrc)
       !$ACC PARALLEL
       !$ACC LOOP GANG
       DO jk = 1,nlev
         !$ACC LOOP VECTOR
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
       !$ACC END DATA

#if defined( _OPENACC )
       CALL warning('GPU:interface_echam_mig','GPU host synchronization should be removed when port is done!')
#endif
       ! satad_v_3D is not yet ported.
       !$ACC UPDATE HOST(xlta, xlqv, xlqc)
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
              & klo      = jkscov                         ,& !> IN
              & kup      = nlev                            & !> IN
              )
#if defined( _OPENACC )
       CALL warning('GPU:interface_echam_mig','GPU device synchronization should be removed when port is done!')
#endif
       !$ACC UPDATE DEVICE(xlta, xlqv, xlqc)
!
    !$ser verbatim call serialize_mig_after_satad1(jg, jb, jcs, jce, nproma, nlev, field,&
    !$ser verbatim        xlta, xlqv, xlqc, xlqi, xlqr, xlqs, xlqg, zqnc,&
    !$ser verbatim        zqrsflux)
       !
       CALL graupel (                            &
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
              & p      =field% presm_old(:,:,jb)    , & !< in:  pressure
              & rho    =field% rho      (:,:,jb)    , & !< in:  density
              & qv     =xlqv            (:,:)       , & !< inout:sp humidity
              & qc     =xlqc            (:,:)       , & !< inout:cloud water
              & qi     =xlqi            (:,:)       , & !< inout: ice
              & qr     =xlqr            (:,:)       , & !< inout:rain
              & qs     =xlqs            (:,:)       , & !< inout: snow
              & qg     =xlqg            (:,:)       , & !< inout: graupel
              & qnc    =zqnc                        , & !< inout: cloud droplet number
              & prr_gsp=field%rain_gsp_rate (:, jb)   , & !< inout precp rate rain
              & prs_gsp=field%snow_gsp_rate (:, jb)   , & !< inout precp rate snow
              & prg_gsp=field%graupel_gsp_rate (:, jb), & !< inout precp rate graupel
              & qrsflux=zqrsflux                    , & !< out  precipitation flux
              & l_cv=.TRUE.                         , &
              & ldiag_ttend = echam_mig_config(jg)%ldiag_ttend , & !< in:  if temp.  tendency shall be diagnosed
              & ldiag_qtend = echam_mig_config(jg)%ldiag_qtend , & !< in:  if moisture tendencies shall be diagnosed
              & ddt_tend_t  = ddt_tend_t (:,:)         , & !< out: tendency temperature
              & ddt_tend_qv = ddt_tend_qv(:,:)         , & !< out: tendency QV
              & ddt_tend_qc = ddt_tend_qc(:,:)         , & !< out: tendency QC
              & ddt_tend_qi = ddt_tend_qi(:,:)         , & !< out: tendency QI
              & ddt_tend_qr = ddt_tend_qr(:,:)         , & !< out: tendency QR
              & ddt_tend_qs = ddt_tend_qs(:,:)         , & !< out: tendency QS
              & ddt_tend_qg = ddt_tend_qg(:,:)         )   !< out: tendency QG

    !$ser verbatim call serialize_mig_before_satad2(jg, jb, jcs, jce, nproma, nlev, field,&
    !$ser verbatim      & xlta, xlqv, xlqc, xlqi, xlqr, xlqs, xlqg, zqnc,&
    !$ser verbatim      & zqrsflux)

#if defined( _OPENACC )
       CALL warning('GPU:interface_echam_mig','GPU host synchronization should be removed when port is done!')
#endif
       !$ACC UPDATE HOST(xlta, xlqv, xlqc)
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
              & klo      = jkscov                         ,& !> IN
              & kup      = nlev                            & !> IN
              )
#if defined( _OPENACC )
       CALL warning('GPU:interface_echam_mig','GPU device synchronization should be removed when port is done!')
#endif
       !$ACC UPDATE DEVICE(xlta, xlqv, xlqc)
       !
       ! Calculate rain and snow
       !
       !$ACC DATA PRESENT( field%rsfl, field%ssfl, field%rain_gsp_rate, field%snow_gsp_rate, &
       !$ACC               field%graupel_gsp_rate )
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG VECTOR
       DO jl = jcs, jce
         field% rsfl(jl,jb) = field%rain_gsp_rate (jl,jb)
         field% ssfl(jl,jb) = field%snow_gsp_rate (jl,jb) + field%graupel_gsp_rate (jl,jb)
       END DO
       !$ACC END PARALLEL
       !$ACC END DATA
       !
       ! Calculate tendencies and convert temperature tendency, as computed for constant volume in satad/graupel
       !  to constant pressure
       !
       !$ACC DATA PRESENT( field%ta, field%qtrc, field% cpair )
       !$ACC PARALLEL DEFAULT(PRESENT)
       !$ACC LOOP GANG
       DO jk = jkscov, nlev
         !$ACC LOOP VECTOR
         DO jl = jcs, jce
         tend_ta_mig(jl,jk) = (xlta(jl,jk)-field% ta(jl,jk,jb))*zdtr                &
                                     & * cvd/field% cpair(jl,jk,jb)
         tend_qtrc_mig(jl,jk,iqv) = MAX(-field% qtrc(jl,jk,jb,iqv)*zdtr,            &
                                     & (xlqv(jl,jk)-field% qtrc(jl,jk,jb,iqv))*zdtr)
         tend_qtrc_mig(jl,jk,iqc) = MAX(-field% qtrc(jl,jk,jb,iqc)*zdtr,            &
                                     & (xlqc(jl,jk)-field% qtrc(jl,jk,jb,iqc))*zdtr)
         tend_qtrc_mig(jl,jk,iqi) = MAX(-field% qtrc(jl,jk,jb,iqi)*zdtr,            &
                                     & (xlqi(jl,jk)-field% qtrc(jl,jk,jb,iqi))*zdtr)
         tend_qtrc_mig(jl,jk,iqr) = MAX(-field% qtrc(jl,jk,jb,iqr)*zdtr,            &
                                     & (xlqr(jl,jk)-field% qtrc(jl,jk,jb,iqr))*zdtr)
         tend_qtrc_mig(jl,jk,iqs) = MAX(-field% qtrc(jl,jk,jb,iqs)*zdtr,            &
                                     & (xlqs(jl,jk)-field% qtrc(jl,jk,jb,iqs))*zdtr)
         tend_qtrc_mig(jl,jk,iqg) = MAX(-field% qtrc(jl,jk,jb,iqg)*zdtr,            &
                                     & (xlqg(jl,jk)-field% qtrc(jl,jk,jb,iqg))*zdtr)
         ENDDO
       ENDDO
       !$ACC END PARALLEL
       !$ACC END DATA
       !
       ! store in memory for output or recycling
       !
       IF (ASSOCIATED(tend% ta_mig)) THEN 
         !$ACC DATA PRESENT( tend%ta_mig, tend_ta_mig )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = jkscov,nlev
           !$ACC LOOP VECTOR
           DO jl = jcs,jce
             tend% ta_mig(jl,jk,jb) = tend_ta_mig(jl,jk)
           ENDDO
         ENDDO
         !$ACC END PARALLEL
         !$ACC END DATA
       ENDIF
       !
       IF (ASSOCIATED(tend% qtrc_mig)) THEN
         !$ACC DATA PRESENT( tend%qtrc_mig, tend_qtrc_mig )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = jkscov,nlev
           !$ACC LOOP VECTOR
           DO jl = jcs,jce
             tend% qtrc_mig(jl,jk,jb,iqv) = tend_qtrc_mig(jl,jk,iqv)
             tend% qtrc_mig(jl,jk,jb,iqc) = tend_qtrc_mig(jl,jk,iqc)
             tend% qtrc_mig(jl,jk,jb,iqi) = tend_qtrc_mig(jl,jk,iqi)
             tend% qtrc_mig(jl,jk,jb,iqr) = tend_qtrc_mig(jl,jk,iqr)
             tend% qtrc_mig(jl,jk,jb,iqs) = tend_qtrc_mig(jl,jk,iqs)
             tend% qtrc_mig(jl,jk,jb,iqg) = tend_qtrc_mig(jl,jk,iqg)
           ENDDO
         ENDDO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       ! store in memory for output in case ldiag_ttend is true
       !
       IF (echam_mig_config(jg)%ldiag_ttend .AND. ASSOCIATED(tend%ddt_tend_t )) THEN
         !$ACC DATA PRESENT( tend%ddt_tend_t, ddt_tend_t )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = jkscov,nlev
           !$ACC LOOP VECTOR
           DO jl = jcs,jce
             tend% ddt_tend_t(jl,jk,jb) = ddt_tend_t(jl,jk)
           ENDDO
         ENDDO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       ! store in memory for output in case ldiag_qtend is true
       !
         IF (echam_mig_config(jg)%ldiag_qtend ) THEN
           IF (ASSOCIATED(tend%ddt_tend_qv)) THEN
             !$ACC DATA PRESENT( tend%ddt_tend_qv, ddt_tend_qv )
             !$ACC PARALLEL DEFAULT(PRESENT)
             !$ACC LOOP GANG
             DO jk = jkscov,nlev
               !$ACC LOOP VECTOR
               DO jl = jcs,jce
                 tend%ddt_tend_qv(jl,jk,jb) = ddt_tend_qv(jl,jk)
               ENDDO
             ENDDO
             !$ACC END PARALLEL
             !$ACC END DATA
           END IF
           !
           IF (ASSOCIATED(tend%ddt_tend_qc)) THEN
             !$ACC DATA PRESENT( tend%ddt_tend_qc, ddt_tend_qc )
             !$ACC PARALLEL DEFAULT(PRESENT)
             !$ACC LOOP GANG
             DO jk = jkscov,nlev
               !$ACC LOOP VECTOR
               DO jl = jcs,jce
                 tend%ddt_tend_qc(jl,jk,jb) = ddt_tend_qc(jl,jk)
               ENDDO
             ENDDO
             !$ACC END PARALLEL
             !$ACC END DATA
           END IF
           !
           IF (ASSOCIATED(tend%ddt_tend_qi)) THEN
             !$ACC DATA PRESENT( tend%ddt_tend_qi, ddt_tend_qi )
             !$ACC PARALLEL DEFAULT(PRESENT)
             !$ACC LOOP GANG
             DO jk = jkscov,nlev
               !$ACC LOOP VECTOR
               DO jl = jcs,jce
                 tend%ddt_tend_qi(jl,jk,jb) = ddt_tend_qi(jl,jk)
               ENDDO
             ENDDO
             !$ACC END PARALLEL
             !$ACC END DATA
           END IF
           !
           IF (ASSOCIATED(tend%ddt_tend_qr)) THEN
             !$ACC DATA PRESENT( tend%ddt_tend_qr, ddt_tend_qr )
             !$ACC PARALLEL DEFAULT(PRESENT)
             !$ACC LOOP GANG
             DO jk = jkscov,nlev
               !$ACC LOOP VECTOR
               DO jl = jcs,jce
                 tend%ddt_tend_qr(jl,jk,jb) = ddt_tend_qr(jl,jk)
               ENDDO
             ENDDO
             !$ACC END PARALLEL
             !$ACC END DATA
           END IF
           !
           IF (ASSOCIATED(tend%ddt_tend_qs)) THEN
             !$ACC DATA PRESENT( tend%ddt_tend_qs, ddt_tend_qs )
             !$ACC PARALLEL DEFAULT(PRESENT)
             !$ACC LOOP GANG
             DO jk = jkscov,nlev
               !$ACC LOOP VECTOR
               DO jl = jcs,jce
                 tend%ddt_tend_qs(jl,jk,jb) = ddt_tend_qs(jl,jk)
               ENDDO
             ENDDO
             !$ACC END PARALLEL
             !$ACC END DATA
           END IF
           !
           IF (ASSOCIATED(tend%ddt_tend_qg)) THEN
             !$ACC DATA PRESENT( tend%ddt_tend_qg, ddt_tend_qg )
             !$ACC PARALLEL DEFAULT(PRESENT)
             !$ACC LOOP GANG
             DO jk = jkscov,nlev
               !$ACC LOOP VECTOR
               DO jl = jcs,jce
                 tend%ddt_tend_qg(jl,jk,jb) = ddt_tend_qg(jl,jk)
               ENDDO
             ENDDO
             !$ACC END PARALLEL
             !$ACC END DATA
           END IF 
         END IF  ! ldiag_qtend
         !
       ELSE    ! is_active
         !
         ! retrieve from memory for recycling
         !
         IF (ASSOCIATED(tend% ta_mig)) THEN 
           !$ACC DATA PRESENT( tend%ta_mig, tend_ta_mig )
           !$ACC PARALLEL DEFAULT(PRESENT)
           !$ACC LOOP GANG
           DO jk = jkscov,nlev
             !$ACC LOOP VECTOR
             DO jl = jcs,jce
               tend_ta_mig(jl,jk) = tend% ta_mig(jl,jk,jb)
             ENDDO
           ENDDO
           !$ACC END PARALLEL
           !$ACC END DATA
         ENDIF
         !
         IF (ASSOCIATED(tend% qtrc_mig)) THEN
           !$ACC DATA PRESENT( tend%qtrc_mig, tend_qtrc_mig )
           !$ACC PARALLEL DEFAULT(PRESENT)
           !$ACC LOOP GANG
           DO jk = jkscov,nlev
             !$ACC LOOP VECTOR
             DO jl = jcs,jce
               tend_qtrc_mig(jl,jk,iqv) = tend% qtrc_mig(jl,jk,jb,iqv)
               tend_qtrc_mig(jl,jk,iqc) = tend% qtrc_mig(jl,jk,jb,iqc)
               tend_qtrc_mig(jl,jk,iqi) = tend% qtrc_mig(jl,jk,jb,iqi)
               tend_qtrc_mig(jl,jk,iqr) = tend% qtrc_mig(jl,jk,jb,iqr)
               tend_qtrc_mig(jl,jk,iqs) = tend% qtrc_mig(jl,jk,jb,iqs)
               tend_qtrc_mig(jl,jk,iqg) = tend% qtrc_mig(jl,jk,jb,iqg)
             ENDDO
           ENDDO
           !$ACC END PARALLEL
           !$ACC END DATA
         END IF
         !
       END IF  ! is_active
       !

       ! accumulate tendencies for later updating the model state
       SELECT CASE(fc_mig)
       CASE(0)
          ! diagnostic, do not use tendency
       CASE(1)
          !$ACC DATA PRESENT( tend%ta_phy, tend%qtrc_phy, tend_ta_mig, tend_qtrc_mig )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG
          DO jk = 1, nlev
            !$ACC LOOP VECTOR
            DO jl = jcs, jce
              ! use tendency to update the model state
              tend%   ta_phy(jl,jk,jb)      = tend%   ta_phy(jl,jk,jb)     + tend_ta_mig  (jl,jk)
              tend% qtrc_phy(jl,jk,jb,iqv)  = tend% qtrc_phy(jl,jk,jb,iqv) + tend_qtrc_mig(jl,jk,iqv)
              tend% qtrc_phy(jl,jk,jb,iqc)  = tend% qtrc_phy(jl,jk,jb,iqc) + tend_qtrc_mig(jl,jk,iqc)
              tend% qtrc_phy(jl,jk,jb,iqi)  = tend% qtrc_phy(jl,jk,jb,iqi) + tend_qtrc_mig(jl,jk,iqi)
              tend% qtrc_phy(jl,jk,jb,iqr)  = tend% qtrc_phy(jl,jk,jb,iqr) + tend_qtrc_mig(jl,jk,iqr)
              tend% qtrc_phy(jl,jk,jb,iqs)  = tend% qtrc_phy(jl,jk,jb,iqs) + tend_qtrc_mig(jl,jk,iqs)
              tend% qtrc_phy(jl,jk,jb,iqg)  = tend% qtrc_phy(jl,jk,jb,iqg) + tend_qtrc_mig(jl,jk,iqg)
            END DO
          END DO
          !$ACC END PARALLEL
          !$ACC END DATA
!!$       CASE(2)
!!$          ! use tendency as forcing in the dynamics
!!$          ...
       END SELECT
       !
       ! update physics state for input to the next physics process
       IF (lparamcpl) THEN
         SELECT CASE(fc_mig)
         CASE(0)
            ! diagnostic, do not use tendency
         CASE(1)
             !$ACC DATA PRESENT( field%ta, field%qtrc, tend_ta_mig, tend_qtrc_mig )
             !$ACC PARALLEL DEFAULT(PRESENT)
             !$ACC LOOP GANG
             DO jk = 1, nlev
               !$ACC LOOP VECTOR
               DO jl = jcs, jce
                 field%   ta(jl,jk,jb)      = field%   ta(jl,jk,jb)      + tend_ta_mig  (jl,jk)    *pdtime
                 field% qtrc(jl,jk,jb,iqv)  = field% qtrc(jl,jk,jb,iqv)  + tend_qtrc_mig(jl,jk,iqv)*pdtime
                 field% qtrc(jl,jk,jb,iqc)  = field% qtrc(jl,jk,jb,iqc)  + tend_qtrc_mig(jl,jk,iqc)*pdtime
                 field% qtrc(jl,jk,jb,iqi)  = field% qtrc(jl,jk,jb,iqi)  + tend_qtrc_mig(jl,jk,iqi)*pdtime
                 field% qtrc(jl,jk,jb,iqr)  = field% qtrc(jl,jk,jb,iqr)  + tend_qtrc_mig(jl,jk,iqr)*pdtime
                 field% qtrc(jl,jk,jb,iqs)  = field% qtrc(jl,jk,jb,iqs)  + tend_qtrc_mig(jl,jk,iqs)*pdtime
                 field% qtrc(jl,jk,jb,iqg)  = field% qtrc(jl,jk,jb,iqg)  + tend_qtrc_mig(jl,jk,iqg)*pdtime
               END DO
             END DO
             !$ACC END PARALLEL
             !$ACC END DATA
             !
!!$         CASE(2)
!!$            ! use tendency as forcing in the dynamics
!!$            ...
         END SELECT
       END IF
       !
    ELSE       ! is_in_sd_ed_interval
       !
       IF (ASSOCIATED(tend% ta_mig)) THEN 
         !$ACC DATA PRESENT( tend%ta_mig )
         !$ACC PARALLEL DEFAULT(PRESENT)
         !$ACC LOOP GANG
         DO jk = 1, nlev
           !$ACC LOOP VECTOR
           DO jl = jcs, jce
             tend% ta_mig(jl,jk,jb) = 0.0_wp
           END DO
         END DO
         !$ACC END PARALLEL
         !$ACC END DATA
       END IF
       !
       IF (ASSOCIATED(tend% qtrc_mig )) THEN
          !$ACC DATA PRESENT( tend%qtrc_mig )
          !$ACC PARALLEL DEFAULT(PRESENT)
          !$ACC LOOP GANG
          DO jk = 1, nlev
            !$ACC LOOP VECTOR
            DO jl = jcs, jce
              tend% qtrc_mig(jl,jk,jb,iqv) = 0.0_wp
              tend% qtrc_mig(jl,jk,jb,iqc) = 0.0_wp
              tend% qtrc_mig(jl,jk,jb,iqi) = 0.0_wp
              tend% qtrc_mig(jl,jk,jb,iqr) = 0.0_wp
              tend% qtrc_mig(jl,jk,jb,iqs) = 0.0_wp
              tend% qtrc_mig(jl,jk,jb,iqg) = 0.0_wp
            END DO
          END DO
          !$ACC END PARALLEL
          !$ACC END DATA
       END IF
       !
    END IF     ! is_in_sd_ed_interval

    !$ACC END DATA

    !$ser verbatim call serialize_mig_output(jg, jb, jcs, jce, nproma, nlev, field)
    ! disassociate pointers

    NULLIFY(field)
    NULLIFY(tend)

    IF (ltimer) call timer_stop(timer_mig)

  END SUBROUTINE interface_echam_mig

END MODULE mo_interface_echam_mig
