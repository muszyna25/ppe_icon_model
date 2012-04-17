!! Contains utilities for diagnose pressure and temperature in nh model
!!
!!
!! @par Revision History
!! Initial release by Almut Gassmann, MPI-M (2009-03-06)
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_diagnose_pres_temp

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_lnd_types,       ONLY: t_lnd_prog
  USE mo_run_config,          ONLY: iqv, iqc, iqi, iqs, iqr, &
    &                               lforcing, iforcing
  USE mo_impl_constants,      ONLY: min_rlcell, MAX_CHAR_LENGTH 
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_physical_constants,  ONLY: rd, grav, vtmpc1, p0ref, rd_o_cpd
  USE mo_timer,               ONLY: timers_level, timer_start, timer_stop, timer_diagnose_pres_temp

  IMPLICIT NONE

  PRIVATE

  REAL(wp), PARAMETER :: cpd_o_rd  = 1._wp / rd_o_cpd
  REAL(wp), PARAMETER :: grav_o_rd = grav / rd

  PUBLIC :: diagnose_pres_temp

  CONTAINS

  ! moved here from mo_nh_stepping to avoid circular dependencies
  !>
  !! diagnose_pres_temp
  !!
  !! Diagnoses pressure and temperature from NH prognostic fields
  !!
  !! @par Revision History
  !! Initial release by Guenther Zaengl (2010-04-15)
  !!
  SUBROUTINE diagnose_pres_temp (p_metrics, pt_prog, pt_prog_rcf, pt_diag, pt_patch, &
    &                            opt_calc_temp, opt_calc_pres, opt_calc_tempv, &
    &                            lnd_prog, opt_calc_temp_ifc, opt_slev, opt_rlend )


    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nh_diagnose_pres_temp:diagnose_pres_temp'

    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),    INTENT(IN)    :: pt_prog      !!the prognostic variables
    TYPE(t_nh_prog),    INTENT(IN)    :: pt_prog_rcf  !!the prognostic variables which are
                                                      !! treated with reduced calling frequency

    TYPE(t_lnd_prog),   INTENT(IN), OPTIONAL :: lnd_prog 
    TYPE(t_nh_diag),    INTENT(INOUT) :: pt_diag      !!the diagnostic variables


    TYPE(t_patch),      INTENT(IN)    :: pt_patch    ! Patch

    LOGICAL, INTENT(IN), OPTIONAL   :: opt_calc_temp, opt_calc_pres, opt_calc_tempv, &
      &                                opt_calc_temp_ifc

    INTEGER, INTENT(IN), OPTIONAL :: opt_slev, opt_rlend 

    INTEGER  :: jb,jk,jc
    INTEGER  :: nlev, nlevp1              !< number of full levels
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
    INTEGER  :: i_rlstart, i_rlend
    INTEGER  :: slev

    LOGICAL  :: l_opt_calc_temp, l_opt_calc_tempv, l_opt_calc_pres, l_opt_calc_temp_ifc

    REAL(wp) :: dz, z_qsum


    IF (timers_level > 2) CALL timer_start(timer_diagnose_pres_temp)

    dz     = 0._wp
    z_qsum = 0._wp
    
    ! Check for optional arguments

    IF ( PRESENT(opt_calc_temp_ifc ) ) THEN
      l_opt_calc_temp_ifc = opt_calc_temp_ifc
    ELSE
      l_opt_calc_temp_ifc = .FALSE.
    ENDIF

    IF ( PRESENT(opt_calc_temp ) ) THEN
      l_opt_calc_temp = opt_calc_temp
    ELSE
      l_opt_calc_temp = .TRUE.
    ENDIF

    IF ( PRESENT(opt_calc_tempv ) ) THEN
      l_opt_calc_tempv = opt_calc_tempv
    ELSE
      l_opt_calc_tempv = .FALSE.
    ENDIF

    IF ( PRESENT(opt_calc_pres ) ) THEN
      l_opt_calc_pres = opt_calc_pres
    ELSE
      l_opt_calc_pres = .TRUE.
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell
    ENDIF

    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    ENDIF

    ! number of child domains
    i_nchdom = MAX(1,pt_patch%n_childdom)

    ! number of vertical levels
    nlev   = pt_patch%nlev
    nlevp1 = pt_patch%nlevp1

    ! Nest boundaries are always included
    i_rlstart = 1

    i_startblk = pt_patch%cells%start_blk(i_rlstart,1)
    i_endblk   = pt_patch%cells%end_blk(i_rlend,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jk, jc, dz, z_qsum) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( pt_patch, jb, i_startblk, i_endblk,      &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )

      IF ( lforcing .AND. iforcing /= 1  ) THEN

        IF ( l_opt_calc_temp .AND. l_opt_calc_tempv) THEN
          DO jk = slev, nlev
            DO jc = i_startidx, i_endidx

              pt_diag%tempv(jc,jk,jb) = pt_prog%theta_v(jc,jk,jb) * pt_prog%exner(jc,jk,jb)

              z_qsum                   =    pt_prog_rcf%tracer (jc,jk,jb,iqc) &
                &                         + pt_prog_rcf%tracer (jc,jk,jb,iqi) &
                &                         + pt_prog_rcf%tracer (jc,jk,jb,iqr) &
                &                         + pt_prog_rcf%tracer (jc,jk,jb,iqs)

              pt_diag%temp(jc,jk,jb)   =  pt_diag%tempv(jc,jk,jb)          &
                &                      /( 1._wp +  vtmpc1                  &
                &                      * pt_prog_rcf%tracer(jc,jk,jb,iqv)  &
                &                      - z_qsum  )
            ENDDO
          ENDDO

        ELSEIF ( l_opt_calc_temp ) THEN
          DO jk = slev, nlev
            DO jc = i_startidx, i_endidx

              z_qsum                   =    pt_prog_rcf%tracer (jc,jk,jb,iqc) &
                &                         + pt_prog_rcf%tracer (jc,jk,jb,iqi) &
                &                         + pt_prog_rcf%tracer (jc,jk,jb,iqr) &
                &                         + pt_prog_rcf%tracer (jc,jk,jb,iqs)

              pt_diag%temp(jc,jk,jb) = pt_prog%theta_v(jc,jk,jb) &
                &                    * pt_prog%exner  (jc,jk,jb) &
                &                    /( 1._wp +  vtmpc1          &
                &                    * pt_prog_rcf%tracer(jc,jk,jb,iqv) &
                &                    - z_qsum )
            ENDDO
          ENDDO

        ELSEIF ( l_opt_calc_tempv ) THEN
          DO jk = slev, nlev
            DO jc = i_startidx, i_endidx
              pt_diag%tempv(jc,jk,jb) = pt_prog%theta_v(jc,jk,jb) * pt_prog%exner(jc,jk,jb)
            ENDDO
          ENDDO
        ENDIF
      ELSE ! .NOT. lforcing or Held-Suarez test forcing
        IF ( l_opt_calc_temp .AND. l_opt_calc_tempv) THEN
          DO jk = slev, nlev
            DO jc = i_startidx, i_endidx
               pt_diag%tempv(jc,jk,jb) = pt_prog%theta_v(jc,jk,jb) * pt_prog%exner(jc,jk,jb)
               pt_diag%temp(jc,jk,jb)  = pt_diag%tempv  (jc,jk,jb)
            ENDDO
          ENDDO
        ELSEIF ( l_opt_calc_temp ) THEN
          DO jk = slev, nlev
            DO jc = i_startidx, i_endidx
              pt_diag%temp(jc,jk,jb) = pt_prog%theta_v(jc,jk,jb) * pt_prog%exner(jc,jk,jb)
            ENDDO
          ENDDO
        ELSEIF ( l_opt_calc_tempv ) THEN
          DO jk = slev, nlev
            DO jc = i_startidx, i_endidx
               pt_diag%tempv(jc,jk,jb) = pt_prog%theta_v(jc,jk,jb) * pt_prog%exner(jc,jk,jb)
            ENDDO
          ENDDO
        ENDIF

      ENDIF

      !-------------------------------------------------------------------------
      !> diagnose pressure on main and interface levels
      !!    and   pressure thickness
      !!
      !-------------------------------------------------------------------------

      IF ( l_opt_calc_pres ) THEN
        
        DO jc = i_startidx, i_endidx
          ! Height difference between surface and first full level
          dz = p_metrics%z_ifc(jc,nlevp1,jb)-p_metrics%z_mc(jc,nlev,jb)

          pt_diag%pres_sfc(jc,jb) =p0ref*EXP(cpd_o_rd*LOG(pt_prog%exner(jc,nlev,jb))) *&
            &                      EXP(-grav_o_rd*dz/pt_diag%temp(jc,nlev,jb))

          pt_diag%pres_ifc(jc,nlevp1,jb) = pt_diag%pres_sfc(jc,jb)
        ENDDO

        !-------------------------------------------------------------------------
        !> diagnose pressure with hydrostatic approximations for Physics 
        !! - at interfaces
        !-------------------------------------------------------------------------

        DO jk = nlev, slev,-1
          DO jc = i_startidx, i_endidx

            pt_diag%pres_ifc(jc,jk,jb) = pt_diag%pres_ifc(jc,jk+1,jb)                 &
              & *EXP(-grav_o_rd*p_metrics%ddqz_z_full(jc,jk,jb)/pt_diag%temp(jc,jk,jb))
          ENDDO
        ENDDO

        !-------------------------------------------------------------------------
        !> - at main levels
        !-------------------------------------------------------------------------

        DO jk = nlev,slev,-1
          DO jc = i_startidx, i_endidx
            pt_diag%pres(jc,jk,jb) = pt_diag%pres_ifc(jc,jk+1,jb)                           &
              & *EXP(-grav_o_rd*0.5_wp*p_metrics%ddqz_z_full(jc,jk,jb)/pt_diag%temp(jc,jk,jb))
          ENDDO
        ENDDO

        !
        !-------------------------------------------------------------------------
        !> Calculate layer thickness for pressure
        !-------------------------------------------------------------------------
        !
        pt_diag%dpres_mc(i_startidx:i_endidx,nlev,jb) =                                          &
          & pt_diag%pres_sfc(i_startidx:i_endidx,jb)-pt_diag%pres_ifc(i_startidx:i_endidx,nlev,jb)

        DO jk = nlev-1, slev,-1
          DO jc =  i_startidx, i_endidx
            pt_diag%dpres_mc(jc,jk,jb) = pt_diag%pres_ifc(jc,jk+1,jb) &
  &                                    - pt_diag%pres_ifc(jc,jk  ,jb)
          ENDDO
        ENDDO
        
      ENDIF ! calc_pres

      !-------------------------------------------------------------------------
      !> diagnose temperature on interface levels
      !-------------------------------------------------------------------------
      
      IF ( l_opt_calc_temp_ifc ) THEN
        
        DO jk = MAX(slev+1,2), nlev
          DO jc =  i_startidx, i_endidx
            pt_diag%temp_ifc(jc,jk,jb) = &
              p_metrics%wgtfac_c(jc,jk,jb)*pt_diag%temp(jc,jk,jb) +      &
              (1._wp-p_metrics%wgtfac_c(jc,jk,jb))*pt_diag%temp(jc,jk-1,jb)
          ENDDO
        ENDDO

        IF ( PRESENT(lnd_prog) ) THEN
          DO jc =  i_startidx, i_endidx
            pt_diag%temp_ifc(jc,     1,jb) = pt_diag%temp (jc,1,jb)
            pt_diag%temp_ifc(jc,nlevp1,jb) = lnd_prog%t_g (jc,jb)
          ENDDO
        ELSE
          DO jc =  i_startidx, i_endidx
            pt_diag%temp_ifc(jc,     1,jb) = pt_diag%temp (jc,1,jb)
            pt_diag%temp_ifc(jc,nlevp1,jb) = pt_diag%temp (jc,nlev,jb)
          ENDDO
        ENDIF

      ENDIF !l_opt_calc_temp_ifc
      
    ENDDO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    
    IF (timers_level > 2) CALL timer_stop(timer_diagnose_pres_temp)

  END SUBROUTINE diagnose_pres_temp

END MODULE  mo_nh_diagnose_pres_temp
