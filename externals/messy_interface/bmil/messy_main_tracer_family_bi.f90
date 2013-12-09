!***********************************************************************
MODULE messy_main_tracer_family_bi
!***********************************************************************

  ! MODULE FOR TRACER FAMILIES (MESSy-SMIL)
  !
  ! Authors: Patrick Joeckel,  MPICH, Jan 2004, June 2007
  !          Astrid Kerkweg,   MPICH, May 2004, June 2007
  !          Joachim Buchholz, MPICH, November 2004
  !

  ! BMIL
#if defined(ECHAM5) || defined(COSMO)
    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR, LGTRSTR
#endif
#ifdef MBM_TRACER
    USE messy_main_tracer_mem_bi, ONLY: S1TRSTR
#endif
#if defined(ECHAM5) || defined(MBM_CMAT)
    USE messy_main_tracer_mem_bi, ONLY: CMATTRSTR
#endif
  ! SMCL
  USE messy_main_tracer_family

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: main_tracer_family_initialize      ! read namelist and initialize
  PUBLIC :: main_tracer_family_new_tracer      ! request family-tracers
  PUBLIC :: main_tracer_family_init_mem        ! attributes: family + fract.
  PUBLIC :: main_tracer_family_init_cpl        ! attributes: tracers
! um_ak_20080709 + 
!!$  PUBLIC :: main_tracer_family_global_start    ! TYPE-1: t2f
  !                                               ! TYPE-2: sum
  ! FOR USE IN BML, BEFORE/AFTER ADVECTION
  PUBLIC :: main_tracer_family_beforeadv
! um_ak_20080709-
  PUBLIC :: main_tracer_family_afteradv
  !
  PUBLIC :: main_tracer_family_free_mem        ! free memory
  !

CONTAINS

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_family_initialize

    ! BMIL
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: start_message_bi, end_message_bi, error_bi
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_family_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i, j

    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL tracer_family_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    CALL start_message_bi(submodstr, 'INITIALISATION', substr)

    IF (p_parallel_io) THEN
       CALL tracfamily_init(status)
    END IF
    CALL p_bcast(status, p_io)
    IF (status /= 0) CALL error_bi('tracfamily_init reported an error', substr)

    ! BROADCAST RESULTS
    CALL p_bcast(l_verbose,   p_io)
    CALL p_bcast(i_diag_pe,   p_io)
    CALL p_bcast(i_diag_jrow, p_io)
    CALL p_bcast(NTF, p_io)
    DO i=1, NMAXTFAM
       CALL p_bcast(XTF(i)%IO%set, p_io)
       CALL p_bcast(XTF(i)%IO%type, p_io)
       CALL p_bcast(XTF(i)%IO%l_rescale, p_io)
       CALL p_bcast(XTF(i)%IO%name, p_io)
       CALL p_bcast(XTF(i)%IO%subname, p_io)
       CALL p_bcast(XTF(i)%fidt, p_io)
       CALL p_bcast(XTF(i)%nt, p_io)
       DO j=1, NMAXTRAC
          CALL p_bcast(XTF(i)%idt(j), p_io)
          CALL p_bcast(XTF(i)%IO%tracer(j), p_io)
          CALL p_bcast(XTF(i)%weight(j), p_io)
       END DO
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) '----> ',NTF,' TRACER FAMILY/IES REQUESTED !'
    END IF

    CALL end_message_bi(submodstr, 'INITIALISATION', substr)    

  END SUBROUTINE main_tracer_family_initialize
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_family_new_tracer

    ! BMIL
    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                      , error_bi
    ! SMCL
    USE messy_main_tracer,        ONLY: tracer_error_str

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'main_tracer_family_new_tracer'
    INTEGER                      :: status

    CALL start_message_bi(submodstr, 'REQUEST TRACERS', substr)

    CALL tracfamily_newtrac(status, p_parallel_io)
    IF (status /=0) CALL error_bi(tracer_error_str(status), substr)

    CALL end_message_bi(submodstr, 'REQUEST TRACERS', substr)

  END SUBROUTINE main_tracer_family_new_tracer
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_family_init_mem

    ! BMIL
    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    USE messy_main_blather_bi,    ONLY: error_bi
    ! SMCL
    USE messy_main_tracer,        ONLY: tracer_error_str

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_family_init_mem'
    INTEGER :: status

    ! flag=1: reset meta information of family-members
    ! flag=2: set meta information of family-members to fraction

#if defined(ECHAM5) || defined(COSMO)
    ! set meta information of family-members to fraction
    CALL tracfamily_meta(status, 2, substr, GPTRSTR, p_parallel_io)
    IF (status /= 0) CALL error_bi(tracer_error_str(status), substr)

!!$    CALL tracfamily_meta(status, 2, substr, LGTRSTR, p_parallel_io)
!!$    IF (status /= 0) CALL error_bi(tracer_error_str(status), substr)
#endif
#ifdef MBM_TRACER
    ! set meta information of family-members to fraction
    CALL tracfamily_meta(status, 2, substr, S1TRSTR, p_parallel_io)
    IF (status /= 0) CALL error_bi( tracer_error_str(status), substr)
#endif

#if defined(ECHAM5) || defined(MBM_CMAT)
    ! mz_ab_20100601+
    ! set meta information of family-members to fraction
    CALL tracfamily_meta(status, 2, substr, CMATTRSTR, p_parallel_io)
    IF (status /= 0) CALL error_bi( tracer_error_str(status), substr)
    ! mz_ab_20100601-
#endif 

    ! ### ADD MORE TRACER SETS HERE

  END SUBROUTINE main_tracer_family_init_mem
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_family_init_cpl

    ! BMIL
    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    USE messy_main_blather_bi,    ONLY: error_bi
    ! SMCL
    USE messy_main_tracer,        ONLY: tracer_error_str

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_family_init_cpl'
    INTEGER :: status

    ! reset meta information of family-members

#if defined(ECHAM5) 
    CALL tracfamily_meta(status, 1, substr, GPTRSTR, p_parallel_io)
    IF (status /= 0) CALL error_bi(tracer_error_str(status),substr)

!!$    CALL tracfamily_meta(status, 1, substr, LGTRSTR, p_parallel_io)
!!$    IF (status /= 0) CALL error_bi(tracer_error_str(status), substr)

    ! NOW STATUS IS: TRACERS ACTIVE
    CALL tracfamily_initmode(4)
#endif

#if defined(COSMO)
    CALL tracfamily_meta(status, 1, substr, GPTRSTR, p_parallel_io)
    IF (status /= 0) CALL error_bi(tracer_error_str(status),substr)

    ! NOW STATUS IS: TRACERS ACTIVE
    CALL tracfamily_initmode(2)
#endif

#ifdef MBM_TRACER
    CALL tracfamily_meta(status, 1, substr, S1TRSTR, p_parallel_io)
    IF (status /= 0) CALL error_bi(tracer_error_str(status), substr)

    ! NOW STATUS IS: TRACERS ACTIVE
    CALL tracfamily_initmode(4)
#endif

#if defined(ECHAM5) || defined(MBM_CMAT)
    ! mz_ab_20100601+
    CALL tracfamily_meta(status, 1, substr, CMATTRSTR, p_parallel_io)
    IF (status /= 0) CALL error_bi(tracer_error_str(status), substr)

    ! NOW STATUS IS: TRACERS ACTIVE
    CALL tracfamily_initmode(4)
    ! mz_ab_20100601-
#endif

    ! ### ADD MORE TRACER SETS HERE

  END SUBROUTINE main_tracer_family_init_cpl
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  ! um_ak_20080709+
  ! subroutine renamed
!!$ SUBROUTINE main_tracer_family_global_start
 SUBROUTINE main_tracer_family_beforeadv
  ! um_ak_20080709-

   ! BMIL
   USE messy_main_data_bi,       ONLY: ngpblks, nproma, npromz
#ifndef MESSYTIMER
   USE messy_main_data_bi,       ONLY: time_step_len
#else
   USE messy_main_timer,         ONLY: time_step_len
#endif
   USE messy_main_mpi_bi,        ONLY: p_pe
   USE messy_main_blather_bi,    ONLY: error_bi
   ! SMCL
   USE messy_main_tracer,        ONLY: tracer_error_str

   IMPLICIT NONE

   ! LOCAL
   CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_family_beforeadv'
   INTEGER :: jjrow, jp, status

   ! NOTE: CALL SUMMATION OF TYPE-2 TRACERS BEFORE TRANSFORMATION OF
   !       TPYE-1 FAMILIES (T2F), I.E. BEFORE THEY ARE
   !       POTENTIALLY CONVERTED TO FRACTIONS (AS MEMBER OF FAMILY TYPE-1);
   !       TYPE-2 FAMILY MEMBERS WITHOUT RESCALING OPTION ARE ALLOWED TO
   !       BE MEMBER OF A TYPE-1 FAMILY

#if defined(ECHAM5) || defined(COSMO)
   DO jjrow=1, ngpblks

      CALL tracfamily_2_sum(GPTRSTR, jjrow)

      IF (jjrow == ngpblks) THEN
         jp = npromz
      ELSE
         jp = nproma
      END IF

      CALL tracfamily_1_t2f(status, substr, p_pe, GPTRSTR, &
           time_step_len, jjrow, jp)
      IF (status /= 0) CALL error_bi(tracer_error_str(status),substr)

   END DO

#endif
#ifdef MBM_TRACER
   DO jjrow=1, ngpblks 

      CALL tracfamily_2_sum(S1TRSTR, jjrow)

      IF (jjrow == ngpblks) THEN
         jp = npromz
      ELSE
         jp = nproma
      END IF

      CALL tracfamily_1_t2f(status, substr, p_pe, S1TRSTR, &
           time_step_len, jjrow, jp)
      IF (status /= 0) CALL error_bi(tracer_error_str(status),substr)

      ! ### ADD MORE TRACER SETS HERE ...

   END DO

   ! ### .. OR HERE

#endif

   ! um_ak_20080709+
!!$ END SUBROUTINE main_tracer_family_global_start
 END SUBROUTINE main_tracer_family_beforeadv
   ! um_ak_20080709-
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
 SUBROUTINE main_tracer_family_afteradv

   ! BMIL
   USE messy_main_data_bi,       ONLY: ngpblks, nproma, npromz
#ifndef MESSYTIMER
   USE messy_main_data_bi,       ONLY: time_step_len
#else
   USE messy_main_timer,         ONLY: time_step_len
#endif
   USE messy_main_mpi_bi,        ONLY: p_pe
   USE messy_main_blather_bi,    ONLY: error_bi
   ! SMCL
   USE messy_main_tracer,        ONLY: tracer_error_str

   IMPLICIT NONE

   ! LOCAL
   CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_family_afteradv'
   INTEGER :: jjrow, jp, status

   ! NOTE: CALL RE-SCALING OF TYPE-2 TRACERS AFTER TRANSFORMATION OF
   !       TPYE-1 FAMILIES (F2T), I.E. AFTER THEY ARE
   !       POTENTIALLY CONVERTED BACK FROM FRACTIONS TO TRACERS
   !       (AS MEMBER OF FAMILY TYPE-1);
   !       TYPE-2 FAMILY MEMBERS WITHOUT RESCALING OPTION ARE ALLOWED TO
   !       BE MEMBER OF A TYPE-1 FAMILY

#if defined(ECHAM5) || defined(COSMO)
   DO jjrow=1, ngpblks
      IF (jjrow == ngpblks) THEN
         jp = npromz
      ELSE
         jp = nproma
      END IF

      CALL tracfamily_1_f2t(status, substr, p_pe, GPTRSTR, &
           time_step_len, jjrow, jp)
      IF (status /= 0) CALL error_bi(tracer_error_str(status),substr)

      CALL tracfamily_2_rsc(GPTRSTR, time_step_len, jjrow)
   END DO

#endif
#ifdef MBM_TRACER
   DO jjrow=1, ngpblks
      IF (jjrow == ngpblks) THEN
         jp = npromz
      ELSE
         jp = nproma
      END IF

      CALL tracfamily_1_f2t(status, substr, p_pe, S1TRSTR, &
           time_step_len, jjrow, jp)
      IF (status /= 0) CALL error_bi(tracer_error_str(status),substr)

      CALL tracfamily_2_rsc(S1TRSTR, time_step_len, jjrow)
   END DO
#endif
   
! bn_ab_201108221 I can't see why this is sensible at the moment
!!$#if defined(ECHAM5) || defined(MBM_CMAT)
!!$   ! mz_ab_20100601+
!!$   DO jjrow=1, ngpblks
!!$      IF (jjrow == ngpblks) THEN
!!$         jp = npromz
!!$      ELSE
!!$         jp = nproma
!!$      END IF
!!$
!!$      CALL tracfamily_1_f2t(status, substr, p_pe, CMATTRSTR, &
!!$           time_step_len, jjrow, jp)
!!$      IF (status /= 0) CALL error_bi(tracer_error_str(status),substr)
!!$
!!$      CALL tracfamily_2_rsc(CMATTRSTR, time_step_len, jjrow)
!!$   END DO
!!$   ! mz_ab_20100601-
!!$#endif   

   ! ### ADD MORE TRACER SETS HERE

 END SUBROUTINE main_tracer_family_afteradv
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
 SUBROUTINE main_tracer_family_free_mem

   IMPLICIT NONE

   CALL tracfamily_freemem

 END SUBROUTINE main_tracer_family_free_mem
! ----------------------------------------------------------------------

!***********************************************************************
END MODULE messy_main_tracer_family_bi
!***********************************************************************
