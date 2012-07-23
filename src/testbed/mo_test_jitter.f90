!----------------------------
#include "omp_definitions.inc"
!----------------------------
!>
!! @brief Main program for the ICON atmospheric model
!!
!! @author
!!  Leonidas Linardakis (MPI-M)
!!  Hui Wan             (MPI-M)
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_test_jitter

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: work_mpi_barrier, my_process_is_stdio
  USE mo_timer,               ONLY: new_timer, timer_start, timer_stop, &
    & print_timer, cleanup_timer

  USE mo_master_control,      ONLY: get_my_process_name, get_my_model_no
  USE mo_icon_testbed_config, ONLY: testbed_iterations, calculate_iterations, &
    & no_of_blocks, no_of_layers

  USE mo_parallel_nml,       ONLY: read_parallel_namelist
  USE mo_parallel_config,    ONLY: nproma

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: test_jitter

CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_jitter(namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    ! 3D variables
    REAL(wp), DIMENSION(nproma,no_of_layers,no_of_blocks) :: a, b, c
        
    INTEGER ::  timer_barrier
        
    CHARACTER(*), PARAMETER :: method_name = "mo_test_jitter:test_jitter"


    !---------------------------------------------------------------------
    CALL read_parallel_namelist(namelist_filename)
    
    CALL message(" ---------------- ", method_name)
    WRITE(message_text,*) "testbed_iterations=", testbed_iterations
    CALL message(" -- ", message_text)
    WRITE(message_text,*) "calculate_iterations=", calculate_iterations
    CALL message(" -- ", message_text)
    WRITE(message_text,*) "nproma=", nproma, " layers=", no_of_layers, " blocks=", no_of_blocks
    CALL message(" -- ", message_text)
   !-------------------------------------------------------------------------
    timer_barrier  = new_timer("mpi_barrier")
    CALL work_mpi_barrier()
    
    !---------------------------------------------------------------------
    ! call some barriers to see how much time it takes
    CALL timer_start(timer_barrier)
    CALL work_mpi_barrier()    
    CALL timer_stop(timer_barrier)
    CALL print_timer()
    CALL cleanup_timer(timer_barrier)
    !---------------------------------------------------------------------
    CALL timer_start(timer_barrier)
    CALL work_mpi_barrier()    
    CALL timer_stop(timer_barrier)
    CALL print_timer()
    CALL cleanup_timer(timer_barrier)
    !---------------------------------------------------------------------
    CALL timer_start(timer_barrier)
    CALL work_mpi_barrier()    
    CALL timer_stop(timer_barrier)
    CALL print_timer()
    CALL cleanup_timer(timer_barrier)
    !---------------------------------------------------------------------
    CALL timer_start(timer_barrier)
    CALL work_mpi_barrier()    
    CALL timer_stop(timer_barrier)
    CALL print_timer()
    CALL cleanup_timer(timer_barrier)

    !---------------------------------------------------------------------
    CALL test_jitter_iter()
    !---------------------------------------------------------------------
  
  END SUBROUTINE test_jitter
  !---------------------------------------------------------------------
    
  !---------------------------------------------------------------------
  !>
  SUBROUTINE test_jitter_iter()

    ! 3D variables
    REAL(wp), DIMENSION(nproma,no_of_layers,no_of_blocks) :: a, b, c
    REAL(wp) :: suma, sumb, sumc
    
    INTEGER :: timer_calculate, timer_barrier
    
    INTEGER :: i, j, k, iter, calculate
    
    CHARACTER(*), PARAMETER :: method_name = "mo_test_jitter:test_jitter_iter"
    
    timer_barrier  = new_timer("mpi_barrier")
    timer_calculate  = new_timer("calculate")
!$OMP PARALLEL
!$OMP DO PRIVATE(i,k,j) ICON_OMP_DEFAULT_SCHEDULE
    DO i = 1, no_of_blocks
      DO k = 1, no_of_layers
        DO j = 1, nproma
          a(j,k,i) = 0.0_wp
          b(j,k,i) = real(k,wp)
          c(j,k,i) = real(k+i,wp)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
    !---------------------------------------------------------------------
    CALL work_mpi_barrier()    
    !---------------------------------------------------------------------

    DO iter=1, testbed_iterations
      !---------------------------------------------------------------------
      ! do some calculations
      CALL timer_start(timer_calculate)
!$OMP PARALLEL
      DO calculate=1,calculate_iterations

!$OMP DO PRIVATE(i,k,j) ICON_OMP_DEFAULT_SCHEDULE
        DO i = 1, no_of_blocks
          DO k = 1, no_of_layers
            DO j = 1, nproma
              a(j,k,i) = b(j,k,i) / c(j,k,i)
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
        
!$OMP DO PRIVATE(i,k,j) ICON_OMP_DEFAULT_SCHEDULE
        DO i = 1, no_of_blocks
          DO k = 1, no_of_layers
            DO j = 1, nproma
              b(j,k,i) = a(j,k,i) + c(j,k,i)
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
        
!$OMP DO PRIVATE(i,k,j) ICON_OMP_DEFAULT_SCHEDULE
        DO i = 1, no_of_blocks
          DO k = 1, no_of_layers
            DO j = 1, nproma
              c(j,k,i) = a(j,k,i) * b(j,k,i)
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO

      ENDDO !calculate=1,calculate_iterations
!$OMP END PARALLEL
            
      CALL timer_stop(timer_calculate)
    
      !---------------------------------------------------------------------
      CALL timer_start(timer_barrier)
      CALL work_mpi_barrier()    
      CALL timer_stop(timer_barrier)
      !---------------------------------------------------------------------
    
      !---------------------------------------------------------------------
      ! print the timers
      CALL print_timer()
      CALL cleanup_timer(timer_calculate)
      CALL cleanup_timer(timer_barrier)
      CALL work_mpi_barrier()    
      !---------------------------------------------------------------------
      
    ENDDO !iter=1, testbed_iterations
             
    !---------------------------------------------------------------------
    ! print something to avoid optimization misfortunes
    suma = SUM(a(:,:,:))
    sumb = SUM(b(:,:,:))
    sumc = SUM(c(:,:,:))
    write(0,*) "sums=", suma, sumb, sumc
    !---------------------------------------------------------------------

  END SUBROUTINE test_jitter_iter
  !-------------------------------------------------------------------------


END MODULE mo_test_jitter

