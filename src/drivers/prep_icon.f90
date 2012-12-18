!>
!! @brief Main program for the prep_icon preprocessor
!!
!! @author
!!  Guenther Zaengl (DWD)
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
PROGRAM prep_icon

  USE mo_prepicon_config,       ONLY: i_oper_mode
  USE mo_exception,             ONLY: message, finish
  USE mo_run_config,            ONLY: iforcing
  USE mo_impl_constants,        ONLY: inwp, MODE_GENERATE_COORDS, MODE_CONVERTIFS, &
    &                                 MODE_VERTINTERP, MODE_REMAP 
  USE mo_extpar_config,         ONLY: itopo
  USE mo_parallel_config,       ONLY: l_test_openmp, num_io_procs, p_test_run
  USE mo_master_nml,            ONLY: lrestart, read_master_namelist
  USE mo_mpi,                   ONLY: start_mpi, p_stop, set_mpi_work_communicators
  USE mo_read_namelists,        ONLY: read_atmo_namelists
  USE mo_nml_crosscheck,        ONLY: atm_crosscheck
  USE mo_run_config,            ONLY: configure_run
  USE mo_remap,                 ONLY: remap_main
!$  USE mo_exception,         ONLY: message_text     ! use only if compiled with OpenMP

#ifdef __INTEL_COMPILER
  USE, INTRINSIC :: ieee_arithmetic
#endif

!-------------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER (LEN=*), PARAMETER :: namelist_filename = TRIM("NAMELIST_PREPICON")

  INTEGER :: error_status

  !declaration of OpenMP Runtime Library Routines:
!$  INTEGER omp_get_max_threads
!$  INTEGER omp_get_num_threads
!$  INTEGER omp_get_num_procs
!$  INTEGER omp_get_thread_num
!$  LOGICAL omp_get_dynamic

!$  INTEGER :: max_threads_omp, num_procs_omp
!$  LOGICAL :: ldynamic_omp

  !--------------------------------------------------------------------
  !BOC

  !print out some information about OpenMP parallelization
!$  max_threads_omp  = omp_get_max_threads()
!$  num_procs_omp    = omp_get_num_procs()
!$  ldynamic_omp     = omp_get_dynamic()
!$  WRITE(message_text,'(A,I3,A,I3)')                &
!$    & "OpenMP:  MAX_THREADS = ", max_threads_omp,  &
!$    & ",  NUM_PROCS = ", num_procs_omp
!$  CALL message('control_model',message_text)
!$  WRITE(message_text,'(A,L3)')  &
!$    & "OpenMP:  DYNAMIC = ", ldynamic_omp
!$  CALL message('control_model',message_text)


#ifdef __INTEL_COMPILER
  ! Important on Intel: disable underflow exceptions:
!  CALL disable_ufl_exception
  CALL ieee_set_halting_mode(ieee_underflow, .FALSE.)
#endif

    !-------------------------------------------------------------------
    ! Initialize MPI, this should aleays be the first call
    CALL start_mpi('PREP_ICON')

    !---------------------------------------------------------------------
    ! 1.1 Read namelists (newly) specified by the user; fill the 
    !     corresponding sections of the configuration states.
    !---------------------------------------------------------------------

    lrestart = .FALSE. ! restarting is not available for prep_icon
    error_status = read_master_namelist("icon_master.namelist")
    CALL read_atmo_namelists(TRIM(namelist_filename),"icon_master.namelist")


    !---------------------------------------------------------------------
    ! 1.2 Cross-check namelist setups
    !---------------------------------------------------------------------

    CALL atm_crosscheck

    !---------------------------------------------------------------------
    ! 2. Call configure_run to finish filling the run_config state.
    !    This needs to be done very early (but anyway after atm_crosscheck)
    !    because some component of the state, e.g., num_lev, may be 
    !    modified in this subroutine which affect the following CALLs.
    !---------------------------------------------------------------------
    CALL configure_run

    ! This is hardcoded for prep_icon because running this program is nonsense otherwise
    itopo = 1
    iforcing = inwp

    !-------------------------------------------------------------------
    ! 3.1 Initialize the mpi work groups
    !-------------------------------------------------------------------
    CALL set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs)

    SELECT CASE(i_oper_mode)
    CASE (MODE_GENERATE_COORDS, MODE_CONVERTIFS, MODE_VERTINTERP)
      CALL finish("prep_icon", "Invalid operation mode!")
    CASE (MODE_REMAP) 
      CALL remap_main(TRIM(namelist_filename))
    CASE DEFAULT
      CALL finish("prep_icon", "Unknown operation mode!")
    END SELECT

    CALL p_stop

END PROGRAM prep_icon
