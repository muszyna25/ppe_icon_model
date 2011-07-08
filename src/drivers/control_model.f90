!>
!! @page pagecontrolmodelf901 ICON Master Program
!!
!!
!! @author
!!     Leonidas Linardakis
!!     (MPI-M)
!!
!! @date 2011-15-6
!!

!>
!! This is the master progam of the ICON model.
!!
!!
!! @par Revision History
!!   
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
PROGRAM control_model


  USE mo_exception,           ONLY: finish
  USE mo_io_units,            ONLY: filename_max
!$ USE mo_exception,          ONLY: message_text, message     ! use only if compiled with OpenMP

  USE mo_mpi,                 ONLY: p_start, p_stop !, p_start_reset
! USE mo_namelist,            ONLY: open_nml,  close_nml, open_nml_output, close_nml_output
! USE mo_output,              ONLY: init_output_files, close_output_files, write_output
  
  USE mo_atmo_model,          ONLY: atmo_model
  USE mo_ocean_model,         ONLY: ocean_model
  USE mo_radiation_model,     ONLY: radiation_model

! USE mo_icon_cpl,            ONLY: ICON_atmos_index, ICON_ocean_index, &
!  &                                comp_id, comp_comm, ierr, &
!  &                                complist
! USE mo_icon_cpl_init,       ONLY: icon_cpl_init
! USE mo_icon_cpl_init_comp,  ONLY: icon_cpl_init_comp
! USE mo_icon_cpl_finalize,   ONLY: icon_cpl_finalize

  USE mo_master_control,      ONLY: init_master_control,                 &
    & get_my_namelist_filename, get_my_process_component, &!is_coupled_run,  &
    & atmo_process, ocean_process,  radiation_process

  IMPLICIT NONE

! INTEGER    :: jg
  INTEGER    :: master_control_status
  
  INTEGER    :: my_process_component
  CHARACTER(len=filename_max) my_namelist_filename

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
  CALL disable_ufl_exception
#endif

  !-------------------------------------------------------------------
  ! Initialize MPI, this should aleays be the first call
  CALL p_start('ICON')

  
  !-------------------------------------------------------------------
  ! Initialize the master control
  master_control_status = init_master_control("icon_master.namelist")
  


!   IF ( is_coupled_run() ) THEN
! 
!      CALL icon_cpl_init
!
!       DO jg = complist(ICON_ocean_index)%min_rank, &
!               complist(ICON_ocean_index)%max_rank, &
!               complist(ICON_ocean_index)%inc_rank
!
!           IF ( p_pe == jg ) THEN
!
!             CALL open_nml(complist(ICON_ocean_index)%nml_name)
!
!             CALL icon_cpl_init_comp ( 'ocean', comp_id, comp_comm, ierr )
!
!           ENDIF
!
!       ENDDO
!
!       DO jg = complist(ICON_atmos_index)%min_rank, &
!               complist(ICON_atmos_index)%max_rank, &
!               complist(ICON_atmos_index)%inc_rank
!
!           IF ( p_pe == jg ) THEN
!
!             CALL open_nml(complist(ICON_atmos_index)%nml_name)
!
!             CALL icon_cpl_init_comp ( 'atmosphere', comp_id, comp_comm, ierr )
!
!           ENDIF
!
!       ENDDO
!
!       CALL p_start_reset ( comp_comm )
!
!     ELSE
!
!       if ( complist(icon_atmos_index)%l_comp_status ) &
!             call open_nml(complist(icon_atmos_index)%nml_name)
!
!       if ( complist(icon_ocean_index)%l_comp_status ) &
!             call open_nml(complist(icon_ocean_index)%nml_name)
!
!     ENDIF

  my_namelist_filename = get_my_namelist_filename()
  my_process_component = get_my_process_component()
  
  SELECT CASE (my_process_component)

  CASE (atmo_process)
    CALL atmo_model(my_namelist_filename)

  CASE (ocean_process)
    CALL ocean_model(my_namelist_filename)

  CASE (radiation_process)
    CALL radiation_model(my_namelist_filename)

  CASE default
    CALL finish("control_model","my_process_component is unkown")
    
  END SELECT
      
!   IF ( is_coupled_run() ) &
!     CALL ICON_cpl_finalize

  ! Shut down MPI
  !
  CALL p_stop

END PROGRAM control_model

