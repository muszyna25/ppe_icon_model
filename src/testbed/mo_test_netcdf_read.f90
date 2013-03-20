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
MODULE mo_test_netcdf_read

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: work_mpi_barrier, my_process_is_stdio
  USE mo_timer,               ONLY: init_timer, ltimer, new_timer, timer_start, timer_stop, &
    & print_timer, activate_sync_timers, timers_level, timer_barrier, timer_radiaton_recv
  USE mo_parallel_config,     ONLY: nproma, icon_comm_method

  USE mo_master_control,      ONLY: get_my_process_name, get_my_model_no
  USE mo_icon_testbed_config, ONLY: testbed_model, o3_filename

  USE mo_model_domain,        ONLY: p_patch
  USE mo_atmo_model,          ONLY: construct_atmo_model, destruct_atmo_model
  USE mo_atmo_hydrostatic,    ONLY: construct_atmo_hydrostatic, destruct_atmo_hydrostatic
  USE mo_netcdf_read

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

PUBLIC :: test_netcdf_read


CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_netcdf_read(namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename
    
    REAL(wp), POINTER :: test_o3(:,:,:,:)  ! is (nproma, vertical_levels, blocks, time (months))
    INTEGER :: return_status

    CHARACTER(*), PARAMETER :: method_name = "test_netcdf_read"

    !---------------------------------------------------------------------
    ltimer = .false.
    timers_level = 0
    activate_sync_timers = .false.
    CALL construct_atmo_model(namelist_filename,shr_namelist_filename)
    CALL construct_atmo_hydrostatic()
    !---------------------------------------------------------------------


    !---------------------------------------------------------------------
    ! test O3
    ! read o3
    CALL message(method_name,   o3_filename)
    NULLIFY(test_o3)
    return_status = netcdf_read_oncells_3D_time(o3_filename, "O3", test_o3, p_patch(1))
    write(0,*) "shape of o3 array:", SHAPE(test_o3)
    !write(0,*) "o3(1,1,1,1)=", test_o3(1,1,1,1)
    !---------------------------------------------------------------------
    ! write the o3 array to output for comparing
    return_status = netcdf_write_oncells_3d_time("O3_test.nc", "O3", test_o3, p_patch(1))
    DEALLOCATE(test_o3)
    !---------------------------------------------------------------------


    !---------------------------------------------------------------------
    ! Carry out the shared clean-up processes
    CALL destruct_atmo_hydrostatic()
    CALL destruct_atmo_model()
     

  END SUBROUTINE test_netcdf_read
  !-------------------------------------------------------------------------
        

END MODULE mo_test_netcdf_read

