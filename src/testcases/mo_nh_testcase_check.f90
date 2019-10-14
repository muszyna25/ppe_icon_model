!>
!! Interface for check of non-hydrostatic testcases.
!!
!! @par Revision History
!! Initial revision by Sebastian Borchert, DWD, (2018-06-21)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

MODULE mo_nh_testcase_check

  USE mo_kind,                     ONLY: wp
  USE mo_impl_constants,           ONLY: MAX_CHAR_LENGTH
  USE mo_nh_testcases_nml,         ONLY: nh_test_name
  USE mo_dynamics_config,          ONLY: ldeepatmo, lcoriolis
  USE mo_nonhydrostatic_config,    ONLY: l_open_ubc, ivctype
  USE mo_sleve_config,             ONLY: min_lay_thckn, top_height
  USE mo_extpar_config,            ONLY: itopo
  USE mo_run_config,               ONLY: iforcing, ltransport, lvert_nest, output_mode
  USE mo_grid_config,              ONLY: grid_sphere_radius, grid_angular_velocity, & 
    &                                    grid_rescale_factor
  USE mo_io_config,                ONLY: inextra_3d
  USE mo_name_list_output_config,  ONLY: first_output_name_list
  USE mo_nh_lahade,                ONLY: check_nh_lahade
  USE mo_upatmo_config,            ONLY: upatmo_dyn_config, upatmo_config

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: check_nh_testcase

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nh_testcase_check'

CONTAINS

  !>
  !! Interface for check of non-hydrostatic testcases.
  !!
  !! @par Revision History
  !! Initial revision by Sebastian Borchert, DWD, (2018-06-21)
  !!
  SUBROUTINE check_nh_testcase()

    ! Local variables
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':check_nh_testcase'
    
    !--------------------------------------------------------------

    ! Note:
    ! - This subroutine is called in 'src/configure_model/mo_nml_crosscheck: atm_crosscheck'
    ! - A query for 'run_nml: ltestcase' encapsulates the call of this subroutine

    SELECT CASE ( TRIM(nh_test_name) )
        
    CASE ('lahade')

      CALL check_nh_lahade(lvert_nest=lvert_nest, l_nml=output_mode%l_nml, ivctype=ivctype, &  !in
        &                  top_height=top_height, grid_sphere_radius=grid_sphere_radius,    &  !in
        &                  grid_rescale_factor=grid_rescale_factor,                         &  !in
        &                  first_output_name_list=first_output_name_list,                   &  !in
        &                  ldeepatmo=ldeepatmo, lcoriolis=lcoriolis, l_open_ubc=l_open_ubc, &  !inout
        &                  ltransport=ltransport, itopo=itopo, iforcing=iforcing,           &  !inout
        &                  inextra_3d=inextra_3d, min_lay_thckn=min_lay_thckn,              &  !inout
        &                  grid_angular_velocity=grid_angular_velocity,                     &  !inout
        &                  upatmo_dyn_config=upatmo_dyn_config, upatmo_config=upatmo_config )  !inout, (opt)in

    END SELECT

  END SUBROUTINE check_nh_testcase

END MODULE mo_nh_testcase_check
