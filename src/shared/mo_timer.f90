!>
!! @author <name, affiliation>
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_timer

  USE mo_real_timer, ONLY: new_timer,                        &
       &                   timer_start,                      &
       &                   timer_stop,                       &
       &                   print_timer   => timer_report,    &
       &                   cleanup_timer => timer_reset_all, &
       &                   delete_timer => del_timer

  USE mo_run_nml,    ONLY: ltimer

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: new_timer, timer_start, timer_stop  !< procedures imported from mo_real_timer
  PUBLIC :: print_timer, cleanup_timer, delete_timer          !< procedures imported and renamed
  PUBLIC :: init_timer                          !< procedure of this module

  PUBLIC :: timer_total                         !< IDs of timers

  PUBLIC :: timer_solve_nh
  PUBLIC :: timer_physics                         !< IDs of timers
  PUBLIC :: timer_radiation, timer_radheat

  PUBLIC :: timer_lrtm_1, timer_lrtm_2

  PUBLIC :: timer_div, timer_grad, timer_gmres  !< IDs of timers
  PUBLIC :: timer_corio, timer_intp             !< IDS of timers
  PUBLIC :: timer_transport
  PUBLIC :: timer_cover, timer_cloud
  PUBLIC :: timer_cucall
  PUBLIC :: timer_vdiff
  PUBLIC :: timer_echam_phy
  PUBLIC :: timer_dyn2phy, timer_phy2dyn

  PUBLIC :: ltimer                              !< if .true., switch on timer

  !-------------------
  ! Module variables
  !-------------------

  ! ID of timer for total model integration time
  INTEGER :: timer_total

  INTEGER :: timer_solve_nh
  INTEGER :: timer_physics

  ! Timer ID's for horizontal operators
  INTEGER :: timer_div
  INTEGER :: timer_grad
  INTEGER :: timer_gmres
  INTEGER :: timer_corio
  INTEGER :: timer_intp

  INTEGER :: timer_transport    ! tracer transport

  ! Timer ID's for ECHAM6 physics
  INTEGER :: timer_cover
  INTEGER :: timer_cloud
  INTEGER :: timer_radiation
  INTEGER :: timer_lrtm_1, timer_lrtm_2
  INTEGER :: timer_radheat
  INTEGER :: timer_cucall
  INTEGER :: timer_vdiff
  INTEGER :: timer_echam_phy

  ! Timer ID's for physics-dynamics coupling

  INTEGER :: timer_dyn2phy
  INTEGER :: timer_phy2dyn

CONTAINS

  SUBROUTINE init_timer

    timer_total     = new_timer("total")
    timer_solve_nh  = new_timer("solve_nh")
    timer_physics   = new_timer("physics")
    timer_radiation = new_timer("radiation")
    timer_lrtm_1    = new_timer("rad_lrtm_1")
    timer_lrtm_2    = new_timer("rad_lrtm_2")

    timer_div       = new_timer("div")
    timer_grad      = new_timer("grad")
    timer_gmres     = new_timer("gmres")
    timer_corio     = new_timer("corio")
    timer_intp      = new_timer("intp")
    timer_transport = new_timer("transport")

    timer_cover     = new_timer("cover")
    timer_cloud     = new_timer("cloud")

    timer_radheat   = new_timer("radheat")

    timer_cucall    = new_timer("cucall")
    timer_vdiff     = new_timer("vdiff")
    timer_echam_phy = new_timer("echam_phy")
    timer_dyn2phy   = new_timer("dyn2phy")
    timer_phy2dyn   = new_timer("phy2dyn")

  END SUBROUTINE init_timer

END MODULE mo_timer








