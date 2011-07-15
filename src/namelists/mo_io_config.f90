!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
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
MODULE mo_io_config

  USE mo_kind,           ONLY: wp
  USE mo_impl_constants, ONLY: MAX_NTRACER, MAX_CHAR_LENGTH, max_dom,&
    &                          SUCCESS
  USE mo_exception,      ONLY: message, finish

  IMPLICIT NONE
  PUBLIC
  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !--------------------------------------------------------------------------
  ! Derived type 
  !--------------------------------------------------------------------------
!  TYPE :: t_io_config

    ! namelist variables

    CHARACTER(len=max_char_length) :: out_expname
    INTEGER :: out_filetype               ! 1 - GRIB1, 2 - netCDF
    LOGICAL :: lkeep_in_sync              ! if .true., sync stream after each timestep
    REAL(wp):: dt_data                    ! output timestep [seconds]
    REAL(wp):: dt_diag                    ! diagnostic output timestep [seconds]
    REAL(wp):: dt_file                    ! timestep [seconds] for triggering new output file
    REAL(wp):: dt_checkpoint              ! timestep [seconds] for triggering new restart file
  
    LOGICAL :: lwrite_vorticity           ! if .true., write out vorticity
    LOGICAL :: lwrite_divergence          ! if .true., write out divergence
    LOGICAL :: lwrite_omega               ! if .true., write out the vertical velocity
    LOGICAL :: lwrite_pres                ! if .true., write out full level pressure
    LOGICAL :: lwrite_z3                  ! if .true., write out geopotential on full levels
    LOGICAL :: lwrite_tracer(max_ntracer) ! for each tracer, if .true. write out
                                          ! tracer on full levels
    LOGICAL :: lwrite_tend_phy            ! if .true., write out physics-induced tendencies
    LOGICAL :: lwrite_radiation           ! if .true., write out fields related to radiation
    LOGICAL :: lwrite_precip              ! if .true., write out precip
    LOGICAL :: lwrite_cloud               ! if .true., write out cloud variables
                                          ! in pressure coordinate
    LOGICAL :: lwrite_tke                 ! if .true., write out TKE
    LOGICAL :: lwrite_surface             ! if .true., write out surface related fields
    LOGICAL :: lwrite_extra               ! if .true., write out extra fields

    ! derived variables

    LOGICAL :: l_outputtime      ! if .true., output is written at the end of the time step.
    LOGICAL :: l_checkpoint_time ! if .true., restart file is written at the end of 
                                 ! the time step.
    LOGICAL :: l_diagtime        ! if .true., diagnostic output is computed and written
                                 ! at the end of the time step.               
                                                                                 
    LOGICAL :: lprepare_output(max_dom) ! if .true., save the prognostic variables
                                        ! to p_prog_out and update p_diag_out.
!  END TYPE t_io_config
  !>
  !!
 ! TYPE(t_io_config):: io_config(max_dom)


!CONTAINS

! SUBROUTINE setup_io_config !(n_dom)
!  END SUBROUTINE setup_io_config

END MODULE mo_io_config
