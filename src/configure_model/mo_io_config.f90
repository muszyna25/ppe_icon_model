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
  USE mo_run_config,     ONLY: dtime, nsteps

  IMPLICIT NONE
  PUBLIC
  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !--------------------------------------------------------------------------
  ! Derived type 
  !--------------------------------------------------------------------------
  !TYPE :: t_io_config

    ! from namelist

    CHARACTER(len=max_char_length) :: out_expname
    INTEGER :: out_filetype               ! 1 - GRIB1, 2 - netCDF
    LOGICAL :: lkeep_in_sync              ! if .true., sync stream after each timestep
    REAL(wp):: dt_data                    ! output timestep [seconds]
    REAL(wp):: dt_diag                    ! diagnostic output timestep [seconds]
    REAL(wp):: dt_file                    ! timestep [seconds] for triggering new output file
    REAL(wp):: dt_checkpoint              ! timestep [seconds] for triggering new restart file

    LOGICAL :: lwrite_initial             ! if .true., write out initial values
    LOGICAL :: lwrite_dblprec             ! if .true. create double precision output
    LOGICAL :: lwrite_decomposition       ! if .true. write field with MPI_RANK

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
    LOGICAL :: lwrite_pzlev               ! if .true. extra output on p- and/or z-levels
    INTEGER :: inextra_2d                 ! number of extra output fields for debugging
    INTEGER :: inextra_3d                 ! number of extra output fields for debugging
    LOGICAL :: lflux_avg                  ! if .FALSE. the output fluxes are accumulated 
                                          !  from the beginning of the run
                                          ! if .TRUE. the output fluxex are average values 
                                          !  from the beginning of the run, except of 
                                          !  TOT_PREC that would be accumulated
    LOGICAL :: lwrite_oce_timestepping    ! if .true. write intermediate ocean variables


    ! derived variables

    LOGICAL :: no_output         ! if .true., no output is written 
    LOGICAL :: l_outputtime      ! if .true., output is written at the end of the time step.
    LOGICAL :: l_diagtime        ! if .true., diagnostic output is computed and written
                                 ! at the end of the time step.

    LOGICAL :: lprepare_output(max_dom) ! if .true., save the prognostic variables
                                        ! to p_prog_out and update p_diag_out.

  !END TYPE t_io_config
  !>
  !!
  !TYPE(t_io_config):: io_config(max_dom)

CONTAINS
!  !>
!  !!
!  subroutine configure_io(n_dom)
!  end subroutine configure_io

   FUNCTION istime4output(sim_time) RESULT(retval)
     REAL(wp), INTENT(IN)  :: sim_time            ! simulation time [s]
     LOGICAL               :: retval
     REAL(wp)              :: nearest_output_time ! nearest output time [s]


     ! get nearest output time
     nearest_output_time = NINT(sim_time/dt_data) * dt_data

     ! write output (true/false)
     retval =( ( (-dtime/2._wp) < (sim_time - nearest_output_time)) .AND.  &
       &       ( (sim_time - nearest_output_time) <= (dtime/2._wp)) )

   END FUNCTION istime4output

!!$   FUNCTION istime4output(current_timestep) RESULT(retval)
!!$     LOGICAL :: retval
!!$     INTEGER, INTENT(IN)  :: current_timestep
!!$
!!$     INTEGER :: n_io
!!$
!!$     n_io    = NINT(dt_data/dtime)        ! write output
!!$
!!$     IF ( (MOD(current_timestep-1,n_io)==0 .AND. current_timestep/=1) ) THEN
!!$       retval = .TRUE.
!!$     ELSE
!!$       retval = .FALSE.
!!$     END IF
!!$   END FUNCTION istime4output

   FUNCTION istime4newoutputfile(current_timestep) RESULT(retval)
     LOGICAL :: retval
     INTEGER, INTENT(IN) :: current_timestep

     INTEGER :: n_file

     n_file  = n_files()        ! trigger new output file

     IF (current_timestep/=1&
       & .AND.(MOD(current_timestep-1,n_file)==0)&
       & .AND.current_timestep/=nsteps) THEN
       retval = .TRUE.
     ELSE
       retval = .FALSE.
     END IF
   END FUNCTION istime4newoutputfile

   FUNCTION n_checkpoints()

     INTEGER :: n_checkpoints

     n_checkpoints = NINT(dt_checkpoint/dtime)  ! write restart files
   END FUNCTION
   FUNCTION n_files()

     INTEGER :: n_files
     n_files  = NINT(dt_file/dtime)        ! trigger new output file
   END FUNCTION
   FUNCTION n_ios()

     INTEGER :: n_ios

     n_ios    = NINT(dt_data/dtime)        ! number of: write output
   END FUNCTION
   FUNCTION n_diags()

     INTEGER :: n_diags

     n_diags  = MAX(1,NINT(dt_diag/dtime)) ! number of: diagnose of total integrals
   END FUNCTION

   FUNCTION is_checkpoint_time(current_step, n_checkpoints, n_steps) RESULT(l_checkpoint)
     INTEGER, INTENT(IN)            :: current_step, n_checkpoints
     INTEGER, INTENT(IN), OPTIONAL  :: n_steps

     LOGICAL              :: l_checkpoint

     IF (PRESENT(n_steps)) THEN
       IF ( MOD(current_step,n_checkpoints)==0 .AND. current_step/=n_steps ) THEN
         l_checkpoint = .TRUE.
       ELSE
         l_checkpoint = .FALSE.
       END IF
     ELSE
       IF ( MOD(current_step,n_checkpoints)==0 ) THEN
         l_checkpoint = .TRUE.
       ELSE
         l_checkpoint = .FALSE.
       END IF
     END IF
   END FUNCTION

END MODULE mo_io_config
