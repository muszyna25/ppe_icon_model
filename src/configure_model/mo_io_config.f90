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
  USE mo_run_config,     ONLY: dtime
  USE mo_io_units,       ONLY: filename_max

  IMPLICIT NONE
  PUBLIC
  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !--------------------------------------------------------------------------
  ! Derived type 
  !--------------------------------------------------------------------------

  ! from namelist
  
  LOGICAL :: lkeep_in_sync              ! if .true., sync stream after each timestep
  REAL(wp):: dt_diag                    ! diagnostic output timestep [seconds]
  REAL(wp):: dt_checkpoint              ! timestep [seconds] for triggering new restart file

  INTEGER :: inextra_2d                 ! number of extra output fields for debugging
  INTEGER :: inextra_3d                 ! number of extra output fields for debugging
  LOGICAL :: lflux_avg                  ! if .FALSE. the output fluxes are accumulated 
  !  from the beginning of the run
  ! if .TRUE. the output fluxex are average values 
  !  from the beginning of the run, except of 
  !  TOT_PREC that would be accumulated
  INTEGER :: itype_pres_msl             ! Specifies method for computation of mean sea level pressure
  INTEGER :: itype_rh                   ! Specifies method for computation of relative humidity

  CHARACTER(LEN=filename_max) :: &
    &        output_nml_dict,    &      !< maps variable names onto the internal ICON names.
    &        netcdf_dict                !< maps internal variable names onto names in output file (NetCDF only).

  LOGICAL :: lzaxis_reference           !< use ZAXIS_REFERENCE instead of ZAXIS_HYBRID for atmospheric 
                                        !  output fields

  ! derived variables

  LOGICAL :: l_outputtime      ! if .true., output is written at the end of the time step.
  LOGICAL :: l_diagtime        ! if .true., diagnostic output is computed and written at the end of the time step.
!  LOGICAL ::  use_set_event_to_simstep = .true.       ! if .true., the set_event_to_simstep routine is activated

  INTEGER, PARAMETER :: read_netcdf_broadcast_method  = 1
  INTEGER, PARAMETER :: read_netcdf_distribute_method = 2
  INTEGER :: default_read_method = read_netcdf_broadcast_method


CONTAINS
  !----------------------------------------------------------------------------------
   FUNCTION n_checkpoints()

     INTEGER :: n_checkpoints

     n_checkpoints = NINT(dt_checkpoint/dtime)  ! write restart files
   END FUNCTION
  !----------------------------------------------------------------------------------
   
  !----------------------------------------------------------------------------------
   FUNCTION n_diags()

     INTEGER :: n_diags

     n_diags  = MAX(1,NINT(dt_diag/dtime)) ! number of: diagnose of total integrals
   END FUNCTION
  !----------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------

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
  !----------------------------------------------------------------------------------
  

END MODULE mo_io_config
