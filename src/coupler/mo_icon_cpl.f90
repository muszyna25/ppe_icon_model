!>
!! ICON coupler version 1.0.0
!!
!! The first version of the ICON coupler is built as
!! a prototype for coupling two ICON grids and an
!! arbitrary number of fields.
!! 
!! Short description of the User Interface
!! =======================================
!!
!!   1.) Initialisation and definition phase
!!
!!       Initialisation of the components
!!       --------------------------------
!!
!!       SUBROUTINE ICON_cpl_init ( comp_name, comp_id, ierror )
!!
!!         CHARACTER(len=*), INTENT(in) :: comp_name
!!         INTEGER, INTENT(out)         :: comp_id
!!         INTEGER, INTENT(out)         :: ierror
!!
!!       Grid Definition
!!       ---------------
!!
!!       SUBROUTINE ICON_cpl_def_grid ( comp_id, shape, grid_glob_index, grid_id, ierror )
!!
!!         INTEGER, INTENT(in)         :: comp_id
!!         INTEGER, INTENT(in)         :: shape(2)
!!         INTEGER, INTENT(in)         :: grid_glob_index(shape(1):shape(2))
!!         INTEGER, INTENT(out)        :: grid_id
!!         INTEGER, INTENT(out)        :: ierror
!!
!!
!!       Definition of coupling fields
!!       -----------------------------
!!
!!       SUBROUTINE ICON_cpl_def_field ( field_name, comp_id, grid_id, field_id, ierror )
!!
!!         CHARACTER(LEN=*), INTENT (in) :: field_name
!!         INTEGER, INTENT(in)           :: comp_id
!!         INTEGER, INTENT(in)           :: grid_id
!!         INTEGER, INTENT(out)          :: field_id
!!         INTEGER, INTENT(out)          :: ierror
!!
!!
!!   2.) Termination of the definition phase
!!
!!       SUBROUTINE ICON_cpl_search
!!
!!
!!   3.) Data exchange
!!
!!       Sending data
!!       ------------
!!
!!       SUBROUTINE ICON_cpl_put ( field_id, shape, send_field, info, ierror )
!!
!!         INTEGER, INTENT(in)    :: field_id 
!!         INTEGER, INTENT(in)    :: field_shape(3)
!!         INTEGER, INTENT(in)    :: info
!!         REAL (wp), INTENT(in)  :: send_field(field_shape(1):field_shape(2),field_shape(3))
!!         INTEGER, INTENT(out)   :: ierror
!!
!!       SUBROUTINE ICON_cpl_get ( field_id, field_shape, recv_field, ierror )
!!
!!         INTEGER, INTENT(in)    :: field_id
!!         INTEGER, INTENT(in)    :: shape(2)
!!         REAL (wp), INTENT(out) :: recv_field (field_shape(1):field_shape(2),field_shape(3))
!!         INTEGER, INTENT(out)   :: ierror
!!
!!
!!   4.) Termination of the coupling
!!
!!       SUBROUTINE ICON_cpl_finalize
!!
!!
!!  Important restrictions and requirements
!!  =======================================
!!
!!  ICON_cpl_init, ICON_cpl_search, and ICON_cpl_finalize are so called collective
!!  calls that have to be called by each ICON process irrespectively whether they
!!  will define a grid later or not. All processes that need to participate in the
!!  data exchange need to announce their local grid and fields. Currently it is
!!  assumed and required that all processes that announce their local grids will later
!!  on exchange the same lists of fields. 
!!
!! @author Rene Redler, Max-Planck Institute for Meteorology, Germany
!!
!! $Id:$
!!
!! @par Revision History
!! first implementation by Rene Redler (2010-02-13)
!!
!! @par Copyright
!! 2010 by MPI-M
!! This software is provided for non-commerncial use only.
!! See the LICENSE and WARRANTY conditions.
!!
!! @par License
!!
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! &ltol>
!! &ltli> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! &ltli> The code may not be re-distributed without the consent of the authors.
!! &ltli> The copyright notice and statement of authorship must appear in all
!!    copies.
!! &ltli> You accept the warranty conditions (see WARRANTY).
!! &ltli> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!!
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_icon_cpl

#ifndef NOMPI
  USE mpi
#endif
  USE mo_kind, ONLY     : wp
  USE mo_datetime, ONLY : t_datetime

  IMPLICIT NONE
  
  PUBLIC

  ! Info for put action

  INTEGER, PARAMETER        :: cpl_field_none = 0
  INTEGER, PARAMETER        :: cpl_field_acc  = 1
  INTEGER, PARAMETER        :: cpl_field_avg  = 2

  ! Maximum allowed length for characters

  INTEGER, PARAMETER        :: maxchar = 132

  ! List of names for coupling fields

  INTEGER, PARAMETER        :: nbr_max_fields = 32
  CHARACTER(len=maxchar)    :: fieldname(nbr_max_fields)

  ! Setting for debug output

  LOGICAL, PARAMETER        :: l_debug     = .true.
  INTEGER, PARAMETER        :: debug_level = 2
  INTEGER, PARAMETER        :: cplout      = 10

  CHARACTER(len=maxchar)    :: filename

  ! Component specific values

  INTEGER                   :: comp_id
  INTEGER                   :: comp_comm
  INTEGER                   :: ierr

  ! MPI Communicator handling

  INTEGER, PARAMETER        :: ICON_Root = 0

  INTEGER                   :: ICON_comm
  INTEGER                   :: ICON_comm_active
  INTEGER                   :: ICON_comp_comm

  INTEGER                   :: ICON_global_size
  INTEGER                   :: ICON_global_rank
  INTEGER                   :: ICON_local_size
  INTEGER                   :: ICON_local_rank

  ! test for cpl own MPI initialisation,
  ! true if ICOM cpl did the MPI initialisation

  LOGICAL                   :: l_MPI_was_initialized

  ! Component and Grid description

  TYPE t_comp
     CHARACTER(len=maxchar) :: comp_name
     CHARACTER(len=maxchar) :: nml_name
     INTEGER                :: comp_process
     LOGICAL                :: l_comp_status
     INTEGER                :: min_rank
     INTEGER                :: max_rank
     INTEGER                :: inc_rank
     LOGICAL                :: L_TIME_AVERAGE
     LOGICAL                :: L_TIME_ACCUMULATION
     INTEGER                :: COUPLING_FREQ
     INTEGER                :: TIME_STEP
  END TYPE t_comp

!>
!! Component Type
!!
!! Data type to store component information containing
!!
!! - Component (short)name
!! - Status of the component: false || true
!!

  TYPE t_grid
     INTEGER                :: comp_id
     INTEGER                :: grid_shape(2)
     LOGICAL                :: l_grid_status
     INTEGER, POINTER       :: grid_glob_index(:) => NULL()
  END TYPE t_grid

!>
!! Grid Type
!!
!! Data type to store grid information containing
!!
!! - Component ID
!! - Shape of array grid_global_index
!! - Array with global gird indices
!! - Status of the grid: false || true
!!

  TYPE t_coupling
     INTEGER                :: time_operation
     INTEGER                :: coupling_freq
     INTEGER                :: time_step
  END TYPE t_coupling

!>
!! Coupling Type
!!
!! Data type to store information about coupling
!!
!! - time_operation (average / accumate / none )
!! - couping_freq coupling frequency in seconds
!!

  TYPE t_field
     INTEGER                :: comp_id
     INTEGER                :: grid_id
     INTEGER                :: global_field_id
     INTEGER                :: event_id
     INTEGER                :: lag
     INTEGER                :: field_shape(3)
     INTEGER                :: accumulation_count
     LOGICAL                :: l_field_status
     CHARACTER(len=maxchar) :: field_name
     REAL(wp), POINTER      :: send_field_acc(:,:) => NULL()
     TYPE(t_coupling)       :: coupling
  END TYPE t_field

!>
!! Field Type
!!
!! Data type to store field information containing
!!
!! - Component ID
!! - Grid ID
!! - Global field ID
!! - Status of the field: false || true
!! - Field (short)name
!!
  TYPE (t_comp),  POINTER    :: comps (:)   => NULL()
  TYPE (t_grid),  POINTER    :: grids (:)   => NULL()
  TYPE (t_field), POINTER    :: fields(:)   => NULL()

  TYPE(t_datetime)           :: initial_date
  TYPE(t_datetime)           :: final_date

  ! Number of active components and grids

  INTEGER                   :: nbr_active_comps
  INTEGER                   :: nbr_active_grids
  INTEGER                   :: nbr_active_fields

  ! Size of increment and allocated data structures

  INTEGER, PARAMETER        :: nbr_ICON_inc = 8

  INTEGER                   :: nbr_ICON_grids
  INTEGER                   :: nbr_ICON_fields
  INTEGER                   :: nbr_ICON_couplings

  ! Marker for the physical component on a process

  INTEGER, PARAMETER        :: nbr_ICON_comps    = 4

  TYPE (t_comp)             :: complist (nbr_ICON_comps)

  ! General MPI function arguments

  INTEGER, PARAMETER        :: initag = 100 ! base tag to mark header
  INTEGER                   :: index        ! returned index from MPI_Waitany
  INTEGER                   :: datatype     ! MPI data type
#ifndef NOMPI
  INTEGER                   :: rstatus(MPI_STATUS_SIZE) ! MPI_Irecv status
  INTEGER                   :: wstatus(MPI_STATUS_SIZE) ! MPI_Wait status
#endif

  ! ===================================================================
  !
  ! PSMILe bsend conversion to MPI datatypes, not really needed
  ! used in psmile_def_datatypes, called from psmile_bsend_init
  ! to initialize the psmile_bsend, a convenient alternative to
  ! MPI_Bsend.

  INTEGER, PARAMETER   :: PRISM_Character        = 1
  INTEGER, PARAMETER   :: PRISM_Integer          = 2
  INTEGER, PARAMETER   :: PRISM_Logical          = 3
  INTEGER, PARAMETER   :: PRISM_Real             = 4
  INTEGER, PARAMETER   :: PRISM_Double_precision = 5
  INTEGER, PARAMETER   :: PRISM_Complex          = 6
  INTEGER, PARAMETER   :: PRISM_Double_complex   = 7
  !
  !   INTEGER, PARAMETER   :: PRISM_Quad_Precision   = 8
  !   INTEGER, PARAMETER   :: PRISM_Double_Quad      = 9
  !
  TYPE t_target_struct
     INTEGER                    :: offset                     ! offset to be considered
     INTEGER                    :: target_rank                ! target partner rank
     INTEGER                    :: source_rank                ! target partner rank
     INTEGER                    :: target_list_len            ! size of target_list
     INTEGER                    :: source_list_len            ! size of source_list
     INTEGER, POINTER           :: target_list(:)  => NULL()  ! subset of target indices
     INTEGER, POINTER           :: source_list(:)  => NULL()  ! found target locations 
  END TYPE t_target_struct

!>
!! Detailed target grid information
!!
!! Data type to store information from the target grid
!! that is collected during the search. For each active
!! relation between two processes this information is
!! instantiated and kept on the target side during the
!! search. This information is later used during the
!! data exchange to extract relevant information from
!! the data that are passed to the coupler. 
!! 
!! - Rank in the global communicator of the target
!! - Length of array target_list
!! - Target list containing the global indices of all
!! - Target points that have a common intersection with
!!   a particular source process
!! - Status of the field: false || true
!! - Rank in the global communicator of the source
!! - Length of array source_list
!! - Source list containing the global indices of all
!!   source points that have a common intersection with
!!   this target process
!!

  TYPE t_source_struct
     INTEGER                    :: target_rank              ! soure partner rank 
     INTEGER                    :: source_list_len          ! size of source_list
     INTEGER, POINTER           :: source_list(:) => NULL() ! found source locations for source
     INTEGER, POINTER           :: target_list(:) => NULL() ! found source locations for target
  END TYPE t_source_struct

!>
!! Detailed source grid information
!!
!! Data type to store information from the source grid
!! and on the source side that is collected during the
!! search. For each active relation between two processes
!! this information is instantiated and kept on the source
!! side during the search. This information is later used
!! during the data exchange to extract relevant information
!! from the data that are passed to the coupler. 
!! 
!! - Rank in the global communicator of the target
!! - Length of array source_list
!! - Source list containing the global indices of all
!!   source points that have a common intersection with
!!   this target process
!!

  TYPE (t_target_struct), POINTER :: target_locs(:) => NULL()
  TYPE (t_source_struct), POINTER :: source_locs(:) => NULL()

  ! Length of header messages

  INTEGER, PARAMETER            :: msg_len = 3

  ! Store and grid extents

  INTEGER, Allocatable          :: all_extents(:,:)

  PRIVATE :: comp_comm

  CONTAINS

    INTEGER FUNCTION set_cpl_local_comm ( comm )

      INTEGER, INTENT (IN) :: comm

      comp_comm = comm

      set_cpl_local_comm = 0

    END FUNCTION set_cpl_local_comm

    INTEGER FUNCTION get_cpl_local_comm ()

#ifndef NOMPI
      IF ( comp_comm /= MPI_COMM_NULL ) THEN
         get_cpl_local_comm = comp_comm
      ELSE
         get_cpl_local_comm = 0

         ! Rene:
         ! TODO: Should we stop here and return an error message
         !       or ignore this case and just provide MPI_COMM_NULL?

      ENDIF
#else
      get_cpl_local_comm = 0
#endif

    END FUNCTION get_cpl_local_comm

END MODULE mo_icon_cpl
