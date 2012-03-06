!>
!! This module contains utilty types
!!
!! @par Revision History
!!  Created by Leonidas Linardakis, MPIM (2012-03-06)
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
!!
MODULE mo_utility_types

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC :: t_subset_range, t_subset_range_index

  !----------------------------------------------------
  !> Defines a subset in a range (in terms of blocks)
  TYPE :: t_subset_range
    
    INTEGER :: start_block
    INTEGER :: start_index
    INTEGER :: end_block
    INTEGER :: end_index
    INTEGER :: block_size   
    
  END TYPE
  !----------------------------------------------------
  

  !----------------------------------------------------
  !> Defines an index for a subset_range 
  TYPE :: t_subset_range_index
  
    INTEGER :: current_block ! the current block in the subset
    INTEGER :: current_index ! the current index in the subset
    
    INTEGER :: current_start_index ! the current start index within the current block,
    INTEGER :: current_end_index   ! the current end index within the current block,

    TYPE(t_subset_range), POINTER :: t_subset_range
  
  END TYPE
  !----------------------------------------------------

END MODULE mo_utility_types

