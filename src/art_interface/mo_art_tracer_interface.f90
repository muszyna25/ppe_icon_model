!>
!! Provides interface to the ART-routine for initializing ART-tracer compounds to the tracer field. 
!!
!! This module provides an interface to the ART-routine art_ini_tracer.
!! The interface is written in such a way, that ICON will compile and run 
!! properly, even if the ART-routines are not available at compile time.
!!
!!
!! @author Kristina Lundgren, KIT
!!
!! @par Revision History
!! Initial revision by Kristina Lundgren (2012-04-03)
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
MODULE mo_art_tracer_interface
    USE mo_model_domain,         ONLY: t_patch
    USE mo_linked_list,          ONLY: t_var_list
    USE mo_nonhydro_types,       ONLY:t_ptr_nh


#ifdef __ICON_ART
    USE mo_art_tracer,       ONLY:art_tracer
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_tracer_interface 


CONTAINS


  !>
  !! Interface for ART-routine art_ini_tracer 
  !!
  !! This interface calls the ART-routine art_ini_tracer, if ICON has been 
  !! built including the ART-package. Otherwise, this is simply a dummy 
  !! routine.
  !!
  !! @par Revision History
  !! Initial revision by Kristina Lundgren, KIT (2012-04-03)
  SUBROUTINE art_tracer_interface(p_patch,p_prog_list,vname_prefix,ptr_arr, &
    &                             timelev,ldims,tlev_source)

!    TYPE(t_var_list), INTENT(INOUT)   :: this_list !< current prognostic state list 
    TYPE(t_patch), TARGET, INTENT(IN) :: & !< current patch
      &  p_patch
    TYPE(t_var_list), INTENT(INOUT)   :: p_prog_list!< current prognostic state list 
    TYPE(t_ptr_nh)      , INTENT(inout)        :: ptr_arr(:)
    INTEGER             , INTENT(in), OPTIONAL :: ldims(3)            ! local dimensions, for checking
    INTEGER             , INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    INTEGER, INTENT(IN) :: timelev

    CHARACTER(len=*), INTENT(IN)      :: & !< list name
      &  vname_prefix
    !-----------------------------------------------------------------------
 
#ifdef __ICON_ART
   
      CALL art_tracer(p_patch,p_prog_list,vname_prefix,ptr_arr,timelev,ldims,tlev_source) 


#endif

  END SUBROUTINE art_tracer_interface 


END MODULE mo_art_tracer_interface

