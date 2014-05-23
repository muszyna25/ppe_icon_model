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
!! Modified by Daniel Rieger, KIT (2014-05-22)
!! - Usage of the generalized ART infrastructure
!! @par Copyright
!! 2002-2010 by DWD, MPI-M, and KIT.
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
!!    an according license agreement with DWD, MPI-M, and KIT.
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
  USE mo_exception,                     ONLY: message
  USE mo_linked_list,                   ONLY: t_var_list
  USE mo_fortran_tools,                 ONLY: t_ptr_2d3d
  USE mo_advection_config,              ONLY: t_advection_config
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_tend
  USE mo_nonhydro_types,                ONLY: t_nh_prog
  USE mo_art_config,                    ONLY: art_config
#ifdef __ICON_ART
  USE mo_art_tracer,                    ONLY: art_tracer
  USE mo_art_init,                      ONLY: art_init
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_tracer_interface 

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_tracer_interface(defcase,jg,nblks_c,this_list,vname_prefix,&
   &                            ptr_arr,advconf,phy_tend, p_prog,&
   &                            timelev,ldims,tlev_source)
!! Interface for ART-routine art_ini_tracer 
!!
!! This interface calls the ART-routine art_ini_tracer, if ICON has been 
!! built including the ART-package. Otherwise, this is simply a dummy 
!! routine.
!!
!! @par Revision History
!! Initial revision by Kristina Lundgren, KIT (2012-04-03)
  INTEGER,INTENT(in)             :: &
    &   jg,                         & !< patch id
    &   nblks_c                       !< patch block
  TYPE(t_var_list),INTENT(INOUT) :: &
    &   this_list                     !< current list: prognostic or phys. tend.    
  TYPE(t_ptr_2d3d),INTENT(inout) :: &
    &   ptr_arr(:)                    !< pointer to each element in list
  TYPE(t_advection_config),INTENT(inout) ::  &
    &   advconf                       !< advection config
      
  TYPE(t_nwp_phy_tend),INTENT(inout),OPTIONAL :: &
    &   phy_tend                      !< physical tendencies
      
  TYPE(t_nh_prog),INTENT(inout),OPTIONAL :: &
    &   p_prog                        !< prognostic variables
      
  INTEGER,INTENT(in), OPTIONAL   :: &
    &   timelev,                    & !< drieg : why is timelevel optional?
    &   ldims(3),                   & !< local dimensions, for checking
    &   tlev_source                   !< actual TL for TL dependent vars

  CHARACTER(len=*), INTENT(IN)   :: & 
    &   vname_prefix                  !< list name
  CHARACTER(len=*), INTENT(IN)   :: & 
    &   defcase                       !< definition of case 

  !-----------------------------------------------------------------------
 
#ifdef __ICON_ART
  IF (art_config(jg)%lart) THEN
    CALL message('','ART: Definition of tracers for defcase: '//TRIM(defcase))
      
    CALL art_tracer(defcase,jg,nblks_c,this_list,vname_prefix,ptr_arr,advconf,phy_tend,p_prog,timelev,ldims,&
      & tlev_source) 
      
    IF (TRIM(defcase) .EQ. 'prog' .AND. timelev .EQ. 1) THEN 
      CALL art_init(jg,this_list)
    END IF
  ENDIF ! lart
#endif

END SUBROUTINE art_tracer_interface 
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_tracer_interface

