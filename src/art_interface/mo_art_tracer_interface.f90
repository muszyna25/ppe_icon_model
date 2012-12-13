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
    USE mo_fortran_tools,        ONLY: t_ptr_2d3d
    USE mo_advection_config,     ONLY: t_advection_config
    USE mo_nwp_phy_types,        ONLY: t_nwp_phy_tend
    USE mo_nonhydro_types,       ONLY: t_nh_prog
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
  SUBROUTINE art_tracer_interface(defcase,jg,nblks_c,this_list,vname_prefix,&
    &                            ptr_arr,advconf,phy_tend, p_prog,&
    &                            timelev,ldims,tlev_source)

    INTEGER                 ,INTENT(in) :: jg                    !< patch id
    INTEGER                 ,INTENT(in) :: nblks_c               !< patch block
    TYPE(t_var_list)        ,INTENT(INOUT) :: this_list          !< current list: prognostic or phys. tend.
    TYPE(t_ptr_2d3d)        ,INTENT(inout) :: ptr_arr(:)         !< pointer to each element in list
    TYPE(t_advection_config),INTENT(inout) :: advconf            !< advection config
    TYPE(t_nwp_phy_tend)    ,INTENT(inout),OPTIONAL :: phy_tend 
    TYPE(t_nh_prog)         ,INTENT(inout),OPTIONAL :: p_prog 
    INTEGER                 ,INTENT(in), OPTIONAL :: timelev 
    INTEGER                 ,INTENT(in), OPTIONAL :: ldims(3)            ! local dimensions, for checking
    INTEGER                 ,INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars

    CHARACTER(len=*), INTENT(IN)      :: & !< list name
      &  vname_prefix
    CHARACTER(len=*), INTENT(IN)      :: & !< definition of case 
      & defcase 

    !-----------------------------------------------------------------------
 
#ifdef __ICON_ART
   
      CALL art_tracer(defcase,jg,nblks_c,this_list,vname_prefix,ptr_arr,advconf,phy_tend,p_prog,timelev,ldims,&
       & tlev_source) 

#endif

  END SUBROUTINE art_tracer_interface 


END MODULE mo_art_tracer_interface

