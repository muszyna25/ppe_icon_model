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
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_art_tracer_interface
  USE mo_exception,                     ONLY: message,finish
  USE mo_linked_list,                   ONLY: t_var_list
  USE mo_fortran_tools,                 ONLY: t_ptr_2d3d
  USE mo_advection_config,              ONLY: t_advection_config
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_tend
  USE mo_nonhydro_types,                ONLY: t_nh_prog
  USE mo_run_config,                    ONLY: lart, iforcing
  USE mo_impl_constants,                ONLY: iecham
  USE mo_radiation_config,              ONLY: irad_o3
  USE mo_echam_rad_config,              ONLY: echam_rad_config
  USE mo_art_config,                    ONLY: art_config
  USE mo_time_config,                   ONLY: time_config
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art_tracInt
#ifdef __ICON_ART
  USE mo_art_tracer,                    ONLY: art_tracer
  USE mo_art_init,                      ONLY: art_init
  USE mo_art_diag_state,                ONLY: art_create_diagnostics
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
  TYPE(t_ptr_2d3d),INTENT(inout),OPTIONAL         :: &
    &   ptr_arr(:)                    !< pointer to each element in list
  TYPE(t_advection_config),INTENT(inout),OPTIONAL ::  &
    &   advconf                       !< advection config
      
  TYPE(t_nwp_phy_tend),INTENT(inout),OPTIONAL     :: &
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
  IF (lart) THEN
    IF (timers_level > 3) CALL timer_start(timer_art_tracInt)

    IF (TRIM(defcase) .NE. 'diag') THEN
    
      IF (.NOT. PRESENT(advconf)) THEN
              CALL finish('mo_art_tracer_interface:art_tracer_interface',  &
         &      'advconf not present at defcase'//TRIM(defcase))
      ENDIF
      IF (.NOT. PRESENT(ptr_arr)) THEN
              CALL finish('mo_art_tracer_interface:art_tracer_interface',  &
         &      'ptr_arr not present at defcase'//TRIM(defcase))
      ENDIF
    
      CALL message('','ART: Definition of tracers for defcase: '//TRIM(defcase))
        
      IF (TRIM(defcase) .EQ. 'prog') THEN 
        CALL art_tracer(defcase,jg,nblks_c,this_list,vname_prefix,ptr_arr,advconf,p_prog=p_prog,timelev=timelev,  &
          & ldims=ldims, tlev_source=tlev_source)
      ELSE
        CALL art_tracer(defcase,jg,nblks_c,this_list,vname_prefix,ptr_arr,advconf,phy_tend=phy_tend,              &
          & ldims=ldims, tlev_source=tlev_source)
      ENDIF
        
      IF (TRIM(defcase) .EQ. 'prog' .AND. timelev .EQ. 1) THEN 

        IF ( iforcing == iecham) THEN
          irad_o3 = echam_rad_config(jg)%irad_o3
        END IF

        CALL art_init(jg, time_config%tc_dt_model, time_config%tc_exp_refdate, &
          &           this_list,tracer=p_prog%tracer)
      ENDIF
    ELSE !defcase is diag
      CALL art_create_diagnostics(jg, this_list)
    ENDIF

    IF (timers_level > 3) CALL timer_stop(timer_art_tracInt)
  ENDIF ! lart
#endif

END SUBROUTINE art_tracer_interface 
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_tracer_interface

