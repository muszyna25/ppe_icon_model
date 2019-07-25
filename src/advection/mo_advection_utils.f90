!>
!! Some utilities which are specific to the transport algorithm.
!!
!! Module contains some functions and procedures which are specifically related
!! to the transport schemes. These subroutines or functions are needed at
!! various places within the transport scheme. Therefore outsourcing these
!! routines protects from possible circular dependencies.
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2010-03-04)
!! Modification by Daniel Reinert, DWD (2010-04-23)
!! - implemented generalized Lax-Friedrich flux function
!!   laxfr_upflux_v, which allows to use the same transport
!!   code for pressure and height based vertical coordinate
!!   systems.
!! Modification by Daniel Reinert, DWD (2016-01-12)
!! - removed obsolete routine tupdate_tracer
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_advection_utils

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: message, message_text
  USE mo_fortran_tools,         ONLY: t_ptr_2d3d
  USE mo_cf_convention,         ONLY: t_cf_var
  USE mo_grib2,                 ONLY: t_grib2_var
  USE mo_var_list,              ONLY: add_ref, find_list_element
  USE mo_linked_list,           ONLY: t_var_list, t_list_element
  USE mo_var_metadata_types,    ONLY: t_var_metadata,                    &
    &                                 t_vert_interp_meta,                &
    &                                 t_hor_interp_meta,                 &
    &                                 MAX_GROUPS,                        &
    &                                 t_post_op_meta
  USE mo_tracer_metadata_types, ONLY: t_tracer_meta

  USE mo_advection_config,      ONLY: t_advection_config

  IMPLICIT NONE

  PRIVATE


  ! functions
  PUBLIC :: laxfr_upflux
  PUBLIC :: laxfr_upflux_v

  ! subroutines
  PUBLIC :: add_tracer_ref            ! add new tracer component

  ! types
  PUBLIC :: t_list2D

  INTERFACE add_tracer_ref
    MODULE PROCEDURE add_var_list_reference_tracer
  END INTERFACE add_tracer_ref

  TYPE t_list2D
    INTEGER, POINTER :: eidx(:,:)
    INTEGER, POINTER :: elev(:,:)
    INTEGER, POINTER :: len(:)
    INTEGER          :: npoints
  END TYPE t_list2D


CONTAINS

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Lax Friedrichs first order upwind flux,.
  !!
  !! Lax Friedrichs first order upwind flux,
  !! used in conservative advection routines.
  !! For passive advection, equivalent to
  !! any other first order upwind flux.
  !!
  !! @par Revision History
  !! Developed  by L.Bonaventura  (2004).
  !!
  FUNCTION laxfr_upflux( p_vn, p_psi1, p_psi2 )  RESULT(p_upflux)
    !

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: p_vn
    REAL(wp), INTENT(in) :: p_psi1, p_psi2

    REAL(wp) :: p_upflux

    !$ACC ROUTINE SEQ

    !-----------------------------------------------------------------------
    p_upflux = 0.5_wp * (        p_vn  *( p_psi1 + p_psi2 )    &
      &                   - ABS( p_vn )*( p_psi2 - p_psi1 ) )

  END FUNCTION laxfr_upflux



  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Lax Friedrichs first order upwind flux for vertical advection,.
  !!
  !! Generalized Lax Friedrichs first order upwind flux,
  !! used in conservative vertical advection routines.
  !! For passive advection, equivalent to any other first
  !! order upwind flux.
  !! Applicable to height based vertical coordinate systems. 
  !!
  !! @par Revision History
  !! Developed  by L.Bonaventura  (2004).
  !! Modification by Daniel Reinert, DWD (2010-04-23)
  !! - generalized for p- and z-based vertical coordinate systems
  !! Modification by Daniel Reinert, DWD (2019-02-18)
  !! - revert generalization for p- and z-based vertical coordinate systems.
  !!   From now on only valid for z-based coordinate systems.
  !!
  FUNCTION laxfr_upflux_v( p_w, p_psi1, p_psi2 )  RESULT(p_upflux)
    !

    IMPLICIT NONE

    REAL(wp), INTENT(in) :: p_w
    REAL(wp), INTENT(in) :: p_psi1, p_psi2

    REAL(wp) :: p_upflux

    !$ACC ROUTINE SEQ

    !-----------------------------------------------------------------------
    p_upflux = 0.5_wp * (        p_w  *( p_psi1 + p_psi2 )    &
      &                   + ABS( p_w )*( p_psi2 - p_psi1 ) )

  END FUNCTION laxfr_upflux_v


  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! reference to an existing pointer to a 3D tracer field
  ! optionally overwrite some default meta data
  !
  SUBROUTINE add_var_list_reference_tracer(this_list, target_name, tracer_name,    &
    &        tracer_idx, ptr_arr, cf, grib2, advconf,jg, ldims, loutput, lrestart, &
    &        isteptype, tlev_source, vert_interp, hor_interp, in_group, post_op,   &
    &        tracer_info)

    TYPE(t_var_list)    , INTENT(inout)        :: this_list
    CHARACTER(len=*)    , INTENT(in)           :: target_name
    CHARACTER(len=*)    , INTENT(in)           :: tracer_name
    INTEGER             , INTENT(inout)        :: tracer_idx       ! index in 4D tracer container
    TYPE(t_ptr_2d3d)    , INTENT(inout)        :: ptr_arr(:)
    TYPE(t_cf_var)      , INTENT(in)           :: cf               ! CF related metadata
    TYPE(t_grib2_var)   , INTENT(in)           :: grib2            ! GRIB2 related metadata
    TYPE(t_advection_config), INTENT(inout)    :: advconf          ! adv configure state
    INTEGER             , INTENT(in), OPTIONAL :: jg               ! patch id
    INTEGER             , INTENT(in), OPTIONAL :: ldims(3)         ! local dimensions, for checking
    LOGICAL             , INTENT(in), OPTIONAL :: loutput          ! output flag
    LOGICAL             , INTENT(in), OPTIONAL :: lrestart         ! restart flag
    INTEGER,              INTENT(in), OPTIONAL :: isteptype        ! type of statistical processing
    INTEGER             , INTENT(in), OPTIONAL :: tlev_source      ! actual TL for TL dependent vars
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp   ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp    ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(MAX_GROUPS)          ! groups to which a variable belongs
    TYPE(t_post_op_meta), INTENT(in), OPTIONAL :: post_op          ! post operation (e.g. scale with const. factor or rho)
    CLASS(t_tracer_meta), INTENT(in), OPTIONAL :: tracer_info      ! tracer meta data

    ! Local variables:
    TYPE(t_list_element), POINTER :: target_element
    TYPE(t_var_metadata), POINTER :: target_info

    INTEGER :: zihadv_tracer, zivadv_tracer

    CHARACTER(*), PARAMETER :: routine = "add_tracer_ref"

  !------------------------------------------------------------------

    ! get pointer to target element (in this case 4D tracer container)
    target_element => find_list_element (this_list, target_name)
    ! get tracer field metadata
    target_info => target_element%field%info

    ! get index of current field in 4D container and set
    ! tracer index accordingly.
    ! Note that this will be repeated for each patch. Otherwise it may happen that
    ! some MPI-processes miss this assignment.
    !
    tracer_idx = target_info%ncontained+1  ! index in 4D tracer container

    WRITE (message_text,'(a,i3,a,a)')                                            &
      & "tracer index ", tracer_idx," assigned to tracer field ", TRIM(tracer_name)
    CALL message(TRIM(routine),message_text)



    ! set default values
    ! If ihadv_tracer or ivadv_tracer are not present, take respective value from the
    ! advection_config state.
    ! If ihadv_tracer or ivadv_tracer are present, take those values, and overwrite
    ! the respective values of the advection_config state.
    !
    IF( PRESENT(tracer_info) ) THEN
      IF ( tracer_info%ihadv_tracer /= 2) THEN
        zihadv_tracer = tracer_info%ihadv_tracer
        ! BE AWARE, THAT ihadv_tracer IS NOT SANITY-CHECKED. THIS OVERWRITES THE
        ! SANITY CHECKED NAMELIST SETTINGS.
        advconf%ihadv_tracer(tracer_idx) = tracer_info%ihadv_tracer
      ELSE
        zihadv_tracer = advconf%ihadv_tracer(tracer_idx)
      ENDIF
    ENDIF

    IF( PRESENT(tracer_info) ) THEN
      IF ( tracer_info%ivadv_tracer /= 3) THEN
        zivadv_tracer = tracer_info%ivadv_tracer
        ! BE AWARE, THAT ivadv_tracer IS NOT SANITY-CHECKED. THIS OVERWRITES THE
        ! SANITY CHECKED NAMELIST SETTINGS.
        advconf%ivadv_tracer(tracer_idx) = tracer_info%ivadv_tracer
      ELSE
        zivadv_tracer = advconf%ivadv_tracer(tracer_idx)
      ENDIF
    ENDIF

    ! create new table entry reference including additional tracer metadata
    CALL add_ref( this_list, target_name, tracer_name, ptr_arr(tracer_idx)%p_3d, &
       &          target_info%hgrid, target_info%vgrid, cf, grib2,               &
       &          ldims=ldims, loutput=loutput, lrestart=lrestart,               &
       &          isteptype=isteptype, tlev_source=tlev_source,                  &
       &          vert_interp=vert_interp, hor_interp=hor_interp,                &
       &          tracer_info=tracer_info, in_group=in_group, post_op=post_op)

  END SUBROUTINE add_var_list_reference_tracer


END MODULE mo_advection_utils
