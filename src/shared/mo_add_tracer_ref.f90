!>
!!               This module is not used
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
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
MODULE mo_add_tracer_ref

  USE mo_kind,             ONLY: wp, i8
  USE mo_exception,        ONLY: message, message_text, finish
  USE mo_fortran_tools,    ONLY: t_ptr_2d3d
  USE mo_cf_convention,    ONLY: t_cf_var
  USE mo_grib2,            ONLY: t_grib2_var
  USE mo_var_list,         ONLY: create_tracer_metadata, add_ref
  USE mo_linked_list,      ONLY: t_var_list, t_list_element,        &
       &                         new_list, delete_list,             &
       &                         append_list_element,               &
       &                         find_list_element,                 &
       &                         delete_list_element
  USE mo_var_metadata,     ONLY: t_var_metadata, t_union_vals,      &
    &                            t_tracer_meta,                     &
    &                            t_vert_interp_meta,                &
    &                            t_hor_interp_meta,                 &
    &                            VARNAME_LEN, VAR_GROUPS,           &
    &                            VINTP_TYPE_LIST,                   &
    &                            t_post_op_meta, POST_OP_NONE
  
  USE mo_advection_config, ONLY: t_advection_config
  USE mo_art_config,       ONLY: t_art_config, art_config

  IMPLICIT NONE

  PRIVATE
  

  PUBLIC :: add_tracer_ref            ! add new tracer component

  INTERFACE add_tracer_ref
    MODULE PROCEDURE add_var_list_reference_tracer
  END INTERFACE add_tracer_ref

  !
CONTAINS


  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! reference to an existing pointer to a 3D tracer field
  ! optionally overwrite some default meta data
  !
  SUBROUTINE add_var_list_reference_tracer(this_list, target_name, tracer_name,  &
    &        tracer_idx, ptr_arr, cf, grib2, advconf,jg, ldims, loutput, lrestart,  &
    &        isteptype, tlev_source, vert_interp, hor_interp, in_group,          &
    &        lis_tracer, ihadv_tracer, ivadv_tracer, lturb_tracer, lsed_tracer,  &
    &        ldep_tracer, lconv_tracer, lwash_tracer, rdiameter_tracer,          &
    &        rrho_tracer)

    TYPE(t_var_list)    , INTENT(inout)        :: this_list
    CHARACTER(len=*)    , INTENT(in)           :: target_name
    CHARACTER(len=*)    , INTENT(in)           :: tracer_name
    INTEGER             , INTENT(inout)        :: tracer_idx     ! index in 4D tracer container
    TYPE(t_ptr_2d3d)    , INTENT(inout)        :: ptr_arr(:)
    TYPE(t_cf_var)      , INTENT(in)           :: cf             ! CF related metadata
    TYPE(t_grib2_var)   , INTENT(in)           :: grib2          ! GRIB2 related metadata
    TYPE(t_advection_config), INTENT(inout)    :: advconf        ! adv configure state
    INTEGER             , INTENT(in), OPTIONAL :: jg             ! patch id
    INTEGER             , INTENT(in), OPTIONAL :: ldims(3)       ! local dimensions, for checking
    LOGICAL             , INTENT(in), OPTIONAL :: loutput        ! output flag
    LOGICAL             , INTENT(in), OPTIONAL :: lrestart       ! restart flag
    INTEGER,              INTENT(in), OPTIONAL :: isteptype      ! type of statistical processing
    INTEGER             , INTENT(in), OPTIONAL :: tlev_source    ! actual TL for TL dependent vars
    TYPE(t_vert_interp_meta),INTENT(in), OPTIONAL :: vert_interp ! vertical interpolation metadata
    TYPE(t_hor_interp_meta), INTENT(in), OPTIONAL :: hor_interp  ! horizontal interpolation metadata
    LOGICAL, INTENT(in), OPTIONAL :: in_group(SIZE(VAR_GROUPS))  ! groups to which a variable belongs
    LOGICAL             , INTENT(in), OPTIONAL :: lis_tracer     ! this is a tracer field (TRUE/FALSE)
    INTEGER             , INTENT(in), OPTIONAL :: ihadv_tracer   ! method for hor. transport
    INTEGER             , INTENT(in), OPTIONAL :: ivadv_tracer   ! method for vert. transport
    LOGICAL             , INTENT(in), OPTIONAL :: lturb_tracer   ! turbulent transport (TRUE/FALSE)
    LOGICAL             , INTENT(in), OPTIONAL :: lsed_tracer    ! sedimentation (TRUE/FALSE)
    LOGICAL             , INTENT(in), OPTIONAL :: ldep_tracer    ! dry deposition (TRUE/FALSE)
    LOGICAL             , INTENT(in), OPTIONAL :: lconv_tracer   ! convection  (TRUE/FALSE)
    LOGICAL             , INTENT(in), OPTIONAL :: lwash_tracer   ! washout (TRUE/FALSE)
    REAL(wp)            , INTENT(in), OPTIONAL :: rdiameter_tracer ! particle diameter in m
    REAL(wp)            , INTENT(in), OPTIONAL :: rrho_tracer    ! particle density in kg m^-3


    ! Local variables:
    TYPE(t_list_element), POINTER :: target_element
    TYPE(t_var_metadata), POINTER :: target_info
    TYPE(t_tracer_meta)           :: tracer_info    ! tracer meta data

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
    IF ( PRESENT( ihadv_tracer )) THEN
      zihadv_tracer = ihadv_tracer
      ! BE AWARE, THAT ivadv_tracer IS NOT SANITY-CHECKED. THIS OVERWRITES THE
      ! SANITY CHECKED NAMELIST SETTINGS.
      advconf%ihadv_tracer(tracer_idx) = ihadv_tracer
    ELSE
      zihadv_tracer = advconf%ihadv_tracer(tracer_idx)
    ENDIF

    IF ( PRESENT( ivadv_tracer )) THEN
      zivadv_tracer = ivadv_tracer
      ! BE AWARE, THAT ivadv_tracer IS NOT SANITY-CHECKED. THIS OVERWRITES THE
      ! SANITY CHECKED NAMELIST SETTINGS.
      advconf%ivadv_tracer(tracer_idx) = ivadv_tracer
    ELSE
      zivadv_tracer = advconf%ivadv_tracer(tracer_idx)
    ENDIF

    ! generate tracer metadata information
    !
    tracer_info = create_tracer_metadata(                                     &
      &                                  lis_tracer       = lis_tracer,       &
      &                                  ihadv_tracer     = zihadv_tracer,    &
      &                                  ivadv_tracer     = zivadv_tracer,    &
      &                                  lturb_tracer     = lturb_tracer,     &
      &                                  lsed_tracer      = lsed_tracer,      &
      &                                  ldep_tracer      = ldep_tracer,      &
      &                                  lconv_tracer     = lconv_tracer,     &
      &                                  lwash_tracer     = lwash_tracer,     &
      &                                  rdiameter_tracer = rdiameter_tracer, &
      &                                  rrho_tracer      = rrho_tracer )



    ! create new table entry reference including additional tracer metadata
    CALL add_ref( this_list, target_name, tracer_name, ptr_arr(tracer_idx)%p_3d, &
       &          target_info%hgrid, target_info%vgrid, cf, grib2,               &
       &          ldims=ldims, loutput=loutput, lrestart=lrestart,               &
       &          isteptype=isteptype, tlev_source=tlev_source,                  &
       &          vert_interp=vert_interp, hor_interp=hor_interp,                &
       &          tracer_info=tracer_info, in_group=in_group                     )

    ! Get the number of convection tracers
    IF (lconv_tracer) THEN
      art_config(jg)%nconv_tracer = art_config(jg)%nconv_tracer + 1
    ENDIF

  END SUBROUTINE add_var_list_reference_tracer


END MODULE mo_add_tracer_ref
