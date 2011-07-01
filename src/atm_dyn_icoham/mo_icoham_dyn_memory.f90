!>
!! Type definition and variable declaration for the dynamical core of ICOHAM.
!!
!! This module contains the variables ("state vectors") used for
!! storing the prognostic and diagnostic variables of the hydrostatic
!! atmospheric dynamical core, as well as the tendency of the prognostic
!! variables.
!! Constructors and destructors for these data structures are also
!! defined here.
!!
!! This module is the substitute and extended version of the old
!! "mo_hydro_state".
!!
!! This module has been changed for new memory allocation infrastructure (29-04-2012)
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
!!     violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!     copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!      an according license agreement with DWD and MPI-M.
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
MODULE mo_icoham_dyn_memory

  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_exception,           ONLY: message,finish
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm, t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_model_domain,        ONLY: t_patch
  USE mo_run_nml,             ONLY: ltheta_dyn
  USE mo_run_nml,             ONLY: ntracer, nproma
  USE mo_advection_nml,       ONLY: ctracer_list
  USE mo_linked_list,         ONLY: t_var_list
  USE mo_var_list,            ONLY: default_var_list_settings, &
                                  & add_var, add_ref,          &
                                  & new_var_list,              &
                                  & delete_var_list
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants 


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: p_hydro_state
  PUBLIC :: construct_icoham_dyn_state, destruct_icoham_dyn_state

  CHARACTER(len=*), PARAMETER :: version = '$Id$'
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_icoham_dyn_memory'

  !!----------------------------------------------------------------------------
  !! Memory buffer
  !!----------------------------------------------------------------------------

  TYPE(t_hydro_atm),TARGET,ALLOCATABLE :: p_hydro_state(:) !< state vector on
                                                           !< different grid levels
                                                           !< shape: (n_dom)

  !!--------------------------------------------------------------------------
  !!                          VARIABLE LISTS
  !!--------------------------------------------------------------------------
  TYPE(t_var_list),PUBLIC,ALLOCATABLE :: hydro_prog_list(:,:)    !< shape: (n_dom,ntimelevel)
  TYPE(t_var_list),PUBLIC,ALLOCATABLE :: hydro_diag_list(:)      !< shape: (n_dom)
  TYPE(t_var_list),PUBLIC,ALLOCATABLE :: hydro_tend_dyn_list(:)  !< shape: (n_dom)
  TYPE(t_var_list),PUBLIC,ALLOCATABLE :: hydro_tend_phy_list(:)  !< shape: (n_dom)
  TYPE(t_var_list),PUBLIC,ALLOCATABLE :: hydro_prog_out_list(:)  !< shape: (n_dom)
  TYPE(t_var_list),PUBLIC,ALLOCATABLE :: hydro_diag_out_list(:)  !< shape: (n_dom)

CONTAINS

  !!----------------------------------------------------------------------------
  !! Subroutines for allocating/deallocating memory
  !!----------------------------------------------------------------------------
  !>
  !! Subroutine that allocates memory for the state vector on ALL grid levels
  !!
  SUBROUTINE construct_icoham_dyn_state( ntimelevel, p_patch )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_icoham_dyn_memory:construct_icoham_dyn_state'

    INTEGER,INTENT(IN)       :: ntimelevel
    TYPE(t_patch),INTENT(IN) :: p_patch (:)

    !local variables
    CHARACTER(len=MAX_CHAR_LENGTH) :: listname, varname_prefix

    INTEGER :: ndomain, jg, jt, ist, nblks_c, nblks_e, nblks_v, nlev

    CALL message(TRIM(routine),'Construction of 3D dynamics state vector started.')

    ndomain = SIZE(p_patch)

    ! Allocate list arrays

    ALLOCATE( hydro_prog_list(ndomain,ntimelevel), &
            & hydro_diag_list(ndomain),            &
            & hydro_tend_dyn_list(ndomain),        &
            & hydro_tend_phy_list(ndomain),        &
            & hydro_prog_out_list(ndomain),        &
            & hydro_diag_out_list(ndomain),        &
            & STAT=ist)

    IF (ist/=SUCCESS) CALL finish(TRIM(thismodule), &
      &'allocation of hydrostatic prog/diag list array failed')

    ! Build a field list and a tendency list for each grid level.
    ! This includes memory allocation. 

    DO jg = 1,ndomain

      nblks_c = p_patch(jg)%nblks_c
      nblks_e = p_patch(jg)%nblks_e
      nblks_v = p_patch(jg)%nblks_v
      nlev    = p_patch(jg)%nlev

      !----------------------------
      ! 1.  For time integration:
      !----------------------------
      ! 1.1 Prognostic variables

      ALLOCATE(p_hydro_state(jg)%prog(1:ntimelevel), STAT=ist)
      IF (ist/=SUCCESS) &
      CALL finish(TRIM(thismodule),'allocation of prognostic state array failed')

      DO jt = 1,ntimelevel

        WRITE(listname,'(a,i2.2,a,i2.2)')  'hydro_prog_D',jg,'_timlev',jt
        WRITE(varname_prefix,'(a,i2.2,a)') 'ha_prog_TL',jt,'_'

        CALL new_hydro_prog_list( nproma, nlev, ntracer, nblks_c, nblks_e, &
                                & TRIM(listname), TRIM(varname_prefix),    &
                                & hydro_prog_list(jg,jt),                  &
                                & p_hydro_state(jg)%prog(jt),              &
                                & lrestart=.TRUE.                          )
      END DO

      ! 1.2 Diagnostic variables

      WRITE(listname,'(a,i2.2)')  'hydro_diag_D',jg
      WRITE(varname_prefix,'(a)') 'ha_diag_'
      CALL new_hydro_diag_list( nproma, nlev, ntracer,                &
                              & nblks_c, nblks_e, nblks_v,            &
                              & TRIM(listname), TRIM(varname_prefix), &
                              & hydro_diag_list(jg),                  &
                              & p_hydro_state(jg)%diag,               &
                              & lrestart=.FALSE.                      )

      ! 1.3 Tendencies

      WRITE(listname,'(a,i2.2)')  'hydro_tend_dyn_D',jg
      WRITE(varname_prefix,'(a)') 'ha_tend_dyn_'
      CALL new_hydro_prog_list( nproma, nlev, ntracer, nblks_c, nblks_e, &
                              & TRIM(listname), TRIM(varname_prefix),    &
                              & hydro_tend_dyn_list(jg),                 &
                              & p_hydro_state(jg)%tend_dyn,              &
                              & lrestart=.TRUE.                          )

      WRITE(listname,'(a,i2.2)')  'hydro_tend_phy_D',jg
      WRITE(varname_prefix,'(a)') 'ha_tend_phy_'
      CALL new_hydro_prog_list( nproma, nlev, ntracer, nblks_c, nblks_e, &
                              & TRIM(listname), TRIM(varname_prefix),    &
                              & hydro_tend_phy_list(jg),                 &
                              & p_hydro_state(jg)%tend_phy,              &
                              & lrestart=.FALSE.                         )

      !----------------------------
      ! 2.  For organizing output
      !----------------------------
      WRITE(listname,'(a,i2.2)')  'hydro_prog_out_D',jg
      WRITE(varname_prefix,'(a)') 'ha_prog_out_'
      CALL new_hydro_prog_list( nproma, nlev, ntracer, nblks_c, nblks_e, &
                              & TRIM(listname), TRIM(varname_prefix),    &
                              & hydro_prog_out_list(jg),                 &
                              & p_hydro_state(jg)%prog_out,              &
                              & lrestart=.FALSE.                         )

      WRITE(listname,'(a,i2.2)')  'hydro_diag_out_D',jg
      WRITE(varname_prefix,'(a)') 'ha_diag_out_'
      CALL new_hydro_diag_list( nproma, nlev, ntracer,                &
                              & nblks_c, nblks_e, nblks_v,            &
                              & TRIM(listname), TRIM(varname_prefix), &
                              & hydro_diag_out_list(jg),              &
                              & p_hydro_state(jg)%diag_out,           &
                              & lrestart=.FALSE.                      )

    ENDDO

    CALL message(TRIM(routine),'Construction of 3D dynamics state vector finished.')

  END SUBROUTINE construct_icoham_dyn_state

  !>
  !! Release memory used by the state variable arrays and list arrays
  !!
  SUBROUTINE destruct_icoham_dyn_state

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_icoham_dyn_memory:destruct_icoham_dyn_state'

    INTEGER :: ntimelevel, ndomain
    INTEGER :: jt    !< time level index
    INTEGER :: jg    !< grid level/domain index
    INTEGER :: ist   !< system status code

    !---

    CALL message(TRIM(routine),'Destruction of 3D dynamics state vector started.')

    ndomain    = SIZE(p_hydro_state)
    ntimelevel = SIZE(p_hydro_state(1)%prog)

    DO jg = 1,ndomain

      ! Prognostic variables
      DO jt = 1,ntimelevel
        CALL delete_var_list( hydro_prog_list(jg,jt) )
      END DO

      ! Diagnostic variables
      CALL delete_var_list( hydro_diag_list(jg) )

      ! Tendencies
      CALL delete_var_list( hydro_tend_dyn_list(jg) )
      CALL delete_var_list( hydro_tend_phy_list(jg) )

      ! Memory used for organizing output
      CALL delete_var_list( hydro_prog_out_list(jg) )
      CALL delete_var_list( hydro_diag_out_list(jg) )

      DEALLOCATE( p_hydro_state(jg)%prog, STAT=ist )
      IF (ist/=SUCCESS) &
      CALL finish(TRIM(thismodule),'deallocation of prognostic state array failed')

    ENDDO

    CALL message(TRIM(routine),'Destruction of 3D dynamics state vector finished.')

  END SUBROUTINE destruct_icoham_dyn_state

  !>
  !!
  !!
  SUBROUTINE new_hydro_prog_list( kproma, klev, ktracer,      &
                                & kblks_c, kblks_e,           &
                                & listname, vname_prefix,     &
                                & field_list, field, lrestart )

    INTEGER,INTENT(IN) :: kproma, klev, ktracer  !< dimension sizes
    INTEGER,INTENT(IN) :: kblks_c, kblks_e       !< dimension sizes
    LOGICAL,INTENT(IN) :: lrestart               !< store in restart file?

    CHARACTER(len=*),INTENT(IN) :: listname, vname_prefix

    TYPE(t_var_list)      ,INTENT(INOUT) :: field_list
    TYPE(t_hydro_atm_prog),INTENT(INOUT) :: field

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d_c(2), shape3d_c(3), shape3d_e(3)
    INTEGER :: ibit, jtrc

    ibit = 16 !size of var in bits

    ! Register a variable list and apply default settings

    CALL new_var_list( field_list, TRIM(listname) )
    CALL default_var_list_settings( field_list,                &
                                  & lrestart=lrestart,         &
                                  & restart_type=FILETYPE_NC2  )

    ! Add variables to the list 

    shape2d_c  = (/kproma,       kblks_c/)
    shape3d_c  = (/kproma, klev, kblks_c/)
    shape3d_e  = (/kproma, klev, kblks_e/)

    cf_desc    = t_cf_var('normal_wind', 'm s-1', 'wind normal to the edge')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( field_list, vname_prefix//'vn', field%vn,  &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HYBRID,      &
                & cf_desc, grib2_desc, ldims=shape3d_e )

    cf_desc    = t_cf_var('surface_pressure', 'Pa', 'surface pressure')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'pres_sfc', field%pres_sfc, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                 &
                & cf_desc, grib2_desc, ldims=shape2d_c )

    cf_desc    = t_cf_var('temperature', 'K', 'absolute temperature')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'temp', field%temp, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,         &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    IF (ltheta_dyn) THEN
      cf_desc    = t_cf_var('potential_temperature', 'K', 'potential temperature')
      grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, vname_prefix//'theta', field%theta, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,           &
                  & cf_desc, grib2_desc, ldims=shape3d_c )
    ENDIF

    IF (ktracer > 0) THEN

      ! Tracer array for (model) internal use

      CALL add_var( field_list, vname_prefix//'tracer', field%tracer,           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                       &
                  & t_cf_var('tracer', 'kg kg-1', 'tracer concentration'),      &
                  & t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL),&
                  & ldims = (/kproma,klev,kblks_c,ktracer/),                    &
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.          )

      ! Reference to individual tracer, for I/O

      ALLOCATE(field%tracer_ptr(ktracer))
      DO jtrc = 1,ktracer

        CALL add_ref( field_list, vname_prefix//'tracer',                          &
                    & vname_prefix//'q'//ctracer_list(jtrc:jtrc),                  &
                    & field%tracer_ptr(jtrc)%p,                                    &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                        &
                    & t_cf_var('tracer_'//ctracer_list(jtrc:jtrc), 'kg kg-1', ''), &
                    & t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL), &
                    & ldims=(/kproma,klev,kblks_c/))

      END DO
    ENDIF ! ktracer > 0

  END SUBROUTINE new_hydro_prog_list

  !>
  !!
  !!
  SUBROUTINE new_hydro_diag_list( kproma, klev, ktracer,      &
                                & kblks_c, kblks_e, kblks_v,  & 
                                & listname, vname_prefix,     &
                                & field_list, field, lrestart )

    INTEGER,INTENT(IN) :: kproma, klev, ktracer      !< dimension sizes
    INTEGER,INTENT(IN) :: kblks_c, kblks_e, kblks_v  !< dimension sizes
    LOGICAL,INTENT(IN) :: lrestart

    CHARACTER(len=*),INTENT(IN) :: listname, vname_prefix

    TYPE(t_var_list)      ,INTENT(INOUT) :: field_list
    TYPE(t_hydro_atm_diag),INTENT(INOUT) :: field

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape3d_c(3), shape3d_e(3), shape3d_v(3), shape4d_e(4)
    INTEGER :: ibit, klevp1, jtrc

    ibit = 16 !size of var in bits

    ! Register a variable list and apply default settings

    CALL new_var_list( field_list, TRIM(listname) )
    CALL default_var_list_settings( field_list, lrestart=lrestart ) 

    !----------------------------
    ! Add variables to the list 
    !----------------------------
    ! Variables on full levels

    shape3d_c  = (/kproma, klev, kblks_c         /)
    shape3d_e  = (/kproma, klev, kblks_e         /)
    shape3d_v  = (/kproma, klev, kblks_v         /)
    shape4d_e  = (/kproma, klev, kblks_e, ktracer /)

    cf_desc    = t_cf_var('condensated_water', 'kg kg-1', 'cloud water + cloud ice')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'qx', field%qx, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,     &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'zonal wind')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'u', field%u, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,   &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('northward_wind', 'm s-1', 'meridional wind')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'v', field%v, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,   &
                & cf_desc, grib2_desc, ldims=shape3d_c  )

    cf_desc    = t_cf_var('tangential_wind', 'm s-1', 'wind tangent to the edge')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( field_list, vname_prefix//'vt', field%vt, &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HYBRID,     &
                & cf_desc, grib2_desc, ldims=shape3d_e  )

    cf_desc    = t_cf_var('vorticity_vertex', 's-1', 'relative vorticity at vertices')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_VERTEX)
    CALL add_var( field_list, vname_prefix//'rel_vort', field%rel_vort, &
                & GRID_UNSTRUCTURED_VERT, ZAXIS_HYBRID,                 &
                & cf_desc, grib2_desc, ldims=shape3d_v )

    cf_desc    = t_cf_var('vorticity_edge', 's-1', 'relative vorticity at edges')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( field_list, vname_prefix//'rel_vort_e', field%rel_vort_e, &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HYBRID,                     &
                & cf_desc, grib2_desc, ldims=shape3d_e )

    cf_desc    = t_cf_var('vorticity_cell', 's-1', 'relative vorticity at cell center')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'rel_vort_c', field%rel_vort_c,&
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                    &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('divergence', 's-1', 'wind divergence')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'div', field%div,         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,               &
                & cf_desc, grib2_desc, ldims=shape3d_c  )

    cf_desc    = t_cf_var('kinetic_energy', 'm2 s-2', 'specific kinetic energy')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'e_kin', field%e_kin,     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,               &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('geopotential_full', 'm2 s-2', 'geopotential at full level')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'geo_mc', field%geo_mc,   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,               &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('vert_vel_full', 'Pa s-1', 'pressure vertical velocity at full level')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'wpres_mc', field%wpres_mc, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                 &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('pressure_full', 'Pa', 'pressure at full level')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'pres_mc', field%pres_mc, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,               &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('delta_pres', 'Pa', 'layer thickness at cell center')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'delp_c', field%delp_c,   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,               &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('delta_pres_new', 'Pa',                           &
                 &        'layer thickness at cell center at next time step')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'delp_c_new', field%delp_c_new, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                     &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('inv_delta_pres', 'Pa-1', 'inverser of layer thickness at cell center')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'rdelp_c', field%rdelp_c, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,               &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('inv_delta_pres_new', 'Pa-1',                          &
                 & 'inverser of layer thickness at cell center at next time step')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'rdelp_c_new', field%rdelp_c_new, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                       &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('delta_pres_e', 'Pa', 'layer thickness at edge center')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( field_list, vname_prefix//'delp_e', field%delp_e, &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HYBRID,             &
                & cf_desc, grib2_desc, ldims=shape3d_e )

    cf_desc    = t_cf_var('delta_pres_v', 'Pa', 'layer thickness at vertex center')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_VERTEX)
    CALL add_var( field_list, vname_prefix//'delp_v', field%delp_v, &
                & GRID_UNSTRUCTURED_VERT, ZAXIS_HYBRID,             &
                & cf_desc, grib2_desc, ldims=shape3d_v )

    cf_desc    = t_cf_var('delta_virtual_temp', 'K', 'virtual temperature increment')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'virt_incr', field%virt_incr, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                   &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('virtual_temperature', 'K', 'virtual temperature')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'tempv', field%tempv,&
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,          &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('rd_log_pr', '', 'Rd * ln(p(k+.5)/p(k-.5)) at cell center')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'rdlnpr_c', field%rdlnpr_c, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                 &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('rd_alpha', '', 'Rd * \alpha at cell center')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'rdalpha_c', field%rdalpha_c, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,                   &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('mass_flux_edge', 'Pa m s-1', 'mass flux at edge')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( field_list, vname_prefix//'mass_flux_e', field%mass_flux_e, &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HYBRID,                       &
                & cf_desc, grib2_desc, ldims=shape3d_e )

    IF(ltheta_dyn)THEN
      cf_desc    = t_cf_var('exner', ' ', 'exner function')
      grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, vname_prefix//'exner', field%exner,&
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,          &
                  & cf_desc, grib2_desc, ldims=shape3d_c )
    END IF

    !--------------------------
    ! Variables on half levels
    !--------------------------

    klevp1     = klev + 1
    shape3d_c  = (/kproma, klevp1, kblks_c/)

    cf_desc    = t_cf_var('geopotential_half', 'm s-1', 'geopotential at half level')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'geo_ic', field%geo_ic,&
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF,       &
                &  cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('vert_vel_half', 'Pa s-1', 'pressure vertical velocity at half level')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'wpres_ic', field%wpres_ic,&
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF,           &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('eta_vert_vel', 'Pa s-1', 'eta vertical velocity times dp_deta')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'weta', field%weta,&
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF,   &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('pressure_half', 'Pa', 'pressure at half level')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'pres_ic', field%pres_ic,&
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF,         &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('pressure_half_new', 'Pa', 'pressure at half level at next time step')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'pres_ic_new', field%pres_ic_new, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF,                  &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('log_pressure', '', 'log of pressure at cell center')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'lnp_ic', field%lnp_ic,&
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF,       &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    !------------------------
    ! Tracers
    !------------------------
    IF (ktracer > 0) THEN

      ! Tracer flux arrays for (model) internal use

      CALL add_var( field_list, vname_prefix//'hfl_tracer', field%hfl_tracer,          &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HYBRID,                              &
                  & t_cf_var('hfl_tracer', 'kg m-1 s-1', 'horizontal flux of tracer'), &
                  & t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_EDGE),       &
                  & ldims = (/kproma,klev,kblks_e,ktracer/),                           &
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.                 )

      CALL add_var( field_list, vname_prefix//'vfl_tracer', field%vfl_tracer,        &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF,                       &
                  & t_cf_var('vfl_tracer', 'kg m-1 s-1', 'vertical flux of tracer'), &
                  & t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL),     &
                  & ldims = (/kproma,klevp1,kblks_c,ktracer/),                       &
                  & lcontainer=.TRUE., lrestart=.FALSE., lpost=.FALSE.               )

      ! Referrence to fluxes individual tracer, for I/O

      ALLOCATE(field%hfl_tracer_ptr(ktracer))
      ALLOCATE(field%vfl_tracer_ptr(ktracer))

      DO jtrc = 1,ktracer

        CALL add_ref( field_list, vname_prefix//'hfl_tracer',                       &
                    & vname_prefix//'hfl_q'//ctracer_list(jtrc:jtrc),               &
                    & field%hfl_tracer_ptr(jtrc)%p,                                 &
                    & GRID_UNSTRUCTURED_EDGE, ZAXIS_HYBRID,                         &
                    & t_cf_var('hfl_q'//ctracer_list(jtrc:jtrc), 'kg m-1 s-1', ''), &
                    & t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_EDGE),  &
                    & ldims = (/kproma,klev,kblks_e/)                               )

        CALL add_ref( field_list, vname_prefix//'vfl_tracer',                       &
                    & vname_prefix//'vfl_q'//ctracer_list(jtrc:jtrc),               &
                    & field%vfl_tracer_ptr(jtrc)%p,                                 &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF,                    &
                    & t_cf_var('vfl_q'//ctracer_list(jtrc:jtrc), 'kg m-1 s-1', ''), &
                    & t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL),  &
                    & ldims = (/kproma,klevp1,kblks_c/)                             )

      END DO
    ENDIF ! ktracer > 0

  END SUBROUTINE new_hydro_diag_list

END module mo_icoham_dyn_memory
