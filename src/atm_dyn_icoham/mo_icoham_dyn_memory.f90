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
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!
MODULE mo_icoham_dyn_memory

  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH, VNAME_LEN, MAX_NTRACER
  USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_EDGE, GRID_UNSTRUCTURED_CELL,   &
    &                               GRID_EDGE, GRID_CELL, GRID_VERTEX,                &
    &                               GRID_UNSTRUCTURED_VERT
  USE mo_exception,           ONLY: message,finish
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm, t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_model_domain,        ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma
  USE mo_advection_config,    ONLY: advection_config
  USE mo_ha_dyn_config,       ONLY: ha_dyn_config
  USE mo_linked_list,         ONLY: t_var_list
  USE mo_var_list,            ONLY: default_var_list_settings, &
                                  & add_var, add_ref,          &
                                  & new_var_list,              &
                                  & delete_var_list
  USE mo_cf_convention,       ONLY: t_cf_var
  USE mo_grib2,               ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                 ONLY: DATATYPE_PACK16, DATATYPE_FLT32, DATATYPE_FLT64, GRID_UNSTRUCTURED
  USE mo_zaxis_type,          ONLY: ZA_SURFACE, ZA_REFERENCE, ZA_REFERENCE_HALF
  USE mo_io_config,           ONLY: lnetcdf_flt64_output

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: p_hydro_state
  PUBLIC :: construct_icoham_dyn_state, destruct_icoham_dyn_state

  !----------------------------------------------------------------------------
  !                          MEMORY BUFFER 
  !----------------------------------------------------------------------------

  TYPE(t_hydro_atm),TARGET,ALLOCATABLE :: p_hydro_state(:) !< state vector on
                                                           !< different grid levels
                                                           !< shape: (n_dom)

  !--------------------------------------------------------------------------
  !                          VARIABLE LISTS
  !--------------------------------------------------------------------------
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
  SUBROUTINE construct_icoham_dyn_state( ntimelevel, ntracer, p_patch )

    INTEGER,      INTENT(IN) :: ntimelevel, ntracer
    TYPE(t_patch),INTENT(IN) :: p_patch(:)

    INTEGER :: ndomain, jg, jt, istat, nblks_c, nblks_e, nblks_v, nlev
    CHARACTER(len=MAX_CHAR_LENGTH) :: listname, varname_prefix
    CHARACTER(len=VNAME_LEN) :: tracer_names(MAX_NTRACER) 

    CHARACTER(len=*),PARAMETER ::  &
             routine = 'mo_icoham_dyn_memory:construct_icoham_dyn_state'

    !---
    CALL message('','')
    CALL message(TRIM(routine),'Construction of 3D dynamics state vector started.')

    ndomain = SIZE(p_patch)

    ! Allocate state array

    ALLOCATE (p_hydro_state(ndomain), stat=istat)
    IF (istat /= success) THEN
      CALL finish(TRIM(routine),'allocation of p_hydro_state failed')
    ENDIF

    ! Allocate list arrays

    ALLOCATE( hydro_prog_list(ndomain,ntimelevel), &
            & hydro_diag_list(ndomain),            &
            & hydro_tend_dyn_list(ndomain),        &
            & hydro_tend_phy_list(ndomain),        &
            & hydro_prog_out_list(ndomain),        &
            & hydro_diag_out_list(ndomain),        &
            & STAT=istat)

    IF (istat/=SUCCESS) CALL finish(TRIM(routine), &
      &'allocation of hydrostatic prog/diag list array failed')

    ! Build a field list and a tendency list for each grid level.
    ! This includes memory allocation. 

    DO jg = 1,ndomain

      nblks_c = p_patch(jg)%nblks_c
      nblks_e = p_patch(jg)%nblks_e
      nblks_v = p_patch(jg)%nblks_v
      nlev    = p_patch(jg)%nlev


      ! get tracer names
      tracer_names = advection_config(jg)%tracer_names

      !----------------------------
      ! 1.  For time integration:
      !----------------------------
      ! 1.1 Prognostic variables

      ALLOCATE(p_hydro_state(jg)%prog(1:ntimelevel), STAT=istat)
      IF (istat/=SUCCESS) &
      CALL finish(TRIM(routine),'allocation of prognostic state array failed')

      DO jt = 1,ntimelevel

        WRITE(listname,'(a,i2.2,a,i2.2)')  'hydro_prog_D',jg,'_timlev',jt
        WRITE(varname_prefix,'(a,i2.2,a)') 'ha_prog_TL',jt,'_'

        CALL new_hydro_prog_list( jg, nproma, nlev, ntracer, tracer_names, &
                                & ha_dyn_config%ltheta_dyn,                &
                                & nblks_c, nblks_e, TRIM(listname),        &
                                & TRIM(varname_prefix),                    &
                                & hydro_prog_list(jg,jt),                  &
                                & p_hydro_state(jg)%prog(jt),              &
                                & store_in_restart=.TRUE.                          )
      END DO

      ! 1.2 Diagnostic variables

      WRITE(listname,'(a,i2.2)')  'hydro_diag_D',jg
      WRITE(varname_prefix,'(a)') 'ha_diag_'
      CALL new_hydro_diag_list( jg, nproma, nlev, ntracer, tracer_names, &
                              & ha_dyn_config%ltheta_dyn,             &
                              & nblks_c, nblks_e, nblks_v,            &
                              & TRIM(listname), TRIM(varname_prefix), &
                              & hydro_diag_list(jg),                  &
                              & p_hydro_state(jg)%diag,               &
                              & store_in_restart=.FALSE.                      )

      ! 1.3 Tendencies

      WRITE(listname,'(a,i2.2)')  'hydro_tend_dyn_D',jg
      WRITE(varname_prefix,'(a)') 'ha_tend_dyn_'
      CALL new_hydro_prog_list( jg, nproma, nlev, ntracer, tracer_names, &
                              & ha_dyn_config%ltheta_dyn,                &
                              & nblks_c, nblks_e, TRIM(listname),        &
                              & TRIM(varname_prefix),                    &
                              & hydro_tend_dyn_list(jg),                 &
                              & p_hydro_state(jg)%tend_dyn,              &
                              & store_in_restart=.TRUE.                          )

      WRITE(listname,'(a,i2.2)')  'hydro_tend_phy_D',jg
      WRITE(varname_prefix,'(a)') 'ha_tend_phy_'
      CALL new_hydro_prog_list( jg, nproma, nlev, ntracer, tracer_names, &
                              & ha_dyn_config%ltheta_dyn,                &
                              & nblks_c, nblks_e, TRIM(listname),        &
                              & TRIM(varname_prefix),                    &
                              & hydro_tend_phy_list(jg),                 &
                              & p_hydro_state(jg)%tend_phy,              &
                              & store_in_restart=.FALSE.                         )

      !----------------------------
      ! 2.  For organizing output
      !----------------------------
      WRITE(listname,'(a,i2.2)')  'hydro_prog_out_D',jg
      WRITE(varname_prefix,'(a)') 'ha_prog_out_'
      CALL new_hydro_prog_list( jg, nproma, nlev, ntracer, tracer_names, &
                              & ha_dyn_config%ltheta_dyn,                &
                              & nblks_c, nblks_e, TRIM(listname),        &
                              & TRIM(varname_prefix),                    &
                              & hydro_prog_out_list(jg),                 &
                              & p_hydro_state(jg)%prog_out,              &
                              & store_in_restart=.FALSE.                         )

      WRITE(listname,'(a,i2.2)')  'hydro_diag_out_D',jg
      WRITE(varname_prefix,'(a)') 'ha_diag_out_'
      CALL new_hydro_diag_list( jg, nproma, nlev, ntracer, tracer_names, &
                              & ha_dyn_config%ltheta_dyn,             &
                              & nblks_c, nblks_e, nblks_v,            &
                              & TRIM(listname), TRIM(varname_prefix), &
                              & hydro_diag_out_list(jg),              &
                              & p_hydro_state(jg)%diag_out,           &
                              & store_in_restart=.FALSE.                      )

    ENDDO

    CALL message(TRIM(routine),'Construction of 3D dynamics state vector finished.')
    CALL message('','')

  END SUBROUTINE construct_icoham_dyn_state

  !>
  !! Release memory used by the state variable arrays and list arrays
  !!
  SUBROUTINE destruct_icoham_dyn_state

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_icoham_dyn_memory:destruct_icoham_dyn_state'

    INTEGER :: ntimelevel, ndomain
    INTEGER :: jt    !< time level index
    INTEGER :: jg    !< grid level/domain index
    INTEGER :: istat !< system status code

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

      DEALLOCATE( p_hydro_state(jg)%prog, STAT=istat )
      IF (istat/=SUCCESS) &
      CALL finish(TRIM(routine),'deallocation of prognostic state array failed')

    ENDDO

    DEALLOCATE (p_hydro_state, STAT=istat)
    IF (istat /= SUCCESS) & 
    CALL finish(TRIM(routine),'deallocation for p_hydro_state failed')

    CALL message(TRIM(routine),'Destruction of 3D dynamics state vector finished.')

  END SUBROUTINE destruct_icoham_dyn_state

  !>
  !!
  !!
  SUBROUTINE new_hydro_prog_list( k_jg, kproma, klev, ktracer,    &
                                & tracer_names, ltheta_dyn,       &
                                & kblks_c, kblks_e,               &
                                & listname, vname_prefix,         &
                                & field_list, field, store_in_restart )

    INTEGER,INTENT(IN) :: k_jg                       !< patch ID
    INTEGER,INTENT(IN) :: kproma, klev, ktracer  !< dimension sizes
    LOGICAL,INTENT(IN) :: ltheta_dyn
    INTEGER,INTENT(IN) :: kblks_c, kblks_e       !< dimension sizes
    LOGICAL,INTENT(IN) :: store_in_restart               !< store in restart file?

    CHARACTER(len=*),INTENT(IN) :: listname, vname_prefix
    CHARACTER(len=VNAME_LEN) :: tracer_names(:)      !< tracer-specific name suffixes

    TYPE(t_var_list)      ,INTENT(INOUT) :: field_list
    TYPE(t_hydro_atm_prog),INTENT(INOUT) :: field

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d_c(2), shape3d_c(3), shape3d_e(3)
    INTEGER :: ibits, jtrc
    INTEGER :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ibits = DATATYPE_PACK16 !size of var in bits

    ! Register a variable list and apply default settings

    CALL new_var_list( field_list, TRIM(listname), patch_id=k_jg )
    CALL default_var_list_settings( field_list,                &
                                  & lrestart=store_in_restart  )

    ! Add variables to the list 

    shape2d_c  = (/kproma,       kblks_c/)
    shape3d_c  = (/kproma, klev, kblks_c/)
    shape3d_e  = (/kproma, klev, kblks_e/)

    cf_desc    = t_cf_var('normal_wind', 'm s-1', 'wind normal to the edge', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( field_list, vname_prefix//'vn', field%vn,  &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,         &
                & cf_desc, grib2_desc, ldims=shape3d_e )

    cf_desc    = t_cf_var('surface_air_pressure', 'Pa', 'surface pressure', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'pres_sfc', field%pres_sfc, &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                   &
                & cf_desc, grib2_desc, ldims=shape2d_c )

    cf_desc    = t_cf_var('air_temperature', 'K', 'absolute temperature', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'temp', field%temp, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,            &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    IF (ltheta_dyn) THEN
      cf_desc    = t_cf_var('potential_temperature', 'K', 'potential temperature', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, vname_prefix//'theta', field%theta, &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,              &
                  & cf_desc, grib2_desc, ldims=shape3d_c )
    ENDIF

    IF (ktracer > 0) THEN

      ! Tracer array for (model) internal use

      CALL add_var( field_list, vname_prefix//'tracer', field%tracer,           &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                  & t_cf_var('tracer', 'kg kg-1', 'tracer concentration', datatype_flt), &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
                  & ldims = (/kproma,klev,kblks_c,ktracer/),                    &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.        )

      ! Reference to individual tracer, for I/O

      ALLOCATE(field%tracer_ptr(ktracer))
      DO jtrc = 1,ktracer

        CALL add_ref( field_list, vname_prefix//'tracer',                          &
                    & vname_prefix//'q'//TRIM(tracer_names(jtrc)),                 &
                    & field%tracer_ptr(jtrc)%p,                                    &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                    & t_cf_var('tracer_'//TRIM(tracer_names(jtrc)), 'kg kg-1', '', datatype_flt),&
                    & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                    & ldims=(/kproma,klev,kblks_c/))

      END DO
    ENDIF ! ktracer > 0

  END SUBROUTINE new_hydro_prog_list

  !>
  !!
  !!
  SUBROUTINE new_hydro_diag_list( k_jg, kproma, klev, ktracer, tracer_names, &
                                & ltheta_dyn,                          &
                                & kblks_c, kblks_e, kblks_v,           & 
                                & listname, vname_prefix,              &
                                & field_list, field, store_in_restart )

    INTEGER,INTENT(IN) :: k_jg                       !< patch ID
    INTEGER,INTENT(IN) :: kproma, klev, ktracer      !< dimension sizes
    LOGICAL,INTENT(IN) :: ltheta_dyn
    INTEGER,INTENT(IN) :: kblks_c, kblks_e, kblks_v  !< dimension sizes
    LOGICAL,INTENT(IN) :: store_in_restart

    CHARACTER(len=*),INTENT(IN) :: listname, vname_prefix
    CHARACTER(len=VNAME_LEN) :: tracer_names(:)      !< tracer-specific name suffixes

    TYPE(t_var_list)      ,INTENT(INOUT) :: field_list
    TYPE(t_hydro_atm_diag),INTENT(INOUT) :: field

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape3d_c(3), shape3d_e(3), shape3d_v(3)
    INTEGER :: ibits, klevp1, jtrc
    INTEGER :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ibits = DATATYPE_PACK16 !size of var in bits

    ! Register a variable list and apply default settings

    CALL new_var_list( field_list, TRIM(listname), patch_id=k_jg )
    CALL default_var_list_settings( field_list, lrestart=store_in_restart ) 

    !----------------------------
    ! Add variables to the list 
    !----------------------------
    ! Variables on full levels

    shape3d_c  = (/kproma, klev, kblks_c         /)
    shape3d_e  = (/kproma, klev, kblks_e         /)
    shape3d_v  = (/kproma, klev, kblks_v         /)

    cf_desc    = t_cf_var('condensated_water', 'kg kg-1', 'cloud water + cloud ice', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'qx', field%qx, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,        &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'zonal wind', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'u', field%u, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,      &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('northward_wind', 'm s-1', 'meridional wind', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'v', field%v, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,      &
                & cf_desc, grib2_desc, ldims=shape3d_c  )

    cf_desc    = t_cf_var('tangential_wind', 'm s-1', 'wind tangent to the edge', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( field_list, vname_prefix//'vt', field%vt, &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,        &
                & cf_desc, grib2_desc, ldims=shape3d_e  )

    cf_desc    = t_cf_var('vorticity_vertex', 's-1', 'relative vorticity at vertices', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_VERTEX)
    CALL add_var( field_list, vname_prefix//'rel_vort', field%rel_vort, &
                & GRID_UNSTRUCTURED_VERT, ZA_REFERENCE,                    &
                & cf_desc, grib2_desc, ldims=shape3d_v )

    cf_desc    = t_cf_var('vorticity_edge', 's-1', 'relative vorticity at edges', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( field_list, vname_prefix//'rel_vort_e', field%rel_vort_e, &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                        &
                & cf_desc, grib2_desc, ldims=shape3d_e )

    cf_desc    = t_cf_var('vorticity_cell', 's-1', 'relative vorticity at cell center', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'rel_vort_c', field%rel_vort_c,&
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                       &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('divergence', 's-1', 'wind divergence', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'div', field%div,         &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                  &
                & cf_desc, grib2_desc, ldims=shape3d_c  )

    cf_desc    = t_cf_var('kinetic_energy', 'm2 s-2', 'specific kinetic energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'e_kin', field%e_kin,     &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                  &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('geopotential', 'm2 s-2', 'geopotential at full level', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'geo_mc', field%geo_mc,   &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                  &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('vert_vel_full', 'Pa s-1', 'pressure vertical velocity at full level', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'wpres_mc', field%wpres_mc, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                    &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('pressure_full', 'Pa', 'pressure at full level', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'pres_mc', field%pres_mc, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                  &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('delta_pres', 'Pa', 'layer thickness at cell center', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'delp_c', field%delp_c,   &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                  &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('delta_pres_new', 'Pa',                           &
                 &        'layer thickness at cell center at next time step', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'delp_c_new', field%delp_c_new, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('inv_delta_pres', 'Pa-1', 'inverser of layer thickness at cell center', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'rdelp_c', field%rdelp_c, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                  &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('inv_delta_pres_new', 'Pa-1',                          &
                 & 'inverser of layer thickness at cell center at next time step', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'rdelp_c_new', field%rdelp_c_new, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('delta_pres_e', 'Pa', 'layer thickness at edge center', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( field_list, vname_prefix//'delp_e', field%delp_e, &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                &
                & cf_desc, grib2_desc, ldims=shape3d_e )

    cf_desc    = t_cf_var('delta_pres_v', 'Pa', 'layer thickness at vertex center', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_VERTEX)
    CALL add_var( field_list, vname_prefix//'delp_v', field%delp_v, &
                & GRID_UNSTRUCTURED_VERT, ZA_REFERENCE,                &
                & cf_desc, grib2_desc, ldims=shape3d_v )

    cf_desc    = t_cf_var('delta_virtual_temp', 'K', 'virtual temperature increment', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'virt_incr', field%virt_incr, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                      &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('virtual_temperature', 'K', 'virtual temperature', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'tempv', field%tempv,&
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,             &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('rd_log_pr', '', 'Rd * ln(p(k+.5)/p(k-.5)) at cell center', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'rdlnpr_c', field%rdlnpr_c, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                    &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('rd_alpha', '', 'Rd * \alpha at cell center', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'rdalpha_c', field%rdalpha_c, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                      &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('mass_flux_edge', 'Pa m s-1', 'mass flux at edge', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( field_list, vname_prefix//'mass_flux_e', field%mass_flux_e, &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                          &
                & cf_desc, grib2_desc, ldims=shape3d_e )

    IF(ltheta_dyn)THEN
      cf_desc    = t_cf_var('exner', ' ', 'exner function', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( field_list, vname_prefix//'exner', field%exner,&
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,             &
                  & cf_desc, grib2_desc, ldims=shape3d_c )
    END IF

    !--------------------------
    ! Variables on half levels
    !--------------------------

    klevp1     = klev + 1
    shape3d_c  = (/kproma, klevp1, kblks_c/)

    cf_desc    = t_cf_var('geopotential', 'm s-1', 'geopotential at half level', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'geo_ic', field%geo_ic,&
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF,          &
                &  cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('vert_vel_half', 'Pa s-1', 'pressure vertical velocity at half level', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'wpres_ic', field%wpres_ic,&
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF,              &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('eta_vert_vel', 'Pa s-1', 'eta vertical velocity times dp_deta', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'weta', field%weta,&
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF,      &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('pressure_half', 'Pa', 'pressure at half level', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'pres_ic', field%pres_ic,&
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF,            &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('pressure_half_new', 'Pa', 'pressure at half level at next time step', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'pres_ic_new', field%pres_ic_new, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF,                     &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('log_pressure', '', 'log of pressure at cell center', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, vname_prefix//'lnp_ic', field%lnp_ic,&
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF,          &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    !------------------------
    ! Tracers
    !------------------------
    IF (ktracer > 0) THEN

      ! Tracer flux arrays for (model) internal use

      CALL add_var( field_list, vname_prefix//'hfl_tracer', field%hfl_tracer,          &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                                 &
                  & t_cf_var('hfl_tracer', 'kg m-1 s-1', 'horizontal flux of tracer', &
                  &          datatype_flt), &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE),     &
                  & ldims = (/kproma,klev,kblks_e,ktracer/),                           &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.               )

      CALL add_var( field_list, vname_prefix//'vfl_tracer', field%vfl_tracer,        &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF,                          &
                  & t_cf_var('vfl_tracer', 'kg m-1 s-1', 'vertical flux of tracer',  &
                  &          datatype_flt), &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
                  & ldims = (/kproma,klevp1,kblks_c,ktracer/),                       &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.             )

      ! Referrence to fluxes individual tracer, for I/O

      ALLOCATE(field%hfl_tracer_ptr(ktracer))
      ALLOCATE(field%vfl_tracer_ptr(ktracer))

      DO jtrc = 1,ktracer

        CALL add_ref( field_list, vname_prefix//'hfl_tracer',                       &
                    & vname_prefix//'hfl_q'//TRIM(tracer_names(jtrc)),              &
                    & field%hfl_tracer_ptr(jtrc)%p,                                 &
                    & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                            &
                    & t_cf_var('hfl_q'//TRIM(tracer_names(jtrc)), 'kg m-1 s-1', '', &
                    &          datatype_flt), &
                    & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE),&
                    & ldims = (/kproma,klev,kblks_e/)                               )

        CALL add_ref( field_list, vname_prefix//'vfl_tracer',                       &
                    & vname_prefix//'vfl_q'//TRIM(tracer_names(jtrc)),              &
                    & field%vfl_tracer_ptr(jtrc)%p,                                 &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF,                       &
                    & t_cf_var('vfl_q'//TRIM(tracer_names(jtrc)), 'kg m-1 s-1', '', &
                    &          datatype_flt), &
                    & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
                    & ldims = (/kproma,klevp1,kblks_c/)                             )

      END DO
    ENDIF ! ktracer > 0

  END SUBROUTINE new_hydro_diag_list

END module mo_icoham_dyn_memory
