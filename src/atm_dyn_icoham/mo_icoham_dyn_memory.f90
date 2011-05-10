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

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_exception,           ONLY: message,finish
  USE mo_icoham_dyn_types,    ONLY: t_hydro_atm, t_hydro_atm_prog, t_hydro_atm_diag
  USE mo_model_domain,        ONLY: t_patch
  USE mo_run_nml,             ONLY: ltheta_dyn
  USE mo_run_nml,             ONLY: nlev, nlevp1, ntracer, nproma
  USE mo_linked_list,         ONLY: t_var_list
  USE mo_var_list,            ONLY: default_var_list_settings, &
                                  & add_var,                   &
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
  TYPE(t_var_list), PUBLIC, POINTER :: hydro_prog_list(:,:)    !< shape: (n_dom,ntimelevel)
  TYPE(t_var_list), PUBLIC, POINTER :: hydro_diag_list (:)     !< shape: (n_dom)
  TYPE(t_var_list), PUBLIC, POINTER :: hydro_tend_dyn_list(:)  !< shape: (n_dom)
  TYPE(t_var_list), PUBLIC, POINTER :: hydro_tend_phy_list (:) !< shape: (n_dom)
  TYPE(t_var_list), PUBLIC, POINTER :: hydro_prog_out_list(:)  !< shape: (n_dom)
  TYPE(t_var_list), PUBLIC, POINTER :: hydro_diag_out_list (:) !< shape: (n_dom)

CONTAINS

  !!----------------------------------------------------------------------------
  !! Subroutines for allocating/deallocating memory
  !!----------------------------------------------------------------------------
  !>
  !! Subroutine that allocates memory for the state vector on ALL grid levels
  !!
  SUBROUTINE construct_icoham_dyn_state( ntimelevel, p_patch )

    INTEGER,INTENT(IN)     :: ntimelevel
    TYPE(t_patch),INTENT(IN) :: p_patch (:)

    !local variables
    CHARACTER(len=MAX_CHAR_LENGTH) :: listname
    TYPE(t_var_list),POINTER :: listptr

    INTEGER :: ndomain, jg, jt, ist, nblks_c, nblks_e, nblks_v, nlev

    CALL message(TRIM(thismodule),'Construction of 3D dynamics state vector started.')

    ndomain = SIZE(p_patch)

    !allocate list array
    ALLOCATE( hydro_prog_list(ndomain,ntimelevel), hydro_diag_list(ndomain), &
      &       hydro_tend_dyn_list(ndomain), hydro_tend_phy_list(ndomain),    &
      &       hydro_prog_out_list(ndomain), hydro_diag_out_list(ndomain), STAT=ist)
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
        WRITE(listname,'(a,i2.2,a,i2.2)') 'hydro_prog_domain_',jg,'_at_timlev_',jt
        listptr => hydro_prog_list(jg,jt)
        CALL make_hydro_prog_list( nproma, nlev, nblks_c, nblks_e, nblks_v, ntracer, &
             &         TRIM(listname), listptr, p_hydro_state(jg)%prog(jt) )
      END DO

      ! 1.2 Diagnostic variables
      WRITE(listname,'(a,i2.2)') 'hydro_diag_of_domain_',jg
      listptr => hydro_diag_list(jg)
      CALL make_hydro_diag_list( nproma, nlev, nblks_c, nblks_e, nblks_v, ntracer,  &
           &          TRIM(listname), listptr, p_hydro_state(jg)%diag )

      ! 1.3 Tendencies
      WRITE(listname,'(a,i2.2)') 'hydro_tend_dyn_of_domain_',jg
      listptr => hydro_tend_dyn_list(jg)
      CALL make_hydro_tend_dyn_list( nproma, nlev, nblks_c, nblks_e, nblks_v, ntracer, &
           &          TRIM(listname), listptr, p_hydro_state(jg)%tend_dyn )

      WRITE(listname,'(a,i2.2)') 'hydro_tend_phy_of_domain_',jg
      listptr => hydro_tend_phy_list(jg)
      CALL make_hydro_tend_phy_list( nproma, nlev, nblks_c, nblks_e, nblks_v, ntracer, &
           &          TRIM(listname), listptr, p_hydro_state(jg)%tend_phy )

      !----------------------------
      ! 2.  For organizing output
      !----------------------------
      WRITE(listname,'(a,i2.2)') 'hydro_prog_out_of_domain_',jg
      listptr => hydro_prog_out_list(jg)
      CALL make_hydro_prog_list( nproma, nlev, nblks_c, nblks_e, nblks_v, ntracer, &
           &          TRIM(listname), listptr, p_hydro_state(jg)%prog_out )

      WRITE(listname,'(a,i2.2)') 'hydro_diag_out_of_domain_',jg
      listptr => hydro_diag_out_list(jg)
      CALL make_hydro_diag_list( nproma, nlev, nblks_c, nblks_e, nblks_v, ntracer, &
           &           TRIM(listname), listptr, p_hydro_state(jg)%diag_out )

    ENDDO

    NULLIFY(listptr)
    CALL message(TRIM(thismodule),'Construction of 3D dynamics state vector finished.')

  END SUBROUTINE construct_icoham_dyn_state

  !>
  !! Release memory used by the state variable arrays and list arrays
  !!
  SUBROUTINE destruct_icoham_dyn_state

    INTEGER :: ntimelevel, ndomain
    INTEGER :: jt    !< time level index
    INTEGER :: jg    !< grid level/domain index
    INTEGER :: ist   !< system status code

    TYPE(t_var_list),POINTER :: listptr
    !---

    CALL message(TRIM(thismodule),'Destruction of 3D dynamics state vector started.')

    ndomain    = SIZE(p_hydro_state)
    ntimelevel = SIZE(p_hydro_state(1)%prog)

    DO jg = 1,ndomain

      ! Prognostic variables
      DO jt = 1,ntimelevel
        listptr => hydro_prog_list(jg,jt)
        CALL delete_var_list( listptr )
      END DO

      ! Diagnostic variables
      listptr => hydro_diag_list(jg)
      CALL delete_var_list( listptr )

      ! Tendencies
      listptr => hydro_tend_dyn_list(jg)
      CALL delete_var_list( listptr )

      listptr => hydro_tend_phy_list(jg)
      CALL delete_var_list( listptr )

      ! Memory used for organizing output
      listptr => hydro_prog_out_list(jg)
      CALL delete_var_list( listptr )

      listptr => hydro_diag_out_list(jg)
      CALL delete_var_list( listptr )

      DEALLOCATE( p_hydro_state(jg)%prog, STAT=ist )
      IF (ist/=SUCCESS) &
      CALL finish(TRIM(thismodule),'deallocation of prognostic state array failed')

    ENDDO

    NULLIFY(listptr)

    CALL message(TRIM(thismodule),'Destruction of 3D dynamics state vector finished.')

  END SUBROUTINE destruct_icoham_dyn_state

  !>
  !!----------------------------------------------------------------
  !! For variables in "hydro_prog_list": of type "t_hydro_atm_prog"
  !!----------------------------------------------------------------
  !!
  SUBROUTINE make_hydro_prog_list( kproma, klev, kblks_c, kblks_e, kblks_v,  & 
             &                     ktracer, listname, field_list, field    )

    INTEGER,INTENT(IN) :: kproma, klev, kblks_c, kblks_e, kblks_v, ktracer !< dimension sizes

    CHARACTER(len=*),INTENT(IN) :: listname

    TYPE(t_var_list),POINTER :: field_list
    TYPE(t_hydro_atm_prog)   :: field

    ! Local variables
    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d_c(2), shape3d_c(3), shape3d_e(3)
    INTEGER :: shape4d_c(4), ibit

    ibit = 16 !size of var in bits

    CALL new_var_list( field_list, TRIM(listname) )
    CALL default_var_list_settings( field_list )
 
    shape2d_c  = (/kproma,       kblks_c         /)
    shape3d_c  = (/kproma, klev, kblks_c         /)
    shape3d_e  = (/kproma, klev, kblks_e         /)
    shape4d_c  = (/kproma, klev, kblks_c, ktracer/)

    ! Register a field list and apply default settings
    !===========================================================
    ! For now using default value 255 for unknwon parameters
    !===========================================================

    !CF:                   long_name      units     std name
    cf_desc    = t_cf_var('normal_wind', 'm s-1', 'wind normal to the edge')

    !GRIB2: discipline, category, parameter, bits, gridtype, subgridtype, leveltype
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_EDGE)

    CALL add_var( field_list, 'vn', field%vn,          &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HYBRID,&
                & cf_desc, grib2_desc, ldims=shape3d_e )

    cf_desc    = t_cf_var('surface_pressure', 'Pa', 'surface pressure')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'pres_sfc', field%pres_sfc, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,   &
                & cf_desc, grib2_desc, ldims=shape2d_c )

    cf_desc    = t_cf_var('temperature', 'K', 'absolute temperature')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'temp', field%temp,      &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,&
                & cf_desc, grib2_desc, ldims=shape3d_c )

    IF (ltheta_dyn) THEN
      cf_desc    = t_cf_var('potential_temperature', 'K', 'potential temperature')
      grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, 'theta', field%theta,    &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,&
                  & cf_desc, grib2_desc, ldims=shape3d_c )
    ENDIF

    IF (ntracer > 0) THEN
      cf_desc    = t_cf_var('tracer', 'kg kg-1', 'tracer concentration')
      grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, 'tracer', field%tracer,  &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,&
                  & cf_desc, grib2_desc, ldims=shape4d_c )
    ENDIF

    !Initialize all fields with zero
    !Comment it so that it can be initialized with NAN 
!    field% pres_sfc = 0.0_wp
!    field% vn       = 0.0_wp
!    field% temp     = 0.0_wp

!    IF (ltheta_dyn)  field% theta  = 0.0_wp
!    IF (ntracer > 0) field% tracer = 0.0_wp

  END SUBROUTINE make_hydro_prog_list

  !>
  !!----------------------------------------------------------------
  !! For variables in "hydro_diag_list": of type "t_hydro_atm_diag"
  !!----------------------------------------------------------------
  !!
  SUBROUTINE make_hydro_diag_list( kproma, klev, kblks_c, kblks_e, kblks_v,  & 
             &                     ktracer, listname, field_list, field   )

    INTEGER,INTENT(IN) :: kproma, klev, kblks_c, kblks_e, kblks_v, ktracer !< dimension sizes

    CHARACTER(len=*),INTENT(IN) :: listname

    TYPE(t_var_list),POINTER :: field_list
    TYPE(t_hydro_atm_diag)   :: field

    ! Local variables
    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape3d_c(3), shape3d_e(3), shape3d_v(3), shape4d_e(4)
    INTEGER :: shape4d_c(4), ibit, klevp1

    ibit = 16 !size of var in bits

    CALL new_var_list( field_list, TRIM(listname) )
    CALL default_var_list_settings( field_list )

    !Variables at full level

    shape3d_c  = (/kproma, klev, kblks_c         /)
    shape3d_e  = (/kproma, klev, kblks_e         /)
    shape3d_v  = (/kproma, klev, kblks_v         /)
    shape4d_e  = (/kproma, klev, kblks_e, ktracer /)

    cf_desc    = t_cf_var('condensated_water', 'kg kg-1', 'cloud water + cloud ice')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'qx', field%qx,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'zonal wind')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'u', field%u,            &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,&
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('northward_wind', 'm s-1', 'meridional wind')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'v', field%v,             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, &
                & cf_desc, grib2_desc, ldims=shape3d_c  )

    cf_desc    = t_cf_var('tangential_wind', 'm s-1', 'wind tangent to the edge')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( field_list, 'vt', field%vt,           &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HYBRID, &
                & cf_desc, grib2_desc, ldims=shape3d_e  )

    cf_desc    = t_cf_var('vorticity_vertex', 's-1', 'relative vorticity at vertices')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_VERTEX)
    CALL add_var( field_list, 'rel_vort', field%rel_vort, &
                & GRID_UNSTRUCTURED_VERT, ZAXIS_HYBRID,   &
                & cf_desc, grib2_desc, ldims=shape3d_v )

    cf_desc    = t_cf_var('vorticity_edge', 's-1', 'relative vorticity at edges')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( field_list, 'rel_vort_e', field%rel_vort_e, &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HYBRID,       &
                & cf_desc, grib2_desc, ldims=shape3d_e )

    cf_desc    = t_cf_var('vorticity_cell', 's-1', 'relative vorticity at cell center')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'rel_vort_c', field%rel_vort_c,&
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,      &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('divergence', 's-1', 'wind divergence')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'div', field%div,         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, &
                & cf_desc, grib2_desc, ldims=shape3d_c  )

    cf_desc    = t_cf_var('kinetic_energy', 'm2 s-2', 'specific kinetic energy')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'e_kin', field%e_kin,     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('geopotential_full', 'm2 s-2', 'geopotential at full level')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'geo_mc', field%geo_mc,   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('vert_vel_full', 'Pa s-1', 'pressure vertical velocity at full level')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'wpres_mc', field%wpres_mc, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,   &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('pressure_full', 'Pa', 'pressure at full level')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'pres_mc', field%pres_mc, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('delta_pres', 'Pa', 'layer thickness at cell center')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'delp_c', field%delp_c,   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('delta_pres_new', 'Pa',                                        &
                 &        'layer thickness at cell center at next time step')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'delp_c_new', field%delp_c_new, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,       &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('inv_delta_pres', 'Pa-1', 'inverser of layer thickness at cell center')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'rdelp_c', field%rdelp_c, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('inv_delta_pres_new', 'Pa-1',                                    &
                 & 'inverser of layer thickness at cell center at next time step')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'rdelp_c_new', field%rdelp_c_new, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,         &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('delta_pres_e', 'Pa', 'layer thickness at edge center')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( field_list, 'delp_e', field%delp_e,  &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HYBRID,&
                & cf_desc, grib2_desc, ldims=shape3d_e )

    cf_desc    = t_cf_var('delta_pres_v', 'Pa', 'layer thickness at vertex center')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_VERTEX)
    CALL add_var( field_list, 'delp_v', field%delp_v,   &
                & GRID_UNSTRUCTURED_VERT, ZAXIS_HYBRID, &
                & cf_desc, grib2_desc, ldims=shape3d_v )

    cf_desc    = t_cf_var('delta_virtual_temp', 'K', 'virtual temperature increment')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'virt_incr', field%virt_incr, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,     &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('virtual_temperature', 'K', 'virtual temperature')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'tempv', field%tempv,    &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,&
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('rd_log_pr', '', 'Rd * ln(p(k+.5)/p(k-.5)) at cell center')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'rdlnpr_c', field%rdlnpr_c, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,   &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('rd_alpha', '', 'Rd * \alpha at cell center')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'rdalpha_c', field%rdalpha_c, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,     &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('mass_flux_edge', 'Pa m s-1', 'mass flux at edge')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( field_list, 'mass_flux_e', field%mass_flux_e, &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HYBRID,         &
                & cf_desc, grib2_desc, ldims=shape3d_e )

    IF(ltheta_dyn)THEN
      cf_desc    = t_cf_var('exner', ' ', 'exner function')
      grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, 'exner', field%exner,    &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,&
                  & cf_desc, grib2_desc, ldims=shape3d_c )
    END IF

    IF (ntracer > 0) THEN
      cf_desc    = t_cf_var('hor_tracer_flux', 'kg m-1 s-1 ',                          &
                   &        'horizontal tracer flux at edges')
      grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( field_list, 'hfl_tracer', field%hfl_tracer, &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HYBRID,       &
                  & cf_desc, grib2_desc, ldims=shape4d_e )
    END IF

    !================================================================
    !Variables at half level
    !================================================================
    klevp1     = nlevp1
    shape3d_c  = (/kproma, klevp1, kblks_c         /)
    shape4d_c  = (/kproma, klevp1, kblks_c, ktracer /)

    cf_desc    = t_cf_var('geopotential_half', 'm s-1', 'geopotential at half level')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'geo_ic', field%geo_ic,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF, &
                &  cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('vert_vel_half', 'Pa s-1', 'pressure vertical velocity at half level')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'wpres_ic', field%wpres_ic,   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF,&
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('eta_vert_vel', 'Pa s-1', 'eta vertical velocity times dp_deta')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'weta', field%weta,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF,&
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('pressure_half', 'Pa', 'pressure at half level')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'pres_ic', field%pres_ic,     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF,&
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('pressure_half_new', 'Pa', 'pressure at half level at next time step')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'pres_ic_new', field%pres_ic_new, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF,    &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    cf_desc    = t_cf_var('log_pressure', '', 'log of pressure at cell center')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'lnp_ic', field%lnp_ic,        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF, &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    IF (ntracer > 0) THEN
      cf_desc    = t_cf_var('vert_tracer_flux', 'kg m-1 s-1 ', 'horizontal tracer flux at cell')
      grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, 'vfl_tracer', field%vfl_tracer, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF,  &
                  & cf_desc, grib2_desc, ldims=shape4d_c )
    END IF


    ! Initialize all components with zero
    ! Comment it so that it can be initialized with NAN 

!    field% qx          = 0.0_wp
!    field% u           = 0.0_wp
!    field% v           = 0.0_wp
!    field% vt          = 0.0_wp
!    field% rel_vort    = 0.0_wp
!    field% rel_vort_e  = 0.0_wp
!    field% rel_vort_c  = 0.0_wp
!    field% div         = 0.0_wp
!    field% e_kin       = 0.0_wp
!    field% geo_ic      = 0.0_wp
!    field% geo_mc      = 0.0_wp
!    field% wpres_mc    = 0.0_wp
!    field% wpres_ic    = 0.0_wp
!    field% weta        = 0.0_wp
!    field% pres_ic     = 0.0_wp
!    field% pres_ic_new = 0.0_wp
!    field% pres_mc     = 0.0_wp
!    field% delp_c      = 0.0_wp
!    field% delp_c_new  = 0.0_wp
!    field% rdelp_c     = 0.0_wp
!    field% rdelp_c_new = 0.0_wp
!    field% delp_e      = 0.0_wp
!    field% delp_v      = 0.0_wp
!    field% virt_incr   = 0.0_wp
!    field% tempv       = 0.0_wp
!    field% rdlnpr_c    = 0.0_wp
!    field% rdalpha_c   = 0.0_wp
!    field% lnp_ic      = 0.0_wp
!    field% mass_flux_e = 0.0_wp
!
!    IF (ltheta_dyn) field% exner =0.0_wp
!
!    IF (ntracer > 0) THEN
!      field% hfl_tracer = 0.0_wp
!      field% vfl_tracer = 0.0_wp
!    END IF

  END SUBROUTINE make_hydro_diag_list

  !>
  !!----------------------------------------------------------------
  !! For variables in "hydro_tend_dyn_list": of type "t_hydro_atm_prog"
  !!----------------------------------------------------------------
  !!
  SUBROUTINE make_hydro_tend_dyn_list( kproma, klev, kblks_c, kblks_e, kblks_v,   &
             &                         ktracer, listname, field_list, field   )

    INTEGER,INTENT(IN) :: kproma, klev, kblks_c, kblks_e, kblks_v, ktracer !< dimension sizes

    CHARACTER(len=*),INTENT(IN) :: listname

    TYPE(t_var_list),POINTER :: field_list
    TYPE(t_hydro_atm_prog)   :: field

    ! Local variables
    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d_c(2), shape3d_c(3), shape3d_e(3)
    INTEGER :: shape4d_c(4), ibit

    ibit = 16 !size of var in bits

    CALL new_var_list( field_list, TRIM(listname) )
    CALL default_var_list_settings( field_list )
 
    shape2d_c  = (/kproma,       kblks_c         /)
    shape3d_c  = (/kproma, klev, kblks_c         /)
    shape3d_e  = (/kproma, klev, kblks_e         /)
    shape4d_c  = (/kproma, klev, kblks_c, ktracer/)

    ! Register a field list and apply default settings

    cf_desc    = t_cf_var('normal_wind_tendency', 'm s-2', 'tendency of wind normal to edge')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( field_list, 'tend_vn', field%vn,      &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HYBRID, &
                & cf_desc, grib2_desc, ldims=shape3d_e )

    cf_desc    = t_cf_var('surface_pressure_tendency', 'Pa s-1', 'surface pressure tendency')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'tend_pres_sfc', field%pres_sfc, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,       &
                & cf_desc, grib2_desc, ldims=shape2d_c )

    cf_desc    = t_cf_var('temperature_tendency', 'K s-1', 'absolute temperature tendency')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'tend_temp', field%temp, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,&
                & cf_desc, grib2_desc, ldims=shape3d_c )

    IF (ltheta_dyn) THEN
      cf_desc    = t_cf_var('potential_temperature_tendency', 'K s-1',  &
                   &        'potential temperature tendency')
      grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, 'tend_theta', field%theta, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,  &
                  & cf_desc, grib2_desc, ldims=shape3d_c )
    ENDIF

    IF (ntracer > 0) THEN
      cf_desc    = t_cf_var('tracer_tendency', 's-1', 'tendency of tracer concentration')
      grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, 'tend_tracer', field%tracer, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,    &
                  & cf_desc, grib2_desc, ldims=shape4d_c )
    ENDIF

    !Initialize all fields with zero
    !Comment it so that it can be initialized with NAN 
!    field% pres_sfc = 0.0_wp
!    field% vn       = 0.0_wp
!    field% temp     = 0.0_wp

!    IF (ltheta_dyn)  field% theta  = 0.0_wp
!    IF (ntracer > 0) field% tracer = 0.0_wp

  END SUBROUTINE make_hydro_tend_dyn_list

  !>
  !!----------------------------------------------------------------
  !! For variables in "hydro_tend_phy_list": of type "t_hydro_atm_prog"
  !!----------------------------------------------------------------
  !!
  SUBROUTINE make_hydro_tend_phy_list( kproma, klev, kblks_c, kblks_e, kblks_v,   &
             &                         ktracer, listname, field_list, field   )

    INTEGER,INTENT(IN) :: kproma, klev, kblks_c, kblks_e, kblks_v, ktracer !< dimension sizes

    CHARACTER(len=*),INTENT(IN) :: listname

    TYPE(t_var_list),POINTER :: field_list
    TYPE(t_hydro_atm_prog)   :: field

    ! Local variables
    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d_c(2), shape3d_c(3), shape3d_e(3)
    INTEGER :: shape4d_c(4), ibit

    ibit = 16 !size of var in bits

    CALL new_var_list( field_list, TRIM(listname) )
    CALL default_var_list_settings( field_list )
 
    shape2d_c  = (/kproma,       kblks_c         /)
    shape3d_c  = (/kproma, klev, kblks_c         /)
    shape3d_e  = (/kproma, klev, kblks_e         /)
    shape4d_c  = (/kproma, klev, kblks_c, ktracer/)

    ! Register a field list and apply default settings

    cf_desc    = t_cf_var('phy_normal_wind_tendency', 'm s-2',      &
                 &        'physical tendency of wind normal to edge')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( field_list, 'phy_tend_vn', field%vn,  &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HYBRID, &
                & cf_desc, grib2_desc, ldims=shape3d_e )

    cf_desc    = t_cf_var('phy_surface_pressure_tendency', 'Pa s-1',&
                 &        'physical surface pressure tendency')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'phy_tend_pres_sfc', field%pres_sfc,&
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,          &
                & cf_desc, grib2_desc, ldims=shape2d_c )

    cf_desc    = t_cf_var('phy_temperature_tendency', 'K s-1', &
                 &        'physical absolute temperature tendency')
    grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
    CALL add_var( field_list, 'phy_tend_temp', field%temp, &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,    &
                & cf_desc, grib2_desc, ldims=shape3d_c )

    IF (ltheta_dyn) THEN
      cf_desc    = t_cf_var('phy_potential_temperature_tendency', 'K s-1',&
                   &        'physical potential temperature tendency')
      grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, 'phy_tend_theta', field%theta, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,      &
                  & cf_desc, grib2_desc, ldims=shape3d_c )
    ENDIF

    IF (ntracer > 0) THEN
      cf_desc    = t_cf_var('phy_tracer_tendency', 's-1', &
                   &       'physical tendency of tracer concentration')
      grib2_desc = t_grib2_var(255, 255, 255, ibit, GRID_REFERENCE, GRID_CELL)
      CALL add_var( field_list, 'phy_tend_tracer', field%tracer, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID,        &
                  & cf_desc, grib2_desc, ldims=shape4d_c )
    ENDIF

    !Initialize all fields with zero
    !Comment it so that it can be initialized with NAN 
!    field% pres_sfc = 0.0_wp
!    field% vn       = 0.0_wp
!    field% temp     = 0.0_wp

!    IF (ltheta_dyn)  field% theta  = 0.0_wp
!    IF (ntracer > 0) field% tracer = 0.0_wp

  END SUBROUTINE make_hydro_tend_phy_list
!------------------------------------------------------------------

END module mo_icoham_dyn_memory
