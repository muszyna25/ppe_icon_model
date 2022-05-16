!NEC$ options "-O1"
!! Data types and variables used by the ECHAM6 physics package
!!  for ND scalingin cloud microphysics via simple plume model.
!!
!! This module is replicated as per the mo_echam_phy_memory module
!!  with associated credit to authors
!!
!! @author Ross Herbert (Uni of Oxford)
!!
!! @par Revision History
!!
MODULE mo_ndscaling_memory

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL,    &
    &                               GRID_CELL
  USE mo_exception,           ONLY: message, finish
  USE mo_parallel_config,     ONLY: nproma
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_model_domain,        ONLY: t_patch
  USE mo_var_list,            ONLY: add_var, t_var_list_ptr
  USE mo_var_list_register,   ONLY: vlr_add, vlr_del
  USE mo_cf_convention,       ONLY: t_cf_var
  USE mo_grib2,               ONLY: t_grib2_var, grib2_var, OPERATOR(+)
  USE mo_cdi,                 ONLY: DATATYPE_PACK16,  &
    &                               DATATYPE_FLT32,  DATATYPE_FLT64,   &
    &                               GRID_UNSTRUCTURED
  USE mo_zaxis_type,          ONLY: ZA_SURFACE
  USE mo_echam_rad_config,    ONLY: echam_rad_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: prm_ndscaling                           !< variables
  PUBLIC :: prm_ndscaling_list                      !< variable lists
  PUBLIC :: construct_ndscaling_list                !< subroutine
  PUBLIC :: destruct_ndscaling_list                 !< subroutines
  PUBLIC :: t_ndscaling

  !!--------------------------------------------------------------------------
  !!                               DATA TYPES
  !!--------------------------------------------------------------------------

  !>
  !! Derived data type: t_ndscaling
  !!
  !! This structure contains components:
  !! <ol>
  !! <li> related to the instantaneous aerosol radiative forcing
  !! </ol>
  !!
  !! All components are arrays of one of the following shapes:
  !! <ol>
  !! <li> (nproma,           nblks_phy)
  !! <li> (nproma, nlev_phy(+1), nblks_phy)
  !! </ol>
  !! Currently the physics grid has the same spatial resolution as the
  !! dynamics grid, but is unstaggered. This means
  !!
  !!    nlev_phy = nlev
  !!   nblks_phy = patch%nblks_c
  !!
  !! In the long run, the physics grid and dynamics grid may differ in
  !! horizontal and/or vertical resolution, or even in shape.


  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_ndscaling_memory'

  TYPE t_ndscaling

  !--- 2D diagnostic fields:

  REAL(wp),  POINTER :: dNovrN(:,:)=>NULL()  !< cloud droplet number concentration scale factor (MACv2-SP)
  REAL(wp),  POINTER :: migCDNC(:,:)=>NULL() !< CDNC in graupel scheme following scaling by dNovrN

  END TYPE t_ndscaling

  !!--------------------------------------------------------------------------
  !!                          STATE VARIABLES
  !!--------------------------------------------------------------------------
  !! The variable names have the prefix "prm_" in order to emphasize that they
  !! are defined for and used in parameterisations.

  TYPE(t_ndscaling),ALLOCATABLE,TARGET :: prm_ndscaling(:)  !< shape: (n_dom)

  !!--------------------------------------------------------------------------
  !!                          VARIABLE LISTS
  !!--------------------------------------------------------------------------
  TYPE(t_var_list_ptr),ALLOCATABLE :: prm_ndscaling_list(:)  !< shape: (n_dom)

CONTAINS


  !!--------------------------------------------------------------------------
  !!                SUBROUTINES FOR BUILDING AND DELETING VARIABLE LISTS
  !!--------------------------------------------------------------------------
  !>
  !! Top-level procedure for building the physics state
  !!
  SUBROUTINE construct_ndscaling_list( patch_array )

    TYPE(t_patch),INTENT(IN) :: patch_array(:)
    CHARACTER(len=MAX_CHAR_LENGTH) :: listname
    INTEGER :: ndomain, jg, ist, nblks, nlev
    CHARACTER(len=*), PARAMETER :: thissubprog = 'construct_ndscaling_list of mo_ndscaling_memory'

    IF (.NOT.(echam_rad_config(1)%irad_aero >= 33)) RETURN
    CALL message(thismodule,'Construction of ndscaling_list started.')

    ! Allocate pointer arrays prm_field and prm_tend,
    ! as well as the corresponding list arrays.

    ndomain = SIZE(patch_array)

    ALLOCATE( prm_ndscaling(ndomain), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(TRIM(thissubprog), &
      &'allocation of prm_ndscaling array failed')

    ALLOCATE( prm_ndscaling_list(ndomain), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(TRIM(thissubprog), &
      &'allocation of prm_ndscaling list array failed')

    ! Build a ndcaling list for each grid level.
    ! This includes memory allocation.

    DO jg = 1,ndomain

      nblks = patch_array(jg)%nblks_c
      nlev  = patch_array(jg)%nlev

      WRITE(listname,'(a,i2.2)') 'prm_ndscaling_D',jg
      CALL new_ndscaling_list( jg, nproma, nlev, nblks,          &
                             & TRIM(listname), '',             &
                             & prm_ndscaling_list(jg),     &
                             & prm_ndscaling(jg)          )

      CALL message(TRIM(thissubprog),'Construction of ndscaling list finished.')

   END DO

  END SUBROUTINE construct_ndscaling_list
  !--------------------------------------------------------------------
  !>
  !! Release memory used by the ndscaling arrays and list arrays
  !!
  SUBROUTINE destruct_ndscaling_list

    INTEGER :: ndomain  !< total # of grid levels/domains
    INTEGER :: jg       !< grid level/domain index
    INTEGER :: ist      !< system status code
    CHARACTER(len=*), PARAMETER :: thissubprog='destruct_ndscaling_list of mo_ndscaling_memory'
    !---
    IF (.NOT.(echam_rad_config(1)% irad_aero >= 33)) RETURN
    CALL message(TRIM(thissubprog),'Destruction of ndscaling_list started.')

    ndomain = SIZE(prm_ndscaling)

    DO jg = 1,ndomain
      CALL vlr_del( prm_ndscaling_list(jg) )
    ENDDO

    DEALLOCATE( prm_ndscaling_list, STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thissubprog), &
      & 'deallocation of ndscaling list array failed')

    DEALLOCATE( prm_ndscaling, STAT=ist )
    IF (ist/=SUCCESS) CALL finish(TRIM(thissubprog), &
      & 'deallocation of ndscaling array failed')

    CALL message(TRIM(thissubprog),'Destruction of ndscaling_list finished.')

  END SUBROUTINE destruct_ndscaling_list
  !--------------------------------------------------------------------

  SUBROUTINE new_ndscaling_list( k_jg, kproma, klev, kblks,          &
                               & listname, prefix,             &
                               & field_list,     &
                               & field          )

    INTEGER,INTENT(IN) :: k_jg !> patch ID
    INTEGER,INTENT(IN) :: kproma, klev, kblks !< dimension sizes

    CHARACTER(len=*),INTENT(IN) :: listname, prefix

    TYPE(t_var_list_ptr),INTENT(INOUT) :: field_list
    TYPE(t_ndscaling),INTENT(INOUT) :: field

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2)
    INTEGER :: ibits
    INTEGER :: datatype_flt

    ibits = DATATYPE_PACK16

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    shape2d  = (/kproma,       kblks/)

    ! Register a field list and apply default settings

    CALL vlr_add(field_list, TRIM(listname), patch_id=k_jg, lrestart=.TRUE.)

    ! 2D diagnostic fields

    cf_desc    = t_cf_var('dNovrN', '1', 'N_PD / N_PI (MACv2-SP) scaling', datatype_flt)
    grib2_desc = grib2_var(0,6,28, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'dNovrN', field%dNovrN,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d, &
                & lrestart = .FALSE. )

    cf_desc    = t_cf_var('migCDNC', 'cm-3', 'CDNC following scaling by dNovrN in graupel scheme', datatype_flt)
    grib2_desc = grib2_var(0,6,28, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'migCDNC', field%migCDNC,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d, &
                & lrestart = .FALSE. )

  END SUBROUTINE new_ndscaling_list

END MODULE mo_ndscaling_memory

