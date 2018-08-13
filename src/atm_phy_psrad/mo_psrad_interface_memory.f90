!>
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!! @author L. Linardakis, MPI-M, Hamburg
!!
!! @remarks
!!  
MODULE mo_psrad_interface_memory

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH,  & 
    &                               VINTP_METHOD_PRES,         &
    &                               VINTP_METHOD_LIN,          &
    &                               VINTP_METHOD_LIN_NLEVP1
  USE mo_exception,           ONLY: message, finish
  USE mo_parallel_config,     ONLY: nproma
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_model_domain,        ONLY: t_patch
  USE mo_alloc_patches,       ONLY: destruct_patches
  USE mtime,                  ONLY: datetime

  USE mo_linked_list,         ONLY: t_var_list
  USE mo_var_list,            ONLY: default_var_list_settings, &
    &                               add_var,                   &
    &                               new_var_list,              &
    &                               delete_var_list
  USE mo_var_metadata,        ONLY: create_vert_interp_metadata, vintp_types
  USE mo_cf_convention,       ONLY: t_cf_var
  USE mo_grib2,               ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                 ONLY: DATATYPE_PACK16, DATATYPE_PACK24,  &
    &                               DATATYPE_FLT32,  DATATYPE_FLT64,   &
    &                               GRID_UNSTRUCTURED,                 &
    &                               TSTEP_INSTANT, TSTEP_CONSTANT,     &
    &                               TSTEP_MIN, TSTEP_MAX,              &
    &                               cdiInqMissval
  USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL
  USE mo_zaxis_type,          ONLY: ZA_REFERENCE, ZA_REFERENCE_HALF,   &
    &                               ZA_REFERENCE_HALF_HHL, ZA_SURFACE

  USE mo_psrad_interface_namelist, ONLY: number_of_levels
 ! USE mo_radiation_config        , ONLY: isolrad, tsi, tsi_radt, ssi_radt, irad_o3, irad_aero
  USE mo_echam_rad_config       , ONLY: echam_rad_config


  USE mo_psrad_general,      ONLY: nbndsw

  USE mo_grid_config,    ONLY: n_dom

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: construct_psrad_interface_memory                !< subroutine
  PUBLIC :: destruct_psrad_interface_memory                 !< subroutines
  PUBLIC :: t_psrad_interface                               !< derived types

  PUBLIC :: psrad_interface_memory                          !< the memory

  !!--------------------------------------------------------------------------
  !!                               DATA TYPES
  !!--------------------------------------------------------------------------
  TYPE t_psrad_interface_const

    INTEGER  ::             &
         irad_aero,         & !< aerosol control
         no_of_levels         !< number of levels
 
    TYPE(t_patch), POINTER  :: patch

    !! communicated once at the beginning of the run
    REAL(WP),POINTER  ::          &
         zf(:,:,:),               & !< geometric height at full level in m
         zh(:,:,:),               & !< geometric height at half level in m
         dz(:,:,:)                  !< geometric height thickness in m

  END TYPE t_psrad_interface_const
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! calculated within ps_rad
  TYPE t_psrad_interface_parameterized
    
    REAL(WP),POINTER  ::          &
         xm_ch4(:,:,:),           & !< ch4 mass in kg/m2
         xm_n2o(:,:,:),           & !< n2o mass in kg/m2
         xm_cfc(:,:,:,:),         & !< cfc mass in kg/m2
         xm_o2(:,:,:),            &  !< o2  mass in kg/m2
         xm_o3(:,:,:)                !< o3  mass in kg/m2, in cases parameterized

  END TYPE t_psrad_interface_parameterized
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! input: communicated each radiation step
  TYPE t_psrad_interface_in

    REAL(wp) :: psctm                         !< orbit and time dependent solar constant for radiation time step
    REAL(wp) :: ssi_factor(nbndsw)            !< fraction of TSI in the 14 RRTM SW bands

    TYPE(datetime), POINTER ::  this_datetime 

    INTEGER,POINTER  ::  convection_type(:,:)

    LOGICAL,POINTER ::              &
         loland(:,:),               & !< land sea mask, land=.true.    note: saved as float
         loglac(:,:)                  !< glacier mask, glacier=.true.  note: saved as float

    REAL(WP),POINTER  ::          &
         pcos_mu0(:,:),           & !< mu0 for solar zenith angle
         daylght_frc(:,:),        & !< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
         alb_vis_dir(:,:),        & !< surface albedo for vis range and dir light
         alb_nir_dir(:,:),        & !< surface albedo for NIR range and dir light
         alb_vis_dif(:,:),        & !< surface albedo for vis range and dif light
         alb_nir_dif(:,:),        & !< surface albedo for NIR range and dif light
         pp_sfc(:,:),             & !< surface pressure in Pa
         pp_fl(:,:,:),            & !< full level pressure in Pa
         tk_sfc(:,:),             & !< surface temperature in K
         tk_fl(:,:,:),            & !< full level temperature in K
         tk_hl(:,:,:),            & !< half level temperature in K
         xm_dry(:,:,:),           & !< dry air     mass in kg/m2
         xm_vap(:,:,:),           & !< water vapor mass in kg/m2
         xm_liq(:,:,:),           & !< cloud water mass in kg/m2
         xm_ice(:,:,:),           & !< cloud ice   mass in kg/m2
         cdnc(:,:,:),             & !< cloud nuclei concentration
         xc_frc(:,:,:),           & !< fractional cloud cover
         xm_co2(:,:,:)              !< co2 mass in kg/m2

  END TYPE t_psrad_interface_in
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! input: calculated each radiation step from input
!   TYPE t_psrad_interface_calc
!   END TYPE t_psrad_interface_calc
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! flux output: communicated each radiation step
  TYPE t_psrad_interface_out

    REAL(wp), POINTER   ::&
      & lw_dnw(:,:,:),    & !< All-sky   downward longwave  at all levels
      & lw_upw(:,:,:),    & !< All-sky   upward   longwave  at all levels
      & sw_dnw(:,:,:),    & !< All-sky   downward shortwave at all levels
      & sw_upw(:,:,:)       !< All-sky   upward   shortwave at all levels

    REAL (wp), POINTER ::           &
        vis_dn_dir_sfc(:,:)       , & !< Direct downward flux surface visible radiation 
        par_dn_dir_sfc(:,:)       , & !< Direct downward flux surface PAR
        nir_dn_dir_sfc(:,:)       , & !< Direct downward flux surface near-infrared radiation
        vis_dn_dff_sfc(:,:)       , & !< Diffuse  downward flux surface visible radiation 
        par_dn_dff_sfc(:,:)       , & !< Diffuse downward flux surface PAR
        nir_dn_dff_sfc(:,:)       , & !< Diffuse downward flux surface near-infrared radiation
        vis_up_sfc    (:,:)       , & !< Upward  flux surface visible radiation 
        par_up_sfc    (:,:)       , & !< Upward  flux surface PAR
        nir_up_sfc    (:,:)           !< Upward  flux surface near-infrared radiation

  END TYPE t_psrad_interface_out
  !--------------------------------------------------------------------------
 
  !--------------------------------------------------------------------------
  !diagnostics: not communicated
  TYPE t_psrad_interface_diagnostics

    REAL(wp), POINTER  :: &
      & lw_dnw_clr(:,:,:),& !< Clear-sky downward longwave  at all levels
      & lw_upw_clr(:,:,:),& !< Clear-sky upward   longwave  at all levels
      & sw_dnw_clr(:,:,:),& !< Clear-sky downward shortwave at all levels
      & sw_upw_clr(:,:,:)   !< Clear-sky upward   shortwave at all levels
 
  END TYPE t_psrad_interface_diagnostics
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  TYPE t_psrad_interface

    TYPE(t_psrad_interface_const)         :: const
    TYPE(t_psrad_interface_parameterized) :: parameterized
    TYPE(t_psrad_interface_in)            :: in
    TYPE(t_psrad_interface_out)           :: out
!     TYPE(t_psrad_interface_calc)          :: calculated
    TYPE(t_psrad_interface_diagnostics)   :: diagnostics

  END TYPE t_psrad_interface
  !--------------------------------------------------------------------------

  !!--------------------------------------------------------------------------
  !!                          variable lists
  TYPE(t_psrad_interface),ALLOCATABLE,TARGET :: psrad_interface_memory(:)  !< shape: (n_dom)

  TYPE(t_var_list),ALLOCATABLE :: psrad_interface_memory_list(:)  !< shape: (n_dom)

  TYPE(t_patch),POINTER  :: patches(:)

  !--------------------------------------------------------------------------
  INTEGER :: number_of_patches = -1

CONTAINS

  !!--------------------------------------------------------------------------
  !!                SUBROUTINES FOR BUILDING AND DELETING VARIABLE LISTS 
  !!--------------------------------------------------------------------------
  !>
  !! Top-level procedure for building the physics state
  !!
  SUBROUTINE construct_psrad_interface_memory( patch_array )

    TYPE(t_patch),TARGET,INTENT(IN) :: patch_array(:)

    CHARACTER(len=MAX_CHAR_LENGTH) :: listname
    INTEGER :: jg, status
    CHARACTER(len=*), PARAMETER :: method_name='construct_psrad_interface_memory'
    
    !---
    CALL message(TRIM(method_name),'Construction of psrad_interface_memory started.')
    number_of_patches = n_dom
    patches => patch_array
    ALLOCATE( psrad_interface_memory(number_of_patches), STAT=status)
    IF (status/=SUCCESS) CALL finish(TRIM(method_name), &
      &'allocation of psrad_interface_memory array failed')

    ALLOCATE( psrad_interface_memory_list(number_of_patches), STAT=status)
    IF (status/=SUCCESS) CALL finish(TRIM(method_name), &
      &'allocation of psrad_interface_memory_list array failed')

    DO jg = 1,number_of_patches
      WRITE(listname,'(a,i2.2)') 'psrad_interface_',jg
      CALL allocate_psrad_interface_memory(                         &
        & TRIM(listname), '',   psrad_interface_memory_list(jg),    &
        & psrad_interface_memory(jg), patches(jg),  number_of_levels(jg)   )
    END DO

    CALL message(TRIM(method_name),'Construction of psrad interface list finished.')

  END SUBROUTINE construct_psrad_interface_memory
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  !>
  !! Release memory used by the psrad_interface arrays and list arrays
  !!
  SUBROUTINE destruct_psrad_interface_memory

    INTEGER :: jg       !< grid level/domain index
    INTEGER :: status      !< system status code
    CHARACTER(len=*), PARAMETER :: method_name='destruct_psrad_interface_memory'
    !---
    CALL message(TRIM(method_name),'Destruction of psrad_interface_memory started.')

    DO jg = 1,number_of_patches
      CALL delete_var_list( psrad_interface_memory_list(jg) )
    ENDDO

    DEALLOCATE( psrad_interface_memory, STAT=status )
    IF (status/=SUCCESS) CALL finish(TRIM(method_name), &
      & 'deallocation of psrad_interface failed')

    DEALLOCATE( psrad_interface_memory_list, STAT=status )
    IF (status/=SUCCESS) CALL finish(TRIM(method_name), &
      & 'deallocation of psrad_interface_memory_list failed')

    CALL destruct_patches( patches )

    CALL message(TRIM(method_name),'Destruction of psrad_interface_memory finished.')

  END SUBROUTINE destruct_psrad_interface_memory
  !--------------------------------------------------------------------

  !--------------------------------------------------------------------
  SUBROUTINE allocate_psrad_interface_memory(listname, prefix,  field_list, psrad_interface_fields, patch, no_of_levels)

    CHARACTER(len=*),      INTENT(IN)       :: listname, prefix
    TYPE(t_var_list),      INTENT(INOUT)    :: field_list
    TYPE(t_psrad_interface),INTENT(INOUT)   :: psrad_interface_fields
    TYPE(t_patch), TARGET                   :: patch
    INTEGER , INTENT(in)                    :: no_of_levels

    ! Local variables
    INTEGER :: jg !> patch ID
    INTEGER :: alloc_cell_blocks

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d(3), shape3d_layer_interfaces(3)
    INTEGER :: ibits
    INTEGER :: datatype_flt

    ibits = DATATYPE_PACK16
 
    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    alloc_cell_blocks = patch%alloc_cell_blocks
    jg         = patch%id

    psrad_interface_fields%const%patch        => patch
    psrad_interface_fields%const%no_of_levels = no_of_levels
    psrad_interface_fields%const%irad_aero    = echam_rad_config(1)%irad_aero

    shape2d  = (/nproma,       alloc_cell_blocks/)
    shape3d  = (/nproma, no_of_levels, alloc_cell_blocks/)
    shape3d_layer_interfaces = (/nproma,no_of_levels+1,alloc_cell_blocks/)

    ! Register a field list and apply default settings
    CALL new_var_list( field_list, TRIM(listname), patch_id=jg )
    CALL default_var_list_settings( field_list,                &
                                  & lrestart=.FALSE.  )

    !----------------------
    ! const variables
    psrad_interface_fields%const%patch => patch

    grib2_desc = grib2_var(0, 3, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    cf_desc    = t_cf_var('geometric_height_at_half_level', 'm',                &
                &         'Geometric height at half level in physics',          &
                &         datatype_flt)
    CALL add_var( field_list, prefix//'zh', psrad_interface_fields%const%zh,    &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_layer_interfaces,                               &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ),                 &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_CONSTANT  )

    cf_desc    = t_cf_var('geometric_height_at_full_level', 'm',                &
                &         'Geometric height at full level in physics',          &
                &         datatype_flt)
    CALL add_var( field_list, prefix//'zf', psrad_interface_fields%const%zf,    &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=shape3d,                                                &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_CONSTANT  )

    cf_desc    = t_cf_var('geometric_height_thickness', 'm',                    &
                &         'Geometric height thickness in physics',              &
                &         datatype_flt)
    CALL add_var( field_list, prefix//'dz', psrad_interface_fields%const%dz,    &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=shape3d,                                                &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_CONSTANT  )


    !----------------------
    ! parameterized variables
    cf_desc    = t_cf_var('ch4', 'kg/m^2', 'ch4', datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'ch4', psrad_interface_fields%parameterized%xm_ch4, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=shape3d,                                                &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_CONSTANT  )

    cf_desc    = t_cf_var('n2o', 'kg/m^2', 'n2o', datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'n2o', psrad_interface_fields%parameterized%xm_n2o, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=shape3d,                                                &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_CONSTANT  )
    
    cf_desc    = t_cf_var('cfc', 'kg/m^2', 'cfc', datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'cfc', psrad_interface_fields%parameterized%xm_cfc, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=(/nproma, no_of_levels, 2, alloc_cell_blocks/),         &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_CONSTANT  )
    
    cf_desc    = t_cf_var('o2', 'kg/m^2', 'o2', datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'o2', psrad_interface_fields%parameterized%xm_o2, &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=shape3d,                                                &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ),                        &
                & lrestart = .FALSE.,                                           &
                & isteptype=TSTEP_CONSTANT  )
 
    cf_desc    = t_cf_var('o3', 'kg/m^2', 'o3', datatype_flt)
    grib2_desc = grib2_var(0,14,1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'o3', psrad_interface_fields%parameterized%xm_o3,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ) )

    !----------------------
    ! in variables
    cf_desc    = t_cf_var('convection_type', '', 'convection_type (0...3)', datatype_flt)
    grib2_desc = grib2_var(0,6,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'convection_type', psrad_interface_fields%in%convection_type, &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart = .FALSE., ldims=shape2d )

    cf_desc    = t_cf_var('loland', '', '', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'loland', psrad_interface_fields%in%loland,  &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                  &
              & cf_desc, grib2_desc, ldims=shape2d, lrestart=.FALSE. )

    cf_desc    = t_cf_var('loglac', '', '', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'loglac', psrad_interface_fields%in%loglac,  &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                  &
              & cf_desc, grib2_desc, ldims=shape2d, lrestart=.FALSE. )

    cf_desc    = t_cf_var('pcos_mu0', '', 'solar zenith angle', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'pcos_mu0', psrad_interface_fields%in%pcos_mu0,  &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                  &
              & cf_desc, grib2_desc, ldims=shape2d, lrestart=.FALSE. )

    cf_desc    = t_cf_var('daylght_frc', '', 'daylight fraction', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'daylght_frc', psrad_interface_fields%in%daylght_frc,  &
              & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                  &
              & cf_desc, grib2_desc, ldims=shape2d, lrestart=.FALSE. )

    ! &       field% albvisdir (nproma,nblks),          &
    cf_desc    = t_cf_var('albvisdir', '', 'albedo VIS direct', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdir', psrad_interface_fields%in%alb_vis_dir,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.FALSE., ldims=shape2d  )

    ! &       field% albvisdif (nproma,nblks),          &
    cf_desc    = t_cf_var('albvisdif', '', 'albedo VIS diffuse', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albvisdif', psrad_interface_fields%in%alb_vis_dif,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.FALSE., ldims=shape2d  )

    ! &       field% albnirdir (nproma,nblks),          &
    cf_desc    = t_cf_var('albnirdir', '', 'albedo NIR direct', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdir', psrad_interface_fields%in%alb_nir_dir,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.FALSE., ldims=shape2d  )

    ! &       field% albnirdif (nproma,nblks),          &
    cf_desc    = t_cf_var('albnirdif', '', 'albedo NIR diffuse', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'albnirdif', psrad_interface_fields%in%alb_nir_dif,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.FALSE., ldims=shape2d  )


    cf_desc    = t_cf_var('surface_pressure', 'Pa', 'surface pressure', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'pp_sfc', psrad_interface_fields%in%pp_sfc,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.FALSE., ldims=shape2d  )

    cf_desc    = t_cf_var('surface_temperature', 'K', 'surface temperature', datatype_flt)
    grib2_desc = grib2_var(0,0,0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'tk_sfc', psrad_interface_fields%in%tk_sfc,        &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
                & lrestart=.FALSE., ldims=shape2d  )

    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at full levels', datatype_flt)
    grib2_desc = grib2_var(0, 3, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'pp_fl', psrad_interface_fields%in%pp_fl,  &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & ldims=shape3d, lrestart = .FALSE.,                           &
                & vert_interp=create_vert_interp_metadata(                     &
                &             vert_intp_type=vintp_types("Z","I"),             &
                &             vert_intp_method=VINTP_METHOD_PRES ) )
    
    cf_desc    = t_cf_var('temperature_fl', 'K', 'temperature at full levels', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'tk_fl', psrad_interface_fields%in%tk_fl,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d, &
                & lrestart = .FALSE.,                                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ) )

    cf_desc    = t_cf_var('temperature_hl', 'K', 'temperature at half levels', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'tk_hl', psrad_interface_fields%in%tk_hl,  &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & ldims=shape3d_layer_interfaces,                               &
                & lrestart = .FALSE.,                                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN ) )
 
    cf_desc    = t_cf_var('dry_air_mass', 'kg/m2', 'dry air mass in layer', &
         &                datatype_flt)
    grib2_desc = grib2_var(0,1,21, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'xm_dry', psrad_interface_fields%in%xm_dry, &
         &        GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
         &        ldims=shape3d, lrestart = .FALSE., initval=1.0_wp,           &
         &        vert_interp=create_vert_interp_metadata(                     &
         &                    vert_intp_type=vintp_types("P","Z","I"),         & 
         &                    vert_intp_method=VINTP_METHOD_LIN,               &
         &                    l_loglin=.FALSE.,                                &
         &                    l_extrapol=.TRUE., l_pd_limit=.FALSE.,           &
         &                    lower_limit=0._wp  ) )


    cf_desc    = t_cf_var('water_vapor', 'kg/m2', 'water vapor', &
         &                datatype_flt)
    grib2_desc = grib2_var(0,1,64, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'xm_vap', psrad_interface_fields%in%xm_vap, &
         &        GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
         &        ldims=shape3d, lrestart = .FALSE., initval=1.0_wp,           &
         &        vert_interp=create_vert_interp_metadata(                     &
         &                    vert_intp_type=vintp_types("P","Z","I"),         & 
         &                    vert_intp_method=VINTP_METHOD_LIN,               &
         &                    l_loglin=.FALSE.,                                &
         &                    l_extrapol=.TRUE., l_pd_limit=.FALSE.,           &
         &                    lower_limit=0._wp  ) )

    cf_desc    = t_cf_var('cloud_water', 'kg/m2', 'cloud water mass', &
         &                datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'xm_liq', psrad_interface_fields%in%xm_liq, &
         &        GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
         &        ldims=shape3d, lrestart = .FALSE., initval=1.0_wp,           &
         &        vert_interp=create_vert_interp_metadata(                     &
         &                    vert_intp_type=vintp_types("P","Z","I"),         & 
         &                    vert_intp_method=VINTP_METHOD_LIN,               &
         &                    l_loglin=.FALSE.,                                &
         &                    l_extrapol=.TRUE., l_pd_limit=.FALSE.,           &
         &                    lower_limit=0._wp  ) )

    cf_desc    = t_cf_var('cloud_ice', 'kg/m2', 'cloud ice mass', &
         &                datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'xm_ice', psrad_interface_fields%in%xm_ice, &
         &        GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
         &        ldims=shape3d, lrestart = .FALSE., initval=1.0_wp,           &
         &        vert_interp=create_vert_interp_metadata(                     &
         &                    vert_intp_type=vintp_types("P","Z","I"),         & 
         &                    vert_intp_method=VINTP_METHOD_LIN,               &
         &                    l_loglin=.FALSE.,                                &
         &                    l_extrapol=.TRUE., l_pd_limit=.FALSE.,           &
         &                    lower_limit=0._wp  ) )

    cf_desc    = t_cf_var('cloud_nuclei', '', 'cloud nuclei concentration', &
         &                datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'cdnc', psrad_interface_fields%in%cdnc, &
         &        GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
         &        ldims=shape3d, lrestart = .FALSE., initval=1.0_wp,           &
         &        vert_interp=create_vert_interp_metadata(                     &
         &                    vert_intp_type=vintp_types("P","Z","I"),         & 
         &                    vert_intp_method=VINTP_METHOD_LIN,               &
         &                    l_loglin=.FALSE.,                                &
         &                    l_extrapol=.TRUE., l_pd_limit=.FALSE.,           &
         &                    lower_limit=0._wp  ) )

    cf_desc    = t_cf_var('fractional_cloud_cover', '', 'fractional cloud cover', &
         &                datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'xc_frc', psrad_interface_fields%in%xc_frc, &
         &        GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
         &        ldims=shape3d, lrestart = .FALSE., initval=1.0_wp,           &
         &        vert_interp=create_vert_interp_metadata(                     &
         &                    vert_intp_type=vintp_types("P","Z","I"),         & 
         &                    vert_intp_method=VINTP_METHOD_LIN,               &
         &                    l_loglin=.FALSE.,                                &
         &                    l_extrapol=.TRUE., l_pd_limit=.FALSE.,           &
         &                    lower_limit=0._wp  ) )

    cf_desc    = t_cf_var('co2', 'kg/m2', ' co2 mass', &
         &                datatype_flt)
    grib2_desc = grib2_var(55,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( field_list, prefix//'xm_co2', psrad_interface_fields%in%xm_co2, &
         &        GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
         &        ldims=shape3d, lrestart = .FALSE., initval=1.0_wp,           &
         &        vert_interp=create_vert_interp_metadata(                     &
         &                    vert_intp_type=vintp_types("P","Z","I"),         & 
         &                    vert_intp_method=VINTP_METHOD_LIN,               &
         &                    l_loglin=.FALSE.,                                &
         &                    l_extrapol=.TRUE., l_pd_limit=.FALSE.,           &
         &                    lower_limit=0._wp  ) )
 
    !----------------------
    ! out variables
     cf_desc    = t_cf_var('downwelling_longwave_flux_in_air', &
         &                'W m-2'                           , &
         &                'downwelling longwave radiation'  , &
         &                datatype_flt                       )
    grib2_desc = grib2_var(0,5,3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'lw_dnw' , psrad_interface_fields%out%lw_dnw     , &
         &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
         &       cf_desc, grib2_desc                          , &
         &       lrestart = .FALSE.                            , &
         &       ldims=shape3d_layer_interfaces               , &
         &       vert_interp=create_vert_interp_metadata        &
         &         (vert_intp_type=vintp_types("P","Z","I") ,   &
         &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1) )

    cf_desc    = t_cf_var('upwelling_longwave_flux_in_air', &
         &                'W m-2'                         , &
         &                'upwelling longwave radiation'  , &
         &                datatype_flt                     )
    grib2_desc = grib2_var(0,5,4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'lw_upw' , psrad_interface_fields%out%lw_upw     , &
         &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
         &       cf_desc, grib2_desc                          , &
         &       lrestart = .FALSE.                            , &
         &       ldims=shape3d_layer_interfaces               , &
         &       vert_interp=create_vert_interp_metadata        &
         &         (vert_intp_type=vintp_types("P","Z","I") ,   &
         &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1) )

    cf_desc    = t_cf_var('upwelling_shortwave_flux_in_air', &
         &                'W m-2'                          , &
         &                'upwelling shortwave radiation'  , &
         &                datatype_flt                     )
    grib2_desc = grib2_var(0,4,8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'sw_upw' , psrad_interface_fields%out%sw_upw , &
         &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF  , &
         &       cf_desc, grib2_desc                         , &
         &       lrestart = .FALSE.                           , &
         &       ldims=shape3d_layer_interfaces              , &
         &       vert_interp=create_vert_interp_metadata       &
         &         (vert_intp_type=vintp_types("P","Z","I") ,  &
         &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1))

    cf_desc    = t_cf_var('downwelling_shortwave_flux_in_air_assuming_clear_sky', &
         &                'W m-2'                                               , &
         &                'downwelling clear-sky shortwave radiation'           , &
         &                datatype_flt                                          )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'sw_dnw' , psrad_interface_fields%out%sw_dnw , &
         &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
         &       cf_desc, grib2_desc                          , &
         &       lrestart = .FALSE.                            , &
         &       ldims=shape3d_layer_interfaces               , &
         &       vert_interp=create_vert_interp_metadata        &
         &         (vert_intp_type=vintp_types("P","Z","I") ,   &
         &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1) )

    ! at the surface
    cf_desc    = t_cf_var('surface_downwelling_direct_visible_flux_in_air_at_rad_time'    , &
         &                'W m-2'                                                         , &
         &                'surface downwelling direct visible radiation at radiation time', &
         &                datatype_flt                                                    )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'vis_dn_dir_sfc', psrad_interface_fields%out%vis_dn_dir_sfc, &
         &       GRID_UNSTRUCTURED_CELL           , ZA_SURFACE       , &
         &       cf_desc, grib2_desc                                 , &
         &       lrestart = .FALSE.,                                   &
         &       ldims=shape2d                                       )

    cf_desc    = t_cf_var('surface_downwelling_direct_par_flux_in_air_at_rad_time'                          , &
         &                'W m-2'                                                                           , &
         &                'surface downwelling direct photosynthetically active radiation at radiation time', &
         &                datatype_flt                                                                      )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'par_dn_dir_sfc', psrad_interface_fields%out%par_dn_dir_sfc, &
         &       GRID_UNSTRUCTURED_CELL           , ZA_SURFACE       , &
         &       cf_desc, grib2_desc                                 , &
         &       lrestart = .FALSE.,                                   &
         &       ldims=shape2d                                       )

    cf_desc    = t_cf_var('surface_downwelling_direct_nearir_flux_in_air_at_rad_time'           , &
         &                'W m-2'                                                               , &
         &                'surface downwelling direct near infrared radiation at radiation time', &
         &                datatype_flt                                                          )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'nir_dn_dir_sfc', psrad_interface_fields%out%nir_dn_dir_sfc, &
         &       GRID_UNSTRUCTURED_CELL           , ZA_SURFACE       , &
         &       cf_desc, grib2_desc                                 , &
         &       lrestart = .FALSE.,                                   &
         &       ldims=shape2d                                       )


    cf_desc    = t_cf_var('surface_downwelling_diffuse_visible_flux_in_air_at_rad_time'    , &
         &                'W m-2'                                                          , &
         &                'surface downwelling diffuse visible radiation at radiation time', &
         &                datatype_flt                                                     )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'vis_dn_dff_sfc', psrad_interface_fields%out%vis_dn_dff_sfc, &
         &       GRID_UNSTRUCTURED_CELL        , ZA_SURFACE          , &
         &       cf_desc, grib2_desc                                 , &
         &       lrestart = .FALSE.,                                   &
         &       ldims=shape2d                                       )

    cf_desc    = t_cf_var('surface_downwelling_diffuse_par_flux_in_air_at_rad_time'                          , &
         &                'W m-2'                                                                            , &
         &                'surface downwelling diffuse photosynthetically active radiation at radiation time', &
         &                datatype_flt                                                                       )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'par_dn_dff_sfc', psrad_interface_fields%out%par_dn_dff_sfc, &
         &       GRID_UNSTRUCTURED_CELL           , ZA_SURFACE       , &
         &       cf_desc, grib2_desc                                 , &
         &       lrestart = .FALSE.,                                   &
         &       ldims=shape2d                                       )

    cf_desc    = t_cf_var('surface_downwelling_diffuse_nearir_flux_in_air_at_rad_time'           , &
         &                'W m-2'                                                                , &
         &                'surface downwelling diffuse near infrared radiation at radiation time', &
         &                datatype_flt                                                           )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'nir_dn_dff_sfc', psrad_interface_fields%out%nir_dn_dff_sfc, &
         &       GRID_UNSTRUCTURED_CELL           , ZA_SURFACE       , &
         &       cf_desc, grib2_desc                                 , &
         &       lrestart = .FALSE.,                                   &
         &       ldims=shape2d                                       )


    cf_desc    = t_cf_var('surface_upwelling_visible_flux_in_air_at_rad_time'    , &
         &                'W m-2'                                                , &
         &                'surface upwelling visible radiation at radiation time', &
         &                datatype_flt                                           )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'vis_up_sfc', psrad_interface_fields%out%vis_up_sfc, &
         &       GRID_UNSTRUCTURED_CELL       , ZA_SURFACE   , &
         &       cf_desc, grib2_desc                         , &
         &       lrestart = .FALSE.                          , &
         &       ldims=shape2d                               )

    cf_desc    = t_cf_var('surface_upwelling_par_flux_in_air_at_rad_time'                          , &
         &                'W m-2'                                                                  , &
         &                'surface upwelling photosynthetically active radiation at radiation time', &
         &                datatype_flt                                                             )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'par_up_sfc', psrad_interface_fields%out%par_up_sfc, &
         &       GRID_UNSTRUCTURED_CELL       , ZA_SURFACE   , &
         &       cf_desc, grib2_desc                         , &
         &       lrestart = .FALSE.                           , &
         &       ldims=shape2d                               )

    cf_desc    = t_cf_var('surface_upwelling_nearir_flux_in_air_at_rad_time'           , &
         &                'W m-2'                                                      , &
         &                'surface upwelling near infrared radiation at radiation time', &
         &                datatype_flt                                                 )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'nir_up_sfc', psrad_interface_fields%out%nir_up_sfc, &
         &       GRID_UNSTRUCTURED_CELL       , ZA_SURFACE   , &
         &       cf_desc, grib2_desc                         , &
         &       lrestart = .FALSE.                          , &
         &       ldims=shape2d                               )


   !----------------------
   ! diagnostics: not communicated

    cf_desc    = t_cf_var('downwelling_longwave_flux_in_air_assuming_clear_sky', &
         &                'W m-2'                                              , &
         &                'downwelling clear-sky longwave radiation'           , &
         &                datatype_flt                                         )
    grib2_desc = grib2_var(0,5,3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'lw_dnw_clr' , psrad_interface_fields%diagnostics%lw_dnw_clr, &
         &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
         &       cf_desc, grib2_desc                          , &
         &       lrestart = .FALSE.                           , &
         &       ldims=shape3d_layer_interfaces               , &
         &       vert_interp=create_vert_interp_metadata        &
         &         (vert_intp_type=vintp_types("P","Z","I") ,   &
         &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1) )

    cf_desc    = t_cf_var('upwelling_longwave_flux_in_air_assuming_clear_sky', &
         &                'W m-2'                                            , &
         &                'upwelling clear-sky longwave radiation'           , &
         &                datatype_flt                                       )
    grib2_desc = grib2_var(0,5,4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'lw_upw_clr' , psrad_interface_fields%diagnostics%lw_upw_clr, &
         &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
         &       cf_desc, grib2_desc                          , &
         &       lrestart = .FALSE.                           , &
         &       ldims=shape3d_layer_interfaces               , &
         &       vert_interp=create_vert_interp_metadata        &
         &         (vert_intp_type=vintp_types("P","Z","I") ,   &
         &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1) )


    cf_desc    = t_cf_var('downwelling_shortwave_flux_in_air_assuming_clear_sky', &
         &                'W m-2'                                               , &
         &                'downwelling clear-sky shortwave radiation'           , &
         &                datatype_flt                                          )
    grib2_desc = grib2_var(0,4,7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'sw_dnw_clr' , psrad_interface_fields%diagnostics%sw_dnw_clr, &
         &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
         &       cf_desc, grib2_desc                          , &
         &       lrestart = .FALSE.                           , &
         &       ldims=shape3d_layer_interfaces               , &
         &       vert_interp=create_vert_interp_metadata        &
         &         (vert_intp_type=vintp_types("P","Z","I") ,   &
         &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1) )

    cf_desc    = t_cf_var('upwelling_shortwave_flux_in_air_assuming_clear_sky', &
         &                'W m-2'                                             , &
         &                'upwelling clear-sky shortwave radiation'           , &
         &                datatype_flt                                        )
    grib2_desc = grib2_var(0,4,8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(field_list, prefix//'sw_upw_clr' , psrad_interface_fields%diagnostics%sw_upw_clr, &
         &       GRID_UNSTRUCTURED_CELL    , ZA_REFERENCE_HALF   , &
         &       cf_desc, grib2_desc                          , &
         &       lrestart = .FALSE.                           , &
         &       ldims=shape3d_layer_interfaces               , &
         &       vert_interp=create_vert_interp_metadata        &
         &         (vert_intp_type=vintp_types("P","Z","I") ,   &
         &          vert_intp_method=VINTP_METHOD_LIN_NLEVP1) )


  END SUBROUTINE allocate_psrad_interface_memory
  !--------------------------------------------------------------------


END MODULE mo_psrad_interface_memory
