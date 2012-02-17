!>
!! Definition, allocation/deallocation and reading of external datasets
!!
!! This module contains the type-declaration for the external datasets,
!! including memory allocation/deallocation and reading.
!!
!! @author Daniel Reinert, DWD
!! @author Hermann Asensio, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2010-07-12)
!! Modification by Hermann Asensio, DWD (2010-07-16)
!!  - add miscellaneous variables for external parameters
!! Modification by Daniel Reinert, DWD (2011-05-03)
!! - Memory allocation method changed from explicit allocation to Luis'
!!   infrastructure
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
MODULE mo_ext_data

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max
  USE mo_parallel_config,    ONLY: nproma
  USE mo_impl_constants,     ONLY: inwp, iecham, ildf_echam, io3_clim, io3_ape, &
    &                              ihs_ocean, ihs_atm_temp, ihs_atm_theta, inh_atmosphere, &
    &                              max_char_length, min_rlcell_int
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_physical_constants, ONLY: ppmv2gg, zemiss_def
  USE mo_run_config,         ONLY: iforcing
  USE mo_ocean_nml,          ONLY: iforc_oce, iforc_omip, iforc_len
  USE mo_extpar_config,      ONLY: itopo, fac_smooth_topo, n_iter_smooth_topo, l_emiss, &
    &                              heightdiff_threshold,                                &
    &                              extpar_filename, generate_filename
  USE mo_dynamics_config,    ONLY: iequations
  USE mo_radiation_config,   ONLY: irad_o3,irad_aero
  USE mo_model_domain,       ONLY: t_patch
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_grid_config,        ONLY: n_dom, nroot, dynamics_grid_filename
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_intp,               ONLY: cells2verts_scalar
  USE mo_math_laplace,       ONLY: nabla2_scalar, nabla4_scalar
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_sync,               ONLY: SYNC_C, SYNC_V, sync_patch_array
  USE mo_mpi,                ONLY: my_process_is_stdio, p_io, p_bcast, &
    &                              p_comm_work_test, p_comm_work
  USE mo_parallel_config,    ONLY: p_test_run
  USE mo_linked_list,        ONLY: t_var_list
  USE mo_var_list,           ONLY: default_var_list_settings, &
    &                              add_var,                   &
    &                              new_var_list,              &
    &                              delete_var_list
  USE mo_master_nml,         ONLY: model_base_dir
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants
  USE mo_util_netcdf,        ONLY: read_netcdf_data, read_netcdf_lu, nf
  USE mo_util_string,        ONLY: t_keyword_list,  &
    &                              associate_keyword, with_keywords


  IMPLICIT NONE

  ! required for testing/reading topography
  INCLUDE 'netcdf.inc'

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  INTEGER::  nlev_o3, nmonths

  CHARACTER(len=6)  :: levelname
  CHARACTER(len=4)  :: zlevelname
  CHARACTER(len=6)  :: cellname
  CHARACTER(len=5)  :: o3name
  CHARACTER(len=20) :: o3unit



  PUBLIC :: t_external_data
  PUBLIC :: t_external_atmos
  PUBLIC :: t_external_atmos_td
  PUBLIC :: t_external_ocean
!  PUBLIC :: external_ocean_td
  PUBLIC :: ext_data
  PUBLIC :: nclass_lu

  PUBLIC :: init_ext_data
  PUBLIC :: destruct_ext_data
  PUBLIC :: smooth_topography 

  PUBLIC :: nlev_o3,nmonths
	


  !>
  !! atmosphere external data class
  !!
  !! atmosphere external data class
  !!
  TYPE :: t_external_atmos

    !
    ! *** Topography ***
    REAL(wp), POINTER ::   &   !< topographic height at cell centers      [m]
      &  topography_c(:,:)     ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< smoothed topographic height at cell centers [m]
      &  topography_smt_c(:,:) ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< topographic height at cell edges        [m]
      &  topography_e(:,:)     ! index1=1,nproma, index2=1,nblks_e

    REAL(wp), POINTER ::   &   !< topographic height at cell vertices     [m]
      &  topography_v(:,:)     ! index1=1,nproma, index2=1,nblks_v

    REAL(wp), POINTER ::   &   !< smoothed topographic height at vertices [m]
      &  topography_smt_v(:,:) ! index1=1,nproma, index2=1,nblks_v

    REAL(wp), POINTER ::   &   !< geopotential (S)                        [m**2/s**2]
      &  fis(:,:)              ! index1=1,nproma, index2=1,nblks_c


    !
    ! *** Land-Sea-Mask ***

    LOGICAL, POINTER  ::   &   !< land-sea-mask for cell centers          [ ]
      &  llsm_atm_c(:,:)       ! .TRUE. if landpoint
                               ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< fraction land in a grid element         [ ]
      &  fr_land(:,:)          ! 0. for water, 1.0 indicates 100% land
                               ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::    &  !< fraction land glacier in a grid element [ ]
      &  fr_glac(:,:)          ! 1.0 indicates 100% glacier
                               ! index1=1,nproma, index2=1,nblks_c    

    REAL(wp), POINTER ::   &   !< fraction sea ice cover in a grid element [ ]
      &  fr_ice(:,:)           ! 1.0 indicates 100% ice
                               ! index1=1,nproma, index2=1,nblks_c 
   
    REAL(wp), POINTER ::   &   !< fraction land in a grid element         [ ]
      &  fr_land_smt(:,:)      !  = smoothed fr_land

    REAL(wp), POINTER ::   &   !< fraction sea ice cover in a grid element [ ]
      &  fr_ice_smt(:,:)       ! = smoothed fr_ice

    REAL(wp), POINTER ::   &   !< fraction land glacier in a grid element [ ]
      &  fr_glac_smt(:,:)      ! = smoothed fr_glac


    !  
    ! *** roughness length ***
    REAL(wp), POINTER ::   &   !< surface roughness                       [m]
      &  z0(:,:)               ! index1=1,nproma, index2=1,nblks_c


    !
    ! *** FLake ***
    REAL(wp), POINTER ::   &   !< fraction of fresh water                 [ ]
      &  fr_lake(:,:)          ! as partition of total area of the
                               ! grid element
                               ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< lake depth                              [m]
      &  depth_lk(:,:)         ! index1=1,nproma, index2=1,nblks_c


    !
    ! *** subgrid scale orography ***
    REAL(wp), POINTER ::   &   !< standard deviation of sub-grid scale orography [m]
      &  sso_stdh(:,:)         ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< anisotropy of sub-grid scale orography  [ ]
      &  sso_gamma(:,:)        ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< angle betw. principal axis of orography and E [rad]
      &  sso_theta(:,:)        ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< mean slope of sub-grid scale orography  [ ]
      &  sso_sigma(:,:)        ! index1=1,nproma, index2=1,nblks_c


    !
    ! *** vegetation parameters ***
    REAL(wp), POINTER ::   &   !< ground fraction covered by plants (vegetation period)  [ ]
      & plcov_mx(:,:)          ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< leaf area index (vegetation period)     [ ]
      &  lai_mx(:,:)           ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< root depth                              [m]
      &  rootdp(:,:)           ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< ground fraction covered by evergreen forest [ ]
      &  for_e(:,:)            ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< ground fraction covered by deciduous forest [ ]
      &  for_d(:,:)            ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< urban area fraction                     [ ]
      &  urban(:,:)            ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< minimum value of stomata resistance     [ s/m ]
      &  rsmin(:,:)            ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< annual maximum NDVI                     [ ]
      &  ndvi_max(:,:)         ! index1=1,nproma, index2=1,nblks_c



    !
    ! *** soil parameters ***
    INTEGER, POINTER  ::   &   !< soil texture, keys 0-9                  [ ]
      &  soiltyp(:,:)          ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< Near surface temperature (climatological mean)  [ K ]
      &  t_cl(:,:)             !  used as climatological layer (deepest layer) of T_SO
                               ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &   !< longwave surface emissivity             [ ]
      &  emis_rad(:,:)         ! index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::  &    !< Landuse class fraction                  [ ]
      & lu_class_fraction(:,:,:) ! index1=1,nproma, index2=1,nblks_c, index3=1,nclass_lu

  END TYPE t_external_atmos


  INTEGER, ALLOCATABLE :: nclass_lu(:)  !< number of landuse classes
                                        !< dim: n_dom


  !>
  !! atmosphere external data class (time dependent)
  !!
  !! Contains auxiliary time dependent versions of some external atmospheric data
  !! fields which are already defined in external_atmos. These fields will
  !! be used to store e.g. montly means from which updated external data can be
  !! derived. The updated interpolated fields can be copied into the time independent
  !! counterparts, which are defined in external_atmos.
  !!
  TYPE :: t_external_atmos_td

    ! *** radiation parameters ***
    REAL(wp), POINTER ::   &   !< aerosol optical thickness of black carbon    [ ]
      &  o3(:,:,:,:)           ! index1=1,nproma, index2=nlev_o3,
                               ! index3=1,nblks_c, index4=1,ntimes

    REAL(wp),POINTER::  &
      &   zf  (:),      &      !full levels of ozone gemometric height
      &   pfoz(:),      &      !full levels of of ozone pressure
      &   phoz(:)              !half levels of ozone pressure field 
    !
    ! *** radiation parameters ***
    REAL(wp), POINTER ::   &   !< aerosol optical thickness of black carbon    [ ]
      &  aer_bc(:,:,:)         ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes

    REAL(wp), POINTER ::   &   !< aerosol optical thickness of ambient aerosol [ ]
      &  aer_dust(:,:,:)       ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes

    REAL(wp), POINTER ::   &   !< aerosol optical thickness of particulate     [ ]
      &  aer_org(:,:,:)        !< organic_matter_ambient_aerosol             
                               ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes

    REAL(wp), POINTER ::   &   !< aerosol optical thickness of sulfate aerosol [ ]
      &  aer_so4(:,:,:)        ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes

    REAL(wp), POINTER ::   &   !< aerosol optical thickness of seasalt aerosol [ ]
      &  aer_ss(:,:,:)         ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes



    !
    ! *** vegetation parameters ***
    REAL(wp), POINTER ::   &   !< normalized difference vegetation index [ ]
      &  ndvi(:,:,:)           !< (monthly mean)
                               ! index1=1,nproma, index2=1,nblks_c, index3=1,ntimes

    REAL(wp), POINTER ::   &   !< (monthly) proportion of actual value/maximum 
      &  ndvi_mrat(:,:)        !< normalized differential vegetation index   [ ]
                               ! index1=1,nproma, index2=1,nblks_c

  END TYPE t_external_atmos_td



  !>
  !! ocean external data class
  !!
  !! ocean external data class
  !!
  TYPE :: t_external_ocean

    ! ocean topography <=> bathymetric height used in the ocean 
    ! cell centers and edges only
    !
    REAL(wp), POINTER ::   &  !<  bathymetric height at cell centers [m]
      &  bathymetry_c(:,:)    !  index1=1,nproma, index2=1,nblks_c

    REAL(wp), POINTER ::   &  !< topographic height at cell edges    [m]
      &  bathymetry_e(:,:)    !  index1=1,nproma, index2=1,nblks_e

  ! REAL(wp), POINTER ::   &  !< topographic height at cell vertices [m]
  !   &  bathymetry_v(:,:)    !  index1=1,nproma, index2=1,nblks_v

    ! *** Land-Sea-Mask ***
    INTEGER, POINTER  ::   &  !< land-sea-mask for cell centers          [ ]
      &  lsm_ctr_c(:,:)       !  index1=1,nproma, index2=1,nblks_c
    INTEGER, POINTER ::    &  !< land-sea-mask for cell edges
      &  lsm_ctr_e(:,:)       !  index1=1,nproma, index2=1,nblks_e

    ! OMIP forcing fluxes on cell centers. no_of_fluxes=3
    !
    REAL(wp), POINTER ::   &       !< omip monthly/daily mean forcing fluxes
      &  omip_forc_mon_c(:,:,:,:)  !  index1=nproma, index2=time, index3=nblks_c, index4=no_of_fluxes

  END TYPE t_external_ocean



!  !>
!  !! ocean external data class (time dependent)
!  !!
!  !! This data type contains auxiliary time dependent versions of
!  !! some external oceanic data fields already defined in external_ocean. These
!  !! fields will be used to store e.g. montly means from which interpolated external
!  !! data can be derived. The updated fields are copied into the time independent
!  !! counterparts which are defined in external_ocean.
!  !!
!  TYPE :: external_ocean_td
!  END TYPE external_ocean_td



  !>
  !! External data class including lists
  !!
  !! External data class including lists
  !!
  TYPE :: t_external_data

    TYPE(t_external_atmos)    :: atm
    TYPE(t_var_list)          :: atm_list

    TYPE(t_external_atmos_td) :: atm_td
    TYPE(t_var_list)          :: atm_td_list

    TYPE(t_external_ocean)    :: oce
    TYPE(t_var_list)          :: oce_list

!    TYPE(t_external_ocean_td) :: oce_td
!    TYPE(t_var_list), POINTER :: oce_td_list

  END TYPE t_external_data


  TYPE(t_external_data),TARGET, ALLOCATABLE :: &
    &  ext_data(:)  ! n_dom

!-------------------------------------------------------------------------

CONTAINS


  !-------------------------------------------------------------------------
  !>
  !! Init external data for atmosphere and ocean
  !!
  !! Init external data for atmosphere and ocean.
  !! 1. Build data structure, including field lists and 
  !!    memory allocation.
  !! 2. External data are read in from netCDF file or set analytically
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-16)
  !!
  SUBROUTINE init_ext_data (p_patch, p_int_state, ext_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_int_state), INTENT(IN)        :: p_int_state(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)

    INTEGER :: jg

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:init_ext_data'

    !-------------------------------------------------------------------------
    CALL message (TRIM(routine), 'Start')

    !-------------------------------------------------------------------------
    !  1.  inquire external files for their data structure
    !-------------------------------------------------------------------------

    ALLOCATE(nclass_lu(n_dom))
    ! Set default value for nclass_lu. Will be overwritten, if external data 
    ! are read from file
    nclass_lu(1:n_dom) = 1

    CALL inquire_external_files(p_patch)

    !------------------------------------------------------------------
    !  2.  construct external fields for the model
    !------------------------------------------------------------------

    ! top-level procedure for building data structures for 
    ! external data.
    CALL construct_ext_data(p_patch, ext_data)

    !-------------------------------------------------------------------------
    !  3.  read the data into the fields
    !-------------------------------------------------------------------------

    ! Check, whether external data should be read from file
    ! 

    !-------------------------------------------------------------------------
    ! Now for atmosphere only:
    !-------------------------------------------------------------------------

    IF (iequations/=ihs_ocean ) THEN

    !-------------------------------------------------------------------------
    !Ozone and aerosols
    !-------------------------------------------------------------------------

      IF((irad_o3 == io3_clim) .OR. (irad_o3 == io3_ape) ) THEN
        CALL message( TRIM(routine),'external ozone data required' )
        CALL read_ext_data_atm (p_patch, ext_data, nlev_o3)   ! read ozone
      ENDIF


    SELECT CASE(itopo)

    CASE(0) ! do not read external data
            ! topography from analytical functions

    !-------------------------------------------------------------------------
    ! surface/vegetation  parameter
    !-------------------------------------------------------------------------

      !
      ! initalize external data with meaningful data, in the case that they 
      ! are not read in from file.
      SELECT CASE ( iforcing )
      CASE ( inwp )
        DO jg = 1, n_dom
          ext_data(jg)%atm%emis_rad(:,:)    = zemiss_def ! longwave surface emissivity
          ext_data(jg)%atm%fr_land(:,:)     = 0._wp      ! land fraction
          ext_data(jg)%atm%fr_land_smt(:,:) = 0._wp      ! land fraction (smoothed)
          ext_data(jg)%atm%fr_glac_smt(:,:) = 0._wp      ! glacier fraction (smoothed)
          ext_data(jg)%atm%llsm_atm_c(:,:)  = .FALSE.    ! land-sea mask
          ext_data(jg)%atm%plcov_mx(:,:)    = 0.5_wp     ! plant cover
          ext_data(jg)%atm%lai_mx(:,:)      = 3._wp      ! max Leaf area index
          ext_data(jg)%atm%rootdp(:,:)      = 1._wp      ! root depth
          ext_data(jg)%atm%rsmin(:,:)       = 150._wp    ! minimal stomata resistence
          ext_data(jg)%atm%soiltyp(:,:)     = 8          ! soil type
          ext_data(jg)%atm%z0(:,:)          = 0.001_wp   ! roughness length
        END DO
      CASE ( iecham, ildf_echam)
        DO jg = 1, n_dom
          ext_data(jg)%atm%emis_rad(:,:)    = zemiss_def ! longwave surface emissivity
        END DO
      END SELECT

      CALL message( TRIM(routine),'Running with analytical topography' )

    CASE(1) ! read external data from netcdf dataset

      CALL message( TRIM(routine),'Start reading external data from file' )

      CALL read_ext_data_atm (p_patch, ext_data, nlev_o3)
      DO jg = 1, n_dom
        CALL smooth_topography (p_patch(jg), p_int_state(jg),  &
                                ext_data(jg)%atm%topography_c, &
                                ext_data(jg)%atm%topography_v, &
                                ext_data(jg)%atm%sso_stdh      )
      ENDDO

      CALL message( TRIM(routine),'Finished reading external data' )

    CASE DEFAULT

      CALL finish( TRIM(routine), 'topography selection not supported' )

    END SELECT

    ! read external ocean data from ocean gridfile (no itopo needed)
    ELSE
      CALL read_ext_data_oce (p_patch, ext_data)
    END IF


  END SUBROUTINE init_ext_data


  !-------------------------------------------------------------------------
  !>
  !! Top-level procedure for building external data structure
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-12)
  !!
  SUBROUTINE construct_ext_data (p_patch, ext_data)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_external_data),  INTENT(INOUT) :: ext_data(:)

    INTEGER :: jg

    CHARACTER(len=MAX_CHAR_LENGTH) :: listname

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:construct_ext_data'

!-------------------------------------------------------------------------


    CALL message (TRIM(routine), 'Construction of data structure for ' // &
      &                          'external data started')

    IF (iequations/=ihs_ocean) THEN  ! atmosphere model ---------------------

      ! Build external data list for constant-in-time fields for the atm model
      write(*,*) 'create new external data list for atmosphere'
      DO jg = 1, n_dom
        WRITE(listname,'(a,i2.2)') 'ext_data_atm_D',jg
        CALL new_ext_data_atm_list(p_patch(jg), ext_data(jg)%atm,       &
          &                        ext_data(jg)%atm_list, TRIM(listname))
      END DO

      ! Build external data list for time-dependent fields
      IF (iforcing > 1 ) THEN ! further distinction is made inside
      DO jg = 1, n_dom
          WRITE(listname,'(a,i2.2)') 'ext_data_atm_td_D',jg
          CALL new_ext_data_atm_td_list(p_patch(jg), ext_data(jg)%atm_td,       &
            &                           ext_data(jg)%atm_td_list, TRIM(listname))
      END DO
      END IF

    ELSE ! iequations==ihs_ocean ------------------------------------------

      write(*,*) 'create new external data list for ocean'
      ! Build external data list for constant-in-time fields for the ocean model
      DO jg = 1, n_dom
        WRITE(listname,'(a,i2.2)') 'ext_data_oce_D',jg
        CALL new_ext_data_oce_list(p_patch(jg), ext_data(jg)%oce,       &
          &                        ext_data(jg)%oce_list, TRIM(listname))
      END DO ! jg = 1,n_dom

      ! Build external data list for time-dependent fields
      ! ### to be done ###

    ENDIF ! atmosphere or ocean ------

    CALL message (TRIM(routine), 'Construction of data structure for ' // &
      &                          'external data finished')

  END SUBROUTINE construct_ext_data

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Allocation of atmospheric external data structure
  !!
  !! Allocation of atmospheric external data structure (constant in time 
  !! elements).
  !!
  !! Initialization of elements with zero.
  !!
  !! @par Revision History
  !! Initial release by Daniel Reinert (2011-05-03)
  !! Statements that assign initial value added by Hui Wan (MPI-M, 2011-05-30)
  !!
  SUBROUTINE new_ext_data_atm_list ( p_patch, p_ext_atm, p_ext_atm_list, &
    &                                listname)
!
    TYPE(t_patch), TARGET, INTENT(IN)     :: & !< current patch
      &  p_patch

    TYPE(t_external_atmos), INTENT(INOUT) :: & !< current external data structure
      &  p_ext_atm 

    TYPE(t_var_list) :: p_ext_atm_list !< current external data list

    CHARACTER(len=*), INTENT(IN)      :: & !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: jg

    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
      &        nblks_e, &    !< number of edge blocks to allocate
      &        nblks_v       !< number of vertex blocks to allocate

    INTEGER :: shape2d_c(2), shape2d_e(2), shape2d_v(2)
    INTEGER :: shape3d_sfc(3)

    INTEGER :: ientr         !< "entropy" of horizontal slice
    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nblks_v = p_patch%nblks_v

    ! get patch ID
    jg = p_patch%id
    ientr = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d_c  = (/ nproma, nblks_c /)
    shape2d_e  = (/ nproma, nblks_e /)
    shape2d_v  = (/ nproma, nblks_v /)
    shape3d_sfc= (/ nproma, nblks_c, nclass_lu(jg) /) 


    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_ext_atm_list, TRIM(listname), patch_id=p_patch%id )
    CALL default_var_list_settings( p_ext_atm_list,            &
                                  & lrestart=.FALSE.,          &
                                  & loutput=.FALSE.,           &
                                  & restart_type=FILETYPE_NC2  )


    ! topography height at cell center
    !
    ! topography_c  p_ext_atm%topography_c(nproma,nblks_c)
    cf_desc    = t_cf_var('surface_height', 'm', &
      &                   'geometric height of the earths surface above sea level')
    grib2_desc = t_grib2_var( 2, 0, 7, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'topography_c', p_ext_atm%topography_c,  &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc,          &
      &           grib2_desc, ldims=shape2d_c )


  SELECT CASE ( iequations )
  CASE ( inh_atmosphere )
    ! smoothed topography height at cell center
    !
    ! topography_smt_c  p_ext_atm%topography_smt_c(nproma,nblks_c)
    cf_desc    = t_cf_var('smoothed_surface_height', 'm', &
      &                   'smoothed geometric height of the earths surface above sea level')
    grib2_desc = t_grib2_var( 2, 0, 7, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'topography_smt_c', p_ext_atm%topography_smt_c, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc,                 &
      &           grib2_desc, ldims=shape2d_c )


    ! topography height at edge midpoint
    !
    ! topography_e  p_ext_atm%topography_e(nproma,nblks_e)
    cf_desc    = t_cf_var('surface_height', 'm', &
      &                   'geometric height of the earths surface above sea level')
    grib2_desc = t_grib2_var( 2, 0, 7, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_ext_atm_list, 'topography_e', p_ext_atm%topography_e, &
      &           GRID_UNSTRUCTURED_EDGE, ZAXIS_SURFACE, cf_desc,         &
      &           grib2_desc, ldims=shape2d_e)


    ! topography height at vertex
    !
    ! topography_v  p_ext_atm%topography_v(nproma,nblks_v)
    cf_desc    = t_cf_var('surface_height', 'm', &
      &                   'geometric height of the earths surface above sea level')
    grib2_desc = t_grib2_var( 2, 0, 7, ientr, GRID_REFERENCE, GRID_VERTEX)
    CALL add_var( p_ext_atm_list, 'topography_v', p_ext_atm%topography_v, &
      &           GRID_UNSTRUCTURED_VERT, ZAXIS_SURFACE, cf_desc,         &
      &           grib2_desc, ldims=shape2d_v )


    ! smoothed topography height at vertex
    !
    ! topography_smt_v  p_ext_atm%topography_smt_v(nproma,nblks_v)
    cf_desc    = t_cf_var('smoothed_surface_height', 'm', &
      &                   'smoothed geometric height of the earths surface above sea level')
    grib2_desc = t_grib2_var( 2, 0, 7, ientr, GRID_REFERENCE, GRID_VERTEX)
    CALL add_var( p_ext_atm_list, 'topography_smt_v', p_ext_atm%topography_smt_v, &
      &           GRID_UNSTRUCTURED_VERT, ZAXIS_SURFACE, cf_desc,                 &
      &           grib2_desc, ldims=shape2d_v )



    ! land sea mask for cells (LOGICAL)
    ! Note: Here "loutput" is set to .FALSE. since the output
    !       scheme operates on REAL model variables only and
    !       throws an error on this.
    !
    ! llsm_atm_c    p_ext_atm%llsm_atm_c(nproma,nblks_c)
    cf_desc    = t_cf_var('land_sea_mask_(cell)', '-', &
      &                   'land sea mask (cell)')
    grib2_desc = t_grib2_var( 2, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'llsm_atm_c', p_ext_atm%llsm_atm_c, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc,     &
      &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )

    ! land fraction
    !
    ! fr_land      p_ext_atm%fr_land(nproma,nblks_c)
    cf_desc    = t_cf_var('land_area_fraction', '-', 'Fraction land')
    grib2_desc = t_grib2_var( 2, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_land', p_ext_atm%fr_land,   &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
      &           grib2_desc, ldims=shape2d_c )


    ! glacier fraction
    !
    ! fr_glac      p_ext_atm%fr_glac(nproma,nblks_c)
    cf_desc    = t_cf_var('glacier_area_fraction', '-', 'Fraction glacier')
    grib2_desc = t_grib2_var( 2, 0, 192, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_glac', p_ext_atm%fr_glac,   &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
      &           grib2_desc, ldims=shape2d_c )


    ! maybe the next three (ice, fr_land_smt, fr_ice_smt)
    ! should be moved into corresponding if block

    ! sea Ice fraction
    !
    ! fr_ice       p_ext_atm%fr_ice(nproma,nblks_c)
    cf_desc    = t_cf_var('Sea_ice_fraction', '-', 'Sea ice fraction')
    grib2_desc = t_grib2_var( 10, 2, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_ice', p_ext_atm%fr_ice,     &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
      &           grib2_desc, ldims=shape2d_c )


    ! land fraction (smoothed)
    !
    ! fr_land_smt  p_ext_atm%fr_land_smt(nproma,nblks_c)
    cf_desc    = t_cf_var('land_area_fraction_(smoothed)', '-', &
      &                   'land area fraction (smoothed)')
    grib2_desc = t_grib2_var( 2, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_land_smt', p_ext_atm%fr_land_smt, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc,       &
      &           grib2_desc, ldims=shape2d_c )



    ! glacier area fraction (smoothed)
    !
    ! fr_glac_smt  p_ext_atm%fr_glac_smt(nproma,nblks_c)
    cf_desc    = t_cf_var('glacier_area_fraction_(smoothed)', '-', &
      &                   'glacier area fraction (smoothed)')
    grib2_desc = t_grib2_var( 2, 0, 192, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_glac_smt', p_ext_atm%fr_glac_smt, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc,       &
      &           grib2_desc, ldims=shape2d_c )



    ! sea Ice fraction (smoothed)
    !
    ! fr_ice_smt  p_ext_atm%fr_glac_smt(nproma,nblks_c)
    cf_desc    = t_cf_var('Sea_ice_fraction (smoothed)', '-', &
      &                   'Sea ice fraction (smoothed)')
    grib2_desc = t_grib2_var( 10, 2, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fr_ice_smt', p_ext_atm%fr_ice_smt, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc,     &
      &           grib2_desc, ldims=shape2d_c )


    ! geopotential (s)
    !
    ! fis          p_ext_atm%fis(nproma,nblks_c)
    cf_desc    = t_cf_var('Geopotential_(s)', 'm2 s-2', &
      &                   'Geopotential (s)')
    grib2_desc = t_grib2_var( 0, 3, 4, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_list, 'fis', p_ext_atm%fis,           &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
      &           grib2_desc, ldims=shape2d_c )



    SELECT CASE ( iforcing )
    CASE ( inwp )
      ! external parameter for NWP forcing

      ! roughness length
      !
      ! z0           p_ext_atm%z0(nproma,nblks_c)
      cf_desc    = t_cf_var('roughtness_length', 'm', 'roughtness length')
      grib2_desc = t_grib2_var( 2, 0, 1, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'z0', p_ext_atm%z0,             &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c )


      ! fraction lake
      !
      ! fr_lake      p_ext_atm%fr_lake(nproma,nblks_c)
      cf_desc    = t_cf_var('fraction_lake', '-', 'fraction lake')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'fr_lake', p_ext_atm%fr_lake,   &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c )


      ! lake depth
      !
      ! depth_lk     p_ext_atm%depth_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('lake_depth', '-', 'lake depth')
      grib2_desc = t_grib2_var( 192, 228, 7, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'depth_lk', p_ext_atm%depth_lk, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c )



      !--------------------------------
      ! sub-gridscale orography
      !--------------------------------

      ! Standard deviation of sub-grid scale orography
      !
      ! sso_stdh     p_ext_atm%sso_stdh(nproma,nblks_c)
      cf_desc    = t_cf_var('standard_deviation_of_height', 'm',&
        &                   'Standard deviation of sub-grid scale orography')
      grib2_desc = t_grib2_var( 0, 3, 20, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_stdh', p_ext_atm%sso_stdh, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c )



      ! Anisotropy of sub-gridscale orography
      !
      ! sso_gamma    p_ext_atm%sso_gamma(nproma,nblks_c)
      cf_desc    = t_cf_var('anisotropy_factor', '-',&
        &                   'Anisotropy of sub-gridscale orography')
      grib2_desc = t_grib2_var( 0, 3, 20, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_gamma', p_ext_atm%sso_gamma, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc,   &
        &           grib2_desc, ldims=shape2d_c )



      ! Angle of sub-gridscale orography
      !
      ! sso_theta    p_ext_atm%sso_theta(nproma,nblks_c)
      cf_desc    = t_cf_var('angle_of_principal_axis', 'radians',&
        &                   'Angle of sub-gridscale orography')
      grib2_desc = t_grib2_var( 0, 3, 21, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_theta', p_ext_atm%sso_theta, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc,   &
        &           grib2_desc, ldims=shape2d_c )



      ! Slope of sub-gridscale orography
      !
      ! sso_sigma    p_ext_atm%sso_sigma(nproma,nblks_c)
      cf_desc    = t_cf_var('slope_of_terrain', '-',&
        &                   'Slope of sub-gridscale orography')
      grib2_desc = t_grib2_var( 0, 3, 22, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'sso_sigma', p_ext_atm%sso_sigma, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc,   &
        &           grib2_desc, ldims=shape2d_c )





      !--------------------------------
      ! vegetation parameters
      !--------------------------------

      ! Plant covering degree in the vegetation phase
      !
      ! plcov_mx     p_ext_atm%plcov_mx(nproma,nblks_c)
      cf_desc    = t_cf_var('vegetation_area_fraction_vegetation_period', '-',&
        &                   'Plant covering degree in the vegetation phase')
      grib2_desc = t_grib2_var( 2, 0, 4, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'plcov_mx', p_ext_atm%plcov_mx, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c )



      ! Max Leaf area index
      !
      ! lai_mx       p_ext_atm%lai_mx(nproma,nblks_c)
      cf_desc    = t_cf_var('leaf_area_index_vegetation_period', '-',&
        &                   'Leaf Area Index Maximum')
      grib2_desc = t_grib2_var( 2, 0, 28, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'lai_mx', p_ext_atm%lai_mx,     &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c )



      ! root depth of vegetation
      !
      ! rootdp      p_ext_atm%rootdp(nproma,nblks_c)
      cf_desc    = t_cf_var('root_depth_of_vegetation', 'm',&
        &                   'root depth of vegetation')
      grib2_desc = t_grib2_var( 2, 0, 32, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'rootdp', p_ext_atm%rootdp,     &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c )



      ! evergreen forest
      !
      ! for_e        p_ext_atm%for_e(nproma,nblks_c)
      cf_desc    = t_cf_var('fraction_of_evergreen_forest_cover', '-',&
        &                   'Fraction of evergreen forest')
      grib2_desc = t_grib2_var( 2, 0, 29, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'for_e', p_ext_atm%for_e,       &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c )



      ! deciduous forest
      !
      ! for_d     p_ext_atm%for_d(nproma,nblks_c)
      cf_desc    = t_cf_var('fraction_of_deciduous_forest_cover', '-',&
        &                   'Fraction of deciduous forest')
      grib2_desc = t_grib2_var( 2, 0, 30, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'for_d', p_ext_atm%for_d,       &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c )



      ! urban area fraction
      !
      ! urban        p_ext_atm%urban(nproma,nblks_c)
      cf_desc    = t_cf_var('fraction_of_urban_areas', '-',&
        &                   'urban area fraction')
      grib2_desc = t_grib2_var( 2, 0, 30, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'urban', p_ext_atm%urban,       &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c )


      ! Minimal stomata resistence
      !
      ! rsmin        p_ext_atm%rsmin(nproma,nblks_c)
      cf_desc    = t_cf_var('RSMIN', 's m-1', 'Minimal stomata resistence')
      grib2_desc = t_grib2_var( 2, 0, 16, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'rsmin', p_ext_atm%rsmin,       &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c )


      ! NDVI yearly maximum
      !
      ! ndvi_max        p_ext_atm%ndvi_max(nproma,nblks_c)
      cf_desc    = t_cf_var('normalized_difference_vegetation_index', '-', &
        &                   'NDVI yearly maximum')
      grib2_desc = t_grib2_var( 2, 0, 31, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'ndvi_max', p_ext_atm%ndvi_max, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c )



      !--------------------------------
      ! soil parameters
      !--------------------------------

      ! soil type
      !
      ! soiltyp      p_ext_atm%soiltyp(nproma,nblks_c)
      cf_desc    = t_cf_var('soil_type', '-','soil type')
      grib2_desc = t_grib2_var( 2, 3, 0, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'soiltyp', p_ext_atm%soiltyp,   &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c, loutput=.FALSE. )


      ! Climat. temperature
      !
      ! t_cl         p_ext_atm%t_cl(nproma,nblks_c)
      cf_desc    = t_cf_var('soil_temperature', 'K',                  &
        &                   'CRU near surface temperature climatology')
      grib2_desc = t_grib2_var( 2, 3, 18, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 't_cl', p_ext_atm%t_cl,         &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c )


      ! longwave surface emissivity
      !
      ! emis_rad     p_ext_atm%emis_rad(nproma,nblks_c)
      cf_desc    = t_cf_var('emis_rad', '-', 'longwave surface emissivity')
      grib2_desc = t_grib2_var( 2, 3, 196, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'emis_rad', p_ext_atm%emis_rad, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c)


      ! landuse class fraction
      !
      ! lu_class_fraction    p_ext_atm%lu_class_fraction(nproma,nblks_c,nclass_lu)
      cf_desc    = t_cf_var('lu_class_fraction', '-', 'landuse class fraction')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'lu_class_fraction', p_ext_atm%lu_class_fraction, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape3d_sfc )

    CASE ( iecham, ildf_echam )

      ! longwave surface emissivity
      !
      ! emis_rad     p_ext_atm%emis_rad(nproma,nblks_c)
      cf_desc    = t_cf_var('emis_rad', '-', 'longwave surface emissivity')
      grib2_desc = t_grib2_var( 2, 3, 196, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'emis_rad', p_ext_atm%emis_rad, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c)

    END SELECT ! iforcing

  CASE ( ihs_atm_temp, ihs_atm_theta )

    SELECT CASE ( iforcing )
    CASE ( iecham, ildf_echam)
      ! longwave surface emissivity
      !
      ! emis_rad     p_ext_atm%emis_rad(nproma,nblks_c)
      cf_desc    = t_cf_var('emis_rad', '-', 'longwave surface emissivity')
      grib2_desc = t_grib2_var( 2, 3, 196, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_list, 'emis_rad', p_ext_atm%emis_rad, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, &
        &           grib2_desc, ldims=shape2d_c)     
    END SELECT
        
  END SELECT ! iequations


  END SUBROUTINE new_ext_data_atm_list



  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Allocation of atmospheric external data structure (time dependent)
  !!
  !! Allocation of atmospheric external data structure (time dependent  
  !! elements).
  !!
  !! Initialization of elements with zero.
  !!
  !! @par Revision History
  !! Initial release by Daniel Reinert (2011-05-03)
  !!
  SUBROUTINE new_ext_data_atm_td_list ( p_patch, p_ext_atm_td, &
    &                               p_ext_atm_td_list, listname)
!
    TYPE(t_patch), TARGET, INTENT(IN)     :: & !< current patch
      &  p_patch

    TYPE(t_external_atmos_td), INTENT(INOUT) :: & !< current external data structure
      &  p_ext_atm_td 

    TYPE(t_var_list) :: p_ext_atm_td_list  !< current external data list

    CHARACTER(len=*), INTENT(IN)      :: & !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c      !< number of cell blocks to allocate

    INTEGER :: shape3d_c(3)
    INTEGER :: shape3d_ape(3)
    INTEGER :: shape4d_c(4)

    INTEGER :: ientr         !< "entropy" of horizontal slice
    INTEGER :: ntimes        !< number of time slices
    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c

    ientr  = 16   ! "entropy" of horizontal slice
    ntimes = 12   ! number of time slices

    ! predefined array shapes
    shape3d_c = (/ nproma, nblks_c, ntimes /)
    shape4d_c = (/ nproma, nlev_o3, nblks_c, nmonths /) 
    shape3d_ape = (/ nproma, nlev_o3, nblks_c/) 

    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_ext_atm_td_list, TRIM(listname), patch_id=p_patch%id )
    CALL default_var_list_settings( p_ext_atm_td_list,         &
                                  & lrestart=.FALSE.,          &
                                  & loutput=.FALSE.,           &
                                  & restart_type=FILETYPE_NC2  )


    !--------------------------------
    ! radiation parameters
    !--------------------------------


    ! ozone on pressure levels
    ! ATTENTION: a GRIB2 number will go to 
    ! the ozone mass mixing ratio...
    !
    IF(irad_o3 == io3_clim .OR. irad_o3 == io3_ape) THEN 

      WRITE(0,*) 'generate ext ozone field'

      ! o3  main height level from read-in file
      cf_desc    = t_cf_var('O3_zf', 'm',   &
        &                   'ozone geometric height level')
      grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_td_list, 'O3_zf', p_ext_atm_td%zf, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, &
        &           grib2_desc, ldims=(/nlev_o3/) )

      ! o3  main pressure level from read-in file
      cf_desc    = t_cf_var('O3_pf', 'Pa',   &
        &                   'ozone main pressure level')
      grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_td_list, 'O3_pf', p_ext_atm_td%pfoz, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE, cf_desc, &
        &           grib2_desc, ldims=(/nlev_o3/) )

      ! o3  intermediate pressure level
      cf_desc    = t_cf_var('O3_ph', 'Pa',   &
        &                   'ozone intermediate pressure level')
      grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_atm_td_list, 'O3_ph', p_ext_atm_td%phoz, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE, cf_desc, &
        &           grib2_desc, ldims=(/nlev_o3+1/) )

        ! o3       p_ext_atm_td%o3(nproma,nlev_o3,nblks_c,nmonths)
        cf_desc    = t_cf_var('O3', TRIM(o3unit),   &
          &                   'mole_fraction_of_ozone_in_air')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_ext_atm_td_list, 'O3', p_ext_atm_td%O3, &
          &           GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE, cf_desc, &
          &           grib2_desc, ldims=shape4d_c )

    END IF ! irad_o3

    IF(iforcing == inwp) THEN

    ! Black carbon aerosol
    !
    ! aer_bc       p_ext_atm%aer_bc(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aerosol optical thickness of black carbon', '-',   &
      &                   'atmosphere_absorption_optical_thickness_due_to_' //&
      &                   'black_carbon_ambient_aerosol')
    grib2_desc = t_grib2_var( 0, 13, 195, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_bc', p_ext_atm_td%aer_bc, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc,  &
      &            grib2_desc, ldims=shape3d_c )


    ! Dust aerosol
    !
    ! aer_dust     p_ext_atm%aer_dust(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aot_dust', '-', &
      &                   'atmosphere absorption optical thickness due '//  &
      &                   'to dust ambient aerosol')
    grib2_desc = t_grib2_var( 0, 13, 193, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_dust', p_ext_atm_td%aer_dust, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
      &           ldims=shape3d_c )


    ! Organic aerosol
    !
    ! aer_org      p_ext_atm%aer_org(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aot_org', '-', &
      &                   'atmosphere absorption optical thickness due '//  &
      &                   'to particulate organic matter ambient aerosol')
    grib2_desc = t_grib2_var( 0, 13, 194, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_org', p_ext_atm_td%aer_org, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c )


    ! Sulfate aerosol
    !
    ! aer_so4      p_ext_atm%aer_so4(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aot_so4', '-', &
      &                   'atmosphere absorption optical thickness due '//  &
      &                   'to sulfate_ambient_aerosol')
    grib2_desc = t_grib2_var( 0, 13, 192, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_so4', p_ext_atm_td%aer_so4, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c )


    ! Seasalt aerosol
    !
    ! aer_ss       p_ext_atm%aer_ss(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('aot_ss', '-', &
      &                   'atmosphere absorption optical thickness due '//  &
      &                   'to seasalt_ambient_aerosol')
    grib2_desc = t_grib2_var( 0, 13, 196, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'aer_ss', p_ext_atm_td%aer_ss, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c )


    !--------------------------------
    ! vegetation parameters
    !--------------------------------

    ! monthly mean normalized difference vegetation index
    !
    ! ndvi         p_ext_atm%ndvi(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('normalized_difference_vegetation_index', '-', &
      &                   'monthly mean NDVI')
    grib2_desc = t_grib2_var( 2, 0, 217, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'ndvi', p_ext_atm_td%ndvi, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c )


    ! (monthly) proportion of actual value/maximum NDVI
    !
    ! ndvi_mrat     p_ext_atm%ndvi_mrat(nproma,nblks_c,ntimes)
    cf_desc    = t_cf_var('normalized_difference_vegetation_index', '-', &
      &                   '(monthly) proportion of actual value/maximum ' // &
      &                   'normalized differential vegetation index')
    grib2_desc = t_grib2_var( 2, 0, 192, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_atm_td_list, 'ndvi_mrat', p_ext_atm_td%ndvi_mrat, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,&
      &           ldims=shape3d_c)

    ENDIF ! inwp

  END SUBROUTINE new_ext_data_atm_td_list



  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Allocation of oceanic external data structure
  !!
  !! Allocation of oceanic external data structure used for elements that are
  !! stationary in time.
  !!
  !! Initialization of elements with zero.
  !!
  !! @par Revision History
  !! Initial release by Daniel Reinert (2011-06-24)
  !!
  SUBROUTINE new_ext_data_oce_list ( p_patch, p_ext_oce, p_ext_oce_list, &
    &                                listname)
!
    TYPE(t_patch), TARGET, INTENT(IN)   :: & !< current patch
      &  p_patch

    TYPE(t_external_ocean), INTENT(INOUT) :: & !< current external data structure
      &  p_ext_oce 

    TYPE(t_var_list) :: p_ext_oce_list !< current external data list

    CHARACTER(len=*), INTENT(IN)  :: & !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
      &        nblks_e       !< number of edge blocks to allocate

    INTEGER :: shape2d_c(2), shape2d_e(2), shape4d_c(4)
    INTEGER :: idim_omip

    INTEGER :: ientr         !< "entropy" of horizontal slice

    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e


    ientr = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d_c = (/ nproma, nblks_c /)
    shape2d_e = (/ nproma, nblks_e /)

    ! OMIP or other flux forcing data on cell centers: 3, 5 or 15 variables, iforc_len data sets
    IF (iforc_omip == 1 ) idim_omip =  3
    IF (iforc_omip == 2 ) idim_omip = 15
    IF (iforc_omip == 3 ) idim_omip =  5
    IF (iforc_omip == 4 ) idim_omip =  9
    IF (iforc_omip == 5 ) idim_omip = 15
    shape4d_c = (/ nproma, iforc_len, nblks_c, idim_omip /)

    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_ext_oce_list, TRIM(listname), patch_id=p_patch%id )
    CALL default_var_list_settings( p_ext_oce_list,            &
                                  & lrestart=.FALSE.,          &
                                  & restart_type=FILETYPE_NC2, &
                                  & model_type='oce'  )

    ! bathymetric height at cell center
    !
    ! bathymetry_c  p_ext_oce%bathymetry_c(nproma,nblks_c)
    cf_desc    = t_cf_var('Model bathymetry at cell center', 'm', &
      &                   'Model bathymetry')
    grib2_desc = t_grib2_var( 192, 140, 219, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ext_oce_list, 'bathymetry_c', p_ext_oce%bathymetry_c,      &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )


    ! bathymetric height at cell edge
    !
    ! bathymetry_e  p_ext_oce%bathymetry_e(nproma,nblks_e)
    cf_desc    = t_cf_var('Model bathymetry at edge', 'm', &
      &                   'Model bathymetry')
    grib2_desc = t_grib2_var( 192, 140, 219, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_ext_oce_list, 'bathymetry_e', p_ext_oce%bathymetry_e,      &
      &           GRID_UNSTRUCTURED_EDGE, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_e )

    ! ocean land-sea-mask at surface on cell centers
    !
    ! lsm_ctr_c  p_ext_oce%lsm_ctr_c(nproma,nblks_c)
    cf_desc    = t_cf_var('Ocean model land-sea-mask at cell center', '-2/-1/1/2', &
      &                   'Ocean model land-sea-mask')
    grib2_desc = t_grib2_var( 192, 140, 219, ientr, GRID_REFERENCE, GRID_CELL)
    !#slo-2011-08-08# does not compile yet?
    CALL add_var( p_ext_oce_list, 'lsm_ctr_c', p_ext_oce%lsm_ctr_c, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_c )

    ! ocean land-sea-mask at surface on cell edge
    !
    cf_desc    = t_cf_var('Ocean model land-sea-mask at cell edge', '-2/0/2', &
      &                   'Ocean model land-sea-mask')
    grib2_desc = t_grib2_var( 192, 140, 219, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_ext_oce_list, 'lsm_ctr_e', p_ext_oce%lsm_ctr_e,      &
      &           GRID_UNSTRUCTURED_EDGE, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape2d_e )

    ! omip forcing data on cell edge
    !
    IF (iforc_oce == 12) THEN
      cf_desc    = t_cf_var('Ocean model OMIP forcing data at cell edge', 'Pa, K', &
        &                   'OMIP forcing data')
      grib2_desc = t_grib2_var( 192, 140, 219, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_ext_oce_list, 'omip_forc_mon_c', p_ext_oce%omip_forc_mon_c,  &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, ldims=shape4d_c )
    END IF

  END SUBROUTINE new_ext_data_oce_list



  !-------------------------------------------------------------------------
  !>
  !! Destruct external data data structure and lists
  !!
  !! Destruct external data data structure and lists
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2011-05-04)
  !!
  SUBROUTINE destruct_ext_data

    INTEGER :: jg, errstat
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_ext_data:destruct_ext_data'
    !-------------------------------------------------------------------------

    CALL message (TRIM(routine), 'Destruction of data structure for ' // &
      &                          'external data started')

    IF (iequations/=ihs_ocean) THEN  ! atmosphere model ------------------

      DO jg = 1,n_dom
        ! Delete list of constant in time atmospheric elements
        CALL delete_var_list( ext_data(jg)%atm_list )
      ENDDO

      IF (iforcing > 1 ) THEN
      DO jg = 1,n_dom
        ! Delete list of time-dependent atmospheric elements
        CALL delete_var_list( ext_data(jg)%atm_td_list )
      ENDDO
      END IF

    ELSE ! iequations==ihs_ocean ------------------------------------------

      DO jg = 1,n_dom
        ! Delete list of constant in time oceanic elements
        CALL delete_var_list( ext_data(jg)%oce_list )
      ENDDO

        ! Delete list of time-dependent oceanic elements
        ! ### to be added if necessary ###

    END IF

    DEALLOCATE(nclass_lu, STAT=errstat)
    IF (errstat /= 0)  &
      CALL finish (TRIM(routine), 'Error in DEALLOCATE operation!')

    CALL message (TRIM(routine), 'Destruction of data structure for ' // &
      &                          'external data finished')

  END SUBROUTINE destruct_ext_data



 !-------------------------------------------------------------------------
 !-------------------------------------------------------------------------

  SUBROUTINE inquire_external_files(p_patch)

    !-------------------------------------------------------
    !
    ! open netcdf files and investigate the data structure  
    ! of the external parameters
    !
    !-------------------------------------------------------

    TYPE(t_patch), INTENT(IN)      :: p_patch(:)

    INTEGER :: jg, mpi_comm
    INTEGER :: no_cells, no_verts
    INTEGER :: ncid, dimid

    LOGICAL :: l_exist

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data: inquire_external_files'

    CHARACTER(filename_max) :: extpar_file !< file name for reading in
    CHARACTER(filename_max) :: ozone_file  !< file name for reading in

!--------------------------------------------------------------------------

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    DO jg= 1,n_dom

      !------------------------------------------------!
      ! 1. Check validity of external parameter file   !
      !------------------------------------------------!
      IF (itopo == 1 .AND. iequations /=ihs_ocean ) THEN

        IF( my_process_is_stdio()) THEN
          !
          ! generate file name
          extpar_file = generate_filename(extpar_filename,                   &
            &                             model_base_dir,                    &
            &                             TRIM(p_patch(jg)%grid_filename))
          CALL message("read_ext_data_atm, extpar_file=",extpar_file)
        
          INQUIRE (FILE=extpar_file, EXIST=l_exist)
          IF (.NOT.l_exist) THEN
            CALL finish(TRIM(routine),'external data file is not found.')
          ENDIF

          !
          ! open file
          !
          CALL nf(nf_open(TRIM(extpar_file), NF_NOWRITE, ncid))

          !
          ! get number of cells and vertices
          !
          CALL nf(nf_inq_dimid(ncid, 'cell', dimid))
          IF (p_patch(jg)%cell_type == 3) THEN ! triangular grid
            CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))
          ELSEIF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid
            CALL nf(nf_inq_dimlen(ncid, dimid, no_verts))
          ENDIF

          CALL nf(nf_inq_dimid(ncid, 'vertex', dimid))
          IF (p_patch(jg)%cell_type == 3) THEN ! triangular grid
            CALL nf(nf_inq_dimlen(ncid, dimid, no_verts))
          ELSEIF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid
            CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))
          ENDIF

          !
          ! check the number of cells and verts
          !
          IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
            CALL finish(TRIM(ROUTINE),&
            & 'Number of patch cells and cells in extpar file do not match.')
          ENDIF
          IF(p_patch(jg)%n_patch_verts_g /= no_verts) THEN
            CALL finish(TRIM(ROUTINE),&
            & 'Number of patch verts and verts in extpar file do not match.')
          ENDIF

          !
          ! get the number of landuse classes
          !
          CALL nf(nf_inq_dimid (ncid, 'nclass_lu', dimid))
          CALL nf(nf_inq_dimlen(ncid, dimid, nclass_lu(jg)))

          !
          ! close file
          !
          CALL nf(nf_close(ncid))

        ENDIF ! my_process_is_stdio()

        ! broadcast from nclass_lu I-Pe to WORK Pes
        CALL p_bcast(nclass_lu(jg), p_io, mpi_comm)

      ENDIF


      !------------------------------------------------!
      ! 2. Check validity of ozone file                !
      !------------------------------------------------!

      ! default values for nlev_o3 and nmonths
      nlev_o3 = 1
      nmonths   = 1

      O3 : IF((irad_o3 == io3_clim) .OR. (irad_o3 == io3_ape )) THEN

        IF(p_test_run) THEN
          mpi_comm = p_comm_work_test
        ELSE
          mpi_comm = p_comm_work
        ENDIF

        IF(irad_o3 == io3_ape ) THEN
          levelname = 'level'
          zlevelname = 'zlev'
          cellname  = 'ncells'
          o3name    = 'O3'
          o3unit    = 'g/g'
        ELSE ! o3_clim
          levelname= 'plev'
          cellname = 'cell'
          o3name   = 'O3'
          o3unit   = 'g/g' !this unit ozon will have after being read out
                           ! and converted from ppmv
        ENDIF

        IF_IO : IF(my_process_is_stdio()) THEN
          
          WRITE(ozone_file,'(a,i2.2,a)') 'o3_icon_DOM',jg,'.nc'

          ! Note resolution assignment is done per script by symbolic links

          INQUIRE (FILE=ozone_file, EXIST=l_exist)
          IF (.NOT.l_exist) THEN
            WRITE(0,*) 'DOMAIN=',jg
            CALL finish(TRIM(routine),'ozone file of domain is not found.')
          ENDIF

          !
          ! open file
          !
          CALL nf(nf_open(TRIM(ozone_file), NF_NOWRITE, ncid))

          WRITE(0,*)'open ozone file'
          ! get number of cells in triangles and hexagons
          !

          !triangles
          IF (p_patch(jg)%cell_type == 3) THEN ! triangular grid
            CALL nf(nf_inq_dimid (ncid, TRIM(cellname), dimid))
            CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))
  
            WRITE(0,*)'number of cells are', no_cells
          !
          ! check the number of cells and verts
          !
            IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
            CALL finish(TRIM(ROUTINE),&
              & 'Number of patch cells and cells in ozone file do not match.')
          ENDIF
        ENDIF
       
          !hexagons
          IF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid
            CALL nf(nf_inq_dimid (ncid, TRIM(cellname), dimid))
            CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))

            WRITE(0,*)'number of hexcells_o3 are', no_cells
            WRITE(0,*)'number of hexverts are', p_patch(jg)%n_patch_verts_g
            WRITE(0,*)'number of hexcellsmo are', p_patch(jg)%n_patch_cells_g 
          !
          ! check the number of cells and verts
          !
            IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
!          IF(p_patch(jg)%n_patch_verts_g /= no_cells) THEN
            CALL finish(TRIM(ROUTINE),&
            & 'Number of patch cells and cells in ozone file do not match.')
          ENDIF
          ENDIF

          ! check the time structure
          CALL nf(nf_inq_dimid (ncid, 'time', dimid))
          CALL nf(nf_inq_dimlen(ncid, dimid, nmonths))

          WRITE(message_text,'(A,I4)')  &
            & 'Number of months in ozone file = ', &
            & nmonths
          CALL message(TRIM(ROUTINE),message_text)

          ! check the vertical structure

!          SELECT CASE (iequations)
!          CASE(ihs_atm_temp,ihs_atm_theta)

            CALL nf(nf_inq_dimid (ncid,TRIM(levelname), dimid))
            CALL nf(nf_inq_dimlen(ncid, dimid, nlev_o3))

            WRITE(message_text,'(A,I4)')  &
              & 'Number of pressure levels in ozone file = ', nlev_o3
            CALL message(TRIM(ROUTINE),message_text)

!          CASE(inh_atmosphere)
!            CALL nf(nf_inq_dimid (ncid,TRIM(zlevelname), dimid))
!            CALL nf(nf_inq_dimlen(ncid, dimid, nlev_o3))
!            WRITE(message_text,'(A,I4)')  &
!              & 'Number of height levels in ozone file = ', nlev_o3
!            CALL message(TRIM(ROUTINE),message_text)
!          END SELECT
  
          !
          ! close file
          !
          CALL nf(nf_close(ncid))

        END IF IF_IO ! pe

        CALL p_bcast(nlev_o3, p_io, mpi_comm)      
        CALL p_bcast(nmonths,   p_io, mpi_comm)      

      END IF O3 !o3

    ENDDO ! ndom

  END SUBROUTINE inquire_external_files

  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Read atmospheric external data
  !!
  !! Read atmospheric external data from netcdf
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2010-07-14)
  !!
  SUBROUTINE read_ext_data_atm (p_patch, ext_data, nlev_o3)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)
    INTEGER, INTENT(IN)                  :: nlev_o3

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:read_ext_data_atm'

    CHARACTER(filename_max) :: extpar_file !< file name for reading in
    CHARACTER(filename_max) :: ozone_file  !< file name for reading in

    INTEGER :: jg, jc, jb, i, mpi_comm
    INTEGER :: i_lev,jk
    INTEGER :: ncid, varid

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk   !> blocks
    INTEGER :: i_startidx, i_endidx   !< slices
    INTEGER :: i_nchdom               !< domain index

    REAL(wp):: zdummy_o3lev(nlev_o3) ! will be used for pressure and height levels

    !----------------------------------------------------------------------

    IF(itopo == 1 ) THEN
      DO jg = 1,n_dom

        i_lev = p_patch(jg)%level

        IF(my_process_is_stdio()) THEN

          !
          ! generate file name and open file
          extpar_file = generate_filename(extpar_filename,                  &
            &                             model_base_dir,                   &
            &                             TRIM(p_patch(jg)%grid_filename))
          CALL nf(nf_open(TRIM(extpar_file), NF_NOWRITE, ncid))

        ENDIF


        !-------------------------------------------------------
        !
        ! Read topography for triangle centers and vertices
        !
        !-------------------------------------------------------

        !
        ! topography
        !
        IF (p_patch(jg)%cell_type == 3) THEN     ! triangular grid

          ! triangle center
          CALL read_netcdf_data (ncid, 'topography_c', p_patch(jg)%n_patch_cells_g,       &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%topography_c)

          ! triangle vertex
          CALL read_netcdf_data (ncid, 'topography_v', p_patch(jg)%n_patch_verts_g,       &
            &                     p_patch(jg)%n_patch_verts, p_patch(jg)%verts%glb_index, &
            &                     ext_data(jg)%atm%topography_v)

        ELSEIF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid

          ! As extpar "knows" only the triangular grid, cells and vertices need to be switched here
          ! triangle center
          CALL read_netcdf_data (ncid, 'topography_c', p_patch(jg)%n_patch_verts_g,       &
            &                     p_patch(jg)%n_patch_verts, p_patch(jg)%verts%glb_index, &
            &                     ext_data(jg)%atm%topography_v)

          ! triangle vertex
          CALL read_netcdf_data (ncid, 'topography_v', p_patch(jg)%n_patch_cells_g,       &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%topography_c)

        ENDIF


        !
        ! other external parameters on triangular grid
        !
        IF (p_patch(jg)%cell_type == 3) THEN     ! triangular grid

          CALL read_netcdf_data (ncid, 'FR_LAND', p_patch(jg)%n_patch_cells_g,              &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index,   &
            &                     ext_data(jg)%atm%fr_land)

          CALL read_netcdf_data (ncid, 'ICE', p_patch(jg)%n_patch_cells_g,                &
            &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
            &                     ext_data(jg)%atm%fr_ice)


          SELECT CASE ( iforcing )
          CASE ( inwp )

            CALL read_netcdf_data (ncid, 'PLCOV_MX', p_patch(jg)%n_patch_cells_g,           &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%plcov_mx)

            CALL read_netcdf_data (ncid, 'LAI_MX', p_patch(jg)%n_patch_cells_g,             &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%lai_mx)

            CALL read_netcdf_data (ncid, 'ROOTDP', p_patch(jg)%n_patch_cells_g,             &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%rootdp)

            CALL read_netcdf_data (ncid, 'RSMIN', p_patch(jg)%n_patch_cells_g,             &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%rsmin)

            CALL read_netcdf_data (ncid, 'FOR_D', p_patch(jg)%n_patch_cells_g,              &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%for_d)

            CALL read_netcdf_data (ncid, 'FOR_E', p_patch(jg)%n_patch_cells_g,              &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%for_e)

            CALL read_netcdf_data (ncid, 'URBAN', p_patch(jg)%n_patch_cells_g,              &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%urban)

            CALL read_netcdf_data (ncid, 'Z0', p_patch(jg)%n_patch_cells_g,                 &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%z0)

            CALL read_netcdf_data (ncid, 'NDVI_MAX', p_patch(jg)%n_patch_cells_g,           &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%ndvi_max)

            CALL read_netcdf_data (ncid, 'SOILTYP', p_patch(jg)%n_patch_cells_g,        &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%soiltyp)

            CALL read_netcdf_lu (ncid, 'LU_CLASS_FRACTION', p_patch(jg)%n_patch_cells_g,  &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     nclass_lu(jg), ext_data(jg)%atm%lu_class_fraction )


            IF ( irad_aero == 6 ) THEN

              CALL read_netcdf_data (ncid, 12, 'AER_SS',           &
                &                    p_patch(jg)%n_patch_cells_g,  &
                &                    p_patch(jg)%n_patch_cells,    &
                &                    p_patch(jg)%cells%glb_index,  & 
                &                    ext_data(jg)%atm_td%aer_ss)

              CALL read_netcdf_data (ncid, 12, 'AER_DUST',         &
                &                    p_patch(jg)%n_patch_cells_g,  &
                &                    p_patch(jg)%n_patch_cells,    &
                &                    p_patch(jg)%cells%glb_index,  & 
                &                    ext_data(jg)%atm_td%aer_dust)

              CALL read_netcdf_data (ncid, 12, 'AER_ORG',          &
                &                    p_patch(jg)%n_patch_cells_g,  &
                &                    p_patch(jg)%n_patch_cells,    &
                &                    p_patch(jg)%cells%glb_index,  & 
                &                    ext_data(jg)%atm_td%aer_org)

              CALL read_netcdf_data (ncid, 12, 'AER_SO4',          &
                &                    p_patch(jg)%n_patch_cells_g,  &
                &                    p_patch(jg)%n_patch_cells,    &
                &                    p_patch(jg)%cells%glb_index,  & 
                &                    ext_data(jg)%atm_td%aer_so4)

              CALL read_netcdf_data (ncid, 12, 'AER_BC',           &
                &                    p_patch(jg)%n_patch_cells_g,  &
                &                    p_patch(jg)%n_patch_cells,    &
                &                    p_patch(jg)%cells%glb_index,  & 
                &                    ext_data(jg)%atm_td%aer_bc)

            ENDIF
           
            IF ( l_emiss ) THEN
              CALL read_netcdf_data (ncid, 'EMIS_RAD', p_patch(jg)%n_patch_cells_g,           &
                &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
                &                     ext_data(jg)%atm%emis_rad)
            ELSE
              ext_data(jg)%atm%emis_rad(:,:)= zemiss_def
            ENDIF

            CALL read_netcdf_data (ncid, 'T_CL', p_patch(jg)%n_patch_cells_g,               &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%t_cl)

            CALL read_netcdf_data (ncid, 'SSO_STDH', p_patch(jg)%n_patch_cells_g,           &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%sso_stdh)

            CALL read_netcdf_data (ncid, 'SSO_THETA', p_patch(jg)%n_patch_cells_g,          &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%sso_theta)

            CALL read_netcdf_data (ncid, 'SSO_GAMMA', p_patch(jg)%n_patch_cells_g,          &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%sso_gamma)

            CALL read_netcdf_data (ncid, 'SSO_SIGMA', p_patch(jg)%n_patch_cells_g,          &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%sso_sigma)

            CALL read_netcdf_data (ncid, 'FR_LAKE', p_patch(jg)%n_patch_cells_g,            &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%fr_lake)

            CALL read_netcdf_data (ncid, 'DEPTH_LK', p_patch(jg)%n_patch_cells_g,           &
              &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
              &                     ext_data(jg)%atm%depth_lk)

          CASE ( iecham, ildf_echam )

            IF ( l_emiss ) THEN
              CALL read_netcdf_data (ncid, 'EMIS_RAD', p_patch(jg)%n_patch_cells_g,           &
                &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
                &                     ext_data(jg)%atm%emis_rad)
            ELSE
              ext_data(jg)%atm%emis_rad(:,:)= zemiss_def
            ENDIF

          END SELECT ! iforcing



          !
          ! derived external parameter fields
          !

          ! land sea mask at cell centers (LOGICAL)
          !
          i_nchdom  = MAX(1,p_patch(jg)%n_childdom)

          rl_start = 1
          rl_end   = min_rlcell_int

          i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
          i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)

          DO jb = i_startblk, i_endblk
            CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
              &                i_startidx, i_endidx, rl_start, rl_end)

            ! Loop starts with 1 instead of i_startidx because the start index is missing in RRTM
            DO jc = 1,i_endidx
              IF (ext_data(jg)%atm%fr_land(jc,jb) > 0.5_wp) THEN
                ext_data(jg)%atm%llsm_atm_c(jc,jb) = .TRUE.  ! land point
              ELSE
                ext_data(jg)%atm%llsm_atm_c(jc,jb) = .FALSE.  ! water point
              ENDIF
            ENDDO
          ENDDO


        ELSEIF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid

          CALL finish(TRIM(ROUTINE),&
            & 'Hexagonal grid is not supported, yet.')

        ENDIF

        !
        ! close file
        !
        IF( my_process_is_stdio()) CALL nf(nf_close(ncid))

      ENDDO

    ENDIF ! itopo


    !-------------------------------------------------------
    ! Read ozone
    !-------------------------------------------------------

    IF((irad_o3 == io3_clim) .OR. (irad_o3 == io3_ape)) THEN

      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test
      ELSE
        mpi_comm = p_comm_work
      ENDIF

      DO jg = 1,n_dom

        IF(my_process_is_stdio()) THEN
          ! open file
          !
          WRITE(ozone_file,'(a,I2.2,a)') 'o3_icon_DOM',jg,'.nc'
          CALL nf(nf_open(TRIM(ozone_file), NF_NOWRITE, ncid))
          WRITE(0,*)'read ozone levels'

          !          SELECT CASE (iequations)
          !          CASE(ihs_atm_temp,ihs_atm_theta)

          CALL nf(nf_inq_varid(ncid, TRIM(levelname), varid))
          CALL nf(nf_get_var_double(ncid, varid,zdummy_o3lev(:) ))

          !          CASE(inh_atmosphere)
          !            CALL nf(nf_inq_varid(ncid, TRIM(zlevelname), varid))
          !            CALL nf(nf_get_var_double(ncid, varid,zdummy_o3lev(:) ))
          !          END SELECT

          !
        ENDIF ! pe

        CALL p_bcast(zdummy_o3lev(:), p_io, mpi_comm)      

        !         SELECT CASE (iequations)
        !         CASE(ihs_atm_temp,ihs_atm_theta)

        DO jk=1,nlev_o3
          ext_data(jg)%atm_td%pfoz(jk)=zdummy_o3lev(jk)
        ENDDO

        ! define half levels of ozone pressure grid
        ! upper boundary: ph =      0.Pa -> extrapolation of uppermost value
        ! lower boundary: ph = 125000.Pa -> extrapolation of lowermost value
        ext_data(jg)%atm_td%phoz(1)           = 0._wp
        ext_data(jg)%atm_td%phoz(2:nlev_o3) = (ext_data(jg)%atm_td%pfoz(1:nlev_o3-1) &
          &                                   +  ext_data(jg)%atm_td%pfoz(2:nlev_o3))*.5_wp
        ext_data(jg)%atm_td%phoz(nlev_o3+1) = 125000._wp

        DO i=1,nlev_o3
          WRITE(0,*) 'full/half level press ozone ', i, ext_data(jg)%atm_td%pfoz(i),&
            &                                           ext_data(jg)%atm_td%phoz(i+1)
        ENDDO

        !         CASE(inh_atmosphere)
        !           DO jk=1,nlev_o3
        !             ext_data(jg)%atm_td%zf(jk)=zdummy_o3lev(jk)
        !           ENDDO
        !       END SELECT

        ! we have 2 different ozone files for hexagons and triangels at the moment
        ! therefore ozone data are both stored at the cell centers

        !         IF (p_patch(jg)%cell_type == 3) THEN     ! triangular grid

        CALL read_netcdf_data (ncid, TRIM(o3name), & ! &
          &                    p_patch(jg)%n_patch_cells_g,  &
          &                    p_patch(jg)%n_patch_cells,    &
          &                    p_patch(jg)%cells%glb_index,  & 
          &                    nlev_o3,  nmonths,          &
          &                    ext_data(jg)%atm_td%O3)

        !        ELSEIF (p_patch(jg)%cell_type == 6) THEN ! hexagonal grid
        !
        !          CALL read_netcdf_data (ncid, TRIM(o3name), & 
        !            &                    p_patch(jg)%n_patch_verts_g,  &
        !            &                    p_patch(jg)%n_patch_verts,    & 
        !            &                    p_patch(jg)%verts%glb_index,  &
        !            &                    nlev_o3, nmonths,           &
        !            &                    ext_data(jg)%atm_td%O3)
        !
        !         ENDIF ! patches


        WRITE(0,*)'MAX/min o3 ppmv',MAXVAL(ext_data(jg)%atm_td%O3(:,:,:,:)),&
          &                        MINVAL(ext_data(jg)%atm_td%O3(:,:,:,:))

        ! convert from ppmv to g/g only in case of APE ozone
        IF(irad_o3 == io3_ape) &
          &         ext_data(jg)%atm_td%O3(:,:,:,:)= ext_data(jg)%atm_td%O3(:,:,:,:)*ppmv2gg

        WRITE(0,*)'MAX/min o3 g/g',MAXVAL(ext_data(jg)%atm_td%O3(:,:,:,:)),&
          &                        MINVAL(ext_data(jg)%atm_td%O3(:,:,:,:))

        ! close file
        IF(my_process_is_stdio()) CALL nf(nf_close(ncid))

      ENDDO ! ndom
    END IF ! irad_o3


  END SUBROUTINE read_ext_data_atm
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Read ocean external data
  !!
  !! Read ocean external data from netcdf
  !!
  !! @par Revision History
  !! Initial revision by Stephan Lorenz, MPI (2011-06-17)
  !!
  SUBROUTINE read_ext_data_oce (p_patch, ext_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)

    CHARACTER(len=max_char_length), PARAMETER :: &
      routine = 'mo_ext_data:read_ext_data_oce'

    CHARACTER(filename_max) :: grid_file   !< file name for reading in
    CHARACTER(filename_max) :: omip_file   !< file name for reading in

    LOGICAL :: l_exist
    INTEGER :: jg, i_lev, i_cell_type, no_cells, no_verts, no_tst
    INTEGER :: ncid, dimid

    !REAL(wp):: z_flux(nproma, 12,p_patch(1)%nblks_c)
    REAL(wp):: z_flux(nproma,iforc_len,p_patch(1)%nblks_c)
    TYPE (t_keyword_list), POINTER :: keywords => NULL()

!-------------------------------------------------------------------------

    CALL message (TRIM(routine), 'start')

!-------------------------------------------------------------------------

    !  READ OCEAN BATHYMETRY

    !-------------------------------------------------------------------------

    !DO jg = 1,n_dom
    jg = 1

    i_lev       = p_patch(jg)%level
    i_cell_type = p_patch(jg)%cell_type

    IF(my_process_is_stdio()) THEN
      !
      ! bathymetry and lsm are read from the general ICON grid file
      ! bathymetry and lsm integrated into grid file by grid generator
      !WRITE (bathy_file,'(a,i0,a,i2.2,a)') 'iconR',nroot,'B',i_lev, '-grid.nc'
      !write(*,*) 'bathy_file = ',TRIM(bathy_file)
      !write(*,*) 'dynamics_grid_filename = ', TRIM(dynamics_grid_filename(1))
      CALL associate_keyword("<path>", TRIM(model_base_dir), keywords)
      grid_file = TRIM(with_keywords(keywords, TRIM(dynamics_grid_filename(jg))))

      INQUIRE (FILE=grid_file, EXIST=l_exist)
      IF (.NOT.l_exist) THEN
        CALL finish(TRIM(routine),'Grid file for reading bathymetry not found.')
      ENDIF

      !
      ! open file
      !
      CALL nf(nf_open(TRIM(grid_file), NF_NOWRITE, ncid))

      !
      ! get number of cells and vertices
      !
      CALL nf(nf_inq_dimid(ncid, 'cell', dimid))
      IF (i_cell_type == 3) THEN ! triangular grid
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))
      ELSEIF (i_cell_type == 6) THEN ! hexagonal grid
        CALL nf(nf_inq_dimlen(ncid, dimid, no_verts))
      ENDIF

      CALL nf(nf_inq_dimid(ncid, 'vertex', dimid))
      IF (i_cell_type == 3) THEN ! triangular grid
        CALL nf(nf_inq_dimlen(ncid, dimid, no_verts))
      ELSEIF (i_cell_type == 6) THEN ! hexagonal grid
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))
      ENDIF

      !
      ! check the number of cells and verts
      !
      IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
        CALL finish(TRIM(ROUTINE),&
        & 'Number of patch cells and cells in bathymetry file do not match.')
      ENDIF
      IF(p_patch(jg)%n_patch_verts_g /= no_verts) THEN
        CALL finish(TRIM(ROUTINE),&
        & 'Number of patch verts and verts in bathymetry file do not match.')
      ENDIF

      WRITE(message_text,'(3(a,i6))') 'No of cells =', no_cells, &
        &                           '  no of edges =', p_patch(jg)%n_patch_edges_g, &
        &                           '  no of verts =', no_verts
      CALL message( TRIM(routine),TRIM(message_text))
    ENDIF


    !-------------------------------------------------------
    !
    ! Read bathymetry for triangle centers and edges
    !
    !-------------------------------------------------------

    ! triangle center and edges

    IF (i_cell_type == 3) THEN     ! triangular grid

      ! These arrays are not included in standard icon-grid, but they are
      ! created by "create_ocean_grid"
      CALL read_netcdf_data (ncid, 'cell_elevation', p_patch(jg)%n_patch_cells_g,     &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     ext_data(jg)%oce%bathymetry_c)

      CALL read_netcdf_data (ncid, 'edge_elevation', p_patch(jg)%n_patch_edges_g,     &
        &                     p_patch(jg)%n_patch_edges, p_patch(jg)%edges%glb_index, &
        &                     ext_data(jg)%oce%bathymetry_e)

      ! get land-sea-mask on cells, integer marks are:
      ! inner sea (-2), boundary sea (-1, cells and vertices), boundary (0, edges),
      ! boundary land (1, cells and vertices), inner land (2)
      CALL read_netcdf_data (ncid, 'cell_sea_land_mask', p_patch(jg)%n_patch_cells_g, &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     ext_data(jg)%oce%lsm_ctr_c)

      CALL read_netcdf_data (ncid, 'edge_sea_land_mask', p_patch(jg)%n_patch_edges_g, &
        &                    p_patch(jg)%n_patch_edges, p_patch(jg)%edges%glb_index,  &
        &                    ext_data(jg)%oce%lsm_ctr_e)

    ENDIF

    !
    ! close file
    !
    IF(my_process_is_stdio()) CALL nf(nf_close(ncid))

    !ENDDO ! jg

    CALL message( TRIM(routine),'Ocean bathymetry for external data read' )

    !-------------------------------------------------------------------------

    !  READ OMIP FORCING

    !-------------------------------------------------------------------------

    IF (iforc_oce == 12) THEN

    !DO jg = 1,n_dom
      jg = 1

      i_lev       = p_patch(jg)%level
      i_cell_type = p_patch(jg)%cell_type

      IF(my_process_is_stdio()) THEN
        !
        WRITE (omip_file,'(a,i0,a,i2.2,a)') 'iconR',nroot,'B',i_lev, '-flux.nc'

        !omip_file=TRIM('/pool/data/ICON/external/iconR2B04-flux.nc')
        CALL message( TRIM(routine),'Ocean OMIP forcing flux file is: '//TRIM(omip_file) )
        INQUIRE (FILE=omip_file, EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(TRIM(routine),'OMIP forcing flux file is not found.')
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(omip_file), NF_NOWRITE, ncid))
        CALL message( TRIM(routine),'Ocean OMIP flux file opened for read' )

        !
        ! get and check number of cells in OMIP data
        !
        CALL nf(nf_inq_dimid(ncid, 'ncells', dimid))
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells))

        IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in OMIP flux file do not match.')
        ENDIF

        !
        ! get number of timesteps
        !
        CALL nf(nf_inq_dimid(ncid, 'time', dimid))
        CALL nf(nf_inq_dimlen(ncid, dimid, no_tst))
        !
        ! check
        !
        WRITE(message_text,'(A,I6,A)')  'Ocean OMIP flux file contains',no_tst,' data sets'
        CALL message( TRIM(routine), TRIM(message_text) )
        IF(no_tst /= iforc_len ) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of forcing timesteps is not equal iforc_len specified in namelist!')
        ENDIF
      ENDIF

      !-------------------------------------------------------
      !
      ! Read complete OMIP data for triangle centers
      !
      !-------------------------------------------------------

      ! provide OMIP fluxes for wind stress forcing
      ! 1:  'stress_x': zonal wind stress       [m/s]
      ! 2:  'stress_y': meridional wind stress  [m/s]
      ! 3:  'SST"     : sea surface temperature [K]

      ! zonal wind stress
      CALL read_netcdf_data (ncid, 'stress_x', p_patch(jg)%n_patch_cells_g,          &
        &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                    no_tst, z_flux(:,:,:))
      ext_data(jg)%oce%omip_forc_mon_c(:,:,:,1) = z_flux(:,:,:)

      ! meridional wind stress
      CALL read_netcdf_data (ncid, 'stress_y', p_patch(jg)%n_patch_cells_g,          &
        &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                    no_tst, z_flux)
      ext_data(jg)%oce%omip_forc_mon_c(:,:,:,2) = z_flux(:,:,:)

      ! SST
      CALL read_netcdf_data (ncid, 'SST', p_patch(jg)%n_patch_cells_g,           &
        &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                    no_tst, z_flux)
      ext_data(jg)%oce%omip_forc_mon_c(:,:,:,3) = z_flux(:,:,:)

      IF (iforc_omip == 2) THEN

      ! Read complete OMIP data sets for focing ocean model
      ! 4:  tafo(:,:),   &  ! 2 m air temperature                              [C]
      ! 5:  ftdew(:,:),  &  ! 2 m dew-point temperature                        [K]
      ! 6:  fu10(:,:) ,  &  ! 10 m wind speed                                  [m/s]
      ! 7:  fclou(:,:),  &  ! Fractional cloud cover
      ! 8:  pao(:,:),    &  ! Surface atmospheric pressure                     [hPa]
      ! 9:  fswr(:,:),   &  ! Incoming surface solar radiation                 [W/m]

        ! 2m-temperature
        CALL read_netcdf_data (ncid, 'temp_2m', p_patch(jg)%n_patch_cells_g,           &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,4) = z_flux(:,:,:)
     
        ! 2m dewpoint temperature
        CALL read_netcdf_data (ncid, 'dpt_temp_2m', p_patch(jg)%n_patch_cells_g,       &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,5) = z_flux(:,:,:)
     
        ! Scalar wind
        CALL read_netcdf_data (ncid, 'scalar_wind', p_patch(jg)%n_patch_cells_g,       &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,6) = z_flux(:,:,:)
     
        ! cloud cover
        CALL read_netcdf_data (ncid, 'cloud', p_patch(jg)%n_patch_cells_g,             &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,7) = z_flux(:,:,:)
     
        ! sea level pressure
        CALL read_netcdf_data (ncid, 'pressure', p_patch(jg)%n_patch_cells_g,          &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,8) = z_flux(:,:,:)
     
        ! total solar radiation
        CALL read_netcdf_data (ncid, 'tot_solar', p_patch(jg)%n_patch_cells_g,         &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,9) = z_flux(:,:,:)
     
        ! precipitation
        CALL read_netcdf_data (ncid, 'precip', p_patch(jg)%n_patch_cells_g,            &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,10) = z_flux(:,:,:)
     
        ! evaporation
        CALL read_netcdf_data (ncid, 'evap', p_patch(jg)%n_patch_cells_g,              &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,11) = z_flux(:,:,:)
     
        ! runoff
        CALL read_netcdf_data (ncid, 'runoff', p_patch(jg)%n_patch_cells_g,            &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,12) = z_flux(:,:,:)

      END IF

      ! provide heat and freshwater flux for focing ocean model
      IF (iforc_omip == 3) THEN

        ! net surface heat flux
        CALL read_netcdf_data (ncid, 'net_hflx', p_patch(jg)%n_patch_cells_g,          &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,4) = z_flux(:,:,:)
     
        ! surface freshwater flux
        CALL read_netcdf_data (ncid, 'net_fflx', p_patch(jg)%n_patch_cells_g,          &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,5) = z_flux(:,:,:)

      END IF

      ! provide 4 parts of heat and 2 parts of freshwater flux for focing ocean model
      IF (iforc_omip == 4) THEN

        ! surface short wave heat flux
        CALL read_netcdf_data (ncid, 'swflxsfc_avg', p_patch(jg)%n_patch_cells_g,      &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,4) = z_flux(:,:,:)
     
        ! surface long wave heat flux
        CALL read_netcdf_data (ncid, 'lwflxsfc_avg', p_patch(jg)%n_patch_cells_g,      &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,5) = z_flux(:,:,:)
     
        ! surface sensible heat flux
        CALL read_netcdf_data (ncid, 'shflx_avg',    p_patch(jg)%n_patch_cells_g,      &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,6) = z_flux(:,:,:)
     
        ! surface latent heat flux
        CALL read_netcdf_data (ncid, 'lhflx_avg',    p_patch(jg)%n_patch_cells_g,      &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,7) = z_flux(:,:,:)
     
        ! total precipiation
        CALL read_netcdf_data (ncid, 'precip', p_patch(jg)%n_patch_cells_g,            &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,8) = z_flux(:,:,:)
     
        ! evaporation
        CALL read_netcdf_data (ncid, 'evap'  , p_patch(jg)%n_patch_cells_g,            &
          &                    p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                    no_tst, z_flux)
        ext_data(jg)%oce%omip_forc_mon_c(:,:,:,9) = z_flux(:,:,:)

      END IF

      !
      ! close file
      !
      IF(my_process_is_stdio()) CALL nf(nf_close(ncid))

    !ENDDO

      CALL message( TRIM(routine),'Ocean OMIP fluxes for external data read' )

    END IF ! iforc_oce=12

  END SUBROUTINE read_ext_data_oce
  !-------------------------------------------------------------------------




  !-------------------------------------------------------------------------

   SUBROUTINE smooth_topography (p_patch, p_int, topography_c, topography_v, sso_stdh)

    TYPE(t_patch), INTENT(IN)       :: p_patch
    TYPE(t_int_state), INTENT(IN)   :: p_int
    REAL(wp), INTENT(INOUT)         :: topography_c(:,:), topography_v(:,:)
    REAL(wp), OPTIONAL, INTENT(INOUT) :: sso_stdh(:,:)

    ! local variables
    INTEGER  :: jg, jb, jc, iter, il, npts, niter
    INTEGER  :: i_startblk, nblks_c, i_startidx, i_endidx
    REAL(wp) :: z_topo(nproma,1,p_patch%nblks_c),z_nabla4_topo(nproma,1,p_patch%nblks_c),    &
      &         z_topo_old(nproma,1,p_patch%nblks_c),z_nabla2_topo(nproma,1,p_patch%nblks_c),&
      &         z_topo_c_sv(nproma,p_patch%nblks_c),z_hdiffmax(nproma,p_patch%nblks_c)
    REAL(wp) :: z_topo_v(nproma,1,p_patch%nblks_v)
    REAL(wp) :: zmaxtop,zmintop,z_topo_new,zdcoeff,z_heightdiff_threshold,zhdiff
    LOGICAL  :: lnabla2_mask(nproma,p_patch%nblks_c)


    jg = p_patch%id

    z_topo_v(:,1,:)  = topography_v(:,:)
    z_topo_c_sv(:,:) = topography_c(:,:)

    z_nabla4_topo(:,1,:) = 0._wp
    z_nabla2_topo(:,1,:) = 0._wp

    nblks_c    = p_patch%nblks_c

    zdcoeff = 0.05_wp ! diffusion coefficient for nabla2 diffusion
    niter   = 20      ! number of iterations for local nabla2 diffusion

    IF (n_iter_smooth_topo(jg) > 0) THEN
      WRITE(message_text,'(a,i3,a,i3)') 'number of topography smoothing steps in domain ', &
        jg, ': ', n_iter_smooth_topo(jg)
      CALL message('', TRIM(message_text))
    ENDIF

    ! Step 1: local nabla2 diffusion at grid points where the height difference to the
    ! neighbors exceeds a certain threshold value

    z_topo(:,1,:)   = topography_c(:,:)

    DO iter = 1, niter ! perform niter iterations
                       ! note: a variable number of iterations (with an exit condition) potentially 
                       ! causes trouble with MPI reproducibility

      z_heightdiff_threshold = heightdiff_threshold(jg)
      IF (iter >= niter-1 .AND. heightdiff_threshold(jg) < 2500._wp) &
        z_heightdiff_threshold = 0.75_wp*heightdiff_threshold(jg)

      z_hdiffmax(:,:)   = 0._wp
      lnabla2_mask(:,:) = .FALSE.

      CALL nabla2_scalar(z_topo, p_patch, p_int, z_nabla2_topo, opt_slev=1, opt_elev=1)

      npts = 0

      i_startblk = p_patch%cells%start_blk(2,1)

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 2)

        DO il=1,p_patch%cell_type
          DO jc = i_startidx, i_endidx
            z_hdiffmax(jc,jb) = MAX(z_hdiffmax(jc,jb), ABS(z_topo(jc,1,jb) - &
              &  z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1,              &
              &         p_patch%cells%neighbor_blk(jc,jb,il))) )

            IF (z_hdiffmax(jc,jb) > z_heightdiff_threshold) lnabla2_mask(jc,jb) = .TRUE.

          ENDDO
        ENDDO

      ENDDO

      ! set diffusion mask also true if one of the neighboring grid points 
      ! fulfills the height difference criterion
      i_startblk = p_patch%cells%start_blk(3,1)

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 3)

        DO il=1,p_patch%cell_type
          DO jc = i_startidx, i_endidx
            IF (z_hdiffmax(p_patch%cells%neighbor_idx(jc,jb,il),                   &
                p_patch%cells%neighbor_blk(jc,jb,il)) > z_heightdiff_threshold) THEN
              lnabla2_mask(jc,jb) = .TRUE.
            ENDIF
          ENDDO
        ENDDO

        DO jc = i_startidx, i_endidx
          IF (lnabla2_mask(jc,jb)) THEN
            npts = npts + 1
            topography_c(jc,jb) = z_topo(jc,1,jb) + zdcoeff *                      &
                                  p_patch%cells%area(jc,jb) * z_nabla2_topo(jc,1,jb)
          ENDIF
        ENDDO

      ENDDO

      z_topo(:,1,:) = topography_c(:,:)
      CALL sync_patch_array(SYNC_C, p_patch, z_topo)
      topography_c(:,:) = z_topo(:,1,:)

    ENDDO ! iteration of local nabla2 smoothing

    i_startblk = p_patch%cells%start_blk(3,1)

    ! Step 2: local nabla4 diffusion with monotonous limiter
    DO iter = 1, n_iter_smooth_topo(jg)
      z_topo(:,1,:)   = topography_c(:,:)
      z_topo_old(:,1,:) = z_topo(:,1,:)

      CALL nabla4_scalar(z_topo, p_patch, p_int, z_nabla4_topo, &
        & opt_slev=1, opt_elev=1, opt_nabla2=z_nabla2_topo )

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 3)

        DO jc = i_startidx, i_endidx

          !Limiter to avoid amplification of local extrema

          ! compute maximum (zmaxtop) and minimum (zmintop) of neighbor cells topography
          !set zmaxtop and zmintop to first neighbor's value
          zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,1),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,1))
          zmintop = zmaxtop
          DO il=2,p_patch%cell_type

            IF ( z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,il)) > &
              & zmaxtop ) THEN
              zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                &          p_patch%cells%neighbor_blk(jc,jb,il) )
            ENDIF
            IF ( z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,il)) < &
              & zmintop ) THEN
              zmintop = z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                &          p_patch%cells%neighbor_blk(jc,jb,il) )
            ENDIF
           !zmaxtop and zmintop are now max resp min of all neighbors
          ENDDO

          z_topo_new = z_topo (jc,1,jb) -  &
              fac_smooth_topo * z_nabla4_topo(jc,1,jb) * &
              p_patch%cells%area(jc,jb)*p_patch%cells%area(jc,jb)

          !If it was a local maximum in the old field, dont make it higher
          IF ( zmaxtop < z_topo (jc,1,jb) ) THEN
             IF ( z_nabla4_topo(jc,1,jb) < 0.0_wp ) CYCLE
          ENDIF
          !If it was a local minimum in the old field, dont make it lower:
          IF ( zmintop > z_topo (jc,1,jb) ) THEN
             IF ( z_nabla4_topo(jc,1,jb) > 0.0_wp ) CYCLE
          ENDIF

          !If it became a local maximum in the new field with regard to old neighbors, avoid it:
          IF (( zmaxtop < z_topo_new  ) .AND. ( z_nabla4_topo(jc,1,jb) < 0.0_wp )) THEN
            topography_c(jc,jb) = MAX(z_topo_old(jc,1,jb),zmaxtop)
            CYCLE
          ENDIF
          !If it became a local minimum in the new field with regard to old neighbors, avoid it:
          IF (( zmintop > z_topo_new  ) .AND. ( z_nabla4_topo(jc,1,jb) > 0.0_wp )) THEN
            topography_c(jc,jb) = MIN(z_topo_old(jc,1,jb),zmintop)
            CYCLE
          ENDIF

          topography_c(jc,jb)=z_topo_new

        ENDDO

      ENDDO

      z_topo(:,1,:)   = topography_c(:,:)
      CALL sync_patch_array(SYNC_C, p_patch, z_topo)
      topography_c(:,:)=z_topo(:,1,:)

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 3)
        DO jc = i_startidx, i_endidx

          !Limiter to avoid amplification of local extrema
          zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,1),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,1))
          zmintop = zmaxtop
          DO il=2,p_patch%cell_type

            IF ( z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,il)) > &
              & zmaxtop ) THEN
              zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                &          p_patch%cells%neighbor_blk(jc,jb,il) )
            ENDIF
            IF ( z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,il)) < &
              & zmintop ) THEN
              zmintop = z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                &          p_patch%cells%neighbor_blk(jc,jb,il) )
            ENDIF
           !zmaxtop and zmintop are now mx resp min of all neighbors
          ENDDO

           !If it became a local maximum in the new field, avoid it
           IF ( ( zmaxtop < z_topo(jc,1,jb) ) .AND. ( z_nabla4_topo(jc,1,jb) < 0.0_wp ) ) THEN
             topography_c(jc,jb) = MAX(z_topo_old(jc,1,jb),zmaxtop)
           !If it became a local minimum in the new field, avoid it
           ELSEIF (( zmintop > z_topo(jc,1,jb) ) .AND. ( z_nabla4_topo(jc,1,jb) > 0.0_wp )) THEN
             topography_c(jc,jb) = MIN(z_topo_old(jc,1,jb),zmintop)
           ENDIF

        ENDDO

      ENDDO

      z_topo(:,1,:)   = topography_c(:,:)
      CALL sync_patch_array(SYNC_C, p_patch, z_topo)
      topography_c(:,:)=z_topo(:,1,:)

    ENDDO !iter


    ! Interpolate smooth topography from cells to vertices
    z_topo(:,1,:)   = topography_c(:,:)
    CALL cells2verts_scalar(z_topo,p_patch,p_int%cells_aw_verts,z_topo_v,1,1)

    CALL sync_patch_array(SYNC_V,p_patch,z_topo_v)

    topography_v(:,:) = z_topo_v(:,1,:)

    ! Apply correction to SSO standard deviation
    IF (PRESENT(sso_stdh)) THEN

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 3)

        DO jc = i_startidx, i_endidx

          zhdiff = topography_c(jc,jb) - z_topo_c_sv(jc,jb)
          sso_stdh(jc,jb) = SQRT(sso_stdh(jc,jb)**2 + zhdiff**2)

        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE smooth_topography

END MODULE mo_ext_data

