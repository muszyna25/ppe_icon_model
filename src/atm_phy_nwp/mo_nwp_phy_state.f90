#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif
MODULE mo_nwp_phy_state
!>
!!  !MODULE:  mo_nwp_phy_state\\
!!
!! Description:  Contains the data structures
!!  to store the physical model state and other auxiliary variables
!!  in order to run the ECHAM physics.
!!  Constructors and destructors for these data structures as well as
!!  initialization of fields are also defined here.
!!  This module should be an analogon to 'mo_hydro_state.f90'

!!  TODO/To think about:
!     - should physics be called before or after dynamics?
!     - allocate fluxes at edges instead at the centers?
!     - horizontal/vertical tracer flux (reconstruct q'v_n' into q'u' and q'v') ?
!     - provide the "virt_inc" with meaning
!     - where to provide the lat/lon info for radiation?
!     - how to implement the echam-modules - rewriting them or "capsulate"?
!     - revision of fields if there are needed or tp be replaced
!     - fill the physics tendency construction/destruction subroutine
!     - later implement already calculated icon gradients for echam physics
!     - think about variables for flexible time steps
!!
!! @author Kristina Froehlich, DWD
!! @author Marco Giorgetta, MPI-M
!!
!!
!! @par Revision History
!! Initial  by Kristina Froehlich (2009-06-10)
!! Memory allocation method changed from explicit allocation to Luis' 
!! infrastructure by Kristina Froehlich (MPI-M, 2011-04-27)
!! Added clch, clcm, clcl, hbas_con, htop_con by Helmut Frank (DWD, 2013-01-17)
!! Added hzerocl, gust10                      by Helmut Frank (DWD, 2013-03-13)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

! !USES:

USE mo_kind,                ONLY: wp, i8
USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag, t_nwp_phy_tend
USE mo_impl_constants,      ONLY: success, max_char_length,           &
  &                               VINTP_METHOD_LIN,VINTP_METHOD_QV,   &
  &                               TASK_COMPUTE_RH, TASK_COMPUTE_PV,   &
  &                               iedmf,                              &
  &                               HINTP_TYPE_LONLAT_NNB,              &
  &                               HINTP_TYPE_LONLAT_BCTR,             &
  &                               HINTP_TYPE_LONLAT_RBF,              &
  &                               nexlevs_rrg_vnest, RTTOV_BT_CL,     &
  &                               RTTOV_RAD_CL, RTTOV_RAD_CS,         &
  &                               TLEV_NNOW_RCF, iss, iorg, ibc, iso4,&
  &                               idu, nclass_aero
USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL,             &
  &                               GRID_CELL
USE mo_parallel_config,     ONLY: nproma
USE mo_run_config,          ONLY: nqtendphy, iqv, iqc, iqi, lart, ldass_lhn
USE mo_exception,           ONLY: message, finish !,message_text
USE mo_model_domain,        ONLY: t_patch, p_patch, p_patch_local_parent
USE mo_grid_config,         ONLY: n_dom, n_dom_start
USE mo_linked_list,         ONLY: t_var_list
USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config, icpl_aero_conv, iprog_aero
USE turb_data,              ONLY: ltkecon
USE mo_extpar_config,       ONLY: itype_vegetation_cycle
USE mo_radiation_config,    ONLY: irad_aero
USE mo_lnd_nwp_config,      ONLY: ntiles_total, ntiles_water, nlev_soil
USE mo_var_list,            ONLY: default_var_list_settings, &
  &                               add_var, add_ref, new_var_list, delete_var_list
USE mo_var_groups,          ONLY: groups, MAX_GROUPS
USE mo_var_metadata_types,  ONLY: POST_OP_SCALE, CLASS_SYNSAT, CLASS_CHEM, VARNAME_LEN
USE mo_var_metadata,        ONLY: create_vert_interp_metadata,  &
  &                               create_hor_interp_metadata,   &
  &                               vintp_types, post_op, &
  &                               new_action, actions
USE mo_nwp_parameters,      ONLY: t_phy_params
USE mo_cf_convention,       ONLY: t_cf_var
USE mo_grib2,               ONLY: t_grib2_var, grib2_var, t_grib2_int_key, OPERATOR(+)
USE mo_io_config,           ONLY: lflux_avg
USE mo_cdi,                 ONLY: TSTEP_MIN, TSTEP_MAX, TSTEP_INSTANT, TSTEP_CONSTANT, &
    &                             TSTEP_AVG, TSTEP_ACCUM, DATATYPE_PACK16,             &
    &                             DATATYPE_FLT32, DATATYPE_FLT64, GRID_UNSTRUCTURED
USE mo_zaxis_type,          ONLY: ZA_REFERENCE, ZA_REFERENCE_HALF,          &
  &                               ZA_SURFACE, ZA_HEIGHT_2M, ZA_HEIGHT_10M,       &
  &                               ZA_HEIGHT_2M_LAYER, ZA_TOA, ZA_DEPTH_BELOW_LAND,   &
  &                               ZA_PRESSURE_0, ZA_PRESSURE_400,&
  &                               ZA_PRESSURE_800, ZA_CLOUD_BASE, ZA_CLOUD_TOP,  &
  &                               ZA_ISOTHERM_ZERO
USE mo_physical_constants,  ONLY: grav
USE mo_ls_forcing_nml,      ONLY: is_ls_forcing

USE mo_advection_config,     ONLY: advection_config
USE mo_synsat_config,        ONLY: lsynsat, num_images, get_synsat_name, num_sensors, &
  &                                total_numchans, get_synsat_grib_triple
USE mo_art_config,           ONLY: nart_tendphy
USE mo_art_tracer_interface, ONLY: art_tracer_interface
USE mo_action,               ONLY: ACTION_RESET
USE mo_les_nml,              ONLY: turb_profile_list, turb_tseries_list
USE mo_io_config,            ONLY: lnetcdf_flt64_output, gust_interval
USE mtime,                   ONLY: max_timedelta_str_len, getPTStringFromMS
USE mo_name_list_output_config, ONLY: &
  &                                first_output_name_list, is_variable_in_output

IMPLICIT NONE
PRIVATE

INCLUDE 'netcdf.inc'


!public interface
!
! subroutines
PUBLIC :: construct_nwp_phy_state
PUBLIC :: destruct_nwp_phy_state
!variables
PUBLIC :: prm_diag 
PUBLIC :: prm_nwp_tend
PUBLIC :: phy_params
PUBLIC :: prm_nwp_diag_list  !< variable lists
PUBLIC :: prm_nwp_tend_list  !< variable lists
!

!!--------------------------------------------------------------------------
!!                          STATE VARIABLES 
!!--------------------------------------------------------------------------
  TYPE(t_nwp_phy_diag), ALLOCATABLE, TARGET :: prm_diag(:) !< shape: (n_dom)
  TYPE(t_nwp_phy_tend), ALLOCATABLE :: prm_nwp_tend(:)     !< shape: (n_dom)
!-------------------------------------------------------------------------
  
!!--------------------------------------------------------------------------
!!                          VARIABLE LISTS
!!--------------------------------------------------------------------------
  TYPE(t_var_list),ALLOCATABLE :: prm_nwp_diag_list(:)  !< shape: (n_dom)
  TYPE(t_var_list),ALLOCATABLE :: prm_nwp_tend_list(:)  !< shape: (n_dom)

!!-------------------------------------------------------------------------
!! Parameters of various physics parameterizations that have to be 
!! domain-dependent (computed during physics initialization phase)
!!-------------------------------------------------------------------------
  TYPE (t_phy_params), ALLOCATABLE :: phy_params(:)  !< shape: (n_dom)


CONTAINS

!-------------------------------------------------------------------------

SUBROUTINE construct_nwp_phy_state( p_patch, l_rh, l_pv)

TYPE(t_patch), TARGET, INTENT(in) :: p_patch(n_dom)
LOGICAL, INTENT(IN) :: l_rh(n_dom), &!< Flag. TRUE if computation of relative humidity desired
                       l_pv(n_dom)   !< Flag. TRUE if computation of potential vorticity desired

CHARACTER(len=max_char_length) :: listname
INTEGER ::  jg,ist, nblks_c, nlev, nlevp1

!-------------------------------------------------------------------------

CALL message('mo_nwp_phy_state:construct_nwp_state', &
  'start to construct 3D state vector')


  ! Allocate pointer arrays prm_diag_nwp and prm_nwp_tend, 
  ! as well as the corresponding list arrays.

  ALLOCATE(prm_diag(n_dom), prm_nwp_diag_list(n_dom),STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_state', &
      'allocation of diagnostic physical array and list failed')
  ENDIF

  ALLOCATE(prm_nwp_tend(n_dom), prm_nwp_tend_list(n_dom), STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_state', &
      'allocation of tendency physical array and list failed')
  ENDIF

  DO jg = 1, n_dom

     !determine size of arrays
     nblks_c = p_patch(jg)%nblks_c
     
     ! number of vertical levels
     nlev   = p_patch(jg)%nlev
     nlevp1 = p_patch(jg)%nlevp1
     
     WRITE(listname,'(a,i2.2)') 'prm_diag_of_domain_',jg

     CALL new_nwp_phy_diag_list( jg, nlev, nlevp1, nblks_c, TRIM(listname),   &
       &                         prm_nwp_diag_list(jg), prm_diag(jg), l_rh(jg), l_pv(jg))
     !
     WRITE(listname,'(a,i2.2)') 'prm_tend_of_domain_',jg
     CALL new_nwp_phy_tend_list ( jg, nlev, nblks_c,&
                                & TRIM(listname), prm_nwp_tend_list(jg), prm_nwp_tend(jg))
  ENDDO


  ! Allocate variable of type t_phy_params containing domain-dependent parameters
  !
  ALLOCATE(phy_params(n_dom), STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_state', &
      'allocation of phy_params array failed')
  ENDIF

  
  CALL message('mo_nwp_phy_state:construct_nwp_state', &
    'construction of state vector finished')

END SUBROUTINE construct_nwp_phy_state

!
SUBROUTINE destruct_nwp_phy_state

  INTEGER :: jg, ist  !< grid level/domain index

  CALL message('mo_nwp_phy_state:destruct_nwp_phy_state', &
  'start to destruct 3D state vector')

  DO jg = 1,n_dom
    CALL delete_var_list( prm_nwp_diag_list(jg) )
    CALL delete_var_list( prm_nwp_tend_list (jg) )

    IF (ASSOCIATED(prm_diag(jg)%buffer_rttov))  DEALLOCATE(prm_diag(jg)%buffer_rttov)  
    IF (ALLOCATED(prm_diag(jg)%synsat_image))   DEALLOCATE(prm_diag(jg)%synsat_image)
  END DO

  DEALLOCATE(prm_diag, prm_nwp_diag_list, STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:destruct_nwp_phy_state', &
       &  'deallocation of NWP physics diagnostic array and list failed')
  ENDIF
 
  DEALLOCATE(prm_nwp_tend, prm_nwp_tend_list, STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:destruct_nwp_phy_state', &
         &' deallocation of NWP physics tendencies array and list failed') 
  ENDIF

  DEALLOCATE(phy_params, STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:destruct_nwp_phy_state', &
         &' deallocation of phy_params array failed') 
  ENDIF

  CALL message('mo_nwp_phy_state:destruct_nwp_phy_state', &
    'destruction of 3D state vector finished')

END SUBROUTINE destruct_nwp_phy_state

     !
SUBROUTINE new_nwp_phy_diag_list( k_jg, klev, klevp1, kblks, &
                     & listname, diag_list, diag, l_rh, l_pv)

    INTEGER,INTENT(IN) :: klev, klevp1, kblks, k_jg !< dimension sizes

    CHARACTER(len=*),INTENT(IN)     :: listname
    CHARACTER(len=max_char_length)  :: vname_prefix
    CHARACTER(LEN=1)                :: csfc
    CHARACTER(LEN=2)                :: caer

    TYPE(t_var_list)    ,INTENT(INOUT) :: diag_list
    TYPE(t_nwp_phy_diag),INTENT(INOUT) :: diag
    LOGICAL, INTENT(IN) :: l_rh, & !< Flag. TRUE if computation of relative humidity desired
                           l_pv    !< Flag. TRUE if computation of potential vorticity desired
    ! Local variables

    INTEGER :: n_updown = 7 !> number of up/downdrafts variables

    TYPE(t_cf_var)    ::    cf_desc, new_cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d(3), shape3dsubs(3), &
      &        shape3dsubsw(3), shape3d_synsat(3),     &
      &        shape2d_synsat(2), shape3d_aero(3)
    INTEGER :: shape3dkp1(3)
    INTEGER :: ibits,  kcloud
    INTEGER :: jsfc, ist
    CHARACTER(len=NF_MAX_NAME) :: long_name
    CHARACTER(len=21) :: name
    CHARACTER(len=3)  :: prefix
    CHARACTER(len=8)  :: meaning
    CHARACTER(len=10) :: varunits  ! variable units, depending on "lflux_avg"
    INTEGER :: a_steptype
    LOGICAL :: lrestart

    LOGICAL :: lradiance, lcloudy
    INTEGER :: ichan, idiscipline, icategory, inumber, &
      &        wave_no, wave_no_scalfac, iimage, isens, k
    CHARACTER(LEN=VARNAME_LEN) :: shortname
    CHARACTER(LEN=128)         :: longname, unit
    CHARACTER(len=max_timedelta_str_len) :: gust_int
    !
    INTEGER :: constituentType                 ! for variable of class 'chem'

    INTEGER :: datatype_flt
    LOGICAL :: in_group(MAX_GROUPS)            ! for adding a variable to one or more groups 

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ibits = DATATYPE_PACK16 ! bits "entropy" of horizontal slice

    shape2d        = (/nproma,               kblks            /)
    shape3d        = (/nproma, klev,         kblks            /)
    shape3dkp1     = (/nproma, klevp1,       kblks            /)
    shape3dsubs    = (/nproma, kblks,        ntiles_total     /)
    shape3dsubsw   = (/nproma, kblks,        ntiles_total+ntiles_water /)
    shape3d_aero   = (/nproma, nclass_aero, kblks             /)

    ! Register a field list and apply default settings

    CALL new_var_list( diag_list, TRIM(listname), patch_id=k_jg )
    CALL default_var_list_settings( diag_list,                 &
                                  & lrestart=.TRUE.  )

   
    !------------------------------
    ! Meteorological quantities
    !------------------------------

    !-------------------
    ! Clouds and precip
    !------------------
    ! 2D and 3D variables

    kcloud= 3


    ! &      diag%rain_gsp_rate(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_gsp_rate', 'kg m-2 s-1', 'gridscale rain rate ', &
      &                   datatype_flt)
    grib2_desc = grib2_var(0, 1, 77, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rain_gsp_rate', diag%rain_gsp_rate,            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & isteptype=TSTEP_INSTANT )


    ! &      diag%snow_gsp_rate(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_gsp_rate', 'kg m-2 s-1', 'gridscale snow rate', &
      &                   datatype_flt)
    grib2_desc = grib2_var(0, 1, 56, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'snow_gsp_rate', diag%snow_gsp_rate,            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & isteptype=TSTEP_INSTANT )

    ! For graupel scheme 
    IF (atm_phy_nwp_config(k_jg)%lhave_graupel) THEN
      
      ! &      diag%graupel_gsp_rate(nproma,nblks_c)
      cf_desc    = t_cf_var('graupel_gsp_rate', 'kg m-2 s-1', 'gridscale graupel rate', &
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 1, 75, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'graupel_gsp_rate', diag%graupel_gsp_rate,            &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d,                                             &
                  & isteptype=TSTEP_INSTANT )
    ENDIF

    !For two moment microphysics
    SELECT CASE (atm_phy_nwp_config(k_jg)%inwp_gscp)
    CASE (4,5,6)

       ! &      diag%ice_gsp_rate(nproma,nblks_c)
      cf_desc    = t_cf_var('ice_gsp_rate', 'kg m-2 s-1', 'gridscale ice rate', &
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 1, 68, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'ice_gsp_rate', diag%ice_gsp_rate,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d, isteptype=TSTEP_INSTANT )
      
       ! &      diag%hail_gsp_rate(nproma,nblks_c)
      cf_desc    = t_cf_var('hail_gsp_rate', 'kg m-2 s-1', 'gridscale hail rate', &
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 1, 73, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'hail_gsp_rate', diag%hail_gsp_rate,            &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d, isteptype=TSTEP_INSTANT )

    END SELECT


    ! &      diag%rain_con_rate(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_con_rate', 'kg m-2 s-1', 'convective rain rate', &
      &                   datatype_flt)
    grib2_desc = grib2_var(0, 1, 76, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rain_con_rate', diag%rain_con_rate,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & isteptype=TSTEP_INSTANT )


    ! &      diag%snow_con_rate(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_con_rate', 'kg m-2 s-1', 'convective snow rate', &
      &                   datatype_flt)
    grib2_desc = grib2_var(0, 1, 55, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'snow_con_rate', diag%snow_con_rate,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & isteptype=TSTEP_INSTANT )

    ! &      diag%rain_con_rate_3d(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('rain_con_rate_3d', 'kg m-2 s-1',                &
      &          '3d convective rain rate', datatype_flt)
    grib2_desc = grib2_var(0, 1, 76, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rain_con_rate_3d', diag%rain_con_rate_3d,       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                & ldims=shape3dkp1,                                           &
                & isteptype=TSTEP_INSTANT )


    ! &      diag%snow_con_rate_3d(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('snow_con_rate_3d', 'kg m-2 s-1',                   &
      &          '3d convective snow rate', datatype_flt)
    grib2_desc = grib2_var(0, 1, 55, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'snow_con_rate_3d', diag%snow_con_rate_3d,       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                & ldims=shape3dkp1,                                           &
                & isteptype=TSTEP_INSTANT )


    IF ( atm_phy_nwp_config(k_jg)%inwp_turb == iedmf ) THEN

      ! &      diag%rain_edmf_rate_3d(nproma,nlevp1,nblks_c)
      cf_desc    = t_cf_var('rain_edmf_rate_3d', 'kg m-2 s-1',                &
        &          '3d EDMF convective rain rate', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'rain_edmf_rate_3d', diag%rain_edmf_rate_3d,       &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                  & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                  & ldims=shape3dkp1,                                           &
                  & isteptype=TSTEP_INSTANT )
      
      
      ! &      diag%snow_edmf_rate_3d(nproma,nlevp1,nblks_c)
      cf_desc    = t_cf_var('snow_edmf_rate_3d', 'kg m-2 s-1',                   &
        &          '3d EDMF convective snow rate', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'snow_edmf_rate_3d', diag%snow_edmf_rate_3d,       &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                  & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                  & ldims=shape3dkp1,                                           &
                  & isteptype=TSTEP_INSTANT )

    ENDIF


    ! &      diag%rain_gsp(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_gsp ', 'kg m-2 ', 'gridscale rain ', datatype_flt)
    grib2_desc = grib2_var(0, 1, 77, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rain_gsp', diag%rain_gsp,                      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d, in_group=groups("precip_vars"),             &
                & isteptype=TSTEP_ACCUM ,                                    &
                & hor_interp=create_hor_interp_metadata(                     &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                   &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB                     &
                & ) )

    ! &      diag%rain_gsp0(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_gsp0', 'kg m-2 ', 'gridscale rain0', datatype_flt)
    grib2_desc = grib2_var(0, 1, 77, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rain_gsp0', diag%rain_gsp0,                    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & loutput=.false., lrestart=.TRUE.,                          &
                & isteptype=TSTEP_ACCUM )


    ! &      diag%snow_gsp(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_gsp', 'kg m-2 ', 'gridscale snow', datatype_flt)
    grib2_desc = grib2_var(0, 1, 56, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'snow_gsp', diag%snow_gsp,                      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & in_group=groups("precip_vars"),                            &
                & isteptype=TSTEP_ACCUM ,                                    &
                & hor_interp=create_hor_interp_metadata(                     &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                   &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB                     &
                & ) )

    ! &      diag%snow_gsp0(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_gsp0', 'kg m-2 ', 'gridscale snow0', datatype_flt)
    grib2_desc = grib2_var(0, 1, 56, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'snow_gsp0', diag%snow_gsp0,                    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & loutput=.false., lrestart=.TRUE.,                          &
                & isteptype=TSTEP_ACCUM )


    !Surface precipitation variables for graupel scheme and two moment microphysics
    IF (atm_phy_nwp_config(k_jg)%lhave_graupel) THEN
       ! &      diag%graupel_gsp(nproma,nblks_c)
      cf_desc    = t_cf_var('graupel_gsp', 'kg m-2', 'gridscale graupel',      &
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 1, 75, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'graupel_gsp', diag%graupel_gsp,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d, in_group=groups("precip_vars"),             &
                  & isteptype=TSTEP_ACCUM,                                     &
                  & hor_interp=create_hor_interp_metadata(                     &
                  &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                   &
                  &    fallback_type=HINTP_TYPE_LONLAT_NNB                     &
                  & ) )
    ENDIF

    SELECT CASE (atm_phy_nwp_config(k_jg)%inwp_gscp)
    CASE (4,5,6)

       ! &      diag%ice_gsp(nproma,nblks_c)
      cf_desc    = t_cf_var('ice_gsp', 'kg m-2', 'gridscale ice', datatype_flt)
      grib2_desc = grib2_var(0, 1, 68, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'ice_gsp', diag%ice_gsp,                        &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d, in_group=groups("precip_vars"),             &
                  & isteptype=TSTEP_ACCUM,                                     &
                  & hor_interp=create_hor_interp_metadata(                     &
                  &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                   &
                  &    fallback_type=HINTP_TYPE_LONLAT_NNB                     &
                  & ) )
      

       ! &      diag%hail_gsp(nproma,nblks_c)
      cf_desc    = t_cf_var('hail_gsp', 'kg m-2', 'gridscale hail', datatype_flt)
      grib2_desc = grib2_var(0, 1, 73, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'hail_gsp', diag%hail_gsp,                      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d, in_group=groups("precip_vars"),             &
                  & isteptype=TSTEP_ACCUM,                                     &
                  & hor_interp=create_hor_interp_metadata(                     &
                  &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                   &
                  &    fallback_type=HINTP_TYPE_LONLAT_NNB                     &
                  & ) )

    END SELECT

    ! &      diag%prec_gsp(nproma,nblks_c)
    cf_desc    = t_cf_var('prec_gsp', 'kg m-2', 'gridscale precip', datatype_flt)
    grib2_desc = grib2_var(0, 1, 54, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'prec_gsp', diag%prec_gsp,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & in_group=groups("precip_vars"),                             &
                & isteptype=TSTEP_ACCUM ,                                     &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB                      &
                & ) )



    ! &      diag%rain_con(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_con', 'kg m-2 ', 'convective rain', datatype_flt)
    grib2_desc = grib2_var(0, 1, 76, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rain_con', diag%rain_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, in_group=groups("precip_vars"),              &
                & isteptype=TSTEP_ACCUM ,                                     &
                & hor_interp=create_hor_interp_metadata(                     &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                   &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB                     &
                & ) )

    ! &      diag%rain_con0(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_con0', 'kg m-2 ', 'convective rain0', datatype_flt)
    grib2_desc = grib2_var(0, 1, 76, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rain_con0', diag%rain_con0,                    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & loutput=.false., lrestart=.TRUE.,                          &
                & isteptype=TSTEP_ACCUM )


    ! &      diag%snow_con(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_con', 'kg m-2', 'convective snow', datatype_flt)
    grib2_desc = grib2_var(0, 1, 55, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'snow_con', diag%snow_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & in_group=groups("precip_vars"),                             &
                & isteptype=TSTEP_ACCUM ,                                     &
                & hor_interp=create_hor_interp_metadata(                     &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                   &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB                     &
                & ) )

    ! &      diag%snow_con0(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_con0', 'kg m-2', 'convective snow0', datatype_flt)
    grib2_desc = grib2_var(0, 1, 55, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'snow_con0', diag%snow_con0,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & loutput=.false., lrestart=.TRUE.,                           &
                & isteptype=TSTEP_ACCUM )



    ! &      diag%prec_con(nproma,nblks_c)
    cf_desc    = t_cf_var('prec_con', 'kg m-2', 'convective precip', datatype_flt)
    grib2_desc = grib2_var(0, 1, 37, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'prec_con', diag%prec_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & in_group=groups("precip_vars"),                             &
                & isteptype=TSTEP_ACCUM ,                                     &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB                      &
                & ) )


    ! &      diag%tot_prec(nproma,nblks_c)
    cf_desc    = t_cf_var('tot_prec', 'kg m-2', 'total precip', datatype_flt)
    grib2_desc = grib2_var(0, 1, 52, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tot_prec', diag%tot_prec,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & in_group=groups("precip_vars"),                             &
                & isteptype=TSTEP_ACCUM ,                                     &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB                      &
                & ) )



    ! &      diag%prec_con_rate_avg(nproma,nblks_c)
    cf_desc    = t_cf_var('prec_con_rate_avg', 'kg m-2 s-1',                  &
      &          'convective precip rate, time average', datatype_flt)
    grib2_desc = grib2_var(0, 1, 37, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'prec_con_rate_avg', diag%prec_con_rate_avg,     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,  lrestart=.FALSE.,                           &
                & in_group=groups("additional_precip_vars"),                  &
                & isteptype=TSTEP_AVG,                                        &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB                      &
                & ) )

    ! &      diag%prec_gsp_rate_avg(nproma,nblks_c)
    cf_desc    = t_cf_var('prec_gsp_rate_avg', 'kg m-2 s-1',                  &
      &          'gridscale precip rate, time average', datatype_flt)
    grib2_desc = grib2_var(0, 1, 54, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'prec_gsp_rate_avg', diag%prec_gsp_rate_avg,     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE.,                            &
                & in_group=groups("additional_precip_vars"),                  &
                & isteptype=TSTEP_AVG,                                        &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB                      &
                & ) )

    ! &      diag%tot_prec_rate_avg(nproma,nblks_c)
    cf_desc    = t_cf_var('tot_prec_rate_avg', 'kg m-2 s-1',                  &
      &          'total precip rate, time average', datatype_flt)
    grib2_desc = grib2_var(0, 1, 52, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tot_prec_rate_avg', diag%tot_prec_rate_avg,     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE.,                            &
                & isteptype=TSTEP_AVG,                                        &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB                      &
                & ) )



    ! &      diag%cape(nproma,nblks_c)
    cf_desc    = t_cf_var('cape', 'J kg-1 ', 'conv avail pot energy', datatype_flt)
    grib2_desc = grib2_var(0, 7, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'cape', diag%cape,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE.,                            &
                & in_group=groups("additional_precip_vars"),                  &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB                      &
                & ) )

    ! &      diag%cape_ml(nproma,nblks_c)
    ! typeOfLevel ZA_SURFACE is changed to 192 in vlistDefVarIntKey
    cf_desc    = t_cf_var('cape_ml', 'J kg-1 ', 'cape of mean surface layer parcel', datatype_flt)
    grib2_desc = grib2_var(0, 7, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)         &
    &           + t_grib2_int_key("typeOfFirstFixedSurface", 192)
    CALL add_var( diag_list, 'cape_ml', diag%cape_ml,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE.,                            &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB) )

    ! &      diag%cin_ml(nproma,nblks_c)
    ! typeOfLevel ZA_SURFACE is changed to 192 in vlistDefVarIntKey
    cf_desc    = t_cf_var('cin_ml', 'J kg-1 ', 'convective inhibition of mean surface layer parcel', datatype_flt)
    grib2_desc = grib2_var(0, 7, 7, ibits, GRID_UNSTRUCTURED, GRID_CELL)         &
    &           + t_grib2_int_key("typeOfFirstFixedSurface", 192)
    CALL add_var( diag_list, 'cin_ml', diag%cin_ml,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE.,                            &
!!$                & lmiss=.TRUE., missval=-999.9_wp,                            &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB) )

    ! &      diag%gust10(nproma,nblks_c)
    CALL getPTStringFromMS(NINT(1000*gust_interval(k_jg), i8), gust_int)
    cf_desc    = t_cf_var('gust10', 'm s-1 ', 'gust at 10 m', datatype_flt)
    grib2_desc = grib2_var( 0, 2, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'gust10', diag%gust10,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,  &
                & ldims=shape2d, lrestart=.TRUE., in_group=groups("pbl_vars"), &
                & isteptype=TSTEP_MAX,                                         &
                & initval=0._wp, resetval=0._wp,                               &
                & action_list=actions(new_action(ACTION_RESET,gust_int)) )

    ! &      diag%dyn_gust(nproma,nblks_c)
    cf_desc    = t_cf_var('dyn_gust', 'm s-1 ', 'dynamical gust', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'dyn_gust', diag%dyn_gust,                        &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,  &
                & ldims=shape2d, lrestart=.FALSE., isteptype=TSTEP_INSTANT,    &
                & hor_interp=create_hor_interp_metadata(                       &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,             &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),             &
                & loutput=.TRUE.                                               )

    ! &      diag%con_gust(nproma,nblks_c)
    cf_desc    = t_cf_var('con_gust', 'm s-1 ', 'convective contribution to wind gust', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'con_gust', diag%con_gust,                        &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,  &
                & ldims=shape2d, lrestart=.TRUE., isteptype=TSTEP_INSTANT,     &
                & loutput=.TRUE. )
   
    ! &      diag%rain_upd(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_upd', 'kg m-2 s-1', 'rain in updroughts', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rain_upd', diag%rain_upd,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE. )

    ! &      diag%con_udd(nproma,nlev,nblks,8)
    cf_desc    = t_cf_var('con_udd', 'unit ', 'convective up/downdraft fields', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'con_udd', diag%con_udd,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                & ldims=(/nproma,klev,kblks,n_updown/),                       &
                & lrestart=.FALSE., loutput=.FALSE. )

    ! &      diag%mbas_con(nproma,nblks_c)
    cf_desc    = t_cf_var('mbas_con', '', 'cloud base level index', datatype_flt)
    grib2_desc = grib2_var(0, 6, 194, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'mbas_con', diag%mbas_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    ! &      diag%mtop_con(nproma,nblks_c)
    cf_desc    = t_cf_var('mtop_con', '', 'cloud top level index', datatype_flt)
    grib2_desc = grib2_var(0, 6, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'mtop_con', diag%mtop_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    ! &      diag%locum(nproma,nblks_c)
    cf_desc    = t_cf_var('locum', '', 'convective activity indicator', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'locum', diag%locum,                             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    ! &      diag%ldshcv(nproma,nblks_c)
    cf_desc    = t_cf_var('ldshcv', '', 'shallow convection indicator', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'ldshcv', diag%ldshcv,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    ! &      diag%ktype(nproma,nblks_c)
    cf_desc    = t_cf_var('ktype', '', 'type of convection', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'ktype', diag%ktype,                              &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,  &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ) )

    ! &      diag%k850(nproma,nblks_c)
    cf_desc    = t_cf_var('k850', '', 'level index corresponding to the HAG of the 850hPa level', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'k850', diag%k850,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,  &
                & isteptype=TSTEP_CONSTANT,                                    &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ) )

    ! &      diag%k950(nproma,nblks_c)
    cf_desc    = t_cf_var('k950', '', 'level index corresponding to the HAG of the 950hPa level', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'k950', diag%k950,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,  &
                & isteptype=TSTEP_CONSTANT,                                    &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ) )

    ! &      diag%k800(nproma,nblks_c)
    cf_desc    = t_cf_var('k800', '', 'level index corresponding to the HAG of the 800hPa level', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'k800', diag%k800,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,  &
                & isteptype=TSTEP_CONSTANT,                                    &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ) )

    ! &      diag%k400(nproma,nblks_c)
    cf_desc    = t_cf_var('k400', '', 'level index corresponding to the HAG of the 400hPa level', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'k400', diag%k400,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,  &
                & isteptype=TSTEP_CONSTANT,                                    &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ) )


    ! &      diag%ktop_envel(nproma,nblks_c)
    cf_desc    = t_cf_var('ktop_envel', '', 'level index of upper boundary of SSO envelope layer', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'ktop_envel', diag%ktop_envel,                    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.TRUE., loutput=.FALSE.  )

    !        diag%snowlmt(nproma,nblks_c)
    cf_desc    = t_cf_var('snowlmt', 'm', 'Height of snow fall limit above MSL', &
      &                   datatype_flt)
    grib2_desc = grib2_var(0, 1, 204, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 101)
    CALL add_var( diag_list, 'snowlmt', diag%snowlmt,                          &
                & GRID_UNSTRUCTURED_CELL, ZA_ISOTHERM_ZERO, cf_desc,           &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE.,                  &
                & loutput=.TRUE.,                                              &
                & lmiss=.TRUE., missval=-999._wp,                              &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ) )


    ! &      diag%clc(nproma,nlev,nblks_c)
    cf_desc      = t_cf_var('clc', '',  'cloud cover', datatype_flt)
    new_cf_desc  = t_cf_var('clc', '%', 'cloud cover', datatype_flt)
    grib2_desc   = grib2_var(0, 6, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'clc', diag%clc,                                 &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,               &
      & ldims=shape3d, lrestart=.FALSE.,                                      &
      & in_group=groups("cloud_diag"),                                        &
      & vert_interp=create_vert_interp_metadata(                              &
      &             vert_intp_type=vintp_types("P","Z","I"),                  &
      &             vert_intp_method=VINTP_METHOD_LIN,                        &
      &             l_loglin=.FALSE.,                                         &
      &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,                   &
      &             lower_limit=0._wp ),                                      &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                          &
      &                 new_cf=new_cf_desc))



    ! &      diag%clct(nproma,nblks_c)
    cf_desc      = t_cf_var('clct', '',  'total cloud cover', datatype_flt)
    new_cf_desc  = t_cf_var('clct', '%', 'total cloud cover', datatype_flt)
    grib2_desc   = grib2_var(0, 6, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'clct', diag%clct,                               &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & in_group=groups("additional_precip_vars"),                            &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                          &
      &                 new_cf=new_cf_desc))

    ! &      diag%clct_mod(nproma,nblks_c)
    cf_desc      = t_cf_var('clct_mod', '', 'modified total cloud cover for media', &
      &                     datatype_flt)
    grib2_desc   = grib2_var(0, 6, 199, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'clct_mod', diag%clct_mod,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & in_group=groups("additional_precip_vars"),                            &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF))

    ! &      diag%clch(nproma,nblks_c)
    cf_desc      = t_cf_var('clch', '', 'high level clouds',  datatype_flt)
    new_cf_desc  = t_cf_var('clch', '%', 'high level clouds', datatype_flt)
    grib2_desc   = grib2_var(0, 6, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'clch', diag%clch,                               &
      & GRID_UNSTRUCTURED_CELL, ZA_PRESSURE_0, cf_desc, grib2_desc,           &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                          &
      &                 new_cf=new_cf_desc))

    ! &      diag%clcm(nproma,nblks_c)
    cf_desc      = t_cf_var('clcm', '',  'mid level clouds', datatype_flt)
    new_cf_desc  = t_cf_var('clcm', '%', 'mid level clouds', datatype_flt)
    grib2_desc   = grib2_var(0, 6, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'clcm', diag%clcm,                               &
      & GRID_UNSTRUCTURED_CELL, ZA_PRESSURE_400, cf_desc, grib2_desc,         &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                          &
      &                 new_cf=new_cf_desc))

    ! &      diag%clcl(nproma,nblks_c)
    cf_desc      = t_cf_var('clcl', '',  'low level clouds', datatype_flt)
    new_cf_desc  = t_cf_var('clcl', '%', 'low level clouds', datatype_flt)
    grib2_desc   = grib2_var(0, 6, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 1)
    CALL add_var( diag_list, 'clcl', diag%clcl,                               &
      & GRID_UNSTRUCTURED_CELL, ZA_PRESSURE_800, cf_desc, grib2_desc,         &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                          &
      &                 new_cf=new_cf_desc))

    ! &      diag%cldepth(nproma,nblks_c)
    cf_desc      = t_cf_var('cldepth', '',  'modified cloud depth for media', datatype_flt)
    grib2_desc   = grib2_var(0, 6, 198, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'cldepth', diag%cldepth,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF))


    ! &      diag%hbas_con(nproma,nblks_c)
    cf_desc    = t_cf_var('hbas_con', 'm', 'height of convective cloud base', datatype_flt)
    grib2_desc = grib2_var(0, 6, 26, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 101)
    CALL add_var( diag_list, 'hbas_con', diag%hbas_con,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE, cf_desc, grib2_desc,           &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
!!$      & lmiss=.TRUE., missval=-500._wp,                                       &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB))

    ! &      diag%htop_con(nproma,nblks_c)
    cf_desc    = t_cf_var('htop_con', 'm', 'height of convective cloud top', datatype_flt)
    grib2_desc = grib2_var(0, 6, 27, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 101)
    CALL add_var( diag_list, 'htop_con', diag%htop_con,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
!!$      & lmiss=.TRUE., missval=-500._wp,                                       &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB))

    ! &      diag%htop_dc(nproma,nblks_c)
    cf_desc    = t_cf_var('htop_dc', 'm', 'height of top of dry convection', datatype_flt)
    grib2_desc = grib2_var(0, 6, 196, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 101)
    CALL add_var( diag_list, 'htop_dc', diag%htop_dc,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
!!$      & lmiss=.TRUE., missval=-999._wp,                                       &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB))

    ! &      diag%acdnc(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('acdnc', 'm-3', 'cloud droplet number concentration', datatype_flt)
    grib2_desc = grib2_var(0, 6, 30, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'acdnc', diag%acdnc,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,               &
      & ldims=shape3d, lrestart=.FALSE.,                                      &
      & isteptype=TSTEP_CONSTANT )

    IF (ldass_lhn) THEN
      ! &      diag%tt_lheat(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('tt_lheat', 'K s-1',                &
        &          '3d latent heat relaese', DATATYPE_FLT32)
      grib2_desc = grib2_var(0, 255, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'tt_lheat', diag%tt_lheat,                       &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                  & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                  & ldims=shape3d,                                              &  
                  & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
                  &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
                  & isteptype=TSTEP_INSTANT )
  
      ! &      diag%qrs_flux(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('qrs_flux', 'kg m-2 s-1',                &
        &          '3d precipitation flux', DATATYPE_FLT32)
      grib2_desc = grib2_var(0, 255, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'qrs_flux', diag%qrs_flux,                       &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                  & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                  & ldims=shape3d,                                              &  
                  & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
                  &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
                  & isteptype=TSTEP_INSTANT )
  
      ! &      diag%lhn_diag(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('lhn_diag', '-',                &
        &          'diagnose of LHN', DATATYPE_FLT32)
      grib2_desc = grib2_var(0, 255, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'lhn_diag', diag%lhn_diag,                       &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                  & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                  & ldims=shape3d,                                              &  
                  & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB), &
                  & isteptype=TSTEP_INSTANT )
  
      ! &      diag%lhn_diag(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('ttend_lhn', 'K s-1',                &
        &          'tempature increment', DATATYPE_FLT32)
      grib2_desc = grib2_var(0, 255, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'ttend_lhn', diag%ttend_lhn,                       &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                  & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                  & ldims=shape3d,                                              &  
                  & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
                  &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
                  & isteptype=TSTEP_INSTANT )
      ! &      diag%lhn_diag(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('qvtend_lhn', 'g/g s-1',                &
        &          'qv increment', DATATYPE_FLT32)
      grib2_desc = grib2_var(0, 255, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'qvtend_lhn', diag%qvtend_lhn,                       &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                  & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                  & ldims=shape3d,                                              &  
                  & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
                  &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
                  & isteptype=TSTEP_INSTANT )
    ELSE
      ALLOCATE (diag%qrs_flux(1,1,kblks))
    ENDIF


    !      diag%tot_cld(nproma,nlev,nblks_c,3)
    cf_desc    = t_cf_var('tot_cld', ' ','total cloud variables (qv,qc,qi)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tot_cld', diag%tot_cld,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=(/nproma,klev,kblks,3/) ,                               &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,         &
                & initval=0.0_wp,                                               &
                & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
                &                                       fallback_type=HINTP_TYPE_LONLAT_RBF))

    ALLOCATE( diag%tot_ptr(kcloud))
    vname_prefix='tot_'

    !QV
    CALL add_ref( diag_list, 'tot_cld',                                            &
                & TRIM(vname_prefix)//'qv_dia', diag%tot_ptr(iqv)%p_3d,            &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                               &
                & t_cf_var(TRIM(vname_prefix)//'qv_dia', 'kg kg-1',                &
                &          'total specific humidity (diagnostic)', datatype_flt),&
                & grib2_var(0, 1, 211, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                & ldims=shape3d,                                                   &
                & vert_interp=create_vert_interp_metadata(                         &
                &             vert_intp_type=vintp_types("P", "Z", "I"),           &
                &             vert_intp_method=VINTP_METHOD_QV,                    &
                &             l_satlimit=.FALSE.,                                  & 
                &             lower_limit=2.5e-7_wp, l_restore_pbldev=.FALSE. ),   &
                & in_group=groups("cloud_diag") )

    !QC
    CALL add_ref( diag_list, 'tot_cld',                                            &
                & TRIM(vname_prefix)//'qc_dia', diag%tot_ptr(iqc)%p_3d,            &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                               &
                & t_cf_var(TRIM(vname_prefix)//'qc_dia', 'kg kg-1',                &
                & 'total specific cloud water content (diagnostic)', datatype_flt),&
                & grib2_var(0, 1, 212, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                & ldims=shape3d,                                                   &
                & vert_interp=create_vert_interp_metadata(                         &
                &             vert_intp_type=vintp_types("P","Z","I"),             &
                &             vert_intp_method=VINTP_METHOD_LIN,                   &
                &             l_loglin=.FALSE.,                                    &
                &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,              &
                &             lower_limit=0._wp ),                                 &
                & in_group=groups("cloud_diag") )

    !QI
    CALL add_ref( diag_list, 'tot_cld',                                            &
                & TRIM(vname_prefix)//'qi_dia', diag%tot_ptr(iqi)%p_3d,            &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                               &
                & t_cf_var(TRIM(vname_prefix)//'qi_dia', 'kg kg-1',                &
                & 'total specific cloud ice content (diagnostic)', datatype_flt),&
                & grib2_var(0, 1, 213, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                & ldims=shape3d,                                                   &
                & vert_interp=create_vert_interp_metadata(                         &
                &             vert_intp_type=vintp_types("P","Z","I"),             &
                &             vert_intp_method=VINTP_METHOD_LIN,                   &
                &             l_loglin=.FALSE.,                                    &
                &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,              &
                &             lower_limit=0._wp ),                                 &
                & in_group=groups("cloud_diag") )

    !      diag%tot_cld_vi(nproma,nblks_c,3)
    cf_desc     = t_cf_var('tot_cld_vi', 'kg m-2','vertical integr total cloud variables', datatype_flt)
    grib2_desc   = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tot_cld_vi', diag%tot_cld_vi,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                &                                 ldims=(/nproma,kblks,3/),   &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.        )

    ! fill the seperate variables belonging to the container tot_cld_vi
    ALLOCATE( diag%tci_ptr(3))
       
    !TQV_DIA
    CALL add_ref( diag_list, 'tot_cld_vi',                        &
      & 'tqv_dia', diag%tci_ptr(iqv)%p_2d,                        &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                       &
      & t_cf_var('tqv_dia', 'kg m**-2',                           &
      & 'total column integrated water vapour (diagnostic)', datatype_flt),   &
      & grib2_var( 0, 1, 214, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("additional_precip_vars"))

    !TQC_DIA
    CALL add_ref( diag_list, 'tot_cld_vi',                         &
      & 'tqc_dia', diag%tci_ptr(iqc)%p_2d,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
      & t_cf_var('tqc_dia', 'kg m**-2',                            &
      & 'total column integrated cloud water (diagnostic)', datatype_flt),    &
      & grib2_var( 0, 1, 215, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("additional_precip_vars"))

    !TQI_DIA
    CALL add_ref( diag_list, 'tot_cld_vi',                         &
      & 'tqi_dia', diag%tci_ptr(iqi)%p_2d,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
      & t_cf_var('tqi_dia', 'kg m**-2',                            &
      & 'total column integrated cloud ice (diagnostic)', datatype_flt),      &
      & grib2_var(0, 1, 216, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("additional_precip_vars"))



    !      diag%tot_cld_vi_avg(nproma,nblks_c,3)
    cf_desc    = t_cf_var('tot_cld_vi_avg', 'unit ','vertical integr total cloud variables', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tot_cld_vi_avg', diag%tot_cld_vi_avg,           &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                &                                 ldims=(/nproma,kblks,3/)  , &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,       &
                & isteptype=TSTEP_AVG )

    ! fill the seperate variables belonging to the container tot_cld_vi_avg
    ALLOCATE( diag%tav_ptr(3))
    vname_prefix='avg_'


    ! TQV_DIA_AVG
    CALL add_ref( diag_list, 'tot_cld_vi_avg',               &
                & TRIM(vname_prefix)//'qv', diag%tav_ptr(iqv)%p_2d,            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var(TRIM(vname_prefix)//'qv', 'kg m-2',                 &
                & 'column integrated water vapour (diagnostic) avg',           &
                & datatype_flt),                                               &
                & grib2_var( 0, 1, 214, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
                & ldims=shape2d, lrestart=.FALSE.,                             &
                & isteptype=TSTEP_AVG,                                         & 
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                     &
                &    fallback_type=HINTP_TYPE_LONLAT_RBF                       &
                & ) )

    ! TQC_DIA_AVG
    CALL add_ref( diag_list, 'tot_cld_vi_avg',      &
                & TRIM(vname_prefix)//'qc', diag%tav_ptr(iqc)%p_2d,            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var(TRIM(vname_prefix)//'qc', '',                       &
                & 'tci specific cloud water content (diagnostic) avg',         &
                & datatype_flt),                                               &
                & grib2_var(0, 1, 215, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
                & ldims=shape2d, lrestart=.FALSE.,                             &
                & isteptype=TSTEP_AVG,                                         &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                     &
                &    fallback_type=HINTP_TYPE_LONLAT_RBF                       &
                & ) )

    ! TQI_DIA_AVG
    CALL add_ref( diag_list, 'tot_cld_vi_avg',          &
                & TRIM(vname_prefix)//'qi', diag%tav_ptr(iqi)%p_2d,            & 
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var(TRIM(vname_prefix)//'qi', '',                       &
                & 'tci specific cloud ice content (diagnostic) avg',           &
                & datatype_flt),                                               &
                & grib2_var(0, 1, 216, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
                & ldims=shape2d, lrestart=.FALSE.,                             &
                & isteptype=TSTEP_AVG,                                         & 
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                     &
                &    fallback_type=HINTP_TYPE_LONLAT_RBF                       &
                & ) )


    ! &      diag%clct_avg(nproma,nblks_c)
    cf_desc      = t_cf_var('clct_avg', '%', 'total cloud cover time avg', &
      &            datatype_flt)
    grib2_desc   = grib2_var(0, 6, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'clct_avg', diag%clct_avg,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF))




    !------------------
    ! Radiation
    !------------------
    ! 2D variables

    !        diag%albdif    (nproma,       nblks),          &
    cf_desc     = t_cf_var('albdif', '', 'Shortwave albedo for diffuse radiation', &
      &                    datatype_flt)
    new_cf_desc = t_cf_var('albdif', '%','Shortwave albedo for diffuse radiation', &
      &                    datatype_flt)
    grib2_desc  = grib2_var(0, 19, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'albdif', diag%albdif,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, in_group=groups("rad_vars"),                         &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                        &
      &                 new_cf=new_cf_desc) )

    !        diag%albvisdif    (nproma,       nblks),          &
    cf_desc     = t_cf_var('albvisdif', '', 'UV visible albedo for diffuse radiation', &
      &                    datatype_flt)
    new_cf_desc = t_cf_var('albvisdif', '%','UV visible albedo for diffuse radiation', &
      &                    datatype_flt)
    grib2_desc  = grib2_var(0, 19, 222, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'albvisdif', diag%albvisdif,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, in_group=groups("rad_vars"),                         &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                        &
      &                 new_cf=new_cf_desc) )

    !        diag%albvisdir    (nproma,       nblks),          &
    cf_desc     = t_cf_var('albvisdir', '', 'UV visible albedo for direct radiation', &
      &                    datatype_flt)
    new_cf_desc = t_cf_var('albvisdir', '%','UV visible albedo for direct radiation', &
      &                    datatype_flt)
    grib2_desc  = grib2_var(192, 128, 15, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'albvisdir', diag%albvisdir,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d,                                                      &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                        &
      &                 new_cf=new_cf_desc) )

    !        diag%albnirdif    (nproma,       nblks),          &
    cf_desc     = t_cf_var('albnirdif', '',  'Near IR albedo for diffuse radiation',&
      &                    datatype_flt)
    new_cf_desc = t_cf_var('albnirdif', '%', 'Near IR albedo for diffuse radiation',&
      &                    datatype_flt)
    grib2_desc  = grib2_var(0, 19, 223, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'albnirdif', diag%albnirdif,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, in_group=groups("rad_vars"),                         &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                        &
      &                 new_cf=new_cf_desc) )

    !        diag%albnirdir    (nproma,       nblks),          &
    cf_desc     = t_cf_var('albnirdir', '',  'Near IR albedo for direct radiation',&
      &                    datatype_flt)
    new_cf_desc = t_cf_var('albnirdir', '%', 'Near IR albedo for direct radiation',&
      &                    datatype_flt)
    grib2_desc  = grib2_var(192, 128, 17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'albnirdir', diag%albnirdir,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d,                                                      &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                        &
      &                 new_cf=new_cf_desc) )


    ! longwave surface emissivity
    !
    ! lw_emiss     diag%lw_emiss(nproma,nblks_c)
    cf_desc    = t_cf_var('lw_emiss', '-', 'longwave surface emissivity', datatype_flt)
    grib2_desc = grib2_var( 2, 3, 199, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lw_emiss', diag%lw_emiss,         &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
      &           grib2_desc, ldims=shape2d, loutput=.TRUE.,    &
      &           lrestart=.TRUE. )


    ! These variables only make sense if the land-surface scheme is switched on.
    IF ( atm_phy_nwp_config(k_jg)%inwp_surface == 1 ) THEN

      !        diag%albdif_t (nproma, nblks, ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('albdif_t', '', &
        &                   'tile-based shortwave albedo for diffusive radiation',&
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 19, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'albdif_t', diag%albdif_t,                &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,       &
        & ldims=shape3dsubsw, lcontainer=.TRUE., lrestart=.FALSE.,       &
        & loutput=.FALSE.)

      ! fill the seperate variables belonging to the container albdif_t
      ALLOCATE(diag%albdif_t_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total+ntiles_water
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( diag_list, 'albdif_t',                               &
           & 'albdif_t_'//TRIM(ADJUSTL(csfc)),                             &
           & diag%albdif_t_ptr(jsfc)%p_2d,                                 &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
           & t_cf_var('albdif_t_'//TRIM(csfc), '', '', datatype_flt),    &
           & grib2_var(0, 19, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
           & ldims=shape2d, lrestart=.TRUE.                                )
      ENDDO

      !        diag%albvisdif_t (nproma, nblks, ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('albvisdif_t', '', &
        &                   'tile-based UV visible albedo for diffusive radiation',&
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 19, 222, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'albvisdif_t', diag%albvisdif_t,               &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
        & ldims=shape3dsubsw, lcontainer=.TRUE., lrestart=.FALSE., &
        & loutput=.FALSE.)

      ! fill the seperate variables belonging to the container albvisdif_t
      ALLOCATE(diag%albvisdif_t_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total+ntiles_water
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( diag_list, 'albvisdif_t',                            &
           & 'albvisdif_t_'//TRIM(ADJUSTL(csfc)),                          &
           & diag%albvisdif_t_ptr(jsfc)%p_2d,                              &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
           & t_cf_var('albvisdif_t_'//TRIM(csfc), '', '', datatype_flt), &
           & grib2_var(0, 19, 222, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
           & ldims=shape2d, lrestart=.TRUE.                                )
      ENDDO

      !        diag%albnirdif_t (nproma, nblks, ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('albnirdif_t', '', &
        &                   'tile-based near IR albedo for diffuse radiation',&
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 19, 223, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'albnirdif_t', diag%albnirdif_t,               &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
        & ldims=shape3dsubsw, lcontainer=.TRUE., lrestart=.FALSE., &
        & loutput=.FALSE.)


      ! fill the seperate variables belonging to the container albnirdif_t
      ALLOCATE(diag%albnirdif_t_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total+ntiles_water
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( diag_list, 'albnirdif_t',                            &
           & 'albnirdif_t_'//TRIM(ADJUSTL(csfc)),                          &
           & diag%albnirdif_t_ptr(jsfc)%p_2d,                              &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
           & t_cf_var('albnirdif_t_'//TRIM(csfc), '', '', datatype_flt), &
           & grib2_var(0, 19, 223, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
           & ldims=shape2d, lrestart=.TRUE.                                )
      ENDDO


      ! &      diag%swflxsfc_t(nproma,nblks_c,ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('sob_s_t', 'W m-2', 'tile-based shortwave net flux at surface', &
           &                datatype_flt)
      grib2_desc = grib2_var(0, 4, 9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'sob_s_t', diag%swflxsfc_t,                 &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
        & ldims=shape3dsubsw, lcontainer=.TRUE., lrestart=.FALSE., &
        & loutput=.FALSE.)

      ! fill the seperate variables belonging to the container swflxsfc_t
      ALLOCATE(diag%swflxsfc_t_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total+ntiles_water
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( diag_list, 'sob_s_t',                                &
           & 'sob_s_t_'//TRIM(ADJUSTL(csfc)),                              &
           & diag%swflxsfc_t_ptr(jsfc)%p_2d,                               &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
           & t_cf_var('swflxsfc_t_'//TRIM(csfc), '', '', datatype_flt),  &
           & grib2_var(0, 4, 9, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
           & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,               &
           & in_group=groups("rad_vars"))
      ENDDO

      ! &      diag%lwflxsfc_t(nproma,nblks_c,ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('thb_s_t', 'W m-2', 'tile-based longwave net flux at surface', &
           &                datatype_flt)
      grib2_desc = grib2_var(0, 5, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'thb_s_t', diag%lwflxsfc_t,                 &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
        & ldims=shape3dsubsw, lcontainer=.TRUE., lrestart=.FALSE.,         &
        & loutput=.FALSE.)

      ! fill the seperate variables belonging to the container lwflxsfc_t
      ALLOCATE(diag%lwflxsfc_t_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total+ntiles_water
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( diag_list, 'thb_s_t',                                &
           & 'thb_s_t_'//TRIM(ADJUSTL(csfc)),                              &
           & diag%lwflxsfc_t_ptr(jsfc)%p_2d,                               &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
           & t_cf_var('lwflxsfc_t_'//TRIM(csfc), '', '', datatype_flt),  &
           & grib2_var(0, 5, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
           & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,               &
           & in_group=groups("rad_vars"))
      ENDDO
    ENDIF

    !        diag%cosmu0    (nproma,nblks)
    cf_desc    = t_cf_var('cosmu0', '-', 'Cosine of solar zenith angle', datatype_flt)
    grib2_desc = grib2_var(192, 214, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'cosmu0', diag%cosmu0,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d)

    ! &      diag%tsfctrad(nproma,nblks_c)
    cf_desc    = t_cf_var('tsfctrad', 'K', 'surface temperature at trad', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tsfctrad', diag%tsfctrad,                     &
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d)

    ! &       diag% flxdwswtoa(nproma,       nblks),          &
    cf_desc    = t_cf_var('sod_t', 'W m-2', 'downward shortwave flux at TOA', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 4, 7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sod_t', diag%flxdwswtoa,              &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,             &
      & ldims=shape2d, lrestart=.FALSE.,in_group=groups("rad_vars"))

    ! &      diag%swflxsfc(nproma,nblks_c)
    cf_desc    = t_cf_var('sob_s', 'W m-2', 'shortwave net flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sob_s', diag%swflxsfc,                        &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d,                                                      &
      & in_group=groups("rad_vars"))

    ! &      diag%lwflxclr_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('thbclr_s', 'W m-2', 'net longwave clear-sky flux at suface', datatype_flt)
    grib2_desc = grib2_var(0, 5, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'thbclr_s', diag%lwflxclr_sfc,             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d)

    ! &      diag%swflxclr_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('sobclr_s', 'W m-2', 'net shortwave clear-sky flux at suface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sobclr_s', diag%swflxclr_sfc,             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,  ldims=shape2d)

    ! &      diag%swflxtoa(nproma,nblks_c)
    cf_desc    = t_cf_var('sob_t', 'W m-2', 'shortwave net flux at TOA', datatype_flt)
    grib2_desc = grib2_var(0, 4, 9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sob_t', diag%swflxtoa,                        &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,                &
      & ldims=shape2d, lrestart=.FALSE.,                                    &
      & in_group=groups("rad_vars"))

    ! &      diag%trsolclr_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('trsolclr_sfc', '', 'shortwave clear-sky transmisivity at suface', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'trsolclr_sfc', diag%trsolclr_sfc,             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, loutput=.FALSE.                                      )

    ! &      diag%trsol_up_toa(nproma,nblks_c)
    cf_desc    = t_cf_var('trsol_up_toa', '', 'shortwave upward transmisivity at TOA', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'trsol_up_toa', diag%trsol_up_toa,             &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,                &
      & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE.                    )

    ! &      diag%trsol_up_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('trsol_up_sfc', '', 'shortwave upward transmisivity at surface', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'trsol_up_sfc', diag%trsol_up_sfc,             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE.                    )

    ! &      diag%trsol_par_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('trsol_par_sfc', '', 'photosynthetically active downward transmisivity at surface', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'trsol_par_sfc', diag%trsol_par_sfc,           &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.TRUE., loutput=.FALSE.                     )

    ! &      diag%trsol_dn_sfc_diff(nproma,nblks_c)
    cf_desc    = t_cf_var('trsol_dn_sfc_diff', '', 'shortwave diffuse downward transmisivity at surface', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'trsol_dn_sfc_diff', diag%trsol_dn_sfc_diff,   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.TRUE., loutput=.FALSE.                    )

    ! &      diag%swflx_up_toa(nproma,nblks_c)
    cf_desc    = t_cf_var('sou_t', 'W m-2', 'shortwave upward flux at TOA', datatype_flt)
    grib2_desc = grib2_var(0, 4, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sou_t', diag%swflx_up_toa,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,                &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("rad_vars")        )

    ! &      diag%swflx_up_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('sou_s', 'W m-2', 'shortwave upward flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sou_s', diag%swflx_up_sfc,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("rad_vars")        )

    ! &      diag%swflx_par_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('swflx_par_sfc', 'W m-2', 'downward photosynthetically active flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'swflx_par_sfc', diag%swflx_par_sfc,           &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.TRUE., in_group=groups("rad_vars")         )

    ! &      diag%swflx_dn_sfc_diff(nproma,nblks_c)
    cf_desc    = t_cf_var('sodifd_s', 'W m-2', 'shortwave diffuse downward flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 199, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sodifd_s', diag%swflx_dn_sfc_diff,              &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("rad_vars")        )

    ! &      diag%lwflxsfc(nproma,nblks_c)
    cf_desc    = t_cf_var('thb_s', 'W m-2', 'longwave net flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 5, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'thb_s', diag%lwflxsfc,                        &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d,                                                      &
      & in_group=groups("rad_vars"))

    ! This is an auxiliary storage field needed because the final diagnostic quantity (below)
    ! is updated each physics time step following the time evolution of the ground temperature
    ! &      diag%lwflx_up_sfc_rs(nproma,nblks_c)
    cf_desc    = t_cf_var('lwflx_up_sfc_rs', 'W m-2', 'longwave upward flux at surface', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lwflx_up_sfc_rs', diag%lwflx_up_sfc_rs,       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d,  lrestart=.TRUE., loutput=.FALSE.                   )

    ! &      diag%lwflx_up_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('thu_s', 'W m-2', 'longwave upward flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 5, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'thu_s', diag%lwflx_up_sfc,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d,  lrestart=.FALSE.,                                   &
      & in_group=groups("rad_vars"))

    IF (lflux_avg ) THEN
        prefix = "a"
        meaning = "mean"
        varunits= "W/m**2"
        a_steptype= TSTEP_AVG
    ELSE
        prefix = "acc"
        meaning = "acc." 
        varunits= "J/m**2"
        a_steptype= TSTEP_ACCUM
    END IF
    WRITE(name,'(A,A5)') TRIM(prefix),"thb_t"
    WRITE(long_name,'(A26,A4,A18)') "TOA net thermal radiation ", meaning, &
                                  & " since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), &
      &                          TRIM(long_name), datatype_flt)
    grib2_desc = grib2_var(0, 5, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%lwflxtoa_a,                  &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,               &
      & ldims=shape2d, isteptype=a_steptype, in_group=groups("rad_vars"),  &
      & hor_interp=create_hor_interp_metadata(                             &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                           &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )


    ! &      diag%swflxtoa_a(nproma,nblks_c)
    WRITE(name,'(A,A5)') TRIM(prefix),"sob_t"
    WRITE(long_name,'(A26,A4,A18)') "TOA net solar radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), &
      &                         TRIM(long_name), datatype_flt)
    grib2_desc = grib2_var(0, 4, 9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name) , diag%swflxtoa_a,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,               &
      & ldims=shape2d,                                                     &
      & isteptype=a_steptype, in_group=groups("rad_vars"),                 &
      & hor_interp=create_hor_interp_metadata(                             &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                           &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )
    
    ! &      diag%lwflxsfc_a(nproma,nblks_c)
    WRITE(name,'(A,A5)') TRIM(prefix),"thb_s"
    WRITE(long_name,'(A30,A4,A18)') "surface net thermal radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), &
      &                      TRIM(long_name), datatype_flt)
    grib2_desc = grib2_var(0, 5, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%lwflxsfc_a,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & isteptype=a_steptype, in_group=groups("rad_vars"),                &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    ! &      diag%lwflxclrsfc_a(nproma,nblks_c)
    WRITE(name,'(A,A8)') TRIM(prefix),"thbclr_s"
    WRITE(long_name,'(A40,A4,A18)') "clear-sky surface net thermal radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), &
      &                      TRIM(long_name), datatype_flt)
    grib2_desc = grib2_var(0, 5, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%lwflxclrsfc_a,              &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, isteptype=a_steptype,                              &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    ! &      diag%swflxsfc_a(nproma,nblks_c)
    WRITE(name,'(A,A5)') TRIM(prefix),"sob_s"
    WRITE(long_name,'(A30,A4,A18)') "Surface net solar radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), &
      &                    TRIM(long_name), datatype_flt)
    grib2_desc = grib2_var(0, 4, 9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%swflxsfc_a,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & isteptype=a_steptype, in_group=groups("rad_vars"),                &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    ! &      diag%swflxclrsfc_a(nproma,nblks_c)
    WRITE(name,'(A,A8)') TRIM(prefix),"sobclr_s"
    WRITE(long_name,'(A40,A4,A18)') "Clear-sky surface net solar radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), &
      &                    TRIM(long_name), datatype_flt)
    grib2_desc = grib2_var(0, 4, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%swflxclrsfc_a,              &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, isteptype=a_steptype,                              &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    ! &       diag%asod_t(nproma,nblks)
    WRITE(name,'(A,A5)') TRIM(prefix),"sod_t"
    WRITE(long_name,'(A30,A4,A18)') "Top down solar radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name),    &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 4, 7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%asod_t,                     &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,              &
      & ldims=shape2d,                                                    &
      & isteptype=a_steptype, in_group=groups("rad_vars"),                &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    ! &       diag%asou_t(nproma,nblks)
    WRITE(name,'(A,A5)') TRIM(prefix),"sou_t"
    WRITE(long_name,'(A30,A4,A18)') "Top up solar radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 4, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%asou_t,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,             &
      & ldims=shape2d,                                                   &
      & isteptype=a_steptype, in_group=groups("rad_vars"),               &
      & hor_interp=create_hor_interp_metadata(                           &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                         &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    ! &       diag%athd_s(nproma,nblks)
    WRITE(name,'(A,A5)') TRIM(prefix),"thd_s"
    WRITE(long_name,'(A30,A4,A18)') "Surface down thermal radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 5, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%athd_s,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
      & ldims=shape2d,                                                   &
      & isteptype=a_steptype, in_group=groups("rad_vars"),               &
      & hor_interp=create_hor_interp_metadata(                           &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                         &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    ! &       diag%athu_s(nproma,nblks)
    WRITE(name,'(A,A5)') TRIM(prefix),"thu_s"
    WRITE(long_name,'(A30,A4,A18)') "Surface up thermal radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 5, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%athu_s,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
      & ldims=shape2d,                                                   &
      & isteptype=a_steptype, in_group=groups("rad_vars"),               &
      & hor_interp=create_hor_interp_metadata(                           &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                         &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    ! &       diag%asodird_s(nproma,nblks)
    WRITE(name,'(A,A8)') TRIM(prefix),"sodird_s"
    WRITE(long_name,'(A30,A4,A18)') "Surface down solar direct rad. ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 4, 198, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%asodird_s,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
      & ldims=shape2d,                                                   &
      & isteptype=a_steptype, in_group=groups("rad_vars"),               &
      & hor_interp=create_hor_interp_metadata(                           &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                         &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    ! &       diag%asodifd_s(nproma,nblks)
    WRITE(name,'(A,A8)') TRIM(prefix),"sodifd_s"
    WRITE(long_name,'(A30,A4,A18)') "Surface down solar diff. rad. ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 4, 199, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%asodifd_s,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
      & ldims=shape2d,                                                   &
      & isteptype=a_steptype, in_group=groups("rad_vars"),               &
      & hor_interp=create_hor_interp_metadata(                           &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                         &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    ! &       diag%asod_s(nproma,nblks)
    WRITE(name,'(A,A5)') TRIM(prefix),"sod_s"
    WRITE(long_name,'(A30,A4,A18)') "Surface down solar rad. ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 4, 7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%asod_s,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
      & ldims=shape2d, isteptype=a_steptype,                             &
      & hor_interp=create_hor_interp_metadata(                           &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                         &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    ! &       diag%asodifu_s(nproma,nblks)
    WRITE(name,'(A,A8)') TRIM(prefix),"sodifu_s"
    WRITE(long_name,'(A30,A4,A18)') "Surface up solar diff. rad. ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 4, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%asodifu_s,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
      & ldims=shape2d,                                                   &
      & isteptype=a_steptype, in_group=groups("rad_vars"),               &
      & hor_interp=create_hor_interp_metadata(                           &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                         &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    ! &      diag%aswflx_par_sfc(nproma,nblks_c)
    WRITE(name,'(A,A13)') TRIM(prefix),"swflx_par_sfc"
    WRITE(long_name,'(A30,A4,A18)') "Downward PAR flux ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), &
      &                   datatype_flt)
    grib2_desc = grib2_var(0, 4, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%aswflx_par_sfc,            &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
      & ldims=shape2d,                                                   & 
      & isteptype=a_steptype, in_group=groups("rad_vars"),               &
      & hor_interp=create_hor_interp_metadata(                           &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                         &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    ! &      diag%vio3(nproma,nblks_c)
    cf_desc    = t_cf_var('vio3', '', 'vertically integrated ozone amount', datatype_flt)
    grib2_desc = grib2_var(0, 14, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'vio3', diag%vio3,                           &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.FALSE. )

    ! &      diag%hmo3(nproma,nblks_c)
    cf_desc    = t_cf_var('hmo3', 'Pa', 'height of O3 maximum (Pa)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'hmo3', diag%hmo3,                           &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.FALSE. )

    ! Reference pressure used for vertical distribution of aerosol optical depths
    ! &      diag%pref_aerdis(nproma,nblks_c)
    cf_desc    = t_cf_var('pref_aerdis', '',                                                             &
      &                   'Reference pressure used for vertical distribution of aerosol optical depths', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'pref_aerdis', diag%pref_aerdis,              &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.FALSE. ) 


    IF (irad_aero == 5) THEN ! Old Tanre aerosol climatology taken over from the COSMO model (to be used with inwp_radiation==2)
      ! &      diag%aersea(nproma,nblks_c)
      cf_desc    = t_cf_var('aersea', '', '', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'aersea', diag%aersea,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape2d, lrestart=.FALSE. ) 

      ! &      diag%aerlan(nproma,nblks_c)
      cf_desc    = t_cf_var('aerlan', '', '', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'aerlan', diag%aerlan,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape2d, lrestart=.FALSE. ) 
    
      ! &      diag%aerurb(nproma,nblks_c)
      cf_desc    = t_cf_var('aerurb', '', '', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'aerurb', diag%aerurb,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape2d, lrestart=.FALSE. ) 

    
      ! &      diag%aerdes(nproma,nblks_c)
      cf_desc    = t_cf_var('aerdes', '', '', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'aerdes', diag%aerdes,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape2d, lrestart=.FALSE. ) 

    ELSE IF (irad_aero == 6 .OR. irad_aero == 9) THEN ! Tegen aerosol climatology, time-interpolated values 
                                  ! (needed as state fields for coupling with microphysics and convection)
      IF (atm_phy_nwp_config(k_jg)%icpl_aero_gscp > 1 .OR. icpl_aero_conv > 1 .OR. iprog_aero > 0) THEN
        lrestart = .TRUE.
      ELSE
        lrestart = .FALSE.
      ENDIF

      ! &      diag%aercl_ss(nproma,nblks_c)
      cf_desc    = t_cf_var('aercl_ss', '', 'sea salt aerosol climatology', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'aercl_ss', diag%aercl_ss,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape2d, lrestart=lrestart ) 

      ! &      diag%aercl_or(nproma,nblks_c)
      cf_desc    = t_cf_var('aercl_or', '', 'organic aerosol climatology', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'aercl_or', diag%aercl_or,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape2d, lrestart=lrestart ) 

      ! &      diag%aercl_bc(nproma,nblks_c)
      cf_desc    = t_cf_var('aercl_bc', '', 'black carbon aerosol climatology', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'aercl_bc', diag%aercl_bc,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape2d, lrestart=lrestart ) 

      ! &      diag%aercl_su(nproma,nblks_c)
      cf_desc    = t_cf_var('aercl_su', '', 'total sulfate aerosol climatology', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'aercl_su', diag%aercl_su,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape2d, lrestart=lrestart ) 

      ! &      diag%aercl_du(nproma,nblks_c)
      cf_desc    = t_cf_var('aercl_du', '', 'total soil dust aerosol climatology', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'aercl_du', diag%aercl_du,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape2d, lrestart=lrestart ) 

      ! &      diag%aerosol(nproma,nclass_aero,nblks_c)
      cf_desc    = t_cf_var('aerosol', '', '', DATATYPE_FLT32)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'aerosol', diag%aerosol,                     &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape3d_aero, lrestart=.FALSE., lcontainer=.TRUE., loutput=.FALSE. ) 

      ALLOCATE(diag%aerosol_ptr(nclass_aero))
      DO k = 1, nclass_aero
        SELECT CASE (k)
        CASE (iss)
          caer='ss'
          constituentType = 62008
          cf_desc    = t_cf_var('aer_'//TRIM(caer), '', &
            &                   'sea salt aerosol', DATATYPE_FLT32)
        CASE (iorg)
          caer='or'
          constituentType = 62010
          cf_desc    = t_cf_var('aer_'//TRIM(caer), '', &
            &                   'organic aerosol', DATATYPE_FLT32)
        CASE (ibc)
          caer='bc'
          constituentType = 62009
          cf_desc    = t_cf_var('aer_'//TRIM(caer), '', &
            &                   'black carbon aerosol', DATATYPE_FLT32)
        CASE (iso4)
          caer='su'
          constituentType = 62006
          cf_desc    = t_cf_var('aer_'//TRIM(caer), '', &
            &                   'total sulfate aerosol', DATATYPE_FLT32)
        CASE (idu)
          caer='du'
          constituentType = 62001
          cf_desc    = t_cf_var('aer_'//TRIM(caer), '', &
            &                   'total soil dust aerosol', DATATYPE_FLT32)
        END SELECT

        grib2_desc = grib2_var(0, 20, 102, ibits, GRID_UNSTRUCTURED, GRID_CELL)   &
          &           + t_grib2_int_key("constituentType", constituentType)
        IF (iprog_aero > 1 .OR. iprog_aero == 1 .AND. k == idu) THEN
          CALL add_ref( diag_list, 'aerosol',                                    &
                & 'aer_'//TRIM(caer), diag%aerosol_ptr(k)%p_2d,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
                & cf_desc,                                                       &
                & grib2_desc,                                                    &
                & ldims=shape2d, lrestart=lrestart, opt_var_ref_pos = 2,         &
                & var_class=CLASS_CHEM,                                          &
                & in_group=groups("dwd_fg_sfc_vars","mode_iau_fg_in",            &
                & "mode_iau_old_fg_in","mode_dwd_fg_in")                         )
        ELSE
          CALL add_ref( diag_list, 'aerosol',                                    &
                & 'aer_'//TRIM(caer), diag%aerosol_ptr(k)%p_2d,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
                & cf_desc,                                                       &
                & grib2_desc,                                                    &
                & ldims=shape2d, lrestart=lrestart, opt_var_ref_pos = 2,         &
                & var_class=CLASS_CHEM                                           )
        ENDIF
      ENDDO

    ENDIF

    IF ( (irad_aero == 6 .OR. irad_aero == 9) .AND.  &
      &  (atm_phy_nwp_config(k_jg)%icpl_aero_gscp == 1 .OR. icpl_aero_conv == 1) ) THEN
      lrestart = .TRUE.
    ELSE
      lrestart = .FALSE.
    ENDIF

    ! &      diag%cloud_num(nproma,nblks_c)
    cf_desc    = t_cf_var('cloud_num', 'm-3', 'cloud droplet number concentration', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'cloud_num', diag%cloud_num,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=lrestart )

    !------------------
    !Radiation 3D variables

    ! &      diag%lwflxall(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('lwflxall', 'W m-2 ', 'longwave net flux', datatype_flt)
    grib2_desc = grib2_var(0, 5, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lwflxall', diag%lwflxall,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,      &
      & ldims=shape3dkp1 )

    ! &      diag%trsolall(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('trsolall', '', 'shortwave net tranmissivity', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'trsolall', diag%trsolall,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, ldims=shape3dkp1)

    !------------------
    !Turbulence 2D variables
    
    ! &      diag%shfl_s(nproma,nblks_c)
    cf_desc    = t_cf_var('shfl_s', 'W m-2 ', 'surface sensible heat flux', datatype_flt)
    grib2_desc = grib2_var(0, 0, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'shfl_s', diag%shfl_s,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & in_group=groups("pbl_vars") )

    WRITE(name,'(A,A6)') TRIM(prefix),"shfl_s"
    WRITE(long_name,'(A27,A4,A18)') "surface sensible heat flux ", meaning, &
                                  & " since model start"
    cf_desc    = t_cf_var(TRIM(name),  TRIM(varunits), TRIM(long_name), datatype_flt)
    grib2_desc = grib2_var(0, 0, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%ashfl_s,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & isteptype=a_steptype, in_group=groups("pbl_vars"),                &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )


    ! &      diag%lhfl_s(nproma,nblks_c)
    cf_desc    = t_cf_var('lhfl_s', 'W m-2 ', 'surface latent heat flux', datatype_flt)
    grib2_desc = grib2_var(0, 0, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lhfl_s', diag%lhfl_s,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & in_group=groups("pbl_vars"))
                
    WRITE(name,'(A,A6)') TRIM(prefix),"lhfl_s"
    WRITE(long_name,'(A27,A4,A18)') "surface latent heat flux ", meaning, &
                                  & " since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), datatype_flt)
    grib2_desc = grib2_var(0, 0, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%alhfl_s,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & isteptype=a_steptype, in_group=groups("pbl_vars"),                &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )


    ! &      diag%lhfl_bs(nproma,nblks_c)
    cf_desc    = t_cf_var('lhfl_bs', 'W m-2 ', 'latent heat flux from bare soil', &
      &          datatype_flt)
    grib2_desc = grib2_var(2, 0, 193, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lhfl_bs', diag%lhfl_bs,                     &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.FALSE.,                                  &
      & in_group=groups("pbl_vars"))

    WRITE(name,'(A,A7)') TRIM(prefix),"lhfl_bs"
    WRITE(long_name,'(A27,A4,A18)') "latent heat flux from bare soil", meaning, &
                                  & " since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), datatype_flt)
    grib2_desc = grib2_var(2, 0, 193, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%alhfl_bs,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.TRUE.,                                   &
      & isteptype=a_steptype,                                             &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    ! &      diag%lhfl_pl(nproma,nlev_soil,nblks_c)
    cf_desc    = t_cf_var('lhfl_pl', 'W m-2 ', 'latent heat flux from plants', &
      &          datatype_flt)
    grib2_desc = grib2_var(2, 0, 194, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lhfl_pl', diag%lhfl_pl,                     &
      & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND, cf_desc, grib2_desc, &
      & ldims=(/nproma,nlev_soil,kblks/), lrestart=.FALSE. )
              
    WRITE(name,'(A,A7)') TRIM(prefix),"lhfl_pl"
    WRITE(long_name,'(A27,A4,A18)') "latent heat flux from plants", meaning, &
                                  & " since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), datatype_flt)
    grib2_desc = grib2_var(2, 0, 194, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%alhfl_pl,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND, cf_desc, grib2_desc, &
      & ldims=(/nproma,nlev_soil,kblks/), lrestart=.TRUE.,                &
      & isteptype=a_steptype,                                             &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    ! &      diag%qhfl_s(nproma,nblks_c)
    cf_desc    = t_cf_var('qhfl_s', 'Kg m-2 s-1', 'surface moisture flux', datatype_flt)
    grib2_desc = grib2_var(2, 0, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'qhfl_s', diag%qhfl_s,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.TRUE.,                                   &
      & in_group=groups("pbl_vars"))

    WRITE(name,'(A,A6)') TRIM(prefix),"qhfl_s"
    WRITE(long_name,'(A23,A4,A18)') "surface moisture flux ", meaning, &
                                  & " since model start"
    IF (lflux_avg ) THEN
      cf_desc = t_cf_var(TRIM(name), 'Kg m-2 s-1', TRIM(long_name), datatype_flt)
    ELSE
      cf_desc = t_cf_var(TRIM(name), 'Kg m-2', TRIM(long_name), datatype_flt)
    ENDIF
    grib2_desc = grib2_var(2, 0, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%aqhfl_s  ,                  &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & isteptype=a_steptype, in_group=groups("pbl_vars"),                &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )


    cf_desc    = t_cf_var('qcfl_s', 'kg m-2 s-1',                         &
      &          'surface cloud water deposition flux due to diffusion',  &
      &          datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'qcfl_s', diag%qcfl_s,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d)

    cf_desc    = t_cf_var('qifl_s', 'kg m-2 s-1',                         &
      &          'surface cloud ice deposition flux due to diffusion',    &
      &          datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'qifl_s', diag%qifl_s,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d)


    ! &      diag%tcm(nproma,nblks_c)
    cf_desc    = t_cf_var('tcm', ' ','turbulent transfer coefficients for momentum', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 2, 29, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tcm', diag%tcm,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & in_group=groups("pbl_vars") )

    ! &      diag%tch(nproma,nblks_c)
    cf_desc    = t_cf_var('tch', ' ','turbulent transfer coefficients for heat', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 0, 19, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tch', diag%tch,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & in_group=groups("pbl_vars") )
    
    ! &      diag%tfm(nproma,nblks_c)
    cf_desc    = t_cf_var('tfm', ' ','factor of laminar transfer of momentum', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tfm', diag%tfm,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d )

    ! &      diag%tfh(nproma,nblks_c)
    cf_desc    = t_cf_var('tfh', ' ','factor of laminar transfer of scalars', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tfh', diag%tfh,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d )

    ! &      diag%tfv(nproma,nblks_c)
    cf_desc    = t_cf_var('tfv', ' ','laminar reduction factor for evaporation', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tfv', diag%tfv,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d )

    ! &      diag%tvm(nproma,nblks_c)
    cf_desc    = t_cf_var('tvm', 'm s-1','turbulent transfer velocity for momentum', &
         &                DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tvm', diag%tvm,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & in_group=groups("pbl_vars") )

    ! &      diag%tvh(nproma,nblks_c)
    cf_desc    = t_cf_var('tvh', 'm s-1 ','turbulent transfer velocity for heat', &
         &                DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tvh', diag%tvh,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & in_group=groups("pbl_vars") )

    ! &      diag%tkr(nproma,nblks_c)
    cf_desc    = t_cf_var('tkr', 'm2 s-1 ','turbulent reference surface diffusion coefficient', &
         &                DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tkr', diag%tkr,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.FALSE.,                                  &
      & in_group=groups("pbl_vars") )

    ! &      diag%tkred_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('tkred_sfc', ' ','reduction factor for minimum diffusion coefficients', &
         &                DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tkred_sfc', diag%tkred_sfc,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.FALSE. )

    ! &      diag%pat_len(nproma,nblks_c)
    cf_desc    = t_cf_var('pat_len', 'm','length scale of sub-grid scale roughness elements', &
         &                DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'pat_len', diag%pat_len,                     &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.FALSE. )

    ! &      diag%gz0(nproma,nblks_c)
    cf_desc     = t_cf_var('gz0', 'm2 s-2 ','roughness length times gravity', datatype_flt)
    new_cf_desc = t_cf_var( 'z0',       'm','roughness length',               datatype_flt)
    grib2_desc = grib2_var(2, 0, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'gz0', diag%gz0,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & post_op=post_op(POST_OP_SCALE, arg1=1._wp/grav,                   &
      &                 new_cf=new_cf_desc),                              &
      & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in",               &
      &                 "mode_iau_fg_in","mode_iau_old_fg_in",            &
      &                 "mode_iniana"),                                   &
      & initval=0.01_wp )

    ! &      diag%t_2m(nproma,nblks_c)
    IF (itype_vegetation_cycle == 3) THEN
      in_group = groups("pbl_vars","dwd_fg_atm_vars","mode_iau_ana_in")
    ELSE
      in_group = groups("pbl_vars","dwd_fg_atm_vars")
    ENDIF
    cf_desc    = t_cf_var('t_2m', 'K ','temperature in 2m', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 't_2m', diag%t_2m,                           &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
      & ldims=shape2d, lrestart=.FALSE.,                                  &
      & in_group=in_group,                                                &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF) )


    ! &      diag%t_2m_land(nproma,nblks_c)
    cf_desc    = t_cf_var('t_2m_land', 'K ','temperature in 2m over land fraction', &
      &          datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL) &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 181)
    CALL add_var( diag_list, 't_2m_land', diag%t_2m_land,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M_LAYER, cf_desc, grib2_desc,  &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars") )

    ! &      diag%tmax_2m(nproma,nblks_c)
    cf_desc    = t_cf_var('tmax_2m', 'K ','Max 2m temperature', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tmax_2m', diag%tmax_2m,                     &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
      & ldims=shape2d, lrestart=.TRUE.,                                   &
      & isteptype=TSTEP_MAX, initval=-999._wp, resetval=-999._wp,         &
      & action_list=actions(new_action(ACTION_RESET,"PT06H")),            &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF) ) 

    ! &      diag%tmin_2m(nproma,nblks_c)
    cf_desc    = t_cf_var('tmin_2m', 'K ','Min 2m temperature', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tmin_2m', diag%tmin_2m,                     &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
      & ldims=shape2d, lrestart=.TRUE.,                                   &
      & isteptype=TSTEP_MIN, initval=999._wp, resetval=999._wp,           &
      & action_list=actions(new_action(ACTION_RESET,"PT06H")),            &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF) )

    ! &      diag%qv_2m(nproma,nblks_c)
    cf_desc    = t_cf_var('qv_2m', 'kg kg-1 ','specific water vapor content in 2m', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 1, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'qv_2m', diag%qv_2m,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars") )

    ! &      diag%rh_2m(nproma,nblks_c)
    cf_desc    = t_cf_var('rh_2m', '%','relative humidity in 2m', datatype_flt)
    grib2_desc = grib2_var(0, 1, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rh_2m', diag%rh_2m,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
      & ldims=shape2d, lrestart=.FALSE. )

    ! &      diag%rh_2m_land(nproma,nblks_c)
    cf_desc    = t_cf_var('rh_2m_land', '%','relative humidity in 2m over land fraction', &
      &          datatype_flt)
    grib2_desc = grib2_var(0, 1, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 181)
    CALL add_var( diag_list, 'rh_2m_land', diag%rh_2m_land,               &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M_LAYER, cf_desc, grib2_desc,  &
      & ldims=shape2d, lrestart=.FALSE. )

    ! &      diag%td_2m(nproma,nblks_c)
    cf_desc    = t_cf_var('td_2m', 'K ','dew-point in 2m', datatype_flt)
    grib2_desc = grib2_var(0, 0, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'td_2m', diag%td_2m,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars","dwd_fg_atm_vars") )

    ! &      diag%td_2m_land(nproma,nblks_c)
    cf_desc    = t_cf_var('td_2m_land', 'K ','dew-point in 2m over land fraction', &
      &           datatype_flt)
    grib2_desc = grib2_var(0, 0, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 181)
    CALL add_var( diag_list, 'td_2m_land', diag%td_2m_land,               &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M_LAYER, cf_desc, grib2_desc,  &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars") )

    ! &      diag%u_10m(nproma,nblks_c)
    cf_desc    = t_cf_var('u_10m', 'm s-1 ','zonal wind in 10m', datatype_flt)
    grib2_desc = grib2_var(0, 2, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'u_10m', diag%u_10m,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars","dwd_fg_atm_vars") )

    ! &      diag%v_10m(nproma,nblks_c)
    cf_desc    = t_cf_var('v_10m', 'm s-1 ','meridional wind in 10m', datatype_flt)
    grib2_desc = grib2_var(0, 2, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'v_10m', diag%v_10m,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars","dwd_fg_atm_vars") )

    ! &      diag%sp_10m(nproma,nblks_c)
    cf_desc    = t_cf_var('sp_10m', 'm s-1 ','wind speed in 10m', datatype_flt)
    grib2_desc = grib2_var(0, 2, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sp_10m', diag%sp_10m,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
      & ldims=shape2d, lrestart=.FALSE. )

    !tiled quantities
    ! &      diag%shfl_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('shfl_s_t', 'W m-2 ', 'tile-based surface sensible heat flux', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 0, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'shfl_s_t', diag%shfl_s_t,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! fill the separate variables belonging to the container shfl_s_t
    ALLOCATE(diag%shfl_s_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'shfl_s_t',                            &
         & 'shfl_s_t_'//TRIM(ADJUSTL(csfc)),                          &
         & diag%shfl_s_t_ptr(jsfc)%p_2d,                              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('shfl_s_t_'//TRIM(csfc), '', '', datatype_flt), &
         & grib2_var(0, 0, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
         & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
    ENDDO

    ! &      diag%lhfl_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('lhfl_s_t', 'W m-2 ', 'tile-based surface latent heat flux', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 0, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lhfl_s_t', diag%lhfl_s_t,                                &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,   &
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! fill the separate variables belonging to the container lhfl_s_t
    ALLOCATE(diag%lhfl_s_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'lhfl_s_t',                            &
         & 'lhfl_s_t_'//TRIM(ADJUSTL(csfc)),                          &
         & diag%lhfl_s_t_ptr(jsfc)%p_2d,                              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        & 
         & t_cf_var('lhfl_s_t_'//TRIM(csfc), '', '', datatype_flt), &
         & grib2_var(0, 0, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL),  & 
         & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
    ENDDO

    ! &      diag%lhfl_bs_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('lhfl_bs_t', 'W m-2 ', 'tile-based latent heat flux from bare soil', &
      &                   datatype_flt)
    grib2_desc = grib2_var(2, 0, 193, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lhfl_bs_t', diag%lhfl_bs_t,                              &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs,    &
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! fill the separate variables belonging to the container lhfl_bs_t
    ALLOCATE(diag%lhfl_bs_t_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'lhfl_bs_t',                           &
         & 'lhfl_bs_t_'//TRIM(ADJUSTL(csfc)),                         &
         & diag%lhfl_bs_t_ptr(jsfc)%p_2d,                             &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        & 
         & t_cf_var('lhfl_bs_t_'//TRIM(csfc), '', '', datatype_flt),&
         & grib2_var(2, 0, 193, ibits, GRID_UNSTRUCTURED, GRID_CELL), & 
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.)
    ENDDO

    ! &      diag%lhfl_pl_t(nproma,nlev_soil,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('lhfl_pl_t', 'W m-2 ', 'tile-based latent heat flux from plants', &
      &                   datatype_flt)
    grib2_desc = grib2_var(2, 0, 194, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lhfl_pl_t', diag%lhfl_pl_t,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND, cf_desc, grib2_desc, &
      & ldims=(/nproma,nlev_soil,kblks,ntiles_total/), lcontainer=.TRUE., &
      & lrestart=.FALSE., loutput=.FALSE.)

    ! fill the separate variables belonging to the container lhfl_pl_t
    ALLOCATE(diag%lhfl_pl_t_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'lhfl_pl_t',                           &
         & 'lhfl_pl_t_'//TRIM(ADJUSTL(csfc)),                         &
         & diag%lhfl_pl_t_ptr(jsfc)%p_3d,                             &
         & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,               & 
         & t_cf_var('lhfl_pl_t_'//TRIM(csfc), '', '', datatype_flt),&
         & grib2_var(2, 0, 194, ibits, GRID_UNSTRUCTURED, GRID_CELL), & 
         & ldims=(/nproma,nlev_soil,kblks/), lrestart=.FALSE.,        &
         & loutput=.TRUE.)
    ENDDO

    ! &      diag%qhfl_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('qhfl_s_t', 'Kg m-2 s-1','tile based surface moisture flux', &
         &                datatype_flt)
    grib2_desc = grib2_var(2, 0, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'qhfl_s_t', diag%qhfl_s_t,                                &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,   &
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! fill the separate variables belonging to the container qhfl_s_t
    ALLOCATE(diag%qhfl_s_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'qhfl_s_t',                            &
         & 'qhfl_s_t_'//TRIM(ADJUSTL(csfc)),                          &
         & diag%qhfl_s_t_ptr(jsfc)%p_2d,                              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        & 
         & t_cf_var('qhfl_s_t_'//TRIM(csfc), '', '', datatype_flt),   &
         & grib2_var(2, 0, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL),   & 
         & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
    ENDDO

    ! &      diag%tcm_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tcm_t', ' ', &
      & 'tile-based turbulent transfer coefficients for momentum', datatype_flt)
    grib2_desc = grib2_var(0, 2, 29, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tcm_t', diag%tcm_t,                                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! fill the separate variables belonging to the container tcm_t
    ALLOCATE(diag%tcm_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tcm_t',                               &
         & 'tcm_t_'//TRIM(ADJUSTL(csfc)),                             &
         & diag%tcm_t_ptr(jsfc)%p_2d,                                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        & 
         & t_cf_var('tcm_t_'//TRIM(csfc), '', '', datatype_flt),    &
         & grib2_var(0, 2, 29, ibits, GRID_UNSTRUCTURED, GRID_CELL),  & 
         & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
    ENDDO

    ! &      diag%tch_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tch_t', ' ', &
         &                'tile-based turbulent transfer coefficients for heat', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 0, 19, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tch_t', diag%tch_t,                                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)                       

    ! fill the separate variables belonging to the container tch_t
    ALLOCATE(diag%tch_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tch_t',                               &
         & 'tch_t_'//TRIM(ADJUSTL(csfc)),                             &
         & diag%tch_t_ptr(jsfc)%p_2d,                                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('tch_t_'//TRIM(csfc), '', '', datatype_flt),    &
         & grib2_var(0, 0, 19, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
         & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
    ENDDO

    ! &      diag%tfv_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tfv_t', ' ', &
         &                'tile-based laminar reduction factor for evaporation', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tfv_t', diag%tfv_t,                                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! fill the separate variables belonging to the container tfv_t
    ALLOCATE(diag%tfv_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tfv_t',                               &
         & 'tfv_t_'//TRIM(ADJUSTL(csfc)),                             &
         & diag%tfv_t_ptr(jsfc)%p_2d,                                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('tfv_t_'//TRIM(csfc), '', '', datatype_flt),    &
         & grib2_var(0, 4, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
         & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
    ENDDO

    ! &      diag%tvm_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tvm_t', 'm s-1', &
      & 'tile-based turbulent transfer velocity for momentum', DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tvm_t', diag%tvm_t,                                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! fill the separate variables belonging to the container tvm_t
    ALLOCATE(diag%tvm_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tvm_t',                               &
         & 'tvm_t_'//TRIM(ADJUSTL(csfc)),                             &
         & diag%tvm_t_ptr(jsfc)%p_2d,                                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('tvm_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),    &
         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
!        & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)  !nec. for rest. only if 'tvm' substitutes 'tcm'
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.)
    ENDDO

    ! &      diag%tvh_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tvh_t', 'm s-1', &
         &                'tile-based turbulent transfer velocity for heat', &
         &                DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tvh_t', diag%tvh_t,                                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! fill the separate variables belonging to the container tvh_t
    ALLOCATE(diag%tvh_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tvh_t',                               &
         & 'tvh_t_'//TRIM(ADJUSTL(csfc)),                             &
         & diag%tvh_t_ptr(jsfc)%p_2d,                                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('tvh_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),    &
         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
!        & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.) !nec. for rest. only if 'tvm' substitutes 'tcm'
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.)
    ENDDO

    ! &      diag%tkr_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tkr_t', 'm2 s-1', &
      & 'tile-based turbulent reference surface diffusion coefficient', DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tkr_t', diag%tkr_t,                                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! fill the separate variables belonging to the container tkr_t
    ALLOCATE(diag%tkr_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tkr_t',                               &
         & 'tkr_t_'//TRIM(ADJUSTL(csfc)),                             &
         & diag%tkr_t_ptr(jsfc)%p_2d,                                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('tkr_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),    &
         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
!        & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.) !nec. for rest. only if 'imode_trancnf>=4'
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.)
    ENDDO


    ! &      diag%gz0_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('gz0_t', 'm2 s-2 ', 'tile-based roughness length times gravity', &
         &                datatype_flt)
    grib2_desc = grib2_var(2, 0, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'gz0_t', diag%gz0_t,                                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw, &
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)



    ! fill the separate variables belonging to the container gz0_t
    ALLOCATE(diag%gz0_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'gz0_t',                               &
         & 'gz0_t_'//TRIM(ADJUSTL(csfc)),                             &
         & diag%gz0_t_ptr(jsfc)%p_2d,                                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('gz0_t_'//TRIM(csfc), '', '', datatype_flt),    &
         & grib2_var(2, 0, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
         & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
    ENDDO


    ! &      diag%tvs_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tvs_s_t', 'm2 s-2 ',                              &
         &                'tile-based turbulence velocity scale at surface', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tvs_s_t', diag%tvs_s_t,                                &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw, &
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)


    ! fill the separate variables belonging to the container tvs_s_t
    ALLOCATE(diag%tvs_s_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tvs_s_t',                                &
         & 'tvs_s_t_'//TRIM(ADJUSTL(csfc)),                              &
         & diag%tvs_s_t_ptr(jsfc)%p_2d,                                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
         & t_cf_var('tvs_s_t_'//TRIM(csfc), '', '', datatype_flt),     &
         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
         & ldims=shape2d, lrestart=.TRUE., loutput=.FALSE.)
    ENDDO



    ! &      diag%tkvm_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tkvm_s_t', 'm s-2 ',                                       &
         &                'tile-based turbulent diff. coeff for momentum at surface', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 2, 31, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tkvm_s_t', diag%tkvm_s_t,                              &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw, &
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)


    ! fill the separate variables belonging to the container tkvm_s_t
    ALLOCATE(diag%tkvm_s_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tkvm_s_t',                                &
         & 'tkvm_s_t_'//TRIM(ADJUSTL(csfc)),                              &
         & diag%tkvm_s_t_ptr(jsfc)%p_2d,                                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
         & t_cf_var('tkvm_s_t_'//TRIM(csfc), '', '', datatype_flt),     &
         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
         & ldims=shape2d, lrestart=.TRUE., loutput=.FALSE.)
    ENDDO


    ! &      diag%tkvh_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tkvh_s_t', 'm s-2 ',                                      &
         &                'tile-based turbulent diff. coeff for heat at surface',    &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 0, 20, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tkvh_s_t', diag%tkvh_s_t,                              &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw, &
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)


    ! fill the separate variables belonging to the container tkvm_s_t
    ALLOCATE(diag%tkvh_s_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tkvh_s_t',                                &
         & 'tkvh_s_t_'//TRIM(ADJUSTL(csfc)),                              &
         & diag%tkvh_s_t_ptr(jsfc)%p_2d,                                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
         & t_cf_var('tkvh_s_t_'//TRIM(csfc), '', '', datatype_flt),     &
         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
         & ldims=shape2d, lrestart=.TRUE., loutput=.FALSE.)
    ENDDO



    ! &      diag%u_10m_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('u_10m_t', 'm s-1 ', 'tile-based zonal wind in 2m', datatype_flt)
    grib2_desc = grib2_var(0, 2, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'u_10m_t', diag%u_10m_t,                                  &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
  
    ! fill the separate variables belonging to the container u_10m_t
    ALLOCATE(diag%u_10m_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'u_10m_t',                             &
         & 'u_10m_t_'//TRIM(ADJUSTL(csfc)),                           &
         & diag%u_10m_t_ptr(jsfc)%p_2d,                               &
         & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M,                     &
         & t_cf_var('u_10m_t_'//TRIM(csfc), '', '', datatype_flt),  &
         & grib2_var(0, 2, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
         & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
    ENDDO


    ! &      diag%v_10m_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('v_10m_t', 'm s-1 ', 'tile-based meridional wind in 2m', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 2, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'v_10m_t', diag%v_10m_t,                                  &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)


    ! fill the separate variables belonging to the container v_10m_t
    ALLOCATE(diag%v_10m_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'v_10m_t',                             &
         & 'v_10m_t_'//TRIM(ADJUSTL(csfc)),                           &
         & diag%v_10m_t_ptr(jsfc)%p_2d,                               &
         & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M,                     &
         & t_cf_var('v_10m_t_'//TRIM(csfc), '', '', datatype_flt),  &
         & grib2_var(0, 2, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
         & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
    ENDDO



    ! &      diag%umfl_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('umfl_s_t', 'N m-2 ', 'u-momentum flux at the surface', datatype_flt)
    grib2_desc = grib2_var(0, 2, 17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'umfl_s_t', diag%umfl_s_t,                                  &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
  
    ! fill the separate variables belonging to the container umfl_s_t
    ALLOCATE(diag%umfl_s_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'umfl_s_t',                            &
         & 'umfl_s_t_'//TRIM(ADJUSTL(csfc)),                          &
         & diag%umfl_s_t_ptr(jsfc)%p_2d,                              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('umfl_s_t_'//TRIM(csfc), '', '', datatype_flt), &
         & grib2_var(0, 2, 17, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
         & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,            &
         & isteptype=TSTEP_INSTANT )
    ENDDO
    !EDMF requires lrestart=.TRUE. 


    ! &      diag%vmfl_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('vmfl_s_t', 'N m-2 ', 'v-momentum flux at the surface', datatype_flt)
    grib2_desc = grib2_var(0, 2, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'vmfl_s_t', diag%vmfl_s_t,                                  &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
  
    ! fill the separate variables belonging to the container vmfl_s_t
    ALLOCATE(diag%vmfl_s_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'vmfl_s_t',                            &
         & 'vmfl_s_t_'//TRIM(ADJUSTL(csfc)),                          &
         & diag%vmfl_s_t_ptr(jsfc)%p_2d,                              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('vmfl_s_t_'//TRIM(csfc), '', '', datatype_flt), &
         & grib2_var(0, 2, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
         & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,            &
         & isteptype=TSTEP_INSTANT )
    ENDDO
    !EDMF requires lrestart=.TRUE. 


    ! &      diag%umfl_s(nproma,nblks_c)
    cf_desc    = t_cf_var('umfl_s', 'N m-2', 'u-momentum flux at the surface', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 2, 17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'umfl_s', diag%umfl_s,                            &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
      & lrestart=.FALSE., loutput=.TRUE.,                                      &
      & isteptype=TSTEP_INSTANT )

    ! &      diag%vmfl_s(nproma,nblks_c)
    cf_desc    = t_cf_var('vmfl_s', 'N m-2', 'v-momentum flux at the surface', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 2, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'vmfl_s', diag%vmfl_s,                            &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
      & lrestart=.FALSE., loutput=.TRUE.,                                      &
      & isteptype=TSTEP_INSTANT )

    IF (lflux_avg ) THEN
        prefix = "a"
        meaning = "mean"
        varunits= "N/m**2"
        a_steptype= TSTEP_AVG 
    ELSE
        prefix = "acc"
        meaning = "acc." 
        varunits= "Ns/m**2"    ! or "kg/(m*s)"
        a_steptype= TSTEP_ACCUM     
    END IF

    WRITE(name,'(A,A6)') TRIM(prefix),"umfl_s"
    WRITE(long_name,'(A26,A4,A18)') "u-momentum flux flux at surface ", meaning, &
                                  & " since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), &
      &                          TRIM(long_name), datatype_flt)
    grib2_desc = grib2_var(0, 2, 17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
!   aumfl_s and avmfl_s are needed for the restart only to get the same output values
!   They are not important to obtain bit identical model results
    CALL add_var( diag_list, TRIM(name), diag%aumfl_s,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
      & isteptype=a_steptype, lrestart=.TRUE., loutput=.TRUE.,                 &
      & hor_interp=create_hor_interp_metadata(                                 &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                               &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    WRITE(name,'(A,A6)') TRIM(prefix),"vmfl_s"
    WRITE(long_name,'(A26,A4,A18)') "v-momentum flux flux at surface ", meaning, &
                                  & " since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), &
      &                          TRIM(long_name), datatype_flt)
    grib2_desc = grib2_var(0, 2, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%avmfl_s,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
      & isteptype=a_steptype, lrestart=.TRUE., loutput=.TRUE.,                 &
      & hor_interp=create_hor_interp_metadata(                                 &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                               &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )


    !------------------------------
    ! SSO surface stress
    !------------------------------

    cf_desc    = t_cf_var('str_u_sso', 'N m-2', 'zonal sso surface stress', datatype_flt)
    grib2_desc = grib2_var(192, 128, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'str_u_sso', diag%str_u_sso,                     &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('str_v_sso', 'N m-2', 'meridional sso surface stress', datatype_flt)
    grib2_desc = grib2_var(192, 128, 196, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'str_v_sso', diag%str_v_sso,                     &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    WRITE(name,'(A,A9)') TRIM(prefix),"str_u_sso"
    WRITE(long_name,'(A25,A4,A18)') "zonal sso surface stress ", meaning, " since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), datatype_flt)
    grib2_desc = grib2_var(192, 128, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%astr_u_sso,                      &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
      & isteptype=a_steptype, lrestart=.TRUE., loutput=.TRUE.,                 &
      & hor_interp=create_hor_interp_metadata(                                 &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                               &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    WRITE(name,'(A,A9)') TRIM(prefix),"str_v_sso"
    WRITE(long_name,'(A30,A4,A18)') "meridional sso surface stress ", meaning, " since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), datatype_flt)
    grib2_desc = grib2_var(192, 128, 196, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%astr_v_sso,                      &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
      & isteptype=a_steptype, lrestart=.TRUE., loutput=.TRUE.,                 &
      & hor_interp=create_hor_interp_metadata(                                 &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                               &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )


    !------------------------------
    ! grid-scale surface stress
    !------------------------------

    cf_desc    = t_cf_var('drag_u_grid', 'N m-2', 'zonal resolved surface stress', datatype_flt)
    grib2_desc = grib2_var(192, 150, 168, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'drag_u_grid', diag%drag_u_grid,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    cf_desc    = t_cf_var('drag_v_grid', 'N m-2', 'meridional resolved surface stress', datatype_flt)
    grib2_desc = grib2_var(192, 150, 169, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'drag_v_grid', diag%drag_v_grid,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    WRITE(name,'(A,A11)') TRIM(prefix),"drag_u_grid"
    WRITE(long_name,'(A30,A4,A18)') "zonal resolved surface stress ", meaning, " since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), datatype_flt)
    grib2_desc = grib2_var(192, 150, 168, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%adrag_u_grid,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
      & isteptype=a_steptype, lrestart=.TRUE., loutput=.TRUE.,                 &
      & hor_interp=create_hor_interp_metadata(                                 &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                               &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )

    WRITE(name,'(A,A11)') TRIM(prefix),"drag_v_grid"
    WRITE(long_name,'(A35,A4,A18)') "meridional resolved surface stress ", meaning, " since model start"
    cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), datatype_flt)
    grib2_desc = grib2_var(192, 150, 169, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%adrag_v_grid,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
      & isteptype=a_steptype, lrestart=.TRUE., loutput=.TRUE.,                 &
      & hor_interp=create_hor_interp_metadata(                                 &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                               &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ) )


    !------------------------------
    ! EDMF
    !------------------------------

    IF( atm_phy_nwp_config(k_jg)%inwp_turb == iedmf) THEN

       ! &      diag%z0m(nproma,nblks_c)
       cf_desc    = t_cf_var('z0m', '', &
            &                'geopotential of the top of the atmospheric boundary layer', &
            &                datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( diag_list, 'z0m', diag%z0m,                             &
         & GRID_UNSTRUCTURED_CELL,ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ENDIF  !inwp_turb == EDMF

    !------------------------------
    ! LES
    !------------------------------

    !Anurag Dipankar, MPI (7 Oct 2013)
    !Diagnostics for LES physics
    IF ( atm_phy_nwp_config(k_jg)%is_les_phy ) THEN

      ! &      diag%z_pbl(nproma,nblks_c)
      cf_desc    = t_cf_var('z_pbl', 'm', 'boundary layer height above sea level', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'z_pbl', diag%z_pbl,                             &
        & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,            &
        & ldims=shape2d, lrestart=.FALSE. )

      ! &      diag%bruvais(nproma,nlev+1,nblks_c)
      cf_desc    = t_cf_var('bruvais', '1/s**2', 'Brunt Vaisala Frequency', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'bruvais', diag%bruvais,                     &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,      &
        & ldims=shape3dkp1, initval=0.0_wp, lrestart=.FALSE. )                                   

      ! &      diag%mech_prod(nproma,nlev+1,nblks_c)
      cf_desc    = t_cf_var('mech_prod', 'm**2/s**3', 'mechanical production term in TKE Eq', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'mech_prod', diag%mech_prod,                     &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,          &
        & ldims=shape3dkp1, initval=0.0_wp, lrestart=.FALSE. )                                   

      !1D and 0D diagnostic variables that can not be part of add_var
      ALLOCATE( diag%turb_diag_1dvar(klevp1,SIZE(turb_profile_list,1)),  &
                diag%turb_diag_0dvar(SIZE(turb_tseries_list,1)), STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:new_nwp_phy_diag_list', &
                    'allocation of 1D and 0D diag var failed!')
      ENDIF
      diag%turb_diag_1dvar = 0._wp
      diag%turb_diag_0dvar = 0._wp

      !  
      !Some diagnostics specific to HDCP2
      !

      ! &      diag%t_cbase(nproma,nblks_c) 
      cf_desc    = t_cf_var('t_cbase', 'K', 'cloud base temperature', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 't_cbase', diag%t_cbase,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,            &
        & ldims=shape2d, lrestart=.FALSE. )
     
      ! &      diag%p_cbase(nproma,nblks_c) 
      cf_desc    = t_cf_var('p_cbase', 'Pa', 'cloud base pressure', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'p_cbase', diag%p_cbase,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,            &
        & ldims=shape2d, lrestart=.FALSE. )

      ! &      diag%t_ctop(nproma,nblks_c) 
      cf_desc    = t_cf_var('t_ctop', 'K', 'cloud top temperature', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 't_ctop', diag%t_ctop,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,          &
        & ldims=shape2d, lrestart=.FALSE. )

      ! &      diag%p_ctop(nproma,nblks_c) 
      cf_desc    = t_cf_var('p_ctop', 'K', 'cloud top pressure', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'p_ctop', diag%p_ctop,                           &
        & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,            &
        & ldims=shape2d, lrestart=.FALSE. )

    END IF  


    !------------------
    !Turbulence 3D variables
    !------------------

   IF (.NOT.atm_phy_nwp_config(k_jg)%is_les_phy) THEN
   ! &      diag%tkvm(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('tkvm', 'm**2/s', ' turbulent diffusion coefficients for momentum', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 2, 31, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tkvm', diag%tkvm,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,        &
      & ldims=shape3dkp1, in_group=groups("pbl_vars") )

   ! &      diag%tkvh(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('tkvh', 'm**2/s', ' turbulent diffusion coefficients for heat', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 0, 20, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tkvh', diag%tkvh,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,        &
      & ldims=shape3dkp1, in_group=groups("pbl_vars") ) 
   ELSE
     ! &      diag%tkvm(nproma,nlevp1,nblks_c)
      cf_desc    = t_cf_var('tkvm', 'kg/(ms)', ' mass weighted turbulent viscosity', &
           &                datatype_flt)
      grib2_desc = grib2_var(0, 2, 31, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'tkvm', diag%tkvm,                             &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,        &
        & ldims=shape3dkp1, lrestart=.FALSE., in_group=groups("pbl_vars") )
  
     ! &      diag%tkvh(nproma,nlevp1,nblks_c)
      cf_desc    = t_cf_var('tkvh', 'kg/(ms)', ' mass weighted turbulent diffusivity', &
           &                datatype_flt)
      grib2_desc = grib2_var(0, 0, 20, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'tkvh', diag%tkvh,                             &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,        &
        & ldims=shape3dkp1, lrestart=.FALSE., in_group=groups("pbl_vars") ) 
    END IF

   ! &      diag%rcld(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('rcld', '', 'standard deviation of the saturation deficit', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rcld', diag%rcld,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,        &
      & ldims=shape3dkp1 )

! SO FAR UNUSED
!!$   ! &      diag%edr(nproma,nlevp1,nblks_c)
!!$    cf_desc    = t_cf_var('edr', '', 'eddy dissipation rate', datatype_flt)
!!$    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
!!$    CALL add_var( diag_list, 'edr', diag%edr,                               &
!!$      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,        &
!!$      & ldims=shape3dkp1, lrestart=.FALSE. ) 


    !------------------
    ! Optional computation of diagnostic fields
    !------------------

    ! Note: These tasks are registered for the post-processing scheduler
    !       which takes care of the regular update:
    ! 
    ! &     relative humidity
    !    
    IF (l_rh) THEN
      cf_desc    = t_cf_var('rh', '%', 'relative humidity', datatype_flt)
      grib2_desc = grib2_var(0, 1, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                       &
                    & "rh", diag%rh,                                                 &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape3d,                                                 &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=vintp_types("P","Z","I"),           &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.FALSE., l_pd_limit=.TRUE.,             &
                    &             lower_limit=0._wp ),                               &
                    & l_pp_scheduler_task=TASK_COMPUTE_RH, lrestart=.FALSE.          )
    END IF
    
    ! &     potential vorticity
    
    IF (l_pv) THEN
      cf_desc    = t_cf_var('pv', 'K m2 kg-1 s-1', 'potential vorticity', DATATYPE_FLT32)
      grib2_desc = grib2_var(0, 2, 14, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                       &
                    & "pv", diag%pv,                                                 &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape3d,                                                 &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=vintp_types("P","Z","I"),           &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.FALSE.),                               &
                    & l_pp_scheduler_task=TASK_COMPUTE_PV, lrestart=.FALSE.          )
    END IF


    !  Height of 0 deg C level
    cf_desc    = t_cf_var('hzerocl', '', 'height of 0 deg C level', datatype_flt)
    grib2_desc = grib2_var(0, 3, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 101)
    CALL add_var( diag_list, 'hzerocl', diag%hzerocl,                                &
      &           GRID_UNSTRUCTURED_CELL, ZA_ISOTHERM_ZERO, cf_desc, grib2_desc,     &
      &           ldims=shape2d, lrestart=.FALSE.,                                   &
!!$      &           lmiss=.TRUE., missval=-999._wp,                                    &
      &           hor_interp=create_hor_interp_metadata(                             &
      &                      hor_intp_type=HINTP_TYPE_LONLAT_NNB) )


    !  significant weather WW
    cf_desc    = t_cf_var('ww', '', 'significant weather', datatype_flt)
    grib2_desc = grib2_var(0, 19, 25, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'ww', diag%iww,                                         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,           &
                & ldims=shape2d, lrestart=.FALSE.,                                   &
                & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB) )
      
    ! mask field to distinguish between tropics and extratropics (for tuning purposes)
    cf_desc    = t_cf_var('tropics_mask', '', 'tropics_mask', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tropics_mask', diag%tropics_mask,                      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,           &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    ! mask field to distinguish between inner tropics and elsewhere (for tuning purposes)
    cf_desc    = t_cf_var('innertropics_mask', '', 'innertropics_mask', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'innertropics_mask', diag%innertropics_mask,            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,           &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    ! buffer field needed for the combination of vertical nesting with a reduced radiation grid
    IF (k_jg > n_dom_start) THEN
      ALLOCATE(diag%buffer_rrg(nproma, 3*nexlevs_rrg_vnest, p_patch_local_parent(k_jg)%nblks_c))
    ENDIF

    ! buffer field needed for the combination of vertical nesting with a reduced radiation grid
    diag%buffer_rttov => NULL()
    IF (lsynsat(k_jg)) THEN
      IF  ((k_jg > n_dom_start) .AND. (p_patch(k_jg)%nshift > 0)) THEN
        ALLOCATE(diag%buffer_rttov(nproma, 5*p_patch(k_jg)%nshift, p_patch_local_parent(k_jg)%nblks_c))
      ENDIF
      
      shape3d_synsat = (/nproma, num_images, p_patch(k_jg)%nblks_c /)
      shape2d_synsat = (/nproma,             p_patch(k_jg)%nblks_c /)

      ! introduce container variable for RTTOV synthetic satellite imagery:
      cf_desc    = t_cf_var('rttov_channels', '', '', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'rttov_channels', diag%synsat_arr,                 &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
        &           ldims=shape3d_synsat ,                                        &
        &           lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,         &
        &           var_class=CLASS_SYNSAT)

      ! add reference variables for the different images:
      ALLOCATE(diag%synsat_image(num_images))

      iimage = 0
      sensor_loop: DO isens = 1, num_sensors

        DO ichan = 1,total_numchans(isens)

          DO k=1,4
            ! the following translation can be derived by gazing at the corresponding RTTOV loop
            lradiance = ((MOD(iimage,4)+1) == RTTOV_RAD_CL) .OR. ((MOD(iimage,4)+1) == RTTOV_RAD_CS)
            lcloudy   = ((MOD(iimage,4)+1) == RTTOV_BT_CL)  .OR. ((MOD(iimage,4)+1) == RTTOV_RAD_CL)
            iimage = iimage + 1
          
            IF (lradiance) THEN
              unit = "mW/cm-1/sr/sq.m"
            ELSE
              unit = "K"
            END IF

            CALL get_synsat_name(lradiance, lcloudy, ichan, shortname, longname)
            CALL get_synsat_grib_triple(lradiance, lcloudy, ichan,       &
              &                         idiscipline, icategory, inumber, &
              &                         wave_no, wave_no_scalfac)
            
            cf_desc    = t_cf_var(TRIM(shortname), TRIM(unit), TRIM(longname), datatype_flt)
            grib2_desc = grib2_var(idiscipline, icategory, inumber, ibits, GRID_UNSTRUCTURED, GRID_CELL)   &
              &           + t_grib2_int_key("scaledValueOfCentralWaveNumber", wave_no)                  &
              &           + t_grib2_int_key("scaleFactorOfCentralWaveNumber", wave_no_scalfac)          &
              &           + t_grib2_int_key("satelliteSeries", 333)                                     &
              &           + t_grib2_int_key("satelliteNumber",  72)                                     &
              &           + t_grib2_int_key("instrumentType",  207)
            CALL add_ref( diag_list, 'rttov_channels', TRIM(shortname),                  &
              &           diag%synsat_image(iimage)%p,                                   &
              &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
              &           cf_desc, grib2_desc, ldims=shape2d_synsat,                     &
              &           opt_var_ref_pos = 2, lrestart=.FALSE., loutput=.TRUE.,         &
              &           in_group=groups("RTTOV"), var_class=CLASS_SYNSAT )
          END DO
        END DO
      END DO sensor_loop

    ENDIF

    CALL message('mo_nwp_phy_state:construct_nwp_phy_diag', &
                 'construction of NWP physical fields finished')  

END SUBROUTINE new_nwp_phy_diag_list



SUBROUTINE new_nwp_phy_tend_list( k_jg, klev,  kblks,   &
                     & listname, phy_tend_list, phy_tend)

    INTEGER,INTENT(IN) :: k_jg, klev, kblks !< dimension sizes

    CHARACTER(len=*),INTENT(IN) :: listname

    TYPE(t_var_list)    ,INTENT(INOUT) :: phy_tend_list
    TYPE(t_nwp_phy_tend),INTENT(INOUT) :: phy_tend

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape3d(3), shape3dkp1(3), shape4d(4), shape4d_conv(4), shape4d_gscp(4)
    INTEGER :: ibits, ktracer, ist, ntr_conv
    LOGICAL :: lrestart
    INTEGER :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ibits = DATATYPE_PACK16 ! "entropy" of horizontal slice

    shape3d    = (/nproma, klev  , kblks            /)
    shape3dkp1 = (/nproma, klev+1, kblks            /)

    IF (lart) THEN
     shape4d    = (/nproma, klev  , kblks, nqtendphy+nart_tendphy /)
    ELSE
     shape4d    = (/nproma, klev  , kblks, nqtendphy /)
    ENDIF 
      
    ! dimension of convective tracer field
    ntr_conv = nqtendphy
    IF (lart)                                        ntr_conv = ntr_conv + nart_tendphy
    IF (atm_phy_nwp_config(k_jg)%ldetrain_conv_prec) ntr_conv = ntr_conv + 2 ! plus qr and qs

    shape4d_conv = (/nproma, klev  , kblks, ntr_conv /)
    shape4d_gscp = (/nproma, klev  , kblks, 5 /)

    NULLIFY(phy_tend%ddt_temp_gscp, phy_tend%ddt_tracer_gscp)

    CALL new_var_list( phy_tend_list, TRIM(listname), patch_id=k_jg )
    CALL default_var_list_settings( phy_tend_list,             &
                                  & lrestart=.TRUE.  )

    
    !------------------------------
    ! Temperature tendencies
    !------------------------------

   ! &      phy_tend%ddt_temp_radsw(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_temp_radsw', 'K s-1', &
         &                            'short wave radiative temperature tendency', datatype_flt)
    grib2_desc = grib2_var(0, 4, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_radsw', phy_tend%ddt_temp_radsw,      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,        &
                & ldims=shape3d, in_group=groups("phys_tendencies") )

   ! &      phy_tend%ddt_temp_radlw(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_temp_radlw', 'K s-1', &
         &                            'long wave radiative temperature tendency', datatype_flt)
    grib2_desc = grib2_var(0, 5, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_radlw', phy_tend%ddt_temp_radlw,      &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,        &
                & ldims=shape3d, in_group=groups("phys_tendencies") )

   ! &      phy_tend%ddt_temp_turb(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_temp_turb', 'K s-1', &
         &                            'turbulence temperature tendency', datatype_flt)
    grib2_desc = grib2_var(192, 162, 121, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_turb', phy_tend%ddt_temp_turb,        &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,        &
                & ldims=shape3d, lrestart=.FALSE., in_group=groups("phys_tendencies"))

    IF ( .NOT. atm_phy_nwp_config(k_jg)%is_les_phy ) THEN
     ! &      phy_tend%ddt_temp_drag(nproma,nlev,nblks)
      cf_desc    = t_cf_var('ddt_temp_drag', 'K s-1', &
           &                'sso + gwdrag temperature tendency', datatype_flt)
      grib2_desc = grib2_var(192, 162, 125, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( phy_tend_list, 'ddt_temp_drag', phy_tend%ddt_temp_drag,        &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,        &
                  & ldims=shape3d, lrestart=.FALSE., in_group=groups("phys_tendencies"))

     ! &      phy_tend%ddt_temp_pconv(nproma,nlev,nblks)
      cf_desc    = t_cf_var('ddt_temp_pconv', 'K s-1', &
           &                            'convective temperature tendency', datatype_flt)
      grib2_desc = grib2_var(0, 0, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( phy_tend_list, 'ddt_temp_pconv', phy_tend%ddt_temp_pconv,        &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d )
    END IF

    IF (atm_phy_nwp_config(k_jg)%is_les_phy .OR. &
        is_variable_in_output(first_output_name_list, var_name="ddt_temp_gscp")) THEN
      ! &      phy_tend%ddt_temp_gscp(nproma,nlev,nblks)
      cf_desc    = t_cf_var('ddt_temp_gscp', 'K s-1', &
           &                            'microphysical temperature tendency', datatype_flt)
      grib2_desc = grib2_var(192, 162, 203, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( phy_tend_list, 'ddt_temp_gscp', phy_tend%ddt_temp_gscp,        &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,        &
                  & ldims=shape3d, lrestart=.FALSE., in_group=groups("phys_tendencies") )
    END IF

    !------------------------------
    ! Zonal Wind tendencies
    !------------------------------

   ! &      phy_tend%ddt_u_turb(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_u_turb', 'm s-2', &
         &                            'turbulence tendency of zonal wind', datatype_flt)
    grib2_desc = grib2_var(192, 162, 119, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_turb', phy_tend%ddt_u_turb,          &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d, lrestart=.FALSE., in_group=groups("phys_tendencies"))

    IF ( .NOT. atm_phy_nwp_config(k_jg)%is_les_phy ) THEN
      ! &      phy_tend%ddt_u_sso(nproma,nlev,nblks)
       cf_desc    = t_cf_var('ddt_u_sso', 'm s-2', &
            &                            'sso tendency of zonal wind', datatype_flt)
       grib2_desc = grib2_var(0, 2, 194, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( phy_tend_list, 'ddt_u_sso', phy_tend%ddt_u_sso,            &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                   & ldims=shape3d, in_group=groups("phys_tendencies") )
   
      ! &      phy_tend%ddt_u_gwd(nproma,nlev,nblks)
       cf_desc    = t_cf_var('ddt_u_gwd', 'm s-2', &
            &                            'GWD tendency of zonal wind', datatype_flt)
       grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( phy_tend_list, 'ddt_u_gwd', phy_tend%ddt_u_gwd,            &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                   & ldims=shape3d, in_group=groups("phys_tendencies") )
   
      ! &      phy_tend%ddt_u_pconv(nproma,nlev,nblks)
       cf_desc    = t_cf_var('ddt_u_pconv', 'm s-2', &
            &                            'convective tendency of zonal wind', datatype_flt)
       grib2_desc = grib2_var(0, 2, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( phy_tend_list, 'ddt_u_pconv', phy_tend%ddt_u_pconv,        &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d)
    END IF


    !------------------------------
    ! Meridional Wind tendencies
    !------------------------------

    ! &      phy_tend%ddt_v_turb(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_v_turb', 'm s-2', &
         &                            'turbulence tendency of meridional wind', datatype_flt)
    grib2_desc = grib2_var(192, 162, 120, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_turb', phy_tend%ddt_v_turb,          &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d, lrestart=.FALSE., in_group=groups("phys_tendencies"))

    IF ( .NOT. atm_phy_nwp_config(k_jg)%is_les_phy ) THEN
      ! &      phy_tend%ddt_v_sso(nproma,nlev,nblks)
      cf_desc    = t_cf_var('ddt_v_sso', 'm s-2', &
           &                            'sso tendency of meridional wind', datatype_flt)
      grib2_desc = grib2_var(0, 2, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( phy_tend_list, 'ddt_v_sso', phy_tend%ddt_v_sso,            &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                  & ldims=shape3d, in_group=groups("phys_tendencies") )
  
      ! &      phy_tend%ddt_v_gwd(nproma,nlev,nblks)
      cf_desc    = t_cf_var('ddt_v_gwd', 'm s-2', &
           &                            'GWD tendency of meridional wind', datatype_flt)
      grib2_desc = grib2_var(192, 128, 221, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( phy_tend_list, 'ddt_v_gwd', phy_tend%ddt_v_gwd,            &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                  & ldims=shape3d, in_group=groups("phys_tendencies") )
  
      ! &      phy_tend%ddt_v_pconv(nproma,nlev,nblks)
      cf_desc    = t_cf_var('ddt_v_pconv', 'm s-2', &
           &                            'convective tendency of meridional wind', datatype_flt)
      grib2_desc = grib2_var(0, 2, 193, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( phy_tend_list, 'ddt_v_pconv', phy_tend%ddt_v_pconv,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape3d )
  
    END IF


    !------------------------------
    ! Vertical Wind tendencies
    !------------------------------

    IF ( atm_phy_nwp_config(k_jg)%is_les_phy ) THEN
      ! &      phy_tend%ddt_w_turb(nproma,nlev+1,nblks)
      cf_desc    = t_cf_var('ddt_w_turb', 'm s-2', &
           &                            'turbulence tendency of vertical wind', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( phy_tend_list, 'ddt_w_turb', phy_tend%ddt_w_turb,          &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                  & ldims=shape3d, lrestart=.FALSE., in_group=groups("phys_tendencies"), &
                  & initval=0._wp)
    END IF


    !------------------------------
    ! Moist tracer tendencies
    !------------------------------

    ! --- Turbulence moist tracer tendencies

    ! &      phy_tend%ddt_tracer_turb(nproma,nlev,nblks,nqtendphy),          &
    cf_desc    = t_cf_var('ddt_tracer_turb', 's-1', &
         &                            'turbulence tendency of tracers', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tracer_turb', phy_tend%ddt_tracer_turb,        &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape4d,&
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    IF (lart) THEN
     ktracer=nqtendphy+nart_tendphy
    ELSE
     ktracer=nqtendphy
    ENDIF
    ALLOCATE( phy_tend%tracer_turb_ptr(ktracer) )

    !qv
    CALL add_ref( phy_tend_list, 'ddt_tracer_turb', &
                & 'ddt_qv_turb', phy_tend%tracer_turb_ptr(1)%p_3d,               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
                & t_cf_var('ddt_qv_turb', 'kg kg**-1 s**-1',                     &
                & 'turbulence tendency of specific humidity', datatype_flt),   &
                & grib2_var(192, 162, 122, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                & ldims=shape3d, lrestart=.FALSE.)
    !qc
    CALL add_ref( phy_tend_list, 'ddt_tracer_turb', &
                & 'ddt_qc_turb', phy_tend%tracer_turb_ptr(2)%p_3d,               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
                & t_cf_var('ddt_qc_turb', 'kg kg**-1 s**-1',                     &
                & 'turbulence tendency of specific cloud water', datatype_flt),&
                & grib2_var(192, 162, 201, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                & ldims=shape3d, lrestart=.FALSE.)
    !qi
    CALL add_ref( phy_tend_list, 'ddt_tracer_turb', &
                & 'ddt_qi_turb', phy_tend%tracer_turb_ptr(3)%p_3d,               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
                & t_cf_var('ddt_qi_turb', 'kg kg**-1 s**-1',                     &
                & 'turbulence tendency of specific cloud ice', datatype_flt),  &
                & grib2_var(192, 162, 202, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                & ldims=shape3d, lrestart=.FALSE.)

    ! art
    IF (lart) THEN
      CALL art_tracer_interface('turb', k_jg, kblks, phy_tend_list,  &
                & 'ddt_', ptr_arr=phy_tend%tracer_turb_ptr,          &
                & advconf=advection_config(k_jg), phy_tend=phy_tend, &
                & ldims=shape3d, tlev_source=TLEV_NNOW_RCF)
    ENDIF

    ! --- Convection moist tracer tendencies

    IF ( .NOT. atm_phy_nwp_config(k_jg)%is_les_phy ) THEN
      cf_desc    = t_cf_var('ddt_tracer_pconv', 'kg m-3 s-1', &
           &                            'convective tendency of tracers', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( phy_tend_list, 'ddt_tracer_pconv', phy_tend%ddt_tracer_pconv,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape4d_conv,&
                  & initval=0._wp, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      IF (lart) THEN
       ktracer=nqtendphy+nart_tendphy 
      ELSE
       ktracer=nqtendphy 
      ENDIF
      IF (atm_phy_nwp_config(k_jg)%ldetrain_conv_prec) ktracer = ktracer+2

      ALLOCATE( phy_tend%tracer_conv_ptr(ktracer) )

      !qv
      CALL add_ref( phy_tend_list, 'ddt_tracer_pconv',                            &
                  & 'ddt_qv_conv', phy_tend%tracer_conv_ptr(1)%p_3d,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                         &
                  & t_cf_var('ddt_qv_conv', 'kg m-3 s-1',                         &
                  & 'convective tendency of absolute humidity', datatype_flt),    &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
                  & ldims=shape3d)
      !qc
      CALL add_ref( phy_tend_list, 'ddt_tracer_pconv',                                &
                  & 'ddt_qc_conv', phy_tend%tracer_conv_ptr(2)%p_3d,                  &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
                  & t_cf_var('ddt_qc_conv', 'kg m-3 s-1',                             &
                  & 'convective tendency of cloud water mass density', datatype_flt), &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
                  & ldims=shape3d)
      !qi
      CALL add_ref( phy_tend_list, 'ddt_tracer_pconv',                              &
                  & 'ddt_qi_conv', phy_tend%tracer_conv_ptr(3)%p_3d,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                  & t_cf_var('ddt_qi_conv', 'kg m-3 s-1',                           &
                  & 'convective tendency of cloud ice mass density', datatype_flt), &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
                  & ldims=shape3d)

      IF (atm_phy_nwp_config(k_jg)%ldetrain_conv_prec) THEN
        !qr
        CALL add_ref( phy_tend_list, 'ddt_tracer_pconv',                           &
                  & 'ddt_qr_conv', phy_tend%tracer_conv_ptr(4)%p_3d,               &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                  & t_cf_var('ddt_qr_conv', 'kg m-3 s-1',                          &
                  & 'convective tendency of rain mass density', datatype_flt),     &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                  & ldims=shape3d)
        !qs
        CALL add_ref( phy_tend_list, 'ddt_tracer_pconv',                           &
                  & 'ddt_qs_conv', phy_tend%tracer_conv_ptr(5)%p_3d,               &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                  & t_cf_var('ddt_qs_conv', 'kg m-3 s-1',                          &
                  & 'convective tendency of snow mass density', datatype_flt),     &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                  & ldims=shape3d)
      ENDIF

      ! art
      IF (lart) THEN
        CALL art_tracer_interface('conv', k_jg, kblks, phy_tend_list,  &
                  & 'ddt_', ptr_arr=phy_tend%tracer_conv_ptr,          &
                  & advconf=advection_config(k_jg), phy_tend=phy_tend, &
                  & ldims=shape3d, tlev_source=TLEV_NNOW_RCF)
      ENDIF

    END IF !.not.is_les_phy

    ! --- Microphysics moist tracer tendencies

    IF ( is_variable_in_output(first_output_name_list, var_name="ddt_qv_gscp") .OR.   &
      &  is_variable_in_output(first_output_name_list, var_name="ddt_qc_gscp") .OR.   &
      &  is_variable_in_output(first_output_name_list, var_name="ddt_qi_gscp") .OR.   &
      &  is_variable_in_output(first_output_name_list, var_name="ddt_qr_gscp") .OR.   &
      &  is_variable_in_output(first_output_name_list, var_name="ddt_qs_gscp") ) THEN
      cf_desc    = t_cf_var('ddt_tracer_gscp', 's-1', &
           &                            'microphysics tendency of tracers', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( phy_tend_list, 'ddt_tracer_gscp', phy_tend%ddt_tracer_gscp,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape4d_gscp,&
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      IF (lart) THEN
       ktracer=5 + nart_tendphy 
      ELSE
       ktracer=5
      ENDIF
      ALLOCATE( phy_tend%tracer_gscp_ptr(ktracer) )

      !qv
      CALL add_ref( phy_tend_list, 'ddt_tracer_gscp',                              &
                  & 'ddt_qv_gscp', phy_tend%tracer_gscp_ptr(1)%p_3d,               &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                  & t_cf_var('ddt_qv_gscp', 'kg kg**-1 s**-1',                     &
                  & 'microphysics tendency of specific humidity', datatype_flt),   &
                  & grib2_var(192, 162, 204, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                  & ldims=shape3d)
      !qc
      CALL add_ref( phy_tend_list, 'ddt_tracer_gscp',                              &
                  & 'ddt_qc_gscp', phy_tend%tracer_gscp_ptr(2)%p_3d,               &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                  & t_cf_var('ddt_qc_gscp', 'kg kg**-1 s**-1',                     &
                  & 'microphysics tendency of specific cloud water', datatype_flt),&
                  & grib2_var(192, 162, 205, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                  & ldims=shape3d)
      !qi
      CALL add_ref( phy_tend_list, 'ddt_tracer_gscp',                              &
                  & 'ddt_qi_gscp', phy_tend%tracer_gscp_ptr(3)%p_3d,               &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                  & t_cf_var('ddt_qi_gscp', 'kg kg**-1 s**-1',                     &
                  & 'microphysics tendency of specific cloud ice', datatype_flt),  &
                  & grib2_var(192, 162, 206, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                  & ldims=shape3d)
      !qr
      CALL add_ref( phy_tend_list, 'ddt_tracer_gscp',                              &
                & 'ddt_qr_gscp', phy_tend%tracer_gscp_ptr(4)%p_3d,                 &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                            &
                & t_cf_var('ddt_qr_gscp', 'kg kg**-1 s**-1',                       &
                & 'microphysics tendency of rain', datatype_flt),                  &
                & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
                & ldims=shape3d)
      !qs
      CALL add_ref( phy_tend_list, 'ddt_tracer_gscp',                              &
                & 'ddt_qs_gscp', phy_tend%tracer_gscp_ptr(5)%p_3d,                 &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                            &
                & t_cf_var('ddt_qs_gscp', 'kg kg**-1 s**-1',                       &
                & 'microphysics tendency of snow', datatype_flt),                  &
                & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
                & ldims=shape3d)

      ! art
      IF (lart) THEN
        CALL art_tracer_interface('gscp', k_jg, kblks, phy_tend_list,  &
                  & 'ddt_', ptr_arr=phy_tend%tracer_gscp_ptr,          &
                  & advconf=advection_config(k_jg), phy_tend=phy_tend, &
                  & ldims=shape3d, tlev_source=TLEV_NNOW_RCF)
      ENDIF

    END IF


    !------------------------------
    ! TKE tendency
    !------------------------------

    IF ( .NOT. atm_phy_nwp_config(k_jg)%is_les_phy ) THEN
      !      phy_tend%ddt_tke(nproma,nlevp1,nblks)
      cf_desc    = t_cf_var('ddt_tke', 'm s-2'          , &
           &                'tendency of turbulent velocity scale', datatype_flt)
      grib2_desc = grib2_var(0, 19, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( phy_tend_list, 'ddt_tke', phy_tend%ddt_tke,             &
                  GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                & ldims=shape3dkp1 )
  
      IF (ltkecon) THEN
        lrestart = .TRUE.
      ELSE
        lrestart = .FALSE.
      ENDIF
  
      !      phy_tend%ddt_tke_pconv(nproma,nlevp1,nblks)
      cf_desc    = t_cf_var('ddt_tke_pconv', 'm**2 s**-3'          , &
           &                'TKE tendency due to sub-grid scale convection', datatype_flt)
      grib2_desc = grib2_var(0, 19, 219, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( phy_tend_list, 'ddt_tke_pconv', phy_tend%ddt_tke_pconv,   &
                  GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                & ldims=shape3dkp1, lrestart=lrestart )  


      !      phy_tend%ddt_tke_hsh(nproma,nlevp1,nblks)
      cf_desc    = t_cf_var('ddt_tke_hsh', 'm**2 s**-3'          , &
           &                'TKE tendency horizonzal shear production', datatype_flt)
      grib2_desc = grib2_var(0, 19, 221, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( phy_tend_list, 'ddt_tke_hsh', phy_tend%ddt_tke_hsh,   &
                  GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                & ldims=shape3dkp1, lrestart=.FALSE.)  
    END IF

   !Anurag Dipankar, MPIM (2013-May-31)
   !Large-scale tendencies for idealized testcases: add_var doesn't work
   !for 1D variables so using ALLOCATE-DEALLOCATE
   !Therefore, these variables can't go into restart/output
   !Initialize them all to 0 
    IF(is_ls_forcing)THEN

      ALLOCATE(phy_tend%ddt_u_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%u_ls failed')
      ENDIF
      phy_tend%ddt_u_ls = 0._wp

      ALLOCATE(phy_tend%ddt_v_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%v_ls failed')
      ENDIF
      phy_tend%ddt_v_ls = 0._wp

      ALLOCATE(phy_tend%ddt_temp_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%temp_ls failed')
      ENDIF
      phy_tend%ddt_temp_ls = 0._wp

      ALLOCATE(phy_tend%ddt_tracer_ls(klev,nqtendphy),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%tracer_ls failed')
      ENDIF
      phy_tend%ddt_tracer_ls = 0._wp
 
      ! Added by Christopher Moseley: 7 output variables for LS tendencies
      
      ALLOCATE(phy_tend%ddt_temp_subs_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%ddt_temp_subs_ls failed')
      ENDIF
      phy_tend%ddt_temp_subs_ls = 0._wp
      
      ALLOCATE(phy_tend%ddt_qv_subs_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%ddt_qv_subs_ls failed')
      ENDIF
      phy_tend%ddt_qv_subs_ls = 0._wp
      
      ALLOCATE(phy_tend%ddt_temp_adv_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%ddt_temp_adv_ls failed')
      ENDIF
      phy_tend%ddt_temp_adv_ls = 0._wp
      
      ALLOCATE(phy_tend%ddt_qv_adv_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%ddt_qv_adv_ls failed')
      ENDIF
      phy_tend%ddt_qv_adv_ls = 0._wp
      
      ALLOCATE(phy_tend%ddt_temp_nud_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%ddt_temp_nud_ls failed')
      ENDIF
      phy_tend%ddt_temp_nud_ls = 0._wp
      
      ALLOCATE(phy_tend%ddt_qv_nud_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%ddt_qv_nud_ls failed')
      ENDIF
      phy_tend%ddt_qv_nud_ls = 0._wp
      
      ALLOCATE(phy_tend%wsub(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%wsub failed')
      ENDIF
      phy_tend%wsub = 0._wp

    END IF

    CALL message('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'construction of NWP physical tendency fields finished')

END SUBROUTINE new_nwp_phy_tend_list


!
!-------------------------------------------------------------------------

END MODULE mo_nwp_phy_state
!<
