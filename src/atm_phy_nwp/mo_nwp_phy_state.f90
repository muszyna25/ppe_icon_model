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

USE mo_kind,                ONLY: wp
USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag, t_nwp_phy_tend
USE mo_impl_constants,      ONLY: success, max_char_length,           &
  &                               VINTP_METHOD_LIN,VINTP_METHOD_QV,   &
  &                               TASK_COMPUTE_RH, iedmf,             &
  &                               HINTP_TYPE_LONLAT_NNB,              &
  &                               HINTP_TYPE_LONLAT_BCTR,             &
  &                               HINTP_TYPE_LONLAT_RBF, nexlevs_rrg_vnest
USE mo_parallel_config,     ONLY: nproma
USE mo_run_config,          ONLY: nqtendphy, iqv, iqc, iqi, lart
USE mo_exception,           ONLY: message, finish !,message_text
USE mo_model_domain,        ONLY: t_patch, p_patch_local_parent
USE mo_grid_config,         ONLY: n_dom, n_dom_start
USE mo_linked_list,         ONLY: t_var_list
USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config, icpl_aero_conv
USE mo_data_turbdiff,       ONLY: ltkecon
USE mo_radiation_config,    ONLY: irad_aero
USE mo_lnd_nwp_config,      ONLY: ntiles_total, ntiles_water, nlev_soil
USE mo_var_list,            ONLY: default_var_list_settings, &
  &                               add_var, add_ref, new_var_list, delete_var_list
USE mo_var_metadata_types,  ONLY: POST_OP_SCALE, CLASS_TILE
USE mo_var_metadata,        ONLY: create_vert_interp_metadata,  &
  &                               create_hor_interp_metadata,   &
  &                               groups, vintp_types, post_op, &
  &                               new_action, actions
USE mo_nwp_parameters,      ONLY: t_phy_params
USE mo_cf_convention,       ONLY: t_cf_var
USE mo_grib2,               ONLY: t_grib2_var
USE mo_io_config,           ONLY: lflux_avg
USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_REFERENCE,        &
  &                               GRID_CELL, ZA_HYBRID, ZA_HYBRID_HALF,          &
  &                               ZA_SURFACE, ZA_HEIGHT_2M, ZA_HEIGHT_10M,       &
  &                               ZA_TOA, ZA_DEPTH_BELOW_LAND, DATATYPE_FLT32,   &
  &                               DATATYPE_PACK16, TSTEP_INSTANT,                &
  &                               TSTEP_ACCUM, TSTEP_AVG, TSTEP_MAX, TSTEP_MIN,  &
  &                               TSTEP_CONSTANT, ZA_PRESSURE_0, ZA_PRESSURE_400,&
  &                               ZA_PRESSURE_800, ZA_CLOUD_BASE, ZA_CLOUD_TOP,  &
  &                               ZA_ISOTHERM_ZERO
USE mo_physical_constants,  ONLY: grav
USE mo_ls_forcing_nml,      ONLY: is_ls_forcing

USE mo_advection_config,     ONLY: advection_config
USE mo_art_config,           ONLY: nart_tendphy
USE mo_art_tracer_interface, ONLY: art_tracer_interface
USE mo_action,               ONLY: ACTION_RESET
USE mo_les_nml,              ONLY: turb_profile_list, turb_tseries_list

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

SUBROUTINE construct_nwp_phy_state( p_patch, l_rh)

TYPE(t_patch), TARGET, INTENT(in) :: p_patch(n_dom)
LOGICAL, INTENT(IN) :: l_rh(n_dom) !< Flag. TRUE if computation of relative humidity desired

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
       &                         prm_nwp_diag_list(jg), prm_diag(jg), l_rh(jg))
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
  ENDDO


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
                     & listname, diag_list, diag, l_rh)

    INTEGER,INTENT(IN) :: klev, klevp1, kblks, k_jg !< dimension sizes

    CHARACTER(len=*),INTENT(IN)     :: listname
    CHARACTER(len=max_char_length)  :: vname_prefix
    CHARACTER(LEN=1)                :: csfc

    TYPE(t_var_list)    ,INTENT(INOUT) :: diag_list
    TYPE(t_nwp_phy_diag),INTENT(INOUT) :: diag
    LOGICAL, INTENT(IN) :: l_rh !< Flag. TRUE if computation of relative humidity desired

    ! Local variables

    INTEGER :: n_updown = 7 !> number of up/downdrafts variables

    TYPE(t_cf_var)    ::    cf_desc, new_cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d(3), shape3dsubs(3), shape3dsubsw(3)
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
 
    ibits = DATATYPE_PACK16 ! bits "entropy" of horizontal slice

    shape2d      = (/nproma,           kblks            /)
    shape3d      = (/nproma, klev,     kblks            /)
    shape3dkp1   = (/nproma, klevp1,   kblks            /)
    shape3dsubs  = (/nproma, kblks,    ntiles_total     /)
    shape3dsubsw = (/nproma, kblks,    ntiles_total+ntiles_water /)

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
      &                   DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 77, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rain_gsp_rate', diag%rain_gsp_rate,            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & isteptype=TSTEP_INSTANT )


    ! &      diag%snow_gsp_rate(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_gsp_rate', 'kg m-2 s-1', 'gridscale snow rate', &
      &                   DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 56, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'snow_gsp_rate', diag%snow_gsp_rate,            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & isteptype=TSTEP_INSTANT )

    ! For graupel scheme 
    SELECT CASE (atm_phy_nwp_config(k_jg)%inwp_gscp)
    CASE (2,4,5,6)
      
      ! &      diag%graupel_gsp_rate(nproma,nblks_c)
      cf_desc    = t_cf_var('graupel_gsp_rate', 'kg m-2 s-1', 'gridscale graupel rate', &
        &                   DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( diag_list, 'graupel_gsp_rate', diag%graupel_gsp_rate,            &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d,                                             &
                  & isteptype=TSTEP_INSTANT )
    END SELECT

    !For two moment microphysics
    SELECT CASE (atm_phy_nwp_config(k_jg)%inwp_gscp)
    CASE (4,5,6)

       ! &      diag%ice_gsp_rate(nproma,nblks_c)
      cf_desc    = t_cf_var('ice_gsp_rate', 'kg m-2 s-1', 'gridscale ice rate', &
        &                   DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( diag_list, 'ice_gsp_rate', diag%ice_gsp_rate,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d, isteptype=TSTEP_INSTANT )
      
       ! &      diag%hail_gsp_rate(nproma,nblks_c)
      cf_desc    = t_cf_var('hail_gsp_rate', 'kg m-2 s-1', 'gridscale hail rate', &
        &                   DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( diag_list, 'hail_gsp_rate', diag%hail_gsp_rate,            &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d, isteptype=TSTEP_INSTANT )

    END SELECT


    ! &      diag%rain_con_rate(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_con_rate', 'kg m-2 s-1', 'convective rain rate', &
      &                   DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 76, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rain_con_rate', diag%rain_con_rate,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & isteptype=TSTEP_INSTANT )


    ! &      diag%snow_con_rate(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_con_rate', 'kg m-2 s-1', 'convective snow rate', &
      &                   DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 55, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'snow_con_rate', diag%snow_con_rate,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & isteptype=TSTEP_INSTANT )

    ! &      diag%rain_con_rate_3d(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('rain_con_rate_3d', 'kg m-2 s-1',                &
      &          '3d convective rain rate', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 76, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rain_con_rate_3d', diag%rain_con_rate_3d,       &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,&
                & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                & ldims=shape3dkp1,                                           &
                & isteptype=TSTEP_INSTANT )


    ! &      diag%snow_con_rate_3d(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('snow_con_rate_3d', 'kg m-2 s-1',                   &
      &          '3d convective snow rate', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 55, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'snow_con_rate_3d', diag%snow_con_rate_3d,       &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,&
                & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                & ldims=shape3dkp1,                                           &
                & isteptype=TSTEP_INSTANT )


    IF ( atm_phy_nwp_config(k_jg)%inwp_turb == iedmf ) THEN

      ! &      diag%rain_edmf_rate_3d(nproma,nlevp1,nblks_c)
      cf_desc    = t_cf_var('rain_edmf_rate_3d', 'kg m-2 s-1',                &
        &          '3d EDMF convective rain rate', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( diag_list, 'rain_edmf_rate_3d', diag%rain_edmf_rate_3d,       &
                  & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,&
                  & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                  & ldims=shape3dkp1,                                           &
                  & isteptype=TSTEP_INSTANT )
      
      
      ! &      diag%snow_edmf_rate_3d(nproma,nlevp1,nblks_c)
      cf_desc    = t_cf_var('snow_edmf_rate_3d', 'kg m-2 s-1',                   &
        &          '3d EDMF convective snow rate', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( diag_list, 'snow_edmf_rate_3d', diag%snow_edmf_rate_3d,       &
                  & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,&
                  & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                  & ldims=shape3dkp1,                                           &
                  & isteptype=TSTEP_INSTANT )

    ENDIF


    ! &      diag%rain_gsp(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_gsp ', 'kg m-2 ', 'gridscale rain ', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 77, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rain_gsp', diag%rain_gsp,                      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d, in_group=groups("precip_vars"),             &
                & isteptype=TSTEP_ACCUM ,                                    &
                & hor_interp=create_hor_interp_metadata(                     &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                   &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB                     &
                & ) )

    ! &      diag%rain_gsp0(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_gsp0', 'kg m-2 ', 'gridscale rain0', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 77, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rain_gsp0', diag%rain_gsp0,                    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & loutput=.false., lrestart=.TRUE.,                          &
                & isteptype=TSTEP_ACCUM )


    ! &      diag%snow_gsp(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_gsp', 'kg m-2 ', 'gridscale snow', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 56, ibits, GRID_REFERENCE, GRID_CELL)
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
    cf_desc    = t_cf_var('snow_gsp0', 'kg m-2 ', 'gridscale snow0', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 56, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'snow_gsp0', diag%snow_gsp0,                    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & loutput=.false., lrestart=.TRUE.,                          &
                & isteptype=TSTEP_ACCUM )


    ! &      diag%rain_con(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_con', 'kg m-2 ', 'convective rain', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 76, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rain_con', diag%rain_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, in_group=groups("precip_vars"),              &
                & isteptype=TSTEP_ACCUM ,                                     &
                & hor_interp=create_hor_interp_metadata(                     &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                   &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB                     &
                & ) )

    ! &      diag%rain_con0(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_con0', 'kg m-2 ', 'convective rain0', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 76, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rain_con0', diag%rain_con0,                    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & loutput=.false., lrestart=.TRUE.,                          &
                & isteptype=TSTEP_ACCUM )


    ! &      diag%snow_con(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_con', 'kg m-2', 'convective snow', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 55, ibits, GRID_REFERENCE, GRID_CELL)
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
    cf_desc    = t_cf_var('snow_con0', 'kg m-2', 'convective snow0', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 55, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'snow_con0', diag%snow_con0,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & loutput=.false., lrestart=.TRUE.,                           &
                & isteptype=TSTEP_ACCUM )


    !Surface precipitation variables for graupel scheme and two moment microphysics
    SELECT CASE (atm_phy_nwp_config(k_jg)%inwp_gscp)
    CASE (2,4,5,6)
       ! &      diag%graupel_gsp(nproma,nblks_c)
      cf_desc    = t_cf_var('graupel_gsp', 'kg m-2', 'gridscale graupel',      &
        &                   DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( diag_list, 'graupel_gsp', diag%graupel_gsp,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d, in_group=groups("precip_vars"),             &
                  & isteptype=TSTEP_ACCUM )
    END SELECT

    SELECT CASE (atm_phy_nwp_config(k_jg)%inwp_gscp)
    CASE (4,5,6)

       ! &      diag%ice_gsp(nproma,nblks_c)
      cf_desc    = t_cf_var('ice_gsp', 'kg m-2', 'gridscale ice', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( diag_list, 'ice_gsp', diag%ice_gsp,                        &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d, in_group=groups("precip_vars"),             &
                  & isteptype=TSTEP_ACCUM )
      

       ! &      diag%hail_gsp(nproma,nblks_c)
      cf_desc    = t_cf_var('hail_gsp', 'kg m-2', 'gridscale hail', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( diag_list, 'hail_gsp', diag%hail_gsp,                      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d, in_group=groups("precip_vars"),             &
                  & isteptype=TSTEP_ACCUM )

    END SELECT


    ! &      diag%tot_prec(nproma,nblks_c)
    cf_desc    = t_cf_var('tot_prec', 'kg m-2', 'total precip', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 52, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tot_prec', diag%tot_prec,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & in_group=groups("precip_vars"),                             &
                & isteptype=TSTEP_ACCUM ,                                     &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB                      &
                & ) )


    ! &      diag%tot_prec_rate_avg(nproma,nblks_c)
    cf_desc    = t_cf_var('tot_prec_rate_avg', 'kg m-2 s-1',                  &
      &          'total precip rate, time average', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 52, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tot_prec_rate_avg', diag%tot_prec_rate_avg,     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE.,                            &
                & isteptype=TSTEP_AVG )


    ! &      diag%con_prec_rate_avg(nproma,nblks_c)
    cf_desc    = t_cf_var('con_prec_rate_avg', 'kg m-2 s-1',                  &
      &          'convective precip rate, time average', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 1, 10, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'con_prec_rate_avg', diag%con_prec_rate_avg,     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,  lrestart=.FALSE.,                           &
                & in_group=groups("additional_precip_vars"),                  &
                & isteptype=TSTEP_AVG )

    ! &      diag%gsp_prec_rate_avg(nproma,nblks_c)
    cf_desc    = t_cf_var('gsp_prec_rate_avg', 'kg m-2 s-1',                  &
      &          'gridscale precip rate, time average', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'gsp_prec_rate_avg', diag%gsp_prec_rate_avg,     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE.,                            &
                & in_group=groups("additional_precip_vars"),                  &
                & isteptype=TSTEP_AVG )

    ! &      diag%cape(nproma,nblks_c)
    cf_desc    = t_cf_var('cape', 'J kg-1 ', 'conv avail pot energy', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 7, 6, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'cape', diag%cape,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE.,                            &
                & in_group=groups("additional_precip_vars"),                  &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB                      &
                & ) )

    ! &      diag%gust10(nproma,nblks_c)
    cf_desc    = t_cf_var('gust10', 'm s-1 ', 'gust at 10 m', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 2, 22, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'gust10', diag%gust10,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,  &
                & ldims=shape2d, lrestart=.TRUE., in_group=groups("pbl_vars"), &
                & isteptype=TSTEP_MAX,                                         &
                & initval_r=0._wp, resetval_r=0._wp,                           &
                & action_list=actions(new_action(ACTION_RESET,'PT03H')) )

    ! &      diag%dyn_gust(nproma,nblks_c)
    cf_desc    = t_cf_var('dyn_gust', 'm s-1 ', 'dynamical gust', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'dyn_gust', diag%dyn_gust,                        &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,  &
                & ldims=shape2d, lrestart=.FALSE., isteptype=TSTEP_INSTANT,    &
                & loutput=.TRUE.                                               )

    ! &      diag%con_gust(nproma,nblks_c)
    cf_desc    = t_cf_var('con_gust', 'm s-1 ', 'convective contribution to wind gust', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'con_gust', diag%con_gust,                        &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,  &
                & ldims=shape2d, lrestart=.TRUE., isteptype=TSTEP_INSTANT,     &
                & loutput=.TRUE. )
   
    ! &      diag%rain_upd(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_upd', 'kg m-2 s-1', 'rain in updroughts', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rain_upd', diag%rain_upd,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE. )

    ! &      diag%con_udd(nproma,nlev,nblks,8)
    cf_desc    = t_cf_var('con_udd', 'unit ', 'convective up/downdraft fields', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'con_udd', diag%con_udd,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,     &
                & ldims=(/nproma,klev,kblks,n_updown/),                       &
                & lrestart=.FALSE., loutput=.FALSE. )

    ! &      diag%mbas_con(nproma,nblks_c)
    cf_desc    = t_cf_var('mbas_con', '', 'cloud base level index', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 6, 194, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'mbas_con', diag%mbas_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    ! &      diag%mtop_con(nproma,nblks_c)
    cf_desc    = t_cf_var('mtop_con', '', 'cloud top level index', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 6, 195, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'mtop_con', diag%mtop_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    ! &      diag%locum(nproma,nblks_c)
    cf_desc    = t_cf_var('locum', '', 'convective activity indicator', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'locum', diag%locum,                             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    ! &      diag%ldshcv(nproma,nblks_c)
    cf_desc    = t_cf_var('ldshcv', '', 'shallow convection indicator', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'ldshcv', diag%ldshcv,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    ! &      diag%ktype(nproma,nblks_c)
    cf_desc    = t_cf_var('ktype', '', 'type of convection', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'ktype', diag%ktype,                              &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,  &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ) )

    ! &      diag%k850(nproma,nblks_c)
    cf_desc    = t_cf_var('k850', '', 'level index corresponding to the HAG of the 850hPa level', &
      &                   DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'k850', diag%k850,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,  &
                & isteptype=TSTEP_CONSTANT,                                    &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ) )

    ! &      diag%k950(nproma,nblks_c)
    cf_desc    = t_cf_var('k950', '', 'level index corresponding to the HAG of the 950hPa level', &
      &                   DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'k950', diag%k950,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,  &
                & isteptype=TSTEP_CONSTANT,                                    &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ) )

    ! &      diag%ktop_envel(nproma,nblks_c)
    cf_desc    = t_cf_var('ktop_envel', '', 'level index of upper boundary of SSO envelope layer', &
      &                   DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'ktop_envel', diag%ktop_envel,                    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.FALSE.  )

    ! &      diag%clc(nproma,nlev,nblks_c)
    cf_desc      = t_cf_var('clc', '',  'cloud cover', DATATYPE_FLT32)
    new_cf_desc  = t_cf_var('clc', '%', 'cloud cover', DATATYPE_FLT32)
    grib2_desc   = t_grib2_var(0, 6, 22, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'clc', diag%clc,                                 &
      & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,               &
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
    cf_desc      = t_cf_var('clct', '',  'total cloud cover', DATATYPE_FLT32)
    new_cf_desc  = t_cf_var('clct', '%', 'total cloud cover', DATATYPE_FLT32)
    grib2_desc   = t_grib2_var(0, 6, 1, ibits, GRID_REFERENCE, GRID_CELL)
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
      &                     DATATYPE_FLT32)
    grib2_desc   = t_grib2_var(0, 6, 199, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'clct_mod', diag%clct_mod,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & in_group=groups("additional_precip_vars"),                            &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF))

    ! &      diag%clch(nproma,nblks_c)
    cf_desc      = t_cf_var('clch', '', 'high_level_clouds',  DATATYPE_FLT32)
    new_cf_desc  = t_cf_var('clch', '%', 'high_level_clouds', DATATYPE_FLT32)
    grib2_desc   = t_grib2_var(0, 6, 22, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'clch', diag%clch,                               &
      & GRID_UNSTRUCTURED_CELL, ZA_PRESSURE_0, cf_desc, grib2_desc,           &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                          &
      &                 new_cf=new_cf_desc))

    ! &      diag%clcm(nproma,nblks_c)
    cf_desc      = t_cf_var('clcm', '',  'mid_level_clouds', DATATYPE_FLT32)
    new_cf_desc  = t_cf_var('clcm', '%', 'mid_level_clouds', DATATYPE_FLT32)
    grib2_desc   = t_grib2_var(0, 6, 22, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'clcm', diag%clcm,                               &
      & GRID_UNSTRUCTURED_CELL, ZA_PRESSURE_400, cf_desc, grib2_desc,         &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                          &
      &                 new_cf=new_cf_desc))

    ! &      diag%clcl(nproma,nblks_c)
    cf_desc      = t_cf_var('clcl', '',  'low_level_clouds', DATATYPE_FLT32)
    new_cf_desc  = t_cf_var('clcl', '%', 'low_level_clouds', DATATYPE_FLT32)
    grib2_desc   = t_grib2_var(0, 6, 22, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'clcl', diag%clcl,                               &
      & GRID_UNSTRUCTURED_CELL, ZA_PRESSURE_800, cf_desc, grib2_desc,         &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                          &
      &                 new_cf=new_cf_desc))

    ! &      diag%cldepth(nproma,nblks_c)
    cf_desc      = t_cf_var('cldepth', '',  'modified cloud depth for media', DATATYPE_FLT32)
    grib2_desc   = t_grib2_var(0, 6, 198, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'cldepth', diag%cldepth,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF))


    ! &      diag%hbas_con(nproma,nblks_c)
    cf_desc    = t_cf_var('hbas_con', 'm', 'height_of_convective_cloud_base', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 6, 26, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'hbas_con', diag%hbas_con,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE, cf_desc, grib2_desc,           &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF))

    ! &      diag%htop_con(nproma,nblks_c)
    cf_desc    = t_cf_var('htop_con', 'm', 'height_of_convective_cloud_top', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 6, 27, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'htop_con', diag%htop_con,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF))

    ! &      diag%htop_dc(nproma,nblks_c)
    cf_desc    = t_cf_var('htop_dc', 'm', 'height_of_top_of_dry_convection', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 6, 196, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'htop_dc', diag%htop_dc,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF))

    ! &      diag%acdnc(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('acdnc', 'm-3', 'cloud droplet number concentration', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 6, 30, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'acdnc', diag%acdnc,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,               &
      & ldims=shape3d, lrestart=.FALSE.,                                      &
      & isteptype=TSTEP_CONSTANT )


    !      diag%tot_cld(nproma,nlev,nblks_c,3)
    cf_desc    = t_cf_var('tot_cld', ' ','total cloud variables (qv,qc,qi)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tot_cld', diag%tot_cld,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,       &
                & ldims=(/nproma,klev,kblks,3/) ,                               &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,         &
                & initval_r=0.0_wp,                                             &
                & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
                &                                       fallback_type=HINTP_TYPE_LONLAT_RBF))

    ALLOCATE( diag%tot_ptr(kcloud))
    vname_prefix='tot_'

           !QV
        CALL add_ref( diag_list, 'tot_cld',                                            &
                    & TRIM(vname_prefix)//'qv_dia', diag%tot_ptr(iqv)%p_3d,            &
                    & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                               &
                    & t_cf_var(TRIM(vname_prefix)//'qv_dia', 'kg kg-1',                &
                    &          'total_specific_humidity_(diagnostic)', DATATYPE_FLT32),&
                    & t_grib2_var(0, 1, 211, ibits, GRID_REFERENCE, GRID_CELL),        &
                    & ldims=shape3d,                                                   &
                    & vert_interp=create_vert_interp_metadata(                         &
                    &             vert_intp_type=vintp_types("P", "Z", "I"),           &
                    &             vert_intp_method=VINTP_METHOD_QV,                    &
                    &             l_satlimit=.FALSE.,                                  & 
                    &             lower_limit=2.5e-6_wp, l_restore_pbldev=.FALSE. ),   &
                    & in_group=groups("cloud_diag") )

           !QC
        CALL add_ref( diag_list, 'tot_cld',                                            &
                    & TRIM(vname_prefix)//'qc_dia', diag%tot_ptr(iqc)%p_3d,            &
                    & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                               &
                    & t_cf_var(TRIM(vname_prefix)//'qc_dia', 'kg kg-1',                &
                    & 'total_specific_cloud_water_content_(diagnostic)', DATATYPE_FLT32),&
                    & t_grib2_var(0, 1, 212, ibits, GRID_REFERENCE, GRID_CELL),        &
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
                    & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                               &
                    & t_cf_var(TRIM(vname_prefix)//'qi_dia', 'kg kg-1',                &
                    & 'total_specific_cloud_ice_content_(diagnostic)', DATATYPE_FLT32),&
                    & t_grib2_var(0, 1, 213, ibits, GRID_REFERENCE, GRID_CELL),        &
                    & ldims=shape3d,                                                   &
                    & vert_interp=create_vert_interp_metadata(                         &
                    &             vert_intp_type=vintp_types("P","Z","I"),             &
                    &             vert_intp_method=VINTP_METHOD_LIN,                   &
                    &             l_loglin=.FALSE.,                                    &
                    &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,              &
                    &             lower_limit=0._wp ),                                 &
                    & in_group=groups("cloud_diag") )






    !      diag%tot_cld_vi(nproma,nblks_c,3)
    cf_desc     = t_cf_var('tot_cld_vi', 'kg m-2','vertical integr total cloud variables', DATATYPE_FLT32)
    grib2_desc   = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
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
      & 'total column integrated water vapour (diagnostic)', DATATYPE_FLT32),   &
      & t_grib2_var( 0, 1, 214, ibits, GRID_REFERENCE, GRID_CELL), &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("additional_precip_vars"))

    !TQC_DIA
    CALL add_ref( diag_list, 'tot_cld_vi',                         &
      & 'tqc_dia', diag%tci_ptr(iqc)%p_2d,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
      & t_cf_var('tqc_dia', 'kg m**-2',                            &
      & 'total column integrated cloud water (diagnostic)', DATATYPE_FLT32),    &
      & t_grib2_var( 0, 1, 215, ibits, GRID_REFERENCE, GRID_CELL), &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("additional_precip_vars"))

    !TQI_DIA
    CALL add_ref( diag_list, 'tot_cld_vi',                         &
      & 'tqi_dia', diag%tci_ptr(iqi)%p_2d,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
      & t_cf_var('tqi_dia', 'kg m**-2',                            &
      & 'total column integrated cloud ice (diagnostic)', DATATYPE_FLT32),      &
      & t_grib2_var(0, 1, 216, ibits, GRID_REFERENCE, GRID_CELL),  &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("additional_precip_vars"))



    !      diag%tot_cld_vi_avg(nproma,nblks_c,3)
    cf_desc    = t_cf_var('tot_cld_vi_avg', 'unit ','vertical integr total cloud variables', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
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
                & 'column integrated water vapour (diagnostic)_avg',           &
                & DATATYPE_FLT32),                                             &
                & t_grib2_var( 0, 1, 214, ibits, GRID_REFERENCE, GRID_CELL),   &
                & ldims=shape2d, lrestart=.FALSE.,                             &
                & isteptype=TSTEP_AVG )

    ! TQC_DIA_AVG
    CALL add_ref( diag_list, 'tot_cld_vi_avg',      &
                & TRIM(vname_prefix)//'qc', diag%tav_ptr(iqc)%p_2d,            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var(TRIM(vname_prefix)//'qc', '',                       &
                & 'tci_specific_cloud_water_content (diagnostic)_avg',         &
                & DATATYPE_FLT32),                                             &
                & t_grib2_var(0, 1, 215, ibits, GRID_REFERENCE, GRID_CELL),    &
                & ldims=shape2d, lrestart=.FALSE.,                             &
                & isteptype=TSTEP_AVG )

    ! TQI_DIA_AVG
    CALL add_ref( diag_list, 'tot_cld_vi_avg',          &
                & TRIM(vname_prefix)//'qi', diag%tav_ptr(iqi)%p_2d,            & 
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                & t_cf_var(TRIM(vname_prefix)//'qi', '',                       &
                & 'tci_specific_cloud_ice_content (diagnostic)_avg',           &
                & DATATYPE_FLT32),                                             &
                & t_grib2_var(0, 1, 216, ibits, GRID_REFERENCE, GRID_CELL),    &
                & ldims=shape2d, lrestart=.FALSE.,                             &
                & isteptype=TSTEP_AVG )


    ! &      diag%clct_avg(nproma,nblks_c)
    cf_desc      = t_cf_var('clct_avg', '%', 'total cloud cover time avg', &
      &            DATATYPE_FLT32)
    grib2_desc   = t_grib2_var(0, 6, 1, ibits, GRID_REFERENCE, GRID_CELL)
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
          &                    DATATYPE_FLT32)
        new_cf_desc = t_cf_var('albdif', '%','Shortwave albedo for diffuse radiation', &
          &                    DATATYPE_FLT32)
        grib2_desc  = t_grib2_var(0, 19, 1, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'albdif', diag%albdif,                         &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
          & ldims=shape2d, in_group=groups("rad_vars"),                         &
          & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                        &
          &                 new_cf=new_cf_desc) )


        !        diag%albvisdif    (nproma,       nblks),          &
        cf_desc     = t_cf_var('albvisdif', '', 'UV visible albedo for diffuse radiation', &
          &                    DATATYPE_FLT32)
        new_cf_desc = t_cf_var('albvisdif', '%','UV visible albedo for diffuse radiation', &
          &                    DATATYPE_FLT32)
        grib2_desc  = t_grib2_var(0, 19, 222, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'albvisdif', diag%albvisdif,                   &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
          & ldims=shape2d, in_group=groups("rad_vars"),                         &
          & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                        &
          &                 new_cf=new_cf_desc) )


        !        diag%albvisdir    (nproma,       nblks),          &
        cf_desc     = t_cf_var('albvisdir', '', 'UV visible albedo for direct radiation', &
          &                    DATATYPE_FLT32)
        new_cf_desc = t_cf_var('albvisdir', '%','UV visible albedo for direct radiation', &
          &                    DATATYPE_FLT32)
        grib2_desc  = t_grib2_var(192, 128, 15, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'albvisdir', diag%albvisdir,                   &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
          & ldims=shape2d,                                                      &
          & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                        &
          &                 new_cf=new_cf_desc) )


        !        diag%albnirdif    (nproma,       nblks),          &
        cf_desc     = t_cf_var('albnirdif', '',  'Near IR albedo for diffuse radiation',&
          &                    DATATYPE_FLT32)
        new_cf_desc = t_cf_var('albnirdif', '%', 'Near IR albedo for diffuse radiation',&
          &                    DATATYPE_FLT32)
        grib2_desc  = t_grib2_var(0, 19, 223, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'albnirdif', diag%albnirdif,                   &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
          & ldims=shape2d, in_group=groups("rad_vars"),                         &
          & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                        &
          &                 new_cf=new_cf_desc) )


        !        diag%albnirdir    (nproma,       nblks),          &
        cf_desc     = t_cf_var('albnirdir', '',  'Near IR albedo for direct radiation',&
          &                    DATATYPE_FLT32)
        new_cf_desc = t_cf_var('albnirdir', '%', 'Near IR albedo for direct radiation',&
          &                    DATATYPE_FLT32)
        grib2_desc  = t_grib2_var(192, 128, 17, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'albnirdir', diag%albnirdir,                   &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
          & ldims=shape2d,                                                      &
          & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                        &
          &                 new_cf=new_cf_desc) )


        ! These variables only make sense if the land-surface scheme is switched on.
        IF ( atm_phy_nwp_config(k_jg)%inwp_surface == 1 ) THEN

          !        diag%albdif_t (nproma, nblks, ntiles_total+ntiles_water)
          cf_desc    = t_cf_var('albdif_t', '', &
            &                   'tile-based shortwave albedo for diffusive radiation',&
            &                   DATATYPE_FLT32)
          grib2_desc = t_grib2_var(0, 19, 1, ibits, GRID_REFERENCE, GRID_CELL)
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
               & t_cf_var('albdif_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),    &
               & t_grib2_var(0, 19, 1, ibits, GRID_REFERENCE, GRID_CELL),      &
               & var_class=CLASS_TILE,                                         &
               & ldims=shape2d, lrestart=.TRUE.                                )
          ENDDO


          !        diag%albvisdif_t (nproma, nblks, ntiles_total+ntiles_water)
          cf_desc    = t_cf_var('albvisdif_t', '', &
            &                   'tile-based UV visible albedo for diffusive radiation',&
            &                   DATATYPE_FLT32)
          grib2_desc = t_grib2_var(0, 19, 222, ibits, GRID_REFERENCE, GRID_CELL)
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
               & t_cf_var('albvisdif_t_'//TRIM(csfc), '', '', DATATYPE_FLT32), &
               & t_grib2_var(0, 19, 222, ibits, GRID_REFERENCE, GRID_CELL),    &
               & var_class=CLASS_TILE,                                         &
               & ldims=shape2d, lrestart=.TRUE.                                )
          ENDDO


          !        diag%albnirdif_t (nproma, nblks, ntiles_total+ntiles_water)
          cf_desc    = t_cf_var('albnirdif_t', '', &
            &                   'tile-based near IR albedo for diffuse radiation',&
            &                   DATATYPE_FLT32)
          grib2_desc = t_grib2_var(0, 19, 223, ibits, GRID_REFERENCE, GRID_CELL)
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
               & t_cf_var('albnirdif_t_'//TRIM(csfc), '', '', DATATYPE_FLT32), &
               & t_grib2_var(0, 19, 223, ibits, GRID_REFERENCE, GRID_CELL),    &
               & var_class=CLASS_TILE,                                         &
               & ldims=shape2d, lrestart=.TRUE.                                )
          ENDDO



          ! &      diag%swflxsfc_t(nproma,nblks_c,ntiles_total+ntiles_water)
          cf_desc    = t_cf_var('sob_s_t', 'W m-2', 'tile-based shortwave net flux at surface', &
               &                DATATYPE_FLT32)
          grib2_desc = t_grib2_var(0, 4, 9, ibits, GRID_REFERENCE, GRID_CELL)
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
               & t_cf_var('swflxsfc_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),  &
               & t_grib2_var(0, 4, 9, ibits, GRID_REFERENCE, GRID_CELL),       &
               & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,               &
               & var_class=CLASS_TILE,                                         &
               & in_group=groups("rad_vars"))
          ENDDO

          ! &      diag%lwflxsfc_t(nproma,nblks_c,ntiles_total+ntiles_water)
          cf_desc    = t_cf_var('thb_s_t', 'W m-2', 'tile_based longwave net flux at surface', &
               &                DATATYPE_FLT32)
          grib2_desc = t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
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
               & t_cf_var('lwflxsfc_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),  &
               & t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL),       &
               & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,               &
               & var_class=CLASS_TILE,                                         &
               & in_group=groups("rad_vars"))
          ENDDO
        ENDIF


        !        diag%cosmu0    (nproma,nblks)
        cf_desc    = t_cf_var('cosmu0', '-', 'Cosine of solar zenith angle', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(192, 214, 1, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'cosmu0', diag%cosmu0,                         &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d)


        ! &      diag%tsfctrad(nproma,nblks_c)
        cf_desc    = t_cf_var('tsfctrad', 'K', 'surface temperature at trad', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tsfctrad', diag%tsfctrad,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d)


        ! &       diag% flxdwswtoa(nproma,       nblks),          &
        cf_desc    = t_cf_var('sod_t', 'W m-2', 'downward shortwave flux at TOA', &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 4, 7, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'sod_t', diag%flxdwswtoa,              &
          & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,             &
          & ldims=shape2d, lrestart=.FALSE.,in_group=groups("rad_vars"))


        ! &      diag%swflxsfc(nproma,nblks_c)
        cf_desc    = t_cf_var('sob_s', 'W m-2', 'shortwave net flux at surface', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 4, 9, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'sob_s', diag%swflxsfc,                        &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
          & ldims=shape2d,                                                      &
          & in_group=groups("rad_vars"))


        ! &      diag%swflxtoa(nproma,nblks_c)
        cf_desc    = t_cf_var('sob_t', 'W m-2', 'shortwave net flux at TOA', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 4, 9, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'sob_t', diag%swflxtoa,                        &
          & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,                &
          & ldims=shape2d, lrestart=.FALSE.,                                    &
          & in_group=groups("rad_vars"))

        ! &      diag%trsolclr_sfc(nproma,nblks_c)
        cf_desc    = t_cf_var('trsolclr_sfc', '', 'shortwave clear-sky transmisivity at suface', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'trsolclr_sfc', diag%trsolclr_sfc,             &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
          & ldims=shape2d, loutput=.FALSE.                                      )

        ! &      diag%trsol_up_toa(nproma,nblks_c)
        cf_desc    = t_cf_var('trsol_up_toa', '', 'shortwave upward transmisivity at TOA', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'trsol_up_toa', diag%trsol_up_toa,             &
          & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,                &
          & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE.                    )

        ! &      diag%trsol_up_sfc(nproma,nblks_c)
        cf_desc    = t_cf_var('trsol_up_sfc', '', 'shortwave upward transmisivity at surface', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'trsol_up_sfc', diag%trsol_up_sfc,             &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
          & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE.                    )

        ! &      diag%trsol_par_sfc(nproma,nblks_c)
        cf_desc    = t_cf_var('trsol_par_sfc', '', 'photosynthetically active downward transmisivity at surface', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'trsol_par_sfc', diag%trsol_par_sfc,           &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
          & ldims=shape2d, lrestart=.TRUE., loutput=.FALSE.                     )

        ! &      diag%trsol_dn_sfc_diff(nproma,nblks_c)
        cf_desc    = t_cf_var('trsol_dn_sfc_diff', '', 'shortwave diffuse downward transmisivity at surface', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'trsol_dn_sfc_diff', diag%trsol_dn_sfc_diff,   &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
          & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE.                    )

        ! &      diag%swflx_up_toa(nproma,nblks_c)
        cf_desc    = t_cf_var('sou_t', 'W m-2', 'shortwave upward flux at TOA', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 4, 8, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'sou_t', diag%swflx_up_toa,                    &
          & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,                &
          & ldims=shape2d, lrestart=.FALSE., in_group=groups("rad_vars")        )

        ! &      diag%swflx_up_sfc(nproma,nblks_c)
        cf_desc    = t_cf_var('sou_s', 'W m-2', 'shortwave upward flux at surface', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 4, 8, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'sou_s', diag%swflx_up_sfc,                    &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
          & ldims=shape2d, lrestart=.FALSE., in_group=groups("rad_vars")        )

        ! &      diag%swflx_par_sfc(nproma,nblks_c)
        cf_desc    = t_cf_var('swflx_par_sfc', 'W m-2', 'downward photosynthetically active flux at surface', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 4, 10, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'swflx_par_sfc', diag%swflx_par_sfc,           &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
          & ldims=shape2d, lrestart=.TRUE., in_group=groups("rad_vars")         )

        ! &      diag%swflx_dn_sfc_diff(nproma,nblks_c)
        cf_desc    = t_cf_var('sodifd_s', 'W m-2', 'shortwave diffuse downward flux at surface', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 4, 199, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'sodifd_s', diag%swflx_dn_sfc_diff,              &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
          & ldims=shape2d, lrestart=.FALSE., in_group=groups("rad_vars")        )

        ! &      diag%lwflxsfc(nproma,nblks_c)
        cf_desc    = t_cf_var('thb_s', 'W m-2', 'longwave net flux at surface', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'thb_s', diag%lwflxsfc,                        &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
          & ldims=shape2d,                                                      &
          & in_group=groups("rad_vars"))

        ! This is an auxiliary storage field needed because the final diagnostic quantity (below)
        ! is updated each physics time step following the time evolution of the ground temperature
        ! &      diag%lwflx_up_sfc_rs(nproma,nblks_c)
        cf_desc    = t_cf_var('lwflx_up_sfc_rs', 'W m-2', 'longwave upward flux at surface', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'lwflx_up_sfc_rs', diag%lwflx_up_sfc_rs,       &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
          & ldims=shape2d,  lrestart=.FALSE., loutput=.FALSE.                   )

        ! &      diag%lwflx_up_sfc(nproma,nblks_c)
        cf_desc    = t_cf_var('thu_s', 'W m-2', 'longwave upward flux at surface', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 5, 4, ibits, GRID_REFERENCE, GRID_CELL)
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
          &                          TRIM(long_name), DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%lwflxtoa_a,                  &
          & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,               &
          & ldims=shape2d, isteptype=a_steptype, in_group=groups("rad_vars") )


        ! &      diag%swflxtoa_a(nproma,nblks_c)
        WRITE(name,'(A,A5)') TRIM(prefix),"sob_t"
        WRITE(long_name,'(A26,A4,A18)') "TOA net solar radiation ", meaning, &
                                      &" since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), &
          &                         TRIM(long_name), DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 4, 9, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name) , diag%swflxtoa_a,                 &
          & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,               &
          & ldims=shape2d,                                                     &
          & isteptype=a_steptype, in_group=groups("rad_vars") )

        
        ! &      diag%lwflxsfc_a(nproma,nblks_c)
        WRITE(name,'(A,A5)') TRIM(prefix),"thb_s"
        WRITE(long_name,'(A30,A4,A18)') "surface net thermal radiation ", meaning, &
                                      &" since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), &
          &                      TRIM(long_name), DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%lwflxsfc_a,                 &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d,                                                    &
          & isteptype=a_steptype, in_group=groups("rad_vars"))


        ! &      diag%swflxsfc_a(nproma,nblks_c)
        WRITE(name,'(A,A5)') TRIM(prefix),"sob_s"
        WRITE(long_name,'(A30,A4,A18)') "Surface net solar radiation ", meaning, &
                                      &" since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), &
          &                    TRIM(long_name), DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 4, 9, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%swflxsfc_a,                 &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d,                                                    &
          & isteptype=a_steptype, in_group=groups("rad_vars"))


        ! &       diag%asod_t(nproma,nblks)
        WRITE(name,'(A,A5)') TRIM(prefix),"sod_t"
        WRITE(long_name,'(A30,A4,A18)') "Top down solar radiation ", meaning, &
                                      &" since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 4, 7, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%asod_t,                    &
          & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,             &
          & ldims=shape2d,                                                   &
          & isteptype=a_steptype, in_group=groups("rad_vars"))


        ! &       diag%asou_t(nproma,nblks)
        WRITE(name,'(A,A5)') TRIM(prefix),"sou_t"
        WRITE(long_name,'(A30,A4,A18)') "Top up solar radiation ", meaning, &
                                      &" since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 4, 8, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%asou_t,                    &
          & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,             &
          & ldims=shape2d,                                                   &
          & isteptype=a_steptype, in_group=groups("rad_vars"))


        ! &       diag%athd_s(nproma,nblks)
        WRITE(name,'(A,A5)') TRIM(prefix),"thd_s"
        WRITE(long_name,'(A30,A4,A18)') "Surface down thermal radiation ", meaning, &
                                      &" since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 5, 3, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%athd_s,                    &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
          & ldims=shape2d,                                                   &
          & isteptype=a_steptype, in_group=groups("rad_vars"))


        ! &       diag%athu_s(nproma,nblks)
        WRITE(name,'(A,A5)') TRIM(prefix),"thu_s"
        WRITE(long_name,'(A30,A4,A18)') "Surface up thermal radiation ", meaning, &
                                      &" since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 5, 4, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%athu_s,                    &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
          & ldims=shape2d,                                                   &
          & isteptype=a_steptype, in_group=groups("rad_vars"))


        ! &       diag%asodird_s(nproma,nblks)
        WRITE(name,'(A,A8)') TRIM(prefix),"sodird_s"
        WRITE(long_name,'(A30,A4,A18)') "Surface down solar direct rad. ", meaning, &
                                      &" since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 4, 198, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%asodird_s,                 &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
          & ldims=shape2d,                                                   &
          & isteptype=a_steptype, in_group=groups("rad_vars"))
   
   
        ! &       diag%asodifd_s(nproma,nblks)
        WRITE(name,'(A,A8)') TRIM(prefix),"sodifd_s"
        WRITE(long_name,'(A30,A4,A18)') "Surface down solar diff. rad. ", meaning, &
                                      &" since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 4, 199, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%asodifd_s,                 &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
          & ldims=shape2d,                                                   &
          & isteptype=a_steptype, in_group=groups("rad_vars"))
   
   
        ! &       diag%asodifu_s(nproma,nblks)
        WRITE(name,'(A,A8)') TRIM(prefix),"sodifu_s"
        WRITE(long_name,'(A30,A4,A18)') "Surface up solar diff. rad. ", meaning, &
                                      &" since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 4, 8, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%asodifu_s,                 &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
          & ldims=shape2d,                                                   &
          & isteptype=a_steptype, in_group=groups("rad_vars"))


        ! &      diag%aswflx_par_sfc(nproma,nblks_c)
        WRITE(name,'(A,A13)') TRIM(prefix),"swflx_par_sfc"
        WRITE(long_name,'(A30,A4,A18)') "Downward PAR flux ", meaning, &
                                      &" since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), &
          &                   DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 4, 10, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%aswflx_par_sfc,            &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
          & ldims=shape2d,                                                   & 
          & isteptype=a_steptype, in_group=groups("rad_vars") )


        ! &      diag%vio3(nproma,nblks_c)
        cf_desc    = t_cf_var('vio3', '', 'vertically integrated ozone amount', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 14, 0, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'vio3', diag%vio3,                           &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d, lrestart=.FALSE. )


        ! &      diag%hmo3(nproma,nblks_c)
        cf_desc    = t_cf_var('hmo3', 'Pa', 'height of O3 maximum (Pa)', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'hmo3', diag%hmo3,                           &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d, lrestart=.FALSE. )


        IF (irad_aero == 5) THEN ! Old Tanre aerosol climatology taken over from the COSMO model (to be used with inwp_radiation==2)
          ! &      diag%aersea(nproma,nblks_c)
          cf_desc    = t_cf_var('aersea', '', '', DATATYPE_FLT32)
          grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
          CALL add_var( diag_list, 'aersea', diag%aersea,                       &
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
            & ldims=shape2d, lrestart=.FALSE. ) 

        
          ! &      diag%aerlan(nproma,nblks_c)
          cf_desc    = t_cf_var('aerlan', '', '', DATATYPE_FLT32)
          grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
          CALL add_var( diag_list, 'aerlan', diag%aerlan,                       &
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
            & ldims=shape2d, lrestart=.FALSE. ) 

        
          ! &      diag%aerurb(nproma,nblks_c)
          cf_desc    = t_cf_var('aerurb', '', '', DATATYPE_FLT32)
          grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
          CALL add_var( diag_list, 'aerurb', diag%aerurb,                       &
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
            & ldims=shape2d, lrestart=.FALSE. ) 

        
          ! &      diag%aerdes(nproma,nblks_c)
          cf_desc    = t_cf_var('aerdes', '', '', DATATYPE_FLT32)
          grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
          CALL add_var( diag_list, 'aerdes', diag%aerdes,                       &
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
            & ldims=shape2d, lrestart=.FALSE. ) 

        ELSE IF (irad_aero == 6) THEN ! Tegen aerosol climatology, time-interpolated values 
                                      ! (needed as state fields for coupling with microphysics and convection)

          IF (atm_phy_nwp_config(k_jg)%icpl_aero_gscp > 1 .OR. icpl_aero_conv > 1) THEN
            lrestart = .TRUE.
          ELSE
            lrestart = .FALSE.
          ENDIF

          ! &      diag%aer_ss(nproma,nblks_c)
          cf_desc    = t_cf_var('aer_ss', '', '', DATATYPE_FLT32)
          grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
          CALL add_var( diag_list, 'aer_ss', diag%aer_ss,                       &
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
            & ldims=shape2d, lrestart=lrestart ) 

          ! &      diag%aer_or(nproma,nblks_c)
          cf_desc    = t_cf_var('aer_or', '', '', DATATYPE_FLT32)
          grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
          CALL add_var( diag_list, 'aer_or', diag%aer_or,                       &
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
            & ldims=shape2d, lrestart=lrestart ) 

          ! &      diag%aer_bc(nproma,nblks_c)
          cf_desc    = t_cf_var('aer_bc', '', '', DATATYPE_FLT32)
          grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
          CALL add_var( diag_list, 'aer_bc', diag%aer_bc,                       &
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
            & ldims=shape2d, lrestart=lrestart ) 

          ! &      diag%aer_su(nproma,nblks_c)
          cf_desc    = t_cf_var('aer_su', '', '', DATATYPE_FLT32)
          grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
          CALL add_var( diag_list, 'aer_su', diag%aer_su,                       &
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
            & ldims=shape2d, lrestart=lrestart ) 

          ! &      diag%aer_du(nproma,nblks_c)
          cf_desc    = t_cf_var('aer_du', '', '', DATATYPE_FLT32)
          grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
          CALL add_var( diag_list, 'aer_du', diag%aer_du,                       &
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
            & ldims=shape2d, lrestart=lrestart ) 

        ENDIF

        IF (irad_aero == 6 .AND. (atm_phy_nwp_config(k_jg)%icpl_aero_gscp == 1 .OR. icpl_aero_conv == 1)) THEN
          lrestart = .TRUE.
        ELSE
          lrestart = .FALSE.
        ENDIF

        ! &      diag%cloud_num(nproma,nblks_c)
        cf_desc    = t_cf_var('cloud_num', 'm-3', 'cloud droplet number concentration', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'cloud_num', diag%cloud_num,                 &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d, lrestart=lrestart )

        !------------------
        !Radiation 3D variables

        ! &      diag%lwflxall(nproma,nlevp1,nblks_c)
        cf_desc    = t_cf_var('lwflxall', 'W m-2 ', 'longwave net flux', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 5, 5, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'lwflxall', diag%lwflxall,                   &
          & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,      &
          & ldims=shape3dkp1 )


        ! &      diag%trsolall(nproma,nlevp1,nblks_c)
        cf_desc    = t_cf_var('trsolall', '', 'shortwave net tranmissivity', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'trsolall', diag%trsolall,                             &
          & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc, ldims=shape3dkp1)


        !------------------
        !Turbulence 2D variables
        
        ! &      diag%shfl_s(nproma,nblks_c)
        cf_desc    = t_cf_var('shfl_s', 'W m-2 ', 'surface sensible heat flux', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 11, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'shfl_s', diag%shfl_s,                       &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d,                                                    &
          & in_group=groups("pbl_vars") )


        WRITE(name,'(A,A6)') TRIM(prefix),"shfl_s"
        WRITE(long_name,'(A27,A4,A18)') "surface sensible heat flux ", meaning, &
                                      & " since model start"
        cf_desc    = t_cf_var(TRIM(name),  TRIM(varunits), TRIM(long_name), DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 11, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%ashfl_s,                    &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d,                                                    &
          & isteptype=a_steptype, in_group=groups("pbl_vars") )


        ! &      diag%lhfl_s(nproma,nblks_c)
        cf_desc    = t_cf_var('lhfl_s', 'W m-2 ', 'surface latent heat flux', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 10, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'lhfl_s', diag%lhfl_s,                       &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d,                                                    &
          & in_group=groups("pbl_vars"))

                    
        WRITE(name,'(A,A6)') TRIM(prefix),"lhfl_s"
        WRITE(long_name,'(A27,A4,A18)') "surface latent heat flux ", meaning, &
                                      & " since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 10, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%alhfl_s,                    &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d,                                                    &
          & isteptype=a_steptype, in_group=groups("pbl_vars") )



        ! &      diag%lhfl_bs(nproma,nblks_c)
        cf_desc    = t_cf_var('lhfl_bs', 'W m-2 ', 'latent heat flux from bare soil', &
          &          DATATYPE_FLT32)
        grib2_desc = t_grib2_var(2, 0, 193, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'lhfl_bs', diag%lhfl_bs,                     &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d, lrestart=.FALSE.,                                  &
          & in_group=groups("pbl_vars"))


                    
        WRITE(name,'(A,A7)') TRIM(prefix),"lhfl_bs"
        WRITE(long_name,'(A27,A4,A18)') "latent heat flux from bare soil", meaning, &
                                      & " since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), DATATYPE_FLT32)
        grib2_desc = t_grib2_var(2, 0, 193, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%alhfl_bs,                   &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d, lrestart=.TRUE.,                                   &
          & isteptype=a_steptype )



        ! &      diag%lhfl_pl(nproma,nlev_soil,nblks_c)
        cf_desc    = t_cf_var('lhfl_pl', 'W m-2 ', 'latent heat flux from plants', &
          &          DATATYPE_FLT32)
        grib2_desc = t_grib2_var(2, 0, 194, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'lhfl_pl', diag%lhfl_pl,                     &
          & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND, cf_desc, grib2_desc, &
          & ldims=(/nproma,nlev_soil,kblks/), lrestart=.FALSE.,               &
          & in_group=groups("pbl_vars"))


                    
        WRITE(name,'(A,A7)') TRIM(prefix),"lhfl_pl"
        WRITE(long_name,'(A27,A4,A18)') "latent heat flux from plants", meaning, &
                                      & " since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), TRIM(long_name), DATATYPE_FLT32)
        grib2_desc = t_grib2_var(2, 0, 194, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%alhfl_pl,                   &
          & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND, cf_desc, grib2_desc, &
          & ldims=(/nproma,nlev_soil,kblks/), lrestart=.TRUE.,                &
          & isteptype=a_steptype )


        ! &      diag%qhfl_s(nproma,nblks_c)
        cf_desc    = t_cf_var('qhfl_s', 'Kg m-2 s-1', 'surface moisture flux', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(2, 0, 6, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'qhfl_s', diag%qhfl_s,                       &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d, lrestart=.TRUE.,                                   &
          & in_group=groups("pbl_vars"))


        WRITE(name,'(A,A6)') TRIM(prefix),"qhfl_s"
        WRITE(long_name,'(A23,A4,A18)') "surface moisture flux ", meaning, &
                                      & " since model start"
        IF (lflux_avg ) THEN
          cf_desc = t_cf_var(TRIM(name), 'Kg m-2 s-1', TRIM(long_name), DATATYPE_FLT32)
        ELSE
          cf_desc = t_cf_var(TRIM(name), 'Kg m-2', TRIM(long_name), DATATYPE_FLT32)
        ENDIF
        grib2_desc = t_grib2_var(2, 0, 6, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%aqhfl_s  ,                  &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d,                                                    &
          & isteptype=a_steptype, in_group=groups("pbl_vars") )


        ! &      diag%tcm(nproma,nblks_c)
        cf_desc    = t_cf_var('tcm', ' ','turbulent transfer coefficients for momentum', &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 2, 29, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tcm', diag%tcm,                             &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d,                                                    &
          & in_group=groups("pbl_vars") )


        ! &      diag%tch(nproma,nblks_c)
        cf_desc    = t_cf_var('tch', ' ','turbulent transfer coefficients for heat', &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 19, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tch', diag%tch,                             &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d,                                                    &
          & in_group=groups("pbl_vars") )

        
        ! &      diag%tfm(nproma,nblks_c)
        cf_desc    = t_cf_var('tfm', ' ','factor of laminar transfer of momentum', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tfm', diag%tfm,                             &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d )


        ! &      diag%tfh(nproma,nblks_c)
        cf_desc    = t_cf_var('tfh', ' ','factor of laminar transfer of scalars', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tfh', diag%tfh,                             &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d )


        ! &      diag%tfv(nproma,nblks_c)
        cf_desc    = t_cf_var('tfv', ' ','laminar reduction factor for evaporation', &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tfv', diag%tfv,                             &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d )


        ! &      diag%gz0(nproma,nblks_c)
        cf_desc     = t_cf_var('gz0', 'm2 s-2 ','roughness length times gravity', DATATYPE_FLT32)
        new_cf_desc = t_cf_var( 'z0',       'm','roughness length',               DATATYPE_FLT32)
        grib2_desc = t_grib2_var(2, 0, 1, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'gz0', diag%gz0,                             &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
          & ldims=shape2d,                                                    &
          & post_op=post_op(POST_OP_SCALE, arg1=1._wp/grav,                   &
          &                 new_cf=new_cf_desc),                              &
          & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in"), &
          & initval_r=0.01_wp )


        ! &      diag%t_2m(nproma,nblks_c)
        cf_desc    = t_cf_var('t_2m', 'K ','temperature in 2m', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 0, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 't_2m', diag%t_2m,                           &
          & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
          & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars","dwd_fg_atm_vars") )


        ! &      diag%tmax_2m(nproma,nblks_c)
        cf_desc    = t_cf_var('tmax_2m', 'K ','Max 2m temperature', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 0, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tmax_2m', diag%tmax_2m,                     &
          & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
          & ldims=shape2d, lrestart=.TRUE.,                                   &
          & isteptype=TSTEP_MAX, initval_r=-999._wp, resetval_r=-999._wp,     &
          & action_list=actions(new_action(ACTION_RESET,'PT03H')),            &
          & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
          &                                       fallback_type=HINTP_TYPE_LONLAT_RBF) ) 

        ! &      diag%tmin_2m(nproma,nblks_c)
        cf_desc    = t_cf_var('tmin_2m', 'K ','Min 2m temperature', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 0, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'tmin_2m', diag%tmin_2m,                     &
          & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
          & ldims=shape2d, lrestart=.TRUE.,                                   &
          & isteptype=TSTEP_MIN, initval_r=999._wp, resetval_r=999._wp,       &
          & action_list=actions(new_action(ACTION_RESET,'PT03H')),            &
          & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
          &                                       fallback_type=HINTP_TYPE_LONLAT_RBF) )


        ! &      diag%qv_2m(nproma,nblks_c)
        cf_desc    = t_cf_var('qv_2m', 'kg kg-1 ','specific water vapor content in 2m', &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 1, 0, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'qv_2m', diag%qv_2m,                         &
          & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
          & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars") )


        ! &      diag%rh_2m(nproma,nblks_c)
        cf_desc    = t_cf_var('rh_2m', '%','relative humidity in 2m', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 1, 1, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'rh_2m', diag%rh_2m,                         &
          & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
          & ldims=shape2d, lrestart=.FALSE. )


        ! &      diag%td_2m(nproma,nblks_c)
        cf_desc    = t_cf_var('td_2m', 'K ','dew-point in 2m', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 6, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'td_2m', diag%td_2m,                         &
          & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
          & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars","dwd_fg_atm_vars") )


        ! &      diag%u_10m(nproma,nblks_c)
        cf_desc    = t_cf_var('u_10m', 'm s-1 ','zonal wind in 10m', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 2, 2, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'u_10m', diag%u_10m,                         &
          & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars","dwd_fg_atm_vars") )


        ! &      diag%v_10m(nproma,nblks_c)
        cf_desc    = t_cf_var('v_10m', 'm s-1 ','meridional wind in 10m', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 2, 3, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'v_10m', diag%v_10m,                         &
          & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars","dwd_fg_atm_vars") )


        ! &      diag%sp_10m(nproma,nblks_c)
        cf_desc    = t_cf_var('sp_10m', 'm s-1 ','wind speed in 10m', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 2, 1, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'sp_10m', diag%sp_10m,                       &
          & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
          & ldims=shape2d, lrestart=.FALSE. )

!tiled quantities
        ! &      diag%shfl_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
        cf_desc    = t_cf_var('shfl_s_t', 'W m-2 ', 'tile-based surface sensible heat flux', &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 11, ibits, GRID_REFERENCE, GRID_CELL)
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
             & t_cf_var('shfl_s_t_'//TRIM(csfc), '', '', DATATYPE_FLT32), &
             & t_grib2_var(0, 0, 11, ibits, GRID_REFERENCE, GRID_CELL),   &
             & var_class=CLASS_TILE,                                      &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%lhfl_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
        cf_desc    = t_cf_var('lhfl_s_t', 'W m-2 ', 'tile-based surface latent heat flux', &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 10, ibits, GRID_REFERENCE, GRID_CELL)
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
             & t_cf_var('lhfl_s_t_'//TRIM(csfc), '', '', DATATYPE_FLT32), &
             & t_grib2_var(0, 0, 10, ibits, GRID_REFERENCE, GRID_CELL),   & 
             & var_class=CLASS_TILE,                                      &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO


        ! &      diag%lhfl_bs_t(nproma,nblks_c,ntiles_total)
        cf_desc    = t_cf_var('lhfl_bs_t', 'W m-2 ', 'tile-based latent heat flux from bare soil', &
          &                   DATATYPE_FLT32)
        grib2_desc = t_grib2_var(2, 0, 193, ibits, GRID_REFERENCE, GRID_CELL)
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
             & t_cf_var('lhfl_bs_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),&
             & t_grib2_var(2, 0, 193, ibits, GRID_REFERENCE, GRID_CELL),  &
             & var_class=CLASS_TILE,                                      &
             & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.)
        ENDDO


        ! &      diag%lhfl_pl_t(nproma,nlev_soil,nblks_c,ntiles_total)
        cf_desc    = t_cf_var('lhfl_pl_t', 'W m-2 ', 'tile-based latent heat flux from plants', &
          &                   DATATYPE_FLT32)
        grib2_desc = t_grib2_var(2, 0, 194, ibits, GRID_REFERENCE, GRID_CELL)
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
             & t_cf_var('lhfl_pl_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),&
             & t_grib2_var(2, 0, 194, ibits, GRID_REFERENCE, GRID_CELL),  & 
             & ldims=(/nproma,nlev_soil,kblks/), lrestart=.FALSE.,        &
             & var_class=CLASS_TILE,                                      &
             & loutput=.TRUE.)
        ENDDO

        ! &      diag%qhfl_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
        cf_desc    = t_cf_var('qhfl_s_t', 'Kg m-2 s-1','tile based surface moisture flux', &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(2, 0, 6, ibits, GRID_REFERENCE, GRID_CELL)
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
             & t_cf_var('qhfl_s_t_'//TRIM(csfc), '', '', DATATYPE_FLT32), &
             & t_grib2_var(0, 0, 10, ibits, GRID_REFERENCE, GRID_CELL),   &
             & var_class=CLASS_TILE,                                      &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO


        ! &      diag%tcm_t(nproma,nblks_c,ntiles_total+ntiles_water)
        cf_desc    = t_cf_var('tcm_t', ' ', &
          & 'tile-based turbulent transfer coefficients for momentum', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 2, 29, ibits, GRID_REFERENCE, GRID_CELL)
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
             & t_cf_var('tcm_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),    &
             & t_grib2_var(0, 2, 29, ibits, GRID_REFERENCE, GRID_CELL),   &
             & var_class=CLASS_TILE,                                      &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%tch_t(nproma,nblks_c,ntiles_total+ntiles_water)
        cf_desc    = t_cf_var('tch_t', ' ', &
             &                'tile-based turbulent transfer coefficients for heat', &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 19, ibits, GRID_REFERENCE, GRID_CELL)
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
             & t_cf_var('tch_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),    &
             & t_grib2_var(0, 0, 19, ibits, GRID_REFERENCE, GRID_CELL),   &
             & var_class=CLASS_TILE,                                      &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO


        ! &      diag%tfv_t(nproma,nblks_c,ntiles_total+ntiles_water)
        cf_desc    = t_cf_var('tfv_t', ' ', &
             &                'tile-based laminar reduction factor for evaporation', &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 11, ibits, GRID_REFERENCE, GRID_CELL)
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
             & t_cf_var('tfv_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),    &
             & t_grib2_var(0, 4, 0, ibits, GRID_REFERENCE, GRID_CELL),    &
             & var_class=CLASS_TILE,                                      &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO


        ! &      diag%gz0_t(nproma,nblks_c,ntiles_total+ntiles_water)
        cf_desc    = t_cf_var('gz0_t', 'm2 s-2 ', 'tile-based roughness length times gravity', &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(2, 0, 1, ibits, GRID_REFERENCE, GRID_CELL)
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
             & t_cf_var('gz0_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),    &
             & t_grib2_var(2, 0, 1, ibits, GRID_REFERENCE, GRID_CELL),    &
             & var_class=CLASS_TILE,                                      &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,            &
             & in_group=groups("dwd_fg_sfc_vars_t") )
        ENDDO


        ! &      diag%tvs_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
        cf_desc    = t_cf_var('tvs_s_t', 'm2 s-2 ',                              &
             &                'tile-based turbulence velocity scale at surface', &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
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
             & t_cf_var('tvs_s_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),     &
             & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL), &
             & var_class=CLASS_TILE,                                         &
             & ldims=shape2d, lrestart=.TRUE., loutput=.FALSE.)
        ENDDO



        ! &      diag%tkvm_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
        cf_desc    = t_cf_var('tkvm_s_t', 'm s-2 ',                                       &
             &                'tile-based turbulent diff. coeff for momentum at surface', &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 2, 31, ibits, GRID_REFERENCE, GRID_CELL)
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
             & t_cf_var('tkvm_s_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),     &
             & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),  &
             & var_class=CLASS_TILE,                                        &
             & ldims=shape2d, lrestart=.TRUE., loutput=.FALSE.)
        ENDDO


        ! &      diag%tkvh_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
        cf_desc    = t_cf_var('tkvh_s_t', 'm s-2 ',                                      &
             &                'tile-based turbulent diff. coeff for heat at surface',    &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 20, ibits, GRID_REFERENCE, GRID_CELL)
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
             & t_cf_var('tkvh_s_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),     &
             & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),  &
             & var_class=CLASS_TILE,                                          &
             & ldims=shape2d, lrestart=.TRUE., loutput=.FALSE.)
        ENDDO


        ! &      diag%u_10m_t(nproma,nblks_c,ntiles_total+ntiles_water)
        cf_desc    = t_cf_var('u_10m_t', 'm s-1 ', 'tile-based zonal wind in 2m', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 11, ibits, GRID_REFERENCE, GRID_CELL)
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
             & t_cf_var('u_10m_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),  &
             & t_grib2_var(0, 4, 0, ibits, GRID_REFERENCE, GRID_CELL),    &
             & var_class=CLASS_TILE,                                      &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO

        ! &      diag%v_10m_t(nproma,nblks_c,ntiles_total+ntiles_water)
        cf_desc    = t_cf_var('v_10m_t', 'm s-1 ', 'tile-based meridional wind in 2m', &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 0, 11, ibits, GRID_REFERENCE, GRID_CELL)
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
             & t_cf_var('v_10m_t_'//TRIM(csfc), '', '', DATATYPE_FLT32),  &
             & t_grib2_var(0, 4, 0, ibits, GRID_REFERENCE, GRID_CELL),    &
             & var_class=CLASS_TILE,                                      &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.)
        ENDDO


        ! &      diag%umfl_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
        cf_desc    = t_cf_var('umfl_s_t', 'N m-2 ', 'u-momentum flux at the surface', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 2, 17, ibits, GRID_REFERENCE, GRID_CELL)
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
             & t_cf_var('umfl_s_t_'//TRIM(csfc), '', '', DATATYPE_FLT32), &
             & t_grib2_var(0, 2, 17, ibits, GRID_REFERENCE, GRID_CELL),   &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,           &
             & isteptype=TSTEP_INSTANT )
        ENDDO
        !EDMF requires lrestart=.TRUE. 


        ! &      diag%vmfl_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
        cf_desc    = t_cf_var('vmfl_s_t', 'N m-2 ', 'v-momentum flux at the surface', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 2, 18, ibits, GRID_REFERENCE, GRID_CELL)
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
             & t_cf_var('vmfl_s_t_'//TRIM(csfc), '', '', DATATYPE_FLT32), &
             & t_grib2_var(0, 2, 18, ibits, GRID_REFERENCE, GRID_CELL),   &
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,           &
             & isteptype=TSTEP_INSTANT )
        ENDDO
        !EDMF requires lrestart=.TRUE. 


        ! &      diag%umfl_s(nproma,nblks_c)
        cf_desc    = t_cf_var('umfl_s', 'N m-2', 'u-momentum flux at the surface', &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 2, 17, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, 'umfl_s', diag%umfl_s,                            &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
          & lrestart=.FALSE., loutput=.TRUE.,                                      &
          & isteptype=TSTEP_INSTANT )

        ! &      diag%vmfl_s(nproma,nblks_c)
        cf_desc    = t_cf_var('vmfl_s', 'N m-2', 'v-momentum flux at the surface', &
             &                DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 2, 18, ibits, GRID_REFERENCE, GRID_CELL)
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
          &                          TRIM(long_name), DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 2, 17, ibits, GRID_REFERENCE, GRID_CELL)
!       aumfl_s and avmfl_s are needed for the restart only to get the same output values
!       They are not important to obtain bit identical model results
        CALL add_var( diag_list, TRIM(name), diag%aumfl_s,                         &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
          & isteptype=a_steptype, lrestart=.TRUE., loutput=.TRUE. )
 
        WRITE(name,'(A,A6)') TRIM(prefix),"vmfl_s"
        WRITE(long_name,'(A26,A4,A18)') "v-momentum flux flux at surface ", meaning, &
                                      & " since model start"
        cf_desc    = t_cf_var(TRIM(name), TRIM(varunits), &
          &                          TRIM(long_name), DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 2, 18, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( diag_list, TRIM(name), diag%avmfl_s,                         &
          & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
          & isteptype=a_steptype, lrestart=.TRUE., loutput=.TRUE. )


  !
  ! EDMF
  !

    IF( atm_phy_nwp_config(k_jg)%inwp_turb == iedmf) THEN

       ! &      diag%z0m(nproma,nblks_c)
       cf_desc    = t_cf_var('z0m', '', &
            &                'geopotential of the top of the atmospheric boundary layer', &
            &                DATATYPE_FLT32)
       grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
       CALL add_var( diag_list, 'z0m', diag%z0m,                             &
         & GRID_UNSTRUCTURED_CELL,ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d )

    ENDIF  !inwp_turb == EDMF

  !
  ! LES
  !

    !Anurag Dipankar, MPI (7 Oct 2013)
    !Diagnostics for LES physics
    IF ( atm_phy_nwp_config(k_jg)%is_les_phy ) THEN

      ! &      diag%z_pbl(nproma,nblks_c)
      cf_desc    = t_cf_var('z_pbl', 'm', 'boundary layer height', &
           &                DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( diag_list, 'z_pbl', diag%z_pbl,                             &
        & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,            &
        & ldims=shape2d, lrestart=.FALSE. )

      ! &      diag%bruvais(nproma,nlev+1,nblks_c)
      cf_desc    = t_cf_var('bruvais', '1/s**2', 'Brunt Vaisala Frequency', &
           &                DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( diag_list, 'bruvais', diag%bruvais,                     &
        & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,      &
        & ldims=shape3dkp1, lrestart=.FALSE. )                                   

      ! &      diag%mech_prod(nproma,nlev+1,nblks_c)
      cf_desc    = t_cf_var('mech_prod', 'm**2/s**3', 'mechanical production term in TKE Eq', &
           &                DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( diag_list, 'mech_prod', diag%mech_prod,                     &
        & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,          &
        & ldims=shape3dkp1, lrestart=.FALSE. )                                   

      !1D and 0D diagnostic variables that can not be part of add_var
      ALLOCATE( diag%turb_diag_1dvar(klevp1,SIZE(turb_profile_list,1)),  &
                diag%turb_diag_0dvar(SIZE(turb_tseries_list,1)), STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:new_nwp_phy_diag_list', &
                    'allocation of 1D and 0D diag var failed!')
      ENDIF
      diag%turb_diag_1dvar = 0._wp
      diag%turb_diag_0dvar = 0._wp

    END IF  


    !------------------
    !Turbulence 3D variables

   ! &      diag%tkvm(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('tkvm', 'm**2/s', ' turbulent diffusion coefficients for momentum', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 31, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tkvm', diag%tkvm,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,        &
      & ldims=shape3dkp1, in_group=groups("pbl_vars") )

   ! &      diag%tkvh(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('tkvh', 'm**2/s', ' turbulent diffusion coefficients for heat', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 0, 20, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tkvh', diag%tkvh,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,        &
      & ldims=shape3dkp1, in_group=groups("pbl_vars") ) 

   ! &      diag%rcld(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('rcld', '', 'standard deviation of the saturation deficit', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'rcld', diag%rcld,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,        &
      & ldims=shape3dkp1 )

! SO FAR UNUSED
!!$   ! &      diag%edr(nproma,nlevp1,nblks_c)
!!$    cf_desc    = t_cf_var('edr', '', 'eddy dissipation rate', DATATYPE_FLT32)
!!$    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
!!$    CALL add_var( diag_list, 'edr', diag%edr,                               &
!!$      & GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,        &
!!$      & ldims=shape3dkp1, lrestart=.FALSE. ) 


    !------------------
    ! Optional computation of diagnostic fields
    !------------------

    ! &     relative humidity
    !
    ! Note: This task is registered for the post-processing scheduler
    !       which takes care of the regular update:
    ! 
    IF (l_rh) THEN
      cf_desc    = t_cf_var('rh', '%', 'relative humidity', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(0, 1, 1, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( diag_list,                                                       &
                    & "rh", diag%rh,                                                 &
                    & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                             &
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


    !  Height of 0 deg C level
    cf_desc    = t_cf_var('hzerocl', '', 'height_of_0_deg_C_level', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 6, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'hzerocl', diag%hzerocl,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_ISOTHERM_ZERO, cf_desc, grib2_desc,        &
      & ldims=shape2d, lrestart=.FALSE.)


    !  significant weather WW
    cf_desc    = t_cf_var('ww', '', 'significant_weather', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 19, 25, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'ww', diag%iww,                                         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,           &
                & ldims=shape2d, lrestart=.FALSE.,                                   &
                & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB) )
      
    ! mask field to distinguish between tropics and extratropics (for tuning purposes)
    cf_desc    = t_cf_var('tropics_mask', '', 'tropics_mask', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( diag_list, 'tropics_mask', diag%tropics_mask,                      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,           &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE. )

    ! buffer field needed for the combination of vertical nesting with a reduced radiation grid
    IF (k_jg > n_dom_start) THEN
      ALLOCATE(diag%buffer_rrg(nproma, 3*nexlevs_rrg_vnest, p_patch_local_parent(k_jg)%nblks_c))
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

    INTEGER :: shape3d(3), shape3dkp1(3), shape4d(4)
    INTEGER :: ibits, ktracer, ist
    LOGICAL :: lrestart

    ibits = DATATYPE_PACK16 ! "entropy" of horizontal slice

    shape3d    = (/nproma, klev  , kblks            /)
    shape3dkp1 = (/nproma, klev+1, kblks            /)

    IF (lart) THEN
     shape4d    = (/nproma, klev  , kblks, nqtendphy+nart_tendphy /)
    ELSE
     shape4d    = (/nproma, klev  , kblks, nqtendphy /)
    ENDIF 

    CALL new_var_list( phy_tend_list, TRIM(listname), patch_id=k_jg )
    CALL default_var_list_settings( phy_tend_list,             &
                                  & lrestart=.TRUE.  )

    
    !------------------------------
    ! Temperature tendencies
    !------------------------------

   ! &      phy_tend%ddt_temp_radsw(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_temp_radsw', 'K s-1', &
         &                            'short wave radiative temperature tendency', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 4, 192, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_radsw', phy_tend%ddt_temp_radsw,      &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,        &
                & ldims=shape3d, in_group=groups("phys_tendencies") )

   ! &      phy_tend%ddt_temp_radlw(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_temp_radlw', 'K s-1', &
         &                            'long wave radiative temperature tendency', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 5, 192, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_radlw', phy_tend%ddt_temp_radlw,      &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,        &
                & ldims=shape3d, in_group=groups("phys_tendencies") )

   ! &      phy_tend%ddt_temp_turb(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_temp_turb', 'K s-1', &
         &                            'turbulence temperature tendency', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(192, 162, 121, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_turb', phy_tend%ddt_temp_turb,        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,        &
                & ldims=shape3d, lrestart=.FALSE., in_group=groups("phys_tendencies"))

   ! &      phy_tend%ddt_temp_drag(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_temp_drag', 'K s-1', &
         &                'sso + gwdrag temperature tendency', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(192, 162, 125, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_drag', phy_tend%ddt_temp_drag,        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,        &
                & ldims=shape3d, lrestart=.FALSE., in_group=groups("phys_tendencies"))

   ! &      phy_tend%ddt_temp_pconv(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_temp_pconv', 'K s-1', &
         &                            'convective temperature tendency', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 0, 192, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_pconv', phy_tend%ddt_temp_pconv,        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )


    !------------------------------
    ! Zonal Wind tendencies
    !------------------------------

   ! &      phy_tend%ddt_u_turb(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_u_turb', 'm s-2', &
         &                            'turbulence tendency of zonal wind', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(192, 162, 119, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_turb', phy_tend%ddt_u_turb,          &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,    &
                & ldims=shape3d, lrestart=.FALSE., in_group=groups("phys_tendencies"))

   ! &      phy_tend%ddt_u_sso(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_u_sso', 'm s-2', &
         &                            'sso tendency of zonal wind', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 194, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_sso', phy_tend%ddt_u_sso,            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,    &
                & ldims=shape3d, in_group=groups("phys_tendencies") )

   ! &      phy_tend%ddt_u_gwd(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_u_gwd', 'm s-2', &
         &                            'GWD tendency of zonal wind', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(192, 128, 220, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_gwd', phy_tend%ddt_u_gwd,            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,    &
                & ldims=shape3d, in_group=groups("phys_tendencies") )

   ! &      phy_tend%ddt_u_pconv(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_u_pconv', 'm s-2', &
         &                            'convective tendency of zonal wind', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 192, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_pconv', phy_tend%ddt_u_pconv,        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d)


    !------------------------------
    ! Meridional Wind tendencies
    !------------------------------

   ! &      phy_tend%ddt_v_turb(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_v_turb', 'm s-2', &
         &                            'turbulence tendency of meridional wind', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(192, 162, 120, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_turb', phy_tend%ddt_v_turb,          &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,    &
                & ldims=shape3d, lrestart=.FALSE., in_group=groups("phys_tendencies"))

   ! &      phy_tend%ddt_v_sso(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_v_sso', 'm s-2', &
         &                            'sso tendency of meridional wind', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 195, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_sso', phy_tend%ddt_v_sso,            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,    &
                & ldims=shape3d, in_group=groups("phys_tendencies") )

   ! &      phy_tend%ddt_v_gwd(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_v_gwd', 'm s-2', &
         &                            'GWD tendency of meridional wind', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(192, 128, 221, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_gwd', phy_tend%ddt_v_gwd,            &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,    &
                & ldims=shape3d, in_group=groups("phys_tendencies") )

   ! &      phy_tend%ddt_v_pconv(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_v_pconv', 'm s-2', &
         &                            'convective tendency of meridional wind', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 193, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_pconv', phy_tend%ddt_v_pconv,                &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape3d )

    !------------------------------
    ! Vertical Wind tendencies
    !------------------------------

    IF ( atm_phy_nwp_config(k_jg)%is_les_phy ) THEN
     ! &      phy_tend%ddt_w_turb(nproma,nlev+1,nblks)
      cf_desc    = t_cf_var('ddt_w_turb', 'm s-2', &
           &                            'turbulence tendency of vertical wind', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( phy_tend_list, 'ddt_w_turb', phy_tend%ddt_w_turb,          &
                  & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,    &
                  & ldims=shape3d, lrestart=.FALSE., in_group=groups("phys_tendencies"), &
                  & initval_r=0._wp)
    END IF

    !------------------------------
    ! Moist tracer tendencies
    !------------------------------

   ! &      phy_tend%ddt_tracer_turb(nproma,nlev,nblks,nqtendphy),          &
    cf_desc    = t_cf_var('ddt_tracer_turb', 's-1', &
         &                            'turbulence tendency of tracers', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tracer_turb', phy_tend%ddt_tracer_turb,        &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape4d,&
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
                    & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                             &
                    & t_cf_var('ddt_qv_turb', 'kg kg**-1 s**-1',                     &
                    & 'turbulence tendency of specific humidity', DATATYPE_FLT32),   &
                    & t_grib2_var(192, 162, 122, ibits, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d, lrestart=.FALSE.)

         !qc
        CALL add_ref( phy_tend_list, 'ddt_tracer_turb', &
                    & 'ddt_qc_turb', phy_tend%tracer_turb_ptr(2)%p_3d,               &
                    & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                             &
                    & t_cf_var('ddt_qc_turb', 'kg kg**-1 s**-1',                     &
                    & 'turbulence tendency of specific cloud_water', DATATYPE_FLT32),&
                    & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d, lrestart=.FALSE.)

         !qi
        CALL add_ref( phy_tend_list, 'ddt_tracer_turb', &
                    & 'ddt_qi_turb', phy_tend%tracer_turb_ptr(3)%p_3d,               &
                    & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                             &
                    & t_cf_var('ddt_qi_turb', 'kg kg**-1 s**-1',                     &
                    & 'turbulence tendency of specific cloud ice', DATATYPE_FLT32),  &
                    & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d, lrestart=.FALSE.)

        ! art
        IF (lart) THEN
          CALL art_tracer_interface('turb', k_jg, kblks, phy_tend_list,  &
                    & 'ddt_', ptr_arr=phy_tend%tracer_turb_ptr,          &
                    & advconf=advection_config(k_jg), phy_tend=phy_tend, &
                    & ldims=shape3d, tlev_source=1)
        ENDIF

    cf_desc    = t_cf_var('ddt_tracer_pconv', 's-1', &
         &                            'convective tendency of tracers', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tracer_pconv', phy_tend%ddt_tracer_pconv,         &
                & GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc, ldims=shape4d,&
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    IF (lart) THEN
     ktracer=nqtendphy+nart_tendphy 
    ELSE
     ktracer=nqtendphy 
    ENDIF
    ALLOCATE( phy_tend%tracer_conv_ptr(ktracer) )

         !qv
        CALL add_ref( phy_tend_list, 'ddt_tracer_pconv', &
                    & 'ddt_qv_conv', phy_tend%tracer_conv_ptr(1)%p_3d,               &
                    & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                             &
                    & t_cf_var('ddt_qv_conv', 'kg kg**-1 s**-1',                     &
                    & 'convective tendency of specific humidity', DATATYPE_FLT32),   &
                    & t_grib2_var(0, 1, 197, ibits, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d)
         !qc
        CALL add_ref( phy_tend_list, 'ddt_tracer_pconv', &
                    & 'ddt_qc_conv', phy_tend%tracer_conv_ptr(2)%p_3d,               &
                    & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                             &
                    & t_cf_var('ddt_qc_conv', 'kg kg**-1 s**-1',                     &
                    & 'convective tendency of specific cloud water', DATATYPE_FLT32),&
                    & t_grib2_var(0, 1, 198, ibits, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d)
         !qi
        CALL add_ref( phy_tend_list, 'ddt_tracer_pconv', &
                    & 'ddt_qi_conv', phy_tend%tracer_conv_ptr(3)%p_3d,               &
                    & GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                             &
                    & t_cf_var('ddt_qi_conv', 'kg kg**-1 s**-1',                     &
                    & 'convective tendency of specific cloud ice', DATATYPE_FLT32),  &
                    & t_grib2_var(0, 1, 199, ibits, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d)

        ! art
        IF (lart) THEN
          CALL art_tracer_interface('conv', k_jg, kblks, phy_tend_list,  &
                    & 'ddt_', ptr_arr=phy_tend%tracer_conv_ptr,          &
                    & advconf=advection_config(k_jg), phy_tend=phy_tend, &
                    & ldims=shape3d, tlev_source=1)
        ENDIF


    !------------------------------
    ! TKE tendency
    !------------------------------

    !      phy_tend%ddt_tke(nproma,nlevp1,nblks)
    cf_desc    = t_cf_var('ddt_tke', 'm s-2'          , &
         &                'tendency of turbulent velocity scale', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 19, 192, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tke', phy_tend%ddt_tke,             &
                GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,&
              & ldims=shape3dkp1 )

    IF (ltkecon) THEN
      lrestart = .TRUE.
    ELSE
      lrestart = .FALSE.
    ENDIF

    !      phy_tend%ddt_tke_pconv(nproma,nlevp1,nblks)
    cf_desc    = t_cf_var('ddt_tke_pconv', 'm**2 s**-3'          , &
         &                'TKE tendency due to sub-grid scale convection', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 19, 219, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tke_pconv', phy_tend%ddt_tke_pconv,   &
                GRID_UNSTRUCTURED_CELL, ZA_HYBRID_HALF, cf_desc, grib2_desc,&
              & ldims=shape3dkp1, lrestart=lrestart )



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
 
    END IF

    CALL message('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'construction of NWP physical tendency fields finished')

END SUBROUTINE new_nwp_phy_tend_list


!
!-------------------------------------------------------------------------

END MODULE mo_nwp_phy_state
!<
