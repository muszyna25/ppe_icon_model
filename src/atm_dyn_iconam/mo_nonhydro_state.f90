#if (defined (__GNUC__) || defined(__SUNPRO_F95) || defined(__SX__))
#define HAVE_F95
#endif
!>
!! Constructs and destructs the state vector of the nonhydrostatic.
!!
!! Constructs and destructs the state vector of the nonhydrostatic
!! model variables. They are subdivided in several classes: prognostics
!! and diagnostics.
!!
!! @author Almut Gassmann (MPI-M)
!! @author Daniel Reinert (DWD-M)
!!
!! @par Revision History
!! Initial release by Almut Gassmann, MPI-M (2009-03-06)
!! Modification by Daniel Reinert, DWD (2011-05-02)
!! - Memory allocation method changed from explicit allocation to Luis'
!!   infrastructure
!! Modification by Daniel Reinert, DWD (2012-02-07)
!! - Moved type definition to new module mo_nonhydro_types, to avoid 
!!   circular dependencies
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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
!!
MODULE mo_nonhydro_state

  USE mo_kind,                 ONLY: wp
  USE mo_impl_constants,       ONLY: SUCCESS, MAX_CHAR_LENGTH, INWP,     &
    &                                VINTP_METHOD_UV, VINTP_TYPE_P_OR_Z, &
    &                                VINTP_METHOD_QV, VINTP_METHOD_PRES, &
    &                                VINTP_METHOD_LIN, VINTP_TYPE_Z,     &
    &                                VINTP_METHOD_LIN_NLEVP1,            &
    &                                TASK_INTP_MSL, HINTP_TYPE_NONE
  USE mo_exception,            ONLY: message, finish, message_text
  USE mo_model_domain,         ONLY: t_patch
  USE mo_nonhydro_types,       ONLY: t_nh_state, t_nh_prog, t_nh_diag,  &
    &                                t_nh_ref, t_nh_metrics, t_buffer_memory
  USE mo_opt_diagnostics,      ONLY: t_nh_diag_pz
  USE mo_grid_config,          ONLY: n_dom, l_limited_area, ifeedback_type
  USE mo_nonhydrostatic_config,ONLY: itime_scheme, l_nest_rcf, igradp_method, iadv_rcf
  USE mo_dynamics_config,      ONLY: nsav1, nsav2
  USE mo_parallel_config,      ONLY: nproma
  USE mo_run_config,           ONLY: iforcing, ntracer,                    &
    &                                iqv, iqc, iqi, iqr, iqs, iqt, iqtvar, &
    &                                nqtendphy, ltestcase 
  USE mo_io_config,            ONLY: lwrite_extra, inextra_2d, inextra_3d
  USE mo_nh_pzlev_config,      ONLY: nh_pzlev_config
  USE mo_advection_config,     ONLY: t_advection_config, advection_config
  USE mo_linked_list,          ONLY: t_var_list
  USE mo_var_list,             ONLY: default_var_list_settings, add_var,     &
    &                                add_ref, new_var_list, delete_var_list, &
    &                                create_tracer_metadata, add_var_list_reference, &
    &                                create_vert_interp_metadata,            &
    &                                create_hor_interp_metadata, groups
  USE mo_linked_list,          ONLY: t_list_element, find_list_element
  USE mo_var_metadata,         ONLY: t_var_metadata, t_tracer_meta
  USE mo_cf_convention,        ONLY: t_cf_var
  USE mo_grib2,                ONLY: t_grib2_var
  USE mo_cdi_constants
  USE mo_art_config,           ONLY: t_art_config,art_config
  USE mo_art_tracer_interface, ONLY: art_tracer_interface
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_art_tracer_interface, ONLY: art_tracer_interface
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


  PUBLIC :: construct_nh_state    ! Constructor for the nonhydrostatic state
  PUBLIC :: destruct_nh_state     ! Destructor for the nonhydrostatic state
  PUBLIC :: duplicate_prog_state  ! Copy the prognostic state

  PUBLIC :: p_nh_state            ! state vector of nonhydrostatic variables (variable)
  PUBLIC :: bufr
  
  TYPE (t_buffer_memory), POINTER :: bufr(:)

  TYPE(t_nh_state), TARGET, ALLOCATABLE :: p_nh_state(:)

  CONTAINS

!-------------------------------------------------------------------------
!!            SUBROUTINES FOR BUILDING AND DELETING VARIABLE LISTS 
!-------------------------------------------------------------------------
!
!
  !>
  !! Constructor for prognostic and diagnostic states.
  !!
  !! Top-level procedure for building the prognostic and diagnostic states.
  !! It calls constructors to single time level prognostic states, and 
  !! diagnostic states.
  !! Initialization of all components with zero.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-03-06)
  !!
  SUBROUTINE construct_nh_state(p_patch, p_nh_state, n_timelevels, l_pres_msl)
!
    TYPE(t_patch),     INTENT(IN)   ::  & ! patch
      &  p_patch(n_dom)
    TYPE(t_nh_state),  INTENT(INOUT)::  & ! nh state at different grid levels
      &  p_nh_state(n_dom)
    INTEGER, OPTIONAL, INTENT(IN)   ::  & ! number of timelevels
      &  n_timelevels    
    LOGICAL :: l_pres_msl(:) !< Flag. TRUE if computation of mean sea level pressure desired

    INTEGER  :: ntl,      &! local number of timelevels
                ntl_pure, &! local number of timelevels (without any extra timelevs)
                ist,      &! status
                jg,       &! grid level counter
                jt         ! time level counter

    LOGICAL  :: l_extra_timelev

    CHARACTER(len=MAX_CHAR_LENGTH) :: listname, varname_prefix

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nonhydro_state:construct_nh_state'
!-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'Construction of NH state started')

    DO jg = 1, n_dom

      IF(PRESENT(n_timelevels))THEN
        ntl      = n_timelevels
        ntl_pure = n_timelevels
      ELSE
        ntl      = 1
        ntl_pure = 1
      ENDIF

      ! If grid nesting is not called at every dynamics time step, an extra time
      ! level is needed for full-field interpolation and boundary-tendency calculation
      IF (l_nest_rcf .AND. n_dom > 1) THEN
        ntl = ntl + 1
        nsav1(jg) = ntl
      ENDIF

      ! In the presence of grid nesting and incremental feedback, another extra time level is needed to
      ! compute the feedback increments
      ! This extra time level is also used to store the driving-model data in the
      ! limited-area mode
      IF (ifeedback_type == 1 .AND. jg > 1 .OR. l_limited_area .AND. jg == 1) ntl = ntl + 1
      nsav2(jg) = ntl

      !
      ! Allocate pointer array p_nh_state(jg)%prog, as well as the 
      ! corresponding list array for each grid level.
      !
      ! create state arrays
      ALLOCATE(p_nh_state(jg)%prog(1:ntl), STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish (TRIM(routine),                               &
          &          'allocation of prognostic state array failed')
      ENDIF

      ! create state list
      ALLOCATE(p_nh_state(jg)%prog_list(1:ntl), STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish (TRIM(routine),                                    &
          &          'allocation of prognostic state list array failed')
      ENDIF

      ! create tracer list (no extra timelevels)
      ALLOCATE(p_nh_state(jg)%tracer_list(1:ntl_pure), STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish (TRIM(routine),                                    &
          &          'allocation of prognostic tracer list array failed')
      ENDIF


      !
      ! Build lists for every timelevel
      ! 
      DO jt = 1, ntl

        ! Tracer fields do not need extra time levels because feedback is not incremental
        ! and the nest-call frequency is always synchronized with the advection time step
        IF (jt > n_timelevels) THEN
          l_extra_timelev = .TRUE.
        ELSE
          l_extra_timelev = .FALSE.
        ENDIF

        WRITE(listname,'(a,i2.2,a,i2.2)') 'nh_state_prog_of_domain_',jg, &
          &                               '_and_timelev_',jt

        ! Build prog state list
        ! includes memory allocation
        !
        ! varname_prefix = 'nh_prog_'
        varname_prefix = ''
        CALL new_nh_state_prog_list(p_patch(jg), p_nh_state(jg)%prog(jt),  &
          &  p_nh_state(jg)%prog_list(jt), listname, TRIM(varname_prefix), &
          &  l_extra_timelev, jt)

        !
        ! Build prog state tracer list
        ! no memory allocation (only references to prog list)
        !
        IF (.NOT. l_extra_timelev) THEN ! not needed for extra timelevel
          WRITE(listname,'(a,i2.2,a,i2.2)') 'nh_state_tracer_of_domain_',jg, &
            &                               '_and_timelev_',jt
          varname_prefix = ''
          CALL new_nh_state_tracer_list(p_patch(jg), p_nh_state(jg)%prog_list(jt), &
            &  p_nh_state(jg)%tracer_list(jt), listname )
        ENDIF
      ENDDO ! jt

      !
      ! Build diag state list
      ! includes memory allocation
      !
      WRITE(listname,'(a,i2.2)') 'nh_state_diag_of_domain_',jg
      CALL new_nh_state_diag_list(p_patch(jg), p_nh_state(jg)%diag, &
        &  p_nh_state(jg)%diag_list, listname, l_pres_msl(jg) )

      !
      ! Build metrics state list
      ! includes memory allocation
      !
      WRITE(listname,'(a,i2.2)') 'nh_state_metrics_of_domain_',jg
      CALL new_nh_metrics_list(p_patch(jg), p_nh_state(jg)%metrics, &
        &  p_nh_state(jg)%metrics_list, listname )

      !
      ! Build ref state list (not needed so far for real case applications)
      ! includes memory allocation
      !
      IF ( ltestcase ) THEN
        WRITE(listname,'(a,i2.2)') 'nh_state_ref_of_domain_',jg
        CALL new_nh_state_ref_list(p_patch(jg), p_nh_state(jg)%ref, &
          &  p_nh_state(jg)%ref_list, listname)
      ENDIF

    ENDDO ! jg

    CALL message (TRIM(routine), 'NH state construction completed')

  END SUBROUTINE construct_nh_state


!-------------------------------------------------------------------------
!
!
  !>
  !! Destructor for prognostic and diagnostic states.
  !!
  !! It calls destructors to
  !! single time level prognostic states, and diagnostic states.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-03-06)
  !!
  SUBROUTINE destruct_nh_state(p_nh_state)
!
    TYPE(t_nh_state), INTENT(INOUT) :: & ! nh state at different grid levels
      &  p_nh_state(n_dom)
                                             
    INTEGER  :: ntl_prog, & ! number of timelevels prog state
                ntl_tra,  & ! number of timelevels 
                ist, &      ! status
                jg,  &      ! grid level counter
                jt          ! time level counter

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nonhydro_state:destruct_nh_state'
!-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'Destruction of NH state started')


    DO jg = 1, n_dom

      ntl_prog = SIZE(p_nh_state(jg)%prog(:))
      IF(ntl_prog==0)THEN
        CALL finish(TRIM(routine), 'prognostic array has no timelevels')
      ENDIF

      ntl_tra = SIZE(p_nh_state(jg)%tracer_list(:))
      IF(ntl_tra==0)THEN
        CALL finish(TRIM(routine), 'tracer list has no timelevels')
      ENDIF

      ! delete reference state list elements
      IF ( ltestcase ) THEN
        CALL delete_var_list( p_nh_state(jg)%ref_list )
      ENDIF

      ! delete diagnostic state list elements
      CALL delete_var_list( p_nh_state(jg)%diag_list )

      ! delete metrics state list elements
      CALL delete_var_list( p_nh_state(jg)%metrics_list )


      ! delete prognostic state list elements
      DO jt = 1, ntl_prog
        CALL delete_var_list( p_nh_state(jg)%prog_list(jt) )
      ENDDO

      ! delete tracer list list elements
      DO jt = 1, ntl_tra
        CALL delete_var_list( p_nh_state(jg)%tracer_list(jt) )
      ENDDO


      ! destruct state lists and arrays
      DEALLOCATE(p_nh_state(jg)%prog_list, STAT=ist )
      IF(ist/=SUCCESS) CALL finish (TRIM(routine),&
        & 'deallocation of prognostic state list array failed')

      DEALLOCATE(p_nh_state(jg)%tracer_list, STAT=ist )
      IF(ist/=SUCCESS) CALL finish (TRIM(routine),&
        & 'deallocation of tracer list array failed')

      DEALLOCATE(p_nh_state(jg)%prog )
      IF(ist/=SUCCESS) CALL finish (TRIM(routine),&
        & 'deallocation of prognostic state array failed')

    ENDDO

    CALL message (TRIM(routine), 'NH state destruction completed')

  END SUBROUTINE destruct_nh_state
  !-------------------------------------------------------------------------
  !
  !
  !>
  !!
  !! duplicate prognostic state 
  !!
  !! @par Revision History
  !! Initial release by P. Ripodas 
  !!
  SUBROUTINE duplicate_prog_state ( p_prog_i, p_prog_d)

      TYPE(t_nh_prog), INTENT(IN)      :: &  !< prognostic state vector to be copied
     &  p_prog_i
      TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< duplicated prognostic state vector
     &  p_prog_d

    !--------------------------------------------------------------
     p_prog_d%w(:,:,:)              = p_prog_i%w(:,:,:)
     p_prog_d%vn(:,:,:)             = p_prog_i%vn(:,:,:)
     p_prog_d%rho(:,:,:)            = p_prog_i%rho(:,:,:)
     p_prog_d%theta_v(:,:,:)        = p_prog_i%theta_v(:,:,:)
     IF (ASSOCIATED(p_prog_i%exner)) &
       p_prog_d%exner(:,:,:)          = p_prog_i%exner(:,:,:)
     IF (ASSOCIATED(p_prog_i%rhotheta_v)) &
       p_prog_d%rhotheta_v(:,:,:)     = p_prog_i%rhotheta_v(:,:,:)
     IF (ASSOCIATED(p_prog_i%tracer)) THEN
      p_prog_d%tracer(:,:,:,:)       = p_prog_i%tracer(:,:,:,:)
     END IF
     IF (ASSOCIATED(p_prog_i%tke )) THEN
      p_prog_d%tke(:,:,:)            = p_prog_i%tke(:,:,:)
     END IF

  END SUBROUTINE duplicate_prog_state

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Allocation of components of prognostic state.
  !!
  !! Initialization of components with zero.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-03-06)
  !!
  SUBROUTINE new_nh_state_prog_list ( p_patch, p_prog, p_prog_list,  &
    &                                 listname, vname_prefix, l_extra_timelev, timelev)
!
    TYPE(t_patch), TARGET, INTENT(IN) :: & !< current patch
      &  p_patch

    TYPE(t_nh_prog),  INTENT(INOUT)   :: & !< current prognostic state
      &  p_prog 

    TYPE(t_var_list), INTENT(INOUT)   :: p_prog_list !< current prognostic state list

    CHARACTER(len=*), INTENT(IN)      :: & !< list name
      &  listname, vname_prefix

    LOGICAL, INTENT(IN) :: l_extra_timelev  !< specifies extra time levels for which 
                                            !< not all variables are allocated

    INTEGER, INTENT(IN) :: timelev


    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
               nblks_e       !< number of edge blocks to allocate

    INTEGER :: nlev, nlevp1, ktracer

    INTEGER :: shape3d_c(3), shape3d_e(3), shape3d_chalf(3), &
      &        shape4d_c(4)

    INTEGER :: ibits         !< "entropy" of horizontal slice

    CHARACTER(len=4) suffix

    TYPE(t_advection_config), POINTER :: advconf
    TYPE(t_art_config), POINTER :: artconf

    INTEGER           :: jt
    CHARACTER(LEN=1)  :: ctracer
    CHARACTER(len=21) :: name

    !**
    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    ! pointer to advection_config(jg) to save some paperwork
    advconf => advection_config(p_patch%id)
    ! ART: pointer to art_config(jg) to save some paperwork
    artconf => art_config(p_patch%id)

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape3d_e     = (/nproma, nlev,   nblks_e  /)
    shape3d_c     = (/nproma, nlev,   nblks_c  /)
    shape3d_chalf = (/nproma, nlevp1, nblks_c  /)
    shape4d_c     = (/nproma, nlev,   nblks_c, ntracer /)

    ! Suffix (mandatory for time level dependent variables)

    WRITE(suffix,'(".TL",i1)') timelev

    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_prog_list, TRIM(listname), patch_id=p_patch%id )
    CALL default_var_list_settings( p_prog_list,               &
                                  & lrestart=.TRUE.,           &
                                  & restart_type=FILETYPE_NC2  )


    !------------------------------
    ! Meteorological quantities
    !------------------------------

    ! vn           p_prog%vn(nproma,nlev,nblks_e)
    cf_desc    = t_cf_var('normal_velocity', 'm s-1', 'velocity normal to edge', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 34, ibits, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_prog_list, TRIM(vname_prefix)//'vn'//suffix, p_prog%vn,     &
      &           GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
      &           ldims=shape3d_e, in_group=groups("nh_prog_vars") )

    ! w            p_prog%w(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('upward_air_velocity', 'm s-1', 'Vertical velocity', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 9, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_prog_list, TRIM(vname_prefix)//'w'//suffix, p_prog%w,   &
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, &
      &          ldims=shape3d_chalf,                                       & 
      &          vert_interp=create_vert_interp_metadata(                   &
      &             vert_intp_type=VINTP_TYPE_P_OR_Z,                       &
      &             vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ),             &
      &          in_group=groups("atmo_ml_vars", "atmo_pl_vars", "atmo_zl_vars") )

    ! rho          p_prog%rho(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('air_density', 'kg m-3', 'density', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 10, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_prog_list, TRIM(vname_prefix)//'rho'//suffix, p_prog%rho,     &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,      &
      &           ldims=shape3d_c, in_group=groups("nh_prog_vars") )

    ! theta_v      p_prog%theta_v(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('virtual_potential_temperature', 'K', &
      &                   'virtual potential temperature', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 0, 15, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_prog_list, TRIM(vname_prefix)//'theta_v'//suffix, p_prog%theta_v, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,          &
      &           ldims=shape3d_c, in_group=groups("nh_prog_vars") )

    ! Initialize pointers that are not always allocated to NULL
    p_prog%exner      => NULL()
    p_prog%rhotheta_v => NULL()
    p_prog%tke        => NULL()
    p_prog%tracer     => NULL()

    IF (.NOT. l_extra_timelev) THEN
      ! exner        p_prog%exner(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('exner_pressure', '-', 'exner pressure', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(0, 3, 26, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_prog_list, TRIM(vname_prefix)//'exner'//suffix, p_prog%exner, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,      &
        &           ldims=shape3d_c, in_group=groups("nh_prog_vars") )

      ! rhotheta_v   p_prog%rhotheta_v(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('rho_virt_pot_temp', 'K', 'rho virt pot temp', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(0, 19, 192, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_prog_list, TRIM(vname_prefix)//'rhotheta_v'//suffix, p_prog%rhotheta_v, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,                &
        &           ldims=shape3d_c )

      ! Tracer array for (model) internal use

      ! tracer         p_prog%tracer(nproma,nlev,nblks_c,ntracer)
      IF (ntracer > 0) THEN
        cf_desc    = t_cf_var('tracer', 'kg kg-1', 'tracer', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_prog_list, 'tracer', p_prog%tracer,                       &
          &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
          &           ldims=shape4d_c ,                                           &
          &           lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      ENDIF

      IF (  iforcing == inwp  ) THEN
        
        ! Reference to individual tracer, for I/O and setting of additional metadata
        ! Note that for qv, qc, qi, qr, qs the corresponding indices iqv, iqc, iqi, 
        ! iqr, iqs are hardcoded. For additional tracers, indices need to be set via 
        ! add_tracer_ref. 
        

        ktracer=ntracer
        ALLOCATE( p_prog%tracer_ptr(ktracer) )
        !QV
        CALL add_ref( p_prog_list, 'tracer',                                         &
                    & TRIM(vname_prefix)//'qv'//suffix, p_prog%tracer_ptr(iqv)%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var('specific_humidity', 'kg kg-1','Specific humidity',   &
                    &          DATATYPE_FLT32),  &
                    & t_grib2_var(192, 201, 28, ibits, GRID_REFERENCE, GRID_CELL),   &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=1,     &              ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                    &             ihadv_tracer=advconf%ihadv_tracer(iqv),            &
                    &             ivadv_tracer=advconf%ivadv_tracer(iqv)),           & 
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=VINTP_TYPE_P_OR_Z,                  &
                    &             vert_intp_method=VINTP_METHOD_QV,                  &
                    &             l_satlimit=.FALSE.,                                &
                    &             lower_limit=2.5e-6_wp, l_restore_pbldev=.FALSE. ), &
                    & in_group=groups("atmo_ml_vars", "atmo_pl_vars", "atmo_zl_vars") )
        !QC
        CALL add_ref( p_prog_list, 'tracer',&
                    & TRIM(vname_prefix)//'qc'//suffix, p_prog%tracer_ptr(iqc)%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var(TRIM(vname_prefix)//'qc',                             &
                    &  'kg kg-1', 'specific_cloud_water_content', DATATYPE_FLT32),   &
                    & t_grib2_var(192, 201, 31, ibits, GRID_REFERENCE, GRID_CELL),   &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=1,     &              ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                    &             ihadv_tracer=advconf%ihadv_tracer(iqc),            &
                    &             ivadv_tracer=advconf%ivadv_tracer(iqc)),           &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=VINTP_TYPE_P_OR_Z,                  &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,             &
                    &             lower_limit=0._wp  ),                              & 
                    & in_group=groups("atmo_ml_vars", "atmo_pl_vars", "atmo_zl_vars") )
        !QI
        CALL add_ref( p_prog_list, 'tracer',                                         &
                    & TRIM(vname_prefix)//'qi'//suffix, p_prog%tracer_ptr(iqi)%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var(TRIM(vname_prefix)//'qi',                             &
                    &  'kg kg-1','specific_cloud_ice_content', DATATYPE_FLT32),      &
                    & t_grib2_var(192, 201, 33, ibits, GRID_REFERENCE, GRID_CELL),   &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=1,     &              ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                    &             ihadv_tracer=advconf%ihadv_tracer(iqi),            &
                    &             ivadv_tracer=advconf%ivadv_tracer(iqi)),           &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=VINTP_TYPE_P_OR_Z,                  &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,             &
                    &             lower_limit=0._wp  ),                              & 
                    & in_group=groups("atmo_ml_vars", "atmo_pl_vars", "atmo_zl_vars")  )
        !QR
        CALL add_ref( p_prog_list, 'tracer',                                         &
                    & TRIM(vname_prefix)//'qr'//suffix, p_prog%tracer_ptr(iqr)%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var(TRIM(vname_prefix)//'qr',                             &
                    &  'kg kg-1','rain_mixing_ratio', DATATYPE_FLT32),               &
                    & t_grib2_var(0, 1, 85, ibits, GRID_REFERENCE, GRID_CELL),       &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=1,     &              ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                    &             ihadv_tracer=advconf%ihadv_tracer(iqr),            &
                    &             ivadv_tracer=advconf%ivadv_tracer(iqr)),           &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=VINTP_TYPE_P_OR_Z,                  &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,             &
                    &             lower_limit=0._wp  ),                              & 
                    & in_group=groups("atmo_ml_vars", "atmo_pl_vars", "atmo_zl_vars")  )
        !QS
        CALL add_ref( p_prog_list, 'tracer',                                         &
                    & TRIM(vname_prefix)//'qs'//suffix, p_prog%tracer_ptr(iqs)%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var(TRIM(vname_prefix)//'qs',                             &
                    &  'kg kg-1','snow_mixing_ratio', DATATYPE_FLT32),               &
                    & t_grib2_var(0, 1, 86, ibits, GRID_REFERENCE, GRID_CELL),       &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=1,     &              ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                    &             ihadv_tracer=advconf%ihadv_tracer(iqs),            &
                    &             ivadv_tracer=advconf%ivadv_tracer(iqs)),           &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=VINTP_TYPE_P_OR_Z,                  &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,             &
                    &             lower_limit=0._wp  ),                              & 
                    & in_group=groups("atmo_ml_vars", "atmo_pl_vars", "atmo_zl_vars")  )
        ! EDMF: total water variance
        IF (atm_phy_nwp_config(p_patch%id)%inwp_turb == 3) THEN
          CALL add_ref( p_prog_list, 'tracer',                                       &
                    & TRIM(vname_prefix)//'qtvar'//suffix, p_prog%tracer_ptr(iqtvar)%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var(TRIM(vname_prefix)//'qtvar',                          &
                    &  'kg2 kg-2','total water variance', DATATYPE_FLT32),           &
                    & t_grib2_var(192, 201, 39, ibits, GRID_REFERENCE, GRID_CELL),   &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=1,     &              ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(),                          &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=VINTP_TYPE_P_OR_Z,                  &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,             &
                    &             lower_limit=0._wp  )  )
        ENDIF

        ! art
        IF (artconf%lart) THEN
          CALL art_tracer_interface(p_patch,p_prog_list,vname_prefix, p_prog%tracer_ptr,&
                    & timelev,ldims=shape3d_c,tlev_source=1)
        ENDIF
   
        ! tke            p_prog%tke(nproma,nlevp1,nblks_c)
        ! for output take field from nnow_rcf slice
        cf_desc    = t_cf_var('specific_kinetic_energy_of_air', 'm2 s-2',         &
          &          'turbulent kinetic energy', DATATYPE_FLT32)
        grib2_desc = t_grib2_var(0, 19, 11, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_prog_list, TRIM(vname_prefix)//'tke'//suffix, p_prog%tke, &
          &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                       &
          &           cf_desc, grib2_desc, ldims=shape3d_chalf,                   &
          &           tlev_source=1,                                              &
          &           vert_interp=create_vert_interp_metadata(                    &
          &             vert_intp_type=VINTP_TYPE_P_OR_Z,                         &
          &             vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ),               &
          &           in_group=groups("atmo_ml_vars", "atmo_pl_vars", "atmo_zl_vars") )
        
      ELSE

        ! add refs to tracers with generic name Q1 ... Qn
        ! (used for example with test cases)
        ALLOCATE( p_prog%tracer_ptr(ntracer) )

        DO jt = 1, ntracer
          ctracer = advconf%ctracer_list(jt:jt)
          WRITE(name,'(A1,A1)') "Q", ctracer
          CALL add_ref( p_prog_list, 'tracer',                                  &
            & TRIM(name)//suffix, p_prog%tracer_ptr(jt)%p_3d,                   &
            & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                             &
            & t_cf_var(TRIM(name), 'kg/kg','Tracer mixing ratio '//TRIM(name),  &
            & DATATYPE_FLT32),                                                  &
            & t_grib2_var(0, 0, 0, ibits, GRID_REFERENCE, GRID_CELL),           &
            & ldims=shape3d_c,                                                  &
            & tlev_source=1,     &              ! output from nnow_rcf slice
            & tracer_info=create_tracer_metadata(),                             &
            & vert_interp=create_vert_interp_metadata(                          &
            &             vert_intp_type=VINTP_TYPE_P_OR_Z,                     & 
            &             vert_intp_method=VINTP_METHOD_LIN,                    &
            &             l_loglin=.FALSE.,                                     &
            &             l_extrapol=.TRUE., l_pd_limit=.FALSE.,                &
            &             lower_limit=0._wp  )  )
        END DO
      ENDIF

    ENDIF ! allocation only if not extra_timelev


  END SUBROUTINE new_nh_state_prog_list



  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Creates tracer var list.
  !!
  !! Creates tracer var list containing references to all prognostic tracer 
  !! fields.
  !!
  !! @par Revision History
  !! Initial release by Daniel Reinert, DWD (2012-02-02)
  !!
  SUBROUTINE new_nh_state_tracer_list ( p_patch, from_var_list, p_tracer_list,  &
    &                                 listname )
!
    TYPE(t_patch), TARGET, INTENT(IN) :: & !< current patch
      &  p_patch

    TYPE(t_var_list), INTENT(IN)      :: & !< source list to be referenced
      &  from_var_list 

    TYPE(t_var_list), INTENT(INOUT)   :: & !< new tracer list (containing all tracers)
      &  p_tracer_list

    CHARACTER(len=*), INTENT(IN)      :: & !< list name
      &  listname

    ! local variables
    TYPE (t_var_metadata), POINTER :: from_info
    TYPE (t_list_element), POINTER :: element
    TYPE (t_list_element), TARGET  :: start_with

    !--------------------------------------------------------------

    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_tracer_list, TRIM(listname), patch_id=p_patch%id )
    CALL default_var_list_settings( p_tracer_list,             &
                                  & lrestart=.FALSE.,          &
                                  & loutput =.FALSE.           )


    !
    ! add references to all tracer fields of the source list (prognostic state)
    !
    element => start_with
    element%next_list_element => from_var_list%p%first_list_element
    !
    for_all_list_elements: DO
      !
      element => element%next_list_element
      IF (.NOT.ASSOCIATED(element)) EXIT
      !
      ! retrieve information from actual linked list element
      !
      from_info => element%field%info

      ! Only add tracer fields to the tracer list
      IF (from_info%tracer%lis_tracer .AND. .NOT. from_info%lcontainer ) THEN

        CALL add_var_list_reference(p_tracer_list, from_info%name, &
          &                         from_var_list%p%name           )

      ENDIF

    ENDDO for_all_list_elements


  END SUBROUTINE new_nh_state_tracer_list


  !-------------------------------------------------------------------------
  !>
  !! Allocation of components of diagnostic state.
  !!
  !! Initialization of components with zero.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-03-06)
  !!
  !! Modification by Kristina Froehlich, DWD, (2010-10-22)
  !! - added pressure on interfaces
  !!
  SUBROUTINE new_nh_state_diag_list ( p_patch, p_diag, p_diag_list,  &
    &                                 listname, l_pres_msl )
!
    TYPE(t_patch), TARGET, INTENT(IN) :: &  !< current patch
      &  p_patch

    TYPE(t_nh_diag),  INTENT(INOUT)   :: &  !< diagnostic state
      &  p_diag 

    TYPE(t_var_list), INTENT(INOUT)   :: &  !< diagnostic state list
      &  p_diag_list
    CHARACTER(len=*), INTENT(IN)      :: &  !< list name
      &  listname
    LOGICAL, INTENT(IN)               :: &  !< Flag. If .TRUE., compute mean sea level pressure
      &  l_pres_msl

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
               nblks_e, &    !< number of edge blocks to allocate
               nblks_v       !< number of vertex blocks to allocate

    INTEGER :: nlev, nlevp1

    INTEGER :: n_timlevs     !< number of time levels for advection 
                             !< tendency fields

    INTEGER :: shape2d_c(2), shape2d_e(2), shape3d_c(3),           &
      &        shape3d_e(3), shape3d_v(3), shape3d_chalf(3),       &
      &        shape3d_ehalf(3), shape4d_chalf(4), shape4d_e(4),   &
      &        shape4d_entl(4), shape4d_chalfntl(4), shape4d_c(4), &
      &        shape3d_ctra(3), shape2d_extra(3), shape3d_extra(4),&
      &        shape3d_c3(3), shape3d_ubcp(3), shape3d_ubcc(3),    &
      &        shape3d_ubcp1(3)
 
    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: jt

    CHARACTER(LEN=2) :: ctrc
    CHARACTER(len=4) suffix
    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nblks_v = p_patch%nblks_v

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    IF (itime_scheme >= 2) THEN
     n_timlevs = 2
    ELSE
     n_timlevs = 1
    ENDIF

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d_c     = (/nproma,          nblks_c    /)
    shape2d_e     = (/nproma,          nblks_e    /)
    shape2d_extra = (/nproma, nblks_c, inextra_2d /)
    shape3d_c     = (/nproma, nlev   , nblks_c    /)
    shape3d_e     = (/nproma, nlev   , nblks_e    /)
    shape3d_v     = (/nproma, nlev   , nblks_v    /)
    shape3d_chalf = (/nproma, nlevp1 , nblks_c    /)
    shape3d_ehalf = (/nproma, nlevp1 , nblks_e    /)
    shape3d_ctra  = (/nproma, nblks_c, ntracer    /)
    shape3d_c3    = (/nproma, nblks_c, nqtendphy  /)
    shape3d_ubcp  = (/nproma, nblks_c, iadv_rcf+2 /)
    shape3d_ubcp1 = (/nproma, nblks_c, iadv_rcf+1 /)
    shape3d_ubcc  = (/nproma, nblks_c, 2  /)
    shape3d_extra = (/nproma, nlev   , nblks_c, inextra_3d  /)
    shape4d_c     = (/nproma, nlev   , nblks_c, ntracer     /)
    shape4d_chalf = (/nproma, nlevp1 , nblks_c, ntracer     /)
    shape4d_e     = (/nproma, nlev   , nblks_e, ntracer     /)
    shape4d_entl  = (/nproma, nlev   , nblks_e, n_timlevs   /)
    shape4d_chalfntl = (/nproma, nlevp1, nblks_c, n_timlevs /)


    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_diag_list, TRIM(listname), patch_id=p_patch%id )
    CALL default_var_list_settings( p_diag_list,               &
                                  & lrestart=.TRUE.,           &
                                  & restart_type=FILETYPE_NC2  )

    ! u           p_diag%u(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'Zonal wind', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 2, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'u', p_diag%u,                                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE.,                            &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=VINTP_TYPE_P_OR_Z,                           &
                &   vert_intp_method=VINTP_METHOD_UV,                           &
                &   l_hires_intp=.FALSE., l_restore_fricred=.FALSE.),           &
                & in_group=groups("atmo_ml_vars", "atmo_pl_vars", "atmo_zl_vars") )

    ! v           p_diag%v(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('northward_wind', 'm s-1', 'Meridional wind', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 3, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'v', p_diag%v,                                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE.,                            &
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=VINTP_TYPE_P_OR_Z,                           &
                &   vert_intp_method=VINTP_METHOD_UV,                           &
                &   l_hires_intp=.FALSE., l_restore_fricred=.FALSE.),           &
                & in_group=groups("atmo_ml_vars", "atmo_pl_vars", "atmo_zl_vars") )

    ! vt           p_diag%vt(nproma,nlev,nblks_e)
    ! *** needs to be saved for restart ***
    cf_desc    = t_cf_var('tangential_wind', 'm s-1', 'tangential-component of wind', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 35, ibits, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_diag_list, 'vt', p_diag%vt,                                 &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_e,                                              &
                & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_NONE ))


    ! omega_z      p_diag%omega_z(nproma,nlev,nblks_v)
    !
    cf_desc    = t_cf_var('atmosphere_relative_vorticity', 'm s-1', 'vertical vorticity', &
      &          DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 10, ibits, GRID_REFERENCE, GRID_VERTEX)
    CALL add_var( p_diag_list, 'omega_z', p_diag%omega_z,                       &
                & GRID_UNSTRUCTURED_VERT, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_v, lrestart=.FALSE.,       &
                & vert_interp=create_vert_interp_metadata( &
                &   vert_intp_type=VINTP_TYPE_P_OR_Z,      &
                &   vert_intp_method=VINTP_METHOD_LIN ),   &
                &   in_group=groups("atmo_derived_vars") )

    ! ddt_vn_phy   p_diag%ddt_vn_phy(nproma,nlev,nblks_e)
    ! *** needs to be saved for restart ***
    cf_desc    = t_cf_var('normal_wind_physical_tendency', 'm s-2',             &
      &                   'normal wind physical tendency', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 199, ibits, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_diag_list, 'ddt_vn_phy', p_diag%ddt_vn_phy,                 &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_e )


    ! ddt_exner    p_diag%ddt_exner(nproma,nlev,nblks_c)
    ! *** needs to be saved for restart ***
    cf_desc    = t_cf_var('exner_pressure_tendency', 's-1', 'exner pressure tendency', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 196, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'ddt_exner', p_diag%ddt_exner,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! ddt_exner_phy  p_diag%ddt_exner_phy(nproma,nlev,nblks_c)
    ! *** needs to be saved for restart ***
    cf_desc    = t_cf_var('exner_pressure_physical_tendency', 's-1',            &
      &                   'exner pressure physical tendency', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 197, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'ddt_exner_phy', p_diag%ddt_exner_phy,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! ddt_temp_dyn  p_diag%ddt_temp_dyn(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('dynamical_temperature_tendency', 'K s-1',            &
      &                   'dynamical temperature tendency', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 197, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'ddt_temp_dyn', p_diag%ddt_temp_dyn,             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE. )

    ! exner_old    p_diag%exner_old(nproma,nlev,nblks_c)
    ! *** needs to be saved for restart ***
    cf_desc    = t_cf_var('old_exner_pressure', '-', 'old exner pressure', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 26, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'exner_old', p_diag%exner_old,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )

    ! exner_dyn_incr    p_diag%exner_dyn_incr(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('exner_dynamics_increment', '-', 'exner dynamics increment', &
         &                DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 196, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'exner_dyn_incr', p_diag%exner_dyn_incr,         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE. )

    ! pres_sfc     p_diag%pres_sfc(nproma,nblks_c)
    !
    cf_desc    = t_cf_var('surface_air_pressure', 'Pa', 'surface pressure', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'pres_sfc', p_diag%pres_sfc,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d_c, lrestart=.FALSE. )

    ! pres_sfc_old     p_diag%pres_sfc_old(nproma,nblks_c)
    !
    cf_desc    = t_cf_var('surface_pressure', 'Pa', 'surface pressure', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'pres_sfc_old', p_diag%pres_sfc_old,             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d_c, lrestart=.FALSE. )

    ! pres_sfc_s6avg     p_diag%pres_sfc_s6avg(nproma,nblks_c)
    !
    cf_desc    = t_cf_var('surface_pressure', 'Pa', 'surface pressure', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'pres_sfc_s6avg', p_diag%pres_sfc_s6avg,         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d_c, lrestart=.FALSE. )

    ! pres_msl           p_diag%pres_msl(nproma,nblks_c)
    !
    ! Note: This task is registered for the post-processing scheduler
    !        which takes care of the regular update.
    !
    IF (l_pres_msl) THEN
      cf_desc    = t_cf_var('mean sea level pressure', 'Pa', &
        &                   'mean sea level pressure', DATATYPE_FLT32)
      grib2_desc = t_grib2_var(0, 3, 1, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'pres_msl', p_diag%pres_msl,                     &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,   &
        &           ldims=shape2d_c, lrestart=.FALSE.,                            &
        &           l_pp_scheduler_task=TASK_INTP_MSL )
    END IF

    ! temp         p_diag%temp(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('air_temperature', 'K', 'Temperature', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 0, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'temp', p_diag%temp,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE.,                            &
                & in_group=groups("atmo_ml_vars", "atmo_pl_vars", "atmo_zl_vars") )

    ! tempv        p_diag%tempv(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('virtual_temperature', 'K', 'Virtual temperature', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 0, 1, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'tempv', p_diag%tempv,                           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE. )


    ! temp_ifc     p_diag%temp_ifc(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('air_temperature', 'K', 'temperature at half level', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 0, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'temp_ifc', p_diag%temp_ifc,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf, lrestart=.FALSE. )


    ! pres         p_diag%pres(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('air_pressure', 'Pa', 'Pressure', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'pres', p_diag%pres,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE. ,                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=VINTP_TYPE_Z,                      &
                &             vert_intp_method=VINTP_METHOD_PRES ),             &
                &             in_group=groups("atmo_ml_vars", "atmo_zl_vars") )

    ! pres_ifc     p_diag%pres_ifc(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('air_pressure', 'Pa', 'pressure at half level', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 3, 0, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'pres_ifc', p_diag%pres_ifc,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf, lrestart=.FALSE. )


    ! dpres_mc     p_diag%dpres_mc(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('pressure_thickness', 'Pa', 'pressure thickness', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'dpres_mc', p_diag%dpres_mc,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE. )


    ! div          p_diag%div(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('divergence_of_wind', 's-1', 'Divergence', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 2, 13, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'div', p_diag%div,                               &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE., in_group=groups("atmo_derived_vars") )


    ! mass_fl_e    p_diag%mass_fl_e(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('horizontal_mass_flux_at_edges', 'kg m-1 s-1',        &
       &         'horizontal mass flux at edges', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_diag_list, 'mass_fl_e', p_diag%mass_fl_e,                   &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_e, lrestart=.FALSE. )


    ! mass_fl_e_sv    p_diag%mass_fl_e_sv(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('storage_field_for_horizontal_mass_flux_at_edges', 'kg m-1 s-1',  &
       &         'storage field for horizontal mass flux at edges', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_diag_list, 'mass_fl_e_sv', p_diag%mass_fl_e_sv,             &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_e, lrestart=.FALSE. )


    ! rho_ic       p_diag%rho_ic(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('density', 'kg m-3', 'density at half level', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 3, 10, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'rho_ic', p_diag%rho_ic,                         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf, lrestart=.FALSE. )


    ! w_concorr_c  p_diag%w_concorr_c(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('contravariant_vertical_correction', 'm s-1',         &
      &                   'contravariant vertical correction', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'w_concorr_c', p_diag%w_concorr_c,               &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf, lrestart=.FALSE.  )


    ! e_kinh       p_diag%e_kinh(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('horizontal specific kinetic energy', 'm2 s-2',       &
      &                   'horizontal specific kinetic energy', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 2, 196, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'e_kinh', p_diag%e_kinh,                         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_c, lrestart=.FALSE. )


    ! theta_v_ic   p_diag%theta_v_ic(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('virtual_potential_temperature_at_half_levels', 'K',&
      &                   'virtual_potential temperature at half levels', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 0, 15, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'theta_v_ic', p_diag%theta_v_ic,               &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                & ldims=shape3d_chalf, lrestart=.FALSE. )


    IF (p_patch%cell_type == 3) THEN

      ! vn_ie        p_diag%vn_ie(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('normal_wind_at_half_level', 'm s-1',               &
        &                   'normal wind at half level', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 0, 2, 197, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'vn_ie', p_diag%vn_ie,                         &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_ehalf, lrestart=.FALSE. )


      ! ddt_vn_adv   p_diag%ddt_vn_adv(nproma,nlev,nblks_e,n_timlevs)
      ! *** needs to be saved for restart (TL nnow)***
      cf_desc    = t_cf_var('advective_normal_wind_tendency', 'm s-2',          &
        &                   'advective normal wind tendency', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 0, 2, 201, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_adv', p_diag%ddt_vn_adv,               &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape4d_entl ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ALLOCATE(p_diag%ddt_vn_adv_ptr(n_timlevs))
      DO jt =1,n_timlevs
        WRITE(suffix,'(".TL",i1)') jt
        CALL add_ref( p_diag_list, 'ddt_vn_adv',                                   &
                    & 'ddt_vn_adv'//suffix, p_diag%ddt_vn_adv_ptr(jt)%p_3d,        &
                    & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT,                        &
                    & t_cf_var('ddt_adv_vn'//suffix, 'm s-2','', DATATYPE_FLT32),  &
                    & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE),&
                    & ldims=shape3d_e )
      ENDDO


      ! ddt_w_adv    p_diag%ddt_w_adv(nproma,nlevp1,nblks_c,n_timlevs)
      ! *** needs to be saved for restart (TL nnow) ***
      cf_desc    = t_cf_var('advective_vertical_wind_tendency', 'm s-2',        &
        &                   'advective vertical wind tendency', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 0, 2, 202, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_w_adv', p_diag%ddt_w_adv,                 &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape4d_chalfntl ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

      ALLOCATE(p_diag%ddt_w_adv_ptr(n_timlevs))
      DO jt =1,n_timlevs
        WRITE(suffix,'(".TL",i1)') jt
        CALL add_ref( p_diag_list, 'ddt_w_adv',                                     &
                    & 'ddt_w_adv'//suffix, p_diag%ddt_w_adv_ptr(jt)%p_3d,           &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                         &
                    & t_cf_var('ddt_adv_w'//suffix, 'm s-2','', DATATYPE_FLT32),    &
                    & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_chalf )
      ENDDO

      ! grf_tend_vn  p_diag%grf_tend_vn(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('normal_wind_tendency', 'm s-2',                    &
        &                   'normal wind tendency (grid refinement)', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 0, 2, 203, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'grf_tend_vn', p_diag%grf_tend_vn,             &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e, lrestart=.FALSE. )


      ! grf_tend_mflx  p_diag%grf_tend_mflx(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('normal_mass_flux_tendency', 'kg m-2 s-2',                    &
        &                   'normal mass flux tendency (grid refinement)', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 0, 2, 203, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'grf_tend_mflx', p_diag%grf_tend_mflx,         &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e, lrestart=.FALSE. )


      ! grf_tend_w  p_diag%grf_tend_w(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('vertical_wind_tendency', 'm s-2',                  &
        &                   'vertical wind tendency (grid refinement)', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 0, 2, 204, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'grf_tend_w', p_diag%grf_tend_w,               &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf, lrestart=.FALSE. )

      ! grf_tend_rho   p_diag%grf_tend_rho(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('density_tendency', 'kg m-3 s-1',                   &
        &                   'density tendency (grid refinement)', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 0, 3, 198, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'grf_tend_rho', p_diag%grf_tend_rho,           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_c, lrestart=.FALSE. )

      ! grf_tend_thv   p_diag%grf_tend_thv(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('virtual_potential_temperature_tendency', 'K s-1',  &
        &                   'virtual potential temperature tendency (grid refinement)', &
        &                   DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'grf_tend_thv', p_diag%grf_tend_thv,           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c, lrestart=.FALSE. )


      ! Storage fields for vertical nesting; the middle index (2) addresses 
      ! the field and its temporal tendency

      ! dvn_ie_int   p_diag%dvn_ie_int(nproma,nblks_e)
      !
      cf_desc    = t_cf_var('normal_velocity_parent_interface_level', 'm s-1',  &
        &                   'normal velocity at parent interface level', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'dvn_ie_int', p_diag%dvn_ie_int,               &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape2d_e, lrestart=.FALSE. )


      ! dvn_ie_ubc   p_diag%dvn_ie_ubc(nproma,nblks_e)
      !
      cf_desc    = t_cf_var('normal_velocity_child_upper_boundary', 'm s-1',    &
        &                   'normal velocity at child upper boundary', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'dvn_ie_ubc', p_diag%dvn_ie_ubc,               &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape2d_e, lrestart=.FALSE. )


      ! mflx_ic_int  p_diag%mflx_ic_int(nproma,nblks_c,iadv_rcf+2)
      !
      cf_desc    = t_cf_var('mass_flux_at_parent_interface_level', 'kg m-3',          &
        &                   'mass flux at parent interface level', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'mflx_ic_int', p_diag%mflx_ic_int,             &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape3d_ubcp, lrestart=.FALSE. )


      ! mflx_ic_ubc  p_diag%mflx_ic_ubc(nproma,nblks_c,2)
      !
      cf_desc    = t_cf_var('mass_flux_at_child_upper_boundary', 'kg m-3',        &
        &                   'mass flux at child upper boundary', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'mflx_ic_ubc', p_diag%mflx_ic_ubc,             &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape3d_ubcc, lrestart=.FALSE. )


      ! dtheta_v_ic_int    p_diag%dtheta_v_ic_int(nproma,nblks_c,iadv_rcf+1)
      !
      cf_desc    = t_cf_var('theta_at_parent_interface_level', 'K',             &
        &                   'potential temperature at parent interface level', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'dtheta_v_ic_int', p_diag%dtheta_v_ic_int,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape3d_ubcp1, lrestart=.FALSE. )


      ! dtheta_v_ic_ubc    p_diag%dtheta_v_ic_ubc(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('theta_at_child_upper_boundary', 'K',               &
        &                   'potential temperature at child upper boundary', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'dtheta_v_ic_ubc', p_diag%dtheta_v_ic_ubc,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape2d_c, lrestart=.FALSE. )


      ! dw_int       p_diag%dw_int(nproma,nblks_c,iadv_rcf+1)
      !
      cf_desc    = t_cf_var('w_at_parent_interface_level', 'm s-1',             &
        &                   'vertical velocity at parent interface level', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'dw_int', p_diag%dw_int,                       &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape3d_ubcp1, lrestart=.FALSE. )


      ! dw_ubc       p_diag%dw_ubc(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('w at child upper boundary', 'm s-1',               &
        &                   'vertical velocity at child upper boundary', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'dw_ubc', p_diag%dw_ubc,                       &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape2d_c, lrestart=.FALSE. )


      ! q_int        p_diag%q_int(nproma,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('q_at_parent_interface_level', 'kg kg-1',           &
        &                   'q at parent interface level', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'q_int', p_diag%q_int,                         &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape3d_ctra ,                                        &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ALLOCATE(p_diag%q_int_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'q_int',                                         &
                    & 'q_int'//ctrc, p_diag%q_int_ptr(jt)%p_2d,                     &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                        &
                    & t_cf_var('q_int'//ctrc, 'kg kg-1','', DATATYPE_FLT32),        &
                    & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape2d_c, lrestart=.FALSE. )
      ENDDO


      ! q_ubc        p_diag%q_ubc(nproma,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('q_at_child_upper_boundary', 'kg kg-1',             &
        &                   'q at child upper boundary', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'q_ubc', p_diag%q_ubc,                         &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape3d_ctra,                                         &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ALLOCATE(p_diag%q_ubc_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'q_ubc',                                         &
                    & 'q_ubc'//ctrc, p_diag%q_ubc_ptr(jt)%p_2d,                     &
                    & GRID_UNSTRUCTURED_CELL,ZAXIS_SURFACE,                         &
                    & t_cf_var('q_ubc'//ctrc, 'kg kg-1','', DATATYPE_FLT32),        &
                    & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape2d_c, lrestart=.FALSE. )
      ENDDO

    ELSE IF (p_patch%cell_type == 6) THEN
 
      ! e_kin        p_diag%e_kin(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('specific_kinetic_energy', 'm2 s-2', 'specific kinetic energy', &
           &                DATATYPE_FLT32)
      grib2_desc = t_grib2_var(0, 2, 196, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'e_kin', p_diag%e_kin,                           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                  & ldims=shape3d_c )

      ! theta_v_impl   p_diag%theta_v_impl(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('(nnow+nnew)/2_from_impl._vert._adv._of_theta_v', 'K', &
        &                   '(nnow+nnew)/2 from impl. vert. adv. of theta_v', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'theta_v_impl', p_diag%theta_v_impl,           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c )

      ! theta_v_ave   p_diag%theta_v_ave(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('time average from horiz. adv. of theta_v', 'K', &
        &                   'time average from horiz. adv. of theta_v', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'theta_v_ave', p_diag%theta_v_ave,           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,&
                  & ldims=shape3d_c )


      ! horpgrad     p_diag%horpgrad(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('covariant_horizontal_pressure_gradient', 'Pa m-1', &
        &                   'covariant horizontal pressure gradient', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'horpgrad', p_diag%horpgrad,                   &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e )


      ! vn_cov       p_diag%vn_cov(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('covariant_normal_wind', 'm s-1',                   &
        &                   'covariant normal wind', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'vn_cov', p_diag%vn_cov,                       &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e )


      ! w_cov        p_diag%w_cov(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('covariant_vertical_wind', 'm s-1',                 &
        &                   'covariant vertical wind at half level', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'w_cov', p_diag%w_cov,                         &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  & 
                  & ldims=shape3d_chalf )

      ! rho_e        p_diag%rho_e(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('density_at_edges_nnow', 'm s-1',                   &
        &                   'density_at_edges_nnow', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'rho_e', p_diag%rho_e,                         &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e )

      ! rho_star_e   p_diag%rho_star_e(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('density_at_edges_estim_step', 'm s-1',              &
        &                   'density_at_edges_estim_step', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'rho_star_e', p_diag%rho_star_e,                &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,   &
                  & ldims=shape3d_e )

      ! omega_t_con  p_diag%omega_t_con(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('tangential_horiz._contravariant_vorticity', 's-1', &
        &                   'tangential horiz. contravariant vorticity', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'omega_t_con', p_diag%omega_t_con,             &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_ehalf )


      ! omega_x      p_diag%omega_x(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('zonal_vorticity', 's-1', 'zonal vorticity', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'omega_x', p_diag%omega_x,                     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c )


      ! omega_y      p_diag%omega_y(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('meridional_vorticity', 's-1', 'meridional vorticity', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'omega_y', p_diag%omega_y,                     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c )


      ! omega_z_con  p_diag%omega_z_con(nproma,nlev,nblks_v)
      !
      cf_desc    = t_cf_var('vertical_vorticity', 's-1',                        &
        &                   'contravariant vertical vorticity', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_VERTEX)
      CALL add_var( p_diag_list, 'omega_z_con', p_diag%omega_z_con,             &
                  & GRID_UNSTRUCTURED_VERT, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_v )


      ! ddt_vn       p_diag%ddt_vn(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('normal_wind_tendency', 'm s-2', 'normal wind tendency', &
           &                DATATYPE_FLT32)
      grib2_desc = t_grib2_var(0, 2, 198, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn', p_diag%ddt_vn,                         &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                  & ldims=shape3d_e )


      ! ddt_vn_vort  p_diag%ddt_vn_vort(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('normal_wind_tendency', 'm s-2',                    &
        &                   'normal wind tendency from vorticity flux term', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_vort', p_diag%ddt_vn_vort,             &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e )


      ! ddt_w        p_diag%ddt_w(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('vertical_wind_tendency', 'm s-2', 'vertical wind tendency', &
           &                DATATYPE_FLT32)
      grib2_desc = t_grib2_var(0, 2, 200, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_w', p_diag%ddt_w,                           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                  & ldims=shape3d_chalf )

 
      ! ddt_w_vort   p_diag%ddt_w_vort(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('vert._wind_tendency', 'm s-2',                     &
        &                   'vert. wind tendency from vorticity flux term', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_w_vort', p_diag%ddt_w_vort,               &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf )

      ! ddt_w_phy   p_diag%ddt_w_phy(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('vert_wind_tendency_phy', 'm s-2',                   &
        &                   'vertical wind tendency from physics', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_w_phy', p_diag%ddt_w_phy,                  &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,   &
                  & ldims=shape3d_chalf )

    ENDIF


    !
    ! tracers
    !
    IF ( ntracer >0 ) THEN

      ! grf_tend_tracer   p_diag%grf_tend_tracer(nproma,nlev,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('tracer_tendency', 'kg kg-1 s-1',                   &
        &                   'tracer_tendency for grid refinement', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'grf_tend_tracer', p_diag%grf_tend_tracer,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape4d_c ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ALLOCATE(p_diag%ddt_grf_trc_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'grf_tend_tracer',                              &
                    & 'ddt_grf_q'//ctrc, p_diag%ddt_grf_trc_ptr(jt)%p_3d,          &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                        &
                    & t_cf_var('ddt_grf_q'//ctrc, 'kg kg-1 s**-1','', DATATYPE_FLT32), &
                    & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape3d_c, lrestart=.FALSE. )
      ENDDO


      ! hfl_tracer   p_diag%hfl_tracer(nproma,nlev,nblks_e,ntracer)
      !
      cf_desc    = t_cf_var('horizontal tracer flux', 'kg m-1 s-1',               &
        &                   'horizontal tracer flux', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'hfl_tracer', p_diag%hfl_tracer,                 &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                  & ldims=shape4d_e ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ALLOCATE(p_diag%hfl_trc_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'hfl_tracer',                                    &
                    & 'hfl_q'//ctrc, p_diag%hfl_trc_ptr(jt)%p_3d,                   &
                    & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT,                         &
                    & t_cf_var('hfl_q'//ctrc, 'kg m-1 s-1','', DATATYPE_FLT32),     &
                    & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE), &
                    & ldims=shape3d_e, lrestart=.FALSE. )
      ENDDO


      ! vfl_tracer   p_diag%vfl_tracer(nproma,nlevp1,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('vertical_tracer_flux', 'kg m-1 s-1',               &
        &                   'vertical tracer flux', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'vfl_tracer', p_diag%vfl_tracer,               &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape4d_chalf ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ALLOCATE(p_diag%vfl_trc_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'vfl_tracer',                                   &
                    & 'vfl_q'//ctrc, p_diag%vfl_trc_ptr(jt)%p_3d,                     &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                        &
                    & t_cf_var('vfl_q'//ctrc, 'kg m-1 s-1','', DATATYPE_FLT32),     &
                    & t_grib2_var(255,255, 255, ibits, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_chalf, lrestart=.FALSE. )
      ENDDO


      ! ddt_tracer_adv   p_diag%ddt_tracer_adv(nproma,nlev,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('advective tracer tendency', 'kg kg-1 s-1',         &
        &                   'advective tracer tendency', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_tracer_adv', p_diag%ddt_tracer_adv,       &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape4d_c ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ALLOCATE(p_diag%ddt_trc_adv_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'ddt_tracer_adv',                               &
                    & 'ddt_adv_q'//ctrc, p_diag%ddt_trc_adv_ptr(jt)%p_3d,             &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                        &
                    & t_cf_var('ddt_adv_q'//ctrc, 'kg kg-1 s-1','', DATATYPE_FLT32),&
                    & t_grib2_var(255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape3d_c, lrestart=.FALSE. )
      ENDDO


      ! tracer_vi(nproma,nblks_c,3), only Q1, Q2, Q3
      cf_desc    = t_cf_var('tracer_vi', '', 'tracer_vi', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'tracer_vi', p_diag%tracer_vi,                  &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c3, lrestart=.FALSE., loutput=.FALSE.,         &
                  & lcontainer=.TRUE.)

      ALLOCATE(p_diag%tracer_vi_ptr(nqtendphy))
      DO jt =1,nqtendphy
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'tracer_vi', 'tracer_vi'//ctrc,               &
          &           p_diag%tracer_vi_ptr(jt)%p_2d,                             &
          &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
          &           cf_desc, grib2_desc, ldims=shape2d_c, lrestart=.FALSE.)
      ENDDO

      ! tracer_vi_avg(nproma,nblks_c,3), only Q1, Q2, Q3
      cf_desc    = t_cf_var('tracer_vi_avg', '', 'tracer_vi_avg', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'tracer_vi_avg', p_diag%tracer_vi_avg,          &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c3, lrestart=.FALSE., loutput=.FALSE.,         &
                  & lcontainer=.TRUE.)

      ALLOCATE(p_diag%tracer_vi_avg_ptr(nqtendphy))
      DO jt =1,nqtendphy
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'tracer_vi_avg', 'tracer_vi_avg'//ctrc,       &
          &           p_diag%tracer_vi_avg_ptr(jt)%p_2d,                         &
          &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                     &
          &           cf_desc, grib2_desc, ldims=shape2d_c, lrestart=.FALSE. )
      ENDDO
    ENDIF



    IF( lwrite_extra) THEN
      WRITE(0,*)'inextra_2d=',inextra_2d

      IF(inextra_2d > 0) THEN

        ! extra_2d   p_diag%extra_2d(nproma,nblks_c,inextra_2d)
        !
        cf_desc    = t_cf_var('extra_field_2D', '-', 'extra field 2D', DATATYPE_FLT32)
        grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_diag_list, 'extra_2d', p_diag%extra_2d,                   &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                    & lcontainer=.TRUE., ldims=shape2d_extra, lrestart=.FALSE. )

        ALLOCATE(p_diag%extra_2d_ptr(inextra_2d))
        DO jt =1,inextra_2d
          WRITE(ctrc,'(I2)')jt
          CALL add_ref( p_diag_list, 'extra_2d', 'extra_2d'//TRIM(ADJUSTL(ctrc)), &
            &           p_diag%extra_2d_ptr(jt)%p_2d,                             &
            &           GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE,                    &
            &           cf_desc, grib2_desc, ldims=shape2d_c, lrestart=.FALSE. )
        ENDDO
      ENDIF


      IF(inextra_3d > 0) THEN

        ! extra_3d   p_diag%extra_3d(nproma,nlev,nblks_c,inextra_3d)
        !
        cf_desc    = t_cf_var('extra_fields_3D', '-', 'extra fields 3D', DATATYPE_FLT32)
        grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_diag_list, 'extra_3d', p_diag%extra_3d,                   &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                    & ldims=shape3d_extra,                                        &
                    & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

        ALLOCATE(p_diag%extra_3d_ptr(inextra_3d))
        DO jt =1,inextra_3d
          WRITE(ctrc,'(I2)')jt
          CALL add_ref( p_diag_list, 'extra_3d', 'extra_3d'//TRIM(ADJUSTL(ctrc)), &
            &           p_diag%extra_3d_ptr(jt)%p_3d,                             &
            &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                     &
            &           cf_desc, grib2_desc, ldims=shape3d_c, lrestart=.FALSE. )
        ENDDO
      ENDIF
    ENDIF

  END SUBROUTINE new_nh_state_diag_list


  !-------------------------------------------------------------------------
  !>
  !! Allocation of components of reference state.
  !!
  !! Initialization of components with zero.
  !!
  !! @par Revision History
  !! Initial release by Daniel Reinert, DWD (2012-06-04)
  !!
  !!
  SUBROUTINE new_nh_state_ref_list ( p_patch, p_ref, p_ref_list,  &
    &                                 listname )
!
    TYPE(t_patch), TARGET, INTENT(IN) :: &  !< current patch
      &  p_patch

    TYPE(t_nh_ref),  INTENT(INOUT)   :: &  !< reference state
      &  p_ref 

    TYPE(t_var_list), INTENT(INOUT)   :: &  !< reference state list
      &  p_ref_list
    CHARACTER(len=*), INTENT(IN)      :: &  !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, nblks_e    !< number of cell/edge blocks to allocate

    INTEGER :: nlev, nlevp1  !< number of vertical full/half levels

    INTEGER :: shape3d_e(3), shape3d_chalf(3)
 
    INTEGER :: ibits         !< "entropy" of horizontal slice

    CHARACTER(len=4) suffix
    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ibits = DATATYPE_PACK   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape3d_e     = (/nproma, nlev   , nblks_e    /)
    shape3d_chalf = (/nproma, nlevp1 , nblks_c    /)


    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_ref_list, TRIM(listname), patch_id=p_patch%id )
    CALL default_var_list_settings( p_ref_list,                &
                                  & lrestart=.FALSE.,          &
                                  & restart_type=FILETYPE_NC2  )

    ! vn_ref     p_ref%vn_ref(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('normal_velocity', 'm s-1', 'velocity normal to edge', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 34, ibits, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_ref_list, 'vn_ref', p_ref%vn_ref,                           &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_e, lrestart=.FALSE.,                            &
                & in_group=groups("atmo_ml_vars", "atmo_pl_vars", "atmo_zl_vars"), &
                & isteptype=TSTEP_CONSTANT )


    ! w_ref      p_ref%w_ref(nproma,nlev+1,nblks_c)
    !
    cf_desc    = t_cf_var('upward_air_velocity', 'm s-1', 'Vertical velocity', DATATYPE_FLT32)
    grib2_desc = t_grib2_var(0, 2, 9, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_ref_list, 'w_ref', p_ref%w_ref,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf, lrestart=.FALSE.,                        &
                & in_group=groups("atmo_ml_vars", "atmo_pl_vars", "atmo_zl_vars"), &
                & isteptype=TSTEP_CONSTANT )

  END SUBROUTINE new_nh_state_ref_list


  !---------------------------------------------------------------------------
  !>
  !! Allocates all metric coefficients defined in type metrics_3d of the patch.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann (2009-04-14)
  !! Modification by Daniel Reinert, DWD, (2010-04-22)
  !! - added geometric height at full levels
  SUBROUTINE new_nh_metrics_list ( p_patch, p_metrics, p_metrics_list,  &
    &                              listname )
!
    TYPE(t_patch), TARGET, INTENT(IN) :: &  !< current patch
      &  p_patch

    TYPE(t_nh_metrics),  INTENT(INOUT):: &  !< diagnostic state
      &  p_metrics 

    TYPE(t_var_list), INTENT(INOUT) :: p_metrics_list   !< diagnostic state list

    CHARACTER(len=*), INTENT(IN)      :: &  !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
               nblks_e, &    !< number of edge blocks to allocate
               nblks_v       !< number of vertex blocks to allocate

    INTEGER :: nlev, nlevp1

    INTEGER :: shape2d_c(2), shape3d_c(3), shape3d_e(3),               &
      &        shape3d_v(3), shape3d_chalf(3), shape3d_ehalf(3),       &
      &        shape2d_ccubed(3), shape2d_ecubed(3), shape3d_vhalf(3), & 
      &        shape3d_esquared(4), shape3d_e8(4)
    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: ist
    !--------------------------------------------------------------

    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nblks_v = p_patch%nblks_v

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d_c        = (/nproma,          nblks_c    /)     
    shape2d_ccubed   = (/nproma, 3      , nblks_c    /)     
    shape2d_ecubed   = (/nproma, 3      , nblks_e    /)     
    shape3d_c        = (/nproma, nlev   , nblks_c    /)     
    shape3d_chalf    = (/nproma, nlevp1 , nblks_c    /)      
    shape3d_e        = (/nproma, nlev   , nblks_e    /)     
    shape3d_ehalf    = (/nproma, nlevp1 , nblks_e    /)     
    shape3d_esquared = (/2     , nproma , nlev   , nblks_e /)
    shape3d_e8       = (/8     , nproma , nlev   , nblks_e /)
    shape3d_v        = (/nproma, nlev   , nblks_v    /)     
    shape3d_vhalf    = (/nproma, nlevp1 , nblks_v    /)


    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_metrics_list, TRIM(listname), patch_id=p_patch%id )
    CALL default_var_list_settings( p_metrics_list,            &
                                  & lrestart=.FALSE.,          &
                                  & restart_type=FILETYPE_NC2  )

    ! geometric height at the vertical interface of cells
    ! z_ifc        p_metrics%z_ifc(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('geometric_height_at_half_level_center', 'm',         &
      &                   'geometric height at half level center', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 3, 6, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'z_ifc', p_metrics%z_ifc,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID_HALF, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf,                                          & 
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=VINTP_TYPE_P_OR_Z,                           &
                &   vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ),                 &
                & in_group=groups("atmo_ml_vars", "atmo_pl_vars"),              &
                & isteptype=TSTEP_CONSTANT )


    ! geometric height at full levels
    ! z_mc         p_metrics%z_mc(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('geometric_height_at_full_level_center', 'm',         &
      &                   'geometric height at full level center', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 3, 6, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'z_mc', p_metrics%z_mc,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HYBRID, cf_desc, grib2_desc,    &
                & ldims=shape3d_c,                                              &
                & isteptype=TSTEP_CONSTANT )


    ! slope of the terrain in normal direction (half level)
    ! ddxn_z_half  p_metrics%ddxn_z_half(nproma,nlevp1,nblks_e)
    !
    cf_desc    = t_cf_var('terrain_slope_in_normal_direction', '-',             &
      &                   'terrain slope in normal direction', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_metrics_list, 'ddxn_z_half', p_metrics%ddxn_z_half,         &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_ehalf,                                          &
                & isteptype=TSTEP_CONSTANT )


    ! slope of the terrain in normal direction (full level)
    ! ddxn_z_full  p_metrics%ddxn_z_full(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('terrain_slope_in_normal_direction', '-',             &
      &                   'terrain slope in normal direction', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_metrics_list, 'ddxn_z_full', p_metrics%ddxn_z_full,         &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_e,                                              &
                & isteptype=TSTEP_CONSTANT )


    ! slope of the terrain in tangential direction (half level)
    ! ddxt_z_half  p_metrics%ddxt_z_half(nproma,nlevp1,nblks_e)
    !
    cf_desc    = t_cf_var('terrain_slope_in_tangential_direction', '-',         &
      &                   'terrain slope in tangential direction', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_metrics_list, 'ddxt_z_half', p_metrics%ddxt_z_half,         &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_ehalf,                                          &
                & isteptype=TSTEP_CONSTANT )

    ! slope of the terrain in tangential direction (full level)
    ! ddxt_z_full  p_metrics%ddxt_z_full(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('terrain_slope_in_tangential_direction', '-',         &
      &                   'terrain slope in tangential direction', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_metrics_list, 'ddxt_z_full', p_metrics%ddxt_z_full,         &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_e,                                              &
                & isteptype=TSTEP_CONSTANT )


    ! functional determinant of the metrics [sqrt(gamma)]
    ! ddqz_z_full  p_metrics%ddqz_z_full(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('metrics_functional_determinant', '-',                &
      &                   'metrics functional determinant', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'ddqz_z_full', p_metrics%ddqz_z_full,         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c,                                              &
                & isteptype=TSTEP_CONSTANT )


    ! functional determinant of the metrics [sqrt(gamma)]
    ! ddqz_z_full_e  p_metrics%ddqz_z_full_e(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('metrics_functional_determinant', '-',                &
      &                   'metrics functional determinant (edge)', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_metrics_list, 'ddqz_z_full_e', p_metrics%ddqz_z_full_e,     &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_e,                                              &
                & isteptype=TSTEP_CONSTANT )


    ! functional determinant of the metrics [sqrt(gamma)]
    ! ddqz_z_half  p_metrics%ddqz_z_half(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('metrics_functional_determinant', '-',                &
      &                   'metrics functional determinant', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'ddqz_z_half', p_metrics%ddqz_z_half,         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf,                                          &
                & isteptype=TSTEP_CONSTANT )


    ! geopotential at full level cell center
    ! geopot       p_metrics%geopot(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('geopotential', 'm2 s-2',                             &
       &                  'geopotential at full level cell centre', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 3, 4, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'geopot', p_metrics%geopot,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c,                                              &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=VINTP_TYPE_P_OR_Z,                 &
                &             vert_intp_method=VINTP_METHOD_LIN,                &
                &             l_extrapol=.TRUE., l_pd_limit=.FALSE.)            )


    ! geopotential above groundlevel at cell center
    ! geopot_agl   p_metrics%geopot_agl(nproma,nlev  ,nblks_c)
    !
    cf_desc    = t_cf_var('geopotential', 'm2 s-2',                             &
      &                   'geopotential above groundlevel at cell center', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 3, 4, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'geopot_agl', p_metrics%geopot_agl,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! geopotential above groundlevel at cell center
    ! geopot_agl_ifc  p_metrics%geopot_agl_ifc(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('geopotential', 'm2 s-2',                             &
      &                   'geopotential above groundlevel at cell center', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 3, 4, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'geopot_agl_ifc', p_metrics%geopot_agl_ifc,   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! geopotential at cell center
    ! dgeopot_mc   p_metrics%dgeopot_mc(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('geopotential', 'm2 s-2',                             &
      &                   'geopotential at cell center', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 0, 3, 4, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'dgeopot_mc', p_metrics%dgeopot_mc,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )

    !------------------------------------------------------------------------------
    !HW: Vertical 1D arrays are not yet supported by add_var. Use allocate for now.
    ! Since this is meant to be a temporary solution, deallocated is not implemented.

      ! Rayleigh damping coefficient for w
      ALLOCATE(p_metrics%rayleigh_w(nlevp1),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for rayleigh_w failed')
      ENDIF

      ! Rayleigh damping coefficient for vn
      ALLOCATE(p_metrics%rayleigh_vn(nlev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for rayleigh_vn failed')
      ENDIF

      ! Background nabla2 diffusion coefficient for upper sponge layer
      ALLOCATE(p_metrics%enhfac_diffu(nlev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nonhydro_state:construct_nh_metrics', &
                    'allocation for enhfac_diffu failed')
      ENDIF
    !----------------------------------------------------------------------------

    ! Explicit weight in vertical wind solver
    ! vwind_expl_wgt   p_metrics%vwind_expl_wgt(nproma,nblks_c)
    !
    cf_desc    = t_cf_var('Explicit_weight_in_vertical_wind_solver', '-',       &
      &                   'Explicit weight in vertical wind solver', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'vwind_expl_wgt', p_metrics%vwind_expl_wgt,   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,        &
                & ldims=shape2d_c )


    ! Implicit weight in vertical wind solver
    ! vwind_impl_wgt  p_metrics%vwind_impl_wgt(nproma,nblks_c)
    !
    cf_desc    = t_cf_var('Implicit_weight_in_vertical_wind_solver', '-',       &
      &                   'Implicit weight in vertical wind solver', DATATYPE_FLT32)
    grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'vwind_impl_wgt', p_metrics%vwind_impl_wgt,   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,        &
                & ldims=shape2d_c )



! These fields are needed for triangles only once the initialization in
! mo_nh_testcases is properly rewritten for hexagons
!    IF (p_patch%cell_type== 3) THEN
      ! weighting factor for interpolation from full to half levels
      ! wgtfac_c     p_metrics%wgtfac_c(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                     'weighting factor for interpolation from full to half levels', &
      &                     DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'wgtfac_c', p_metrics%wgtfac_c,             &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf,                                        &
                  & isteptype=TSTEP_CONSTANT )


      ! weighting factor for interpolation from full to half levels
      ! wgtfac_e     p_metrics%wgtfac_e(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                     'weighting factor for interpolation from full to half levels', &
      &                     DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'wgtfac_e', p_metrics%wgtfac_e,             &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_ehalf,                                        &
                  & isteptype=TSTEP_CONSTANT )


      ! weighting factor for quadratic interpolation to surface
      ! wgtfacq_c    p_metrics%wgtfacq_c(nproma,3,nblks_c)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                     'weighting factor for quadratic interpolation to surface', &
      &                     DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'wgtfacq_c', p_metrics%wgtfacq_c,           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape2d_ccubed, loutput=.FALSE.,                      &
                  & isteptype=TSTEP_CONSTANT )


      ! weighting factor for quadratic interpolation to surface
      ! wgtfacq_e    p_metrics%wgtfacq_e(nproma,3,nblks_e)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                     'weighting factor for quadratic interpolation to surface', &
      &                     DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'wgtfacq_e', p_metrics%wgtfacq_e,           &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape2d_ecubed, loutput=.FALSE.,                      &
                  & isteptype=TSTEP_CONSTANT )


      ! weighting factor for quadratic interpolation to model top
      ! wgtfacq1_c    p_metrics%wgtfacq1_c(nproma,3,nblks_c)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                     'weighting factor for quadratic interpolation to model top', &
      &                     DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'wgtfacq1_c', p_metrics%wgtfacq1_c,         &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape2d_ccubed, loutput=.FALSE.,                      &
                  & isteptype=TSTEP_CONSTANT )


      ! weighting factor for quadratic interpolation to model top
      ! wgtfacq1_e   p_metrics%wgtfacq1_e(nproma,3,nblks_e)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                     'weighting factor for quadratic interpolation to model top', &
      &                     DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'wgtfacq1_e', p_metrics%wgtfacq1_e,         &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape2d_ecubed, loutput=.FALSE.,                      &
                  & isteptype=TSTEP_CONSTANT )


      ! Inverse layer thickness of full levels
      ! inv_ddqz_z_full   p_metrics%inv_ddqz_z_full(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Inverse_layer_thickness', 'm-1',                   &
      &                     'Inverse layer thickness of full levels', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'inv_ddqz_z_full', p_metrics%inv_ddqz_z_full, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                  & ldims=shape3d_c,                                              &
                  & isteptype=TSTEP_CONSTANT )


      ! Inverse distance between full levels jk+1 and jk-1
      ! inv_ddqz_z_half2  p_metrics%inv_ddqz_z_half2(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Inverse_distance_between_full_levels', 'm-1',      &
      &                     'Inverse distance between full levels jk+1 and jk-1', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'inv_ddqz_z_half2', p_metrics%inv_ddqz_z_half2, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,      &
                  & ldims=shape3d_c,                                                &
                  & isteptype=TSTEP_CONSTANT )

      ! Reference atmosphere field exner
      ! exner_ref_mc  p_metrics%exner_ref_mc(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_exner', '-',            &
      &                     'Reference atmosphere field exner', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'exner_ref_mc', p_metrics%exner_ref_mc,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c,                                            &
                  & isteptype=TSTEP_CONSTANT )


    IF (p_patch%cell_type== 3) THEN

      ! Vertical index of neighbor points needed for Taylor-expansion-based pressure gradient
      ! vertidx_gradp  p_metrics%vertidx_gradp(2,nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('Vertical_index', '-',                              &
      &                     'Vertical index', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'vertidx_gradp', p_metrics%vertidx_gradp,   &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_esquared, loutput=.FALSE.,                    &
                  & isteptype=TSTEP_CONSTANT )


      IF (igradp_method <= 3) THEN
        ! Height differences between local edge point and neighbor cell points used for
        ! pressure gradient computation
        ! zdiff_gradp  p_metrics%zdiff_gradp(2,nproma,nlev,nblks_e)
        !
        cf_desc    = t_cf_var('Height_differences', 'm',                          &
        &                     'Height differences', DATATYPE_FLT32)
        grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
        CALL add_var( p_metrics_list, 'zdiff_gradp', p_metrics%zdiff_gradp,       &
                    & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                    & ldims=shape3d_esquared, loutput=.FALSE.,                    &
                    & isteptype=TSTEP_CONSTANT  )
      ELSE
        ! Coefficients for cubic interpolation of Exner pressure
        ! coeff_gradp  p_metrics%coeff_gradp(8,nproma,nlev,nblks_e)
        !
        cf_desc    = t_cf_var('Interpolation_coefficients', '-',                  &
        &                     'Interpolation coefficients', DATATYPE_FLT32)
        grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
        CALL add_var( p_metrics_list, 'coeff_gradp', p_metrics%coeff_gradp,       &
                    & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                    & ldims=shape3d_e8, loutput=.FALSE.,                          &
                    & isteptype=TSTEP_CONSTANT  )
      ENDIF


      ! Extrapolation factor for Exner pressure
      ! exner_exfac  p_metrics%exner_exfac(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Extrapolation_factor_for_Exner_pressure', '-',     &
      &                     'Extrapolation factor for Exner pressure', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'exner_exfac', p_metrics%exner_exfac,       &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c,                                            &
                  & isteptype=TSTEP_CONSTANT )


      ! Reference atmosphere field theta
      ! theta_ref_mc  p_metrics%theta_ref_mc(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_theta', 'K',            &
      &                     'Reference atmosphere field theta', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'theta_ref_mc', p_metrics%theta_ref_mc,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c,                                            &
                  & isteptype=TSTEP_CONSTANT )


      ! Reference atmosphere field theta (edges)
      ! theta_ref_me  p_metrics%theta_ref_me(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_theta', 'K',            &
      &                     'Reference atmosphere field theta', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'theta_ref_me', p_metrics%theta_ref_me,     &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e,                                            &
                  & isteptype=TSTEP_CONSTANT )


      ! Reference atmosphere field theta
      ! theta_ref_ic  p_metrics%theta_ref_ic(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_theta', 'K',            &
      &                     'Reference atmosphere field theta', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'theta_ref_ic', p_metrics%theta_ref_ic,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf,                                        &
                  & isteptype=TSTEP_CONSTANT )


      ! Reference surface temperature
      ! tsfc_ref  p_metrics%tsfc_ref(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_surface_temperature', 'K',               &
      &                     'Reference surface temperature', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'tsfc_ref', p_metrics%tsfc_ref,             &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape2d_c,                                            &
                  & isteptype=TSTEP_CONSTANT )


      ! Reference atmosphere field density
      ! rho_ref_mc  p_metrics%rho_ref_mc(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_density', '-',          &
      &                     'Reference atmosphere field density', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'rho_ref_mc', p_metrics%rho_ref_mc,         &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c,                                            &
                  & isteptype=TSTEP_CONSTANT )


      ! Reference atmosphere field density (edges)
      ! rho_ref_me  p_metrics%rho_ref_me(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_density', '-',          &
      &                     'Reference atmosphere field density', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'rho_ref_me', p_metrics%rho_ref_me,         &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e,                                            &
                  & isteptype=TSTEP_CONSTANT )


      ! Reference atmosphere field exner
      ! d_exner_dz_ref_ic  p_metrics%d_exner_dz_ref_ic(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_exner', '-',            &
      &                     'Reference atmosphere field exner', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'd_exner_dz_ref_ic', p_metrics%d_exner_dz_ref_ic, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf,                                        &
                  & isteptype=TSTEP_CONSTANT )


      IF (igradp_method <= 3) THEN
        ! Reference atmosphere field exner
        ! d2dexdz2_fac1_mc  p_metrics%d2dexdz2_fac1_mc(nproma,nlev,nblks_c)
        !
        cf_desc    = t_cf_var('Reference_atmosphere_field_exner', '-',            &
        &                     'Reference atmosphere field exner', DATATYPE_FLT32)
        grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_metrics_list, 'd2dexdz2_fac1_mc', p_metrics%d2dexdz2_fac1_mc, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,      &
                    & ldims=shape3d_c,                                            &
                    & isteptype=TSTEP_CONSTANT )


        ! Reference atmosphere field exner
        ! d2dexdz2_fac2_mc  p_metrics%d2dexdz2_fac2_mc(nproma,nlev,nblks_c)
        !
        cf_desc    = t_cf_var('Reference_atmosphere_field_exner', '-',            &
        &                     'Reference atmosphere field exner', DATATYPE_FLT32)
        grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_metrics_list, 'd2dexdz2_fac2_mc', p_metrics%d2dexdz2_fac2_mc, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,      &
                    & ldims=shape3d_c,                                            &
                    & isteptype=TSTEP_CONSTANT )
      ENDIF

      ! mask field that excludes boundary halo points
      ! mask_prog_halo_c  p_metrics%mask_prog_halo_c(nproma,nblks_c)
      ! Note: Here "loutput" is set to .FALSE. since the output
      !       scheme operates on REAL model variables only and
      !       throws an error on this.
      !
      cf_desc    = t_cf_var('mask_field', '-',                                  &
      &                     'mask field that excludes boundary halo points', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'mask_prog_halo_c', p_metrics%mask_prog_halo_c, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,     &
                  & ldims=shape2d_c, loutput=.FALSE.,                           &
                  & isteptype=TSTEP_CONSTANT )


    ELSE IF (p_patch%cell_type== 6) THEN


      ! geometric height at full level edges
      ! z_mc_e       p_metrics%z_mc_e(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('geometric_height_at_full_level_edge', 'm',           &
        &                   'geometric height at full level edge', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'z_mc_e', p_metrics%z_mc_e,                   &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                  & ldims=shape3d_e,                                              &
                  & isteptype=TSTEP_CONSTANT )


      ! slope of the coordinate lines in Northern direction
      ! ddnorth_z    p_metrics%ddnorth_z(nproma,nlev,nblks_v)
      !
      cf_desc    = t_cf_var('slope_of_coordinate_lines_in_Northern_direction', '-', &
      &                     'slope of coordinate lines in Northern direction', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_VERTEX)
      CALL add_var( p_metrics_list, 'ddnorth_z', p_metrics%ddnorth_z,               &
                  & GRID_UNSTRUCTURED_VERT, ZAXIS_HEIGHT, cf_desc, grib2_desc,      &
                  & ldims=shape3d_v,                                                &
                  & isteptype=TSTEP_CONSTANT )


      ! slope of the coordinate lines in Eastern direction
      ! ddeast_z     p_metrics%ddeast_z(nproma,nlev,nblks_v)
      !
      cf_desc    = t_cf_var('slope_of_coordinate_lines_in_Eastern_direction', '-', &
      &                     'slope of coordinate lines in Eastern direction', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_VERTEX)
      CALL add_var( p_metrics_list, 'ddeast_z', p_metrics%ddeast_z,                &
                  & GRID_UNSTRUCTURED_VERT, ZAXIS_HEIGHT, cf_desc, grib2_desc,     &
                  & ldims=shape3d_v,                                               &
                  & isteptype=TSTEP_CONSTANT )


      ! functional determinant of the metrics [sqrt(gamma)]
      ! ddqz_z_full_v   p_metrics%ddqz_z_full_v(nproma,nlev,nblks_v)
      !
      cf_desc    = t_cf_var('metrics_functional_determinant', '-',              &
      &                     'metrics functional determinant', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_VERTEX)
      CALL add_var( p_metrics_list, 'ddqz_z_full_v', p_metrics%ddqz_z_full_v,   &
                  & GRID_UNSTRUCTURED_VERT, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_v,                                            &
                  & isteptype=TSTEP_CONSTANT )


      ! functional determinant of the metrics [sqrt(gamma)]
      ! ddqz_z_full_r   p_metrics%ddqz_z_full_r(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('metrics_functional_determinant', '-',              &
      &                     'metrics functional determinant', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'ddqz_z_full_r', p_metrics%ddqz_z_full_r,   &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e,                                            &
                  & isteptype=TSTEP_CONSTANT )


      ! functional determinant of the metrics [sqrt(gamma)]
      ! ddqz_z_half_e   p_metrics%ddqz_z_half_e(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('metrics_functional_determinant', '-',              &
      &                     'metrics functional determinant', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'ddqz_z_half_e', p_metrics%ddqz_z_half_e,   &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_ehalf,                                        &
                  & isteptype=TSTEP_CONSTANT )


      ! 1/dz(k-1)-1/dz(k)
      ! diff_1_o_dz  p_metrics%diff_1_o_dz(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('difference_1_over_dz', 'm-1', 'difference 1 over dz', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'diff_1_o_dz', p_metrics%diff_1_o_dz,       &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf,                                        &
                  & isteptype=TSTEP_CONSTANT )
               

      ! 1/(dz(k-1)*dz(k))
      ! mult_1_o_dz  p_metrics%mult_1_o_dz(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('mult_1_over_dz', 'm-2', 'mult 1 over dz', DATATYPE_FLT32)
      grib2_desc = t_grib2_var( 255, 255, 255, ibits, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'mult_1_o_dz', p_metrics%mult_1_o_dz,       &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf,                                        &
                  & isteptype=TSTEP_CONSTANT )

    ENDIF

  END SUBROUTINE new_nh_metrics_list

END MODULE mo_nonhydro_state





