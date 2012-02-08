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
  USE mo_impl_constants,       ONLY: SUCCESS, MAX_CHAR_LENGTH, INWP
  USE mo_exception,            ONLY: message, finish
  USE mo_model_domain,         ONLY: t_patch
  USE mo_nonhydro_types,       ONLY: t_nh_state, t_nh_prog, t_nh_diag,        &
    &                                t_nh_metrics, t_ptr_nh, t_nh_diag_pz,    &
    &                                t_buffer_memory
  USE mo_model_domain_import,  ONLY: n_dom, l_limited_area
  USE mo_nonhydrostatic_config,ONLY: itime_scheme, l_nest_rcf
  USE mo_dynamics_config,      ONLY: nsav1, nsav2
  USE mo_parallel_config,      ONLY: nproma
  USE mo_run_config,           ONLY: iforcing, &!ltransport,    &
    &                                ntracer, ntracer_static,      &
    &                                iqv, iqc, iqi, iqr, iqs, io3, &
    &                                nqtendphy
  USE mo_radiation_config,     ONLY: irad_o3
  USE mo_io_config,            ONLY: lwrite_extra, inextra_2d, inextra_3d, &
    &                                lwrite_pzlev
  USE mo_nh_pzlev_config,      ONLY: nh_pzlev_config
  USE mo_linked_list,          ONLY: t_var_list
  USE mo_var_list,             ONLY: default_var_list_settings, add_var,     &
    &                                add_ref, new_var_list, delete_var_list, &
    &                                create_tracer_metadata, add_var_list_reference
  USE mo_linked_list,          ONLY: t_list_element, find_list_element
  USE mo_var_metadata,         ONLY: t_var_metadata, t_tracer_meta
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants


  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


  PUBLIC :: construct_nh_state    ! Constructor for the nonhydrostatic state
  PUBLIC :: destruct_nh_state     ! Destructor for the nonhydrostatic state
  PUBLIC :: duplicate_prog_state  ! Copy the prognostic state

  PUBLIC :: p_nh_state            ! state vector of nonhydrostatic variables (variable)
  PUBLIC :: bufr


!!$  INTERFACE add_tracer_ref
!!$    MODULE PROCEDURE add_var_list_reference_tracer
!!$  END INTERFACE add_tracer_ref


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
  SUBROUTINE construct_nh_state(p_patch, p_nh_state, n_timelevels)
!
    TYPE(t_patch),     INTENT(IN)   ::  & ! patch
      &  p_patch(n_dom)

    TYPE(t_nh_state),  INTENT(INOUT)::  & ! nh state at different grid levels
      &  p_nh_state(n_dom)

    INTEGER, OPTIONAL, INTENT(IN)   ::  & ! number of timelevels
      &  n_timelevels    

    INTEGER  :: ntl, &    ! local number of timelevels
                ist, &    ! status
                jg,  &    ! grid level counter
                jt        ! time level counter

    LOGICAL  :: l_extra_timelev

    CHARACTER(len=MAX_CHAR_LENGTH) :: listname, varname_prefix

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nonhydro_state:construct_nh_state'
!-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'Construction of NH state started')

    DO jg = 1, n_dom

      IF(PRESENT(n_timelevels))THEN
        ntl = n_timelevels
      ELSE
        ntl = 1
      ENDIF

      ! If grid nesting is not called at every dynamics time step, an extra time
      ! level is needed for full-field interpolation and boundary-tendency calculation
      IF (l_nest_rcf .AND. n_dom > 1) THEN
        ntl = ntl + 1
        nsav1(jg) = ntl
      ENDIF

      ! In the presence of grid nesting, another extra time level is needed to save
      ! the feedback increments
      ! This extra time level is also used to store the driving-model data in the
      ! limited-area mode
      IF (l_limited_area .OR. jg > 1) ntl = ntl + 1
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

      ! create tracer list
      ALLOCATE(p_nh_state(jg)%tracer_list(1:ntl), STAT=ist)
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
      ENDDO

      !
      ! Build diag state list
      ! includes memory allocation
      !
      WRITE(listname,'(a,i2.2)') 'nh_state_diag_of_domain_',jg
      CALL new_nh_state_diag_list(p_patch(jg), p_nh_state(jg)%diag, &
        &  p_nh_state(jg)%diag_list, listname)


      ! Build diag_p and diag_z state lists
      ! includes memory allocation
      !
      IF (lwrite_pzlev) THEN
        IF (nh_pzlev_config(jg)%lwrite_plev) THEN
          WRITE(listname,'(a,i2.2)') 'nh_state_diag_p_of_domain_',jg
          CALL new_nh_state_diag_p_list(p_patch(jg), p_nh_state(jg)%diag_p, &
            &  p_nh_state(jg)%diag_p_list, listname)
        ENDIF

        WRITE(listname,'(a,i2.2)') 'nh_state_diag_z_of_domain_',jg
        CALL new_nh_state_diag_z_list(p_patch(jg), p_nh_state(jg)%diag_z, &
          &  p_nh_state(jg)%diag_z_list, listname)
      ENDIF


      !
      ! Build metrics state list
      ! includes memory allocation
      !
      WRITE(listname,'(a,i2.2)') 'nh_state_metrics_of_domain_',jg
      CALL new_nh_metrics_list(p_patch(jg), p_nh_state(jg)%metrics, &
        &  p_nh_state(jg)%metrics_list, listname )

    ENDDO

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
                                             
    INTEGER  :: ntl, &    ! local number of timelevels
                ist, &    ! status
                jg,  &    ! grid level counter
                jt        ! time level counter

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nonhydro_state:destruct_nh_state'
!-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'Destruction of NH state started')


    DO jg = 1, n_dom

      ntl = SIZE(p_nh_state(jg)%prog(:))
      IF(ntl==0)THEN
        CALL finish(TRIM(routine), 'prognostic array has no timelevels')
      ENDIF

      ! delete diagnostic state list elements
      CALL delete_var_list( p_nh_state(jg)%diag_list )

      IF (lwrite_pzlev) THEN
        CALL delete_var_list( p_nh_state(jg)%diag_z_list )

        IF (nh_pzlev_config(jg)%lwrite_plev) THEN
          CALL delete_var_list( p_nh_state(jg)%diag_p_list )
        ENDIF
      ENDIF


      ! delete metrics state list elements
      CALL delete_var_list( p_nh_state(jg)%metrics_list )


      ! delete prognostic state list elements
      DO jt = 1, ntl
        CALL delete_var_list( p_nh_state(jg)%prog_list(jt) )
      ENDDO

      ! delete tracer list list elements
      DO jt = 1, 2  ! no extra time level  !! quick fix !!
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

    LOGICAL, INTENT(IN) :: l_extra_timelev  !< specifies extra time levels for which not all variables are allocated

    INTEGER, INTENT(IN) :: timelev


    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
               nblks_e       !< number of edge blocks to allocate

    INTEGER :: nlev, nlevp1, ktracer

    INTEGER :: shape3d_c(3), shape3d_e(3), shape3d_chalf(3), &
      &        shape4d_c(4)

    INTEGER :: ientr         !< "entropy" of horizontal slice

    CHARACTER(len=4) suffix

    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ientr = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape3d_e     = (/nproma, nlev,   nblks_e  /)
    shape3d_c     = (/nproma, nlev,   nblks_c  /)
    shape3d_chalf = (/nproma, nlevp1, nblks_c  /)
    shape4d_c     = (/nproma, nlev,   nblks_c, ntracer+ntracer_static/)

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
    cf_desc    = t_cf_var('normal_velocity', 'm s-1', 'velocity normal to edge')
    grib2_desc = t_grib2_var(0, 2, 197, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_prog_list, TRIM(vname_prefix)//'vn'//suffix, p_prog%vn,     &
      &           GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
      &           ldims=shape3d_e )

    ! w            p_prog%w(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('upward air velocity', 'm s-1', 'upward air velocity')
    grib2_desc = t_grib2_var(0, 2, 9, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_prog_list, TRIM(vname_prefix)//'w'//suffix, p_prog%w,   &
      &          GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc, &
      &          ldims=shape3d_chalf )

    ! rho          p_prog%rho(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('density', 'kg m-3', 'density')
    grib2_desc = t_grib2_var(0, 3, 10, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_prog_list, TRIM(vname_prefix)//'rho'//suffix, p_prog%rho,     &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,      &
      &           ldims=shape3d_c )

    ! theta_v      p_prog%theta_v(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('virtual_potential_temperature', 'K', &
      &                   'virtual potential temperature')
    grib2_desc = t_grib2_var(0, 0, 1, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_prog_list, TRIM(vname_prefix)//'theta_v'//suffix, p_prog%theta_v, &
      &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,          &
      &           ldims=shape3d_c )

    ! Initialize pointers that are not always allocated to NULL
    p_prog%exner      => NULL()
    p_prog%rhotheta_v => NULL()
    p_prog%tke        => NULL()
    p_prog%tracer     => NULL()

    IF (.NOT. l_extra_timelev) THEN
      ! exner        p_prog%exner(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('exner_pressure', '-', 'exner pressure')
      grib2_desc = t_grib2_var(0, 3, 195, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_prog_list, TRIM(vname_prefix)//'exner'//suffix, p_prog%exner, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,      &
        &           ldims=shape3d_c )

      ! rhotheta_v   p_prog%rhotheta_v(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('rho_virt_pot_temp', 'K', 'rho virt pot temp')
      grib2_desc = t_grib2_var(0, 19, 192, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_prog_list, TRIM(vname_prefix)//'rhotheta_v'//suffix, p_prog%rhotheta_v, &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,                &
        &           ldims=shape3d_c )

      ! Tracer array for (model) internal use

      ! tracer         p_prog%tracer(nproma,nlev,nblks_c,ntracer+ntracer_static)
      IF (ntracer > 0) THEN
        cf_desc    = t_cf_var('tracer', 'kg kg-1', 'tracer')
        grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_prog_list, 'tracer', p_prog%tracer,                       &
          &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
          &           ldims=shape4d_c ,                                           &
          &           lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)
      ENDIF



      IF (  iforcing == inwp  ) THEN

      ! Reference to individual tracer, for I/O

        ktracer=ntracer+ntracer_static
        ALLOCATE( p_prog%tracer_ptr(ktracer) )

           !QV
        CALL add_ref( p_prog_list, 'tracer',                                         &
                    & TRIM(vname_prefix)//'qv'//suffix, p_prog%tracer_ptr(iqv)%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var(TRIM(vname_prefix)//'qv',                             &
                    &  'kg kg-1','specific_humidity'),                               &
                    & t_grib2_var(0, 1, 0, ientr, GRID_REFERENCE, GRID_CELL),        &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=1,     &              ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.) )
           !QC
        CALL add_ref( p_prog_list, 'tracer',&
                    & TRIM(vname_prefix)//'qc'//suffix, p_prog%tracer_ptr(iqc)%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var(TRIM(vname_prefix)//'qc',                             &
                    &  'kg kg-1', 'specific_cloud_water_content'),                   &
                    & t_grib2_var(192, 201, 31, ientr, GRID_REFERENCE, GRID_CELL),   &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=1,     &              ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.) )
           !QI
        CALL add_ref( p_prog_list, 'tracer',                                         &
                    & TRIM(vname_prefix)//'qi'//suffix, p_prog%tracer_ptr(iqi)%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var(TRIM(vname_prefix)//'qi',                             &
                    &  'kg kg-1','specific_cloud_ice_content'),                      &
                    & t_grib2_var(192, 201, 33, ientr, GRID_REFERENCE, GRID_CELL),   &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=1,     &              ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.) )
           !QR
        CALL add_ref( p_prog_list, 'tracer',                                         &
                    & TRIM(vname_prefix)//'qr'//suffix, p_prog%tracer_ptr(iqr)%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var(TRIM(vname_prefix)//'qr',                             &
                    &  'kg kg-1','rain_mixing_ratio'),                               &
                    & t_grib2_var(0, 1, 24, ientr, GRID_REFERENCE, GRID_CELL),       &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=1,     &              ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.) )
           !QS
        CALL add_ref( p_prog_list, 'tracer',                                         &
                    & TRIM(vname_prefix)//'qs'//suffix, p_prog%tracer_ptr(iqs)%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
                    & t_cf_var(TRIM(vname_prefix)//'qs',                             &
                    &  'kg kg-1','snow_mixing_ratio'),                               &
                    & t_grib2_var(0, 1, 25, ientr, GRID_REFERENCE, GRID_CELL),       &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=1,     &              ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.) )


!!$           !QV
!!$        CALL add_tracer_ref( p_prog_list, 'tracer',                   &
!!$                    & TRIM(vname_prefix)//'qv'//suffix, iqv, p_prog%tracer_ptr,      &
!!$                    & t_cf_var(TRIM(vname_prefix)//'qv',                             &
!!$                    &  'kg kg-1','specific_humidity'),                               &
!!$                    & t_grib2_var(0, 1, 0, ientr, GRID_REFERENCE, GRID_CELL),        &
!!$                    & ldims=shape3d_c,                                               &
!!$                    & tlev_source=1,     &              ! output from nnow_rcf slice
!!$                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.) )
!!$
!!$           !QC
!!$        CALL add_tracer_ref( p_prog_list, 'tracer',&
!!$                    & TRIM(vname_prefix)//'qc'//suffix, iqc, p_prog%tracer_ptr,      &
!!$                    & t_cf_var(TRIM(vname_prefix)//'qc',                             &
!!$                    &  'kg kg-1', 'specific_cloud_water_content'),                   &
!!$                    & t_grib2_var(192, 201, 31, ientr, GRID_REFERENCE, GRID_CELL),   &
!!$                    & ldims=shape3d_c,                                               &
!!$                    & tlev_source=1,     &              ! output from nnow_rcf slice
!!$                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.) )
!!$           !QI
!!$        CALL add_tracer_ref( p_prog_list, 'tracer',                                  &
!!$                    & TRIM(vname_prefix)//'qi'//suffix, iqi, p_prog%tracer_ptr,      &
!!$                    & t_cf_var(TRIM(vname_prefix)//'qi',                             &
!!$                    &  'kg kg-1','specific_cloud_ice_content'),                      &
!!$                    & t_grib2_var(192, 201, 33, ientr, GRID_REFERENCE, GRID_CELL),   &
!!$                    & ldims=shape3d_c,                                               &
!!$                    & tlev_source=1,     &              ! output from nnow_rcf slice
!!$                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.) )
!!$           !QR
!!$        CALL add_tracer_ref( p_prog_list, 'tracer',                                  &
!!$                    & TRIM(vname_prefix)//'qr'//suffix, iqr, p_prog%tracer_ptr,      &
!!$                    & t_cf_var(TRIM(vname_prefix)//'qr',                             &
!!$                    &  'kg kg-1','rain_mixing_ratio'),                               &
!!$                    & t_grib2_var(0, 1, 24, ientr, GRID_REFERENCE, GRID_CELL),       &
!!$                    & ldims=shape3d_c,                                               &
!!$                    & tlev_source=1,     &              ! output from nnow_rcf slice
!!$                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.) )
!!$           !QS
!!$        CALL add_tracer_ref( p_prog_list, 'tracer',                                  &
!!$                    & TRIM(vname_prefix)//'qs'//suffix, iqs, p_prog%tracer_ptr,      &
!!$                    & t_cf_var(TRIM(vname_prefix)//'qs',                             &
!!$                    &  'kg kg-1','snow_mixing_ratio'),                               &
!!$                    & t_grib2_var(0, 1, 25, ientr, GRID_REFERENCE, GRID_CELL),       &
!!$                    & ldims=shape3d_c,                                               &
!!$                    & tlev_source=1,     &              ! output from nnow_rcf slice
!!$                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.) )


        IF( irad_o3 == 4 .OR. irad_o3 == 6 .OR. irad_o3 == 7 ) THEN
           !O3
          CALL add_ref( p_prog_list, 'tracer',                               &
            & TRIM(vname_prefix)//'O3'//suffix, p_prog%tracer_ptr(io3)%p_3d, &
            & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                          &
            & t_cf_var(TRIM(vname_prefix)//'O3',                             &
            &  'kg kg-1','ozone_mass_mixing_ratio'),                         &
            & t_grib2_var(0, 14, 1, ientr, GRID_REFERENCE, GRID_CELL),       &
            & ldims=shape3d_c,                                               &
            & tlev_source=0,     &              ! output from nnow_rcf slice
            & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.) )
        ENDIF

        ! tke            p_prog%tke(nproma,nlevp1,nblks_c)
        cf_desc    = t_cf_var('turbulent_kinetic_energy', 'm2 s-2', 'turbulent kinetic energy')
        grib2_desc = t_grib2_var(0, 19, 11, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_prog_list, TRIM(vname_prefix)//'tke'//suffix, p_prog%tke, &
          &           GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                       &
          &           cf_desc, grib2_desc, ldims=shape3d_chalf,                   &
          &           tlev_source=1 ) ! for output take field from nnow_rcf slice
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

!DR set loutput=.FALSE. for each individual field until a possible bug in the new 
!DR output scheme is fixed.
        CALL add_var_list_reference(p_tracer_list, from_info%name, &
          &                         from_var_list%p%name, loutput=.FALSE.)

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
    &                                 listname )
!
    TYPE(t_patch), TARGET, INTENT(IN) :: &  !< current patch
      &  p_patch

    TYPE(t_nh_diag),  INTENT(INOUT)   :: &  !< diagnostic state
      &  p_diag 

    TYPE(t_var_list), INTENT(INOUT)   :: &  !< diagnostic state list
      &  p_diag_list
    CHARACTER(len=*), INTENT(IN)      :: &  !< list name
      &  listname

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
      &        shape3d_c3(3)
 
    INTEGER :: ientr         !< "entropy" of horizontal slice
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

    IF (itime_scheme == 4 .OR. itime_scheme == 6) THEN
     n_timlevs = 2
    ELSE
     n_timlevs = 1
    ENDIF

    ientr = 16   ! "entropy" of horizontal slice

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
    cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'u-component of wind')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'u', p_diag%u,                                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE. )


    ! v           p_diag%v(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('northward_wind', 'm s-1', 'v-component of wind')
    grib2_desc = t_grib2_var(0, 2, 3, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'v', p_diag%v,                                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE. )

    ! vt           p_diag%vt(nproma,nlev,nblks_e)
    ! *** needs to be saved for restart ***
    cf_desc    = t_cf_var('tangential_wind', 'm s-1', 'tangential-component of wind')
    grib2_desc = t_grib2_var(0, 2, 198, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_diag_list, 'vt', p_diag%vt,                                 &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_e )


    ! omega_z      p_diag%omega_z(nproma,nlev,nblks_v)
    !
    cf_desc    = t_cf_var('vertical_vorticity', 'm s-1', 'vertical vorticity')
    grib2_desc = t_grib2_var(0, 2, 197, ientr, GRID_REFERENCE, GRID_VERTEX)
    CALL add_var( p_diag_list, 'omega_z', p_diag%omega_z,                       &
                & GRID_UNSTRUCTURED_VERT, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_v, lrestart=.FALSE., loutput=.FALSE. )

    ! omega_z_c    p_diag%omega_z_c(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('vertical_vorticity', 'm s-1', 'vertical vorticity')
    grib2_desc = t_grib2_var(0, 2, 197, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'omega_z', p_diag%omega_z_c,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE. )

    ! ddt_vn_phy   p_diag%ddt_vn_phy(nproma,nlev,nblks_e)
    ! *** needs to be saved for restart ***
    cf_desc    = t_cf_var('normal_wind_physical_tendency', 'm s-2',             &
      &                   'normal wind physical tendency')
    grib2_desc = t_grib2_var(0, 2, 199, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_diag_list, 'ddt_vn_phy', p_diag%ddt_vn_phy,                 &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_e )


    ! ddt_exner    p_diag%ddt_exner(nproma,nlev,nblks_c)
    ! *** needs to be saved for restart ***
    cf_desc    = t_cf_var('exner_pressure_tendency', 's-1', 'exner pressure tendency')
    grib2_desc = t_grib2_var(0, 3, 196, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'ddt_exner', p_diag%ddt_exner,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! ddt_exner_phy  p_diag%ddt_exner_phy(nproma,nlev,nblks_c)
    ! *** needs to be saved for restart ***
    cf_desc    = t_cf_var('exner_pressure_physical_tendency', 's-1',            &
      &                   'exner pressure physical tendency')
    grib2_desc = t_grib2_var(0, 3, 197, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'ddt_exner_phy', p_diag%ddt_exner_phy,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! ddt_temp_dyn  p_diag%ddt_temp_dyn(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('dynamical_temperature_tendency', 'K s-1',            &
      &                   'dynamical temperature tendency')
    grib2_desc = t_grib2_var(0, 3, 197, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'ddt_temp_dyn', p_diag%ddt_temp_dyn,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE. )

    ! exner_old    p_diag%exner_old(nproma,nlev,nblks_c)
    ! *** needs to be saved for restart ***
    cf_desc    = t_cf_var('old_exner_pressure', '-', 'old exner pressure')
    grib2_desc = t_grib2_var(0, 3, 196, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'exner_old', p_diag%exner_old,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )

    ! exner_dyn_incr    p_diag%exner_dyn_incr(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('exner_dynamics_increment', '-', 'exner dynamics increment')
    grib2_desc = t_grib2_var(0, 3, 196, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'exner_dyn_incr', p_diag%exner_dyn_incr,       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE. )

    ! pres_sfc     p_diag%pres_sfc(nproma,nblks_c)
    !
    cf_desc    = t_cf_var('surface_pressure', 'Pa', 'surface pressure')
    grib2_desc = t_grib2_var(0, 3, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'pres_sfc', p_diag%pres_sfc,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d_c, lrestart=.FALSE. )

    ! pres_sfc_old     p_diag%pres_sfc_old(nproma,nblks_c)
    !
    cf_desc    = t_cf_var('surface_pressure', 'Pa', 'surface pressure')
    grib2_desc = t_grib2_var(0, 3, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'pres_sfc_old', p_diag%pres_sfc_old,             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d_c, lrestart=.FALSE. )

    ! pres_sfc_s6avg     p_diag%pres_sfc_s6avg(nproma,nblks_c)
    !
    cf_desc    = t_cf_var('surface_pressure', 'Pa', 'surface pressure')
    grib2_desc = t_grib2_var(0, 3, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'pres_sfc_s6avg', p_diag%pres_sfc_s6avg,         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d_c, lrestart=.FALSE. )

    ! temp         p_diag%temp(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('temperature', 'K', 'temperature')
    grib2_desc = t_grib2_var(0, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'temp', p_diag%temp,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE. )


    ! tempv        p_diag%tempv(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('virtual_temperature', 'K', 'virtual temperature')
    grib2_desc = t_grib2_var(0, 0, 192, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'tempv', p_diag%tempv,                           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE. )


    ! temp_ifc     p_diag%temp_ifc(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('temperature', 'K', 'temperature at half level')
    grib2_desc = t_grib2_var(0, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'temp_ifc', p_diag%temp_ifc,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf, lrestart=.FALSE. )


    ! pres         p_diag%pres(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure')
    grib2_desc = t_grib2_var(0, 3, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'pres', p_diag%pres,                             &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE. )


    ! pres_ifc     p_diag%pres_ifc(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure at half level')
    grib2_desc = t_grib2_var(0, 3, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'pres_ifc', p_diag%pres_ifc,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf, lrestart=.FALSE. )


    ! dpres_mc     p_diag%dpres_mc(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('pressure_thickness', 'Pa', 'pressure thickness')
    grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'dpres_mc', p_diag%dpres_mc,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_c, lrestart=.FALSE. )


    ! div          p_diag%div(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('divergence', 's-1', 'divergence')
    grib2_desc = t_grib2_var( 0, 2, 13, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'div', p_diag%div,                               &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_c, lrestart=.FALSE. )


    ! mass_fl_e    p_diag%mass_fl_e(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('horizontal_mass_flux_at_edges', 'kg m-1 s-1',        &
       &         'horizontal mass flux at edges')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_diag_list, 'mass_fl_e', p_diag%mass_fl_e,                   &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_e, lrestart=.FALSE. )


    ! rho_ic       p_diag%rho_ic(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('density', 'kg m-3', 'density at half level')
    grib2_desc = t_grib2_var( 0, 3, 10, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'rho_ic', p_diag%rho_ic,                         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf, lrestart=.FALSE. )


    ! w_concorr_c  p_diag%w_concorr_c(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('contravariant_vertical_correction', 'm s-1',         &
      &                   'contravariant vertical correction')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'w_concorr_c', p_diag%w_concorr_c,               &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf, lrestart=.FALSE.  )


    ! e_kinh       p_diag%e_kinh(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('horizontal specific kinetic energy', 'm2 s-2',       &
      &                   'horizontal specific kinetic energy')
    grib2_desc = t_grib2_var( 0, 2, 196, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'e_kinh', p_diag%e_kinh,                         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_c, lrestart=.FALSE. )


    ! theta_v_ic   p_diag%theta_v_ic(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('virtual_potential_temperature_at_half_levels', 'K',&
      &                   'virtual_potential temperature at half levels')
    grib2_desc = t_grib2_var( 0, 0, 1, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_list, 'theta_v_ic', p_diag%theta_v_ic,               &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                & ldims=shape3d_chalf, lrestart=.FALSE. )


    IF (p_patch%cell_type == 3) THEN

      ! vn_ie        p_diag%vn_ie(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('normal_wind_at_half_level', 'm s-1',               &
        &                   'normal wind at half level')
      grib2_desc = t_grib2_var( 0, 2, 197, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'vn_ie', p_diag%vn_ie,                         &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_ehalf, lrestart=.FALSE. )


      ! ddt_vn_adv   p_diag%ddt_vn_adv(nproma,nlev,nblks_e,n_timlevs)
      ! *** needs to be saved for restart (TL nnow)***
      cf_desc    = t_cf_var('advective_normal_wind_tendency', 'm s-2',          &
        &                   'advective normal wind tendency')
      grib2_desc = t_grib2_var( 0, 2, 201, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_adv', p_diag%ddt_vn_adv,               &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape4d_entl ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ALLOCATE(p_diag%ddt_vn_adv_ptr(n_timlevs))
      DO jt =1,n_timlevs
        WRITE(suffix,'(".TL",i1)') jt
        CALL add_ref( p_diag_list, 'ddt_vn_adv',                                   &
                    & 'ddt_vn_adv'//suffix, p_diag%ddt_vn_adv_ptr(jt)%p_3d,        &
                    & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT,                        &
                    & t_cf_var('ddt_adv_vn'//suffix, 'm s-2',''),                  &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape3d_e )
      ENDDO


      ! ddt_w_adv    p_diag%ddt_w_adv(nproma,nlevp1,nblks_c,n_timlevs)
      ! *** needs to be saved for restart (TL nnow) ***
      cf_desc    = t_cf_var('advective_vertical_wind_tendency', 'm s-2',        &
        &                   'advective vertical wind tendency')
      grib2_desc = t_grib2_var( 0, 2, 202, ientr, GRID_REFERENCE, GRID_CELL)
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
                    & t_cf_var('ddt_adv_w'//suffix, 'm s-2',''),                    &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_chalf )
      ENDDO

      ! grf_tend_vn  p_diag%grf_tend_vn(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('normal_wind_tendency', 'm s-2',                    &
        &                   'normal wind tendency (grid refinement)')
      grib2_desc = t_grib2_var( 0, 2, 203, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'grf_tend_vn', p_diag%grf_tend_vn,             &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_e, lrestart=.FALSE. )


      ! grf_tend_w  p_diag%grf_tend_w(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('vertical_wind_tendency', 'm s-2',                  &
        &                   'vertical wind tendency (grid refinement)')
      grib2_desc = t_grib2_var( 0, 2, 204, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'grf_tend_w', p_diag%grf_tend_w,               &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf, lrestart=.FALSE. )

      ! grf_tend_rho   p_diag%grf_tend_rho(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('density_tendency', 'kg m-3 s-1',                   &
        &                   'density tendency (grid refinement)')
      grib2_desc = t_grib2_var( 0, 3, 198, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'grf_tend_rho', p_diag%grf_tend_rho,           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_c, lrestart=.FALSE. )

      ! grf_tend_thv   p_diag%grf_tend_thv(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('virtual_potential_temperature_tendency', 'K s-1',  &
        &                   'virtual potential temperature tendency (grid refinement)')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'grf_tend_thv', p_diag%grf_tend_thv,           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_c, lrestart=.FALSE. )


      ! Storage fields for vertical nesting; the middle index (2) addresses 
      ! the field and its temporal tendency

      ! dvn_ie_int   p_diag%dvn_ie_int(nproma,nblks_e)
      !
      cf_desc    = t_cf_var('normal_velocity_parent_interface_level', 'm s-1',  &
        &                   'normal velocity at parent interface level')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'dvn_ie_int', p_diag%dvn_ie_int,               &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_e, lrestart=.FALSE. )


      ! dvn_ie_ubc   p_diag%dvn_ie_ubc(nproma,nblks_e)
      !
      cf_desc    = t_cf_var('normal_velocity_child_upper_boundary', 'm s-1',    &
        &                   'normal velocity at child upper boundary')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'dvn_ie_ubc', p_diag%dvn_ie_ubc,               &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_e, lrestart=.FALSE. )


      ! drho_ic_int  p_diag%drho_ic_int(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('rho_at_parent_interface_level', 'kg m-3',          &
        &                   'rho at parent interface level')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'drho_ic_int', p_diag%drho_ic_int,             &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_c, lrestart=.FALSE. )


      ! drho_ic_ubc  p_diag%drho_ic_ubc(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('density_at_child_upper_boundary', 'kg m-3',        &
        &                   'density at child upper boundary')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'drho_ic_ubc', p_diag%drho_ic_ubc,             &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_c, lrestart=.FALSE. )


      ! dtheta_v_ic_int    p_diag%dtheta_v_ic_int(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('theta_at_parent_interface_level', 'K',             &
        &                   'potential temperature at parent interface level')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'dtheta_v_ic_int', p_diag%dtheta_v_ic_int,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_c, lrestart=.FALSE. )


      ! dtheta_v_ic_ubc    p_diag%dtheta_v_ic_ubc(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('theta_at_child_upper_boundary', 'K',               &
        &                   'potential temperature at child upper boundary')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'dtheta_v_ic_ubc', p_diag%dtheta_v_ic_ubc,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_c, lrestart=.FALSE. )


      ! dw_int       p_diag%dw_int(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('w_at_parent_interface_level', 'm s-1',             &
        &                   'vertical velocity at parent interface level')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'dw_int', p_diag%dw_int,                       &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape2d_c, lrestart=.FALSE. )


      ! dw_ubc       p_diag%dw_ubc(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('w at child upper boundary', 'm s-1',               &
        &                   'vertical velocity at child upper boundary')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'dw_ubc', p_diag%dw_ubc,                       &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                  & ldims=shape2d_c, lrestart=.FALSE. )


      ! q_int        p_diag%q_int(nproma,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('q_at_parent_interface_level', 'kg kg-1',           &
        &                   'q at parent interface level')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
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
                    & t_cf_var('q_int'//ctrc, 'kg kg-1',''),                        &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape2d_c, lrestart=.FALSE. )
      ENDDO


      ! q_ubc        p_diag%q_ubc(nproma,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('q_at_child_upper_boundary', 'kg kg-1',             &
        &                   'q at child upper boundary')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
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
                    & t_cf_var('q_ubc'//ctrc, 'kg kg-1',''),                        &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape2d_c, lrestart=.FALSE. )
      ENDDO

    ELSE IF (p_patch%cell_type == 6) THEN
 
      ! e_kin        p_diag%e_kin(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('specific_kinetic_energy', 'm2 s-2', 'specific kinetic energy')
      grib2_desc = t_grib2_var(0, 2, 196, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'e_kin', p_diag%e_kin,                           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                  & ldims=shape3d_c )

      ! theta_v_impl   p_diag%theta_v_impl(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('(nnow+nnew)/2_from_impl._vert._adv._of_theta_v', 'K', &
        &                   '(nnow+nnew)/2 from impl. vert. adv. of theta_v')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'theta_v_impl', p_diag%theta_v_impl,           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c )

      ! theta_v_ave   p_diag%theta_v_ave(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('time average from horiz. adv. of theta_v', 'K', &
        &                   'time average from horiz. adv. of theta_v')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'theta_v_ave', p_diag%theta_v_ave,           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,&
                  & ldims=shape3d_c )


      ! horpgrad     p_diag%horpgrad(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('covariant_horizontal_pressure_gradient', 'Pa m-1', &
        &                   'covariant horizontal pressure gradient')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'horpgrad', p_diag%horpgrad,                   &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e )


      ! vn_cov       p_diag%vn_cov(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('covariant_normal_wind', 'm s-1',                   &
        &                   'covariant normal wind')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'vn_cov', p_diag%vn_cov,                       &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e )


      ! w_cov        p_diag%w_cov(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('covariant_vertical_wind', 'm s-1',                 &
        &                   'covariant vertical wind at half level')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'w_cov', p_diag%w_cov,                         &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  & 
                  & ldims=shape3d_chalf )

      ! rho_e        p_diag%rho_e(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('density_at_edges_nnow', 'm s-1',                   &
        &                   'density_at_edges_nnow')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'rho_e', p_diag%rho_e,                         &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e )

      ! rho_star_e   p_diag%rho_star_e(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('density_at_edges_estim_step', 'm s-1',              &
        &                   'density_at_edges_estim_step')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'rho_star_e', p_diag%rho_star_e,                &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,   &
                  & ldims=shape3d_e )

      ! omega_t_con  p_diag%omega_t_con(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('tangential_horiz._contravariant_vorticity', 's-1', &
        &                   'tangential horiz. contravariant vorticity')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'omega_t_con', p_diag%omega_t_con,             &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_ehalf )


      ! omega_x      p_diag%omega_x(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('zonal_vorticity', 's-1', 'zonal vorticity')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'omega_x', p_diag%omega_x,                     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c )


      ! omega_y      p_diag%omega_y(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('meridional_vorticity', 's-1', 'meridional vorticity')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'omega_y', p_diag%omega_y,                     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c )


      ! omega_z_con  p_diag%omega_z_con(nproma,nlev,nblks_v)
      !
      cf_desc    = t_cf_var('vertical_vorticity', 's-1',                        &
        &                   'contravariant vertical vorticity')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_VERTEX)
      CALL add_var( p_diag_list, 'omega_z_con', p_diag%omega_z_con,             &
                  & GRID_UNSTRUCTURED_VERT, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_v )


      ! ddt_vn       p_diag%ddt_vn(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('normal_wind_tendency', 'm s-2', 'normal wind tendency')
      grib2_desc = t_grib2_var(0, 2, 198, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn', p_diag%ddt_vn,                         &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                  & ldims=shape3d_e )


      ! ddt_vn_vort  p_diag%ddt_vn_vort(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('normal_wind_tendency', 'm s-2',                    &
        &                   'normal wind tendency from vorticity flux term')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_vort', p_diag%ddt_vn_vort,             &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e )


      ! ddt_w        p_diag%ddt_w(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('vertical_wind_tendency', 'm s-2', 'vertical wind tendency')
      grib2_desc = t_grib2_var(0, 2, 200, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_w', p_diag%ddt_w,                           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                  & ldims=shape3d_chalf )

 
      ! ddt_w_vort   p_diag%ddt_w_vort(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('vert._wind_tendency', 'm s-2',                     &
        &                   'vert. wind tendency from vorticity flux term')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_w_vort', p_diag%ddt_w_vort,               &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf )

      ! ddt_w_phy   p_diag%ddt_w_phy(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('vert_wind_tendency_phy', 'm s-2',                   &
        &                   'vertical wind tendency from physics')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_w_phy', p_diag%ddt_w_phy,                  &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,   &
                  & ldims=shape3d_chalf )

    ENDIF


    !
    ! tracers
    !
!    IF ( ltransport ) THEN
    IF ( ntracer >0 ) THEN

      ! grf_tend_tracer   p_diag%grf_tend_tracer(nproma,nlev,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('tracer_tendency', 'kg kg-1 s-1',                   &
        &                   'tracer_tendency for grid refinement')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_list, 'grf_tend_tracer', p_diag%grf_tend_tracer,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape4d_c ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

      ALLOCATE(p_diag%ddt_grf_trc_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'grf_tend_tracer',                              &
                    & 'ddt_grf_q'//ctrc, p_diag%ddt_grf_trc_ptr(jt)%p_3d,             &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT,                        &
                    & t_cf_var('ddt_grf_q'//ctrc, 'kg kg-1 s**-1',''),             &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape3d_c, lrestart=.FALSE. )
      ENDDO


      ! hfl_tracer   p_diag%hfl_tracer(nproma,nlev,nblks_e,ntracer)
      !
      cf_desc    = t_cf_var('horizontal tracer flux', 'kg m-1 s-1',               &
        &                   'horizontal tracer flux')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
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
                    & t_cf_var('hfl_q'//ctrc, 'kg m-1 s-1',''),                     &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_e, lrestart=.FALSE. )
      ENDDO


      ! vfl_tracer   p_diag%vfl_tracer(nproma,nlevp1,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('vertical_tracer_flux', 'kg m-1 s-1',               &
        &                   'vertical tracer flux')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
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
                    & t_cf_var('vfl_q'//ctrc, 'kg m-1 s-1',''),                    &
                    & t_grib2_var(255,255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_chalf, lrestart=.FALSE. )
      ENDDO


      ! ddt_tracer_adv   p_diag%ddt_tracer_adv(nproma,nlev,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('advective tracer tendency', 'kg kg-1 s-1',         &
        &                   'advective tracer tendency')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
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
                    & t_cf_var('ddt_adv_q'//ctrc, 'kg kg-1 s-1',''),               &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL),&
                    & ldims=shape3d_c, lrestart=.FALSE. )
      ENDDO


      ! tracer_vi(nproma,nblks_c,3), only Q1, Q2, Q3
      cf_desc    = t_cf_var('tracer_vi', '', 'tracer_vi')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
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
          &           cf_desc, grib2_desc, ldims=shape3d_c3, lrestart=.FALSE.)
      ENDDO

      ! tracer_vi_avg(nproma,nblks_c,3), only Q1, Q2, Q3
      cf_desc    = t_cf_var('tracer_vi_avg', '', 'tracer_vi_avg')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
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
          &           cf_desc, grib2_desc, ldims=shape3d_c3, lrestart=.FALSE. )
      ENDDO
    ENDIF



    IF( lwrite_extra) THEN
      WRITE(0,*)'inextra_2d=',inextra_2d

      IF(inextra_2d > 0) THEN

        ! extra_2d   p_diag%extra_2d(nproma,nblks_c,inextra_2d)
        !
        cf_desc    = t_cf_var('extra_field_2D', '-', 'extra field 2D')
        grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_diag_list, 'extra_2d', p_diag%extra_2d,                   &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc, &
                    & ldims=shape2d_extra, lrestart=.FALSE. )
      ENDIF

      IF(inextra_3d > 0) THEN

        ! extra_3d   p_diag%extra_3d(nproma,nlev,nblks_c,inextra_3d)
        !
        cf_desc    = t_cf_var('extra_fields_3D', '-', 'extra fields 3D')
        grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_diag_list, 'extra_3d', p_diag%extra_3d,                   &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                    & ldims=shape3d_extra,                                        &
                    & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

      ENDIF
    ENDIF

  END SUBROUTINE new_nh_state_diag_list




  !-------------------------------------------------------------------------
  !>
  !! Allocation of components of diagnostic z-state.
  !!
  !! Initialization of components with zero.
  !!
  !! @par Revision History
  !! Initial release by Daniel Reinert, DWD (2011-09-05)
  !!
  SUBROUTINE new_nh_state_diag_z_list ( p_patch, p_diag_z, p_diag_z_list,  &
    &                                   listname )
!
    TYPE(t_patch), TARGET, INTENT(IN) :: &  !< current patch
      &  p_patch

    TYPE(t_nh_diag_pz), INTENT(INOUT) :: &  !< diagnostic state
      &  p_diag_z 

    TYPE(t_var_list), INTENT(INOUT)   :: &  !< diagnostic state list
      &  p_diag_z_list

    CHARACTER(len=*), INTENT(IN)      :: &  !< list name
      &  listname

    CHARACTER(len=max_char_length)    :: vname_prefix

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c       !< number of cell blocks to allocate

    INTEGER :: nlev          !< number of vertical levels

    INTEGER :: jg            !< patch ID
 
    INTEGER :: shape3d_c(3), shape4d_c(4)
 
    INTEGER :: ientr         !< "entropy" of horizontal slice

    INTEGER :: kcloud, ktracer

    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c

    ! determine patch ID
    jg = p_patch%id

    ! number of vertical levels
    nlev   = nh_pzlev_config(jg)%nzlev

    ientr = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape3d_c     = (/nproma, nlev, nblks_c /)
    shape4d_c     = (/nproma, nlev, nblks_c, ntracer+ntracer_static/)

    kcloud= 4

    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_diag_z_list, TRIM(listname), patch_id=p_patch%id, level_type=3 )
    CALL default_var_list_settings( p_diag_z_list,             &
                                  & lrestart=.FALSE.,          &
                                  & restart_type=FILETYPE_NC2  )


    ! u           p_diag_z%u(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'u-component of wind')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_z_list, 'u', p_diag_z%u,                              &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE, cf_desc, grib2_desc, &
                & ldims=shape3d_c )


    ! v           p_diag_z%v(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('northward_wind', 'm s-1', 'v-component of wind')
    grib2_desc = t_grib2_var(0, 2, 3, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_z_list, 'v', p_diag_z%v,                              &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE, cf_desc, grib2_desc, &
                & ldims=shape3d_c )


    ! pres         p_diag_z%pres(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure')
    grib2_desc = t_grib2_var(0, 3, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_z_list, 'pres', p_diag_z%pres,                        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE, cf_desc, grib2_desc, &
                & ldims=shape3d_c )


    ! temp         p_diag_z%temp(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('temperature', 'K', 'temperature')
    grib2_desc = t_grib2_var(0, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_z_list, 'temp', p_diag_z%temp,                        &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE, cf_desc, grib2_desc, &
                & ldims=shape3d_c )



    ! Tracer array for (model) internal use

    ! tracer        p_diag_z%tracer(nproma,nlev,nblks_c,ntracer+ntracer_static)
    IF ( ntracer > 0 ) THEN

      cf_desc    = t_cf_var('tracer', 'kg kg-1', 'tracer')
      grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_z_list, 'tracer', p_diag_z%tracer,                   &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE, cf_desc, grib2_desc,&
        &           ldims=shape4d_c ,                                           &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)


      IF (  iforcing == inwp  ) THEN

      ! Reference to individual tracer, for I/O

        ktracer=ntracer+ntracer_static
        ALLOCATE( p_diag_z%tracer_ptr(ktracer) )
        vname_prefix='tracer_'

           !QV
        CALL add_ref( p_diag_z_list, 'tracer',                                 &
                    & TRIM(vname_prefix)//'qv', p_diag_z%tracer_ptr(iqv)%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,                  &
                    & t_cf_var(TRIM(vname_prefix)//'qv',                       &
                    &  'kg kg-1','specific_humidity'),                         &
                    & t_grib2_var(0, 1, 0, ientr, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d_c)
           !QC
        CALL add_ref( p_diag_z_list, 'tracer',                                     &
                    & TRIM(vname_prefix)//'qc', p_diag_z%tracer_ptr(iqc)%p_3d,     &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,                      &
                    & t_cf_var(TRIM(vname_prefix)//'qc',                           &
                    &  'kg kg-1', 'specific_cloud_water_content'),                 &
                    & t_grib2_var(192, 201, 31, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)
           !QI
        CALL add_ref( p_diag_z_list, 'tracer',                                     &
                    & TRIM(vname_prefix)//'qi', p_diag_z%tracer_ptr(iqi)%p_3d,     &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,                      &
                    & t_cf_var(TRIM(vname_prefix)//'qi',                           & 
                    &  'kg kg-1','specific_cloud_ice_content'),                    &
                    & t_grib2_var(192, 201, 33, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)
           !QR
        CALL add_ref( p_diag_z_list, 'tracer',                                  &
                    & TRIM(vname_prefix)//'qr', p_diag_z%tracer_ptr(iqr)%p_3d,  &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,                   &
                    & t_cf_var(TRIM(vname_prefix)//'qr',                        &       
                    &  'kg kg-1','rain_mixing_ratio'),                          &
                    & t_grib2_var(0, 1, 24, ientr, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d_c)
           !QS
        CALL add_ref( p_diag_z_list, 'tracer',                                 &
                    & TRIM(vname_prefix)//'qs', p_diag_z%tracer_ptr(iqs)%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,                  &
                    & t_cf_var(TRIM(vname_prefix)//'qs',                       &
                    &  'kg kg-1','snow_mixing_ratio'),                         &
                    & t_grib2_var(0, 1, 25, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)

        IF(irad_o3 == 3) THEN
           !O3
          CALL add_ref( p_diag_z_list, 'tracer',                       &
            & TRIM(vname_prefix)//'O3', p_diag_z%tracer_ptr(io3)%p_3d, &
            & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,                  &
            & t_cf_var(TRIM(vname_prefix)//'O3',                       &
            &  'kg kg-1','ozone_mass_mixing_ratio'),                   &
            & t_grib2_var(0, 14, 1, ientr, GRID_REFERENCE, GRID_CELL), &
            & ldims=shape3d_c)
        ENDIF




        ! &      p_diag_z%tot_cld(nproma,nlev,nblks_c,4)
        cf_desc    = t_cf_var('tot_cld', 'kg kg-1','total cloud variables (cc,qv,qc,qi)')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_diag_z_list, 'tot_cld', p_diag_z%tot_cld,                  &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE, cf_desc, grib2_desc, &
                    &                         ldims=(/nproma,nlev,nblks_c,kcloud/),&
                    & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)


        ALLOCATE( p_diag_z%tot_ptr(kcloud))
        vname_prefix='tot_'

           !CC
        CALL add_ref( p_diag_z_list, 'tot_cld',                                     &
                    & TRIM(vname_prefix)//'cc', p_diag_z%tot_ptr(1)%p_3d,           &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,                       &
                    & t_cf_var(TRIM(vname_prefix)//'cc', '','total_cloud_cover'),   &
                    & t_grib2_var(192, 128, 164, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)

           !QV
        CALL add_ref( p_diag_z_list, 'tot_cld',                                        &
                    & TRIM(vname_prefix)//'qv', p_diag_z%tot_ptr(2)%p_3d,              &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,                          &
                    & t_cf_var(TRIM(vname_prefix)//'qv', '','total_specific_humidity'),&
                    & t_grib2_var( 0, 1, 0, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)

           !QC
        CALL add_ref( p_diag_z_list, 'tot_cld',                                     &
                    & TRIM(vname_prefix)//'qc', p_diag_z%tot_ptr(3)%p_3d,           &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,                       &
                    & t_cf_var(TRIM(vname_prefix)//'qc', '',                        &
                    & 'total_specific_cloud_water_content'),                        &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)

           !QI
        CALL add_ref( p_diag_z_list, 'tot_cld',                                     &
                    & TRIM(vname_prefix)//'qi', p_diag_z%tot_ptr(4)%p_3d,           &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_ALTITUDE,                       &
                    & t_cf_var(TRIM(vname_prefix)//'qi', '',                        &
                    & 'total_specific_cloud_ice_content'),                          &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)


      ENDIF ! iforcing
    ENDIF ! ntracer



  END SUBROUTINE new_nh_state_diag_z_list



  !-------------------------------------------------------------------------
  !>
  !! Allocation of components of diagnostic p-state.
  !!
  !! Initialization of components with zero.
  !!
  !! @par Revision History
  !! Initial release by Daniel Reinert, DWD (2011-09-05)
  !!
  SUBROUTINE new_nh_state_diag_p_list ( p_patch, p_diag_p, p_diag_p_list,  &
    &                                 listname )
!
    TYPE(t_patch), TARGET, INTENT(IN) :: &  !< current patch
      &  p_patch

    TYPE(t_nh_diag_pz), INTENT(INOUT) :: &  !< diagnostic state
      &  p_diag_p 

    TYPE(t_var_list), INTENT(INOUT)   :: &  !< diagnostic state list
      &  p_diag_p_list

    CHARACTER(len=*), INTENT(IN)      :: &  !< list name
      &  listname

    CHARACTER(len=max_char_length)    :: vname_prefix

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c       !< number of cell blocks to allocate

    INTEGER :: nlev          !< number of vertical levels

    INTEGER :: jg            !< patch ID

    INTEGER :: shape3d_c(3), shape4d_c(4)
 
    INTEGER :: ientr         !< "entropy" of horizontal slice

    INTEGER :: kcloud, ktracer

    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c

    ! determine patch ID
    jg = p_patch%id

    ! number of vertical levels
    nlev   =  nh_pzlev_config(jg)%nplev

    ientr = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape3d_c     = (/nproma, nlev, nblks_c /)
    shape4d_c     = (/nproma, nlev, nblks_c, ntracer+ntracer_static/)

    kcloud= 4


    !
    ! Register a field list and apply default settings
    !
    CALL new_var_list( p_diag_p_list, TRIM(listname), patch_id=p_patch%id, level_type=2 )
    CALL default_var_list_settings( p_diag_p_list,             &
                                  & lrestart=.FALSE.,          &
                                  & restart_type=FILETYPE_NC2  )


    ! u           p_diag_p%u(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'u-component of wind')
    grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_p_list, 'u', p_diag_p%u,                                 &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! v           p_diag_p%v(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('northward_wind', 'm s-1', 'v-component of wind')
    grib2_desc = t_grib2_var(0, 2, 3, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_p_list, 'v', p_diag_p%v,                                 &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! pres         p_diag_p%pres(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('pressure', 'Pa', 'pressure')
    grib2_desc = t_grib2_var(0, 3, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_p_list, 'pres', p_diag_p%pres,                           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! temp         p_diag_p%temp(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('temperature', 'K', 'temperature')
    grib2_desc = t_grib2_var(0, 0, 0, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_p_list, 'temp', p_diag_p%temp,                           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! z           p_diag_p%geopot(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('z', 'm2 s-2', 'geopotential')
    grib2_desc = t_grib2_var(0, 3, 4, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_diag_p_list, 'z', p_diag_p%geopot,                            &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )



    ! Tracer array for (model) internal use

    ! tracer         p_diag_p%tracer(nproma,nlev,nblks_c,ntracer+ntracer_static)
    IF ( ntracer > 0 ) THEN

      cf_desc    = t_cf_var('tracer', 'kg kg-1', 'tracer')
      grib2_desc = t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_diag_p_list, 'tracer', p_diag_p%tracer,                   &
        &           GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE, cf_desc, grib2_desc,&
        &           ldims=shape4d_c ,                                           &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)


      IF (  iforcing == inwp  ) THEN

      ! Reference to individual tracer, for I/O

        ktracer=ntracer+ntracer_static
        ALLOCATE( p_diag_p%tracer_ptr(ktracer) )
        vname_prefix='tracer_'

           !QV
        CALL add_ref( p_diag_p_list, 'tracer',                                 &
                    & TRIM(vname_prefix)//'qv', p_diag_p%tracer_ptr(iqv)%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE,                  &
                    & t_cf_var(TRIM(vname_prefix)//'qv',                       &
                    &  'kg kg-1','specific_humidity'),                         &
                    & t_grib2_var(0, 1, 0, ientr, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d_c)
           !QC
        CALL add_ref( p_diag_p_list, 'tracer',                                     &
                    & TRIM(vname_prefix)//'qc', p_diag_p%tracer_ptr(iqc)%p_3d,     &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE,                      &
                    & t_cf_var(TRIM(vname_prefix)//'qc',                           &
                    &  'kg kg-1', 'specific_cloud_water_content'),                 &
                    & t_grib2_var(192, 201, 31, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)
           !QI
        CALL add_ref( p_diag_p_list, 'tracer',                                     &
                    & TRIM(vname_prefix)//'qi', p_diag_p%tracer_ptr(iqi)%p_3d,     &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE,                      &
                    & t_cf_var(TRIM(vname_prefix)//'qi',                           & 
                    &  'kg kg-1','specific_cloud_ice_content'),                    &
                    & t_grib2_var(192, 201, 33, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)
           !QR
        CALL add_ref( p_diag_p_list, 'tracer',                                  &
                    & TRIM(vname_prefix)//'qr', p_diag_p%tracer_ptr(iqr)%p_3d,  &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE,                   &
                    & t_cf_var(TRIM(vname_prefix)//'qr',                        &       
                    &  'kg kg-1','rain_mixing_ratio'),                          &
                    & t_grib2_var(0, 1, 24, ientr, GRID_REFERENCE, GRID_CELL),  &
                    & ldims=shape3d_c)
           !QS
        CALL add_ref( p_diag_p_list, 'tracer',                                 &
                    & TRIM(vname_prefix)//'qs', p_diag_p%tracer_ptr(iqs)%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE,                  &
                    & t_cf_var(TRIM(vname_prefix)//'qs',                       &
                    &  'kg kg-1','snow_mixing_ratio'),                         &
                    & t_grib2_var(0, 1, 25, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)

        IF(irad_o3 == 3) THEN
           !O3
          CALL add_ref( p_diag_p_list, 'tracer',                       &
            & TRIM(vname_prefix)//'O3', p_diag_p%tracer_ptr(io3)%p_3d, &
            & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE,                  &
            & t_cf_var(TRIM(vname_prefix)//'O3',                       &
            &  'kg kg-1','ozone_mass_mixing_ratio'),                   &
            & t_grib2_var(0, 14, 1, ientr, GRID_REFERENCE, GRID_CELL), &
            & ldims=shape3d_c)
        ENDIF



        ! &      p_diag_p%tot_cld(nproma,nlev,nblks_c,4)
        cf_desc    = t_cf_var('tot_cld', 'kg kg-1','total cloud variables (cc,qv,qc,qi)')
        grib2_desc = t_grib2_var(0, 2, 2, ientr, GRID_REFERENCE, GRID_CELL)
        CALL add_var( p_diag_p_list, 'tot_cld', p_diag_p%tot_cld,                  &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE, cf_desc, grib2_desc, &
                    &                         ldims=(/nproma,nlev,nblks_c,kcloud/),&
                    & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)


        ALLOCATE( p_diag_p%tot_ptr(kcloud))
        vname_prefix='tot_'

           !CC
        CALL add_ref( p_diag_p_list, 'tot_cld',                                     &
                    & TRIM(vname_prefix)//'cc', p_diag_p%tot_ptr(1)%p_3d,           &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE,                       &
                    & t_cf_var(TRIM(vname_prefix)//'cc', '','total_cloud_cover'),   &
                    & t_grib2_var(192, 128, 164, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)

           !QV
        CALL add_ref( p_diag_p_list, 'tot_cld',                                        &
                    & TRIM(vname_prefix)//'qv', p_diag_p%tot_ptr(2)%p_3d,              &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE,                          &
                    & t_cf_var(TRIM(vname_prefix)//'qv', '','total_specific_humidity'),&
                    & t_grib2_var( 0, 1, 0, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)

           !QC
        CALL add_ref( p_diag_p_list, 'tot_cld',                                     &
                    & TRIM(vname_prefix)//'qc', p_diag_p%tot_ptr(3)%p_3d,           &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE,                       &
                    & t_cf_var(TRIM(vname_prefix)//'qc', '',                        &
                    & 'total_specific_cloud_water_content'),                        &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)

           !QI
        CALL add_ref( p_diag_p_list, 'tot_cld',                                     &
                    & TRIM(vname_prefix)//'qi', p_diag_p%tot_ptr(4)%p_3d,           &
                    & GRID_UNSTRUCTURED_CELL, ZAXIS_PRESSURE,                       &
                    & t_cf_var(TRIM(vname_prefix)//'qi', '',                        &
                    & 'total_specific_cloud_ice_content'),                          &
                    & t_grib2_var(255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL), &
                    & ldims=shape3d_c)


      ENDIF ! iforcing
    ENDIF ! ktracer


  END SUBROUTINE new_nh_state_diag_p_list


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
      &        shape3d_esquared(4) 
    INTEGER :: ientr         !< "entropy" of horizontal slice
    INTEGER :: ist
    !--------------------------------------------------------------

    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nblks_v = p_patch%nblks_v

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ientr = 16   ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d_c        = (/nproma,          nblks_c    /)     
    shape2d_ccubed   = (/nproma, 3      , nblks_c    /)     
    shape2d_ecubed   = (/nproma, 3      , nblks_e    /)     
    shape3d_c        = (/nproma, nlev   , nblks_c    /)     
    shape3d_chalf    = (/nproma, nlevp1 , nblks_c    /)      
    shape3d_e        = (/nproma, nlev   , nblks_e    /)     
    shape3d_ehalf    = (/nproma, nlevp1 , nblks_e    /)     
    shape3d_esquared = (/2     , nproma , nlev   , nblks_e /)
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
      &                   'geometric height at half level center')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'z_ifc', p_metrics%z_ifc,                     &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! geometric height at full levels
    ! z_mc         p_metrics%z_mc(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('geometric_height_at_full_level_center', 'm',         &
      &                   'geometric height at full level center')
    grib2_desc = t_grib2_var( 0, 3, 6, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'z_mc', p_metrics%z_mc,                       &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! slope of the terrain in normal direction (half level)
    ! ddxn_z_half  p_metrics%ddxn_z_half(nproma,nlevp1,nblks_e)
    !
    cf_desc    = t_cf_var('terrain_slope_in_normal_direction', '-',             &
      &                   'terrain slope in normal direction')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_metrics_list, 'ddxn_z_half', p_metrics%ddxn_z_half,         &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_ehalf )


    ! slope of the terrain in normal direction (full level)
    ! ddxn_z_full  p_metrics%ddxn_z_full(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('terrain_slope_in_normal_direction', '-',             &
      &                   'terrain slope in normal direction')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_metrics_list, 'ddxn_z_full', p_metrics%ddxn_z_full,         &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_e )


    ! slope of the terrain in tangential direction (half level)
    ! ddxt_z_half  p_metrics%ddxt_z_half(nproma,nlevp1,nblks_e)
    !
    cf_desc    = t_cf_var('terrain_slope_in_tangential_direction', '-',         &
      &                   'terrain slope in tangential direction')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_metrics_list, 'ddxt_z_half', p_metrics%ddxt_z_half,         &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_ehalf )


    ! functional determinant of the metrics [sqrt(gamma)]
    ! ddqz_z_full  p_metrics%ddqz_z_full(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('metrics_functional_determinant', '-',                &
      &                   'metrics functional determinant')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'ddqz_z_full', p_metrics%ddqz_z_full,         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_c )


    ! functional determinant of the metrics [sqrt(gamma)]
    ! ddqz_z_full_e  p_metrics%ddqz_z_full_e(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('metrics_functional_determinant', '-',                &
      &                   'metrics functional determinant (edge)')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
    CALL add_var( p_metrics_list, 'ddqz_z_full_e', p_metrics%ddqz_z_full_e,     &
                & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,         &
                & ldims=shape3d_e )


    ! functional determinant of the metrics [sqrt(gamma)]
    ! ddqz_z_half  p_metrics%ddqz_z_half(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('metrics_functional_determinant', '-',                &
      &                   'metrics functional determinant')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'ddqz_z_half', p_metrics%ddqz_z_half,         &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! geopotential at full level cell center
    ! geopot       p_metrics%geopot(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('geopotential', 'm2 s-2',                             &
       &                  'geopotential at full level cell centre')
    grib2_desc = t_grib2_var( 0, 3, 4, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'geopot', p_metrics%geopot,                   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! geopotential above groundlevel at cell center
    ! geopot_agl   p_metrics%geopot_agl(nproma,nlev  ,nblks_c)
    !
    cf_desc    = t_cf_var('geopotential', 'm2 s-2',                             &
      &                   'geopotential above groundlevel at cell center')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'geopot_agl', p_metrics%geopot_agl,           &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_c )


    ! geopotential above groundlevel at cell center
    ! geopot_agl_ifc  p_metrics%geopot_agl_ifc(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('geopotential', 'm2 s-2',                             &
      &                   'geopotential above groundlevel at cell center')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'geopot_agl_ifc', p_metrics%geopot_agl_ifc,   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                & ldims=shape3d_chalf )


    ! geopotential at cell center
    ! dgeopot_mc   p_metrics%dgeopot_mc(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('geopotential', 'm2 s-2',                             &
      &                   'geopotential at cell center')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
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
      &                   'Explicit weight in vertical wind solver')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
    CALL add_var( p_metrics_list, 'vwind_expl_wgt', p_metrics%vwind_expl_wgt,   &
                & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,        &
                & ldims=shape2d_c )


    ! Implicit weight in vertical wind solver
    ! vwind_impl_wgt  p_metrics%vwind_impl_wgt(nproma,nblks_c)
    !
    cf_desc    = t_cf_var('Implicit_weight_in_vertical_wind_solver', '-',       &
      &                   'Implicit weight in vertical wind solver')
    grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
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
      &                   'weighting factor for interpolation from full to half levels')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'wgtfac_c', p_metrics%wgtfac_c,             &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf )


      ! weighting factor for interpolation from full to half levels
      ! wgtfac_e     p_metrics%wgtfac_e(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                   'weighting factor for interpolation from full to half levels')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'wgtfac_e', p_metrics%wgtfac_e,             &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_ehalf )


      ! weighting factor for quadratic interpolation to surface
      ! wgtfacq_c    p_metrics%wgtfacq_c(nproma,3,nblks_c)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                   'weighting factor for quadratic interpolation to surface')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'wgtfacq_c', p_metrics%wgtfacq_c,           &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape2d_ccubed, loutput=.FALSE. )


      ! weighting factor for quadratic interpolation to surface
      ! wgtfacq_e    p_metrics%wgtfacq_e(nproma,3,nblks_e)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                   'weighting factor for quadratic interpolation to surface')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'wgtfacq_e', p_metrics%wgtfacq_e,           &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape2d_ecubed, loutput=.FALSE. )


      ! weighting factor for quadratic interpolation to model top
      ! wgtfacq1_c    p_metrics%wgtfacq1_c(nproma,3,nblks_c)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                   'weighting factor for quadratic interpolation to model top')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'wgtfacq1_c', p_metrics%wgtfacq1_c,         &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape2d_ccubed, loutput=.FALSE. )


      ! weighting factor for quadratic interpolation to model top
      ! wgtfacq1_e   p_metrics%wgtfacq1_e(nproma,3,nblks_e)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                   'weighting factor for quadratic interpolation to model top')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'wgtfacq1_e', p_metrics%wgtfacq1_e,         &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape2d_ecubed, loutput=.FALSE. )


      ! Inverse layer thickness of full levels
      ! inv_ddqz_z_full   p_metrics%inv_ddqz_z_full(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Inverse_layer_thickness', 'm-1',                   &
      &                     'Inverse layer thickness of full levels')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'inv_ddqz_z_full', p_metrics%inv_ddqz_z_full, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_c )


      ! Inverse distance between full levels jk+1 and jk-1
      ! inv_ddqz_z_half2  p_metrics%inv_ddqz_z_half2(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Inverse_distance_between_full_levels', 'm-1',      &
      &                     'Inverse distance between full levels jk+1 and jk-1')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'inv_ddqz_z_half2', p_metrics%inv_ddqz_z_half2, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_c )



    IF (p_patch%cell_type== 3) THEN

      ! Vertical index of neighbor points needed for Taylor-expansion-based pressure gradient
      ! vertidx_gradp  p_metrics%vertidx_gradp(2,nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('Vertical_index', '-',                              &
      &                     'Vertical index')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'vertidx_gradp', p_metrics%vertidx_gradp,   &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_esquared, loutput=.FALSE. )


      ! Height differences between local edge point and neighbor cell points used for
      ! pressure gradient computation
      ! zdiff_gradp  p_metrics%zdiff_gradp(2,nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('Height_differences', 'm',                          &
      &                     'Height differences')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'zdiff_gradp', p_metrics%zdiff_gradp,       &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_esquared, loutput=.FALSE.  )


      ! Extrapolation factor for Exner pressure
      ! exner_exfac  p_metrics%exner_exfac(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Extrapolation_factor_for_Exner_pressure', '-',     &
      &                     'Extrapolation factor for Exner pressure')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'exner_exfac', p_metrics%exner_exfac,       &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c )


      ! Reference atmosphere field theta
      ! theta_ref_mc  p_metrics%theta_ref_mc(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_theta', 'K',            &
      &                     'Reference atmosphere field theta')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'theta_ref_mc', p_metrics%theta_ref_mc,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c )


      ! Reference atmosphere field theta
      ! theta_ref_ic  p_metrics%theta_ref_ic(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_theta', 'K',            &
      &                     'Reference atmosphere field theta')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'theta_ref_ic', p_metrics%theta_ref_ic,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf )


      ! Reference atmosphere field exner
      ! exner_ref_mc  p_metrics%exner_ref_mc(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_exner', '-',            &
      &                     'Reference atmosphere field exner')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'exner_ref_mc', p_metrics%exner_ref_mc,     &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,       &
                  & ldims=shape3d_c )


      ! Reference atmosphere field density
      ! rho_ref_mc  p_metrics%rho_ref_mc(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_density', '-',            &
      &                     'Reference atmosphere field density')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'rho_ref_mc', p_metrics%rho_ref_mc,         &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c )


      ! Reference atmosphere field exner
      ! d_exner_dz_ref_ic  p_metrics%d_exner_dz_ref_ic(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_exner', '-',            &
      &                     'Reference atmosphere field exner')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'd_exner_dz_ref_ic', p_metrics%d_exner_dz_ref_ic, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_chalf )


      ! Reference atmosphere field exner
      ! d2dexdz2_fac1_mc  p_metrics%d2dexdz2_fac1_mc(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_exner', '-',            &
      &                     'Reference atmosphere field exner')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'd2dexdz2_fac1_mc', p_metrics%d2dexdz2_fac1_mc, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,      &
                  & ldims=shape3d_c )


      ! Reference atmosphere field exner
      ! d2dexdz2_fac2_mc  p_metrics%d2dexdz2_fac2_mc(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_exner', '-',            &
      &                     'Reference atmosphere field exner')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'd2dexdz2_fac2_mc', p_metrics%d2dexdz2_fac2_mc, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,      &
                  & ldims=shape3d_c )


      ! mask field that excludes boundary halo points
      ! mask_prog_halo_c  p_metrics%mask_prog_halo_c(nproma,nblks_c)
      ! Note: Here "loutput" is set to .FALSE. since the output
      !       scheme operates on REAL model variables only and
      !       throws an error on this.
      !
      cf_desc    = t_cf_var('mask_field', '-',                                  &
      &                     'mask field that excludes boundary halo points')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'mask_prog_halo_c', p_metrics%mask_prog_halo_c, &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_SURFACE, cf_desc, grib2_desc,     &
                  & ldims=shape2d_c, loutput=.FALSE. )


    ELSE IF (p_patch%cell_type== 6) THEN


      ! geometric height at full level edges
      ! z_mc_e       p_metrics%z_mc_e(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('geometric_height_at_full_level_edge', 'm',           &
        &                   'geometric height at full level edge')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'z_mc_e', p_metrics%z_mc_e,                   &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                  & ldims=shape3d_e )


      ! slope of the coordinate lines in Northern direction
      ! ddnorth_z    p_metrics%ddnorth_z(nproma,nlev,nblks_v)
      !
      cf_desc    = t_cf_var('slope_of_coordinate_lines_in_Northern_direction', '-', &
      &                     'slope of coordinate lines in Northern direction')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_VERTEX)
      CALL add_var( p_metrics_list, 'ddnorth_z', p_metrics%ddnorth_z,               &
                  & GRID_UNSTRUCTURED_VERT, ZAXIS_HEIGHT, cf_desc, grib2_desc,      &
                  & ldims=shape3d_v )


      ! slope of the coordinate lines in Eastern direction
      ! ddeast_z     p_metrics%ddeast_z(nproma,nlev,nblks_v)
      !
      cf_desc    = t_cf_var('slope_of_coordinate_lines_in_Eastern_direction', '-', &
      &                     'slope of coordinate lines in Eastern direction')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_VERTEX)
      CALL add_var( p_metrics_list, 'ddeast_z', p_metrics%ddeast_z,                &
                  & GRID_UNSTRUCTURED_VERT, ZAXIS_HEIGHT, cf_desc, grib2_desc,     &
                  & ldims=shape3d_v )


      ! functional determinant of the metrics [sqrt(gamma)]
      ! ddqz_z_full_v   p_metrics%ddqz_z_full_v(nproma,nlev,nblks_v)
      !
      cf_desc    = t_cf_var('metrics_functional_determinant', '-',              &
      &                     'metrics functional determinant')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_VERTEX)
      CALL add_var( p_metrics_list, 'ddqz_z_full_v', p_metrics%ddqz_z_full_v,   &
                  & GRID_UNSTRUCTURED_VERT, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_v )

      ! functional determinant of the metrics [sqrt(gamma)]
      ! ddqz_z_full_r   p_metrics%ddqz_z_full_r(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('metrics_functional_determinant', '-',              &
      &                     'metrics functional determinant')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'ddqz_z_full_r', p_metrics%ddqz_z_full_r,   &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e )

      ! functional determinant of the metrics [sqrt(gamma)]
      ! ddqz_z_half_e   p_metrics%ddqz_z_half_e(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('metrics_functional_determinant', '-',              &
      &                     'metrics functional determinant')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_EDGE)
      CALL add_var( p_metrics_list, 'ddqz_z_half_e', p_metrics%ddqz_z_half_e,   &
                  & GRID_UNSTRUCTURED_EDGE, ZAXIS_HEIGHT, cf_desc, grib2_desc,  &
                  & ldims=shape3d_ehalf )

      ! 1/dz(k-1)-1/dz(k)
      ! diff_1_o_dz  p_metrics%diff_1_o_dz(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('difference_1_over_dz', 'm-1', 'difference 1 over dz')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'diff_1_o_dz', p_metrics%diff_1_o_dz,         &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                  & ldims=shape3d_chalf )

      ! 1/(dz(k-1)*dz(k))
      ! mult_1_o_dz  p_metrics%mult_1_o_dz(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('mult_1_over_dz', 'm-2', 'mult 1 over dz')
      grib2_desc = t_grib2_var( 255, 255, 255, ientr, GRID_REFERENCE, GRID_CELL)
      CALL add_var( p_metrics_list, 'mult_1_o_dz', p_metrics%mult_1_o_dz,         &
                  & GRID_UNSTRUCTURED_CELL, ZAXIS_HEIGHT, cf_desc, grib2_desc,    &
                  & ldims=shape3d_chalf )

    ENDIF

  END SUBROUTINE new_nh_metrics_list


  !------------------------------------------------------------------------------------------------
  !
  ! create (allocate) a new table entry
  ! reference to an existing pointer to a 3D tracer field
  ! optionally overwrite some default meta data 
  !
  SUBROUTINE add_var_list_reference_tracer(this_list, target_name, tracer_name, &
    &        tracer_idx, ptr_arr, cf, grib2, ldims, loutput, lrestart,          &
    &        tlev_source, tracer_info)

    TYPE(t_var_list)    , INTENT(inout)        :: this_list
    CHARACTER(len=*)    , INTENT(in)           :: target_name
    CHARACTER(len=*)    , INTENT(in)           :: tracer_name
    INTEGER             , INTENT(out)          :: tracer_idx          ! index in 4D tracer container
    TYPE(t_ptr_nh)      , INTENT(inout)        :: ptr_arr(:)
    TYPE(t_cf_var)      , INTENT(in)           :: cf                  ! CF related metadata
    TYPE(t_grib2_var)   , INTENT(in)           :: grib2               ! GRIB2 related metadata
    INTEGER             , INTENT(in), OPTIONAL :: ldims(3)            ! local dimensions, for checking
    LOGICAL             , INTENT(in), OPTIONAL :: loutput             ! output flag
    LOGICAL             , INTENT(in), OPTIONAL :: lrestart            ! restart flag
    INTEGER             , INTENT(in), OPTIONAL :: tlev_source         ! actual TL for TL dependent vars
    TYPE(t_tracer_meta) , INTENT(in), OPTIONAL :: tracer_info         ! tracer meta data

    ! Local variables:
    TYPE(t_list_element), POINTER :: target_element  
    TYPE(t_var_metadata), POINTER :: target_info

    ! get pointer to target element (in this case 4D tracer container)
    target_element => find_list_element (this_list, target_name)
    ! get tracer field metadata
    target_info => target_element%field%info

    ! get index in 4D
    tracer_idx = target_info%ncontained+1  ! index in 4D tracer container

    ! create new table entry reference
    CALL add_ref( this_list, target_name, tracer_name, ptr_arr(tracer_idx)%p_3d, &
       &          target_info%hgrid, target_info%vgrid, cf, grib2,               &
       &          ldims=ldims, loutput=loutput, lrestart=lrestart,               &
       &          tlev_source=tlev_source, tracer_info=tracer_info )


  END SUBROUTINE add_var_list_reference_tracer

END MODULE mo_nonhydro_state



