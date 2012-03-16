!>
!! Definition of optional (diagnostic) model variables
!!
!! In the metadata of each model variable ("add_var") it is specified
!! *how* certain post-processing tasks, e.g. vertical interpolation
!! onto p/z-levels, are treated. In the namelist, users can then
!! specify *if* computations for these variables are performed.
!!
!! If so, the resulting model fields are appended to the list of
!! internal post-processing tasks (each field forms its own task). As
!! we do not know in advance the contents of this list, we call them
!! "optional diagnostics".
!!
!! @author F. Prill (DWD)
!!
!! @par Revision History
!! Initial implementation by F. Prill, DWD (2012-03-07)
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
MODULE mo_opt_diagnostics

  USE mo_kind,                 ONLY: wp
  USE mo_parallel_config,      ONLY: nproma
  USE mo_linked_list,          ONLY: t_var_list
  USE mo_model_domain,         ONLY: t_patch
  USE mo_impl_constants,       ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_exception,            ONLY: finish
  USE mo_grid_config,          ONLY: n_dom
  USE mo_var_list,             ONLY: default_var_list_settings,     &
    &                                new_var_list, delete_var_list
  USE mo_var_list_element,     ONLY: level_type_pl,level_type_hl
  USE mo_cdi_constants,        ONLY: FILETYPE_NC2


  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = &
    & '$Id$'

  ! data types
  PUBLIC :: t_nh_opt_diag         ! optional diagnostic variables (data type)
  PUBLIC :: p_nh_opt_diag         ! state vector of optional diagnostic variables
                                  ! e.g. variables on p- and/or z-levels
  PUBLIC :: t_nh_diag_pz
  PUBLIC :: t_vcoeff
  ! subroutines
  PUBLIC :: vcoeff_allocate, vcoeff_deallocate
  PUBLIC :: construct_opt_diag
  PUBLIC :: destruct_opt_diag


  ! Derived type containing coefficient tables for vertical
  ! interpolation. There exist to different kinds of coefficients: For
  ! p- and for z-level-interpolation.
  TYPE t_vcoeff
    LOGICAL :: &
      & l_initialized = .FALSE.,  &
      & l_allocated   = .FALSE.

    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::  &
      &  wfac_lin, coef1, coef2, coef3
    INTEGER,  ALLOCATABLE, DIMENSION(:,:,:) ::  &
      &  idx0_lin, idx0_cub
    INTEGER,  ALLOCATABLE, DIMENSION(:,:)   ::  &
      &  bot_idx_lin, bot_idx_cub

    ! The following entries are only allocated for z-level
    ! interpolation:
    INTEGER,  ALLOCATABLE, DIMENSION(:,:)   ::  &
      &  kpbl1, kpbl2
    REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::  &
      &  wfacpbl1, wfacpbl2
  END TYPE t_vcoeff


  ! State vector for diagnostic variables on p- and/or z-levels
  ! 
  ! @note The pointers which are collected in this derived type
  !       constitute only the minimum set of fields that are required
  !       for pz-level interpolation. All other variables are stored
  !       inside the "opt_diag_list_p", "opt_diag_list_z" variable
  !       lists.
  TYPE t_nh_diag_pz

    REAL(wp), POINTER ::    &
      ! fields that are essential for pz-level interpolation:
      &  z_temp(:,:,:),        & ! temperature (nproma,nlev,nblks_c)          [K]
      &  z_pres(:,:,:),        & ! pressure (nproma,nlev,nblks_c)             [Pa]
      &  z_tracer_iqv(:,:,:),  & ! tracer concentration (i.e. prognostic cloud 
      &  z_tot_cld_iqv(:,:,:), & ! total cloud variables (cc,qv,qc,qi)        [kg/kg]
      ! fields that are essential for p-level interpolation only:
      &  p_geopot(:,:,:),      & ! geopotential (nproma,nlev,nblks_c)         [m2/s2]
      &  p_temp(:,:,:)           ! temperature (nproma,nlev,nblks_c)          [K]

    ! coefficient tables for vertical interpolation. There exist to
    ! different kinds of coefficients: For p- and for
    ! z-level-interpolation.
    TYPE(t_vcoeff) :: vcoeff_z, vcoeff_p

  END TYPE t_nh_diag_pz


  ! List of optional diagnostics + necessary meta data
  TYPE t_nh_opt_diag

    ! diag_pz: data structure containing coefficient tables and
    ! pointers to a few number of fields which are required for
    ! interpolation of model variables to p/z-levels
    TYPE(t_nh_diag_pz) :: diag_pz

    ! diag_pz_list: this contains all variables that have been
    ! interpolated onto p/z-levels
    TYPE(t_var_list)   :: opt_diag_list_p, opt_diag_list_z

  END TYPE t_nh_opt_diag


  ! Actual instantiation of optional diagnostics type "t_nh_opt_diag"
  TYPE(t_nh_opt_diag), TARGET, ALLOCATABLE :: p_nh_opt_diag(:)


CONTAINS


  !-------------
  !
  !> Add optional diagnostic variable lists (might remain empty)
  !
  SUBROUTINE construct_opt_diag(p_patch)
    TYPE(t_patch),                 INTENT(IN)   :: p_patch(n_dom)
    ! local variables
    CHARACTER(*), PARAMETER :: routine =  &
      &  TRIM("mo_opt_diagnostics:construct_opt_diag")
    INTEGER                            :: jg, ist
    CHARACTER(len=MAX_CHAR_LENGTH)     :: listname

    ! initialize data structure for optional diagnostics
    ALLOCATE(p_nh_opt_diag(n_dom), STAT=ist)
    IF (ist /= SUCCESS) &
      CALL finish (TRIM(routine), 'Allocation of optional diagnostics failed')

    DO jg = 1, n_dom
      WRITE(listname,'(a,i2.2)') 'nh_state_opt_diag_z_of_domain_',jg
      CALL new_var_list( p_nh_opt_diag(jg)%opt_diag_list_z, TRIM(listname), &
        & patch_id=p_patch(jg)%id, level_type=level_type_hl )
      CALL default_var_list_settings( p_nh_opt_diag(jg)%opt_diag_list_z,    &
        & lrestart=.FALSE., restart_type=FILETYPE_NC2  )

      WRITE(listname,'(a,i2.2)') 'nh_state_opt_diag_p_of_domain_',jg
      CALL new_var_list( p_nh_opt_diag(jg)%opt_diag_list_p, TRIM(listname), &
        & patch_id=p_patch(jg)%id, level_type=level_type_pl )
      CALL default_var_list_settings( p_nh_opt_diag(jg)%opt_diag_list_p,    &
        & lrestart=.FALSE., restart_type=FILETYPE_NC2  )
    ENDDO ! jg

  END SUBROUTINE construct_opt_diag


  !-------------
  !
  !> Clear optional diagnostic variable lists
  !
  SUBROUTINE destruct_opt_diag()
    ! local variables
    CHARACTER(*), PARAMETER :: routine =  &
      &  TRIM("mo_opt_diagnostics:destruct_opt_diag")
    INTEGER                            :: jg, ist

    DO jg = 1, n_dom
      CALL delete_var_list( p_nh_opt_diag(jg)%opt_diag_list_z )
      CALL delete_var_list( p_nh_opt_diag(jg)%opt_diag_list_p )
    ENDDO ! jg

    ! Delete optional diagnostics
    DEALLOCATE(p_nh_opt_diag, STAT=ist)
    IF (ist /= SUCCESS) &
      CALL finish(TRIM(routine),'Deallocation for optional diagnostics failed.')

  END SUBROUTINE destruct_opt_diag


  !-------------
  !>  
  ! Initialize a variable containing coefficient tables for vertical
  ! interpolation. There exist to different kinds of coefficients: For
  ! p- and for z-level-interpolation.
  SUBROUTINE vcoeff_allocate(p_patch, nlev, vcoeff)
    TYPE(t_patch),       TARGET,       INTENT(IN)    :: p_patch    
    INTEGER,                           INTENT(IN)    :: nlev
    TYPE(t_vcoeff),                    INTENT(INOUT) :: vcoeff

    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_nh_vert_interp:vcoeff_allocate")
    INTEGER :: ierrstat

    IF (.NOT. vcoeff%l_allocated) THEN
      ! real(wp)
      ALLOCATE( &
        &  vcoeff%wfac_lin(nproma,nlev,p_patch%nblks_c), &
        &  vcoeff%coef1(nproma,nlev,p_patch%nblks_c),    &
        &  vcoeff%coef2(nproma,nlev,p_patch%nblks_c),    &
        &  vcoeff%coef3(nproma,nlev,p_patch%nblks_c),    &
        &  STAT=ierrstat )
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

      ! integer
      ALLOCATE( &
        &  vcoeff%idx0_lin(nproma,nlev,p_patch%nblks_c), &
        &  vcoeff%idx0_cub(nproma,nlev,p_patch%nblks_c), &
        &  vcoeff%bot_idx_lin(nproma,p_patch%nblks_c),    &
        &  vcoeff%bot_idx_cub(nproma,p_patch%nblks_c),    &
        &  STAT=ierrstat )
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

      ! The following entries are only allocated for z-level
      ! interpolation:
      ALLOCATE( &
        &  vcoeff%wfacpbl1(nproma,p_patch%nblks_c),            &
        &  vcoeff%wfacpbl2(nproma,p_patch%nblks_c),            &
        &  STAT=ierrstat )
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      ALLOCATE( &
        &  vcoeff%kpbl1(nproma,p_patch%nblks_c),               &
        &  vcoeff%kpbl2(nproma,p_patch%nblks_c),               &
        &  STAT=ierrstat )
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

      vcoeff%l_allocated = .TRUE.
    END IF
  END SUBROUTINE vcoeff_allocate


  !-------------
  !>  
  ! Clear a variable containing coefficient tables for vertical
  ! interpolation. There exist to different kinds of coefficients: For
  ! p- and for z-level-interpolation.
  SUBROUTINE vcoeff_deallocate(vcoeff)
    TYPE(t_vcoeff),                    INTENT(INOUT) :: vcoeff

    CHARACTER(*), PARAMETER :: routine = &
      &  TRIM("mo_nh_vert_interp:vcoeff_deallocate")
    INTEGER :: ierrstat

    ! deallocate coefficient tables:
    IF (vcoeff%l_allocated) THEN
      ! real(wp)
      DEALLOCATE( &
        &  vcoeff%wfac_lin, &
        &  vcoeff%coef1,    &
        &  vcoeff%coef2,    &
        &  vcoeff%coef3,    &
        &  STAT=ierrstat )
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
      ! integer
      DEALLOCATE( &
        &  vcoeff%idx0_lin,       &
        &  vcoeff%idx0_cub,       &
        &  vcoeff%bot_idx_lin,    &
        &  vcoeff%bot_idx_cub,    &
        &  STAT=ierrstat )
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

      DEALLOCATE( vcoeff%wfacpbl1, vcoeff%wfacpbl2,      &
        &         STAT=ierrstat )
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
      DEALLOCATE( vcoeff%kpbl1, vcoeff%kpbl2,               &
        &         STAT=ierrstat )
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

      vcoeff%l_allocated = .FALSE.
    END IF

    vcoeff%l_initialized = .FALSE.
  END SUBROUTINE vcoeff_deallocate


END MODULE mo_opt_diagnostics



