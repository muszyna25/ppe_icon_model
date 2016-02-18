!>
!! Contains the setup of variables related to horizontal diffusion
!!        
!! @par Revision History
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
MODULE mo_diffusion_nml

  USE mo_diffusion_config,    ONLY: diffusion_config
  USE mo_dynamics_config,     ONLY: iequations
  USE mo_kind,                ONLY: wp
  USE mo_mpi,                 ONLY: my_process_is_stdio 
  USE mo_exception,           ONLY: message, finish
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_diffusion_namelist

  !-------------------------------------------------------------------------
  ! 1.0 Namelist variables and auxiliary parameters setting up the
  !     configuration of the dynamical core
  !-------------------------------------------------------------------------
  INTEGER :: hdiff_order  ! order of horizontal diffusion
                          ! -1: no diffusion
                          ! 2: 2nd order linear diffusion on all vertical levels 
                          ! 3: Smagorinsky diffusion without background diffusion
                          ! 4: 4th order linear diffusion on all vertical levels 
                          ! 5: Smagorinsky diffusion with fourth-order background diffusion
                          ! 24 or 42: 2nd order linear diffusion for upper levels,
                          !           4th order for lower levels

  REAL(wp) :: k2_pres_max ! (relevant only when hdiff_order = 24 or 42)
                          ! pressure (in Pa) specified by the user
                          ! to determine the lowest vertical level 
                          ! to which 2nd order linear diffusion is applied.
                          ! For the levels with pressure > k2_pres_max, 
                          ! 4th order linear diffusion is applied. 

  INTEGER  :: k2_klev_max ! (relevant only when hdiff_order = 24 or 42)
                          ! vertical level index specified by the user
                          ! to determine the lowest vertical level 
                          ! to which 2nd order linear diffusion is applied.
                          ! For the levels with k > k2_klev_max, 
                          ! 4th order linear diffusion is applied. 

  REAL(wp) :: hdiff_efdt_ratio      ! ratio of e-folding time to (2*)time step
  REAL(wp) :: hdiff_w_efdt_ratio    ! ratio of e-folding time to time step for w diffusion (NH only)
  REAL(wp) :: hdiff_min_efdt_ratio  ! minimum value of hdiff_efdt_ratio 
                                    ! (for upper sponge layer)
  REAL(wp) :: hdiff_tv_ratio        ! the ratio of diffusion coefficient: temp:mom
  REAL(wp) :: hdiff_smag_fac        ! scaling factor for Smagorinsky diffusion
  REAL(wp) :: hdiff_multfac         ! multiplication factor of normalized diffusion
                                    ! coefficient for nested domains
  INTEGER  :: itype_vn_diffu        ! options for discretizing the Smagorinsky momentum diffusion
  INTEGER  :: itype_t_diffu         ! options for discretizing the Smagorinsky temperature diffusion
  LOGICAL :: lhdiff_temp   ! if .TRUE., apply horizontal diffusion to temp.
  LOGICAL :: lhdiff_vn     ! if .TRUE., apply horizontal diffusion to horizontal momentum.
  LOGICAL :: lhdiff_w      ! if .TRUE., apply horizontal diffusion to vertical momentum.
  LOGICAL :: lsmag_3d      ! if .TRUE., compute 3D Smagorinsky diffusion coefficient.

  NAMELIST/diffusion_nml/ hdiff_order, k2_klev_max, k2_pres_max,              &
                          hdiff_efdt_ratio, hdiff_min_efdt_ratio,             &
                          hdiff_tv_ratio, hdiff_smag_fac, hdiff_multfac,      &
                          lhdiff_temp, lhdiff_vn, itype_vn_diffu,             &
                          itype_t_diffu, hdiff_w_efdt_ratio, lhdiff_w, lsmag_3d

CONTAINS
  !-------------------------------------------------------------------------
  !>
  !! Read Namelist for horizontal diffusion. 
  !!
  !! This subroutine 
  !! - reads the Namelist for diffusion
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)  
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-07-06)
  !!
  SUBROUTINE read_diffusion_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename 
    INTEGER :: istat, funit, iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_diffusion_nml: read_diffusion_namelist'

    !-----------------------
    ! 1. default settings
    !-----------------------
    lhdiff_temp          = .TRUE.
    lhdiff_vn            = .TRUE.
    lhdiff_w             = .TRUE.
    lsmag_3d             = .FALSE.

    IF (iequations == 3) THEN
      hdiff_order          = 5
      hdiff_efdt_ratio     = 36.0_wp
      hdiff_smag_fac       = 0.015_wp
    ELSE
      hdiff_order          = 4
      hdiff_efdt_ratio     = 1.0_wp
      hdiff_smag_fac       = 0.15_wp
    ENDIF

    hdiff_min_efdt_ratio = 1.0_wp
    hdiff_w_efdt_ratio   = 15.0_wp
    hdiff_multfac        = 1.0_wp
    hdiff_tv_ratio       = 1.0_wp
    itype_vn_diffu       = 1
    itype_t_diffu        = 2

    k2_pres_max          = -99.0_wp                                                    
    k2_klev_max          = 0

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('diffusion_nml')
      READ(funit,NML=diffusion_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('diffusion_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, diffusion_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, diffusion_nml)                                       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, diffusion_nml)   ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------
    SELECT CASE( hdiff_order)
    CASE(-1)
      CALL message(TRIM(routine),'Horizontal diffusion switched off.')
      lhdiff_temp = .FALSE.
      lhdiff_vn   = .FALSE.
      lhdiff_w    = .FALSE.

    CASE(2,3,4,5,24,42)

      IF ((.NOT.lhdiff_temp).AND.(.NOT.lhdiff_vn)) THEN
        CALL message('','')
        CALL message('','lhdiff_temp and lhdiff_vn both set to .FALSE. by user.')
        CALL message('','Horizontal diffusion is thus switched off and '// &
                        'hdiff_order reset to -1')
        CALL message('','')
        hdiff_order = -1
      END IF

    CASE DEFAULT
      CALL finish(TRIM(routine),                         &
        & 'Error: Invalid choice of  hdiff_order. '// &                
        & 'Choose from -1, 2, 3, 4, 5, 24, and 42.')
    END SELECT

    IF ( hdiff_efdt_ratio<=0._wp) THEN
      CALL message(TRIM(routine),'No horizontal background diffusion is used')
    ENDIF

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------
    diffusion_config(:)% lhdiff_temp          =  lhdiff_temp
    diffusion_config(:)% lhdiff_vn            =  lhdiff_vn
    diffusion_config(:)% lhdiff_w             =  lhdiff_w
    diffusion_config(:)% lsmag_3d             =  lsmag_3d
    diffusion_config(:)% hdiff_order          =  hdiff_order
    diffusion_config(:)% k2_klev_max          =  k2_klev_max
    diffusion_config(:)% k2_pres_max          =  k2_pres_max
    diffusion_config(:)% hdiff_efdt_ratio     =  hdiff_efdt_ratio
    diffusion_config(:)% hdiff_w_efdt_ratio   =  hdiff_w_efdt_ratio
    diffusion_config(:)% hdiff_min_efdt_ratio =  hdiff_min_efdt_ratio
    diffusion_config(:)% hdiff_smag_fac       =  hdiff_smag_fac
    diffusion_config(:)% hdiff_multfac        =  hdiff_multfac
    diffusion_config(:)% hdiff_tv_ratio       =  hdiff_tv_ratio 
    diffusion_config(:)%itype_vn_diffu        =  itype_vn_diffu
    diffusion_config(:)%itype_t_diffu         =  itype_t_diffu 

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=diffusion_nml)                    
      CALL store_and_close_namelist(funit,'diffusion_nml') 
    ENDIF
    ! 7. write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=diffusion_nml)

  END SUBROUTINE read_diffusion_namelist

END MODULE mo_diffusion_nml
