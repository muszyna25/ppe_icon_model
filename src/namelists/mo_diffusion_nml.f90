!>
!! Contains the setup of variables related to horizontal diffusion
!!        
!! @par Revision History
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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
MODULE mo_diffusion_nml

  USE mo_diffusion_config,    ONLY: diffusion_config
  USE mo_kind,                ONLY: wp
  USE mo_mpi,                 ONLY: my_process_is_stdio 
  USE mo_exception,           ONLY: message, finish
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_diffusion_namelist

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !-------------------------------------------------------------------------
  ! 1.0 Namelist variables and auxiliary parameters setting up the
  !     configuration of the dynamical core
  !-------------------------------------------------------------------------
  INTEGER :: hdiff_order  ! order of horizontal diffusion
                          ! -1: no diffusion
                          ! 2: 2nd order linear diffusion on all vertical levels 
                          ! 3: Smagorinsky diffusion for hexagonal model
                          ! 4: 4th order linear diffusion on all vertical levels 
                          ! 5: Smagorinsky diffusion for triangular model
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
  REAL(wp) :: hdiff_min_efdt_ratio  ! minimum value of hdiff_efdt_ratio 
                                    ! (for upper sponge layer)
  REAL(wp) :: hdiff_tv_ratio        ! the ratio of diffusion coefficient: temp:mom
  REAL(wp) :: hdiff_smag_fac        ! scaling factor for Smagorinsky diffusion
  REAL(wp) :: hdiff_multfac         ! multiplication factor of normalized diffusion
                                    ! coefficient for nested domains
  LOGICAL :: lhdiff_temp   ! if .TRUE., apply horizontal diffusion to temp.
  LOGICAL :: lhdiff_vn     ! if .TRUE., apply horizontal diffusion to momentum.

  NAMELIST/diffusion_nml/ hdiff_order, k2_klev_max, k2_pres_max,         &
                          hdiff_efdt_ratio, hdiff_min_efdt_ratio,        &
                          hdiff_tv_ratio, hdiff_smag_fac, hdiff_multfac, &
                          lhdiff_temp, lhdiff_vn

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
    INTEGER :: istat, funit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_diffusion_nml: read_diffusion_namelist'

    !-----------------------
    ! 1. default settings
    !-----------------------
    lhdiff_temp          = .TRUE.
    lhdiff_vn            = .TRUE.

    hdiff_order          = 4
    hdiff_efdt_ratio     = 1.0_wp
    hdiff_min_efdt_ratio = 1.0_wp
    hdiff_multfac        = 1.0_wp
    hdiff_smag_fac       = 0.15_wp
    hdiff_tv_ratio       = 1.0_wp

    k2_pres_max          = -99.0_wp                                                    
    k2_klev_max          = 0

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('diffusion_nml')
      READ(funit,NML=diffusion_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('diffusion_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, diffusion_nml)
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
    diffusion_config(:)% hdiff_order          =  hdiff_order
    diffusion_config(:)% k2_klev_max          =  k2_klev_max
    diffusion_config(:)% k2_pres_max          =  k2_pres_max
    diffusion_config(:)% hdiff_efdt_ratio     =  hdiff_efdt_ratio
    diffusion_config(:)% hdiff_min_efdt_ratio =  hdiff_min_efdt_ratio
    diffusion_config(:)% hdiff_smag_fac       =  hdiff_smag_fac
    diffusion_config(:)% hdiff_multfac        =  hdiff_multfac
    diffusion_config(:)% hdiff_tv_ratio       =  hdiff_tv_ratio 

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
