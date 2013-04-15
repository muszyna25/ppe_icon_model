!>
!! Contains the setup of variables related to large eddy simulation setup
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
MODULE mo_les_nml

  USE mo_les_config,          ONLY: les_config
  USE mo_dynamics_config,     ONLY: iequations
  USE mo_kind,                ONLY: wp
  USE mo_mpi,                 ONLY: my_process_is_stdio 
  USE mo_exception,           ONLY: message, finish
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
                                  & open_and_restore_namelist, close_tmpfile
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_dom

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_les_namelist

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  REAL(wp) :: sst        ! prescribed SST
  REAL(wp) :: shflx      ! prescribed sensible heat flux (Km/s)
  REAL(wp) :: lhflx      ! prescribed latent heat flux   (Km/s)
  INTEGER  :: isrfc_type ! 1=fixed sst, 2=fixed flux

  REAL(wp) :: ugeo(2)    ! ugeo(1)=constant, ugeo(2)=gradient
  REAL(wp) :: vgeo(2)    ! vgeo(1)=constant, vgeo(2)=gradient
  REAL(wp) :: umean(2)   ! umean(1)=constant, umean(2)=gradient
  REAL(wp) :: vmean(2)   ! vmean(1)=constant, vmean(2)=gradient
  REAL(wp) :: ufric      ! friction velocity
 
  LOGICAL  :: is_dry_cbl  !special case for CBL testcase
  LOGICAL  :: set_geowind !TRUE is geostrophic wind is set
 
  !Some parameters
  REAL(wp) :: karman_constant
  REAL(wp) :: rkarman_constant  !inverse karman constant
  REAL(wp) :: smag_constant
  REAL(wp) :: turb_prandtl 
  REAL(wp) :: rturb_prandtl     !inverse turbulent prandtl number
 
  NAMELIST/les_nml/ sst, shflx, lhflx, isrfc_type, ugeo, vgeo, umean,       &
                    vmean, ufric, is_dry_cbl, set_geowind, karman_constant, &
                    rkarman_constant, smag_constant, turb_prandtl,          &
                    rturb_prandtl

CONTAINS
  !-------------------------------------------------------------------------
  !>
  !! Read Namelist for LES
  !!
  !! This subroutine 
  !! - reads the Namelist 
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state 
  !!
  !! @par Revision History
  !!  by Anurag Dipankar, MPIM (2013-04)
  !!
  SUBROUTINE read_les_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename 
    INTEGER :: istat, funit, jg

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_les_nml: read_les_namelist'

    !-----------------------
    ! 1. default settings
    !-----------------------
    sst          = 300._wp
    ugeo(1:2)    = 0._wp 
    vgeo(1:2)    = 0._wp 
    umean(1:2)   = 0._wp
    vmean(1:2)   = 0._wp
    shflx        = 0._wp 
    lhflx        = 0._wp 
    isrfc_type   = 2      !fixed SST
    ufric        = -1._wp 

    is_dry_cbl   = .FALSE.
    set_geowind  = .FALSE.

    !parameters
    karman_constant  = 0.4_wp
    rkarman_constant = 2.5_wp
    smag_constant    = 0.23_wp
    turb_prandtl     = 0.33333333333_wp
    rturb_prandtl    = 3.0_wp

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('les_nml')
      READ(funit,NML=les_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('les_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, les_nml)
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------
    IF(isrfc_type==1)THEN
       shflx = 0._wp   
       lhflx = 0._wp   
    END IF
    IF(is_dry_cbl)THEN
       lhflx = 0._wp
    END IF

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------
    DO jg = 1 , max_dom
      les_config(jg)% sst          =  sst
      les_config(jg)% ugeo(:)    =  ugeo(:)
      les_config(jg)% vgeo(:)    =  vgeo(:)
      les_config(jg)% umean(:)   =  umean(:)
      les_config(jg)% vmean(:)   =  vmean(:)
      les_config(jg)% shflx        =  shflx
      les_config(jg)% lhflx        =  lhflx
      les_config(jg)% isrfc_type   =  isrfc_type
      les_config(jg)% ufric        =  ufric
      les_config(jg)% is_dry_cbl   =  is_dry_cbl
      les_config(jg)% set_geowind  =  set_geowind
      les_config(jg)% karman_constant   =  karman_constant
      les_config(jg)% rkarman_constant  =  rkarman_constant
      les_config(jg)% smag_constant     =  smag_constant
      les_config(jg)% turb_prandtl      =  turb_prandtl
      les_config(jg)% rturb_prandtl     =  rturb_prandtl
    END DO

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=les_nml)                    
      CALL store_and_close_namelist(funit,'les_nml') 
    ENDIF
    ! 7. write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=les_nml)

  END SUBROUTINE read_les_namelist

END MODULE mo_les_nml
