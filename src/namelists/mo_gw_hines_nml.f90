!>
!! Namelist for the configuration of the vertical diffusion
!!
!!
!! @par Revision History
!! Revision history in mo_echam_vdiff_params.f90 (r4300)
!! Modification by Constantin Junk, MPI-M (2011-05-05)
!! - moved echam_vdiff namelist variables and subroutine setup_vdiff
!!   from mo_echam_vdiff_params to namelists/mo_echam_vdiff_nml
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
MODULE mo_gw_hines_nml

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_dom
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_exception,           ONLY: message, print_value
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_master_nml,          ONLY: lrestart
  USE mo_gw_hines_config,     ONLY: gw_hines_config
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist, &     
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_mpi,                ONLY: p_pe, p_io

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_gw_hines_namelist
  PUBLIC :: gw_hines_nml            !< namelist for Hines gravity wave parameterization

  PUBLIC :: lheatcal, emiss_lev, rmscon, kstar, m_min
!!$  PUBLIC :: lfront, rms_front, front_thres
!!$  PUBLIC :: lozpr, pcrit, pcons
!!$  PUBLIC :: lrmscon_lat, lat_rmscon_lo, lat_rmscon_hi, rmscon_lo, rmscon_hi

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !-----------------------------------!
  ! gw_hines_nml namelist variables   !
  !-----------------------------------!

  LOGICAL  :: lheatcal      !< true : compute momentum flux dep., heating and diffusion coefficient
                            !< false: compute only momentum flux deposition

  INTEGER  :: emiss_lev     !< root mean square gravity wave wind at lowest level (m/s)
  REAL(wp) :: rmscon        !< number of levels above the ground at which gw are emitted
  REAL(wp) :: kstar         !< typical gravity wave horizontal wavenumber (1/m)
  REAL(wp) :: m_min         !< minimum bound in  vertical wavenumber (1/m)

!!$  LOGICAL  :: lfront        !< true: compute gw sources emerging from fronts and background
!!$                            !< (Charron and Manzini, 2002)
!!$  REAL(wp) :: rms_front     !< rms frontal gw wind at source level  (m/s)
!!$  REAL(wp) :: front_thres   !< minimum value of the frontogenesis function, for which
!!$                            !< gravity waves are emitted from fronts [(K/m)^2/hour]
!!$
!!$  LOGICAL  :: lozpr         !< true: for background enhancement associated with precipitation
!!$                            !< (Manzini et al., 1997)
!!$  REAL(wp) :: pcrit         !< critical precipitation value (mm/d), above which 
!!$                            !< gravity wave rms wind enhancement is applied
!!$  REAL(wp) :: pcons         !< adimensional factor for background enhancement 
!!$                            !< associated with precipitation
!!$
!!$  LOGICAL  :: lrmscon_lat   !< true:   use latitude dependent rmscon as defined
!!$                            !< through rmscon_lo, rmscon_hi, lat_rmscon_lo, and lat_rmscon_hi
!!$                            !< false:  use uniform rmscon
!!$                            !< attention: may be overwritten if lfront or lozpr is true
!!$  REAL(wp) :: lat_rmscon_lo !< rmscon_lo is used equatorward of this latitude (degN)
!!$  REAL(wp) :: lat_rmscon_hi !< rmscon_hi is used poleward of this latitude (degN)
!!$  REAL(wp) :: rmscon_lo     !< rmscon used equatorward of lat_rmscon_lo
!!$  REAL(wp) :: rmscon_hi     !< rmscon used poleward of lat_rmscon_hi


NAMELIST /gw_hines_nml/ &
  & lheatcal, rmscon, emiss_lev, kstar, m_min !!$,                    &
!!$  & lfront, rms_front, front_thres,                                &
!!$  & lozpr, pcrit, pcons,                                           &
!!$  & lrmscon_lat, lat_rmscon_lo, lat_rmscon_hi, rmscon_lo, rmscon_hi

CONTAINS

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for Hines gw parameterization. 
  !!
  !! This subroutine 
  !! - reads the Namelist for Hines gw parameterization
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)    
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-06-07)
  !!
  SUBROUTINE read_gw_hines_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg 

    !0!CHARACTER(len=*), PARAMETER :: routine = 'mo_gw_hines_nml:read_gw_hines_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------
    lheatcal = .TRUE.

    emiss_lev = 10          ! is correct for L31 and L47
    rmscon    = 1.0_wp      ! default value used in ECHAM5
    kstar     = 5.0e-5_wp   ! = 2*pi/(126000 m)
    m_min     = 0.0_wp

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('gw_hines_nml')
      READ(funit,NML=gw_hines_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('gw_hines_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, gw_hines_nml)
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------

    DO jg = 1,max_dom
      gw_hines_config(jg)%lheatcal  = lheatcal
      gw_hines_config(jg)%emiss_lev = emiss_lev
      gw_hines_config(jg)%rmscon    = rmscon
      gw_hines_config(jg)%kstar     = kstar
      gw_hines_config(jg)%m_min     = m_min
    ENDDO

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=gw_hines_nml)                    
    CALL store_and_close_namelist(funit, 'gw_hines_nml') 

    ! 6. write the contents of the namelist to an ASCII file
    !
    IF(p_pe == p_io) WRITE(nnml_output,nml=gw_hines_nml)

  END SUBROUTINE read_gw_hines_namelist


END MODULE mo_gw_hines_nml
