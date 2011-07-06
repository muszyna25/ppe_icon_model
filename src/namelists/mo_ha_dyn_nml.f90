!>
!! Namelist module for the hydro atm dynamical core.
!! 
!! @par Revision History
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_ha_dyn_nml

  USE mo_kind,                  ONLY: wp
  USE mo_io_units,              ONLY: nnml, nnml_output
  USE mo_namelist,              ONLY: position_nml, positioned
  USE mo_master_nml,            ONLY: lrestart
  USE mo_ha_dyn_config,         ONLY: ha_dyn_config
  USE mo_io_restart_attributes, ONLY: get_restart_attribute
  USE mo_io_restart_namelist,   ONLY: open_tmpfile, store_and_close_namelist, &       
                                      open_and_restore_namelist, close_tmpfile

  IMPLICIT NONE

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  PUBLIC

  !---------------------------------------------------------------
  ! Namelist variables and auxiliary parameters setting up the
  ! configuration of the dynamical core
  !---------------------------------------------------------------

  INTEGER  :: ileapfrog_startup  ! choice of the first time step in
                                 ! a leapfrog time stepping scheme
                                 ! 1 = Euler forward
                                 ! 2 = several sub-steps

  REAL(wp) :: asselin_coeff      ! parameter used in Asselin filter

  INTEGER  :: si_expl_scheme     ! scheme for the explicit part of the
                                 ! 2-time-level semi-implicit time integration.
                                 ! See mo_impl_constants for the options.

  REAL(wp) :: si_coeff           !  = 0 : explicit scheme(for *d*,*t*,*alps*).
                                 !  = 1 : semi implicit scheme.
                                 !  in (0,1): a weighted scheme

  REAL(wp) :: si_offctr          ! weighting parameter used in calculating the
                                 ! second temporal derivatives in the semi-implicit
                                 ! correction scheme. The value read from namelist are
                                 ! assumed to be the offcentering (i.e. between 0 and 1).

  LOGICAL  :: lsi_3d             ! if .true., solve the 3D equation

  REAL(wp) :: si_rtol            ! relative tolerance

  REAL(wp) :: si_cmin            ! min. phase speed of the decomposed modes to be
                                 ! solved by the semi-implicit correction scheme

  REAL(wp) :: sw_ref_height      ! reference height to linearize around if using
                                 ! lshallow_water and semi-implicit correction
  
  LOGICAL :: ldry_dycore ! if .TRUE., ignore the effact of water vapor,
                         ! cloud liquid and cloud ice on virtual temperature.
  LOGICAL :: lref_temp   ! if .TRUE., involve the reference temperature profile
                         ! in the calculation of pressure gradient force.

  NAMELIST/ha_dyn_nml/ ileapfrog_startup, asselin_coeff, &
                       si_coeff, si_offctr, si_rtol,     &
                       si_cmin, lsi_3d, si_expl_scheme,  &
                       ldry_dycore, lref_temp

CONTAINS
  !>
  !!
  SUBROUTINE read_ha_dyn_namelist()

    INTEGER :: istat, funit

    !------------------------------------------------------------
    ! 1. Set up the default values for dynamics_ctl
    !------------------------------------------------------------
 
    ileapfrog_startup = 1
    asselin_coeff     = 0.1_wp
 
    si_expl_scheme    = AB2 
 
    si_coeff          = 1.0_wp
    si_offctr         = 0.7_wp
    lsi_3d            = .FALSE.
    si_rtol           = 1.e-3_wp
    si_cmin           = 30._wp
 
    ldry_dycore       = .TRUE.
    lref_temp         = .FALSE.

    !------------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above by 
    !    values in the restart file
    !------------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('dynamics_ctl')
      READ(funit,NML=dynamics_ctl)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! 3. Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL position_nml ('dynamics_ctl', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, dynamics_ctl)
    END SELECT

    !-----------------------------------------------------
    ! 4. Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=dynamics_ctl)
    CALL store_and_close_namelist(funit, 'dynamics_ctl')

  ! ! write the contents of the namelist to an ASCII file
  ! IF(p_pe == p_io) WRITE(nnml_output,nml=dynamics_ctl)

    !-----------------------------------------------------
    ! 5. Fill configuration state
    !-----------------------------------------------------

    ha_dyn_config(:)% ileapfrog_startup = ileapfrog_startup
    ha_dyn_config(:)% asselin_coeff     = asselin_coeff
    ha_dyn_config(:)% si_expl_scheme    = si_expl_scheme
    ha_dyn_config(:)% si_coeff          = si_coeff
    ha_dyn_config(:)% si_offctr         = si_offctr
    ha_dyn_config(:)% si_cmin           = si_cmin
    ha_dyn_config(:)% si_rtol           = si_rtol
    ha_dyn_config(:)% lsi_3d            = lsi_3d
    ha_dyn_config(:)% ldry_dycore       = ldry_dycore
    ha_dyn_config(:)% lref_temp         = lref_temp

  END SUBROUTINE read_ha_dyn_namelist

END MODULE mo_ha_dyn_nml
