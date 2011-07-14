!>
!! Namelist specific to the hydro atm dynamical core. 
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

  USE mo_ha_dyn_config,         ONLY: ha_dyn_config
  USE mo_kind,                  ONLY: wp
  USE mo_mpi,                   ONLY: p_pe, p_io
  USE mo_io_units,              ONLY: nnml, nnml_output
  USE mo_exception,             ONLY: message, message_text, finish
  USE mo_namelist,              ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_master_nml,            ONLY: lrestart
  USE mo_io_restart_attributes, ONLY: get_restart_attribute
  USE mo_io_restart_namelist,   ONLY: open_tmpfile, store_and_close_namelist, &       
                                      open_and_restore_namelist, close_tmpfile
  USE mo_impl_constants,        ONLY:  AB2

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_ha_dyn_namelist

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !---------------------
  ! namelist variables
  !---------------------

  INTEGER  :: nml_ileapfrog_startup  ! choice of the first time step in
                                     ! a leapfrog time stepping scheme
                                     ! 1 = Euler forward
                                     ! 2 = several sub-steps

  REAL(wp) :: nml_asselin_coeff      ! parameter used in Asselin filter

  INTEGER  :: nml_si_expl_scheme     ! scheme for the explicit part of the
                                     ! 2-time-level semi-implicit time integration.
                                     ! See mo_impl_constants for the options.
  REAL(wp) :: nml_si_2tls

  REAL(wp) :: nml_si_rtol     ! relative tolerance

  REAL(wp) :: nml_si_coeff    !  = 0 : explicit scheme(for *d*,*t*,*alps*).
                              !  = 1 : semi implicit scheme.
                              !  in (0,1): a weighted scheme

  REAL(wp) :: nml_si_offctr   ! weighting parameter used in calculating the
                              ! second temporal derivatives in the semi-implicit
                              ! correction scheme. The value read from namelist are
                              ! assumed to be the offcentering (i.e. between 0 and 1).



  REAL(wp) :: nml_si_cmin     ! min. phase speed of the decomposed modes to be
                              ! solved by the semi-implicit correction scheme
  LOGICAL  :: nml_lsi_3d      ! if .true., solve the 3D equation

  LOGICAL  :: nml_lref_temp   ! if .TRUE., involve the reference temperature profile
                              ! in the calculation of pressure gradient force.

  NAMELIST/ha_dyn_nml/ nml_ileapfrog_startup, nml_asselin_coeff,     &
                       nml_si_expl_scheme, nml_si_2tls, nml_si_rtol, &
                       nml_si_coeff, nml_si_offctr, nml_si_cmin,     &
                       nml_lsi_3d, nml_lref_temp

CONTAINS
  !>
  !!
  SUBROUTINE read_ha_dyn_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    CHARACTER(len=*),PARAMETER :: routine = 'mo_ha_dyn_nml:read_ha_dyn_namelist'

    !------------------------------------------------------------
    ! Default values
    !------------------------------------------------------------
    nml_ileapfrog_startup = 1
    nml_asselin_coeff     = 0.1_wp
 
    nml_si_expl_scheme    = AB2 
    nml_si_2tls           = 0.6_wp 
    nml_si_rtol           = 1.e-3_wp

    nml_si_coeff          = 1.0_wp
    nml_si_offctr         = 0.7_wp
    nml_si_cmin           = 30._wp
    nml_lsi_3d            = .FALSE.
 
    nml_lref_temp         = .FALSE.

    !------------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above by 
    ! values in the restart file
    !------------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('ha_dyn_nml')
      READ(funit,NML=ha_dyn_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processors)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('ha_dyn_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, ha_dyn_nml)
    END SELECT
    CALL close_nml

    !-----------------------------------------------------
    ! Sanity Check
    !-----------------------------------------------------
    IF (nml_asselin_coeff<0._wp) CALL finish( TRIM(routine), &
      'wrong (negative) coefficient of Asselin filter')

    IF (nml_si_offctr>1._wp.OR.nml_si_offctr<0._wp) CALL finish( TRIM(routine), &
      'Invalid value for parameter nml_si_offctr. Valid range is [0,1].')

    IF (nml_si_2tls<0.5_wp.OR.nml_si_2tls>1._wp) CALL finish( TRIM(routine), &
      'Improper value for parameter nml_si_2tls. Should be in [0.5,1]')

    IF (nml_lref_temp) THEN
      CALL message(TRIM(routine),                            &
          'use of reference temperature switched ON in ' //  &
          'calculation of pressure gradient force.')
    ELSE                                                                    
      CALL message(TRIM(routine),                            &
           'use of reference temperature switched OFF in '// &
           'calculation of pressure gradient force.')
    ENDIF 

    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=ha_dyn_nml)
    CALL store_and_close_namelist(funit, 'ha_dyn_nml')

    ! write the contents of the namelist to an ASCII file
    IF(p_pe == p_io) WRITE(nnml_output,nml=ha_dyn_nml)

    !-----------------------------------------------------
    ! Fill configuration state
    !-----------------------------------------------------
    ha_dyn_config% ileapfrog_startup = nml_ileapfrog_startup
    ha_dyn_config% asselin_coeff     = nml_asselin_coeff
    ha_dyn_config% si_expl_scheme    = nml_si_expl_scheme
    ha_dyn_config% si_2tls           = nml_si_2tls
    ha_dyn_config% si_rtol           = nml_si_rtol
    ha_dyn_config% si_coeff          = nml_si_coeff
    ha_dyn_config% si_offctr         = nml_si_offctr
    ha_dyn_config% si_cmin           = nml_si_cmin
    ha_dyn_config% lsi_3d            = nml_lsi_3d
    ha_dyn_config% lref_temp         = nml_lref_temp

  END SUBROUTINE read_ha_dyn_namelist

END MODULE mo_ha_dyn_nml
