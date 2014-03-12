!>
!! Namelist for turbulent diffusion (turbdiff)
!!
!! These subroutines are called by read_atmo_namelists and do the turbulent
!! diffusion setup (for turbdiff).
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2011-09-19)
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
MODULE mo_turbdiff_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_control,      ONLY: is_restart_run
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_dom
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_turbdiff_config,     ONLY: turbdiff_config 

  USE mo_data_turbdiff,       ONLY: &
    & itype_tran, itype_sher, itype_wcld, itype_synd, &
    & imode_tran, imode_turb, icldm_tran, icldm_turb, &
    & ltkesso, ltkecon, lexpcor, ltmpcor, lprfcor, lnonloc, lcpfluc, lsflcnd, limpltkediff, &
    & tur_len, pat_len, a_stab, tkhmin, tkmmin, c_diff, &
    & rlam_heat, rlam_mom, rat_sea, tkesmot, impl_s, impl_t
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_turbdiff_namelist


  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  !----------------------------------!
  ! additional turbdiff_nml namelist variables  !
  !----------------------------------!

  LOGICAL :: lconst_z0   ! TRUE: horizontally homogeneous roughness length 
                         ! (for idealized testcases)

  REAL(wp) :: const_z0   ! horizontally homogeneous roughness length 
                         ! (for idealized testcases)

  NAMELIST/turbdiff_nml/ &
    & itype_tran, itype_sher, itype_wcld, itype_synd, &
    & imode_tran, imode_turb, icldm_tran, icldm_turb, &
    & ltkesso, ltkecon, lexpcor, ltmpcor, lprfcor, lnonloc, lcpfluc, lsflcnd, limpltkediff, &
    & tur_len, pat_len, a_stab, tkhmin, tkmmin, c_diff, &
    & rlam_heat, rlam_mom, rat_sea, tkesmot, impl_s, impl_t, &
!   additional namelist parameters:
    & lconst_z0, const_z0

CONTAINS

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for turbulent diffusion. 
  !!
  !! This subroutine 
  !! - reads the Namelist for turbulent diffusion
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)   
  !!
  !! @par Revision History
  !!  by Daniel Reinert, DWD (2011-09-19)
  !!
  SUBROUTINE read_turbdiff_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg          !< patch loop index
    INTEGER :: iunit

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_turbdiff_nml: read_turbdiff_nml'

    ! 0. default settings of internal turbdiff namelist variables
    !    is done by initialization in MODULE 'mo_data_turbdiff'

    !-----------------------
    ! 1. default settings of additional namelist variables:
    !-----------------------

    lconst_z0    =.FALSE. ! horizontally homogeneous roughness length 
                          ! (for idealized testcases)

    const_z0     = 0.001_wp ! horizontally homogeneous roughness length
                            ! (for idealized testcases)

    !-----------------------
    ! 1.b temorarily overwrite default settings of namelist variables:
    !-----------------------

!DR This is a dirty hack, since it implies that the namelist nwp_phy_nml
!DR must be read before turbdiff_nml. Any cross dependency between namelist-modules should 
!DR be avoided!! 
    IF ( ANY( (/10,11/)==atm_phy_nwp_config(1)%inwp_turb ) ) THEN
      imode_tran = 0       ! mode of surface-atmosphere transfer
      icldm_tran = 2       ! mode of cloud representation in transfer parametr.
    END IF
  
    IF ( ANY( (/10,12/)==atm_phy_nwp_config(1)%inwp_turb ) ) THEN
      imode_turb = 1       ! mode of turbulent diffusion parametrization
    END IF
  
    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (is_restart_run()) THEN
      funit = open_and_restore_namelist('turbdiff_nml')
      READ(funit,NML=turbdiff_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('turbdiff_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults() 
      WRITE(iunit, turbdiff_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, turbdiff_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, turbdiff_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------





    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg= 0,max_dom
      turbdiff_config(jg)%itype_tran   = itype_tran
      turbdiff_config(jg)%imode_tran   = imode_tran
      turbdiff_config(jg)%icldm_tran   = icldm_tran
      turbdiff_config(jg)%imode_turb   = imode_turb
      turbdiff_config(jg)%icldm_turb   = icldm_turb
      turbdiff_config(jg)%itype_sher   = itype_sher
      turbdiff_config(jg)%ltkesso      = ltkesso
      turbdiff_config(jg)%ltkecon      = ltkecon
      turbdiff_config(jg)%lexpcor      = lexpcor
      turbdiff_config(jg)%ltmpcor      = ltmpcor
      turbdiff_config(jg)%lprfcor      = lprfcor
      turbdiff_config(jg)%lnonloc      = lnonloc
      turbdiff_config(jg)%lcpfluc      = lcpfluc
      turbdiff_config(jg)%lsflcnd      = lsflcnd
      turbdiff_config(jg)%limpltkediff = limpltkediff
      turbdiff_config(jg)%itype_wcld   = itype_wcld
      turbdiff_config(jg)%itype_synd   = itype_synd
      turbdiff_config(jg)%tur_len      = tur_len
      turbdiff_config(jg)%pat_len      = pat_len
      turbdiff_config(jg)%a_stab       = a_stab
      turbdiff_config(jg)%tkhmin       = tkhmin
      turbdiff_config(jg)%tkmmin       = tkmmin
      turbdiff_config(jg)%c_diff       = c_diff
      turbdiff_config(jg)%rlam_heat    = rlam_heat
      turbdiff_config(jg)%rlam_mom     = rlam_mom
      turbdiff_config(jg)%rat_sea      = rat_sea
      turbdiff_config(jg)%tkesmot      = tkesmot
      turbdiff_config(jg)%impl_s       = impl_s
      turbdiff_config(jg)%impl_t       = impl_t

      turbdiff_config(jg)%lconst_z0    = lconst_z0
      turbdiff_config(jg)%const_z0     = const_z0
    ENDDO

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=turbdiff_nml)                    
      CALL store_and_close_namelist(funit, 'turbdiff_nml')             
    ENDIF

    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=turbdiff_nml)


  END SUBROUTINE read_turbdiff_namelist

END MODULE mo_turbdiff_nml
