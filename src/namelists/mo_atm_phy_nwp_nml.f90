!>
!!  Namelist for NWP physics
!!
!!  these Subroutines are called by control model and construct the
!!  physics composition
!!
!! @author <Kristina Froehlich, DWD>
!!
!!
!! @par Revision History
!! First implementation by Kristina Froehlich, DWD (2010-06-20>)
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
MODULE mo_atm_phy_nwp_nml

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: max_dom,MAX_CHAR_LENGTH
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: p_pe, p_io
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_nml,          ONLY: lrestart

  USE mo_io_restart_namelist, ONLY: open_tmpfile, store_and_close_namelist,  &
    &                               open_and_restore_namelist, close_tmpfile

  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_data_turbdiff,       ONLY: imode_turb,                              &
    &                               limpltkediff, ltkesso, lexpcor,          &
    &                               tur_len, pat_len, a_stab,                &
    &                               tkhmin, tkmmin, c_diff,                  &
    &                               itype_wcld, icldm_turb,                  &
    &                               itype_tran, rlam_heat, rlam_mom, rat_sea,&
    &                               llake, lseaice

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_nwp_phy_namelist 

   !
   ! user defined calling intervals
   !
  REAL(wp) :: dt_conv(max_dom)   !> field element for convection
  REAL(wp) :: dt_ccov(max_dom)   !! field element for subscale cloud cover
  REAL(wp) :: dt_rad(max_dom)    !! "-"                     radiation
  REAL(wp) :: dt_radheat(max_dom)!! "-" rad. heating from radiative fluxes with updated cosmu0 
  REAL(wp) :: dt_sso(max_dom)    !! "-"  for subscale orographic gravity waves
  REAL(wp) :: dt_gwd(max_dom)    !! "-"  for subscale gravity waves
  REAL(wp) :: dt_gscp(max_dom)   !! field element for microphysics
  REAL(wp) :: dt_satad(max_dom)  !! field element for sat. adjustment
  REAL(wp) :: dt_turb(max_dom)   !! field element for turbulence
  REAL(wp) :: dt_sfc(max_dom)    !! field element for surface
  REAL(wp) :: dt_update(max_dom) !! field element for tracer phys update


  INTEGER  :: inwp_convection    !! convection
  INTEGER  :: inwp_cldcover      !! cloud cover
  INTEGER  :: inwp_radiation     !! radiation
  INTEGER  :: inwp_sso           !! sso
  INTEGER  :: inwp_gwd           !! non-orographic gravity wave drag
  INTEGER  :: inwp_gscp          !! microphysics
  INTEGER  :: inwp_satad         !! saturation adjustment
  INTEGER  :: inwp_turb          !! turbulence
  INTEGER  :: inwp_surface       !! surface including soil, ocean, ice,lake

!  INTEGER  :: imode_turb, itype_wcld, icldm_turb, itype_tran
!   
!  LOGICAL  :: limpltkediff, ltkesso, lexpcor
!
!  REAL(wp) :: tur_len, pat_len, a_stab,                &
!    &         tkhmin, tkmmin, c_diff,                  &
!    &         rlam_heat, rlam_mom, rat_sea
!
!  LOGICAL  :: lseaice  !> forecast with sea ice model
!  LOGICAL  :: llake    !! forecst with lake model FLake
!
!  REAL(wp)::  qi0, qc0 !! variables for hydci_pp

  NAMELIST /nwp_phy_nml/ inwp_convection, inwp_cldcover,           &
    &                    inwp_radiation, inwp_sso, inwp_gwd,       &
    &                    inwp_gscp, inwp_satad,                    &
    &                    inwp_turb, inwp_surface,                  &
    &                    dt_conv, dt_ccov,                         &
    &                    dt_rad, dt_radheat,                       &
    &                    dt_sso, dt_gwd,                           &
    &                    dt_gscp, dt_satad,                        &
    &                    dt_turb, dt_sfc ! ,                          & 
!    &                    imode_turb,                               &
!    &                    limpltkediff, ltkesso, lexpcor,           &
!    &                    tur_len, pat_len, a_stab,                 &
!    &                    tkhmin, tkmmin, c_diff,                   &
!    &                    itype_wcld, icldm_turb,                   &
!    &                    itype_tran, rlam_heat, rlam_mom, rat_sea, &
!    &                    qi0, qc0



CONTAINS
  !-------------------------------------------------------------------------
  !>
  !! Read physics Namelist
  !!
  !! @par Revision History
  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
  !!
 SUBROUTINE read_inwp_nml

  INTEGER :: i_status

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &  routine = 'mo_atm_phy_nwp_nml:read_inwp_nml'

    CALL position_nml ('nwp_phy_nml', status=i_status)
    !
    SELECT CASE (i_status)
    CASE (POSITIONED)
      READ (nnml, nwp_phy_nml)
    END SELECT
  !  write the contents of the namelist to an ASCII file
    IF(p_pe == p_io) WRITE(nnml_output,nml=nwp_phy_nml)

 END SUBROUTINE read_inwp_nml



  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Read Namelist for NWP physics. 
  !!
  !! This subroutine 
  !! - reads the Namelist for NWP physics
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
  SUBROUTINE read_nwp_phy_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit, jg
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_atm_phy_nwp_nml:read_nwp_phy_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------
    inwp_gscp       = 0           !> 0 = no microphysics
    inwp_satad      = 0           !> 1 = saturation adjustment on
    inwp_convection = 0           !> 0 = no convection
    inwp_radiation  = 0           !> 0 = no radiation
    inwp_sso        = 0           !> 0 = no sso
    inwp_gwd        = 0           !> 0 = no gwd, 1= IFS gwd scheme
    inwp_cldcover   = 1           !> 1 = use grid-scale clouds for radiation
    inwp_turb       = 0           !> 0 = no turbulence,1= cosmo/turbdiff,2=echam/vdiff
    inwp_surface    = 0           !> 0 = no surface, 1 =  cosmo surface

    DO jg=1, max_dom
      dt_conv (jg) = 600._wp      !seconds
      dt_ccov (jg) = dt_conv(jg)  !presently not used; cloud cover is synchronized with radiation
      dt_rad  (jg) = 1800._wp     !seconds
      dt_sso  (jg) = 3600._wp     !seconds
      dt_gwd  (jg) = 3600._wp     !seconds
      dt_gscp (jg) = 100._wp      !seconds
      dt_turb (jg) = 100._wp      !seconds
      dt_sfc  (jg) = 100._wp      !seconds
      dt_satad(jg) = 100._wp      !seconds
      dt_update(jg) =  dt_satad (jg)
      dt_radheat(jg)=  dt_update(jg)
    ENDDO

!    qi0 = 0.0_wp 
!    qc0 = 0.0_wp 

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('nwp_phy_nml')
      READ(funit,NML=nwp_phy_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('nwp_phy_nml', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, nwp_phy_nml)
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------
    DO jg=1,max_dom
      atm_phy_nwp_config(jg)%inwp_convection = inwp_convection
      atm_phy_nwp_config(jg)%inwp_cldcover   = inwp_cldcover
      atm_phy_nwp_config(jg)%inwp_radiation  = inwp_radiation
      atm_phy_nwp_config(jg)%inwp_sso        = inwp_sso
      atm_phy_nwp_config(jg)%inwp_gwd        = inwp_gwd     
      atm_phy_nwp_config(jg)%inwp_gscp       = inwp_gscp 
      atm_phy_nwp_config(jg)%inwp_satad      = inwp_satad
      atm_phy_nwp_config(jg)%inwp_turb       = inwp_turb
      atm_phy_nwp_config(jg)%inwp_surface    = inwp_surface
      
      atm_phy_nwp_config(jg)% dt_conv        = dt_conv (jg) 
      atm_phy_nwp_config(jg)% dt_ccov        = dt_ccov (jg)
      atm_phy_nwp_config(jg)% dt_rad         = dt_rad  (jg)
      atm_phy_nwp_config(jg)% dt_radheat     = dt_radheat(jg)
      atm_phy_nwp_config(jg)% dt_sso         = dt_sso   (jg)
      atm_phy_nwp_config(jg)% dt_gwd         = dt_gwd   (jg)
      atm_phy_nwp_config(jg)% dt_gscp        = dt_gscp  (jg)
      atm_phy_nwp_config(jg)% dt_satad       = dt_satad (jg)
      atm_phy_nwp_config(jg)% dt_turb        = dt_turb  (jg)
      atm_phy_nwp_config(jg)% dt_sfc         = dt_sfc   (jg)
      atm_phy_nwp_config(jg)% dt_update      = dt_update(jg)

!      atm_phy_nwp_config(jg)%lseaice         = lseaice 
!      atm_phy_nwp_config(jg)%llake           = llake
!      atm_phy_nwp_config(jg)%imode_turb      = imode_turb  
!      atm_phy_nwp_config(jg)%limpltkediff    = limpltkediff
!      atm_phy_nwp_config(jg)%ltkesso         = ltkesso
!      atm_phy_nwp_config(jg)%lexpcor         = lexpcor
!      atm_phy_nwp_config(jg)%tur_len         = tur_len
!      atm_phy_nwp_config(jg)%pat_len         = pat_len
!      atm_phy_nwp_config(jg)%a_stab          = a_stab
!      atm_phy_nwp_config(jg)%tkhmin          = tkhmin
!      atm_phy_nwp_config(jg)%tkmmin          = tkmmin
!      atm_phy_nwp_config(jg)%c_diff          = c_diff
!      atm_phy_nwp_config(jg)%itype_wcld      = itype_wcld
!      atm_phy_nwp_config(jg)%icldm_turb      = icldm_turb
!      atm_phy_nwp_config(jg)%qi0             = qi0 
!      atm_phy_nwp_config(jg)%qc0             = qc0 

    ENDDO

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=nwp_phy_nml)                    
    CALL store_and_close_namelist(funit, 'nwp_phy_nml') 

    ! 6. write the contents of the namelist to an ASCII file
    !
    IF(p_pe == p_io) WRITE(nnml_output,nml=nwp_phy_nml)

  END SUBROUTINE read_nwp_phy_namelist


END MODULE mo_atm_phy_nwp_nml

