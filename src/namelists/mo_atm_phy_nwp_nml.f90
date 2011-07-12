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
  USE mo_namelist,            ONLY: position_nml, POSITIONED
  USE mo_mpi,                 ONLY: p_pe, p_io
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_master_nml,          ONLY: lrestart
 
  USE mo_model_domain,        ONLY: t_patch

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

!  PRIVATE

   !
   ! user defined calling intervals
   !
  REAL(wp) :: nml_dt_conv(max_dom)   !> field element for convection
  REAL(wp) :: nml_dt_ccov(max_dom)   !! field element for subscale cloud cover
  REAL(wp) :: nml_dt_rad(max_dom)    !! "-"                     radiation
  REAL(wp) :: nml_dt_radheat(max_dom)!! "-" rad. heating from radiative fluxes with updated cosmu0 
  REAL(wp) :: nml_dt_sso(max_dom)    !! "-"  for subscale orographic gravity waves
  REAL(wp) :: nml_dt_gwd(max_dom)    !! "-"  for subscale gravity waves
  REAL(wp) :: nml_dt_gscp(max_dom)   !! field element for microphysics
  REAL(wp) :: nml_dt_turb(max_dom)   !! field element for turbulence
  REAL(wp) :: nml_dt_sfc(max_dom)    !! field element for surface
  REAL(wp) :: nml_dt_satad(max_dom)  !! field element for sat. adjustment
  REAL(wp) :: nml_dt_update(max_dom) !! field element for tracer phys update


  INTEGER ::  nml_inwp_gscp        !> microphysics
  INTEGER ::  nml_inwp_satad       !! saturation adjustment
  INTEGER ::  nml_inwp_convection  !! convection
  INTEGER ::  nml_inwp_radiation   !! radiation
  INTEGER ::  nml_inwp_sso         !! sso
  INTEGER ::  nml_inwp_gwd         !! non-orographic gravity wave drag
  INTEGER ::  nml_inwp_cldcover    !! cloud cover
  INTEGER ::  nml_inwp_turb        !! turbulence
  INTEGER ::  nml_inwp_surface     !! surface including soil, ocean, ice,lake

  INTEGER :: jg


  INTEGER :: nml_imode_turb, nml_itype_wcld, nml_icldm_turb,nml_itype_tran
   
  LOGICAL :: nml_limpltkediff, nml_ltkesso, nml_lexpcor

  REAL(wp):: nml_tur_len, nml_pat_len, nml_a_stab,                &
    &        nml_tkhmin, nml_tkmmin, nml_c_diff,                  &
    &        nml_rlam_heat, nml_rlam_mom, nml_rat_sea

  LOGICAL ::       &
   &    nml_lseaice,    & !> forecast with sea ice model
   &    nml_llake,      & !! forecst with lake model FLake
   &    nml_lmelt     , & !! soil model with melting process
   &    nml_lmelt_var     !! freezing temperature dependent on water content

!> Variables for hydci_pp
! --------------------------------------

  REAL(wp)::  nml_qi0, nml_qc0

  NAMELIST /nwp_phy_ctl/ nml_inwp_gscp, nml_inwp_satad, nml_inwp_convection,  &
    &                  nml_inwp_radiation, nml_inwp_sso, nml_inwp_cldcover, &
    &                  nml_inwp_gwd,                                        &
    &                  nml_inwp_turb, nml_inwp_surface,                     &
    &                  nml_dt_conv, nml_dt_ccov, nml_dt_rad,                &
    &                  nml_dt_radheat,                                      &
    &                  nml_dt_sso, nml_dt_gscp, nml_dt_satad,               &
    &                  nml_dt_turb, nml_dt_sfc, nml_dt_gwd,                 & 
    &                  nml_imode_turb,                                      &
    &                  nml_limpltkediff, nml_ltkesso, nml_lexpcor,          &
    &                  nml_tur_len, nml_pat_len, nml_a_stab,                &
    &                  nml_tkhmin, nml_tkmmin, nml_c_diff,                  &
    &                  nml_itype_wcld, nml_icldm_turb,                      &
    &                  nml_itype_tran, nml_rlam_heat, nml_rlam_mom, nml_rat_sea,&
    &                  nml_qi0, nml_qc0


   PUBLIC :: read_nwp_phy_namelist 

 CONTAINS


  !-------------------------------------------------------------------------
  !
  !>
  !! Read physics Namelist
  !!
  !!
  !! @par Revision History
  !! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
  !!
 SUBROUTINE read_inwp_nml

  INTEGER :: jg
  INTEGER :: i_status

  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: routine =  &
                              'mo_atm_phy_nwp_nml'
  !-----------------------------------------------------------------------

    CALL position_nml ('nwp_phy_ctl', status=i_status)
    !
    SELECT CASE (i_status)
    CASE (POSITIONED)
      READ (nnml, nwp_phy_ctl)
    END SELECT
  !  write the contents of the namelist to an ASCII file
    IF(p_pe == p_io) WRITE(nnml_output,nml=nwp_phy_ctl)

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
  SUBROUTINE read_nwp_phy_namelist (p_patch)
    !
    TYPE(t_patch), OPTIONAL, INTENT(IN) :: p_patch(:)
    INTEGER :: istat, funit

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_atm_phy_nwp_nml: read_nwp_phy_namelist'

    !-----------------------------------------------------------------------

    !-----------------------!
    ! 1. default settings   !
    !-----------------------!
    nml_inwp_gscp       = 0           !> 0 = no microphysics
    nml_inwp_satad      = 0           !> 1 = saturation adjustment on
    nml_inwp_convection = 0           !> 0 = no convection
    nml_inwp_radiation  = 0           !> 0 = no radiation
    nml_inwp_sso        = 0           !> 0 = no sso
    nml_inwp_gwd        = 0           !> 0 = no gwd, 1= IFS gwd scheme
    nml_inwp_cldcover   = 1           !> 1 = use grid-scale clouds for radiation
    nml_inwp_turb       = 0           !> 0 = no turbulence,1= cosmo/turbdiff,2=echam/vdiff
    nml_inwp_surface    = 0           !> 0 = no surface, 1 =  cosmo surface

    DO jg=1, max_dom
      nml_dt_conv (jg) = 600._wp      !seconds
      nml_dt_ccov (jg) = nml_dt_conv(jg)  !presently not used; cloud cover is synchronized with radiation
      nml_dt_rad  (jg) = 1800._wp     !seconds
      nml_dt_sso  (jg) = 3600._wp     !seconds
      nml_dt_gwd  (jg) = 3600._wp     !seconds
      nml_dt_gscp (jg) = 100._wp      !seconds
      nml_dt_turb (jg) = 100._wp      !seconds
      nml_dt_sfc  (jg) = 100._wp      !seconds
      nml_dt_satad(jg) = 100._wp      !seconds
      nml_dt_update(jg) =  nml_dt_satad (jg)
      nml_dt_radheat(jg)=  nml_dt_update(jg)
    ENDDO

  !> KF  current settings to get NWP turbulence running
    nml_lseaice    = lseaice
    nml_llake      = llake

    nml_imode_turb   = imode_turb 
    nml_limpltkediff = limpltkediff
    nml_ltkesso      = ltkesso
    nml_lexpcor      = lexpcor
    nml_tur_len      = tur_len
    nml_pat_len      = pat_len
    nml_a_stab       = a_stab
    nml_tkhmin       = tkhmin
    nml_tkmmin       = tkmmin
    nml_c_diff       = c_diff
    nml_itype_wcld   = itype_wcld 
    nml_icldm_turb   = icldm_turb
    nml_itype_tran   = itype_tran
    nml_rlam_heat    = rlam_heat 
    nml_rlam_mom     = rlam_mom
    nml_rat_sea      = rat_sea

    nml_qi0 = 0.0_wp 
    nml_qc0 = 0.0_wp 

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (lrestart) THEN
      funit = open_and_restore_namelist('nwp_phy_ctl')
      READ(funit,NML=nwp_phy_ctl)
      CALL close_tmpfile(funit)
    END IF


    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL position_nml ('nwp_phy_ctl', status=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, nwp_phy_ctl)
    END SELECT


    !----------------------------------------------------
    ! 4. Fill the configuration state
    !----------------------------------------------------
    DO jg= 1,max_dom
      atm_phy_nwp_config(jg)%inwp_gscp       = nml_inwp_gscp 
      atm_phy_nwp_config(jg)%inwp_satad      = nml_inwp_satad
      atm_phy_nwp_config(jg)%inwp_convection = nml_inwp_convection
      atm_phy_nwp_config(jg)%inwp_radiation  = nml_inwp_radiation
      atm_phy_nwp_config(jg)%inwp_sso        = nml_inwp_sso
      atm_phy_nwp_config(jg)%inwp_gwd        = nml_inwp_gwd     
      atm_phy_nwp_config(jg)%inwp_cldcover   = nml_inwp_cldcover
      atm_phy_nwp_config(jg)%inwp_turb       = nml_inwp_turb
      atm_phy_nwp_config(jg)%inwp_surface    = nml_inwp_surface
      atm_phy_nwp_config(jg)% dt_conv        = nml_dt_conv (jg) 
      atm_phy_nwp_config(jg)% dt_ccov        = nml_dt_ccov (jg)
      atm_phy_nwp_config(jg)% dt_rad         = nml_dt_rad  (jg)
      atm_phy_nwp_config(jg)% dt_radheat     = nml_dt_radheat(jg)
      atm_phy_nwp_config(jg)% dt_sso         = nml_dt_sso   (jg)
      atm_phy_nwp_config(jg)% dt_gwd         = nml_dt_gwd   (jg)
      atm_phy_nwp_config(jg)% dt_gscp        = nml_dt_gscp  (jg)
      atm_phy_nwp_config(jg)% dt_turb        = nml_dt_turb  (jg)
      atm_phy_nwp_config(jg)% dt_sfc         = nml_dt_sfc   (jg)
      atm_phy_nwp_config(jg)% dt_satad       = nml_dt_satad (jg)
      atm_phy_nwp_config(jg)% dt_update      = nml_dt_update(jg)

      atm_phy_nwp_config(jg)%lseaice         = nml_lseaice 
      atm_phy_nwp_config(jg)%llake           = nml_llake
      atm_phy_nwp_config(jg)%imode_turb      = nml_imode_turb  
      atm_phy_nwp_config(jg)%limpltkediff    = nml_limpltkediff
      atm_phy_nwp_config(jg)%ltkesso         = nml_ltkesso
      atm_phy_nwp_config(jg)%lexpcor         = nml_lexpcor
      atm_phy_nwp_config(jg)%tur_len         = nml_tur_len
      atm_phy_nwp_config(jg)%pat_len         = nml_pat_len
      atm_phy_nwp_config(jg)%a_stab          = nml_a_stab
      atm_phy_nwp_config(jg)%tkhmin          = nml_tkhmin
      atm_phy_nwp_config(jg)%tkmmin          = nml_tkmmin
      atm_phy_nwp_config(jg)%c_diff          = nml_c_diff
      atm_phy_nwp_config(jg)%itype_wcld      = nml_itype_wcld
      atm_phy_nwp_config(jg)%icldm_turb      = nml_icldm_turb
      atm_phy_nwp_config(jg)%qi0             = nml_qi0 
      atm_phy_nwp_config(jg)%qc0             = nml_qc0 

    ENDDO

    !-----------------------------------------------------
    ! 5. Store the namelist for restart
    !-----------------------------------------------------
    funit = open_tmpfile()
    WRITE(funit,NML=nwp_phy_ctl)                    
    CALL store_and_close_namelist(funit, 'nwp_phy_ctl') 


    ! 6. write the contents of the namelist to an ASCII file
    !
    IF(p_pe == p_io) WRITE(nnml_output,nml=nwp_phy_ctl)


  END SUBROUTINE read_nwp_phy_namelist


END MODULE mo_atm_phy_nwp_nml

