!>
!! This module is the interface between nwp_nh_interface to the radiation schemes
!! (RRTM or Ritter-Geleyn).
!!
!! @author Thorsten Reinhardt, AGeoBw, Offenbach
!!
!! @par Revision History
!! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
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

MODULE mo_nwp_rad_interface

  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_datetime,             ONLY: t_datetime
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_parallel_config,      ONLY: nproma, parallel_radiation_mode
  USE mo_impl_constants,       ONLY: max_char_length, MODIS 
  USE mo_kind,                 ONLY: wp
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_model_domain,         ONLY: t_patch
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_radiation_config,     ONLY: albedo_type
  USE mo_radiation,            ONLY: pre_radiation_nwp_steps
  USE mo_nwp_rrtm_interface,   ONLY: nwp_rrtm_radiation,             &
    &                                nwp_rrtm_radiation_reduced,     &
    &                                nwp_rrtm_radiation_repartition, &
    &                                nwp_rrtm_ozon_aerosol
  USE mo_nwp_rg_interface,     ONLY: nwp_rg_radiation,               &
    &                                nwp_rg_radiation_reduced
  USE mo_albedo,               ONLY: sfc_albedo, sfc_albedo_modis

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER:: &
    &  version = '$Id$'


  PUBLIC :: nwp_radiation
  

 CONTAINS
  
  !---------------------------------------------------------------------------------------
  !>
  !! This subroutine is the interface between nwp_nh_interface to the radiation schemes.
  !! Depending on inwp_radiation, it can call RRTM (1) or Ritter-Geleyn (2).
  !!
  !! @par Revision History
  !! Initial release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-13)
  !!
  SUBROUTINE nwp_radiation ( lredgrid, p_sim_time, datetime, pt_patch,pt_par_patch, &
    & ext_data, lnd_diag, pt_prog, pt_diag, prm_diag, lnd_prog, wtr_prog, p_metrics )

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER::  &
      &  routine = 'mo_nwp_rad_interface:nwp_radiation'
    
    LOGICAL, INTENT(in)         :: lredgrid        !< use reduced grid for radiation

    REAL(wp),INTENT(in)         :: p_sim_time

    TYPE(t_datetime)            ,INTENT(in) :: datetime
    TYPE(t_patch)       , TARGET,INTENT(in) :: pt_patch     !<grid/patch info.
    TYPE(t_patch)       , TARGET,INTENT(in) :: pt_par_patch !<grid/patch info (parent grid)
    TYPE(t_external_data)       ,INTENT(inout):: ext_data
    TYPE(t_lnd_diag)            ,INTENT(in)   :: lnd_diag   !<diag vars for sfc
    TYPE(t_nh_prog)     , TARGET,INTENT(inout):: pt_prog    !<the prognostic variables
    TYPE(t_nh_diag)     , TARGET,INTENT(inout):: pt_diag    !<the diagnostic variables
    TYPE(t_nh_metrics)          ,INTENT(in)   :: p_metrics
    TYPE(t_nwp_phy_diag)        ,INTENT(inout):: prm_diag
    TYPE(t_lnd_prog)            ,INTENT(inout):: lnd_prog   ! time level new
    TYPE(t_wtr_prog)            ,INTENT(   in):: wtr_prog   ! time level new

    REAL(wp) :: &
      & zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)

    
    INTEGER :: jg

    REAL(wp):: zsct        ! solar constant (at time of year)
    REAL(wp):: cosmu0_dark ! minimum cosmu0, for smaller values no shortwave calculations



    ! patch ID
    jg = pt_patch%id



    !-------------------------------------------------------------------------
    !> Radiation setup
    !-------------------------------------------------------------------------

    IF (atm_phy_nwp_config(jg)%inwp_radiation == 1 ) THEN
      ! RRTM
      ! In radiative transfer routine RRTM skips all points with cosmu0<=0. That's why 
      ! points to be skipped need to be marked with a value <=0
      cosmu0_dark = -1.e-9_wp  ! minimum cosmu0, for smaller values no shortwave calculations
    ELSE
      ! Ritter-Geleyn
      ! Skipping of points is performed on block- rather than cell-level. I.e. if a block 
      ! contains at least 1 point with cosmu0>1.E-8, radiatve transfer is computed for 
      ! the entire block. Therefore cosmu0_dark = -1.e-9_wp does not work here (crashes).
      ! For all points cosmu0 must be <0.
      cosmu0_dark =  1.e-9_wp   ! minimum cosmu0, for smaller values no shortwave calculations
    ENDIF


    ! Calculation of zenith angle optimal during dt_rad.
    ! (For radheat, actual zenith angle is calculated separately.)
    CALL pre_radiation_nwp_steps (                        &
      & kbdim        = nproma,                            & !in
      & cosmu0_dark  = cosmu0_dark,                       & !in
      & p_inc_rad    = atm_phy_nwp_config(jg)%dt_rad,     & !in
      & p_inc_radheat= atm_phy_nwp_config(jg)%dt_fastphy, & !in
      & p_sim_time   = p_sim_time,                        & !in
      & pt_patch     = pt_patch,                          & !in
      & zsmu0        = prm_diag%cosmu0(:,:),              & !out
      & zsct         = zsct                               ) !out, optional




    ! Compute tile-based and aggregated surface-albedo
    !
    IF ( albedo_type == MODIS ) THEN
      ! MODIS albedo
      CALL sfc_albedo_modis(pt_patch, ext_data, lnd_prog, wtr_prog, lnd_diag, prm_diag)
    ELSE
      ! albedo based on tabulated bare soil values
      CALL sfc_albedo(pt_patch, ext_data, lnd_prog, wtr_prog, lnd_diag, prm_diag)
    ENDIF



    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------
    ! RRTM
    !
    IF (atm_phy_nwp_config(jg)%inwp_radiation == 1 ) THEN
       
      CALL nwp_rrtm_ozon_aerosol ( p_sim_time, datetime, pt_patch, ext_data, &
        & pt_diag,prm_diag,zaeq1,zaeq2,zaeq3,zaeq4,zaeq5 )
    
      IF ( .NOT. lredgrid ) THEN

        SELECT CASE(parallel_radiation_mode(pt_patch%id))
        CASE(1) 
          CALL nwp_rrtm_radiation_repartition ( pt_patch, ext_data, &
            & zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,                    &
            & pt_diag, prm_diag, lnd_prog, p_metrics )
!         CASE(2)
!           CALL nwp_omp_rrtm_interface ( pt_patch, ext_data, &
!              &  lnd_diag, pt_diag, prm_diag, lnd_prog )
          
        CASE default
          CALL nwp_rrtm_radiation ( pt_patch, ext_data, &
            & zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,        &
            & pt_diag, prm_diag, lnd_prog, p_metrics )

       END SELECT
       
      ELSE 

        CALL nwp_rrtm_radiation_reduced ( pt_patch,pt_par_patch, ext_data, &
          & zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,                             &
          & pt_diag, prm_diag, lnd_prog, p_metrics )
          
      ENDIF

      RETURN
    ENDIF !inwp_radiation = 1 (RRTM)


    ! Ritter-Geleyn
    !
    IF ( atm_phy_nwp_config(jg)%inwp_radiation == 2 .AND. .NOT. lredgrid) THEN
    
      CALL nwp_rg_radiation ( p_sim_time, datetime, pt_patch, &
        & ext_data,pt_prog,pt_diag,prm_diag, lnd_prog, zsct )


    ELSEIF ( atm_phy_nwp_config(jg)%inwp_radiation == 2 .AND. lredgrid) THEN

      CALL nwp_rg_radiation_reduced ( p_sim_time, datetime, pt_patch,pt_par_patch, &
        & ext_data, pt_prog, pt_diag, prm_diag, lnd_prog, zsct )


    ENDIF !inwp_radiation = 2 (Ritter-Geleyn)

  END SUBROUTINE nwp_radiation


END MODULE mo_nwp_rad_interface

