!>
!! The main program unit of the sea-ice parameterization scheme for NWP. 
!! It contains two procedures, viz., 
!! SUBROUTINE seaice_init_nwp
!! that initializes the scheme and performs some consistency checks, 
!! and 
!! SUBROUTINE seaice_timestep_nwp
!! that advances prognostic variables of the sea-ice scheme one time step.
!!
!! The present sea-ice parameterization scheme is a bulk thermodynamic (no rheology) scheme 
!! intended for use in NWP and similar applications.  
!! The scheme is based on a self-similar parametric representation (assumed shape) of the
!! evolving temperature profile within the ice and on the integral heat budget of the ice slab. 
!! The scheme carries ordinary differential equations (in time) 
!! for the ice surface temperature and the ice thickness. 
!! An explicit Euler scheme is used for time advance.
!! In the current configuration of the scheme, snow over sea ice is not treated explicitly. 
!! The effect of snow above the ice is accounted for implicitely (parametrically) through 
!! an empirical temperature dependence of the ice surface albedo with respect to solar radiation. 
!! For the "sea water" type ICON grid boxes, the snow thickness is set to zero and
!! the snow surface temperature is set equal to the ice surface temperature
!! (both temperatures are set equal to the fresh-water freezing point if the ice is absent).
!! Prognostic equations for the ice thickness and the ice surface temperature are solved 
!! for the ICON grid boxes with the ice fraction 
!! (area fraction of a given model grid box of the type "sea water" that is covered by ice)
!! that exceeds a threshold value of 0.03. Otherwise, the grid box is treated as ice-free.
!! The ice fraction is determined on the basis of observational data
!! by the data assimilation scheme and is kept constant over the entire model forecast period. 
!! However, if the ice melts out during the forecast, the ice fraction is reset to zero. 
!! This is done within 
!! SUBROUTINE <name>, 
!! where the grid-box-mean (aggregated) surface temperature 
!! and grid-box-mean surface fluxes are computed.
!! If the ICON grid box is set ice-free during the initialization, 
!! no ice is created over the forecast period. 
!! If observational data indicate open water conditions for a given ICON grid box,
!! residual ice from the previous model run is removed, 
!! i.e. the ice thickness is set to zero and 
!! the ice surface temperature is set to the fresh-water freezing point. 
!! The newly formed ice has the surface temperature equal to the salt-water freezing point 
!! and the thickness of 0.5 m. The new ice is formed instantaneously 
!! if the data assimilation scheme indicates the presence of ice in a given ICON grid box
!! but there was no ice in that grid box during the previous model run. 
!! Prognostic ice thickness is limited by a maximum value of 3 m and a minimum value of 0.05 m. 
!! Constant values of the density, molecular heat conductivity, specific heat of ice, the latent
!! heat of fusion, and the salt-water freezing point are used.
!!
!! A detailed description of the sea ice scheme is given in
!! Mironov, D., B. Ritter, J.-P. Schulz, M. Buchhold, M. Lange, and E. Machulskaya, 2012:
!! Parameterization of sea and lake ice in numerical weather prediction models
!! of the German Weather Service.
!! Tellus A, 64, 17330. doi:10.3402/tellusa.v64i0.17330
!!
!! The present seas ice scheme (with minor modifications) 
!! is also implemented into the NWP models GME and COSMO 
!! (see Mironov et al. 2012, for details).
!!
!!
!! @author Dmitrii Mironov, DWD. 
!!
!! @par Revision History
!! Initial release by Dmitrii Mironov, DWD (2012-07-24)
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

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

! Lines embraced with "!_tmp>" and "!_tmp<" contain temporary parts of the code.
! Lines embraced/marked with "!_dev>" and "!_dev<" may be replaced
! as improved formulations are developed and tested.
! Lines embraced/marked with "!_cdm>" and "!_cdm<" are DM's comments that may be helpful to a user.
! Lines embraced/marked with "!_dbg>" and "!_dbg<" are used for debugging purposes only.
! Lines starting with "!_nu" are not used. 

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

MODULE mo_seaice_nwp

  USE mo_kind, ONLY:      &
                   &  wp  !< KIND-type parameter for real variables

  USE mo_exception, ONLY:           &          
                        &  finish , &  !< external procedure, finishes model run and reports the reason  
                        &  message     !< external procedure, sends a message (error, warning, etc.)

!_cdm>
! Note that ki is equal to 2.03 in ICON, but is 2.29 in COSMO and GME. 
!_cdm<
  USE mo_physical_constants, ONLY:                      &
                                 &  tf_fresh => tmelt , &  !< fresh-water freezing point [K]
                                 &              alf   , &  !< latent heat of fusion [J/kg]
                                 &              rhoi  , &  !< density of ice [kg/m^3]
                                 &              ci    , &  !< specific heat of ice [J/(kg K)]
                                 &              ki         !< molecular heat conductivity of ice [J/(m s K)]  

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

!_cdm>
! The value of the salt-water freezing point is the same as in GME and COSMO (-1.7 dgr C).
! Note that a different value (-1.8 dgr C) is defined in "mo_physical_constants".
!_cdm<
  REAL (wp), PARAMETER ::                             &
                       &  tf_salt      = 271.45_wp  , &  !< salt-water freezing point [K] 
                       &  frsi_min     = 0.03_wp    , &  !< minimum sea-ice fraction [-]
                       &  hice_min     = 0.05_wp    , &  !< minimum sea-ice thickness [m]
                       &  hice_max     = 3.0_wp     , &  !< maximum sea-ice thickness [m]
                       &  hice_new     = 0.5_wp     , &  !< thickness of the newly formed sea ice [m]
                       &  csi_lin      = 0.5_wp     , &  !< shape factor for linear temperature profile [-]
                       &  phiipr0_lin  = 1.0_wp     , &  !< derivative (at zeta=0) of the linear 
                                                         !< temperature profile shape function [-]
                       &  csidp_nlin   = 2.0_wp     , &  !< disposable parameter for non-linear 
                                                         !< temperature profile shape factor [-]
                       &  csidp_nlin_d =              &  !< derived disposable parameter for non-linear 
                                                         !< temperature profile shape factor [-]
                       &  (1._wp+csidp_nlin)/12._wp , &
                       &  cmaxearg     = 1.0E+02_wp , &  !< maximum value of the EXP function argument (security constant) [-]
                       &  csmall       = 1.0E-05_wp      !< small number (security constant) [-]
 
!_nu  <type>, PARAMETER :: <parameter> !<  <Parameter description>

  !>
  !! Optical characteristics of sea ice. 
  !!
  !! A storage for an n-band approximation of the exponential decay law 
  !! for the flux of solar radiation is allocated.
  !! A maximum value of the wave-length bands is currently set equal to two.
  !!
  
  INTEGER, PARAMETER ::                      &
                     &  nband_optic_max = 2  !< maximum number of wave-length bands in the decay law
                                             !< for the solar radiation flux [-]

  TYPE opticpar_seaice
    INTEGER ::                                       &
            &  nband_optic                              !< number of wave-length bands [-]
    REAL (wp) ::                                     &
              &  frac_optic(nband_optic_max)       , &  !< fractions of total solar radiation flux for different bands [-]
              &  extincoef_optic(nband_optic_max)       !< extinction coefficients for different bands [-]
  END TYPE opticpar_seaice

  ! One-band approximation for opaque sea ice.
  ! The use of large extinction coefficient prevents the penetration of solar radiation
  ! into the ice interior, i.e. the volumetric character of the solar radiation heating is ignored.
  TYPE (opticpar_seaice), PARAMETER ::                             &
                                    &  opticpar_seaice_opaque =    & 
                                    &  opticpar_seaice(1,          &
                                    &  (/1._wp, 0._wp/),           &
                                    &  (/1.0E+07_wp, 1.E+10_wp/)) 

  ! Minimum values of the sea-ice fraction and of the sea-ice thickness are used outside 
  ! "mo_seaice_nwp" to compose an index list of grid boxes where sea ice is present. 
  PUBLIC ::                       &
         &  frsi_min            , & ! parameter
         &  hice_min            , & ! parameter 
         &  seaice_init_nwp     , & ! procedure     
         &  seaice_timestep_nwp     ! procedure  


!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

CONTAINS

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

  !>
  !! Prognostic variables of the sea-ice parameterization scheme are initialized 
  !! and some consistency checks are performed. 
  !!
  !! The procedure arguments are arrays (vectors) of the sea-ice scheme prognostic variables
  !! and of the sea-ice fraction. The vector length is equal to the number of ICON grid boxes
  !! (within a given block)
  !! where sea ice is present, i.e. where the sea-ice fraction exceeds its minimum value.
  !! First, the sea-ice fraction is checked. 
  !! If a value less than a minimum threshold value is found 
  !! (indicating that a grid box is declared as partially ice-covered but no sea ice 
  !! should be present), an error message is sent and the model abort is called.
  !! Next, "new" ice is formed in the grid boxes where "old" value of the sea-ice thickness
  !! is less than a minimum threshold value (i.e. there was no ice at the end 
  !! of the the previous model run). The newly formed ice has the surface temperature equal 
  !! to the salt-water freezing point and the thickness of 0.5 m. 
  !! Then, the ice thickness is limited from above and below 
  !! by the maximum and minimum threshold values, 
  !! and the ice surface temperature is limited from above by the fresh-water freezing point. 
  !! These security measures (taken for each grid box where sea ice is present)
  !! are required to avoid non-allowable values of prognostic variables 
  !! that may occure due to a loss of accuracy during the model IO 
  !! (e.g. due to GRIB encoding and decoding).  
  !! Finally, the snow thickness is set to zero, and
  !! the snow surface temperature is set equal to the ice surface temperature 
  !! (recall that snow over sea ice is not treated explicitly).
  !! 
  !! 
  !! @par Revision History
  !! Initial release by Dmitrii Mironov, DWD (2012-07-24)
  !! Modification by Daniel Reinert, DWD (2012-11-05)
  !! - modified initialization procedure for the case that the sea-ice thickness 
  !!   field is not provided as input (i.e. when starting from IFS analysis)
  !!

  SUBROUTINE seaice_init_nwp (                                  & 
                          &  nsigb,                             &
                          &  frsi,                              &
                          &  t_seasfc,                          &
                          &  l_hice_in,                         &
                          &  tice_p, hice_p, tsnow_p, hsnow_p,  &
                          &  tice_n, hice_n, tsnow_n, hsnow_n   &
                          &  )

    IMPLICIT NONE

    ! Procedure arguments 

    INTEGER, INTENT(IN) ::        &
                        &  nsigb  !< Array (vector) dimension
                                  !< (equal to the number of grid boxes within a block 
                                  !< where the sea ice is present) 

    REAL(wp), DIMENSION(:), INTENT(IN)    ::       &
                                           &  frsi      !< sea-ice fraction [-]
 
    REAL(wp), DIMENSION(:), INTENT(IN)    ::       &
                                           &  t_seasfc  !< sea surface temperature (including sea-ice) [K]

    LOGICAL, INTENT(IN)                   ::       &
                                           &  l_hice_in !< Logical switch, if sea-ice thickness field is 
                                                        !< provided as input                     

    REAL(wp), DIMENSION(:), INTENT(INOUT) ::               &
                                              &  tice_p  , &  !< temperature of ice upper surface at previous time level [K] 
                                              &  hice_p  , &  !< ice thickness at previous time level [m] 
                                              &  tsnow_p , &  !< temperature of snow upper surface at previous time level [K] 
                                              &  hsnow_p , &  !< snow thickness at previous time level [m] 
                                              &  tice_n  , &  !< temperature of ice upper surface at new time level [K] 

                                              &  hice_n  , &  !< ice thickness at new time level [m] 
                                              &  tsnow_n , &  !< temperature of snow upper surface at new time level [K] 
                                              &  hsnow_n      !< snow thickness at new time level [m] 

    ! Local variables 

    INTEGER ::      &
            &  isi  !< DO loop index

    CHARACTER(len=256) ::           &
                       &  nameerr , &  !< name of procedure where an error occurs
                       &  texterr      !< error/warning message text

    LOGICAL ::             &
            &  lcallabort  !< logical switch, set .TRUE. if errors are encountered 
                           !< (used to call modell abort outside a DO loop)


    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------

    ! Logical switch, default value is .FALSE. 
    lcallabort = .FALSE.  

    ! Loop over grid boxes where sea ice is present
    GridBoxesWithSeaIce: DO isi=1, nsigb

      ! Check sea-ice fraction
      IF( frsi(isi) < frsi_min ) THEN 
        ! Sea-ice fraction less than a minimum threshold value is found
        ! Set logical switch 
        lcallabort = .TRUE.  
        ! Exit DO loop to call model abort
        EXIT GridBoxesWithSeaIce 
      END IF 
 

      IF ( l_hice_in ) THEN ! sea-ice thickness field provided as input
        ! Create new ice 
        IF( hice_p(isi) < (hice_min-csmall) ) THEN 
          hice_p(isi) = hice_new
          tice_p(isi) = tf_salt
        END IF
      ELSE  ! sea-ice thickness field NOT provided as input
        hice_p(isi) = hice_new
        tice_p(isi) = t_seasfc(isi)  ! use sea surface temperature
      ENDIF 


      ! Take security measures 
      hice_p(isi) = MAX(MIN(hice_p(isi), hice_max), hice_min) 
      tice_p(isi) = MIN(tice_p(isi), tf_fresh)

      ! Set temperature of snow upper surface and snow thickness
      tsnow_p(isi) = tice_p(isi)
      hsnow_p(isi) = 0._wp

      ! Set variables at new time level
      tice_n(isi)  = tice_p(isi)     
      hice_n(isi)  = hice_p(isi)     
      tsnow_n(isi) = tsnow_p(isi)    
      hsnow_n(isi) = hsnow_p(isi)    

    END DO GridBoxesWithSeaIce  

    ! Call model abort if errors are encountered
    IF( lcallabort ) THEN 
      ! Send an error message 
      WRITE(nameerr,*) "MODULE mo_seaice_nwp, SUBROUTINE seaice_init_nwp" 
      WRITE(texterr,*) "Sea-ice fraction ", frsi(isi),                        & 
                    &  " is less than a minimum threshold value ", frsi_min,  & 
                    &  " Call model abort."
      CALL message(TRIM(nameerr), TRIM(texterr))

      ! Call model abort
      WRITE(nameerr,*) "mo_seaice_nwp:seaice_init_nwp" 
      WRITE(texterr,*) "error in sea-ice fraction"
      CALL finish(TRIM(nameerr), TRIM(texterr))
    END IF 

    !-----------------------------------------------------------------------------------------------
    !  End calculations
    !===============================================================================================

  END SUBROUTINE seaice_init_nwp

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

  !>
  !! Prognostic variables of the sea-ice scheme are advanced one time step.
  !!
  !! Ordinary differential equations (in time) for the ice surface temperature and 
  !! the ice thickness are solved using an explicit Euler scheme for time advance.
  !! The shape factor for the temperature profile within the ice and the derivative of the 
  !! temperature profile shape function at the underside of the ice are functions of the ice 
  !! thickness (see Mironov et al. 2012, for details).  
  !! Optionally, constant values of the shape factor and of the shape-function derivative
  !! corresponding to the linear temperature profile within the ice 
  !! (cf. the sea ice schemes of GME and COSMO) can be used 
  !! (as there is no logical switch to activate this option, changes in the code should be made).
  !! In the regime of ice growth or melting from below, the solution may become numerically 
  !! unstable if the ice is thin. In such a case, a quasi-steady heat transfer through the ice 
  !! is assumed and the differential equation for the ice surface temperature 
  !! is reduced to an algebraic relation. 
  !! In the current configuration of the sea ice scheme, snow over ice is not treated explicitly. 
  !! The effect of snow is accounted for implicitely through an empirical temperature dependence 
  !! of the ice surface albedo with respect to solar radiation. 
  !! For the "sea water" type grid boxes, the snow thickness is set to zero and
  !! the snow surface temperature is set equal to the ice surface temperature.
  !! As the coupling between the sea ice and the sea water beneath is not considered, 
  !! the heat flux from water to ice is neglected.
  !! Prognostic ice thickness is limited by a maximum value of 3 m and a minimum value of 0.05 m. 
  !! No ice is created during the forecast period. 
  !! If the ice melts out during the forecast (i.e. the ice becomes thinner than 0.05 m), 
  !! the ice tickness is set to zero and 
  !! the ice surface temperatures is set to the fresh-water freezing point.
  !! The procedure arguments are arrays (vectors) of the sea-ice scheme prognostic variables
  !! and of the components of the heat balance at the ice upper surface, 
  !! i.e. the fluxes of sensible and latent heat, the net flux of long-wave radiation,
  !! and the net flux of solar radiation (with due regard for the ice surface albedo).
  !! Fluxes are positive when directed downward.
  !! The vector length is equal to the number of model grid boxes 
  !! (within a given block) where sea ice is present  
  !! (i.e. where the sea-ice fraction exceeds its minimum threshold value).
  !! The time tendecies of the ice thickness, the ice surface temperature, 
  !! the snow thickness and the snow surface temperature are also computed. 
  !! These are optional arguments of the procedure. 
  !!
  !!
  !! @par Revision History
  !! Initial release by Dmitrii Mironov, DWD (2012-07-24)
  !!

  SUBROUTINE seaice_timestep_nwp (                                      &
                              &  dtime,                                 &
                              &  nproma, nsigb,                         &
                              &  qsen, qlat, qlwrnet, qsolnet,          &
                              &  tice_p, hice_p, tsnow_p, hsnow_p,      &
                              &  tice_n, hice_n, tsnow_n, hsnow_n,      &
                              &  opt_dticedt, opt_dhicedt, opt_dtsnowdt,&
                              &  opt_dhsnowdt                           )

    IMPLICIT NONE

    ! Procedure arguments 

    REAL(wp), INTENT(IN) ::        & 
                         &  dtime  !< model time step [s]

    INTEGER, INTENT(IN) ::        &
                        &  nproma !< Array (vector) dimension

    INTEGER, INTENT(IN) ::        &
                        &  nsigb  !< number of grid boxes within a block 
                                  !< where the sea ice is present 

    REAL(wp), DIMENSION(:), INTENT(IN) ::           & 
                                           &  qsen    , &  !< sensible heat flux at the surface [W/m^2]
                                           &  qlat    , &  !< latent heat flux at the surface [W/m^2]
                                           &  qlwrnet , &  !< net long-wave radiation flux at the surface [W/m^2] 
                                           &  qsolnet      !< net solar radiation flux at the surface [W/m^2] 

    REAL(wp), DIMENSION(:), INTENT(IN) ::           &
                                           &  tice_p  , &  !< temperature of ice upper surface at previous time level [K] 
                                           &  hice_p  , &  !< ice thickness at previous time level [m] 
                                           &  tsnow_p , &  !< temperature of snow upper surface at previous time level [K] 
                                           &  hsnow_p      !< snow thickness at previous time level [m] 

    REAL(wp), DIMENSION(:), INTENT(OUT) ::           &
                                            &  tice_n  , &  !< temperature of ice upper surface at new time level [K] 
                                            &  hice_n  , &  !< ice thickness at new time level [m] 
                                            &  tsnow_n , &  !< temperature of snow upper surface at new time level [K] 
                                            &  hsnow_n      !< snow thickness at new time level [m] 

    REAL(wp), DIMENSION(:), INTENT(OUT), OPTIONAL ::         &
                                            &  opt_dticedt , &  !< time tendency of ice surface temperature [K/s] 
                                            &  opt_dhicedt , &  !< time tendency of ice thickness [m/s] 
                                            &  opt_dtsnowdt, &  !< time tendency of snow surface temperature [K/s] 
                                            &  opt_dhsnowdt     !< time tendency of snow thickness [m/s] 

    ! Derived parameters 
    ! (combinations of physical constants encountered several times in the code)

    REAL (wp), PARAMETER ::                                   &  
                         &  r_rhoici     = 1._wp/(rhoi*ci)  , &  
                         &  ki_o_rhoici  = ki*r_rhoici      , & 
                         &  r_rhoialf    = 1._wp/(rhoi*alf) , &
                         &  ki_o_rhoialf = ki*r_rhoialf     , & 
                         &  ci_o_alf     = ci/alf        

    ! Local variables 
    REAL(wp), DIMENSION(nproma) :: &
                                &  dticedt , &  !< time tendency of ice surface temperature [K/s] 
                                &  dhicedt , &  !< time tendency of ice thickness [m/s] 
                                &  dtsnowdt, &  !< time tendency of snow surface temperature [K/s] 
                                &  dhsnowdt     !< time tendency of snow thickness [m/s] 

    INTEGER ::      &
            &  isi  !< DO loop index

    REAL (wp) ::                &
              &  qsoliw       , &  !< solar radiation flux at the ice-water interface (positive downward) [W/m^2]
              &  qatm         , &  !< "total atmospheric heat flux" for the ice slab  (positive downward) [W/m^2]
                                   !< (sum of the sensible heat flux, latent heat flux and the net flux 
                                   !< of long-waver radiation at the ice upper surface, 
                                   !< and the difference of solar radiation fluxes 
                                   !< at the upper and lower surfaces of the ice slab)
              &  qwat         , &  !< heat flux at the ice-water interface (positive downward) [W/m^2]
              &  csice        , &  !< shape factor for the temperature profile within the ice [-]
              &  phiipr0      , &  !< derivative (at zeta=0) of the temperature profile shape function [-]
              &  hice_thrshld , &  !< threshold value of the ice tickness to switch between a quasi-equilibrium
                                   !< model of heat transfer through the ice and a complete model [m]
              &  rti          , &  !< dimensionless parameter [-]
              &  r_dtime      , &  !< reciprocal of the time step [1/s]
              &  strg_1       , &  !< temporary storage variable 
              &  strg_2            !< temporary storage variable 

    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------


    ! Reciprocal of the time step 
    r_dtime = 1._wp/dtime

    ! Loop over grid boxes where sea ice is present
    GridBoxesWithSeaIce: DO isi=1, nsigb

      ! Compute solar radiation flux at the ice-water interface (positive downward) 
      qsoliw = qsolnet(isi)*(                                                                 &
             & opticpar_seaice_opaque%frac_optic(1)                                           &
             & *EXP(-MIN(opticpar_seaice_opaque%extincoef_optic(1)*hice_p(isi), cmaxearg)) +  &
             & opticpar_seaice_opaque%frac_optic(2)                                           &
             & *EXP(-MIN(opticpar_seaice_opaque%extincoef_optic(2)*hice_p(isi), cmaxearg)) )

      ! Compute total atmospheric heat flux for the ice slab  (positive downward) 
      qatm = qsen(isi) + qlat(isi) + qlwrnet(isi) + qsolnet(isi) - qsoliw 

!_dev>
      ! Provision is made to account for the heat flux from water to ice 
      ! Currently the heat flux from water to ice is set to zero 
      qwat = 0._wp
!_dev<

      ! Compute temperature profile shape factor and temporary storage variables
      ! (to recover linear temperature profile, set csice=csi_lin)
      csice = csi_lin - csidp_nlin_d*hice_p(isi)/hice_max
      rti = ci_o_alf*(tice_p(isi)-tf_salt)
      strg_1 = rti*(1._wp-csice)
      strg_2 = 1._wp + strg_1 

      FreezingMeltingRegime: IF( tice_p(isi)>=(tf_fresh-csmall) .AND. qatm>0._wp ) THEN 

        ! Melting from above

        ! Set the ice surface temperature equal to the fresh-water freezing point
        tice_n(isi) = tf_fresh

        ! Compute the rate of ice melting (note the sign of heat fluxes)
        dhicedt(isi) = -(qatm-qwat)*r_rhoialf/strg_2 
        ! Update the ice thickness
        hice_n(isi) = hice_p(isi) + dtime*dhicedt(isi)

      ELSE FreezingMeltingRegime 

        ! Freezing or melting from below

        ! Derivative (at zeta=0) of the temperature profile shape function 
        ! (to recover linear temperature profile, set phiipr0=1)
        phiipr0 = 1._wp - hice_p(isi)/hice_max

        ! Compute threshold value of the ice thickness  
        ! Note that an expression in parentheses should be multiplied with 
        ! MAX(1._wp, ABS(2._wp*csice*rti)) (cf. the code of the lake parameterization scheme FLake).
        ! However, |2*csice*(ci/alf)*(tice-tf_salt)| < 1 
        ! at all conceivable values of the ice surface temperature.
        hice_thrshld = SQRT(phiipr0*ki_o_rhoici*dtime/csice)

        IF( hice_p(isi)<hice_thrshld) THEN 

          ! Use a quasi-equilibrium model of heat transfer through the ice 

          ! Compute the time-rate-of-change of the ice thickness (note the sign of heat fluxes)
          dhicedt(isi) = -(qatm-qwat*strg_2)*r_rhoialf
          ! Update the ice thickness
          hice_n(isi) = hice_p(isi) + dtime*dhicedt(isi)

          ! Compute the (updated) ice surface temperature (note the sign of heat fluxes)
          tice_n(isi) = tf_salt + (qatm-qwat*strg_1)*hice_n(isi)/(phiipr0*ki)

        ELSE

          ! Use a complete model of heat transfer through the ice 

          ! Use dhsnowdt as a temporary storage
          dhsnowdt(isi) = phiipr0*(tice_p(isi)-tf_salt)/hice_p(isi)

          ! Compute the time-rate-of-change of the ice surface temperature 
          ! (note the sign of heat fluxes)
          dticedt(isi) = ( (qatm-strg_1*qwat)*r_rhoici - ki_o_rhoici*dhsnowdt(isi)*strg_2 )  & 
                     & /(csice*hice_p(isi))  
          ! Update the ice surface temperature
          tice_n(isi) = tice_p(isi) + dtime*dticedt(isi)

          ! Compute the time-rate-of-change of the ice thickness (note the sign of heat fluxes)
          dhicedt(isi) = -ki_o_rhoialf*dhsnowdt(isi) + qwat*r_rhoialf
          ! Update the ice thickness
          hice_n(isi) = hice_p(isi) + dtime*dhicedt(isi)

        END IF 

      END IF FreezingMeltingRegime 

      ! Remove too thin ice or impose security constraints
      IF( hice_n(isi)<hice_min ) THEN
        ! Remove too thin ice
        hice_n(isi) = 0._wp
        ! Set the ice surface temperature equal to the fresh water freezing point
        tice_n(isi) = tf_fresh
      ELSE 
        ! Limit the ice thickness from above (security)
        hice_n(isi) = MIN(hice_n(isi), hice_max)
        ! Limit the ice surface temperature from above (security)
        tice_n(isi) = MIN(tice_n(isi), tf_fresh)
      END IF 

      ! Currently, snow over sea ice is not treated explicitly
      ! Set the snow thickness to zero 
      hsnow_n(isi) = 0._wp
      ! Set the snow surface temperature equal to the ice surface temperature
      tsnow_n(isi) = tice_n(isi)

      ! Compute tendencies (for eventual use outside the sea-ice scheme program units)
      dticedt(isi)  = (tice_n(isi)-tice_p(isi))*r_dtime
      dhicedt(isi)  = (hice_n(isi)-hice_p(isi))*r_dtime
      dtsnowdt(isi) = (tsnow_n(isi)-tsnow_p(isi))*r_dtime
      dhsnowdt(isi) = (hsnow_n(isi)-hsnow_p(isi))*r_dtime

    END DO GridBoxesWithSeaIce


    IF (PRESENT(opt_dticedt)) THEN
      opt_dticedt(:) = dticedt(:)
    ENDIF
    IF (PRESENT(opt_dhicedt)) THEN
      opt_dhicedt(:) = dhicedt(:)
    ENDIF
    IF (PRESENT(opt_dtsnowdt)) THEN
      opt_dtsnowdt(:) = dtsnowdt(:)
    ENDIF
    IF (PRESENT(opt_dhsnowdt)) THEN
      opt_dhsnowdt(:) = dhsnowdt(:)
    ENDIF
 
    !-----------------------------------------------------------------------------------------------
    !  End calculations
    !===============================================================================================

  END SUBROUTINE seaice_timestep_nwp

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

END MODULE mo_seaice_nwp

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890
