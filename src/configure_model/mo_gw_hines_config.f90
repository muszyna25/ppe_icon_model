!>
!! <Short description of module for listings and indices>
!!
!! <Describe the concepts of the procedures and algorithms used in the module.>
!! <Details of procedures are documented below with their definitions.>
!! <Include any applicable external references inline as module::procedure,>
!! <external_procedure(), or by using @see.>
!! <Don't forget references to literature.>
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_gw_hines_config

  USE mo_kind,           ONLY: wp
  USE mo_impl_constants, ONLY: max_dom

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_gw_hines_config, gw_hines_config

  !--------------------------------------------------------------------------
  ! Basic configuration setup for Hines gravity wave drag scheme
  !--------------------------------------------------------------------------
  TYPE t_gw_hines_config

  LOGICAL  :: lheatcal      !< true : compute momentum flux dep., heating and diffusion coefficient
                            !< false: compute only momentum flux deposition

  INTEGER  :: emiss_lev     !< number of levels above the ground at which gw are emitted
  REAL(wp) :: rmscon        !< [m/s] root mean square gravity wave wind at emission level
  REAL(wp) :: kstar         !< [1/m] typical gravity wave horizontal wavenumber
  REAL(wp) :: m_min         !< [1/m] minimum bound in vertical wavenumber

!!$  LOGICAL  :: lfront        !< true: compute gw sources emerging from fronts and background
!!$                            !< (Charron and Manzini, 2002)
!!$  REAL(wp) :: rms_front     !< [m/s] rms frontal gw wind at source level
!!$  REAL(wp) :: front_thres   !< [(K/m)^2/hour] minimum value of the frontogenesis function,
!!$                            !< for whichgravity waves are emitted from fronts
!!$
!!$  LOGICAL  :: lozpr         !< true: for background enhancement associated with precipitation
!!$                            !< (Manzini et al., 1997)
!!$  REAL(wp) :: pcrit         !< [mm/d] critical precipitation value, above which 
!!$                            !< gravity wave rms wind enhancement is applied
!!$  REAL(wp) :: pcons         !< [] adimensional factor for background enhancement 
!!$                            !< associated with precipitation
!!$
  LOGICAL  :: lrmscon_lat   !< true:  use latitude dependent rmscon
                            !< - |latitude| >= lat_rmscon:
                            !<      use rmscon
                            !< - |latitude| <= lat_rmscon_eq:
                            !<      use rmscon_eq
                            !< - lat_rmscon_eq < |latitude| < lat_rmscon: 
                            !<      use linear interpolation between rmscon_eq and rmscon
                            !< false: use rmscon for all latitudes
                            !< attention: may be overwritten if lfront or lozpr is true
  REAL(wp) :: lat_rmscon_eq !< [degN] rmscon_eq is used equatorward of this latitude
  REAL(wp) :: lat_rmscon    !< [degN] rmscon is used poleward of this latitude
  REAL(wp) :: rmscon_eq     !< [m/s]  rms constant used equatorward of lat_rmscon_eq

  END TYPE t_gw_hines_config

  !>
  !!
  TYPE(t_gw_hines_config) :: gw_hines_config(max_dom)

END MODULE mo_gw_hines_config
