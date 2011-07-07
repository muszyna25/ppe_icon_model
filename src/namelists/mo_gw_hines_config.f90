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
MODULE mo_gw_hines_config

  USE mo_kind,           ONLY: wp
  USE mo_impl_constants, ONLY: MAX_NTRACER, MAX_CHAR_LENGTH, max_dom

  IMPLICIT NONE
  PUBLIC
  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'


  !!--------------------------------------------------------------------------
  !! Basic configuration setup for GW Hines
  !!--------------------------------------------------------------------------


  TYPE t_gw_hines_config

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

  END TYPE t_gw_hines_config

  !>
  !!
  TYPE(t_gw_hines_config) :: gw_hines_config(max_dom)

END MODULE mo_gw_hines_config
