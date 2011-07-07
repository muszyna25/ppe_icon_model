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
MODULE mo_diffusion_config

  USE mo_kind,           ONLY: wp
  USE mo_impl_constants, ONLY: MAX_NTRACER, MAX_CHAR_LENGTH, max_dom

  IMPLICIT NONE
  PUBLIC
  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'


  !!--------------------------------------------------------------------------
  !! Basic configuration setup for diffusion
  !!--------------------------------------------------------------------------


  TYPE t_diffusion_config

  INTEGER :: hdiff_order  ! order of horizontal diffusion
                          ! 2: 2nd order linear diffusion on all vertical levels 
                          ! 3: Smagorinsky diffusion for hexagonal model
                          ! 4: 4th order linear diffusion on all vertical levels 
                          ! 5: Smagorinsky diffusion for triangular model
                          ! 24 or 42: 2nd order linear diffusion for upper levels,
                          !           4th order for lower levels

  REAL(wp) :: k2_pres_max  ! (relevant only when hdiff_order = 24 or 42)
                           ! pressure (in Pa) specified by the user
                           ! to determine the lowest vertical level 
                           ! to which 2nd order linear diffusion is applied.
                           ! For the levels with pressure > k2_pres_max, 
                           ! 4th order linear diffusion is applied. 

  INTEGER  :: k2_klev_max  ! (relevant only when hdiff_order = 24 or 42)
                           ! vertical level index specified by the user
                           ! to determine the lowest vertical level 
                           ! to which 2nd order linear diffusion is applied.
                           ! For the levels with k > k2_klev_max, 
                           ! 4th order linear diffusion is applied. 

  REAL(wp) ::           &
    & hdiff_efdt_ratio, &! ratio of e-folding time to (2*)time step
    & hdiff_min_efdt_ratio, &! minimum value of hdiff_efdt_ratio (for upper sponge layer)
    & hdiff_tv_ratio,   &! the ratio of diffusion coefficient: temp:mom
    & hdiff_smag_fac,   &! scaling factor for Smagorinsky diffusion
    & hdiff_multfac      ! multiplication factor of normalized diffusion coefficient
                         ! for nested domains

  REAL(wp), ALLOCATABLE, DIMENSION(:) :: &
    & k6, k4, k2       ! numerical diffusion coefficients
                       ! Values for these parameters are not directly
                       ! specified by the user, but derived from the ratio 
                       ! between the e-folding time and the model time step
                       ! (hdiff_efdt_ratio above), and the horizontal 
                       ! resolution of the model

  INTEGER k2s, k2e, k4s, k4e  ! indices defining to which vertical levels
                              ! 2nd and 4th linear diffusion are applied.
                              ! The values are not specified by the user via namelist,
                              ! but determined from k2_klev_max, k2_pres_max
                              ! and the configuration of the vertical coordinate

  LOGICAL ::          &
    & lhdiff_temp,    &! if .TRUE., apply horizontal diffusion to temp.
    & lhdiff_vn        ! if .TRUE., apply horizontal diffusion to momentum.


  END TYPE t_diffusion_config
  !>
  !!
  TYPE(t_diffusion_config) :: diffusion_config(max_dom)

END MODULE mo_diffusion_config
