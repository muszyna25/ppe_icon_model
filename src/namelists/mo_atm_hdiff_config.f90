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
MODULE mo_atm_hdiff_config

  USE mo_kind, ONLY: wp

  IMPLICIT NONE

  PUBLIC ! all "get" functions; "fill/kill" subroutines

  PRIVATE :: t_hdiff_config, atm_hdiff_config

  CHARACTER(len=*),PARAMETER,PRIVATE :: version = '$Id$'

  !>
  !!--------------------------------------------------------------------------
  !! Type definition 
  !!--------------------------------------------------------------------------
  TYPE :: t_hdiff_config

    LOGICAL :: lhdiff_temp  ! if .TRUE., apply horizontal diffusion to thermodynamic variable.
    LOGICAL :: lhdiff_vn    ! if .TRUE., apply horizontal diffusion to momentum.

    INTEGER :: hdiff_order  ! order of horizontal diffusion
                            ! 2: 2nd order linear diffusion on all vertical levels 
                            ! 3: Smagorinsky diffusion for hexagonal model
                            ! 4: 4th order linear diffusion on all vertical levels 
                            ! 5: Smagorinsky diffusion for triangular model
                            ! 24 or 42: 2nd order linear diffusion for upper levels,
                            !           4th order for lower levels
  
    REAL(wp) ::               &
      & hdiff_efdt_ratio,     &! ratio of e-folding time to (2*)time step
      & hdiff_min_efdt_ratio, &! minimum value of hdiff_efdt_ratio (for upper sponge layer)
      & hdiff_tv_ratio,       &! the ratio of diffusion coefficient: temp:mom
      & hdiff_smag_fac,       &! scaling factor for Smagorinsky diffusion
      & hdiff_multfac          ! multiplication factor of normalized diffusion coefficient
                               ! for nested domains
    REAL(wp) ::              &
      & k6, k4, k2           ! numerical diffusion coefficients
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
  END TYPE t_hdiff_config

  !>
  !!
  TYPE(t_hdiff_config),ALLOCATABLE :: atm_hdiff_config(:) !< shape: (ndom)

CONTAINS
  !----------------------------------------------------------------
  ! Subroutines for filling and killing the configuration state
  !----------------------------------------------------------------
  !>
  !!
  SUBROUTINE fill_atm_hdiff_config(ndom)

    INTEGER,INTENT(IN) :: ndom

    ALLOCATE(atm_hdiff_config(ndom))

  END SUBROUTINE fill_atm_hdiff_config
  !----------------------------------------------------------------
  !>
  !!
  SUBROUTINE kill_atm_hdiff_config
    DEALLOCATE(atm_hdiff_config)
  END SUBROUTINE kill_atm_hdiff_config

  !----------------------------------------------------------------
  ! "Get" functions
  !----------------------------------------------------------------
  !>
  !!
  LOGICAL FUNCTION get_lhdiff_temp(jg)
    INTEGER,INTENT(IN) :: jg
    get_lhdiff_temp = atm_hdiff_config(jg)%lhdiff_temp
  END FUNCTION get_lhdiff_temp
  !----------------------------------------------------------------
  !>
  !!
  LOGICAL FUNCTION get_lhdiff_vn(jg)
    INTEGER,INTENT(IN) :: jg
    get_lhdiff_vn = atm_hdiff_config(jg)%lhdiff_vn
  END FUNCTION get_lhdiff_vn
  !----------------------------------------------------------------
  !>
  !!
  REAL(wp) FUNCTION get_k2(jg)
    INTEGER,INTENT(IN) :: jg
    get_k2 = atm_hdiff_config(jg)%k2
  END FUNCTION get_k2
  !----------------------------------------------------------------
  !>
  !!
  REAL(wp) FUNCTION get_k4(jg)
    INTEGER,INTENT(IN) :: jg
    get_k4 = atm_hdiff_config(jg)%k4
  END FUNCTION get_k4
  !----------------------------------------------------------------
  !>
  !!
  INTEGER FUNCTION get_k2s(jg)
    INTEGER,INTENT(IN) :: jg
    get_k2s = atm_hdiff_config(jg)%k2s
  END FUNCTION get_k2s
  !----------------------------------------------------------------
  !>
  !!
  INTEGER FUNCTION get_k2e(jg)
    INTEGER,INTENT(IN) :: jg
    get_k2e = atm_hdiff_config(jg)%k2e
  END FUNCTION get_k2e
  !----------------------------------------------------------------
  !>
  !!
  INTEGER FUNCTION get_k4s(jg)
    INTEGER,INTENT(IN) :: jg
    get_k4s = atm_hdiff_config(jg)%k4s
  END FUNCTION get_k4s
  !----------------------------------------------------------------
  !>
  !!
  INTEGER FUNCTION get_k4e(jg)
    INTEGER,INTENT(IN) :: jg
    get_k4e = atm_hdiff_config(jg)%k4e
  END FUNCTION get_k4e
  !----------------------------------------------------------------

END MODULE mo_atm_hdiff_config

