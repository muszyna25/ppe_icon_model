!>
!! @brief variables for marine biology.
!!
!! Definition of variables and allocation of memory
!!
!! @author S.Legutke (MPI-M)
!!
!! @par Revision History
!!
!! First version by S.Legutke    (MPI-M)    Oct 31, 2001
!!
!! @par Copyright
!! 2002-2013 by MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
MODULE mo_biomod

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC

  INTEGER :: n90depth,n1000depth,n2000depth
  INTEGER, DIMENSION (:), ALLOCATABLE  :: kbo   !< k-index of bottom layer (2d)
  REAL(wp), DIMENSION (:), ALLOCATABLE :: bolay !<  height of bottom cell

  REAL(wp), DIMENSION (:), ALLOCATABLE :: expoor
  REAL(wp), DIMENSION (:), ALLOCATABLE :: expoca
  REAL(wp), DIMENSION (:), ALLOCATABLE :: exposi

  REAL(wp), DIMENSION (:), ALLOCATABLE :: strahl

  REAL(wp), DIMENSION (:), ALLOCATABLE :: alar1max
  REAL(wp), DIMENSION (:), ALLOCATABLE :: TSFmax
  REAL(wp), DIMENSION (:), ALLOCATABLE :: TMFmax

 
  REAL(wp), DIMENSION(6):: dmsp

  REAL(wp) :: phytomi, grami, grazra, rrrcl
  REAL(wp) :: remido, dyphy, zinges, epsher, spemor, gammap, gammaz, ecan
  REAL(wp) :: ro2ut, rcar, rnit, rnoi, rnit23, rnit13, rcalc, ropal, rn2, p2gtc
  REAL(wp) :: bkphy, bkzoo, bkopal, bifr13, bifr14, plafr13, plafr14
  REAL(wp) :: drempoc, dremdoc, dremcalc, dremn2o
  REAL(wp) :: remido_cya, dremdoc_cya
  REAL(wp) :: thresh_aerob ! O2 threshold for aerob POC remineralization
  REAL(wp) :: thresh_sred ! O2 threshold for sulfate reduction
  REAL(wp) :: thresh_o2, prodn2o 

  REAL(wp) :: mc_fac, mc_eu_speed, mc_depth
  REAL(wp) :: n2_fixation
  REAL(wp) :: sulfate_reduction
  REAL(wp) :: denitrification


  REAL(wp) :: wopal ! daily sinking speed of opal (namelist parameter)
  REAL(wp) :: wcal  ! daily sinking speed of cal (namelist parameter)
  REAL(wp) :: wcya  ! daily boyancy speed of cyanos

  REAL(wp) :: dremopal, calmax, gutc
  REAL(wp) :: psedi, csedi, ssedi
  REAL(wp) :: perc_diron, riron, fesoly, relaxfe, wdust,bolaymin
  REAL(wp) :: pi_alpha, fPAR
  REAL(wp) :: nitdem, n2prod, ro2nitri
  REAL(wp) :: dremn3o, ro2bal

  REAL(wp) :: dustd1, dustd2, dustd3, dustsink

  REAL(wp) :: cycdec, pi_alpha_cya,cya_growth_max,Topt_cya,T1_cya,T2_cya ! (namelist parameter)
  REAL(wp) :: bkcya_P, bkcya_fe, bkcya_N, doccya_fac           ! (namelist parameter)
  REAL(wp) :: buoyancyspeed_cya                                          ! daily buyancy speed of cyano (namelist parameter)
! ---------


CONTAINS

  SUBROUTINE ALLOC_MEM_BIOMOD

    USE mo_control_bgc
    USE mo_param1_bgc


    ALLOCATE (expoor(bgc_nproma))
    expoor(:)=0._wp
    ALLOCATE (expoca(bgc_nproma))
    expoca(:)=0._wp
    ALLOCATE (exposi(bgc_nproma))
    exposi(:)=0._wp
    ALLOCATE (kbo(bgc_nproma))
    ALLOCATE (strahl(bgc_nproma))
    strahl=200._wp
    ALLOCATE (bolay(bgc_nproma))
    ALLOCATE (alar1max(bgc_nproma))
    ALLOCATE (TSFmax(bgc_nproma))
    ALLOCATE (TMFmax(bgc_nproma))

  END SUBROUTINE ALLOC_MEM_BIOMOD

END MODULE mo_biomod
