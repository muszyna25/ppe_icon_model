!>
!! @brief Parameters for the boundary layer parameterization. 
!!
!! @author Hui Wan, MPI-M
!!
!! @par Revision History
!! First version by Hui Wan, MPI-M (2010-08-31)
!! Modification by Constantin Junk, MPI-M (2011-05-05)
!! - moved echam_vdiff_ctl and setup_vdiff to namelists/mo_echam_vdiff_nml
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
MODULE mo_echam_vdiff_params

  USE mo_kind,              ONLY: wp
  USE mo_physical_constants,ONLY: grav
  USE mo_impl_constants,    ONLY: cvdifts

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ckap, cchar, cb, cc                         !< parameters
  PUBLIC :: eps_shear, eps_corio, totte_min             !< parameters
  PUBLIC :: chneu, shn, smn, da1                        !< parameters
  PUBLIC :: cons5                                       !< parameters
  PUBLIC :: tpfac1, tpfac2, tpfac3, tpfac4     !< parameters
  PUBLIC :: itop                                        !< parameters
  
  !-------------------
  ! Module parameters
  !-------------------

  REAL(wp),PARAMETER :: ckap    = 0.4_wp       !< karman constant.
  REAL(wp),PARAMETER :: cchar   = 0.018_wp     !< charnock constant.
  REAL(wp),PARAMETER :: cb      = 5._wp        !< stability parameter near neutrality.
  REAL(wp),PARAMETER :: cc      = 5._wp        !< stability parameter for unstable cases.

  REAL(wp),PARAMETER :: lmix_max  = 150._wp    !< maximum mixing length in neutral and stable conditions
  REAL(wp),PARAMETER :: eps_shear = 1.e-5_wp   !< zepshr in sbr. vdiff of ECHAM6
  REAL(wp),PARAMETER :: eps_corio = 5.e-5_wp   !< zepcor in sbr. vdiff of ECHAM6
  REAL(wp),PARAMETER :: totte_min = 1.e-10_wp  !< minimum total turbulent energy 

  REAL(wp),PARAMETER :: chneu  = 0.3_wp

  REAL(wp),PARAMETER :: cons5  = 3._wp*cb*cc

  ! Parameters related to time step weighting in *rhs* of *vdiff* and *scv*

  REAL(wp),PARAMETER :: tpfac1  = cvdifts
  REAL(wp),PARAMETER :: tpfac2  = 1._wp / tpfac1
  REAL(wp),PARAMETER :: tpfac3  = 1._wp - tpfac2
  REAL(wp),PARAMETER :: tpfac4  = 1._wp + tpfac3

  REAL(wp),PARAMETER :: shn = 2.22_wp*0.22_wp*SQRT(2._wp)
  REAL(wp),PARAMETER :: smn = shn*1.24_wp*2.37_wp/3.69_wp
  REAL(wp),PARAMETER :: da1 = 1._wp/smn**3

  INTEGER ,PARAMETER :: itop = 1

END MODULE mo_echam_vdiff_params

