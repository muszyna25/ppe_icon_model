!>
!! Switches and parameters of the ECHAM-HAM aerosol module.
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
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
MODULE mo_ham_aerosol_params

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE



  INTEGER :: ncdnc = 0  !< CDNC activation is off
  INTEGER :: nicnc = 0  !< ICNC activation is off

  !------------------
  ! log-normal modes
  !------------------

  INTEGER, PARAMETER :: nmod=7  !< # of modes considered by the aerosol microphysics module
  INTEGER, PARAMETER :: inucs=1,  iaits=2,  iaccs=3,  icoas=4,  iaiti=5,  iacci=6,  icoai=7
  !                     nucl.   | aitk.   | acc.    | coar.   | aitk.   | acc.    | coar.   |
  !                     soluble | soluble | soluble | soluble | insol.  | insol.  | insol.  |

  TYPE t_vmem3d
    REAL(wp), POINTER  :: ptr(:,:,:)
  END TYPE t_vmem3d

  TYPE(t_vmem3d) :: rwet(nmod)    !< wet radius (?)

  !------------------

  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_ham_aerosol_params'


  PUBLIC  :: ncdnc, nicnc
  PUBLIC  :: iaiti, iacci, icoai
  PUBLIC  :: t_vmem3d
  PUBLIC  :: rwet

END MODULE mo_ham_aerosol_params

