!>
!! Switches and parameters of the ECHAM-HAM aerosol module.
!!
!! @author <name, affiliation>
!! @author <name, affiliation>
!!
!! @par Revision History
!! <Description of activity> by <name, affiliation> (<YYYY-MM-DD>)
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
MODULE mo_ham_aerosol_params

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: ncdnc, nicnc
  PUBLIC  :: iaiti, iacci, icoai
  PUBLIC  :: t_vmem3d
  PUBLIC  :: rwet

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

  CHARACTER(len=*), PARAMETER :: version = '$Id$'
  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_ham_aerosol_params'

END MODULE mo_ham_aerosol_params

