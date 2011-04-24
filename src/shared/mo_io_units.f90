!>
!!                 Sets the standard I/O units.
!!
!!
!! @par Revision History
!!  Initial version  by Luis Kornblueh (2004).
!!  Modification by Thomas Heinze (2005-08-29):
!!   - no more comments in ifdef declarations
!!   - included VERSION CONTROL
!!  Modification by Thomas Heinze (2006-02-21):
!!  - renamed m_modules to mo_modules
!!  Modification by Hui Wan (2007-02):
!!  - parameters <i>input1</i> and <i>input3</i> deleted.
!!    (They specified the input namelist files of the grid generator,
!!    but were not used at all in the shallow water model.)
!!  - 'ini_file' renamed as 'namelist_file' and set to 'NAMELIST_ICON'
!!  - 'root_dir' defined, which is used for setting the default path to the
!!     files providing grid/patch/topography data.
!!  Modification by Thomas Heinze, DWD, (2007-03-05):
!!  - reintroduced nlog
!!  Modification by Thomas Heinze, DWD, (2007-05-10):
!!  - moved unit1, unit2 and output1 from mo_io_graph here to avoid
!!    cyclic dependencies in Makefile
!!  Modifications by Luis Kornblueh, MPI-M, 2008-04-21:
!!  - going back to original design
!!  - added some predefined units
!!  - added function to find next free unit
!!
!! @par Copyright
!! 2002-2007 by DWD and MPI-M
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
!!
MODULE mo_io_units
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: version = '$Id$'

  PUBLIC

! This paramter is taken from /usr/include/stdio.h (ANSI C standard). If problems
! with filename length appear, check the before mentioned file.

  INTEGER, PARAMETER :: filename_max = 1024

! Standard I/O-units

#ifdef hpux
  INTEGER, PARAMETER :: nerr  = 7     ! error output
#else
  INTEGER, PARAMETER :: nerr  = 0     ! error output
#endif
  INTEGER, PARAMETER :: nlog  = 1     ! standard log file unit
  INTEGER, PARAMETER :: nnml  = 2     ! standard namelist file unit
  INTEGER, PARAMETER :: nstat = 3     ! standard statistics file unit
  INTEGER, PARAMETER :: ngmt  = 4     ! standard GMT output file unit
  INTEGER, PARAMETER :: nin   = 5     ! standard input
  INTEGER, PARAMETER :: nout  = 6     ! standard output

  INTEGER, PARAMETER :: nan   = -1    ! unit given back, when nothing
                                      ! in the allowed range is available

  INTEGER :: nnml_output ! unit of the ASCII output that contains the
                         ! namelist variables and their actual values.
!-------------------------------------------------------------------------

CONTAINS

  FUNCTION find_next_free_unit(istart,istop) RESULT(iunit)
    INTEGER :: iunit
    INTEGER, INTENT(in) :: istart, istop
    INTEGER :: kstart, kstop
    LOGICAL :: found, opened
    INTEGER :: i
    CHARACTER(len=32) :: info

    found = .FALSE.
    kstart = istart
    kstop  = istop
    IF (kstart < 10) kstart = 10
    IF (kstop <= kstart) kstop = kstart+10
    DO i = kstart, kstop
      INQUIRE(unit=i, opened=opened)
      IF (.NOT. opened) THEN
        iunit = i
        found = .TRUE.
        EXIT
      END IF
    END DO

    IF (.NOT. found) THEN
      WRITE(info,'(a,i0,a,i0,a)') &
           'No unit in range <', kstart, ':', kstop, '> free.'
      WRITE(nerr,'(a,a,a)') 'find_next_free_unit', ': ',TRIM(info)
      iunit = nan
    END IF

  END FUNCTION find_next_free_unit

END MODULE mo_io_units














