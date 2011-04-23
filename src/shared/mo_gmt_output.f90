!>
!! Subroutien gmt_out writes out an array into an ASCII file
!! which can be visualized using the "pscontour" command of GMT.
!! This is mainly ment to serve as a debugging tool.
!!
!! @author Hui Wan
!!
!! @par Revision History
!! Initial version by Hui Wan (2010-08)
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
MODULE mo_gmt_output

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: ngmt, nnml, nnml_output
  USE mo_mpi,                ONLY: p_pe, p_io
  USE mo_namelist,           ONLY: position_nml, POSITIONED
  USE mo_exception,          ONLY: finish
  USE mo_math_constants,     ONLY: rad2deg
  USE mo_atmo_control,       ONLY: p_patch

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: lgmt_output, igmt_output_lev    !< namelist variables
  PUBLIC  :: setup_gmt_output                !< subroutine
  PUBLIC  :: gmt_out                         !< subroutine (interface)

  INTERFACE gmt_out
    MODULE PROCEDURE gmt_out_real  !< for variables of type REAL
    MODULE PROCEDURE gmt_out_int   !< for variables of type INTEGER
  END INTERFACE gmt_out

  LOGICAL :: lgmt_output           !< .TRUE. for writting GMT files
  INTEGER :: igmt_output_lev       !< the vertical level to write out

  NAMELIST/gmt_ctl/ lgmt_output, igmt_output_lev

CONTAINS
  !>
  !! Read namelist "gmt_ctl" to set up the GMT output
  !!
  !! @par Revision History
  !! Initial version by Hui Wan, MPI-M (2010-08-04)
  !! Modification by Daniel Reinert, DWD (2011-02-08)
  !! - adapted to vertical nesting. Default-nlev is set to  
  !!   the minimum of all patch-nlevs.  
  !!
  SUBROUTINE setup_gmt_output(nlev)

    INTEGER, INTENT(IN) :: nlev  !< number of full levels
    INTEGER :: ist

    ! Set default values

    lgmt_output     = .FALSE.
    igmt_output_lev = nlev

    ! Read namelist (every CPU does this)

    CALL position_nml ('gmt_ctl', STATUS=ist)
    SELECT CASE (ist)
    CASE (POSITIONED)
      READ (nnml, gmt_ctl)
    END SELECT

    ! Sanity check, then write out the namelist

    IF (igmt_output_lev < 1 .OR. igmt_output_lev > nlev ) THEN
       CALL finish('setup_gmt_output','invalid value for igmt_output_lev')
    END IF
    IF(p_pe == p_io) WRITE(nnml_output,nml=gmt_ctl)

  END SUBROUTINE setup_gmt_output
  !-------------
  !>
  !! @brief Write out a floating point array
  !!
  !! @par Revision History
  !! Initial version by Hui Wan, MPI-M (2010-08-03)
  !!
  SUBROUTINE gmt_out_real( jstep,jg,jb,jce,nbdim, pos,varname,var )

    INTEGER,INTENT(IN)          :: jstep, jg, jb, jce, nbdim
    CHARACTER(len=*),INTENT(IN) :: pos
    CHARACTER(len=*),INTENT(IN) :: varname
    REAL(wp),INTENT(IN)         :: var(nbdim)

    CHARACTER(len=100) :: strtmp
    INTEGER :: jc

    WRITE(strtmp,'(a,i3.3,a,i3.3,5a,i3.3,a,i5.5,a)')             &
         & 'step',jstep,'_lev',igmt_output_lev,'_',              &
         & TRIM(pos),'_',TRIM(varname),'_PE',p_pe,'_b',jb,'.gmt'

    OPEN(ngmt,FILE=TRIM(strtmp),FORM='formatted',STATUS='unknown')
    DO jc = 1,jce
       WRITE(ngmt,'(3f20.10)') p_patch(jg)%cells%center(jc,jb)%lon*rad2deg, &
                               p_patch(jg)%cells%center(jc,jb)%lat*rad2deg, &
                               var(jc)
    ENDDO
    CLOSE(ngmt)

  END SUBROUTINE gmt_out_real
  !-------------
  !>
  !! @brief Write out an integer array
  !!
  !! @par Revision History
  !! Initial version by Hui Wan, MPI-M (2010-08-03)
  !!
  SUBROUTINE gmt_out_int( jstep,jg,jb,jce,nbdim, pos,varname,var )

    INTEGER,INTENT(IN)          :: jstep, jg, jb, jce, nbdim
    CHARACTER(len=*),INTENT(IN) :: pos
    CHARACTER(len=*),INTENT(IN) :: varname
    INTEGER,INTENT(IN)          :: var(nbdim)

    CHARACTER(len=100) :: strtmp
    INTEGER :: jc

    WRITE(strtmp,'(a,i3.3,a,i3.3,5a,i3.3,a,i5.5,a)')             &
         & 'step',jstep,'_lev',igmt_output_lev,'_',              &
         & TRIM(pos),'_',TRIM(varname),'_PE',p_pe,'_b',jb,'.gmt'

    OPEN(ngmt,FILE=TRIM(strtmp),FORM='formatted',STATUS='unknown')
    DO jc = 1,jce
       WRITE(ngmt,'(2f20.10,i15)') p_patch(jg)%cells%center(jc,jb)%lon*rad2deg, &
                                   p_patch(jg)%cells%center(jc,jb)%lat*rad2deg, &
                                   var(jc)
    ENDDO
    CLOSE(ngmt)

  END SUBROUTINE gmt_out_int
  !-------------

END MODULE mo_gmt_output
