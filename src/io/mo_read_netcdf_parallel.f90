!>
!!               This module provides wrappers for netcdf functions.
!!
!!               This module provides wrappers for netcdf functions
!! for reading a NetCDF file in a parallel run.
!! These wrappers have the same interface as the corresponding
!! NetCDF routines, but only one processor actually reads the file
!! and broadcasts the data to the others.
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Added p_comm_input_bcast by Rainer Johanni, Oct 2010
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
MODULE mo_read_netcdf_parallel
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!

USE mo_kind
USE mo_mpi
USE mo_mpi, ONLY: p_comm_input_bcast

IMPLICIT NONE

PRIVATE

INCLUDE 'netcdf.inc'

INTEGER, PARAMETER :: nf_read = NF_NOWRITE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

!modules interface-------------------------------------------
!subroutines
PUBLIC :: p_nf_open
PUBLIC :: p_nf_close
PUBLIC :: p_nf_inq_dimid
PUBLIC :: p_nf_inq_dimlen
PUBLIC :: p_nf_inq_varid
PUBLIC :: p_nf_get_att_text
PUBLIC :: p_nf_get_att_int
PUBLIC :: p_nf_get_att_double
PUBLIC :: p_nf_get_var_int
PUBLIC :: p_nf_get_var_double
PUBLIC :: p_nf_get_vara_double

! constants
PUBLIC :: nf_read

! make some names from netcdf.inc also global

PUBLIC :: nf_nowrite, nf_global, nf_noerr, nf_strerror

INTERFACE p_nf_get_att_int
   MODULE PROCEDURE p_nf_get_att_int_0
   MODULE PROCEDURE p_nf_get_att_int_1
END INTERFACE

INTERFACE p_nf_get_att_double
   MODULE PROCEDURE p_nf_get_att_double_single
   MODULE PROCEDURE p_nf_get_att_double_array
END INTERFACE

!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!
!

!>
!!               Wrapper for nf_open.
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Added p_comm_input_bcast by Rainer Johanni, Oct 2010
!!
INTEGER FUNCTION p_nf_open(path, omode, ncid)

!
   CHARACTER*(*), INTENT(IN) :: path
   INTEGER, INTENT(IN) :: omode
   INTEGER, INTENT(OUT) :: ncid

   INTEGER :: res


!-----------------------------------------------------------------------

   IF(p_pe == p_io) THEN
      res = nf_open(path, omode, ncid)
   ELSE
      ncid = -1 ! set it to an invalid value
   ENDIF

   CALL p_bcast(res, p_io, p_comm_input_bcast)
   p_nf_open = res

END FUNCTION p_nf_open

!-------------------------------------------------------------------------
!
!

!>
!!               Wrapper for nf_close.
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Added p_comm_input_bcast by Rainer Johanni, Oct 2010
!!
INTEGER FUNCTION p_nf_close(ncid)

!
   INTEGER, INTENT(IN) :: ncid

   INTEGER :: res


!-----------------------------------------------------------------------

   IF(p_pe == p_io) THEN
      res = nf_close(ncid)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_input_bcast)
   p_nf_close = res

END FUNCTION p_nf_close

!-------------------------------------------------------------------------
!
!

!>
!!               Wrapper for nf_inq_dimid.
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Added p_comm_input_bcast by Rainer Johanni, Oct 2010
!!
INTEGER FUNCTION p_nf_inq_dimid(ncid, name, dimid)

!
   INTEGER, INTENT(IN) :: ncid
   CHARACTER*(*), INTENT(IN) :: name
   INTEGER, INTENT(OUT) :: dimid

   INTEGER :: res


!-----------------------------------------------------------------------

   IF(p_pe == p_io) THEN
      res = nf_inq_dimid(ncid, name, dimid)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_input_bcast)
   p_nf_inq_dimid = res

   CALL p_bcast(dimid, p_io, p_comm_input_bcast)

END FUNCTION p_nf_inq_dimid

!-------------------------------------------------------------------------
!
!

!>
!!               Wrapper for nf_inq_dimlen.
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Added p_comm_input_bcast by Rainer Johanni, Oct 2010
!!
INTEGER FUNCTION p_nf_inq_dimlen(ncid, dimid, len)

!
   INTEGER, INTENT(IN) :: ncid, dimid
   INTEGER, INTENT(OUT) :: len

   INTEGER :: res


!-----------------------------------------------------------------------

   IF(p_pe == p_io) THEN
      res = nf_inq_dimlen(ncid, dimid, len)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_input_bcast)
   p_nf_inq_dimlen = res

   CALL p_bcast(len, p_io, p_comm_input_bcast)

END FUNCTION p_nf_inq_dimlen

!-------------------------------------------------------------------------
!
!

!>
!!               Wrapper for nf_inq_varid.
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Added p_comm_input_bcast by Rainer Johanni, Oct 2010
!!
INTEGER FUNCTION p_nf_inq_varid(ncid, name, varid)

!
   INTEGER, INTENT(IN) :: ncid
   CHARACTER*(*), INTENT(IN) :: name
   INTEGER, INTENT(OUT) :: varid

   INTEGER :: res


!-----------------------------------------------------------------------

   IF(p_pe == p_io) THEN
      res = nf_inq_varid(ncid, name, varid)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_input_bcast)
   p_nf_inq_varid = res

   CALL p_bcast(varid, p_io, p_comm_input_bcast)

END FUNCTION p_nf_inq_varid

!-------------------------------------------------------------------------
!
!

!>
!!               Wrapper for nf_get_att_text.
!!
!!
!! @par Revision History
!! Initial version by Luis Kornblueh, Jan 2011
!!
INTEGER FUNCTION p_nf_get_att_text(ncid, varid, name, tval)
  INTEGER,          INTENT(in)  :: ncid, varid
  CHARACTER(len=*), INTENT(in)  :: name
  CHARACTER(len=*), INTENT(out) :: tval
  
  INTEGER :: res
  
  !-----------------------------------------------------------------------

  IF (p_pe == p_io) THEN
    res = nf_get_att_text(ncid, varid, name, tval)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_input_bcast)
  p_nf_get_att_text = res
  
  CALL p_bcast(tval, p_io, p_comm_input_bcast)

END FUNCTION p_nf_get_att_text
!-----------------------------------------------------------

!-----------------------------------------------------------
!>
!!               Wrapper for nf_get_att_double
!!
!! @par Revision History
!! Initial version by Leonidas Linardakis, May 2012
!!
INTEGER FUNCTION p_nf_get_att_double_single(ncid, varid, name, dvalue)

!
   INTEGER, INTENT(IN) :: ncid, varid
   CHARACTER*(*), INTENT(IN) :: name
   REAL(dp), INTENT(OUT) :: dvalue

   INTEGER :: res

   IF(p_pe == p_io) THEN
      res = nf_get_att_double(ncid, varid, name, dvalue)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_input_bcast)
   p_nf_get_att_double_single = res

   CALL p_bcast(dvalue, p_io, p_comm_input_bcast)

END FUNCTION p_nf_get_att_double_single
!-----------------------------------------------------------------------

!-----------------------------------------------------------
!>
!!               Wrapper for nf_get_att_double
!!
!! @par Revision History
!! Initial version by Leonidas Linardakis, May 2012
!!
INTEGER FUNCTION p_nf_get_att_double_array(ncid, varid, name, dvalue)

!
   INTEGER, INTENT(IN) :: ncid, varid
   CHARACTER*(*), INTENT(IN) :: name
   REAL(dp), INTENT(OUT) :: dvalue(:)

   INTEGER :: res

   IF(p_pe == p_io) THEN
      res = nf_get_att_double(ncid, varid, name, dvalue)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_input_bcast)
   p_nf_get_att_double_array = res

   CALL p_bcast(dvalue, p_io, p_comm_input_bcast)

END FUNCTION p_nf_get_att_double_array
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!>
!!               Wrapper for nf_get_att_int.
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Added p_comm_input_bcast by Rainer Johanni, Oct 2010
!!
INTEGER FUNCTION p_nf_get_att_int_0(ncid, varid, name, ivals)

!
   INTEGER, INTENT(IN) :: ncid, varid
   CHARACTER*(*), INTENT(IN) :: name
   INTEGER, INTENT(OUT) :: ivals

   INTEGER :: res


!-----------------------------------------------------------------------

   IF(p_pe == p_io) THEN
      res = nf_get_att_int(ncid, varid, name, ivals)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_input_bcast)
   p_nf_get_att_int_0 = res

   CALL p_bcast(ivals, p_io, p_comm_input_bcast)

END FUNCTION p_nf_get_att_int_0

!-------------------------------------------------------------------------
!
!

!>
!!               Wrapper for nf_get_att_int.
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Added p_comm_input_bcast by Rainer Johanni, Oct 2010
!!
INTEGER FUNCTION p_nf_get_att_int_1(ncid, varid, name, ivals)

!
   INTEGER, INTENT(IN) :: ncid, varid
   CHARACTER*(*), INTENT(IN) :: name
   INTEGER, INTENT(OUT) :: ivals(:)

   INTEGER :: res, len


!-----------------------------------------------------------------------

   IF(p_pe == p_io) THEN
      ! First get the length of the attribute
      res = nf_inq_attlen (ncid, varid, name, len)
      IF(res == nf_noerr) &
         res = nf_get_att_int(ncid, varid, name, ivals)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_input_bcast)
   p_nf_get_att_int_1 = res

   ! If there was an error, don't try to broadcast the values

   IF(res /= nf_noerr) return

   ! Broadcast number of values and values themselves

   CALL p_bcast(len, p_io, p_comm_input_bcast)
   CALL p_bcast(ivals(1:len), p_io, p_comm_input_bcast)

END FUNCTION p_nf_get_att_int_1

!-------------------------------------------------------------------------
!
!

!>
!!               Wrapper for nf_get_var_int.
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Added p_comm_input_bcast by Rainer Johanni, Oct 2010
!!
INTEGER FUNCTION p_nf_get_var_int(ncid, varid, ivals)

!
   INTEGER, INTENT(IN) :: ncid, varid
   INTEGER, INTENT(OUT) :: ivals(*)

   INTEGER :: res, len, ndims, dimids(NF_MAX_VAR_DIMS), dimlen, i


!-----------------------------------------------------------------------

   IF(p_pe == p_io) THEN

      ! First get the length of the array

      res = nf_inq_varndims(ncid, varid, ndims)
      IF(res /= nf_noerr) GOTO 9999
      res = nf_inq_vardimid(ncid, varid, dimids)
      IF(res /= nf_noerr) GOTO 9999

      len = 1
      DO i = 1, ndims
         res = nf_inq_dimlen(ncid, dimids(i), dimlen)
         IF(res /= nf_noerr) GOTO 9999
         len = len * dimlen
      ENDDO

      res = nf_get_var_int(ncid, varid, ivals)

   ENDIF

9999 CONTINUE

   CALL p_bcast(res, p_io, p_comm_input_bcast)
   p_nf_get_var_int = res

   ! If there was an error, don't try to broadcast the values

   IF(res /= nf_noerr) return

   ! Broadcast number of values and values themselves

   CALL p_bcast(len, p_io, p_comm_input_bcast)
   CALL p_bcast(ivals(1:len), p_io, p_comm_input_bcast)

END FUNCTION p_nf_get_var_int

!-------------------------------------------------------------------------
!
!

!>
!!               Wrapper for nf_get_var_double.
!!
!!
!! @par Revision History
!! Initial version by Rainer Johanni, Nov 2009
!! Added p_comm_input_bcast by Rainer Johanni, Oct 2010
!!
INTEGER FUNCTION p_nf_get_var_double(ncid, varid, dvals)

!
   INTEGER,  INTENT(IN)  :: ncid, varid
   REAL(dp), INTENT(OUT) :: dvals(*)

   INTEGER :: res, len, ndims, dimids(NF_MAX_VAR_DIMS), dimlen, i


!-----------------------------------------------------------------------

   IF(p_pe == p_io) THEN

      ! First get the length of the array

      res = nf_inq_varndims(ncid, varid, ndims)
      IF(res /= nf_noerr) GOTO 9999
      res = nf_inq_vardimid(ncid, varid, dimids)
      IF(res /= nf_noerr) GOTO 9999

      len = 1
      DO i = 1, ndims
         res = nf_inq_dimlen(ncid, dimids(i), dimlen)
         IF(res /= nf_noerr) GOTO 9999
         len = len * dimlen
      ENDDO

      res = nf_get_var_double(ncid, varid, dvals)

   ENDIF

9999 CONTINUE

   CALL p_bcast(res, p_io, p_comm_input_bcast)
   p_nf_get_var_double = res

   ! If there was an error, don't try to broadcast the values

   IF(res /= nf_noerr) return

   ! Broadcast number of values and values themselves

   CALL p_bcast(len, p_io, p_comm_input_bcast)
   CALL p_bcast(dvals(1:len), p_io, p_comm_input_bcast)

END FUNCTION p_nf_get_var_double

!-------------------------------------------------------------------------
!
!

!>
!!               Wrapper for nf_get_vara_double.
!!
!!
!! @par Revision History
!! Initial version by Marco Giorgetta, Sept 2010
!! Added p_comm_input_bcast by Rainer Johanni, Oct 2010
!!
INTEGER FUNCTION p_nf_get_vara_double(ncid, varid, start, count, dvals)

!
   INTEGER,  INTENT(IN)  :: ncid, varid, start(*), count(*)
   REAL(dp), INTENT(OUT) :: dvals(*)

   INTEGER :: res, len, ndims, dimids(NF_MAX_VAR_DIMS), i


!-----------------------------------------------------------------------

   IF(p_pe == p_io) THEN

      ! First get the length of the array

      res = nf_inq_varndims(ncid, varid, ndims)
      IF(res /= nf_noerr) GOTO 9999
      res = nf_inq_vardimid(ncid, varid, dimids)
      IF(res /= nf_noerr) GOTO 9999

      len = 1
      DO i = 1, ndims
         len = len * count(i)
      ENDDO

      res = nf_get_vara_double(ncid, varid, start, count, dvals)

   ENDIF

9999 CONTINUE

   CALL p_bcast(res, p_io, p_comm_input_bcast)
   p_nf_get_vara_double = res

   ! If there was an error, don't try to broadcast the values

   IF(res /= nf_noerr) return

   ! Broadcast number of values and values themselves

   CALL p_bcast(len, p_io, p_comm_input_bcast)
   CALL p_bcast(dvals(1:len), p_io, p_comm_input_bcast)

END FUNCTION p_nf_get_vara_double

!-------------------------------------------------------------------------

END MODULE mo_read_netcdf_parallel
