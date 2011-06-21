!-----------------------------------------------------------------------
! Copyright 2006-2010, NEC Europe Ltd., London, UK.
! All rights reserved. Use is subject to OASIS4 license terms.
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: PSMILe_Char2buf
!
! !INTERFACE:
!
      SUBROUTINE psmile_char2buf (ilubuf, ndibuf, ipos, string)
!
! !USES:
!
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      CHARACTER (len=*), INTENT(in) :: string

!     Character string, to be written to buffer.

      INTEGER, INTENT(in)           :: ndibuf

!     Dimension of buffer 'ILUBUF'.

!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER, INTENT(inout)        :: ipos

!     Pointer to the last used index in 'ILUBUF'
!     Length of 'ILUBUF' currently used.

      INTEGER, INTENT(inout)        :: ilubuf (ndibuf)

!     Integer buffer.
!
! !LOCAL VARIABLES
!
      CHARACTER (len=12)      :: form
      INTEGER                 :: i, irest, lav, lenstr
      INTEGER                 :: length, length_of_integer
!
! !DESCRIPTION:
!
!  Write Character String to Buffer
!
!  Subroutine "PSMILe_Char2buf" writes the character string 'STRING'
!  to buffer 'ILUBUF(IPOS+1:)' and 'IPOS' is updated to point to the
!  last index used in buffer 'ILUBUF'.
!
!  If character string 'STRING' is too long to be stored entirely in
!  the buffer, the storable length is stored in vector 'ILUBUF' and
!  'IPOS' points to the last or first index which corresponds to
!  the index in 'ILUBUF' necessary to store the string entirely.
!
! !REVISION HISTORY:
!
!   Date      Programmer   Description
! ----------  ----------   -----------
! 01.12.03    H. Ritzdorf  created
!
!EOP
!-----------------------------------------------------------------------
!
!  $Id$
!  $Autor$

!-----------------------------------------------------------------------
!
!===> Initialization
!
      length_of_integer = BIT_SIZE (i) / 8
      lenstr = LEN (string)
      length = (lenstr-1) / length_of_integer
!
      lav = ndibuf - ipos
!
!===> Create format ``form''
!
      irest = lenstr - length*length_of_integer
!
      IF (length .GT. 0) THEN
         WRITE (form, 9990) length, length_of_integer, irest
      ELSE
         WRITE (form, 9980) irest
      ENDIF

write ( * , * ) form

!
!===> Write to buffer ``ilubuf''
!
      READ (string, form) (ilubuf(ipos+i), i= 1, MIN (length+1, lav))
!
      ipos = ipos + (length+1)
!
!  Formats:
!
9990  FORMAT ('(', i4, 'a', i1, ', a', i1, ')')
9980  FORMAT ('(a', i1, ')')

      END SUBROUTINE psmile_char2buf
