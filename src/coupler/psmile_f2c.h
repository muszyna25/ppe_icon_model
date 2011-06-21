/**********************************************************************/

/* Copyright 2006, C&C Research Laboratories, NEC Europe Ltd., St. Augustin, Germany. */
/* All rights reserved. Use is subject to OASIS4 license terms. */

/**********************************************************************

!----------------------------------------------------------------------
!BOP
!
! !INCLUDE:  PSMILe_f2c.h
!
! !DESCRIPTION:
! 
  Include file of C-Interface of PSMILe Library:

  This include-file contains the installation dependent definitions
  of the FORTRAN to C-Interface of the PSMILe Library.

! !REVISION HISTORY:
!
!   Date      Programmer   Description
! ----------  ----------   -----------
! 01.12.03    H. Ritzdorf  created
!
!EOP
!-----------------------------------------------------------------------

  ----------------------------------------------------------------------

  $Id: PSMILe_f2c.h,v 1.8.2.2 2010-01-20 13:25:35 m300083 Exp $

 ***********************************************************************/

/* =====================================================================
   Definitions of Variables used in the Interface between Fortran and C

   INTEGER    : Corresponding type to Fortran type integer
   REAL       : Corresponding type to Fortran type real

   LONG_INT_C : C type which corresponds to the Fortran type
                LONG_INT_F.
                LONG_INT_F must be large enough to store relative
                addresses which are pointers to allocated buffers.
   ===================================================================== */

#ifndef INTEGER
#  if defined ( PRISM_EXTENDED_WIDTH )
/*    Fortran INTEGER corresponds to C's long long */
#     define INTEGER              long long
#  else
#     define INTEGER              int
#  endif
#endif /* INTEGER */

#ifndef REAL
#  if defined ( PRISM_EXTENDED_WIDTH )
/*    Fortran REAL corresponds to C's double */
#     define REAL                 double
#  else
#     define REAL                 float
#  endif
#endif /* REAL */

#ifndef LONG_INT_C
#   if defined (POINTER_64_BITS)
#      define LONG_INT_C           long long
#   else
#      define LONG_INT_C           INTEGER
#   endif
#endif

#ifndef MALLOC

#   define MALLOC(size) malloc ((size_t) size)
#   define FREE(ptr)    free (ptr)

#endif /* MALLOC */

#ifdef PRISM_ASSERTION
#   undef ASSERT
#ifdef __ANSI_CPP__
#   define ASSERT(c) \
if (!(c)) {\
   fprintf(stderr, "### Assertion violation: %s in %s:%d\n",\
           #c, __FILE__, __LINE__);\
   abort ();\
}
#else
#   define ASSERT(c) \
if (!(c)) {\
   fprintf(stderr, "### Assertion violation: %s in %s:%d\n",\
           "expression", __FILE__, __LINE__);\
   abort ();\
}
#endif /* __ANSI_CPP__ */

#else

#   define ASSERT(c) 

#endif /* PRISM_ASSERTION */

/* =====================================================================
   Names of C-Routines which are called by FORTRAN routines
   ===================================================================== */

#if defined (FORTRANCAPS)

#  define psmile_bsend       PSMILE_BSEND
#  define psmile_bsend_init  PSMILE_BSEND_INIT
#  define psmile_redirstdout PSMILE_REDIRSTDOUT

#elif defined (FORTRANDOUBLEUNDERSCORE)

#  define psmile_bsend       psmile_bsend__
#  define psmile_bsend_init  psmile_bsend_init__
#  define psmile_redirstdout psmile_redirstdout__

#else

#  define psmile_bsend       psmile_bsend_
#  define psmile_bsend_init  psmile_bsend_init_
#  define psmile_redirstdout psmile_redirstdout_

#endif
