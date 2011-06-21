/* --------------------------------------------------------------------- */
/* $COPYRIGHT_CCRL */
/* --------------------------------------------------------------------- */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#ifndef NOMPI
#include "mpi.h"
#endif
#include "psmile_f2c.h"

#ifndef PSMILE_BSEND

#define PRISM_Error_Alloc 13

#define MAX(a,b) (a) > (b) ? (a) : (b)

#   define ASSERT2(c, a, b) \
if (!(c)) {\
   fprintf(stderr, "### Assertion violation: %s (%s = %d, %s = %d) in %s:%d\n", #c, #a, a, #b, b, __FILE__, __LINE__);\
   abort();\
}

/* -----------------------------------------------------------------------
   Local parameters :
 
   n_to_free = Number of message buffers from which is freed
   dfree     = Length of buffers allocated from which is freed
   N_HASH    = Number of hash values
   ----------------------------------------------------------------------- */

#ifdef DEBUG
#   define n_to_free 0
#else
#   define n_to_free 64
#endif /* DEBUG */

#define dfree  64*2048
#define N_HASH 16 /* Test these value with OpenMPI */

/* -----------------------------------------------------------------------
   Local variables :
 
   requests  = List of MPI_Isend requests
               Dimension: request [max_alloc]
   buffers   = List of MPI_Isend buffers
               Dimension: buffers [max_alloc]
   lengths   = Lengths of MPI_Isend buffers
               Dimension: lengths [max_alloc]

   max_alloc = Dimension of "requests", "buffers" and "lengths"
   nalloc    = Number of currently used requests in "requests", "buffers"
               and "lengths".

   dalloc    = Total number of allocated data in Bytes.

   min_type  = Minimum index of Fortran datatypes
   max_type  = Maximum index of Fortran datatypes
   dbsend    = Lengths of Fortran datatypes in bytes if use_hash == 0.
               Dimension: dbsend [max_type-mintype+1]
               Set in psmile_bsend_init ().

   use_hash  = Use hash table to translate Fortran datatypes
   dbtype    = Hashed list   of Fortran datatypes
   dblen     = Hashed length of Fortran datatypes
   ----------------------------------------------------------------------- */

static double dalloc = 0.0;

#ifndef NOMPI
static MPI_Request   *requests = NULL;
#endif
static void         **buffers = NULL;
static long long     *lengths = NULL;
static int           nalloc = 0, max_alloc = 0;

static int           use_hash = 0;
static int           max_type = -1, min_type = 0;
static int           *dbsend = NULL;
static int           *dbtype = NULL;
static int           *dblen  = NULL;

/* -----------------------------------------------------------------------
   Internal functions
   ----------------------------------------------------------------------- */

	void psmile_bsend_init (INTEGER *ftypes, INTEGER *flengths,
                                INTEGER *number_of_ftypes, INTEGER *ierror);

        void psmile_bsend (void *buf, INTEGER *lenbuf, INTEGER *dtype,
                           INTEGER *dest, INTEGER *tag, INTEGER *comm,
		           INTEGER *ierror);

static	void	error_in_testany (int imin, int ind, int error);
static	int 	free_buffers (long long len_bytes, int free_memory);
static 	int  	release_buffer (int i, long long len_bytes, int *imin,
                                int free_memory);

/* -----------------------------------------------------------------------
//BOP
 
! !ROUTINE: psmile_bsend
 
  !INTERFACE:
*/
 
#ifdef ANSI_C

void psmile_bsend (void *buf, INTEGER *lenbuf, INTEGER *dtype,
                   INTEGER *dest, INTEGER *tag, INTEGER *comm,
		   INTEGER *ierror)

#else

void psmile_bsend (buf, lenbuf, dtype,
                   dest, tag, comm, ierror)

void    *buf;
INTEGER *lenbuf, *dtype, *dest, *tag, *comm, *ierror;

#endif
/*
 
  !INPUT PARAMETERS:
 
      void, Intent (In)               :: buf (*)
 
      Message buffer to be sent
 
      Integer, Intent (In)            :: lenbuf
 
      Length of the buffer (of type "dtype") to be sent.
 
      Integer, Intent (In)            :: dtype
 
      Fortran Datatype of buffer
 
      Integer, Intent (In)            :: dest
 
      Destination of message.
 
      Integer, Intent (In)            :: tag
 
      Message tag to be used
 
      Integer, Intent (In)            :: comm
 
      Communicator
 
  !RETURN VALUES:
  
      Integer, Intent (Out)           :: ierror
 
      Returns the error code of psmile_bsend:
              ierror = 0 : No error
              ierror > 0 : Severe error
 
  !DESCRIPTION:
 
   Taken from the OASIS4 distribution.
   Subroutine "psmile_bsend" performs a buffered send.
   It allocates a message buffer of the length "lenbuf",
      copies buf(1:lenbuf) into this message buffer and
      sends  this message to the destination process using
      non-blocking sends.

  !FILES USED:

         <stdio.h>
         <stddef.h>
         <stdlib.h>
         <string.h>

         "mpi.h"
         "PSMILe_f2c.h"

  !REVISION HISTORY:

    Date      Programmer   Description
  ----------  ----------   -----------
  06.07.03    H. Ritzdorf  created
 
//EOP

 ----------------------------------------------------------------------
  $Id: psmile_bsend.c,v 1.16.2.5 2009-10-05 09:00:14 m300083 Exp $
  $Author: m300083 $
 ---------------------------------------------------------------------- */

{
#ifndef NOMPI
   long long len_bytes; /* Message length in Bytes */
   MPI_Datatype dtypec = MPI_Type_f2c(*dtype);
   MPI_Comm commc  = MPI_Comm_f2c(*comm);
   
   /* Internal control */

   ASSERT (*dtype >= min_type && *dtype <= max_type)

   if (use_hash) {
      register int i, j = (*dtype) % N_HASH;
#pragma vdir vector
      for (i = dbsend[j]; i < dbsend[j+1]; j++)
         if (dbtype[i] == (*dtype)) break;

      ASSERT (i < dbsend[j+1]);
      len_bytes = *lenbuf * dblen[i];
   }
   else {
      ASSERT (dbsend [*dtype-min_type] > 0)
      len_bytes = *lenbuf * dbsend [*dtype-min_type];
   }

#ifdef VERBOSE_COMM
#define VERBOSE_COMM
   fprintf (stdout, "-> %d; tag %d comm %d lenbuf %d %lld\n", 
	            *dest, *tag, *comm, *lenbuf, len_bytes);
   fflush (stdout);
#else /* VERBOSE_COMM */
 
#ifdef DEBUG
   if (*lenbuf < 1) {
      int rank;
      MPI_Comm_rank (commc, &rank);
      fprintf (stderr, "(%d)-> psmile_bsend: dest = %d, tag = %d, comm %d, lenbuf %d\n",
	       rank, *dest, *tag, *comm, *lenbuf);
   }
#endif /* DEBUG */
#endif /* VERBOSE_COMM */

/* =======================================================================
   Special case : Can the message be sent by mpi_send ?
   ======================================================================= */
 
#ifdef PSMILE_SEND_LENBUF
 
   if (len_bytes <= PSMILE_SEND_LENBUF) {
      *ierror = MPI_Send (buf, (int) *lenbuf, dtypec,
	                       (int) *dest,   (int) *tag,   commc);

      return;
   }
 
#else  /* PSMILE_SEND_LENBUF */
/* =======================================================================
   Special case : Empty message
   ======================================================================= */

   if (*lenbuf <= 0) {
      MPI_Request lrequest;

      *ierror = MPI_Isend (buf, (int) *lenbuf, dtypec,
	                   (int) *dest, (int) *tag, commc, &lrequest);
      if (*ierror != MPI_SUCCESS) return;
 
      *ierror = MPI_Request_free (&lrequest);
      return;
   }
#endif /* PSMILE_SEND_LENBUF */

/* =======================================================================
   Send message using non blocking send

   ======================================================================= */
 
   /* Look for a free message buffer.
      (i) Control message buffers allocated
 
      buffers (i) = reference address
      lengths (i) = size of buffer allocated
   */
 
   {
      int imin = -1;
      int free_memory = nalloc > n_to_free || dalloc > dfree;

      if (free_memory) {
	 imin = free_buffers (len_bytes, free_memory);
         ASSERT2 (imin < nalloc, imin, nalloc);
      }

   /* ==============================================================
      Send message
 
      imin = -1 : Allocate message buffer
                  Otherwise, use old message buffer
      ============================================================== */
 
      if (imin == -1) {
	 if (nalloc == max_alloc) {
	    register int i;
	    max_alloc += 32;
	    if (nalloc == 0) {
	       requests = (MPI_Request *) MALLOC (max_alloc * sizeof(MPI_Request));
	       buffers  = (void **)     MALLOC (max_alloc * sizeof(void *));
	       lengths  = (long long *) MALLOC (max_alloc * sizeof(long long));
	    }
	    else {
	       requests = (MPI_Request *) realloc (requests, max_alloc * sizeof(MPI_Request));
	       buffers  = (void **)       realloc (buffers, max_alloc * sizeof(void *));
	       lengths  = (long long *)   realloc (lengths, max_alloc * sizeof(long long));
	    }

	    if (! requests || ! buffers || ! lengths) {
	       fprintf (stderr, "pmsile_bsend: Cannot allocate %d bytes in order to allocate requests\n",
		       max_alloc*(sizeof(MPI_Request)+sizeof(void *)+sizeof(long long)));
               *ierror = MAX (1307, MPI_ERR_LASTCODE+10);
	       return;
	    }

            /* Initialize newly allocated parts of control vectors */
#pragma vdir vector
	    for (i=nalloc; i < max_alloc; i++) {
	       requests[i] = MPI_REQUEST_NULL;
	       buffers[i] = (void *) NULL;
	       lengths[i] = 0;
	    }
	 }

         /* Allocate a new buffer */
 
allocate_buffer:
	 buffers [nalloc] = (void *) MALLOC (len_bytes);
	 if (buffers [nalloc] == NULL) {
            if (! free_memory) {
	       /* Try to free buffers */
               free_memory = 1;
	       imin = free_buffers (len_bytes, free_memory);
               if (imin == -1) goto allocate_buffer;
               ASSERT2 (imin < nalloc, imin, nalloc);
	    }
	    else {
	       fprintf (stderr, "pmsile_bsend: Cannot allocate %lld bytes in order to send buffer\n",
		       len_bytes);
               *ierror = MAX (1307, MPI_ERR_LASTCODE+10);
	       return;
	    }
	 }
	 else {
            /* Store data of buffer newly allocated */

            imin = nalloc;
            nalloc ++;
            dalloc += len_bytes;
            lengths [ imin] = len_bytes;
	 }
      }

      /* Internal Control */

      ASSERT2 (imin < nalloc, imin, nalloc);
      ASSERT (imin >= 0);
      ASSERT (lengths[imin] >= len_bytes);
      ASSERT (requests[imin] == MPI_REQUEST_NULL);
 
      /* Copy and send message */

      memcpy (buffers[imin], buf, (size_t) len_bytes);
 
      *ierror = MPI_Isend (buffers[imin], *lenbuf, dtypec,
                           *dest, *tag, commc, &requests[imin]);

      return;
   }
}

/* =========================================================================
   Initialize datatype information for "psmile_bsend"

  !PARAMETERS:

   ftypes          : Vector containing the FORTRAN MPI Datatypes
                     Dimension: ftypes (number_of_ftypes)
   flengths        : Vector containing the size of the datatype in bytes
                     Dimension: flengths (number_of_ftypes)
   number_of_ftypes: Number of FORTRAN MPI Datatypes
                     Dimension: flengths (number_of_ftypes)

 
  !RETURN VALUES:
  
      Integer, Intent (Out)           :: ierror
 
      Returns the error code of psmile_bsend_init:
              ierror = 0 : No error
              ierror > 0 : Severe error
   ========================================================================= */
 
#ifdef ANSI_C

void psmile_bsend_init (INTEGER *ftypes, INTEGER *flengths,
                        INTEGER *number_of_ftypes, INTEGER *ierror)

#else

void psmile_bsend_init (ftypes, flengths, number_of_ftypes,
                        ierror)

INTEGER *ftypes, *flengths, *number_of_ftypes, *ierror;

#endif

{
   register int i, j, k, ndbsnd;
   register int n_ftypes = *number_of_ftypes;
   int *num = NULL;

   max_type = min_type = ftypes[0];
#pragma vdir vector
   for (i=1; i < n_ftypes; i++) 
      max_type = ftypes[i] > max_type ? ftypes [i] : max_type;

#pragma vdir vector
   for (i=1; i < n_ftypes; i++) 
      min_type = ftypes[i] < min_type ? ftypes [i] : min_type;

   use_hash = max_type - min_type > 64;
   if (! use_hash) {
      /* Small datatype range */
      ndbsnd = max_type - min_type + 1;
      dbsend = (int *) MALLOC (ndbsnd*sizeof(int));
      if (! dbsend) {
         fprintf (stderr, "Error in psmile_bsend_init: Cannot allocate %d bytes\n",
	                  ndbsnd*sizeof(int));
         *ierror = PRISM_Error_Alloc;
         return;
      }

#pragma vdir vector
      for (i = 0; i < ndbsnd; i++) dbsend [i] = 0;

#pragma vdir vector
      for (i = 0; i < n_ftypes; i++)
         dbsend [ftypes[i]-min_type] = flengths[i];
   }
   else {
      /* Large datatype range; generate simple hash table */
      ndbsnd = N_HASH + 1;

      dbsend = (int *) MALLOC ((ndbsnd+2*n_ftypes)*sizeof(int));
      if (! dbsend) {
         fprintf (stderr, "Error in psmile_bsend_init: Cannot allocate %d bytes\n",
	                  (ndbsnd+2*n_ftypes)*sizeof(int));
         *ierror = PRISM_Error_Alloc;
         return;
      }

      num = (int *) MALLOC (ndbsnd*sizeof(int));
      if (! num) {
         fprintf (stderr, "Error in psmile_bsend_init: Cannot allocate %d bytes\n",
	                  ndbsnd*sizeof(int));
         *ierror = PRISM_Error_Alloc;
         return;
      }

      dbtype = dbsend + ndbsnd;
      dblen  = dbtype + n_ftypes;

#pragma vdir vector
      for (j = 0; j < ndbsnd; j++) 
         dbsend [j] = 0;

#pragma vdir vector
      for (i = 0; i < n_ftypes; i++) {
         j = ftypes[i] % N_HASH;
         dbsend[j+1] ++;
      }

      /* dbsend[j] = start index of ftypes, where ftypes[i]%N_HASH == j
                     end index is dbsend[j+1]-1 */

#pragma vdir vector
      for (j = 1; j < ndbsnd; j++) {
         dbsend[j] += dbsend[j-1];
      }

#pragma vdir vector
      for (j = 0; j < ndbsnd; j++) {
         num[j] = 0;
      }

#pragma vdir vector
      for (i = 0; i < n_ftypes; i++) {
         j = ftypes[i] % N_HASH;
	 k = dbsend[j] + num[j]; /* Index in dbtype and dblen */
         num[j] ++;

         dbtype[k] = ftypes[i];
         dblen [k] = flengths[i];
      }

#ifdef DEBUGXX
      fprintf (stderr, "ftypes:");
      for (i = 0; i < n_ftypes; i++) fprintf (stderr, "%d ", ftypes[i]);
      fprintf (stderr, "\n");

      fprintf (stderr, "flengths:");
      for (i = 0; i < n_ftypes; i++) fprintf (stderr, "%d ", flengths[i]);
      fprintf (stderr, "\n");

      fprintf (stderr, "dbsend:");
      for (i = 0; i < ndbsnd; i++) fprintf (stderr, "%d ", dbsend[i]);
      fprintf (stderr, "\n");

      fprintf (stderr, "dbtype:");
      for (i = 0; i < n_ftypes; i++) fprintf (stderr, "%d ", dbtype[i]);
      fprintf (stderr, "\n");

      fprintf (stderr, "dblen:");
      for (i = 0; i < n_ftypes; i++) fprintf (stderr, "%d ", dblen[i]);
      fprintf (stderr, "\n");
#endif /* DEBUGXX */

      /* Free temporary space */

      FREE (num);
   }

   *ierror = 0;
   return;
}

/* =========================================================================
   Look for an appropriate buffer and free buffers.

   Input Arguments:

   len_bytes   = (Minimal) Length of an appropriate buffer.
   free_memory = Should the memory be freed ?

   Function "free_buffers" returns the index of an appropriate buffer
   if such an buffer is found.
   Otherwise, -1 is returned.
   ========================================================================= */

int free_buffers (long long len_bytes, int free_memory)
{
   int imin = -1;
   int ibeg = 0;
   int error, flag;
   MPI_Status lstatus;
 
   /* --------------------------------------------------------------
      Look within the requests which are already fulfilled
      -------------------------------------------------------------- */

   while (ibeg < nalloc) {
      register int i;

#pragma vdir vector
      for (i=ibeg; i < nalloc; i++) {
         if (requests [i] == MPI_REQUEST_NULL) break;
      }

      if (i >= nalloc) break;


      ibeg = release_buffer (i, len_bytes, &imin, free_memory);
      if (imin >= 0 && lengths[imin] == len_bytes) {
	 return imin;
      }

   } /* end while */

   /* All buffers are controlled
      If a valid buffer IMIN was found, send message */
 
   if (imin != -1) return imin;
 
   /* ------------------------------------------------
      Look in requests which are not already fulfilled
      ------------------------------------------------ */
 
   ibeg = 0;

   while (ibeg < nalloc) {
      int ind;
      error = MPI_Testany (nalloc-ibeg, requests+ibeg, &ind, &flag,
                           &lstatus);
      if (error != MPI_SUCCESS) {
         error_in_testany (ibeg, ind, error);
         return imin;
      }

      /* no request fullfilled ? */
 

      if (! flag || ind == MPI_UNDEFINED) return imin;

      /* Request ibeg+ind was finished */
 
      ind = ibeg + ind;

      /* Control MPI function */

      ASSERT (requests[ind] == MPI_REQUEST_NULL)
 
      ibeg = release_buffer (ind, len_bytes, &imin, free_memory);
      /* if (lengths[imin] == len_bytes && ! free_memory) return imin; */
      if (imin >= 0 && lengths[imin] == len_bytes) return imin;
   }

   return imin;
}

/* =========================================================================
   Release buffer with index I

   Arguments:

   i           = Index of buffer to be freed.
   len_bytes   = Length of buffer required to be sent.
   imin        = Index of the buffer which is free and 
                 whose size if greater than "len_bytes".
   free_memory = Should the memory be freed ?

   Returns next index to be controlled; i.e. it returns
    (*) I if buffer with index I was removed and
	     the index was filled by the last buffer
    (*) I++ otherwise.
   ========================================================================= */

static int release_buffer (int i, long long len_bytes, int *imin,
                           int free_memory)
{
   ASSERT (0 <= i && i < nalloc)
   ASSERT (i != *imin)

   if (lengths [i] < len_bytes ||
       (*imin >= 0 && lengths[i] >= lengths[*imin])) {
 
      /* Buffer is too small or too large */
 
      if (free_memory) {
         FREE (buffers[i]);
 
         dalloc -= lengths [i];
         nalloc --;

         buffers [i] = buffers [nalloc];
         lengths [i] = lengths [nalloc];
 
         requests [i] = requests [nalloc];

#ifdef PRISM_ASSERTION
	 /* for the Assertions */
         requests [nalloc] = MPI_REQUEST_NULL;
#endif
      }
      else i++;
   }
   else if (free_memory && *imin >= 0) {
      /* ... Buffer "I" is smaller than buffer "*IMIN".
             Release buffer "*IMIN" */
 
      ASSERT2 (*imin < nalloc, *imin, nalloc)
      ASSERT (*imin >= 0)
      ASSERT (lengths[*imin] >= len_bytes)
      ASSERT (requests [*imin] == MPI_REQUEST_NULL)
      ASSERT (requests [i]    == MPI_REQUEST_NULL)
 
      FREE (buffers[*imin]);

      dalloc -= lengths [*imin];
      nalloc --;

      /* Replace entry *imin */

      buffers [*imin] = buffers [i];
      lengths [*imin] = lengths [i];

      /* Move last entry */
 
      buffers [i] = buffers [nalloc];
      lengths [i] = lengths [nalloc];
 
      requests [i] = requests [nalloc];

#ifdef PRISM_ASSERTION
      /* for the Assertions */
      requests [nalloc] = MPI_REQUEST_NULL;
#endif
   }
   else {
      *imin = i;
      i++;
   }

   ASSERT2 (*imin < nalloc, *imin, nalloc)

   return i;
}

/* =========================================================================
   Error in MPI_Testany
   Used for internal debugging and 
   some MPI implementations had problems with MPI_Testany ()
   ========================================================================= */

static void error_in_testany (int ibeg, int ind, int error)
{
   register int i;
   fprintf (stderr, "\nError in MPI_Testany () called form psmile_bsend: error %d\n",
	           error);
   fprintf (stderr, "ind %d request %d nalloc %d ibeg %d:\n",
	             ind, requests[ind], nalloc, ibeg);
 
   for (i=0; i < nalloc; i++) {
      fprintf (stderr, "i = %d, start 0x%p len %lld, request %d\n",
	               i, buffers [i], lengths [i], requests[i]);
   }
 
   /* Set ``requests [i]'' to NULL to avoid indefinite loop */
 
   requests [ind] = MPI_REQUEST_NULL;
#endif  /* ! defined NOMPI */
}
#endif /* ! defined psmile_bsend */
