/* --------------------------------------------------------------------- */
/* Copyright 2006-2010, NEC Europe Ltd., London, UK. */
/* All rights reserved. Use is subject to OASIS4 license terms. */
/* --------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "psmile_f2c.h"

/*

//BOP
 
! !ROUTINE: psmile_redirstdout
 
  !INTERFACE:

*/

void psmile_redirstdout (char *filestem, INTEGER *lenstr,
                         INTEGER *parallel, INTEGER *my_pe,
                         INTEGER *npes, INTEGER *ierror )
/*

  !INPUT PARAMETERS:

   filestem = Basename of stdout/stderr files
   lenstr   = Significant length (without \0) of "filestem"
   
   parallel = Flag for direction of output
              parallel = 1 : Create single output files for each process
              Otherwise    : Create output files only for process 0.
                             The output of other processes is not redirected.

   my_pe    = Rank

   npes     = number of processes

   !OUTPUT PARAMETERS:

   ierror   = Return code
              = 0 : No error
              = 1 : Error in malloc

  !DESCRIPTION:

   Redirect standard output

  !FILES USED:

         <math.h>
         <stdio.h>
         <stddef.h>
         "psmile_f2c.h"

  !REVISION HISTORY:

    Date      Programmer   Description
  ----------  ----------   -----------
  01.12.03    R. Redler    created
  01.12.05    H. Ritzdorf  revised

//EOP

 ----------------------------------------------------------------------
  $Id
  $Author
 ---------------------------------------------------------------------- */

{
   size_t len_alloc = *lenstr;
   register int n_digits, num;
   register int redirect = 1;

   char *sname, *ename;
   register char *p_err, *p_std;

   ASSERT (*lenstr > 0)
   ASSERT (*my_pe  >= 0)

   *ierror = 0;

   /* Get number of digits in rank */

   num = (*npes > 0) ? *npes : 1;
   n_digits = (int) (log10((double)(num) + (double) 0.5)) + 1;
   ASSERT (n_digits > 0)

   /* allocate file names */

#define NUMBER_OF_EXTENSIONS 2
#define LENGTH_OF_ERR 3

   len_alloc += n_digits + NUMBER_OF_EXTENSIONS + LENGTH_OF_ERR + 1;

   sname = (char *) MALLOC (len_alloc);
   ename = (char *) MALLOC (len_alloc);

   if (!sname || !ename) {
      fprintf (stderr, "PSMILe_redirstdout: Cannot allocate memory. %d bytes\n", len_alloc);
      *ierror = 1;
      return;
   }

   /* Copy base name to standard output and standard error file name */

   memcpy (sname, filestem, (size_t) *lenstr);
   memcpy (ename, filestem, (size_t) *lenstr);

   p_std = sname + (size_t) *lenstr;
   p_err = ename + (size_t) *lenstr;

   if(* parallel == 1) {
     /* all processes write into their own file */

     sprintf (p_std, ".%0*d",     n_digits, *my_pe);
     sprintf (p_err, ".err.%0*d", n_digits, *my_pe);
   }
   else {
     /* only local root redirectes stdout */
     redirect = (*my_pe == 0);

     if( redirect ) {
       sprintf (p_std, ".log");
       sprintf (p_err, ".err");
     }
   }

   /* Redirect STDOUT/STDERR */

   if (redirect) {
      freopen (sname, "w", stdout);
      freopen (ename, "w", stderr);
   }

   /* Free file names */

   FREE (ename);
   FREE (sname);
}
