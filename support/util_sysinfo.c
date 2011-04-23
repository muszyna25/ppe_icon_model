/*
 * This file contains a set of routines to provide operating system information
 * for f90 on Unix machines. Unfortunatelly there is no implemented set of
 * such functions available. To restrict the f90 source code to f90 standard
 * only this mostly POSIX C compliant implementations are developed. One
 * basic convention is to name all routines with a leading util_ to make
 * this visible to the developer.
 *
 * Authors: L. Kornblueh      Max-Planck-Institute for Meteorology, Hamburg
 *          U. Schulzweida    Max-Planck-Institute for Meteorology, Hamburg
 *
 * Date:    7.5.1998
 *
 *
 */

#include "config.h"

#include <sys/param.h>
#include <sys/utsname.h>
#ifndef __XT3__
#include <netdb.h>
#endif

#include <sys/types.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#elif HAVE_SYS_UNISTD_H
#include <sys/unistd.h>
#endif
 
#include <pwd.h>
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#ifdef __uxp__
#include <sys/systeminfo.h>
#endif

#include "cfortran.h"


/****************************************************************************/
/* cfortran prototypes:                                                     */

void cf_util_os_system(char *name, int *actual_len);
FCALLSCSUB2(cf_util_os_system, UTIL_OS_SYSTEM, util_os_system,
	    PSTRING, PINT)
void cf_util_user_name(char *name, int *actual_len);
FCALLSCSUB2(cf_util_user_name, UTIL_USER_NAME, util_user_name,
	    PSTRING, PINT)
void cf_util_node_name(char *name, int *actual_len);
FCALLSCSUB2(cf_util_node_name, UTIL_NODE_NAME, util_node_name,
	    PSTRING, PINT)

/****************************************************************************/

/* funcion implemetations */ 

void cf_util_os_system(char *name, int *actual_len)
{
  struct utsname utname;
 
  uname( &utname );

  /* Cray doesn't handle operating system information proper */

#ifdef _UNICOSMP
  strcpy(name, "UNICOS");
#else
  strcpy(name, utname.sysname);
#endif
  strcat(name, " ");
  strcat(name, utname.release);
  strcat(name, " on ");
  strcat(name, utname.machine);

  *actual_len = strlen(name);
}

void cf_util_user_name(char *name, int *actual_len)
{
  struct passwd *current;
 
#ifdef __XT3__
  current = NULL;
#else
  current = getpwuid(getuid());
#endif
  if (current == NULL) {
#ifdef __XT3__
    strcpy(name, "user name not available");
#else
    strcpy(name, "unknown user name");
#endif
  } else {
    if (strlen(current->pw_name) == 0) {
      strcpy(name, "unknown user name");
    } else {
#ifdef __hpux
      strncpy(name, current->pw_gecos, strcspn(current->pw_gecos, ","));
#else      
      strcpy(name, current->pw_gecos);
#endif
      strcat(name, " (");
      strcat(name, current->pw_name);
      strcat(name, ")");
    }
  }    

  *actual_len = strlen(name);
}

void cf_util_node_name(char *name, int *actual_len)
{
  char *hostname;

  if ((hostname = getenv ("HOST")) == NULL) {
    if ((hostname = getenv ("HOSTNAME")) == NULL) {
      strcpy(name, "unknown");
    } else {
      strcpy(name, hostname);
    }
  } else {
    strcpy (name, hostname);
  }


  *actual_len = strlen(name);
}






