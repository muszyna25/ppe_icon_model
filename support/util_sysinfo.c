#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <pwd.h>

#include <sys/types.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#elif HAVE_SYS_UNISTD_H
#include <sys/unistd.h>
#endif

#include <sys/utsname.h>

#ifdef _AIX
#include <sys/systemcfg.h>
#endif

/* funcion implemetations */ 

void util_user_name(char *name, int *actual_len)
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
      strcpy(name, current->pw_gecos);
      strcat(name, " (");
      strcat(name, current->pw_name);
      strcat(name, ")");
    }
  }    

  *actual_len = strlen(name);
}

void util_os_system(char *name, int *actual_len)
{
  struct utsname system_info;
  int iret;
#ifdef _AIX
  int imp;
#endif

  iret = uname(&system_info);

  /* Cray doesn't handle operating system information proper */

#ifdef _UNICOSMP
  strcpy(name, "UNICOS");
#else
  strcpy(name, system_info.sysname);
#endif
  strcat(name, " ");
#ifdef _AIX
  strcat(name, system_info.version);
  strcat(name, ".");
  strcat(name, system_info.release);
#else
  strcat(name, system_info.release);
#endif
  strcat(name, " ");
#ifdef _AIX
  imp = _system_configuration.implementation;
  if (imp == POWER_4)
    strcat(name, "Power4");
  else if (imp == POWER_5)
    strcat(name, "Power5");
  else if (imp == POWER_6)
    strcat(name, "Power6");
  else if (imp == POWER_7)
    strcat(name, "Power7");
  else {
    strcat(name, "unknown");
  }
#else
  strcat(name, system_info.machine);
#endif

  *actual_len = strlen(name);
}

void util_node_name(char *name, int *actual_len)
{
  struct utsname host_info;
  int iret;
  char *fqdn;

  iret = uname(&host_info);
  if (iret != EFAULT) 
    {
      strcpy (name, host_info.nodename);
      
    } 
  else 
    {
      strcpy(name, "unknown");
    }

  /* if FQDN is returned remove domainname */
  if ((fqdn = strchr(name,'.')) != NULL) 
    {
      *fqdn = '\0';
    }
  
  *actual_len = strlen(name);

  return;
}

