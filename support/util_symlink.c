#define _POSIX_C_SOURCE 200112L 

#include "config.h"

#include <errno.h>

#include <sys/types.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#elif HAVE_SYS_UNISTD_H
#include <sys/unistd.h>
#endif

int util_islink(char *path)
{
#ifdef PATH_MAX
  char buf[PATH_MAX];  
#else
  char buf[1024];
#endif
  ssize_t len;
  
  if ((len = readlink(path, buf, sizeof(buf)-1)) != -1)
    {
      buf[len] = '\0';
      /* file is a link ..., return C true */
      return 1;
    }
  else
    {
      if (errno == EINVAL)
	{
	  /* file is not a link ... */
	  return 0;
	}
    }
  
  /* something else went wrong ... */
  return -1;
}
