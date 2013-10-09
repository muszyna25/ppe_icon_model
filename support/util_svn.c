#include <stdio.h>

/* function implementations */ 

/* print expanded compiler macro __SVN_VERSION containing the current
   source code revision number. */
void print_svn_version()
{
#ifdef __SVN_VERSION
#define xstr(s) str(s)
#define str(s) #s
  fprintf(stderr, "\n ----------------------\n", xstr(__SVN_VERSION));
  fprintf(stderr, " ICON revision #%s\n", xstr(__SVN_VERSION));
  fprintf(stderr, " ----------------------\n\n", xstr(__SVN_VERSION));
#endif
}

