/* -------------------------------------------------------------------
 *
 * ICON_GRID_GET
 *
 * Utility for finding ICON grid files from the ICON XML grid table.
 *
 * Search procedure:
 *
 * 1. Extract basename from full grid URI and try to locate a file with this
 *    name in the local directory.
 *
 * 2. Similar to 1. but traversing a directory given by the environment
 *    variable "ICON_GRID_DIR".
 *    
 * 3. Try to retrieve from full URI (not yet implemented).
 *
 *
 * @todo Implementation of network-based access to URI (preliminary code
 *       available).
 * 
 * @todo Error handling, i.e. non-existing directory, check if result is
 *       actually a file (not a directory), ...
 *
 * @todo Proper usage information and option parsing.
 *
 *
 * @author L. Kornblueh, F. Prill (2013-04-05)
 *
 * -------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <ctype.h>
#include "xml_scan.h"


const char *ICON_GRID_DIR       = "ICON_GRID_DIR";
const char *ICON_XML_GRID_TABLE = "ICON_XML_GRID_TABLE";



/*
 * Scans given directory for file @p basename.
 */
int scan_dir(const char *directory, const char *basename)
{
    DIR *dp;
    struct dirent *ep;     
  
    dp = opendir (directory);
  
    if (dp != NULL)
    {
	while ((ep = readdir (dp)) != NULL) 
	{
	    if (strcmp(ep->d_name, basename) == 0) return 1;
	}
	(void) closedir (dp);
    }
    else
    {
	perror ("Couldn't open the directory");
    }

    return 0;
}


void trim_printf(size_t len, const char *str)
{
  if(len == 0) return;
  const char *end;
  size_t out_size;
  // Trim leading space
  while(isspace(*str)) str++;
  if(*str == 0)  // All spaces?
  {
    return;
  }
  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;
  end++;
  // Set output size to minimum of trimmed string length and buffer size minus 1
  out_size = (end - str) < len-1 ? (end - str) : len-1;
  printf("%.*s\n", (int) out_size, str);
}


/* --------------------------------------------------------------
 * Main program
 * -------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  int iret = 0;

  if (argc < 5) {
      printf("Usage: %s  number_of_grid_used centre subcentre type\n", argv[0]);
      return 0;
  }

  char *icon_grid_dir       = getenv(ICON_GRID_DIR);
  char *icon_xml_grid_table = getenv(ICON_XML_GRID_TABLE);

  if (icon_xml_grid_table == NULL) {
      printf("Environment variable %s missing!\n", ICON_XML_GRID_TABLE);
      return 0;
  }

  struct t_grid_attr grid_attr = {.number_of_grid_used = atoi(argv[1]), 
				  .centre              = atoi(argv[2]),
				  .subcentre           = atoi(argv[3]),
				  .type                = argv[4]};

  /* read XML file and look for a given grid URI and description. */
  parseXML(icon_xml_grid_table);

  /* extract file name from URI */
  const char *uri         = get_grid_uri(&grid_attr);
  const char *extpar      = get_grid_extpar(&grid_attr);
  const char *description = get_grid_description(&grid_attr) ;

  printf("URI:    %s\n", uri);
  if (description != NULL) {
    printf("ExtPar: %s\n", extpar);
  }
  if (description != NULL) {
    printf("\n        ");
    trim_printf(strlen(description), description);
    printf("\n");
  }

  char *basename = strrchr(uri, '/') + 1;
  printf("search for \"%s\"\n", basename);
  
  /* scan local directory for file */
  iret = scan_dir("./", basename);
  if (iret) {
      printf("File found in\n  ./%s\n", basename);
      return 0;
  }

  /* scan directory given by environment variable */
  if (icon_grid_dir != NULL) 
    {
      iret = scan_dir(icon_grid_dir, basename);
      if (iret) {
	  printf("File found in\n  %s/%s\n", icon_grid_dir, basename);
	  return 0;
      }
    }

  printf("File not found!\n");
  return 0;
}



