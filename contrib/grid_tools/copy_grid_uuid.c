/* -------------------------------------------------------------------
 *
 * ICON_GRID_GET
 *
 * Utility for finding ICON grid files from the ICON XML grid table.
 *
 * @author F. Prill (2013-04-30)
 *
 * -------------------------------------------------------------------
 */
#include <stdio.h>
#include <stdlib.h>
#include "netcdf.h"


/* --------------------------------------------------------------
 * Main program
 * -------------------------------------------------------------- */

int main(int argc, char *argv[])
{
    if (argc < 3) {
	printf("Usage: %s  infile.nc outfile.nc\n", argv[0]);
	return 0;
    }

    int  ncid1, ncid2;         /* netCDF IDs */

    ncid1 = ncopen(argv[1], NC_NOWRITE);
    ncid2 = ncopen(argv[2], NC_WRITE);

    ncredef(ncid2);           /* enter define mode */

    ncattcopy(ncid1, NC_GLOBAL, "uuidOfHGrid",         ncid2, NC_GLOBAL);
    ncattcopy(ncid1, NC_GLOBAL, "number_of_grid_used", ncid2, NC_GLOBAL);
    ncattcopy(ncid1, NC_GLOBAL, "ICON_grid_file_uri",  ncid2, NC_GLOBAL);
    ncattcopy(ncid1, NC_GLOBAL, "centre",              ncid2, NC_GLOBAL);
    ncattcopy(ncid1, NC_GLOBAL, "subcentre",           ncid2, NC_GLOBAL);
    
    ncendef(ncid2);           /* leave define mode */
    ncclose(ncid1);
    ncclose(ncid2);
    return 0;
}
