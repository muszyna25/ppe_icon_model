HPC_GNU_NETCDF_INCLUDE = -I/e/uhome/hanlauf/opt/gfortran47/include/
HPC_GNU_LIBS           = -L/e/uhome/hanlauf/opt/gfortran47/lib -lnetcdf

INCLUDE                = ${HPC_GNU_NETCDF_INCLUDE}
LIBS                   = ${HPC_GNU_LIBS} 


all:
# build "icon_grid_get" utility
	ragel xml_scan.rl -o xml_scan.h
	gcc icon_grid_get.c -o icon_grid_get -std=gnu99 -O0 -g -Wall
# build "copy_grid_uuid" utility
	gcc ${INCLUDE} copy_grid_uuid.c -o copy_grid_uuid -std=gnu99 -O0 -g -Wall ${LIBS}