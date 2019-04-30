# This file contains dependencies that can not be detected automatically.
# The file is not included by the Makefile.main but is only parsed
# with './utils/mkhelper/deplist.py'.

# Fortran to C dependencies:
src/io/restart/mo_c_restart_util.o: support/util_multifile_restart.o
src/io/shared/mo_util_file.o: support/util_file.o
src/io/shared/mo_util_nml.o: support/nml_annotate.o
src/shared/mo_expression.o: support/util_arithmetic_expr.o
src/shared/mo_util_backtrace.o: support/util_backtrace.o
src/shared/mo_util_hash.o: support/util_hash.o
src/shared/mo_util_stride.o: support/util_stride.o
src/shared/mo_util_string_parse.o: support/util_string_parse.o
src/shared/mo_util_sysinfo.o: support/util_sysinfo.o
src/shared/mo_util_system.o: support/util_system.o
src/shared/mo_util_timer.o: support/util_timer.o
src/shared/mo_util_uuid.o: support/util_uuid.o
src/shared/mo_util_vcs.o: version.o

# Undetectable Fortran dependencies:
src/atm_chem_cariolle/lcariolle_do3dt.o:      \
  src/atm_chem_cariolle/lcariolle_o3_column.o

src/atm_phy_echam/mo_echam_phy_init.o:          \
  src/atm_chem_cariolle/lcariolle_init.o         \
  src/atm_chem_cariolle/lcariolle_init_o3.o      \
  src/atm_chem_cariolle/lcariolle_lat_intp_li.o  \
  src/atm_chem_cariolle/lcariolle_pres_intp_li.o

src/atm_phy_echam/mo_interface_echam_car.o:     \
  src/atm_chem_cariolle/lcariolle_do3dt.o        \
  src/atm_chem_cariolle/lcariolle_lat_intp_li.o  \
  src/atm_chem_cariolle/lcariolle_pres_intp_li.o

src/hamocc/icon_specific/bgc_icon.o: \
  src/hamocc/common/chemcon.o         \
  src/hamocc/common/ocprod.o          \
  src/hamocc/common/sedshi.o          \
  src/hamocc/common/swr_absorption.o

src/ocean/drivers/mo_hydro_ocean_run.o:  \
  src/hamocc/icon_specific/bgc_icon.o     \
  src/hamocc/icon_specific/ini_bgc_icon.o

