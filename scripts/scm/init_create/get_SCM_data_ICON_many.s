#!/bin/bash
# ----------------------------------------------------------------------------
# prepare SCM input data for multiple points
#
# attention: lon between -180 and 180!!
#
# work flow SCM from ICON input:
#  - ICON ini:   read_icon_ana_oper_mem1_40km.s
#  - ICON run:   run_ICON_4_SCMini
#  - SCM ini:    get_SCM_data_ICON.py
#  - SCM extpar: create_SCM_extpar_ICON.py
#  - SCM run:    run_SCM_ICONini
#  - plot SCM:   plot-scm-*.py
#
# 2021-02-16 Martin Koehler
# ----------------------------------------------------------------------------

for ((ii=1; ii<=17; ii+=1)) ; do

  lat_scm='-60'
  lon_scm=$((ii*20-180))
  echo '--- processing lat/lon: ' $lat_scm $lon_scm ' ---'

  python3 get_SCM_data_ICON.py $lat_scm $lon_scm

done
