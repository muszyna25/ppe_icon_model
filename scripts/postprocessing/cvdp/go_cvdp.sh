#!/bin/bash -x

\rm cvdp/namelist
\rm cvdp/namelist_obs

ystart=2800
yend=2926

for expid in slo1016
do


expdir=/work/mh0287/users/stephan/Icon/Git_repos/icon.oes.mergecpl1710/experiments/${expid}

rm -rf data_${expid}_${ystart}-${yend}
mkdir data_${expid}_${ystart}-${yend}


cdo -O -f nc cat $expdir/*atm_2d_ml*  BOT.nc
REMAP="-remapnn,r360x180"
SELYEAR="-selyear,$ystart/$yend"
cdo $REMAP $SELYEAR  -selvar,tas  BOT.nc data_${expid}_${ystart}-${yend}/tas_${expid}_icon_BOT_mm_${ystart}01_${yend}12.nc
cdo $REMAP $SELYEAR  -selvar,pr    BOT.nc data_${expid}_${ystart}-${yend}/pr_${expid}_icon_BOT_mm_${ystart}01_${yend}12.nc
cdo $REMAP $SELYEAR  -selvar,psl  BOT.nc data_${expid}_${ystart}-${yend}/psl_${expid}_icon_BOT_mm_${ystart}01_${yend}12.nc
cdo $REMAP $SELYEAR  -selvar,ts  BOT.nc data_${expid}_${ystart}-${yend}/ts_${expid}_icon_BOT_mm_${ystart}01_${yend}12.nc
##cdo $REMAP $SELYEAR  -selvar,snd  BOT.nc data_$expid/snd_${expid}_icon_BOT_mm_${ystart}01_${yend}12.nc

\rm BOT.nc

cat >>cvdp/namelist<<EOF
$expid | /work/mh0033/m211054/projects/icon/cvdp/data_${expid}_${ystart}-${yend}/ | ${ystart} | ${yend}
EOF

done


cd cvdp
ncl driver.ncl


