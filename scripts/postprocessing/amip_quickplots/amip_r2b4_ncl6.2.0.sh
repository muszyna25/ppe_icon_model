#!/bin/sh
#
set -e
#####################################################
#
# Please adjust the following variables in the script
# explanation see below
#


TYP=ANN


EXP=mbe0662

YY1=1979
YY2=2008

NAME=ml_${EXP}_${YY1}-${YY2}_${TYP}.nc

COMMENT='amip r2b4'


SINGLE=1
PAGE=0

ATM_3d=0
ATM_2d=1
TAB=0


WORKDIR=/home/zmaw/m214091/contrib_quick/test
MODELDIR=/scratch/local1/m214091/icon-aes-bko
#######################################################
#
QUELLE=${MODELDIR}/scripts/postprocessing/amip_quickplots/r2b4_ncl6.2.0
export QUELLE
echo QUELLE path $QUELLE
echo $QUELLE

# cell=filled triangles cont=filled contour 
# default= 1: model=cell, ERAinterim=cont, model-ERAinterim=cont
default=1
# otherwise default = 0
# if cell=1 all plots are filled triangles
# if cell=0 all plots are filled contour
  cell=1
if [ "$default" = "1" ];then
   cell=-1
fi

date
#
# ERAinterim (1979-2008)
TO=2008
# era_RES= atmospheric grid resolution 63 
#          (used for ERAinterim-data and
era_RES=63
# atm_RES= ICON-grid resolution 
atm_RES=r2b4
# oce_RES= ocean grid resolution GR15 GR30 TP04 TP10
oce_RES=GR15
# Lev 47 possible
LEV=47
FIXCODE=/pool/data/ECHAM6/post/FixCodes


GrdDir=/pool/data/ICON/grids/private/r2b4_amip/
   GrdInfoFile=${GrdDir}/r2b4_amip.nc
export GrdInfoFile
   GrdInfoCellFile=${GrdDir}/remapdis_r2b4_amip_cell_t63grid.nc
export GrdInfoCellFile
   ERAinDir=/pool/data/ICON/post/r2b4_amip/ERA-Interim/
export ERAinDir

  PLTDIR=${WORKDIR}/${EXP}_${TYP}
  export PLTDIR
  if [ ! -d ${PLTDIR} ] ; then
    mkdir ${PLTDIR}
    echo ${PLTDIR}
  fi

cd ${PLTDIR}
pwd



 . /sw/share/Modules/init/ksh 
 module unload ncl
  module load ncl/6.2.0-precompiled 
 
which ncl



#
#--------------- TABLE --------------------------
if [ "$TAB" = "1" ]
then
cp ${FIXCODE}/F${era_RES}${oce_RES}_LAND F_LAND
cp ${FIXCODE}/F${era_RES}_GLACIER F_GLACIER

${QUELLE}/TABLEr2b4_job $TYP $NAME $EXP $YY1 $YY2  $WORKDIR 
fi
#
#--------------- QUICKPLOTS ---------------------
#
if [ "$ATM_2d" = "1" ]
then

MEANTIME="(${YY1}-${YY2})"
cat >var.txt << eof00
$NAME
$TYP
$EXP
$MEANTIME
$COMMENT
$PLTDIR
eof00
#
#

${QUELLE}/PREPAREatm_2d_r2b4 $TYP $NAME $atm_RES $TO $WORKDIR

if [ "$PAGE" = "1" ]
then
  nclsh  ${QUELLE}/atm_2d_r2b4_page.ncl -default=${default} -cell=${cell}

fi
if [ "$SINGLE" = "1" ]
then
  nclsh  ${QUELLE}/atm_2d_r2b4_single.ncl -default=${default} -cell=${cell}
fi

set +e
rm Ubusy_*.nc   var.txt  sea.nc


set -e

echo '####################################################'
echo  you find your plots in
echo ${PLTDIR}
echo '#####################################################'

fi

############################################################
if [ "$ATM_3d" = "1" ]
then
echo ATM_dyn

MEANTIME="(${YY1}-${YY2})"
cat >var.txt << eof00
$NAME
$TYP
$EXP
$MEANTIME
$COMMENT
$PLTDIR
eof00

${QUELLE}/PREPAREatm_3d_logp_r2b4 $TYP $NAME $atm_RES $era_RES $TO $LEV $WORKDIR
if [ "$PAGE" = "1" ]
then
  ncl   ${QUELLE}/atm_3d_logp_page.ncl
fi
if [ "$SINGLE" = "1" ]
then
  ncl   ${QUELLE}/atm_3d_logp_single.ncl
fi
rm -f Ubusy_*.nc  Uatm_dyn_pl_log Uatm_dyn_pl

${QUELLE}/PREPAREatm_3d_r2b4 $TYP $NAME $atm_RES $era_RES $TO $WORKDIR

if [ "$PAGE" = "1" ]
then
  nclsh  ${QUELLE}/atm_3d_linp_page.ncl 
  nclsh  ${QUELLE}/atm_3d_map_page.ncl -default=${default} -cell=${cell}
fi
if [ "$SINGLE" = "1" ]
then
  nclsh  ${QUELLE}/atm_3d_linp_single.ncl 
  nclsh  ${QUELLE}/atm_3d_map_single.ncl -default=${default} -cell=${cell}
fi

rm -f Ubusy_*.nc   var.txt var1.txt Uatm_dyn_pl_log Uatm_dyn_pl


echo '####################################################'
echo  you find your plots in
echo ${PLTDIR}
echo '#####################################################'

fi
exit

#Please adjust the following variables in the script:
#
# EXP= experiment number, appears in the caption of the plots
#
# COMMENT= the comment appears in the subtitle of the plots
#          maximum length 20 characters 
# PRINTER = 1 Name of the black and white printer for the summary table 
#           0 ghostview  the plot
# PRINTERC= 1 Name of the color printer
#           0 ghostview  the plot
# TYP= average to compare with ERAinterim-data(1979-1999)or (1979-2008)
#      ANN(annual), DJF(Dec-Feb), MAM(mar-may)  JJA(jul-aug), SON(sep-nov),
#      JAN ... DEC
#
# YY1= start date, appears in the caption of the plots
# YY2= end date, appears in the caption of the plots
#                                
#      
# NAME= XXX name of data files (maximum length 10 characters)
# WORKDIR= working directory (containing the input data atm_dyn_XXX and atm_phy_XXX)
#
# MODELDIR= model directory
#
#
# ATM_dyn= 1 plot atmosphere data
#          0 no plot of atmospheric data
# ATM_phy= 1 plot surface data
#          0 no plot of surface data
#
#       the plot program expects the following two files:
#                 atm_phy_XXX (surface data, containing at least:
#                           variable: 
#                                     clwvi Liquid water + ice content
#                                     clt   total cloud cover     
#                                     psl   sea level pressure    
#                                     tas   2 m temperature       
#                                     ts    surface temperature   
#                                     tauu  zonal wind stress     
#                                     prw   column water vapor    
#                                           vertical integral of cloud liquid water
#                                     pr    total precipitation   
#                                                          
#       the interpolation from model level to pressure level computes this programm automatically 
# 
#                atm_3d_XXX (atmosphere data, with the following pressure levels 
#                         in hPa:  1000,925,850,775,700,600,500,400,300,250,
#                                   200,150,100,70,50,30,10
#                         containing at least:
#                             variable: ta  temperature           
#                                       ua  zonal wind            
#                                       va  meridional wind       
#                                       hus specific humidity     
#                                           velocity potential
#                                       clw cloud liquid water    
#                                       cli cloud ice             
#                                       zg  geopotential height
#                                       hur relative humidity
#                                       cl  cloud cover
#
#                atm_3d_XXX (atmosphere data, 
#                         e.g. with the following pressure 47 levels in hPa:
#                         100900,99500,97100,93900,90200,86100,81700,77200,
#                         72500,67900,63300,58800,54300,49900,45700,41600,
#                         37700,33900,30402,27015,23833,20867,18116,15578,
#                         13239,11066,9102,7406,5964,4752,3743,2914,2235,1685,
#                         1245,901,637,440,296,193,122,74,43,23,11,4,1
#
#                         
#                         containing at least:
#                             variable: ta  temperature
#                                       ua  zonal wind
#                                       va  meridional wind
#        

