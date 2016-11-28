#!/bin/bash
#########################################################
# Part 7 in MJO-Analysis Suite
#   * Loops through all variables
#   * Calls NCL-Skript
#     "mjo_analysis_clivar_rmm-index.ncl"
#     and provides required information on
#     dimensions etc.
#   * NCL Script
#     - Computes and plots EOFs
#-------------------------------------------------------
# DWD, FE 13, Julia Keller, 02/2016
########################################################


for variable in  OLR TOT_PREC U200 U850 V850 
do
  fileinpart=${variable}"_"${dataext}"anom"
  echo "Process "${variable} 
  ncl 'var="'${variable}'"' \
      'infile="'${filepath}${dataset}'_'${fileinpart}'.nc"' \
      'plotdir="'${plotpath}'"' \
      'datainfo="'${dataset}'"' \
      ymdStrt=${StrtDate}   \
      ymdLast=${LastDate}   \
      bpmin=${bandpassmin}  \
      bpmax=${bandpassmax}  \
      bpwgt=${bandpasswgt}  \
      spd=${samperday}      \
      latmin=-30             \
      latmax=30             \
      neof=3                \
      mjo_analysis_clivar_eofs.ncl
done


