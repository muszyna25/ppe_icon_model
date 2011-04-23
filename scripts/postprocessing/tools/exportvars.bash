#! /bin/bash
#------------------------------------
# Marco Giorgetta (MPI-M, 2009-04-03)
#------------------------------------
# This script defines environment variables that are used
# by NCL plotting scripts, e.g. by "plotvar.ncl".
#
# Run this script by ". exportvar.sh" before starting NCL
# scripts reading these parameters.
#
# (1) ICON directory
#------------------------------------------------------------------
export icon_dir=~/svn/icon/branches/icon-1.0.6_RESTRUCT

# (2) Variable
#------------------------------------------------------------------
export icon_varname=P     # variable {VOR,DIV,U,V,OMEGA,P,PS,...}

# (3.1) Source 1 used by plotvar.ncl and plotdiff.ncl
#------------------------------------------------------------------
export icon_expdir_1=${icon_dir}/experiments/IIIEEEETTTT
export icon_expfile_1=${icon_expdir_1}/IIIEEEETTTT_R2B04L31_0001.nc
export icon_istep_1=20    # in [0,nsteps] 0=initial state
export icon_ilev_1=18     # in [1,nlev]   1=uppermost level
#
# (3.2) Source 2 used by plotdiff.ncl for difference 2-1
#------------------------------------------------------------------
export icon_expdir_2=${icon_dir}/experiments/IIIEEEETTTT
export icon_expfile_2=${icon_expdir_2}/IIIEEEETTTT_R2B04L31_0001.nc
export icon_istep_2=24    # in [0,nsteps] 0=initial state
export icon_ilev_2=18     # in [1,nlev]   1=uppermost level
#
# (4) NCL plots
#------------------------------------------------------------------
# Output type ("NCL workstation")
export icon_nclwks="X11"  # "X11" -> graphics window
                          # "ps"  -> ps  file in $plot_path
                          # "eps" -> eps file in $plot_path
                          # "pdf" -> pdf file in $plot_path
#
# directory for plot files
export icon_plotdir=${icon_expdir_1}
#
# plot file name used by plotvar.ncl : VAR_IT1_IL1
export icon_plotvarname=${icon_varname}_IT${icon_istep_1}_IL${icon_ilev_1}
export icon_plotvarfile=$icon_plotdir/${icon_plotvarname}
#
# plot file name used by plotdiff.ncl: VAR_IT2_IL2-IT1_IL1
export icon_plotdiffname=${icon_varname}_IT${icon_istep_2}_IL${icon_ilev_2}-IT${icon_istep_1}_IL${icon_ilev_1}
export icon_plotdifffile=$icon_plotdir/${icon_plotdiffname}
#
