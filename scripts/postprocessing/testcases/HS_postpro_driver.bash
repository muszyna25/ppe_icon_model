#!/bin/bash
#
#==========================================================================
#      Driver script for postprocessing and visulization for the 
#                        Held-Suarez test 
#==========================================================================
#
# History: 
# Initial version by Hui Wan (MPI-M, 2010-07)
# Pilar Ripodas (DWD, 2010-10), options for the Non Hydrostatic model
#
# Short description:
# This script 
#
# Software needed:
# - CDO (Climate Data Operators, www.mpimet.mpg.de/cdo) 
#   for interplation and for the spectral transform;
# - NCL (NCAR Command Language, www.ncl.ucar.edu)
#   for visualization.
#
#==============================================================================
#
check_error()
{

# Check if the first parameter (return status) is not OK 

  if [ $1 -ne 0 ] 
  then

# Stop running this script and return the error status
    echo "ERROR: $2"
    exit $1
  fi
}

#==========================================================================
#                          USER'S SPECIFICATIONS 
#--------------------------------------------------------------------------
# 0. About the dynamical core ( Hydrostatic or Non-Hydrostatic )
#--------------------------------------------------------------------------
# The variable Model can be ICOHDC (Hydrostatic) or ICONAM (Non-Hydrostatic)
# If set to ICOHDC the evolution of zonal temperature variance is 
# calculated and ploted at 750 hPa and the vertical variable in the plots 
# of the zonal mean climate is eta
# If set to ICONAM, the evolution of zonal temperature variance is 
# calculated and ploted at model level 24 and the vertical variable in the plots 
# of the zonal mean climate is the height
# 

Model="ICOHDC"
#Model="ICONAM"

#--------------------------------------------------------------------------
# 1. About the model output
#--------------------------------------------------------------------------
# 1.1 What is the experiment ID? 
# The experiment ID is the first part of the file name of your model output, 
# which is followed by information on hor. and vert. resolution (cf. 1.5),
# and a file index (cf. 1.6) if the output has been split on several files.
# You may have used something like "A0001" just for simplicity.

EXP="hat_heldsuarez"  
#EXP="nh35_heldsuarez"
#EXP="nh35_heldsuarez_div1_open_ubc_vwind_off_0p15_ex_expol_0p66"
#EXP="nh35_hex_heldsuarezB3"

# 1.2 In which directory are the data files located? (Don't forget the trailing "/"!)

model_data_path="somewhere/"
#model_data_path='/scratch/work/mh0287/hwan/exps/icon-dev/output_201007_port_phy/hat_hs_2tl/'
model_data_path="/e/uscratch/mripodas/icon-dev/experiments/${EXP}/"

# The ICON model developers have an automatic testing system which uses
# the following convention:
 
#dir=`pwd -P`
#icon_path=${dir%%scripts/ncl_scripts}
#model_data_path="${icon_path}experiments/$EXP/"

# 1.3 The experiment identifier that will appear in the plots, e.g.,
# "brk R2B5L31" 

config_string="HS R2B04L35 spr0.90"

# 1.4 Shape of control volume (3 = triangle, 6 = hexagon/pentagon) 

export cell_type=3

# 1.5 Spatial resolution
# (This script assumes that a string indicating the spatial resolution,
#  e.g., "_R2B04L31", has been added after ${EXP} as part of the name
#  of the model output.)

horizontal_resolution="R2B04"
vertical_resolution="L31"        #ICOHDC
#vertical_resolution="L35"        #ICONAM  

#--------------------------------------------------------------------------
# 2. Decide what to do
#--------------------------------------------------------------------------
# A complete set of the contour plots includes
#
# a. evolution of zonal temperature variance at 750 hPa;
# b. the zonal mean "climate" computed for the user-specified period, including
#    > time and zonal mean temperature, u-wind, and v-wind,
#    > meridional flux of heat and zonal momentum caused by eddy activities,
#    > eddy kinetic energy.
# c. temporal evolution of the surface pressure glomal mean
#
# If you want to have data or plot of (a), set the following variable to 1

diag_evolution=1

# If you want to have data or plot of (b), set the following variable to 1

diag_climate=1

# If you want to have data or plot of (c), set the following variable to 1

diag_evolution_ps=1

# Since the Held-Suarez test is an idealized climate simulation,
# this script assumes that the output has been split into a series of files. 
# The first file is named ${EXP}"_RxBxxLxx_0001",the second one "_0002",
# and so on. Now specify the starting and ending file indices for plot (a)
# (see above)

evol_istart=1
evol_iend=30 

# Specify the starting and ending file indices for plot (b):

clim_istart=21
clim_iend=30

# Specify the starting and ending file indices for plot (c):

evol_ps_istart=1
evol_ps_iend=30

# If you want to obtain plots, set the following variable to 1 

make_plot=1

# What format do you want, pdf, eps or ps?

export plot_file_format="ps"

# Before plotting, interpolation and statistics computation is needed. 
# If you have already done that and stored the data, turn this step off.

do_computation=1    # (1=ON,0=OFF)

# To diagnose the zonal mean "climate", the ICON model output 
# on the geodesic grid needs to be interpolated to a Gaussian grid. 
# Choose the Gaussian grid by specifying the triangular truncation:

export trunc=85

# For the horizontal interpolation from the ICON grid to Gaussian grid,
# a large part of the time will be spent on calculating the remapping 
# weights. This may take very long when the resolution is high.
# Therefore we suggest calculating the weights only once and store them 
# for later use. 

compute_remap_weights=0   # (1=ON,0=OFF)

# The remapping weights file generated by this script will be named, 
# e.g., icon_R2B04_spr0.90_tri_cell_to_T159.nc 
# If the weights are already available, specify the location:
# (Don't forget the trailing "/")

#remap_weights_path="${model_data_path}remap_weights/"
remap_weights_path='/e/uwork/mripodas/weights/'

# For the remapping, we also need to know which optimization was used 
# to generate the ICON grid. For example,
#   "ori"     : original icosahedral grid without optimization
#   "hro"     : Heikes-Randall optimization
#   "spr0.90" : spring dynamics, with spring coefficient 0.90

grid_optimization="spr0.90"

# Remove these files after finishing the diagnoses? 
# (1=REMOVE,0=SAVE FOR LATER USE)

rm_tmp_files=0

# Where should the plot files be located? Don't forget the trailing "/".

plot_file_path="${model_data_path}plots/"

#--------------------------------------------------------------------------
# Now specify the directory in which the intermediate files 
# (excluding the remapping weights) should be placed. 
# Don't forget the trailing "/".

tmp_data_path="${model_data_path}tmp/"

#--------------------------------------------------------------------------
# Do you want CDO to run in silence mode, or to report everything
# it is doing?

cdo_silence=1   #( 1 = silence mode; 0 = detailed report )

#--------------------------------------------------------------------------
#                    END OF USER'S SPECIFICATIONS 
#==========================================================================

echo
echo "**********************************************************"
echo "***       ICON Tool Kit for Idealized Test Cases       ***"
echo "**********************************************************"
echo 
echo "=== Postprocessing started for the Held-Suarez test."

# Temporary variables

#export script_path=`pwd`"/HS_postpro_scripts/"
#export script_path='/e/uscratch/mripodas/icon-dev/scripts/ncl_scripts/'
export script_path='/e/uhome/mripodas/icon-dev/scripts/ncl_scripts/'
#
#directory with the colomap files
export NCARG_COLORMAPS=${script_path}color_map:${NCARG_COLORMAPS}
#
#echo $script_path
#echo "\n"
resolution=${horizontal_resolution}${vertical_resolution}
fori=${model_data_path}${EXP}_${resolution}
ftmp=${tmp_data_path}${EXP}_${resolution}

if [ $cdo_silence -eq 1 ]; then
   silence='-s'
else
   silence=''
fi

# The directories for intermediate data and plots will be created, if 
# not already there

if [ ! -d ${plot_file_path} ]; then
   mkdir -p ${plot_file_path} 
fi
if [ ! -d ${tmp_data_path} ]; then
   mkdir -p ${tmp_data_path} 
fi

#==========================================================================
# Prepare remappping weights
#==========================================================================
weights=${remap_weights_path}"icon_"${horizontal_resolution}

if [ ${cell_type} -eq 3 ]; then
   weights=${weights}"_"${grid_optimization}"_cell_to_T"${trunc}".nc"
elif [ ${cell_type} -eq 6 ]; then
   weights=${weights}"_"${grid_optimization}"_vert_to_T"${trunc}".nc"
else
   echo "Wrong choice of cell_type. Should be 3 or 6 !"
   exit 1
fi

# compute weights

if [ ${compute_remap_weights} -eq 1 ]; then

 if [ ! -d ${remap_weights_path} ]; then
    mkdir -p ${remap_weights_path} 
 fi

 echo
 echo "=== Computing remapping weights (ICON to Gaussian) ..."

 label=$(printf "%04d" ${evol_istart})
 cdo $silence gendis,t${trunc}grid -seltimestep,1 -selname,PS ${fori}"_"${label}".nc" ${weights}

 check_error $? "In script HS_postpro_driver.bash: part 'Computing remapping weights (ICON to Gaussian)'"
 echo "=== Done."
fi


#====================================================================================================
# Compute the zonal temperature variance on 750 hPa (ICOHDC) or at model level 24 (ICONAM); make plot.
#====================================================================================================
if [ ${diag_evolution} -eq 1 ]; then

      if [ $Model == "ICOHDC" ]; then
       templ="temp750"
      fi
      if [ $Model == "ICONAM" ]; then
       templ="templ24"
      fi

   if [ ${do_computation} -eq 1 ]; then

      #-------------------------
      # vertical interpolation
   
      echo
      if [ $Model == "ICOHDC" ]; then
       echo "=== Interpolating model temperature to 750 hPa ..."
       level="750 hPa"
      fi
      if [ $Model == "ICONAM" ]; then
       echo "=== Selecting temperature at model level 24 ..."
       level="level 24"
      fi
   
      ifile=$evol_istart
      while [ $ifile -le $evol_iend ] ; do
   
        label=$(printf "%04d" $ifile)
        if [ $Model == "ICOHDC" ]; then
         cdo $silence ml2pl,75000 -selname,T,PS,PHIS ${fori}"_"$label".nc" ${ftmp}"_"$label"."${templ}".nc"
         check_error $? "In script HS_postpro_driver.bash: part 'Interpolating model temperature to 750 hPa', file $label"
        fi
        if [ $Model == "ICONAM" ]; then
         cdo $silence sellevel,24 -selname,T ${fori}"_"$label".nc" ${ftmp}"_"$label"."${templ}".nc"
         check_error $? "In script HS_postpro_driver.bash: part 'Selecting model temperature at model level 24', file $label"
        fi
        echo "  File" $ifile "done"
        ifile=` expr $ifile + 1 `
      done
      echo "=== Done."
   
      #--------------------------
      # merge to a single file 
   
      cdo $silence -r copy ${ftmp}"_"[0-9][0-9][0-9][0-9]"."${templ}".nc" \
                           ${ftmp}"_"${evol_istart}-${evol_iend}"."${templ}".nc"
   
      check_error $? "In script HS_postpro_driver.bash: file merging"
   
      rm ${ftmp}"_"[0-9][0-9][0-9][0-9]"."${templ}".nc"
   
      #-----------------------------
      # interpolation to Gaus grid
   
      echo
      echo "=== Computing zonal variance of temperature at "${level}"..."
   
      cdo $silence -r sqr -zonstd -remap,t${trunc}grid,${weights} \
                                         ${ftmp}"_"${evol_istart}-${evol_iend}"."${templ}".nc" \
                                         ${ftmp}"_"${evol_istart}-${evol_iend}"."${templ}".T"${trunc}".zonvar.nc"
   
      check_error $? "In script HS_postpro_driver.bash: part 'Computing zonal variance of temperature at '"${level}
      echo "=== Done."
   fi

   if [ ${make_plot} -eq 1 ]; then

      export Model
      export DataPath=${tmp_data_path}
      export DataID=${EXP}_${resolution}"_"${evol_istart}-${evol_iend}"."${templ}".T"${trunc}".zonvar"
      export PlotPath=${plot_file_path}
      export ConfigStr=${config_string}
      export VarName="T"
      export Resolution=${horizontal_resolution}${vertical_resolution}

      echo "=== Plotting temperature variance ..."
      ncl ${script_path}HS_plot_T-var_time_series.ncl 
      check_error $? "In scripte HS_postpro_driver.bash: part 'Plotting temperature variance'"
      echo "=== Done."
   fi
fi

#========================================================================
# Diagnose the zonal mean climate; make plot
#========================================================================
if [ $diag_climate -eq 1 ]; then

   if [ $do_computation -eq 1 ]; then

      echo
      echo "=== Computing statistics ..."

      #------------------
      # Select variables 
      #------------------
      ifile=$clim_istart
      while [ $ifile -le $clim_iend ] ; do

        label=$(printf "%04d" $ifile)
        cdo $silence selname,T ${fori}"_"$label".nc" ${ftmp}"_"$label".T.nc"
        cdo $silence selname,U ${fori}"_"$label".nc" ${ftmp}"_"$label".U.nc"
        cdo $silence selname,V ${fori}"_"$label".nc" ${ftmp}"_"$label".V.nc"

        check_error $? "In script HS_postpro_driver.bash: part 'selname', file $label"
        echo "  File" $ifile "done"
        ifile=` expr $ifile + 1 `
      done

      #-----------------------------------------
      # Merge all time steps into a single file
      #-----------------------------------------
      for varname in T U V ; do
          cdo $silence -r copy  ${ftmp}"_"[0-9][0-9][0-9][0-9]"."$varname".nc" \
                                ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname".nc"

          rm ${ftmp}"_"[0-9][0-9][0-9][0-9]"."$varname".nc"
          check_error $? "In script HS_postpro_driver.bash: file merging"
      done

      #-----------------------------------------
      # Compute mean and variances
      #-----------------------------------------
      for varname in T U V ; do
       for opr in timavg timvar ; do
           cdo $silence $opr ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname".nc"     \
                             ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname"."$opr".nc"
       done
      done
      check_error $? "In script HS_postpro_driver.bash: part 'mean and variance'"

      #-----------------------------------------
      # Eddy kinetic energy
      #-----------------------------------------

      cdo $silence mulc,0.5 -add ${ftmp}"_"$clim_istart"-"$clim_iend".U.timvar.nc" \
                                 ${ftmp}"_"$clim_istart"-"$clim_iend".V.timvar.nc" \
                                 ${ftmp}"_"$clim_istart"-"$clim_iend".eddy-KE.nc"

      rm ${ftmp}"_"$clim_istart"-"$clim_iend".U.timvar.nc"
      rm ${ftmp}"_"$clim_istart"-"$clim_iend".V.timvar.nc"

      check_error $? "In script HS_postpro_driver.bash: eddy kinetic energy"

      #-----------------------------------------
      # Eddy momentum/heat flux
      #-----------------------------------------

      for varname in T U ; do

          cdo $silence timavg -mul ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname".nc" \
                                   ${ftmp}"_"$clim_istart"-"$clim_iend".V.nc"          \
                                   ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname"V.timavg.nc"

          cdo $silence mul ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname".timavg.nc" \
                           ${ftmp}"_"$clim_istart"-"$clim_iend".V.timavg.nc"          \
                           ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname"aVa.nc"

          cdo $silence sub ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname"V.timavg.nc" \
                           ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname"aVa.nc"      \
                           ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname"pVp.nc"

          rm ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname"V.timavg.nc"
          rm ${ftmp}"_"$clim_istart"-"$clim_iend"."$varname"aVa.nc"

          check_error $? "In script HS_postpro_driver.bash: eddy flux of $varname"
      done

      
      #-------------------------------------------------
      # Interpolate to Gaus grid and compute zonal mean
      #-------------------------------------------------

      varlabel1=("U.timavg" "V.timavg" "T.timavg" "UpVp" "TpVp" "T.timvar" "eddy-KE")
      varlabel2=("U"        "V"        "T"        "UpVp" "TpVp" "TpTp"     "eddy-KE")

      nvar=${#varlabel1[*]} 
      ivar=0

      while [ $ivar -lt $nvar ]; do

        cdo $silence zonavg -remap,t${trunc}grid,${weights}    \
                                  ${ftmp}"_"$clim_istart"-"$clim_iend"."${varlabel1[$ivar]}".nc" \
                                  ${ftmp}"_"$clim_istart"-"$clim_iend"."${varlabel2[$ivar]}".zmta.T"$trunc".nc"

        check_error $? "In script HS_postpro_driver.bash: remap and zonal mean"
        rm ${ftmp}"_"$clim_istart"-"$clim_iend"."${varlabel1[$ivar]}".nc"


      ivar=` expr $ivar + 1 `
      done
      echo "=== Done."
   fi

   if [ $make_plot -eq 1 ]; then

      export Model
      export DataPath=${tmp_data_path}
      export DataID=${EXP}_${resolution}"_"${clim_istart}-${clim_iend}
      export PlotPath=${plot_file_path}
      export ConfigStr=${config_string}
      export Resolution=${horizontal_resolution}${vertical_resolution}
      export trunc

      ncl ${script_path}HS_plot_climate_colored.ncl
   fi
fi

#====================================================================================================
# Compute surface pressure global mean ; make plot.
#====================================================================================================
if [ ${diag_evolution_ps} -eq 1 ]; then

   if [ ${do_computation} -eq 1 ]; then

      #----------------------------------
      # select PS and calculate global mean

      ifile=$evol_ps_istart
      while [ $ifile -le $evol_ps_iend ] ; do

        label=$(printf "%04d" $ifile)

        cdo $silence fldmean  -selname,PS ${fori}"_"$label".nc" ${ftmp}"_"$label".PSmean.nc"       

        echo "  File" $ifile "done"
        ifile=` expr $ifile + 1 `

      done
      echo "=== Done."
   
      #--------------------------
      # merge to a single file 

      cdo $silence -r copy ${ftmp}"_"[0-9][0-9][0-9][0-9]".PSmean.nc" \
                           ${ftmp}"_"${evol_ps_istart}-${evol_ps_iend}".PSmean.nc"
   
      check_error $? "In script HS_postpro_driver.bash: file merging"
   
      rm ${ftmp}"_"[0-9][0-9][0-9][0-9]".PSmean.nc"
   fi

   if [ ${make_plot} -eq 1 ]; then

      export Model
      export DataPath=${tmp_data_path}
      export DataID=${EXP}_${resolution}"_"${evol_ps_istart}-${evol_ps_iend}".PSmean"
      export PlotPath=${plot_file_path}
      export ConfigStr=${config_string}
      export VarName="PS"
      export Resolution=${horizontal_resolution}${vertical_resolution}
      export EXP

      echo "=== Plotting Surface Pressure evolution ..."
      ncl ${script_path}HS_plot_PS-mean_time_series.ncl 
      check_error $? "In scripte HS_postpro_driver.bash: part 'Plotting Surface Pressure global mean'"
      echo "=== Done."
   fi

fi

#------------------------------------------------------------------------
# Clean up
#------------------------------------------------------------------------

echo 
echo "=== Postprocessing finished for the Held-Suarez test."
echo "=== The plots can be found in "${plot_file_path}

if [ $rm_tmp_files -eq 1 ]; then

   rm ${tmp_data_path}/${EXP}_${resolution}*.nc 

   if [ `ls ${tmp_data_path} |wc -l` -eq 0 ]; then
      rm -rf ${tmp_data_path}
   fi 
   echo "=== The temporary data have been removed."

else
   echo "=== The pressure level data and spectral coefficients can be found in "${tmp_data_path}
fi
