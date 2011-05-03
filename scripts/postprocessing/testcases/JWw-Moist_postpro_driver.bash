#!/bin/bash
#
#==========================================================================
#      Driver script for postprocessing and visulization for the 
#           Jablonowski-Williamson baroclinic wave test 
#                      with moist processes
#==========================================================================
#
# History: 
# Initial version for dry test case by Hui Wan (MPI-M, 2009-04)
# Modified for the moist case by Hui Wan (MPI-M, 2010-07)
#
# Short description:
# This script interpolates the orginal ICOHDC output (on hybrid 
# vertical levels) to selected pressure levels, then makes
# contour plots.
#
# Software needed:
# - CDO (Climate Data Operators, www.mpimet.mpg.de/cdo) 
#   for interplation and for the spectral transform;
# - NCL (NCAR Command Language, www.ncl.ucar.edu)
#   for visualization.
#
#==============================================================================
#
#
#==============================================================================

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
# Check if a parameter is given

if [ "x$1" != "x" ]
then
  set_env=$1 
else
  set_env=/null
fi

# 1. About the model output
#--------------------------------------------------------------------------
# 1.0 Which model is it, ECHAM or ICOHAM?

export Model="ICOHAM"     # or: export Model="ECHAM"

# 1.1 The directory in which the model output can be found. 
# Don't forget the trailing "/".

# for automatic testing
dir=`pwd -P`
icon_path=${dir%%scripts/postprocessing/testcases}

# Output frequency expressed as number of time slices per day. 
# (The value you give here must be ASCII representations of integer values, 
# such as "1" for daily output, and "4" for 6-hourly data.)

export output_frequency=4

# 1.2 The data file name is constructed in the same way as in the model.
# The name is composed of experiment name, hor. and vert. resolution (cf. 1.5),
# and a file index (cf. 1.6), if the output has been split on several files.
# The experiment name is the first part of the file name of your model output. 
# You may have used something like "A0001" just for simplicity.

EXP="hat_jww_echam"  
#DOMAIN="DOM01"

# For ECHAM data, specify the year and month that show up in the data file name
#month=197801

# 1.3 The experiment identifier that will appear in the plots, e.g.,
# "brk R2B5L31" 

config_string="JWw + ECHAM physics, R2B4L31" 

# 1.4 Shape of control volume (3 = triangle, 6 = hexagon/pentagon) 

export cell_type=3

# 1.5 Spatial resolution
# (This script assumes that a string indicating the spatial resolution,
#  e.g., "_R2B04L31", has been added after ${EXP} as part of the name
#  of the model output.)

horizontal_resolution="R2B04"
vertical_resolution="L31"

# 1.6 Whether the model output is in a single NetCDF file or 
# has been split into separate files. (Here we assume that
# in case of a single file, the string "_0001" has been added to
# the name of the model output after ${EXP}"_RxBxxLxx".
# In case of output splitting, the first file is labeled "_0001",
# the second one "_0002", and so on.)

data_file_split=1     # 0 = single file; 1 = multiple files
export post_start=1
export post_end=10

# Note that for the 850 hPa contour plots, only the first 10 days 
# of model output are needed. For the kinetic energy spectra
# and for calculation of the difference norms, we need all the 30 days.

#--------------------------------------------------------------------------
# 2. Decide what to do
#--------------------------------------------------------------------------
# It is assumed that the simulation was performed with moist processes,
# (i.e., phyiscs = 1 or 2)

physics=2  #( 0 = no moist processes; 2 = ECHAM physics; 3 = NWP physics)

# this scripts visualizes 
# - (large scale and convective) precipitation rate at the Earth's surface,
#   overlaid with surface pressure;
# - 850 hPa cloud cover overlaid with vertical velocity.
#
# You can let the script produce the whole set by setting the following 
# variable to 1, or switch off the contour plotting completely (e.g., when
# you only want to do vertical interpolation) by setting the 
# variable to zero,

cn_plot_option=2

# or set it to any other integer, and switch on or off each
# variable separately (1=ON,0=OFF):

plot_precip=1
plot_cloud=1

# in which format (pdf, eps or ps) should the plots be drawn?

export plot_file_format="pdf"

# Should the files contain multiple or single panel per page? 
# (Use single for Buildbot!)

export multi_panel=1   # 1=multi,0=single

#--------------------------------------------------------------------------
# Before plotting, vertical interpolation is needed. If you have already
# done this and stored the data, turn this step off.
# Like for plotting, you can choose to interpolate all the variables
# listed above (except for surface pressure) by setting the following 
# variable to 1, to turn off interpolation completely by setting it to 0, 

interp_option=2    # (1=ON,0=OFF)

# or switch on or off individual variables (1=ON,0=OFF):
interp_precip=1    # Needed only for ECHAM. In fact the precipitation
                   # rates are directly available in the model output.
                   # We just need to call afterburner to get the gridpoint
                   # values of surface pressure
interp_cloud=1

# Remove these files after finishing the diagnoses? 
# (1=REMOVE,0=SAVE FOR LATER USE)

rm_tmp_files=1

if [ -f ${set_env} ] 
then 

  echo " "
  echo " !!!!! Use setting from ./${set_env}"
  echo " "
  source ./${set_env}
fi

# Where should the plot files be located? Don't forget the trailing "/".

model_data_path="${icon_path}experiments/$EXP/"

plot_file_path="${model_data_path}plots/"

#--------------------------------------------------------------------------
# Now specify the directory in which the pressure level data and other 
# intermediate files (excluding the remapping weights) should be placed. 
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
echo "=== Postprocessing started for the moist JW wave test."

export script_path=`pwd`"/JWw-Moist_postpro_scripts/"
cd ${script_path}

# The directories for intermediate data and plots will be created, if 
# not already there

if [ ! -d ${plot_file_path} ]; then
   mkdir -p ${plot_file_path} 
fi
if [ ! -d ${tmp_data_path} ]; then
   mkdir -p ${tmp_data_path} 
fi

#==========================================================================
# Perform vertical interpolation and make contour plots
#==========================================================================
if [ ${interp_option} -eq 1 ]; then
   interp_precip=1
   interp_cloud=1
elif [ ${interp_option} -eq 0 ]; then
   interp_precip=0
   interp_cloud=0
fi

if [ ${cn_plot_option} -eq 1 ]; then
   plot_precip=1
   plot_cloud=1
elif [ ${cn_plot_option} -eq 0 ]; then
   plot_precip=0
   plot_cloud=0
fi

# For the NCL plotting scripts

resolution=${horizontal_resolution}${vertical_resolution}
export PlotPath=${plot_file_path}
export Resolution=${resolution}
export ConfigStr=${config_string}
export ExpName=${EXP}_${resolution}

# Temporary variables
if [ $cdo_silence -eq 1 ]; then
   silence='-s'
else
   silence=''
fi

#------------------------------------------------------------------------
# for ECHAM, get grid point values of ps for the precipitation plot
#------------------------------------------------------------------------
if [ $Model == "ECHAM" -a $interp_precip -eq 1 ];then

cat > $Model_$EXP.nml <<EOF
&SELECT
 TYPE = 20, FORMAT = 2,
 CODE = 134,142,143,
&END
134 = aps
142 = aprl
143 = aprc
EOF

   echo "=== Getting ECHAM data in grid point space ..."
   fid=`expr $post_start + 1`
   while [ $fid -le `expr $post_end + 1 ` ]; do

     tid2=$(printf "%02d" $fid)
     tid4=$(printf "%04d" $fid)
     fostem=${tmp_data_path}${EXP}_${resolution}_${tid4}

     fecham=${model_data_path}${EXP}_${month}.${tid2}"_echam.nc"
     fetmp=${fostem}"_etmp.nc"
     after $fecham $fetmp < $Model_$EXP.nml 

     fsurf=${model_data_path}${EXP}_${month}.${tid2}"_surf.nc"
     fstmp=${fostem}"_stmp.nc"
     cdo $silence selvar,jrsfl,jssfl,jrsfc,jssfc $fsurf $fstmp
     cdo $silence merge $fetmp $fstmp  $fostem".nc"

     rm $fetmp $fstmp
     echo " file "$tid4" done"
     fid=`expr $fid + 1`
   done
   echo "=== Done."
fi
#------------------------------------------------------------------------
# plot precipitation at selected time steps (see the ncl script)
#------------------------------------------------------------------------     
if [ $plot_precip -eq 1 ]; then

   if [ $Model == "ECHAM" ]; then
       export DataPath=${tmp_data_path}
   else
       export DataPath=${model_data_path}
   fi
   echo
   echo "=== Plotting precipitation ..."
   ncl JWw-Moist_plot_PS-Precip.ncl 
   check_error $? "In script JWw-Moist_postpro_driver.bash: part 'Plotting precipitation'"
   echo "=== Done."
fi
#------------------------------------------------------------------------
# interpolate and plot vertical velocity, cloud cover and hydrometeor 
#------------------------------------------------------------------------
if [ $interp_cloud -eq 1 ];then

   if [ $Model == "ECHAM" ];then

      cat > $Model_$EXP.nml <<EOF
      &SELECT
       TYPE = 30, FORMAT = 2,
       CODE = 134,135,162,153,154, 
       LEVEL = 90000,
      &END
      LEVEL = 90000,85000,80000,70000
      134 = aps
      135 = omega
      162 = aclc
      153 = xl
      154 = xi
EOF

      echo "=== Getting ECHAM data on pressure levels ..."
      fid=`expr $post_start + 1`
      while [ $fid -le `expr $post_end + 1 ` ]; do
   
        tid2=$(printf "%02d" $fid)
        tid4=$(printf "%04d" $fid)
        fecham=${model_data_path}${EXP}_${month}.${tid2}"_echam.nc"
        fostem=${tmp_data_path}${EXP}_${resolution}_${tid4}
        after $fecham ${fostem}"_cloud_p.nc" < $Model_$EXP.nml 
   
        echo " file "$tid4" done"
        fid=`expr $fid + 1`
      done
      echo "=== Done."

   else #---- ICOHAM -----
      echo
      echo "=== Interpolating omega and cloud cover ..."
      fid=$post_start
      while [ $fid -le $post_end ]; do
   
        tid4=$(printf "%04d" $fid)
        fori=${model_data_path}${EXP}${DOMAIN}_${resolution}_${tid4}".nc"
        ftmp=${tmp_data_path}${EXP}_${resolution}_${tid4}"_cloud_p.nc"
       #cdo $silence ml2pl,90000,85000,80000,70000 \
        cdo $silence ml2pl,90000 \
                     -selname,OMEGA_PHY,ACLC,Qw,Qi,PS    ${fori} ${ftmp}
        check_error $? "In script JWw-Moist_postpro_driver.bash: part 'Interpolating cloud over'"
        echo " File $tid4 done"
   
        fid=`expr $fid + 1`
      done
      echo "=== Done."
   fi
fi

if [ $plot_cloud -eq 1 ]; then

   export DataPath=${tmp_data_path}
   echo
   echo "=== Plotting omega and cloud cover ..."
   ncl JWw-Moist_plot_omega-cloud.ncl 
   check_error $? "In script JWw-Moist_postpro_driver.bash: part 'Plotting omega and cloud cover'"
   echo "=== Done."
fi

#------------------------------------------------------------------------
# Clean up
#------------------------------------------------------------------------
echo 
echo "=== Postprocessing finished for the moist JW wave test."
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
