#!/bin/ksh
#
# Function for checking if the return status is OK 

check_error()
{
  if [ $1 -eq 0 ]; then  # No problem

    echo $2: done.

  else # Stop running this script and return the error status

    echo $2: ERROR ENCOUNTERED!
    exit $1
  fi
}

# argument list processing

while getopts "v:e:" option
do
  case $option in
  v ) varname=$OPTARG ;;
  e ) set_env=$OPTARG ;;
  esac
done

#==========================================================================
#                          USER'S SPECIFICATIONS 
#--------------------------------------------------------------------------
# 0. Model name 
#--------------------------------------------------------------------------

Model="ICOHAM"   # or "ECHAM", "ICONAM"

#--------------------------------------------------------------------------
# 1. About the model output
#--------------------------------------------------------------------------
# 1.1 What is the experiment ID? 
# The experiment ID is the first part of the file name of your model output, 
# which is followed by information on hor. and vert. resolution (cf. 1.5),
# and a file index (cf. 1.6) if the output has been split on several files.
# You may have used something like "A0001" just for simplicity.

EXP="hat_ape_qobs_moist_RH00-clip"

# 1.2 In which directory are the data files located? (Don't forget the trailing "/"!)

model_data_path="somewhere/"

# The ICON model developers have an automatic testing system which uses
# the following convention:
 
dir=`pwd -P`
icon_path=${dir%%scripts/postprocessing/testcases}
model_data_path="${icon_path}experiments/$EXP/"

# 1.3 Shape of control volume (3 = triangle, 6 = hexagon/pentagon) 

cell_type=3

# 1.4 Spatial resolution
# (This script assumes that a string indicating the spatial resolution,
#  e.g., "_R2B04L31", has been added after ${EXP} as part of the name
#  of the model output.)

horizontal_resolution="R2B04"
vertical_resolution="L31"

# 1.5 The experiment identifier that will appear in the plots, e.g.,
# "brk R2B5L31" 

export exp_config="$Model $horizontal_resolution$vertical_resolution"

#--------------------------------------------------------------------------
# 2. Decide what to do
#--------------------------------------------------------------------------
# This script assumes that the output has been split into a series of files. 
# The first file is named ${EXP}"_RxBxxLxx_0001",the second one "_0002",
# and so on. Now specify the starting and ending file indices for the 
# evolution plot.

evol_istart=1
evol_iend=10

# If you want to obtain plots, set the next variable to 1 
make_plot=1

# What format do you want, pdf, eps or ps?
export plot_file_format="pdf"

# Set the orientation to "landscape" if you explicitly set it so.
# (No specification means NCL default - "portrait".)
#export wkOrientation="landscape"

# Before plotting, some data processing is needed. 
# If you have already done that and stored the data, turn this step off.

do_computation=1   # (1=ON,0=OFF)

# Select vertical levels: 
#   - plev_evol for pressure levels (unit: Pa)
#   - hlev_evol for height levels   (unit: km)
#   - mlev_evol for model levels (i.e., level indices 1, 2, etc.)
# If the user assigns a meaningful (i.e., not null) string to 
# one of the three variables, the corresponding vertical levels
# will be selected; If the user does not set any of the three variables
# or sets all of them to null (""), the default pressure levels
# (85000,50000,20000,10000 Pa) will be selected; If more than one 
# of the three variables are assigned meaningful strings, the order
# of precedence is: plev_evol > hlev_evol > mlev_evol.

plev_evol=""

# To diagnose the zonal mean circulation, the ICON model output 
# on the geodesic grid needs to be interpolated to a Gaussian grid. 
# Choose the Gaussian grid by specifying the triangular truncation:

trunc=63

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

remap_weights_path="${model_data_path}remap_weights/"

# For the remapping, we also need to know which optimization was used 
# to generate the ICON grid. For example,
#   "ori"     : original icosahedral grid without optimization
#   "hro"     : Heikes-Randall optimization
#   "spr0.90" : spring dynamics, with spring coefficient 0.90

grid_optimization="spr0.90"

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

if [ -f ${set_env} ] ; then
  source ./${set_env}
fi

# Check user's specification of vertical levels

if [ ${plev_evol:-0} -eq 0 ]; then   # User didn't specify any pressure level
 if [ ${hlev_evol:-0} -eq 0 ]; then  # User didn't specify any height level
  if [ ${mlev_evol:-0} -eq 0 ]; then # User didn't specify any model level
   
   # User didn't specify any kind of levels. Will use default pressure values.
   levtype="p"
   cmd="ml2pl"
   levs="85000,50000,20000,10000"    # default pressure levels

  else # Will select model levels
   levtype="m"
   cmd="sellevel"
   levs=${mlev_evol}
  fi
 else  # Will interpolate to height levels
   levtype="h"
   cmd="ml2hl"
   levs=${hlev_evol}
 fi
else   # Will interpolate to user-specified pessure values
   levtype="p"
   cmd="ml2pl"
   levs=${plev_evol} 
fi
level_list=`echo $levs | sed 's/,/ /g'`

echo
echo "**********************************************************"
echo "***      ICON Tool Kit for APE: evolution plot         ***"
echo "**********************************************************"
echo 

# Temporary variables

script_path=${icon_path}'/scripts/postprocessing/'

resolution=${horizontal_resolution}${vertical_resolution}
fori=${model_data_path}${EXP}_${resolution}
ftmp=${tmp_data_path}${EXP}_${resolution}

if [ $cdo_silence -eq 1 ]; then
   silence='-s'
else
   silence=''
fi

if [ $Model == "ECHAM" ]; then
   compute_remap_weights=0
fi

# The directories for intermediate data and plots will be created, if 
# not already there

if [ ! -d ${plot_file_path} ]; then
   mkdir -p ${plot_file_path} 
fi
if [ ! -d ${tmp_data_path} ]; then
   mkdir -p ${tmp_data_path} 
fi

# Create a directory for soft links. This is used later for 
# reading zonal statistics from multiple files.

datetime=`date +%F-%H%M%S-%N`
lnkdir=${tmp_data_path}"lnk_"$datetime
if [ -d $lnkdir ]; then
   rm -rf $lnkdir
fi
mkdir $lnkdir

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

if [[ ${compute_remap_weights} -eq 1 && ! -a ${weights} ]]; then

 if [ ! -d ${remap_weights_path} ]; then
    mkdir -p ${remap_weights_path} 
 fi

 echo
 echo "=== Computing remapping weights (ICON to Gaussian) ..."

 label=$(printf "%04d" ${evol_istart})
 cdo $silence gendis,t${trunc}grid \
     -seltimestep,1 -selname,PS ${fori}"_"${label}".nc" ${weights}

 check_error $? "Computing remapping weights"

fi

#========================================================================
# Calculate zonal statistics
#========================================================================

    timerange=${evol_istart}"-"${evol_iend}

    # Look up the user-specified variable name in a registry in order to
    # inquire necessary information for data postprocessing and plotting.
    # The script lookup_variable.ksh returns the entry ID (as variable "ie"), 
    # all cause exit if the variable is not found in the registry.

    export varname
    . ./lookup_variable.ksh $varname

    # Pre-plot data processing

    if [ $do_computation -eq 1 ]; then

      echo
      echo "=== Computing zonal statistics for variable $varname ..."

      ifile=$evol_istart
      while [ $ifile -le $evol_iend ] ; do  #loooooooooooooooooooooooooooooooop
        label=$(printf "%04d" $ifile)
        echo File $label

        #------------------------------------------------------
        # Select variables (if this has not been done before)
        #------------------------------------------------------
        if [ ! -f ${ftmp}"_"$label"_"${varname}".nc" ]; then

           case $Model in 
           "ICOHAM" | "ICONAM")

              cdo $silence -r selname,${varname} ${fori}"_"$label".nc" \
                                                 ${ftmp}"_"$label"_"${varname}".nc"
           ;;
           "ECHAM")

             if [ ${afterbn[$ie]} -eq 1 ]; then  # dynamics variables

                cat >$Model_$EXP_$varname.nml <<EOF
                &SELECT
                 TYPE = 20, FORMAT = 2,
                 CODE = ${varcode[$ie]}
                &END
EOF
                after ${fori}"_"$label".nc" \
                      ${ftmp}"_"$label"_"${varcode[$ie]}".nc" \
                      < $Model_$EXP_$varname.nml >after.out

                cdo $silence -r setname,${varname} \
                    ${ftmp}"_"$label"_"${varcode[$ie]}".nc" \
                    ${ftmp}"_"$label"_"${varname}".nc"
  
                rm  ${ftmp}"_"$label"_"${varcode[$ie]}".nc"
                rm  after.out $Model_$EXP_$varname.nml 

             else # physics variables

                cdo $silence -r setname,$varname -selname,${tmpname[$ie]} \
                    ${fori}"_"$label".nc" \
                    ${ftmp}"_"$label"_"$varname".nc"
             fi
       
           ;; 
           *)
             echo "Wrong model name! Abort."
             exit
           esac 
           check_error $? ' - Set/selname'
        fi

        #----------------------------------------------------------------------------
        # Interpolate from model levels to pressure/height levels
        # or simply select the model levels the user specified.
        #----------------------------------------------------------------------------
        case $varname in 
        "PS" | "PHIS")
           echo " - Vertical interpolation/level selection skipped for $varname" 
        ;;
        *) 
          if [[ $levtype == "p" || $levtype == "h" ]]; then
          
            filein=${ftmp}"_"$label"_"$varname"_ml2phl_tmp.nc"

            # Merge PS and PHIS with $varname into one file
            if [ ! -a $filein ]; then
              cdo $silence merge \
                  ${ftmp}"_"$label"_"$varname".nc" \
                  ${ftmp}"_"$label"_"PS".nc"       \
                  ${ftmp}"_"$label"_"PHIS".nc"     \
                  $filein 
              check_error $? " - Merge $varname with PS and PHIS"
            fi

          else
            filein=${ftmp}"_"$label"_"$varname".nc"
          fi

          # Vertical interpolation or level selection

          for lev in $level_list ; do
            fileout=${ftmp}"_"$label"_"$varname"_"$levtype$lev".nc"

            if [ ! -a $fileout ]; then
              cdo $silence $cmd,$lev $filein $fileout
              check_error $? " - ${cmd}, ${lev}"
            fi 
          done
        ;;
        esac #$varname


        #-------------------------------------------------------------------
        # Compute zonal statistics for each time step. Note that
        #  1. For the ICON models, interpolation to Gaus grid is necessary. 
        #  2. We do not merge all time steps into a single file but rather
        #     let the NCL plotting script read data from multiple files. 
        #  3. Soft links to the zonal mean data files are created to facilitate 
        #     data loading in NCL.
        #-------------------------------------------------------------------
        case $varname in 
        "PS" | "PHIS")
           echo " - Computation of zonal statistics skipped for $varname" ;;
        *) 

          case $Model in 
          "ICOHAM" | "ICONAM") intp="-remap,t${trunc}grid,${weights}" ;;
          "ECHAM")             intp="" ;;
          esac

          for lev in $level_list ; do
            filein=${ftmp}"_"$label"_"$varname"_"$levtype$lev".nc"

           #for stat in zonavg zonvar ; do
            for stat in zonavg ; do
              fileout=${ftmp}"_"$label"_"$varname"_"$levtype$lev"_T"$trunc"_"$stat".nc"

              if [ ! -a $fileout ]; then
                cdo $silence $stat $intp $filein $fileout 
                check_error $? " - $stat"
              fi

              flnk=$lnkdir"/"${EXP}_${resolution}
              flnk=$flnk"_"$label"_"$varname"_"$levtype$lev"_T"$trunc"_"$stat".nc"
              if [ ! -a $flnk ]; then
                ln -s $fileout $flnk
              fi
            done
          done

        ;;
        esac
        #----

        ifile=` expr $ifile + 1 `
      done  #loooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooop
    fi # do_computation -eq 1

   #-----------------------------
   # Plotting
   #-----------------------------
   if [ $make_plot -eq 1 ] ; then

     echo
     if [ $varname == "PS" ] || [ $varname == "PHIS" ] ; then
       echo "=== Plotting skipped for ${varname}."
     else

       export resolution
       export timerange
       export evol_istart evol_iend
       export LongName=${longname[$ie]}
       export Scale=${plotscale[$ie]}
       export Min=${plotmin[$ie]}
       export Max=${plotmax[$ie]}
       export Int=${plotint[$ie]}
       export ColorMap=${colormap[$ie]}
       export ColorStart=${colorstart[$ie]}
       export ColorEnd=${colorend[$ie]}
       export DataFNameP1=${lnkdir}"/"${EXP}_${resolution}"_"

       for lev in $level_list ; do
      #for stat in zonavg zonvar ; do
       for stat in zonavg ; do

         export DataFNameP2="_"$varname"_"$levtype$lev"_T"$trunc"_"$stat".nc"
         export PlotFile=${plot_file_path}${EXP}"_"${resolution}"_"$timerange
         export PlotFile=${PlotFile}"_"$varname"_"$levtype$lev"_T"$trunc"_"$stat"-evol"
        
         echo "=== Plotting time evolution ($levtype$lev, $stat) ..."  
         ncl plot_lat-time.ncl >ncl_output_${datetime}.log

         r1=$?
         r2=`grep -i fatal ncl_output_${datetime}.log` ;  echo "$r2"
         check_error $r1 "=== Plotting time evolution ($levtype$lev, $stat) using plot_lat-time.ncl"
         rm ncl_output_${datetime}.log

       done
       done
     fi
   fi #if make_plot -eq 1

#========================================================================
# Clean up
#========================================================================

if [ -d $lnkdir ]; then
   rm -rf $lnkdir
fi
