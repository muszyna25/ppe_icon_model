#! /bin/ksh

# Specify the experiment to be processed

Model=ICOHAM
EXP=hat_ape_qobs_moist_RH00-clip  # String will show up in the plot
exp_config="$Model APE w/o O3"

cell_type=3
horizontal_resolution=R2B04
vertical_resolution=L31

model_data_path=/scratch/work/mh0287/hwan/exps/icon-dev/output_201007_port_phy/$EXP/
remap_weights_path=$model_data_path/remap_weights/
plot_file_path=$model_data_path/plots/
tmp_data_path=$model_data_path/tmp/

compute_remap_weights=1
grid_optimization="spr0.90"
trunc=63
diag_climate=1
clim_istart=171
clim_iend=172
do_computation=1
make_plot=1
plot_file_format=pdf
TopHeight=35
plevs=""

# Specify the reference experiment for making the difference plot

ref_config="same run, days 151-220"  # String will show up in the difference plot
ref_exp_name=$EXP
ref_resolution=${horizontal_resolution}${vertical_resolution}
ref_timerange="151-220"
ref_data_path=$tmp_data_path

# define individual set-env file
set_env="Post.${EXP}.set_env"

#
cd ./APE_postpro

# Remove old ${set_env} file
 
if [ -f ${set_env} ]; then
  rm -f ${set_env}
fi

# Write new ${set_env} file
 
echo Model=$Model                                 >  ${set_env}
echo EXP=$EXP                                     >> ${set_env}
echo export exp_config=\"$exp_config\"            >> ${set_env}

echo cell_type=$cell_type                         >> ${set_env}
echo horizontal_resolution=$horizontal_resolution >> ${set_env}
echo vertical_resolution=$vertical_resolution     >> ${set_env}

echo model_data_path=$model_data_path             >> ${set_env}
echo remap_weights_path=$remap_weights_path       >> ${set_env}
echo plot_file_path=$plot_file_path               >> ${set_env}
echo tmp_data_path=$tmp_data_path                 >> ${set_env}

echo compute_remap_weights=$compute_remap_weights >> ${set_env}
echo grid_optimization=$grid_optimization         >> ${set_env}
echo trunc=$trunc                                 >> ${set_env}
echo diag_climate=$diag_climate                   >> ${set_env}
echo clim_istart=$clim_istart                     >> ${set_env}
echo clim_iend=$clim_iend                         >> ${set_env}
echo do_computation=$do_computation               >> ${set_env}
echo make_plot=$make_plot                         >> ${set_env}
echo export plot_file_format=$plot_file_format    >> ${set_env}
echo export TopHeight=$TopHeight                  >> ${set_env}
echo plevs=\"${plevs}\"                           >> ${set_env}

echo export ref_config=\"$ref_config\"            >> ${set_env}
echo ref_exp_name=$ref_exp_name                   >> ${set_env}
echo ref_resolution=$ref_resolution               >> ${set_env}
echo ref_data_path=$ref_data_path                 >> ${set_env}
echo ref_timerange=$ref_timerange                 >> ${set_env}

# Remove old "ncl_output.log"  file

if [ -f ncl_output.log ]; then
  rm -f ncl_output.log
fi

# Run script
#
# Note: Make sure that
#   1. all variable listed below exist in model output and 
#      have been registered in "lookup_variable.ksh";
#   2. PS and PHIS are processed before any other variable
#      (because they are needed for the vertical 
#      interpolation from model levels to pressure levels).
 
for var in PS PHIS T U V OMEGA Qv Qw Qi ACLC ; do
  ./zonal_clim.ksh $var

  # Switch off repeated computation
  echo compute_remap_weights=0 >> ${set_env}
done

exit
