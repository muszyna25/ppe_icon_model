#!/bin/bash

# Slurm job options (job-name, compute nodes, job time)

#SBATCH --job-name=test_ICON_cycle2
#SBATCH --account=n02-RECAP
#SBATCH --mem=0

# standard Q
##SBATCH --time=24:00:0
##SBATCH --nodes=8
##SBATCH --tasks-per-node=128
##SBATCH --cpus-per-task=1
##SBATCH --partition=standard
##SBATCH --qos=standard

# short Q
#SBATCH --time=0:10:00
#SBATCH --nodes=4
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --reservation=shortqos

export icon_dir=/work/n02/n02/rherbe/ICON/nextgems_cycle2

. $icon_dir/modules_cray
module load libfabric
module load atp

export ATP_ENABLED=1
export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH

cd $icon_dir/run
./TEST_cr2b4_n2_OXF_5day.run
