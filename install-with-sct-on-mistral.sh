#! /bin/bash

ICONROOT=$(pwd)

cd externals

if [[ -e sct ]]
then
    rm -rf sct
fi

prefix=$(pwd)

git clone git@git.mpimet.mpg.de:sct.git

cd sct

module swap intel intel/16.0
module load fca/2.5.2431
module load mxm/3.4.3082
module load bullxmpi_mlx/bullxmpi_mlx-1.2.9.2

HDF5ROOT=/sw/rhel6-x64/hdf5/hdf5-1.8.14-parallel-bullxmpi-intel14

./autogen.sh

CC=icc FC=ifort MPICC=mpicc MPIFC=mpif90 ./configure \
  --prefix=$prefix \
  --enable-openmp \
  --enable-mpi \
  --enable-hdf5 \
  --with-libhdf5-prefix=$HDF5ROOT

make -j 4
make install

cd $ICONROOT

./configure --with-sct=$(pwd)/externals --with-fortran=intel --with-openmmp 

make -j 12

./make_runscripts 

cat << 'EOF'

Add the following lines to your run script:

sct_base_filename="icon-timer"
sct_file_suffix="h5"
sct_file_counter=0

if [[ $restart =~ false ]]
then
    rm -rf ${sct_base_filename}*
else
    sct_last_file=$(ls -1 ${sct_base_filename}* | tail -1)
    t1=${sct_last_file%.${sct_file_suffix}}
    t2=${t1##*-}
    (( sct_file_counter = 10#$t2 + 1 ))
fi

export SCT_OUT="hdf5"
export SCT_FILENAME=$(printf "%s-%04d.%s" $sct_base_filename $sct_file_counter $sct_file_suffix)

if [[ ! -z ${SLURM_JOBID:-} ]]
then
    export SLURM_TIME_FORMAT=standard
    export SCT_JOB_ID="$SLURM_JOB_ID"
    export SCT_JOB_NAME="$SLURM_JOB_NAME"
    export SCT_SUBMIT_DATE=$(sacct -j $SLURM_JOBID -n -o Submit | sed 's/T/ /')
elif [[ ! -z ${PBS_JOBID:-} ]]
then
    export SCT_JOB_ID="$PBS_JOBID"
    export SCT_JOB_NAME="$PBS_JOBNAME"
    export SCT_SUBMIT_DATE=$(qstat -f 4881315.xcepbs00 | awk -F'=' '/qtime/{print $2}' | date +"%Y-%m-%d %H:%M:%S" -f -)
else
    export SCT_JOB_ID="$$"
    export SCT_JOB_NAME="$0"
    export SCT_SUBMIT_DATE=$(date +"%Y-%m-%d %H:%M:%S")
fi

EOF
