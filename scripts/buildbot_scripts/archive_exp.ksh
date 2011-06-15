#!/bin/ksh


DATE=$1
builder=$2
builder_ID=$3

target_server=squall.zmaw.de
target_dir=/scratch/mpi/mh0287/data/archive/${DATE}/buildbot/${builder}

ssh ${target_server} mkdir -p ${target_dir}
scp -r experiments ${target_server}:${target_dir}/${builder_ID}

