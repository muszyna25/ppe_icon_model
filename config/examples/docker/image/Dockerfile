FROM nvidia/cuda:10.1-devel-ubuntu16.04

RUN useradd -ms /bin/bash cscs

RUN apt-get update && \
    apt-get install \
    ant \
    bison \
    default-jdk \
    environment-modules \
    flex \
    git \
    ksh \
    less \
    libcurl4-openssl-dev \
    libxml2-dev \
    python-pip \
    python2.7 \
    rsync \
    ssh \
    vim \
    wget \
    --yes --no-install-recommends && \
    apt-get clean

RUN pip install --upgrade pip setuptools

WORKDIR /home/cscs

USER cscs

RUN pip install --user easybuild

RUN git clone https://github.com/lucamar/production.git && \
    cd production && \
    git checkout tsa76

ENV PATH=$PATH:/home/cscs/.local/bin \
    EASYBUILD_MODULES_TOOL=EnvironmentModulesC \
    EASYBUILD_MODULE_SYNTAX=Tcl \
    EASYBUILD_ROBOT_PATHS=/home/cscs/production/easybuild/easyconfigs:/home/cscs/.local/easybuild/easyconfigs

RUN eb CMake-3.14.1.eb -r

RUN mkdir -p /home/cscs/.local/easybuild/sources/PGI

COPY pgilinux-2019-1910-x86-64.tar.gz /home/cscs/.local/easybuild/sources/PGI/
COPY PGI-19.10-GCC-7.3.0-2.30.eb /home/cscs/production/easybuild/easyconfigs/p/PGI/
RUN eb PGI-19.10-GCC-7.3.0-2.30.eb -r

COPY Szip-2.1.1-PGI-19.10-GCC-7.3.0-2.30.eb /home/cscs/production/easybuild/easyconfigs/s/Szip/
COPY HDF5-1.10.5-PGI-19.10-GCC-7.3.0-2.30.eb /home/cscs/production/easybuild/easyconfigs/H/HDF5/
RUN eb HDF5-1.10.5-PGI-19.10-GCC-7.3.0-2.30.eb -r

COPY netCDF-4.6.1-PGI-19.10-GCC-7.3.0-2.30.eb /home/cscs/production/easybuild/easyconfigs/n/netCDF/
COPY netCDF-Fortran-4.4.4-PGI-19.10-GCC-7.3.0-2.30.eb /home/cscs/production/easybuild/easyconfigs/n/netCDF-Fortran/
RUN eb netCDF-Fortran-4.4.4-PGI-19.10-GCC-7.3.0-2.30.eb -r

COPY MPICH-3.1.4-PGI-19.10-GCC-7.3.0-2.30.eb /home/cscs/production/easybuild/easyconfigs/m/MPICH/
RUN eb MPICH-3.1.4-PGI-19.10-GCC-7.3.0-2.30.eb -r

COPY JasPer-2.0.14-PGI-19.10-GCC-7.3.0-2.30.eb /home/cscs/production/easybuild/easyconfigs/j/JasPer
RUN eb JasPer-2.0.14-PGI-19.10-GCC-7.3.0-2.30.eb -r

COPY ecCodes-2.13.0-PGI-19.10-GCC-7.3.0-2.30.eb /home/cscs/production/easybuild/easyconfigs/c/ecCodes/
RUN eb ecCodes-2.13.0-PGI-19.10-GCC-7.3.0-2.30.eb

RUN git clone https://github.com/claw-project/claw-compiler.git && \
    cd claw-compiler && \
    git submodule init && \
    git submodule update --remote && \
    mkdir build && \
    cd build && \
    . /usr/share/modules/init/sh && \
    module use /home/cscs/.local/easybuild/modules/all && \
    module load PGI MPICH CMake && \
    FC=mpif90 CC=mpicc CXX=mpicxx cmake -DCMAKE_INSTALL_PREFIX=/home/cscs/.local/claw -DOMNI_MPI_FC="MPI_FC=mpif90" -DOMNI_MPI_CC="MPI_CC=mpicc" .. && \
    make && \
    make install

ENV ICON_DOCKER=1

RUN echo  '. /usr/share/modules/init/bash' >> /home/cscs/.bashrc && \
    echo  'module use /home/cscs/.local/easybuild/modules/all' >> /home/cscs/.bashrc

