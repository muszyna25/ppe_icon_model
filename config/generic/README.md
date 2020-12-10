<!--
This file is written using Markdown language, which might make it difficult to
read it in a plain text editor. Please, visit ICON project page on DKRZ GitLab
(https://gitlab.dkrz.de/icon/icon/-/tree/master/config/generic) to see this file
rendered or use a Markdown viewer of your choice
(https://www.google.com/search?q=markdown+viewer).
-->

# Introduction

This directory contains generic
[configure wrappers](../../README.md##configuration-wrappers) for ICON. You can
use them either inside the [docker containers](#docker-containers) or
[natively](#native-building) on your system.

# Docker containers

The easiest way to build and run ICON on your personal machine is to use Docker
images from the
[Gitlab container registry](https://gitlab.dkrz.de/icon/icon-cimd/container_registry)
associated with the [icon-cimd project](https://gitlab.dkrz.de/icon/icon-cimd).

First, you will need to [install Docker](https://docs.docker.com/get-docker/) on
your machine. Once that is done, you need to log in to the registry with the
same credentials you are using to access ICON source code repository at
[gitlab.dkrz.de](https://gitlab.dkrz.de/):
```bash
docker login registry.gitlab.dkrz.de
```

The recommended workflow is to build and run ICON inside the container and
manage and edit the source code in a separate terminal, using the tools
available on your machine. This scenario implies that the directory with ICON
source code is located on the host machine and mounted to the container. You
will also need to mount a so-called `pool` directory containing ICON input
files, e.g. grid files. The contents and the layout of the directory depends on
the experiment you want to run and are not covered in this document.

Run the container in the interactive mode as follows:
```bash
docker run -it -v /path/to/icon-src:/home/icon/icon -v /path/to/pool:/home/icon/pool registry.gitlab.dkrz.de/icon/icon-cimd/icon-dev
```
where `/path/to/icon-src` and `/path/to/pool` are paths to ICON source and
`pool` directories on your machine, and `/home/icon/icon` and `/home/icon/pool`
are respective mount points of the directories inside the container.

As a result of the previous command you will get an interactive command prompt
of the container. You can now configure, build and run ICON using the following
commands as a reference:
```console
icon@dev-gcc-6.5.0$ cd ./icon
icon@dev-gcc-6.5.0$ ./config/generic/gcc
icon@dev-gcc-6.5.0$ make -j4
icon@dev-gcc-6.5.0$ ./make_runscripts -s atm_ampi_test
icon@dev-gcc-6.5.0$ cd ./run
icon@dev-gcc-6.5.0$ ./exp.atm_amip_test.run
```
> **_NOTE:_** To be able to run ICON inside the container, you might need to
increase the amount of RAM available to Docker (Preferences->Resources->Memory).

# Native building

The generic wrappers in this directory are written with the assumption that the
required [software libraries](#software-libraries) are installed to the same
prefix. The prefix defaults to `/opt/local` on macOS and to `/usr` on other
platforms. The default values can be overridden by setting the environment
variable `ICON_SW_PREFIX`:
```bash
export ICON_SW_PREFIX='/path/to/icon/prerequisites'
```

## Prerequisites

This section provides a list of software required for building and running
ICON. The users can build and install (to the same prefix) the listed packages
manually or use a
[package managers](https://en.wikipedia.org/wiki/Package_manager) available for
their platform. The basic instructions on how to do it on several popular
platforms are provided in section [Tested platforms](#tested-platforms).

### Building tools

- [GNU Make](https://www.gnu.org/software/make) v3.81+
- [Python](https://www.python.org) v2.6+ or v3.5+
- [Perl](https://www.perl.org) v5.10+
- Interoperable C and Fortran compilers (see also known
[compiler issues](https://gitlab.dkrz.de/icon/icon/-/boards/189))

### Software libraries

- [MPICH](https://www.mpich.org), [OpenMPI](https://www.open-mpi.org) or any
other [MPI](https://www.mpi-forum.org) implementation that provides compiler
wrappers `mpicc` and `mpif90` for C and Fortran, repsectively, as well as the
job launcher `mpiexec`.
- [HDF5](https://support.hdfgroup.org/HDF5) with high-level interface (for
<a href="#netcdf-c">NetCDF-С</a>), thread-safety (for <a href="#cdo">CDO</a>),
and szlib filtering support (only C interface required, not a direct dependency
of ICON)
- <a name="netcdf-c"/> [NetCDF-C](https://www.unidata.ucar.edu/software/netcdf/docs)
with NetCDF-4 support
- [NetCDF-Fortran](https://www.unidata.ucar.edu/software/netcdf/docs-fortran)
- [ecCodes](https://confluence.ecmwf.int/display/ECC) with JPEG2000 and AEC
support (only C interface required)
- [BLAS](http://www.netlib.org/blas)
- [LAPACK](http://www.netlib.org/lapack)
- [Libxml2](http://www.xmlsoft.org)

See section [ICON dependencies](../../README.md#icon-dependencies) for more
details.

### Optional tools

- [KornShell](http://www.kornshell.com) for the
[generated runscripts](../../README.md#running)
- <a name="cdo"/> [CDO](https://code.mpimet.mpg.de/projects/cdo) for pre- and
post-processing, also used by some of the
[generated runscripts](../../README.md#running)
- [rsync](https://rsync.samba.org/) for the
[generated runscipts](../../README.md#running) in the case of
[out-of-source building](../../README.md#out-of-source-configuration-building)

## Tested platforms and tools

This section provides basic instructions on how to install most commonly
required subset of [ICON dependencies](../../README.md#icon-dependencies) on
different operating systems using relevant
[package managers](https://en.wikipedia.org/wiki/Package_manager).

> **_NOTE:_** The current recommended version of [GCC](https://gcc.gnu.org/) is
**6.x**. This is why most of the examples below are based on this version.

### macOS with [MacPorts](https://www.macports.org)

**Tested on `macOS Catalina 10.15.5`.**

Most of the required software packages are either already available on the
system or installed together with [Xcode](https://developer.apple.com/xcode) and
[Command Line Tools for Xcode](https://developer.apple.com/download/more/),
which are [prerequsites for MacPorts](https://www.macports.org/install.php). The
rest of the required software can be installed by running the following
commands:

```bash
# Install building tools and ICON dependencies:
sudo port -N install       \
  gcc6                     \
  mpich-gcc6               \
  hdf5 +hl+threadsafe+szip \
  netcdf                   \
  netcdf-fortran +gcc6     \
  eccodes                  \
  libxml2

# Select the compiler and MPI compiler wrappers:
sudo port select --set gcc mp-gcc6
sudo port select --set mpi mpich-gcc6-fortran
hash -r

# Install optional tools:
sudo port -N install cdo +netcdf
```

### Ubuntu with [Apt](https://wiki.debian.org/Apt)

**Tested on `Ubuntu Bionic Beaver 18.04.4 LTS`.**

```bash
# Install building tools and ICON dependencies:
sudo apt install -y \
  build-essential   \
  python3           \
  gcc-6             \
  gfortran-6        \
  libmpich-dev      \
  libhdf5-dev       \
  libnetcdf-dev     \
  libnetcdff-dev    \
  libeccodes-dev    \
  libblas-dev       \
  liblapack-dev     \
  libxml2-dev

# Select the compiler:
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-6 50 \
  --slave /usr/bin/gfortran gfortran /usr/bin/gfortran-6

# The command above can be reverted as follows:
# sudo update-alternatives --remove gcc /usr/bin/gcc-6

# Select MPI libraries and compiler wrappers
sudo update-alternatives --set mpi /usr/include/mpich
sudo update-alternatives --set mpirun /usr/bin/mpirun.mpich

# If the two non-interactive commands above do not work,
# try the interactive analogues:
# sudo update-alternatives --config mpirun
# sudo update-alternatives --config mpi

# Install optional tools:
sudo apt install -y ksh cdo
```

### Arch Linux with [Pacman](https://wiki.archlinux.org/index.php/pacman)

**Tested on `Arch Linux 2020.07.01`.**

```bash
# Install building tools and ICON dependencies:
sudo pacman -S --noconfirm \
  base-devel               \
  git                      \
  python                   \
  gcc-fortran              \
  openmpi                  \
  hdf5                     \
  netcdf                   \
  netcdf-fortran           \
  blas                     \
  lapack                   \
  libxml2

# Install ecCodes from the Arch User Repository:
git clone https://aur.archlinux.org/openjpeg.git && cd openjpeg && makepkg -csi --noconfirm && cd ..
git clone https://aur.archlinux.org/eccodes.git && cd eccodes && makepkg -csi --noconfirm && cd ..

# Install optional tools:
sudo pacman -S --noconfirm rsync ksh
```

> **_NOTE:_** [ecCodes package](https://aur.archlinux.org/packages/eccodes) is
currently broken and its `PKGBUILD` file might need to be extended with an
additional CMake argument `-DCMAKE_Fortran_FLAGS=-fallow-argument-mismatch`.

### macOS/Linux with [Spack](https://spack.io)

**Tested on `Ubuntu Bionic Beaver 18.04.4 LTS`.**

> **_NOTE:_** Spack has its own list of
[prerequisites](https://spack.readthedocs.io/en/latest/getting_started.html#prerequisites).
Most of them are available out-of-the box on most Unix systems. On Ubuntu
systems, however, the users might need to install them separetly:
>```bash
>sudo apt install -y build-essential python3
>```

> **_NOTE:_** Although it is possible to install Fortran and C compilers with
Spack, the current recommendation is to do it using system's default package
manager. For example, on Ubuntu system, the users can run:
>```bash
>sudo apt install -y gcc-6 g++-6 gfortran-6
>```

> **_NOTE:_** Currently, Spack does not provide `KornShell`. The users are
recommended to install it using system's default package manager. For example,
on Ubuntu system, the users can run:
>```bash
>sudo apt install -y ksh
>```

```bash
# Install Spack:
git clone https://github.com/spack/spack.git
. ./spack/share/spack/setup-env.sh

# Find the compilers;
spack compiler find

# Ignore the unneeded compilers:
spack compiler remove --all gcc@7:

# Install ICON dependencies:
spack install openmpi &&
  spack install netcdf-fortran ^hdf5+hl+szip+threadsafe &&
  spack install eccodes &&
  spack install netlib-lapack &&
  spack install libxml2

# Symlink the dependencies to a single prefix (e.g. to $HOME/icon-sw):
export ICON_SW_PREFIX="$HOME/icon-sw"
spack view symlink -i "$ICON_SW_PREFIX" \
  openmpi                               \
  netcdf-fortran                        \
  eccodes                               \
  netlib-lapack                         \
  libxml2

# Install optional tools:
spack install cdo ^hdf5+hl+szip+threadsafe
spack load cdo
```