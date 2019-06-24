# New configuration script

This is a branch of ICON v2.6.0-rc with a new building system.

## Dependencies
Fig. 1 shows a partial dependency tree of the model. A dependency can be
mandatory or optional and some of the dependencies are provided together
with the source code of ICON (bundled) as git submodules (see
[.gitmodules](/.gitmodules)), which might make it easier to build the model.

<a name="icon-depgraph"/>![ICON dependency tree](/doc/config/icon-config-doc-depgraph.png)

*Fig. 1. ICON dependency tree*</a>

Run the following command to see the list of possible configuration options and
their default values:
```console
    $ ./configure --help
```

The list of the required libraries (packages) depends on the selected options
(see [Table 1](#icon-deptable)). To make the configuration process more
transparent, the configure script is implemented without arguments
`--with-package` that would accept paths to the installation directories of the
packages, which would be used to extend the corresponding compiler flags.
Instead, all the libraries need to be provided as `FCFLAGS`, `LDFLAGS` and
`LIBS`. Some of the packages required by ICON have common dependencies, which
might make it difficult to specify the libraries in the right topological
order. For that reason, we provide the following table. We recommend to specify
the `LIBS` argument according to the order indicated in the table, and the rest
of the compiler flag arguments (i.e. `FCFLAGS`, `CPPFLAGS`, `LDFLAGS`) in the
reversed order. You can find configuration examples in
[/config/examples](/config/examples).

Upon successful configuration, the building can be started with the standard
`make` command:
```console
    $ make -j8
```

<a name="icon-deptable">*Table 1. ICON dependencies*</a>

| Linking order | Package | Repo | Dependency condition | Required flags |
| :---: | :---: | :---: | :---: | :---: |
|  1 | [SELF](https://code.mpimet.mpg.de/projects/self-standard-extendible-library-for-fortran) | git@git.mpimet.mpg.de:libself.git (branch: `libself-icon_v0.2-config`) | `--with-external-self` | `FCFLAGS='-I/path/to/libself/include' LDFLAGS='-L/path/to/libself/lib' LIBS='-lself'` |
|  2 | ICON-TIXI (a modified version of [TIXI](https://github.com/DLR-SC/tixi)) | git@gitlab.dkrz.de:icon-new-config/icon-tixi.git (branch: `new-config-release`) | `--enable-art --with-external-tixi` | `FCFLAGS='-I/path/to/tixi/include' LDFLAGS='-L/path/to/tixi/lib' LIBS='-licon_tixi'` |
|  3 | [YAC](https://doc.redmine.dkrz.de/YAC/html/) | git@gitlab.dkrz.de:icon-new-config/YAC-dev.git (branch: `new-config-release`) | `--enable-coupling --with-external-yac` | `FCFLAGS='-I/path/to/yac/include' LDFLAGS='-L/path/to/yac/lib -L/path/to/libxml2/lib -L/path/to/netcdf/lib' LIBS='-lyac -lxml2 -lnetcdf'` |
|  4 | [XML2](http://www.xmlsoft.org/) | | `--enable-coupling` or `--enable-art --without-external-tixi` | `CPPFLAGS='-I/path/to/libxml2/include/libxml2' LDFLAGS='-L/path/to/libxml2/lib' LIBS='-lxml2'` |
|  5 | [LAPACK](http://www.netlib.org/lapack/) (or analogue) | | mandatory | `LDFLAGS='-L/path/to/lapack/lib' LIBS='-llapack'` (depends on the implementation) |
|  6 | [BLAS](http://www.netlib.org/blas/) (or analogue) | | mandatory | `LDFLAGS='-L/path/to/blas/lib' LIBS='-lblas'` (depends on the implementation) |
|  7 | [MTIME](https://code.mpimet.mpg.de/projects/mtime) | git@git.mpimet.mpg.de:libmtime.git (branch: `m300488/icon-new-config`) | `--with-external-mtime` or `--without-external-mtime --enable-coupling --with-external-yac` | `FCFLAGS='-I/path/to/mtime/include' CPPFLAGS='-I/path/to/mtime/include' LDFLAGS='-L/path/to/mtime/lib' LIBS='-lmtime'` |
|  8 | [SERIALBOX2](https://github.com/eth-cscs/serialbox2) | | `--enable-serialization` | `FCFLAGS='-I/path/to/serialbox2/include' LDFLAGS='-L/path/to/serialbox2/lib' LIBS='-lSerialboxFortranShared'` |
|  9 | [CDI](https://code.mpimet.mpg.de/projects/cdi/) | git@git.mpimet.mpg.de:libcdi.git (branch: `cdi-1.8.x`) | `--with-external-cdi` | `FCFLAGS='-I/path/to/libcdi/include' LDFLAGS='-L/path/to/libcdi/lib' LIBS='-lcdi_f2003 -lcdi'` |
| 10 | [ECCODES](https://confluence.ecmwf.int/display/ECC) or [GRIB-API](https://confluence.ecmwf.int/display/GRIB/Home) | | `--enable-grib2 --without-external-cdi` or `--enable-yaxt --without-external-yaxt --with-external-cdi` (as a workaround for YAXT, which helps only if `--enable-rpaths`, otherwise, requires setting RPATH flags manually) | `CPPFLAGS='-I/path/to/eccodes/include' LDFLAGS='-L/path/to/eccodes/lib' LIBS='-leccodes'` (or `LIBS='-lgrib_api'`) |
| 11 | [YAXT](https://www.dkrz.de/redmine/projects/yaxt/wiki) | git@gitlab.dkrz.de:icon-new-config/yaxt.git (branch: `release-0.7.0-patched`) | `--enable-yaxt --with-external-yaxt` | `FCFLAGS='-I/path/to/yaxt/include' LDFLAGS='-L/path/to/yaxt/lib' LIBS='-lyaxt'` |
| 12 | [SCT](https://code.mpimet.mpg.de/projects/performance-monitoring/wiki/Access_of_stored_performance_data) | git@gitlab.dkrz.de:icon-new-config/sct.git (branch: `master-release`) | `--enable-sct --with-external-sct` | `FCFLAGS='-I/path/to/sct/include' LDFLAGS='-L/path/to/sct/lib' LIBS='-lsct'` |
| 13 | RTTOV (a modified? version of [RTTOV](https://www.nwpsaf.eu/site/software/rttov/)) | | `--enable-rttov` | `FCFLAGS='-I/path/to/rttov/include' LDFLAGS='-L/path/to/rttov/lib' LIBS='-lradiance -lrttov10.2'` |
| 14 | [ECRAD](https://confluence.ecmwf.int/display/ECRAD/ECMWF+Radiation+Scheme+Home) | git@gitlab.dkrz.de:icon-new-config/ecrad.git (branch: `master`) | `--enable-ecrad --with-external-ecrad` | `FCFLAGS='-I/path/to/ecrad/include' LDFLAGS='-L/path/to/ecrad/lib' LIBS='-lradiation -lifsrrtm -lutilities -lifsaux'` |
| 15 | [RTE+RRTMGP](https://github.com/RobertPincus/rte-rrtmgp) | git@github.com:skosukhin/rte-rrtmgp.git (branch: `icon`) | `--enable-rte-rrtmgp --with-external-rte-rrtmgp` | `FCFLAGS='-I/path/to/rte-rrtmgp/include' LDFLAGS='-L/path/to/rte-rrtmgp/lib' LIBS='-lrrtmgp -lrte'` |
| 16 | [NetCDF-Fortran](https://www.unidata.ucar.edu/software/netcdf/docs-fortran/) | | mandatory | `FCFLAGS='-I/path/to/netcdf-fortran/include' LDFLAGS='-L/path/to/netcdf-fortran/lib' LIBS='-lnetcdff'` |
| 17 | [NetCDF-C](https://www.unidata.ucar.edu/software/netcdf/docs/) | | `--without-external-cdi` or `--enable-coupling` | `CPPFLAGS='-I/path/to/netcdf/include' LDFLAGS='-L/path/to/netcdf/lib' LIBS='-lnetcdf'` |
| 18 | [HDF5](https://support.hdfgroup.org/HDF5/) | | `--enable-sct --without-external-sct` | `CPPFLAGS='-I/path/to/hdf5/include' LDFLAGS='-L/path/to/hdf5/lib' LIBS='-lhdf5'` |
| 19 | [AEC](https://gitlab.dkrz.de/k202009/libaec) or [SZIP](https://support.hdfgroup.org/doc_resource/SZIP/) | | static linking | `LDFLAGS='-L/path/to/aec/lib' LIBS='-laec'` (or `LIBS='-lsz'`) |
| 20 | [MPI](https://www.mpi-forum.org/) (Fortran interface) | | `--enable-mpi` or `--enable-yaxt --without-external-yaxt` or `--enable-coupling --without-external-yac` (both mpif.h and mpi.mod interfaces are required) | `FC='/path/to/mpi/bin/mpif90'` or `FCFLAGS='-I/path/to/mpi/include' LDFLAGS='-L/path/to/mpi/lib' LIBS='-lmpifort -lmpi'` (depends on the implementation) |
| 21 | [MPI](https://www.mpi-forum.org/) (C interface) | | `--enable-mpi --enable-coupling --without-external-yac` or `--enable-yaxt --without-external-yaxt` or `--enable-mpi --enable-sct --without-external-sct` | `CC=/path/to/mpi/bin/mpicc` or `CPPFLAGS=-I/path/to/mpi/include LDFLAGS='-L/path/to/mpi/lib' LIBS='-lmpi'` (depends on the implementation) |


