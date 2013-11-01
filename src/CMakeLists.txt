file(GLOB_RECURSE Fortran_sourcefiles "*.?90")
file(GLOB_RECURSE C_sourcefiles "*.c")


add_executable (control_model ${Fortran_sourcefiles} ${C_sourcefiles})
add_executable (grid_generator ${Fortran_sourcefiles} ${C_sourcefiles})
add_executable (partitioner ${Fortran_sourcefiles} ${C_sourcefiles})

target_link_libraries (control_model ${MPI_Fortran_LIBRARIES} ${NETCDF_LIBRARIES} mtime support lapack blas)

