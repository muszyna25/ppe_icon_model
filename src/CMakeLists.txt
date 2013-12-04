file (GLOB_RECURSE Fortran_sourcefiles "*.?90")
file (GLOB_RECURSE C_sourcefiles "*.c")

file (GLOB_RECURSE control_model_source "control_model.f90") 
file (GLOB_RECURSE grid_command_source "grid_command.f90") 
file (GLOB_RECURSE partitioner_source "partitioner.f90") 

message (STATUS "control_model: ${control_model_source}")
message (STATUS "grid_command : ${grid_command_source}")
message (STATUS "partitioner  : ${partitioner_source}")

list (REMOVE_ITEM Fortran_sourcefiles 
  ${control_model_source} ${grid_command_source} ${partitioner_source})

add_executable (control_model ${control_model_source} 
  ${Fortran_sourcefiles} 
  ${C_sourcefiles}
  version.c)
add_custom_command (OUTPUT version.c 
  COMMAND ${CMAKE_HOME_DIRECTORY}/config/pvcs.pl --srcdir ${CMAKE_HOME_DIRECTORY})
add_custom_target (version.c)
target_link_libraries (control_model 
  ${MPI_Fortran_LIBRARIES} 
  ${NETCDF_LIBRARIES} 
  ${HDF5_LIBRARIES} 
  mtime support lapack blas)

add_executable (grid_command ${grid_command_source} 
  ${Fortran_sourcefiles} 
  ${C_sourcefiles})
target_link_libraries (grid_command
  ${MPI_Fortran_LIBRARIES} 
  ${NETCDF_LIBRARIES} 
  ${HDF5_LIBRARIES} 
  mtime support lapack blas)

add_executable (partitioner ${partitioner_source} 
  ${Fortran_sourcefiles} 
  ${C_sourcefiles})
target_link_libraries (partitioner
  ${MPI_Fortran_LIBRARIES} 
  ${NETCDF_LIBRARIES} 
  ${HDF5_LIBRARIES} 
  mtime support lapack blas)


