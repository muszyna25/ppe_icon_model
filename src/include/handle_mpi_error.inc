#define HANDLE_MPI_ERROR_USE use mo_mpi, only: handle_mpi_error, mpi_success
#define HANDLE_MPI_ERROR(ierror, mpi_routine) IF (ierror /= mpi_success) \
call handle_mpi_error(ierror, routine, __LINE__, mpi_routine)
