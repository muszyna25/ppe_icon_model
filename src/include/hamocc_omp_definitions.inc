!-----------------------------------------------
! define default omp parameters
#define HAMOCC_OMP_PARALLEL $omp parallel
#define HAMOCC_OMP_END_PARALLEL $omp end parallel
#define HAMOCC_OMP_DO $omp do
#define HAMOCC_OMP_END_DO $omp end do
#define HAMOCC_OMP_END_DO_NOWAIT $omp end do nowait
#define HAMOCC_OMP $omp
#define HAMOCC_OMP_PARALLEL_DO $omp parallel do
#define HAMOCC_OMP_END_PARALLEL_DO $omp end parallel do
#define HAMOCC_OMP_BARRIER $omp barrier
#define HAMOCC_OMP_SINGLE $omp single
#define HAMOCC_OMP_END_SINGLE $omp end single
#define HAMOCC_OMP_MASTER $omp master
#define HAMOCC_OMP_END_MASTER $omp end master
#define HAMOCC_OMP_WORKSHARE $omp workshare
#define HAMOCC_OMP_END_WORKSHARE $omp end workshare
!-----------------------------------------------

#ifdef _CRAYFTN
#define HAMOCC_OMP_DEFAULT_SCHEDULE SCHEDULE(guided)
#define HAMOCC_OMP_RUNTIME_SCHEDULE SCHEDULE(runtime)
#define HAMOCC_OMP_GUIDED_SCHEDULE  SCHEDULE(guided)
#else 
#define HAMOCC_OMP_DEFAULT_SCHEDULE SCHEDULE(dynamic,1)
#define HAMOCC_OMP_RUNTIME_SCHEDULE SCHEDULE(runtime)
#define HAMOCC_OMP_GUIDED_SCHEDULE  SCHEDULE(guided)
#endif


!-----------------------------------------------

