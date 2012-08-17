!-----------------------------------------------
! define default omp parameters
#define ICON_OMP_PARALLEL $omp parallel
#define ICON_OMP_END_PARALLEL $omp end parallel
#define ICON_OMP_END_DO $omp end do
#define ICON_OMP_END_DO_NOWAIT $omp end do nowait

#ifdef __SX__
#define ICON_OMP_DO $omp do
#define ICON_OMP_DEFAULT_SCHEDULE
#define ICON_OMP_RUNTIME_SCHEDULE SCHEDULE(runtime)
#define ICON_OMP_GUIDED_SCHEDULE  SCHEDULE(guided)
#else
#define ICON_OMP_DO $omp do schedule(dynamic,1)
#define ICON_OMP_DEFAULT_SCHEDULE SCHEDULE(dynamic,1)
#define ICON_OMP_RUNTIME_SCHEDULE SCHEDULE(dynamic,1)
#define ICON_OMP_GUIDED_SCHEDULE  SCHEDULE(dynamic,1)
#endif

!-----------------------------------------------

