!>
!! @brief Subroutine interface_echam_par calls a generic parameterization "echam_par".
!!  The parameterization is called in a loop over all rows of a block containing all
!!  columns of a given patch of the grid. The loop is embedded in an OMP region.
!!
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version 2017-12
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_interface_echam_par

  USE mo_kind           ,ONLY: wp
  USE mtime             ,ONLY: datetime
  USE mo_model_domain   ,ONLY: t_patch
  USE mo_loopindices    ,ONLY: get_indices_c
  USE mo_impl_constants ,ONLY: min_rlcell_int, grf_bdywidth_c

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_par

CONTAINS

  SUBROUTINE interface_echam_par(echam_par,            &
       &                         is_in_sd_ed_interval, &
       &                         is_active,            &
       &                         patch,                &
       &                         datetime_old,         &
       &                         pdtime                )

    INTERFACE
       !
       SUBROUTINE echam_par(is_in_sd_ed_interval, &
            &               is_active,            &
            &               jg, jb, jcs, jce,     &
            &               datetime_old,         &
            &               pdtime                )
         !
         IMPORT :: wp, datetime
         !
         LOGICAL            ,INTENT(in) :: is_in_sd_ed_interval
         LOGICAL            ,INTENT(in) :: is_active
         INTEGER            ,INTENT(in) :: jg, jb, jcs, jce
         TYPE(datetime)     ,POINTER    :: datetime_old
         REAL(wp)           ,INTENT(in) :: pdtime
         !
       END SUBROUTINE echam_par
       !
    END INTERFACE

    LOGICAL                 ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in) :: is_active
    TYPE(t_patch)   ,TARGET ,INTENT(in) :: patch
    TYPE(datetime)          ,POINTER    :: datetime_old
    REAL(wp)                ,INTENT(in) :: pdtime

    INTEGER  :: jg          !< grid index
    INTEGER  :: ncd         !< number of child domains of grid jg
    INTEGER  :: rls, rle    !< row limits, start and end

    INTEGER  :: jb          !< index of block loop
    INTEGER  :: jbs, jbe    !< start and end indices of block   loop
    INTEGER  :: jcs, jce    !< start and end indices of columns loop

    jg  = patch%id

    rls = grf_bdywidth_c+1
    rle = min_rlcell_int

    ncd = MAX(1,patch%n_childdom)
    jbs = patch%cells%start_blk(rls,  1)
    jbe = patch%cells%  end_blk(rle,ncd)

!$OMP PARALLEL DO PRIVATE(jcs,jce)
    DO jb = jbs,jbe
       !
       CALL get_indices_c(patch, jb, jbs, jbe, jcs, jce, rls, rle)
       !
       CALL echam_par(is_in_sd_ed_interval, &
            &         is_active,            &
            &         jg, jb, jcs, jce,     &
            &         datetime_old,         &
            &         pdtime                )
       !
    END DO
!$OMP END PARALLEL DO 

  END SUBROUTINE interface_echam_par

END MODULE mo_interface_echam_par
