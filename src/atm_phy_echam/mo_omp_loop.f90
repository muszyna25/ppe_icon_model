!>
!! @brief Subroutine omp_loop_cell/edge_par/dia calls a generic subroutine in an
!!  OMP parallelized loop over all blocks of a given patch of the grid.
!!
!!  omp_loop_cell/edge_par passes indices and time control arguments.
!!  omp_loop_cell/edge_dia passes only indices.
!!
!!$!!  omp_loop_cell_par/dia works on cells of a patch
!!$!!  omp_loop_edge_par/dia works on edges of a patch
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

MODULE mo_omp_loop

  USE mo_kind               ,ONLY: wp
  USE mtime                 ,ONLY: datetime
  USE mo_parallel_config    ,ONLY: nproma
  USE mo_model_domain       ,ONLY: t_patch
  USE mo_loopindices        ,ONLY: get_indices_c
  USE mo_impl_constants     ,ONLY: min_rlcell_int
  USE mo_impl_constants_grf ,ONLY: grf_bdywidth_c

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: omp_loop_cell, omp_loop_cell_tc

CONTAINS


  SUBROUTINE omp_loop_cell_tc(patch                ,&
       &                      routine              ,&
       &                      is_in_sd_ed_interval ,&
       &                      is_active            ,&
       &                      datetime_old         ,&
       &                      pdtime               )

    ! Arguments
    !
    TYPE(t_patch)   ,TARGET ,INTENT(in) :: patch
    !
    INTERFACE
       !
       SUBROUTINE routine(jg,jb,jcs,jce        ,&
            &             nproma,nlev          ,& 
            &             is_in_sd_ed_interval ,&
            &             is_active            ,&
            &             datetime_old         ,&
            &             pdtime               )
         !
         IMPORT :: wp, datetime
         !
         INTEGER        ,INTENT(in) :: jg,jb,jcs,jce
         INTEGER        ,INTENT(in) :: nproma,nlev
         LOGICAL        ,INTENT(in) :: is_in_sd_ed_interval
         LOGICAL        ,INTENT(in) :: is_active
         TYPE(datetime) ,POINTER    :: datetime_old
         REAL(wp)       ,INTENT(in) :: pdtime
         !
       END SUBROUTINE routine
       !
    END INTERFACE
    !
    LOGICAL                 ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL                 ,INTENT(in) :: is_active
    TYPE(datetime)          ,POINTER    :: datetime_old
    REAL(wp)                ,INTENT(in) :: pdtime

    ! Local variables
    !
    INTEGER  :: jg          !< grid index
    INTEGER  :: ncd         !< number of child domains of grid jg
    INTEGER  :: rls, rle    !< row limits, start and end
    INTEGER  :: jb          !< index of block loop
    INTEGER  :: jbs, jbe    !< start and end indices of block   loop
    INTEGER  :: jcs, jce    !< start and end indices of columns loop
    INTEGER  :: nlev        !< number of levels

    jg  = patch%id

    rls = grf_bdywidth_c+1
    rle = min_rlcell_int

    ncd = MAX(1,patch%n_childdom)
    jbs = patch%cells%start_blk(rls,  1)
    jbe = patch%cells%  end_blk(rle,ncd)

    nlev= patch%nlev

!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = jbs,jbe
       !
       CALL get_indices_c(patch,jb,jbs,jbe,jcs,jce,rls,rle)
       IF (jcs>jce) CYCLE
       !
       CALL routine(jg,jb,jcs,jce        ,&
            &       nproma,nlev          ,&
            &       is_in_sd_ed_interval ,&
            &       is_active            ,&
            &       datetime_old         ,&
            &       pdtime               )
       !
    END DO
!$OMP END PARALLEL DO 

  END SUBROUTINE omp_loop_cell_tc


  SUBROUTINE omp_loop_cell(patch  ,&
       &                   routine)

    ! Arguments
    !
    TYPE(t_patch)   ,TARGET ,INTENT(in) :: patch
    !
    INTERFACE
       !
       SUBROUTINE routine(jg,jb,jcs,jce ,&
            &             nproma,nlev   )
         !
         INTEGER        ,INTENT(in) :: jg,jb,jcs,jce
         INTEGER        ,INTENT(in) :: nproma,nlev
         !
       END SUBROUTINE routine
       !
    END INTERFACE

    ! Local variables
    !
    INTEGER  :: jg          !< grid index
    INTEGER  :: ncd         !< number of child domains of grid jg
    INTEGER  :: rls, rle    !< row limits, start and end
    INTEGER  :: jb          !< index of block loop
    INTEGER  :: jbs, jbe    !< start and end indices of block   loop
    INTEGER  :: jcs, jce    !< start and end indices of columns loop
    INTEGER  :: nlev        !< number of levels

    jg  = patch%id

    rls = grf_bdywidth_c+1
    rle = min_rlcell_int

    ncd = MAX(1,patch%n_childdom)
    jbs = patch%cells%start_blk(rls,  1)
    jbe = patch%cells%  end_blk(rle,ncd)

    nlev= patch%nlev

!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = jbs,jbe
       !
       CALL get_indices_c(patch,jb,jbs,jbe,jcs,jce,rls,rle)
       IF (jcs>jce) CYCLE
       !
       CALL routine(jg,jb,jcs,jce ,&
            &       nproma,nlev   )
       !
    END DO
!$OMP END PARALLEL DO 

  END SUBROUTINE omp_loop_cell

END MODULE mo_omp_loop
