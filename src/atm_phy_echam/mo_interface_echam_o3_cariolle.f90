!>
!! @brief Subroutine echam_phy_main calls all the parameterization schemes
!!
!! @author Hui Wan, MPI-M
!! @author Marco Giorgetta, MPI-M
!!
!! @par Revision History
!!  Original version from ECHAM6 (revision 2028)
!!  Modified for ICOHAM by Hui Wan and Marco Giorgetta (2010)
!!  Modified for ICONAM by Marco Giorgetta (2014)
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

MODULE mo_interface_echam_o3_cariolle

  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: pi
  USE mo_impl_constants      ,ONLY: min_rlcell_int, grf_bdywidth_c
  USE mo_run_config,          ONLY: nlev, io3
  USE mo_physical_constants,  ONLY:  amd, amo3
  USE mo_echam_phy_memory,    ONLY: t_echam_phy_field,     &
    &                               t_echam_phy_tend
  USE mtime,                  ONLY: datetime
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights, &
    &                                  calculate_time_interpolation_weights

  USE mo_lcariolle_types,     ONLY: t_avi, t_time_interpolation

  USE mo_parallel_config     ,ONLY: nproma
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_model_domain        ,ONLY: t_patch
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: interface_echam_o3_cariolle

CONTAINS

  !----------------------------------------------------------------
  !---------------------------------------------------------------------
  SUBROUTINE interface_echam_o3_cariolle(patch, rl_start, rl_end, field, tend, this_datetime)
    TYPE(t_patch)   ,INTENT(in), TARGET :: patch           !< grid/patch info
    INTEGER         ,INTENT(in) :: rl_start, rl_end
    TYPE(t_echam_phy_field),   POINTER :: field
    TYPE(t_echam_phy_tend) ,   POINTER :: tend
    TYPE(datetime),  POINTER     :: this_datetime  !< time step

     ! Temporary variables for Cariolle scheme (ozone)
    REAL(wp)    :: do3dt(nproma,nlev)
    TYPE(t_time_interpolation) :: time_interpolation
    EXTERNAL       lcariolle_lat_intp_li, lcariolle_pres_intp_li
    TYPE(t_time_interpolation_weights) :: current_time_interpolation_weights
    TYPE(t_avi) :: avi

    INTEGER  :: i_nchdom
    INTEGER  :: i_startblk,i_endblk
    INTEGER  :: jb             !< block index
    INTEGER  :: jcs, jce       !< start/end column index within this block
 
    i_nchdom   = MAX(1,patch%n_childdom)
    i_startblk = patch%cells%start_blk(rl_start,1)
    i_endblk   = patch%cells%end_blk(rl_end,i_nchdom)

    current_time_interpolation_weights = calculate_time_interpolation_weights(this_datetime)
    time_interpolation%imonth1=current_time_interpolation_weights%month1_index
    time_interpolation%imonth2=current_time_interpolation_weights%month2_index
    time_interpolation%weight1=current_time_interpolation_weights%weight1
    time_interpolation%weight2=current_time_interpolation_weights%weight2

!  NOTE: something is wrong with the avi dimensions; cannot be run in OpenMP
!$OMP PARALLEL PRIVATE(avi)
    ALLOCATE(avi%o3_vmr(nproma,nlev), avi%vmr2molm2(nproma,nlev), avi%cell_center_lat(nproma), &
        & avi%lday(nproma))
!$OMP DO PRIVATE(jcs,jce,do3dt)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb,i_startblk,i_endblk, jcs,jce, rl_start, rl_end)
 
      avi%ldown=.TRUE.
      avi%o3_vmr(jcs:jce,:)        = field%qtrc(jcs:jce,:,jb,io3)*amd/amo3
      avi%tmprt                    => field%ta(:,:,jb)
      avi%vmr2molm2(jcs:jce,:)     = field%mdry(jcs:jce,:,jb) / amd * 1.e3_wp
      avi%pres                     => field%presm_old(jcs:jce,:,jb)
      avi%cell_center_lat(jcs:jce) = patch%cells%center(jcs:jce,jb)%lat
      avi%lday(jcs:jce)            = field%cosmu0(jcs:jce,jb) > 1.e-3_wp

      CALL lcariolle_do3dt(                                                    &
          & jcs,                    jce,                nproma,                  &
          & nlev,                   time_interpolation, lcariolle_lat_intp_li,  &
          & lcariolle_pres_intp_li, avi,                do3dt                   )
      tend% qtrc(jcs:jce,:,jb,io3) = tend% qtrc(jcs:jce,:,jb,io3) + do3dt(jcs:jce,:)*amo3/amd
    ENDDO
!$OMP END DO 
    DEALLOCATE(avi%o3_vmr, avi%vmr2molm2, avi%cell_center_lat, avi%lday)
!$OMP END PARALLEL

  END SUBROUTINE interface_echam_o3_cariolle
  !---------------------------------------------------------------------

END MODULE mo_interface_echam_o3_cariolle
