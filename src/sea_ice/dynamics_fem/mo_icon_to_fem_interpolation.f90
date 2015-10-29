!>
!! Contains the implementation of interpolation from ICON grid to FEM grid for the sea-ice dynamics
!!
!! @par Revision History
!! Based on: *list modules which were used for this development*
!! Developed  by Einar (2014), Vladimir (2015)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_icon_to_fem_interpolation
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: min_rlcell, min_rledge, min_rlvert
!  USE mo_grid_config,         ONLY: l_limited_area
  USE mo_model_domain,        ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: ltimer
  USE mo_loopindices,         ONLY: get_indices_v!, get_indices_c, get_indices_e
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_intp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cells2verts_scalar_seaice

CONTAINS

!-----------------------------------------------------------------------
!
!  ! averaging and interpolation routines and
!  ! routines needed to compute the coefficients therein
!
!-----------------------------------------------------------------------

!------------------------------------------------------------------------
!>
!!  Computes  average of scalar fields from centers of cells to vertices.
!!
!!
!! @par Revision History
!! Developed  by Almut Gassmann, MPI-M (2009-01-28)
!!
SUBROUTINE cells2verts_scalar_seaice( p_cell_in, ptr_patch, c_int, p_vert_out,  &
  &                            opt_slev, opt_elev, opt_rlstart, opt_rlend )
!

TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! cell based scalar input field
REAL(wp), INTENT(in) :: p_cell_in(:,:,:)   ! dim: (nproma,nlev,nblks_c)

! coefficients for interpolation
REAL(wp), INTENT(in) :: c_int(:,:,:)       ! dim: (nproma,9-cell_type,nblks_v)

INTEGER, INTENT(in), OPTIONAL :: opt_slev  ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL :: opt_elev  ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL ::  opt_rlstart, opt_rlend

! vertex based scalar output field
REAL(wp), INTENT(inout) :: p_vert_out(:,:,:) ! dim: (nproma,nlev,nblks_v)

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jv, jk, jb, ji
INTEGER :: rl_start, rl_end
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER :: cell_index, cell_block

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

!-----------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(p_cell_in,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlvert
END IF

iidx => ptr_patch%verts%cell_idx
iblk => ptr_patch%verts%cell_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%verts%start_blk(rl_start,1)
i_endblk   = ptr_patch%verts%end_blk(rl_end,i_nchdom)


IF (ltimer) CALL timer_start(timer_intp)

IF (ptr_patch%geometry_info%cell_type == 6) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)


#ifdef __LOOP_EXCHANGE
    DO jv = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=6
    DO jk = slev, elev
      DO jv = i_startidx, i_endidx
#endif

         p_vert_out(jv,jk,jb) =                       &
           c_int(jv,1,jb) * p_cell_in(iidx(jv,jb,1),jk,iblk(jv,jb,1)) + &
           c_int(jv,2,jb) * p_cell_in(iidx(jv,jb,2),jk,iblk(jv,jb,2)) + &
           c_int(jv,3,jb) * p_cell_in(iidx(jv,jb,3),jk,iblk(jv,jb,3))

      ENDDO
    ENDDO

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
ELSE IF (ptr_patch%geometry_info%cell_type == 3) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jv,jk) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_v(ptr_patch, jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, rl_start, rl_end)


#ifdef __LOOP_EXCHANGE
    DO jv = i_startidx, i_endidx
      DO jk = slev, elev
#else
!CDIR UNROLL=6
    DO jk = slev, elev
      DO jv = i_startidx, i_endidx
#endif

    ! Sea-ice module calls this funciton to interpolate from ICON grid to FEM P1P1 grid
    ! An indexing error used to occur because of zeros in iidx(jv,jb,:) at the pentagons centers and boundary vertices

    ! In the atmosphere module this is done by calling this function somewhere at the initialization step
    !    CALL move_dummies_to_end_idxblk( &
    !        ptr_patch%verts%cell_idx(:,:,1:max_verts_connectivity), &
    !        ptr_patch%n_patch_verts, max_verts_connectivity, &
    !        use_duplicated_connectivity)

!    if ( (iidx(jv,jb,6)==0) .or. (iblk(jv,jb,6)==0) ) then
!        print *, 'lat, lon =', ptr_patch%verts%vertex(jv,jb)
!    !    print *, '6 interp coeff =', c_int(jv,1:6,jb)
!        print *, 'iidx(jv,jb,:)', iidx(jv,jb,:)
!        print *, 'iblk(jv,jb,:)', iblk(jv,jb,:)
!        print *, 'max_connectivity', ptr_patch%verts%max_connectivity
!        print *, 'ptr_patch%verts%num_edges(jv,jb)', ptr_patch%verts%num_edges(jv,jb)
!    endif

        p_vert_out(jv,jk,jb) = 0.0_wp

        DO ji = 1, 6
          cell_index = iidx(jv,jb,ji)
          cell_block = iblk(jv,jb,ji)
          IF (cell_index > 0)                                      &
            & p_vert_out(jv,jk,jb) = p_vert_out(jv,jk,jb) +                                        &
            &      c_int(jv,ji,jb) * p_cell_in(cell_index,jk,cell_block)
        ENDDO

!       Old version of the code that relies on non-zero values in iidx => ptr_patch%verts%cell_idx
!         p_vert_out(jv,jk,jb) =                       &
!           c_int(jv,1,jb) * p_cell_in(iidx(jv,jb,1),jk,iblk(jv,jb,1)) + &
!           c_int(jv,2,jb) * p_cell_in(iidx(jv,jb,2),jk,iblk(jv,jb,2)) + &
!           c_int(jv,3,jb) * p_cell_in(iidx(jv,jb,3),jk,iblk(jv,jb,3)) + &
!           c_int(jv,4,jb) * p_cell_in(iidx(jv,jb,4),jk,iblk(jv,jb,4)) + &
!           c_int(jv,5,jb) * p_cell_in(iidx(jv,jb,5),jk,iblk(jv,jb,5)) + &
!           c_int(jv,6,jb) * p_cell_in(iidx(jv,jb,6),jk,iblk(jv,jb,6))

      ENDDO
    ENDDO

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
ENDIF

IF (ltimer) CALL timer_stop(timer_intp)


END SUBROUTINE cells2verts_scalar_seaice
!------------------------------------------------------------------------

END MODULE mo_icon_to_fem_interpolation
