!>
!! Updates dynamical fields with slow physics tendencies in case ldynamics=F
!!
!! Updates dynamical fields with slow physics tendencies in the special case
!! that dynamics are switched off.
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial revision by Daniel Reinert, DWD (2013-11-28)
!!
!! @par Copyright
!! 2002-2013 by DWD and MPI-M
!! This software is provided for non-commercial use only.
!! See the LICENSE and the WARRANTY conditions.
!!
!! @par License
!! The use of ICON is hereby granted free of charge for an unlimited time,
!! provided the following rules are accepted and applied:
!! <ol>
!! <li> You may use or modify this code for your own non commercial and non
!!    violent purposes.
!! <li> The code may not be re-distributed without the consent of the authors.
!! <li> The copyright notice and statement of authorship must appear in all
!!    copies.
!! <li> You accept the warranty conditions (see WARRANTY).
!! <li> In case you intend to use the code commercially, we oblige you to sign
!!    an according license agreement with DWD and MPI-M.
!! </ol>
!!
!! @par Warranty
!! This code has been tested up to a certain level. Defects and weaknesses,
!! which may be included in the code, do not establish any warranties by the
!! authors.
!! The authors do not make any warranty, express or implied, or assume any
!! liability or responsibility for the use, acquisition or application of this
!! software.
!!
MODULE mo_update_dyn

  USE mo_kind,               ONLY: wp
  USE mo_loopindices,        ONLY: get_indices_c, get_indices_e
  USE mo_model_domain,       ONLY: t_patch
  USE mo_nonhydro_types,     ONLY: t_nh_state
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_physical_constants, ONLY: p0ref, rd, cvd_o_rd
  USE mo_impl_constants,     ONLY: min_rlcell_int, min_rledge_int
  USE mo_sync,               ONLY: SYNC_E, SYNC_C, sync_patch_array

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: add_slowphys

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


CONTAINS

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Updates dynamical fields with slow physics tendencies
  !!
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2013-11-28)
  !! 
  !!
  SUBROUTINE add_slowphys(p_nh, p_patch, nnow, nnew, dtime)

    TYPE(t_nh_state),  TARGET, INTENT(INOUT) :: p_nh
    TYPE(t_patch),     TARGET, INTENT(IN)    :: p_patch

    ! Time levels
    INTEGER,                   INTENT(IN)    :: nnow, nnew
    ! Time step
    REAL(wp),                  INTENT(IN)    :: dtime

    INTEGER :: nlev                 ! number of vertical (full) levels

    INTEGER :: jc, je, jk, jb        ! loop indices
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk  ! start and end block
    INTEGER :: i_startidx, i_endidx  ! start and end indices
    INTEGER :: i_nchdom

   !-----------------------------------------------------------------------

    ! number of vertical levels
    nlev   = p_patch%nlev

    i_nchdom = MAX(1,p_patch%n_childdom)

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
    rl_start = 3
    rl_end = min_rlcell_int - 1

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx

          ! rho (simply copy)
          p_nh%prog(nnew)%rho(jc,jk,jb) = p_nh%prog(nnow)%rho(jc,jk,jb)

          ! w (simply copy)
          p_nh%prog(nnew)%w(jc,jk,jb) = p_nh%prog(nnow)%w(jc,jk,jb)

          ! exner (update)
          p_nh%prog(nnew)%exner(jc,jk,jb) = p_nh%prog(nnow)%exner(jc,jk,jb)  &
             &                            + dtime*p_nh%diag%ddt_exner_phy(jc,jk,jb)

          ! diagnose theta_v from updated exner
          p_nh%prog(nnew)%theta_v(jc,jk,jb) = (p0ref/(rd*p_nh%prog(nnew)%rho(jc,jk,jb)))  &
             &                              * p_nh%prog(nnew)%exner(jc,jk,jb)**cvd_o_rd

        ENDDO  ! jc
      ENDDO  ! jk
    ENDDO  ! jb
!$OMP ENDDO NOWAIT


    rl_start = 7
    rl_end   = min_rledge_int-1

    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_endblk   = p_patch%edges%end_blk  (rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)

      DO jk = 1, nlev
        DO je = i_startidx, i_endidx

          ! vn (update)
          p_nh%prog(nnew)%vn(je,jk,jb) = p_nh%prog(nnow)%vn(je,jk,jb)  &
            &                          + dtime * p_nh%diag%ddt_vn_phy(je,jk,jb)
        ENDDO  ! je
      ENDDO  ! jk
    ENDDO  ! jb
!$OMP ENDDO
!$OMP END PARALLEL


    ! Synchronize updated prognostic variables
    CALL sync_patch_array(SYNC_C,p_patch,p_nh%prog(nnew)%exner)
    CALL sync_patch_array(SYNC_E,p_patch,p_nh%prog(nnew)%vn)

  END SUBROUTINE add_slowphys 


END MODULE mo_update_dyn



