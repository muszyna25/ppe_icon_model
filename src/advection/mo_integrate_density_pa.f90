!>
!! Integrates density for idealized test cases
!!
!! Integrates density for idealized test cases, i.e. when the dynamical core 
!! is switched off.
!!
!! <Don't forget references to literature.>
!!
!! @author Daniel Reinert, DWD
!!
!!
!! @par Revision History
!! Initial version by Daniel Reinert (2012-07-02)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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
MODULE mo_integrate_density_pa

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, min_rledge, min_rlcell,  &
    &                               MIURA, MIURA3
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e
  USE mo_advection_config,    ONLY: advection_config 
  USE mo_advection_hflux,     ONLY: upwind_hflux_miura, upwind_hflux_miura3
  USE mo_math_divrot,         ONLY: div 


  IMPLICIT NONE

  PRIVATE


  CHARACTER(len=*), PARAMETER :: version = '$Id$'


  PUBLIC :: get_nh_df_mflx_rho


CONTAINS

  !---------------------------------------------------------------------------
  !>
  !! SUBROUTINE get_nh_df_mflx_rho
  !!
  !! Short description:
  !! computes updated massflux and density for the deformational flow 
  !! test case, if necessary
  !!
  !!
  !! @par Revision History
  !! Initial revision by Daniel Reinert (2011-03-11)
  !!
  SUBROUTINE get_nh_df_mflx_rho( p_patch, p_int, p_prog_now, p_prog_new, &
    &                            p_metrics, p_diag, p_dtime )

    !INPUT PARAMETERS:
    TYPE(t_patch),      INTENT(IN)    :: p_patch
    TYPE(t_int_state),  INTENT(IN)    :: p_int
    TYPE(t_nh_prog),    INTENT(IN)    :: p_prog_now
    TYPE(t_nh_prog),    INTENT(INOUT) :: p_prog_new
    TYPE(t_nh_metrics), INTENT(IN)    :: p_metrics
    TYPE(t_nh_diag),    INTENT(INOUT) :: p_diag

    REAL(wp), INTENT(IN) :: p_dtime

    REAL(wp) ::   &                   !< density edge value
      &  z_rho_e(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) ::   &                   !< time averaged normal velocities
      &  z_vn_traj(nproma,p_patch%nlev,p_patch%nblks_e)
    REAL(wp) ::  &                    !< flux divergence at cell center
      &  z_fluxdiv_rho(nproma,p_patch%nlev,p_patch%nblks_c)

    INTEGER  :: nlev                !< number of full levels
    INTEGER  :: i_startblk, i_startidx, i_endblk, i_endidx
    INTEGER  :: i_rlstart, i_rlend, i_nchdom

    INTEGER  :: jb,je,jc,jk         !< loop indices

    LOGICAL  :: lcompute, lcleanup
    LOGICAL  :: lcoupled_rho        !< grid_sphere_radius-integrate mass continuity
                                    !< equation (.TRUE.)
    INTEGER  :: pid                 !< patch ID
    !---------------------------------------------------------------------------

    lcoupled_rho=.FALSE. ! grid_sphere_radius-integrate mass continuity equation (.TRUE.)

    ! get patch ID
    pid = p_patch%id

    ! number of vertical levels
    nlev = p_patch%nlev

    ! get new mass flux
    i_nchdom   = MAX(1,p_patch%n_childdom)

    lcompute =.TRUE.
    lcleanup =.TRUE.

    i_rlstart = 1
    i_rlend   = min_rledge

    i_startblk = p_patch%edges%start_blk(i_rlstart,1)
    i_endblk   = p_patch%edges%end_blk(i_rlend,i_nchdom)


    IF (lcoupled_rho) THEN
    ! explicit (grid_sphere_radius)-integration of mass continuity equation.

      !
      ! get time-averaged velocity field for mass flux and trajectory 
      ! computation.
!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, i_rlstart, i_rlend)
        DO jk = 1,nlev
          DO je = i_startidx, i_endidx

            z_vn_traj(je,jk,jb) = 0.5_wp * (p_prog_now%vn(je,jk,jb) &
              &                           + p_prog_new%vn(je,jk,jb)) 

          ENDDO
        ENDDO
      ENDDO
!!$OMP END DO
!!$OMP END PARALLEL

      !
      ! get new Density 'edge-value'
      !

      ! Select desired flux calculation method
      !
      ! Note: It is implicitly assumed, that only one tracer is advected.
      !       The namelist settings of the first tracer are applied to 
      !       the density as well.
      SELECT CASE( advection_config(pid)%ihadv_tracer(1) )
      CASE( MIURA )

        CALL upwind_hflux_miura(p_patch, p_prog_now%rho, p_prog_new%vn,        &
          &                     z_vn_traj, p_dtime, p_int, lcompute, lcleanup, &
          &                     advection_config(pid)%igrad_c_miura,           &
          &                     advection_config(pid)%itype_hlimit(1),         &
          &                     advection_config(pid)%iord_backtraj,           &
          &                     z_rho_e, opt_lout_edge=.TRUE.   )

      CASE( MIURA3 )

        CALL upwind_hflux_miura3(p_patch, p_prog_now%rho, p_prog_new%vn,       &
          &                      z_vn_traj, p_dtime, p_int, lcompute, lcleanup,&
          &                      advection_config(pid)%itype_hlimit(1),        &
          &                      z_rho_e, opt_lout_edge=.TRUE.                 )
      END SELECT


      !
      ! massflux at edges
      !
!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, i_rlstart, i_rlend)

        DO jk = 1,nlev
          DO je = i_startidx, i_endidx
            p_diag%mass_fl_e(je,jk,jb) = z_vn_traj(je,jk,jb)               & 
              &                         * p_metrics%ddqz_z_full_e(je,jk,jb) &
              &                         * z_rho_e(je,jk,jb)   
          ENDDO
        ENDDO
      ENDDO
!!$OMP END DO
!!$OMP END PARALLEL

      !
      ! Compute updated density field
      !
      CALL div( p_diag%mass_fl_e(:,:,:), p_patch, &! in
        &       p_int, z_fluxdiv_rho(:,:,:)       )! in,inout


      i_rlstart = 1
      i_rlend   = min_rlcell

      i_startblk = p_patch%cells%start_blk(i_rlstart,1)
      i_endblk   = p_patch%cells%end_blk(i_rlend,i_nchdom)

!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, i_rlstart, i_rlend)

        ! New density field
        DO jk = 1,nlev
          DO jc = i_startidx, i_endidx
            p_prog_new%rho(jc,jk,jb)= p_prog_now%rho(jc,jk,jb)          &
              &                    - p_dtime * z_fluxdiv_rho(jc,jk,jb)  &
              &                    * p_metrics%inv_ddqz_z_full(jc,jk,jb)

          ENDDO
        ENDDO
      ENDDO
!!$OMP END DO
!!$OMP END PARALLEL

    ELSE
    ! mass continuity equation will not be (grid_sphere_radius)-integrated.
    ! Only the updated mass flux is provided, assuming a homogeneuos,
    ! constant in time density field with \rho=1 everywhere.

      !
      ! massflux at edges
      !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx)
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, i_rlstart, i_rlend)

        DO jk = 1,nlev
          DO je = i_startidx, i_endidx
            p_diag%mass_fl_e(je,jk,jb) = 0.5_wp * (p_prog_now%vn(je,jk,jb) &
              &                        + p_prog_new%vn(je,jk,jb))          &
              &                        * p_metrics%ddqz_z_full_e(je,jk,jb)
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ENDIF

  END SUBROUTINE get_nh_df_mflx_rho



END MODULE mo_integrate_density_pa

