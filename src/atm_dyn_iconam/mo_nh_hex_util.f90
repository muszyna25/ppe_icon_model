!>
!! Initializes and controls the time stepping in the nonhydrostatic model.
!!
!!
!! @par Revision History
!! Initial release by Almut Gassmann, MPI-M (2009-02-06)
!!
!! @par Copyright
!! 2002-2009 by DWD and MPI-M
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
!!
MODULE mo_nh_hex_util

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_parallel_config,     ONLY: nproma
  USE mo_nonhydrostatic_config,ONLY: ltheta_up_vert
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_state,      ONLY: t_nh_state, t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_exception,           ONLY: message, finish
  USE mo_physical_constants,  ONLY: cvd_o_rd
  USE mo_interpolation,       ONLY: t_int_state, &
                                    verts2edges_scalar,&
                                    cells2verts_scalar,&
                                    cells2edges_scalar, &
                                    sick_a, sick_o
  USE mo_vector_operations,   ONLY: covariant_velocities,       &
                                    contravariant_vorticities,  &
                                    contravariant_vert_mass_flux,&
                                    orthogonal_vorticities,     &
                                    vorticity_tendencies,       &
                                    kinetic_energy,             &
                                    impl_vert_adv_vn
  USE mo_divergent_modes,     ONLY: impl_vert_adv_theta, impl_vert_adv_theta_5diag
  USE mo_math_operators,      ONLY: nabla2_scalar, nabla2_vec
  USE mo_sync,                ONLY: SYNC_C, sync_patch_array

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: forcing_straka, momentum_adv

CONTAINS

  !-----------------------------------------------------------------------------
  !>
  !! forcing_straka
  !!
  !! Computes the tendency from an (artificial) forcing term.
  !! Currently that contains the Straka test case.
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-10.16)
  !!
  SUBROUTINE  forcing_straka(p_nh_prog,p_patch,p_int_state,p_metrics,p_nh_diag)

    TYPE(t_patch),TARGET, INTENT(in) :: p_patch    !< single patch
    TYPE(t_int_state),INTENT(in)  :: p_int_state!< single interpolation state
    TYPE(t_nh_metrics),INTENT(in) :: p_metrics  !< single metrics state
    TYPE(t_nh_prog), INTENT(in)   :: p_nh_prog  !< single nh prognostic state
    TYPE(t_nh_diag), INTENT(inout):: p_nh_diag  !< single nh diagnostic state

    REAL(wp), PARAMETER :: z_ny=75.0_wp !< Straka test diffusion coefficient
    REAL(wp) :: z_lapl_theta(nproma,p_patch%nlev  ,p_patch%nblks_c),  &
                z_lapl_w    (nproma,p_patch%nlevp1,p_patch%nblks_c),  &
                z_lapl_vn   (nproma,p_patch%nlev  ,p_patch%nblks_e)
    INTEGER :: jb, jk, nlen, nblks_e, npromz_e, nblks_c, npromz_c, &
               kp1, km1
    INTEGER :: nlev, nlevp1         !< number of full and half levels
    !--------------------------------------------------------------------------

    nblks_c   = p_patch%nblks_int_c
    npromz_c  = p_patch%npromz_int_c
    nblks_e   = p_patch%nblks_int_e
    npromz_e  = p_patch%npromz_int_e

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! Horizontal part
    CALL nabla2_scalar(p_nh_prog%theta_v, p_patch, p_int_state, z_lapl_theta)
    CALL nabla2_scalar(p_nh_prog%w, p_patch, p_int_state, z_lapl_w, 2, nlev)
    CALL nabla2_vec(p_nh_prog%vn, p_patch, p_int_state, z_lapl_vn)

    ! Vertical part
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk,kp1,km1)
    DO jb = 1, nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      DO jk = 1, nlev
        kp1=MIN(jk+1,nlev)
        km1=MAX(jk-1,1)
        p_nh_diag%ddt_vn(1:nlen,jk,jb) = z_ny*(z_lapl_vn(1:nlen,jk,jb)+ &
          (p_nh_prog%vn(1:nlen,km1,jb)-2.0_wp*p_nh_prog%vn(1:nlen,jk,jb)&
           +p_nh_prog%vn(1:nlen,kp1,jb))&
           /(p_metrics%ddqz_z_full_e(1:nlen,jk,jb)**2))
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, jk,kp1,km1)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = npromz_c
      ENDIF
      DO jk = 1, nlev
        kp1=MIN(jk+1,nlev)
        km1=MAX(jk-1,1)
        p_nh_diag%ddt_exner(1:nlen,jk,jb) = z_ny*&
          (z_lapl_theta(1:nlen,jk,jb)+(p_nh_prog%theta_v(1:nlen,km1,jb)&
          -2.0_wp*p_nh_prog%theta_v(1:nlen,jk,jb)&
                 +p_nh_prog%theta_v(1:nlen,kp1,jb))&
          /(p_metrics%ddqz_z_full(1:nlen,jk,jb)**2))&
          /cvd_o_rd*p_nh_prog%exner(1:nlen,jk,jb)&
          /p_nh_prog%theta_v(1:nlen,jk,jb)
      ENDDO
      DO jk = 2, nlev ! bottom and top tendencies do not exist
        kp1=MIN(jk+1,nlevp1)
        km1=MAX(jk-1,1)
        p_nh_diag%ddt_w(1:nlen,jk,jb) = z_ny*&
          (z_lapl_w(1:nlen,jk,jb)+(p_nh_prog%w(1:nlen,km1,jb)&
          -2.0_wp*p_nh_prog%w(1:nlen,jk,jb)+p_nh_prog%w(1:nlen,kp1,jb))&
          /(p_metrics%ddqz_z_half(1:nlen,jk,jb)**2))
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE forcing_straka
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  !>
  !! momentum_adv
  !!
  !! Computes the tendency from the velocity advection
  !!
  !! @par Revision History
  !! Initial release by Almut Gassmann, MPI-M (2009-04.17)
  !! Modification by Almut Gassman, MPI-M (2010-01.11)
  !! - adjustment to new time stepping scheme
  !! Modification by Almut Gassman, MPI-M (2011-01)
  !! - RK2 stepping and implicit vertical advection
  !! Complete rewriting by Almut Gassmann, MPI-M (2011-04)
  !! - avoid double evaluations in the dynamical core, improved consistencies
  !!
  SUBROUTINE momentum_adv(p_nh,p_patch,p_int,know,knew,l_predictor)

    IMPLICIT NONE

    TYPE(t_patch), TARGET, INTENT(in)    :: p_patch    !< single patch
    TYPE(t_int_state),TARGET,INTENT(in)  :: p_int      !< single interpolation state
    TYPE(t_nh_state),TARGET,INTENT(inout):: p_nh       !< nh state
    INTEGER, INTENT(in)                  :: know, knew !< time levels involved
    !INTEGER, INTENT(in)                  :: k_call_no  !< counter for time stepping state
    LOGICAL, INTENT(in)                  :: l_predictor  !< counter for time stepping state

    REAL(wp),POINTER :: p_w(:,:,:)
    REAL(wp),POINTER :: p_vn(:,:,:)
    REAL(wp),POINTER :: p_rho(:,:,:)
    REAL(wp),POINTER :: p_rho_e(:,:,:)
    REAL(wp):: z_vn_impl    (nproma,p_patch%nlev  ,p_patch%nblks_e), &
    &          z_rho_l      (nproma,p_patch%nlevp1,p_patch%nblks_c), &
    &          z_rho_e_l    (nproma,p_patch%nlevp1,p_patch%nblks_e), &
    &          z_rho_e_lnt  (nproma,p_patch%nlevp1,p_patch%nblks_e), &
    &          z_rho_v      (nproma,p_patch%nlev  ,p_patch%nblks_v), &
    &          z_rho_a      (nproma,p_patch%nlev  ,p_patch%nblks_e), &
    &          z_wmflx_ort  (nproma,p_patch%nlevp1,p_patch%nblks_c), &
    &          z_wmflx_ort_e(nproma,p_patch%nlevp1,p_patch%nblks_e), &
    &          z_wmflx_con  (nproma,p_patch%nlevp1,p_patch%nblks_c), &
    &          z_wmflx_con_e(nproma,p_patch%nlevp1,p_patch%nblks_e), &
    &          z_tmp_v_l    (nproma,p_patch%nlevp1,p_patch%nblks_v), &
    &          z_tmp_e_l    (nproma,p_patch%nlevp1,p_patch%nblks_e), &
    &          z_theta_v_impl(nproma,p_patch%nlev ,p_patch%nblks_c)
    INTEGER :: nblks_c, npromz_c, nblks_e, npromz_e, nlen, jk, jb, nlev, nlevp1

    !--------------------------------------------------------------------------

    nblks_c  = p_patch%nblks_int_c
    npromz_c = p_patch%npromz_int_c
    nblks_e  = p_patch%nblks_int_e
    npromz_e = p_patch%npromz_int_e
    nlev     = p_patch%nlev
    nlevp1   = p_patch%nlevp1

    IF (l_predictor) THEN ! first call

      ! density at edges
      CALL cells2edges_scalar(p_nh%prog(know)%rho,p_patch,p_int%c_lin_e,p_nh%diag%rho_e)
      CALL cells2verts_scalar(p_nh%prog(know)%rho,p_patch,p_int%cells_aw_verts,z_rho_v)
      CALL verts2edges_scalar(z_rho_v,p_patch,p_int%v_1o2_e,z_rho_a)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen)
      DO jb = 1,nblks_e
        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF
        DO jk = 1,nlev
          p_nh%diag%rho_e(1:nlen,jk,jb) = sick_a*z_rho_a(1:nlen,jk,jb)&
          &                      +sick_o*p_nh%diag%rho_e(1:nlen,jk,jb)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      p_rho_e => p_nh%diag%rho_e
      p_vn    => p_nh%prog(know)%vn
      p_rho   => p_nh%prog(know)%rho
      p_w     => p_nh%prog(know)%w

    ELSEIF(.NOT.l_predictor) THEN ! second call

      ! density at edges
      CALL cells2edges_scalar(p_nh%prog(knew)%rho,p_patch,p_int%c_lin_e,p_nh%diag%rho_star_e)
      CALL cells2verts_scalar(p_nh%prog(knew)%rho,p_patch,p_int%cells_aw_verts,z_rho_v)
      CALL verts2edges_scalar(z_rho_v,p_patch,p_int%v_1o2_e,z_rho_a)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen)
      DO jb = 1,nblks_e
        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF
        DO jk = 1,nlev
          p_nh%diag%rho_star_e(1:nlen,jk,jb) = sick_a*z_rho_a(1:nlen,jk,jb)&
          &                      +sick_o*p_nh%diag%rho_star_e(1:nlen,jk,jb)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      p_rho_e => p_nh%diag%rho_star_e
      p_vn    => p_nh%prog(knew)%vn
      p_rho   => p_nh%prog(knew)%rho
      p_w     => p_nh%prog(knew)%w

    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk)
    DO jb = 1,p_patch%nblks_c
      ! density at half levels
      z_rho_l(1:nproma,1,jb)     = p_rho(1:nproma,1,jb)
      DO jk=2,nlev
        z_rho_l(1:nproma,jk,jb)  = 0.5_wp*(p_rho(1:nproma,jk-1,jb)+p_rho(1:nproma,jk,jb))
      ENDDO
      z_rho_l(1:nproma,nlevp1,jb)= p_rho(1:nproma,nlev,jb)
      ! orthogonal vertical mass flux
      DO jk = 1, nlevp1
        z_wmflx_ort(1:nproma,jk,jb)= p_w(1:nproma,jk,jb)*z_rho_l(1:nproma,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk)
    DO jb = 1,p_patch%nblks_e
      ! density at half levels and edges with tilde
      z_rho_e_l(1:nproma,1,jb)     = p_rho_e(1:nproma,1,jb)
      DO jk=2,nlev
        z_rho_e_l(1:nproma,jk,jb)  = 0.5_wp*(p_rho_e(1:nproma,jk-1,jb)+p_rho_e(1:nproma,jk,jb))
      ENDDO
      z_rho_e_l(1:nproma,nlevp1,jb)= p_rho_e(1:nproma,nlev,jb)
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! unaveraged rho at half level edges (needed for horiz. pot. vorticity)
    CALL cells2edges_scalar(z_rho_l,p_patch,p_int%c_lin_e,z_rho_e_lnt,1,nlevp1)

    ! contravariant vertical mass flux (does not use the tilde averages rho*w)
    CALL contravariant_vert_mass_flux(p_vn,p_rho_e,z_wmflx_ort,p_patch,p_int,&
    &                                 p_nh%metrics,z_wmflx_con,.TRUE.)
    CALL sync_patch_array(SYNC_C,p_patch,z_wmflx_con)

    ! implicit vertical advection for theta with the above contrav. vert. mass flux
    IF (ltheta_up_vert) THEN
      CALL impl_vert_adv_theta_5diag(p_nh%prog(know)%theta_v,z_wmflx_con,p_rho,p_patch,&
      &                        p_nh%metrics,z_theta_v_impl)
    ELSE
      CALL impl_vert_adv_theta(p_nh%prog(know)%theta_v,z_wmflx_con,p_rho,p_patch,&
      &                        p_nh%metrics,z_theta_v_impl)
    ENDIF
    IF (l_predictor) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen)
      DO jb = 1,nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        p_nh%diag%theta_v_impl(1:nlen,:,jb)=z_theta_v_impl(1:nlen,:,jb)
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ELSE
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen)
      DO jb = 1,nblks_c
        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF
        p_nh%diag%theta_v_impl(1:nlen,:,jb) = &
          & 0.5_wp*( p_nh%diag%theta_v_impl(1:nlen,:,jb)  + z_theta_v_impl(1:nlen,:,jb) )
        p_nh%diag%theta_v_ave(1:nlen,:,jb)  = &
          & 0.5_wp*( p_nh%prog(know)%theta_v(1:nlen,:,jb) +p_nh%prog(knew)%theta_v(1:nlen,:,jb) )
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ENDIF

    ! vertical mass flux at edges
    CALL cells2edges_scalar(z_wmflx_con,p_patch,p_int%c_lin_e,z_wmflx_con_e,1,nlevp1)
    CALL cells2edges_scalar(z_wmflx_ort,p_patch,p_int%c_lin_e,z_wmflx_ort_e,1,nlevp1)

    ! contravariant vertical mass flux (does use the tilde averages rho*w)
    CALL cells2verts_scalar(z_wmflx_ort,p_patch,p_int%cells_aw_verts,z_tmp_v_l,1,nlevp1)
    CALL verts2edges_scalar(z_tmp_v_l,p_patch,p_int%v_1o2_e,z_tmp_e_l,1,nlevp1)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen)
    DO jb = 1,nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      DO jk = 2, nlev
        ! here is only the metric correction stored which is not averaged
        z_wmflx_con_e(1:nlen,jk,jb) = z_wmflx_con_e(1:nlen,jk,jb)-z_wmflx_ort_e(1:nlen,jk,jb)
      ENDDO
      DO jk = 1, nlevp1
        ! now the vertical flux is averaged
        z_wmflx_ort_e(1:nlen,jk,jb) = sick_a*z_tmp_e_l(1:nlen,jk,jb) &
        &                            +sick_o*z_wmflx_ort_e(1:nlen,jk,jb)
      ENDDO
      DO jk = 2, nlev
        ! the full contravariant values is there again
        z_wmflx_con_e(1:nlen,jk,jb) = z_wmflx_con_e(1:nlen,jk,jb)+z_wmflx_ort_e(1:nlen,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen)
    DO jb = 1,nblks_e
      IF (jb /= nblks_e) THEN
        nlen = nproma
      ELSE
        nlen = npromz_e
      ENDIF
      DO jk = 2, nlev
        ! the rho cancel here (from horiz. pot. vorticity)
        ! z_rho_e_l has not tilde e average
        z_wmflx_con_e(1:nlen,jk,jb) = z_wmflx_con_e(1:nlen,jk,jb)/z_rho_e_lnt(1:nlen,jk,jb)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! implicit vertical advection for horizontal velocity
    CALL impl_vert_adv_vn(p_nh%prog(know)%vn,z_wmflx_con_e,p_patch,p_nh%metrics,z_vn_impl)

    ! covariant velocities
    CALL covariant_velocities(p_vn,p_w,p_patch,p_int,p_nh%metrics,p_nh%diag)

    ! contravariant vorticities
    CALL contravariant_vorticities(p_vn,z_vn_impl,p_patch,p_int,p_nh%metrics,p_nh%diag)

    ! orthogonal vorticities
    CALL orthogonal_vorticities(p_patch,p_int,p_nh%metrics,p_nh%diag)

    ! vorticity flux term
    CALL vorticity_tendencies(p_vn,z_vn_impl,z_wmflx_ort_e,p_rho,p_rho_e,z_rho_e_l,&
    &                         p_patch,p_int,p_nh%metrics,p_nh%diag,l_predictor)

    ! kinetic energy
    CALL kinetic_energy(p_vn,z_vn_impl,p_w,p_patch,p_int,p_nh%metrics,p_nh%diag,2,&
    & l_predictor)

  END SUBROUTINE momentum_adv
  


END MODULE mo_nh_hex_util


