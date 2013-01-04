!>
!!  Subroutine to initialize the CBL case for HDCP2
!!
!!
!! @par Revision History
!! - first version by Anurag Dipankar , MPIM, (2012-12-12)
!! @par Literature
!! -
!!
!! @par Copyright
!! 2002-2008 by DWD and MPI-M
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
MODULE mo_nh_torus_exp
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2008
!
!-------------------------------------------------------------------------
!
!
!

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, min_rledge, min_rlcell
  USE mo_physical_constants,  ONLY: rd, cpd, p0ref, rd_o_cpd, grav
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_math_constants,      ONLY: pi, pi2
  USE mo_math_utilities,      ONLY: gnomonic_proj, t_geographical_coordinates, &
     &                              t_cartesian_coordinates, gc2cc, az_eqdist_proj
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_math_utilities,      ONLY: plane_torus_distance, arc_length
  USE mo_grid_config,         ONLY: grid_sphere_radius, is_plane_torus
  USE mo_sync,                ONLY: global_sum_array, sync_patch_array, SYNC_C, SYNC_E
  USE mo_nh_init_utils,       ONLY: init_w

   IMPLICIT NONE


  !DEFINED PARAMETERS:
  REAL(wp), PARAMETER :: zp0     = 100000._wp !< surface pressure
  REAL(wp), PARAMETER :: zt0     = 300._wp    !< temperature 

  REAL(wp):: waveno_x, waveno_y, u0, v0

  PUBLIC :: init_nh_state_prog_advtest, waveno_x, waveno_y, u0, v0

!--------------------------------------------------------------------

   CONTAINS
!-------------------------------------------------------------------------
!
!
  !>
  !! Initialization of prognostic state vector for the nh CBL test case 
  !!  without moisture
  !!
  !! @par Revision History
  !!
  !!
  SUBROUTINE init_nh_state_prog_advtest( ptr_patch, ptr_nh_prog, ptr_nh_diag,  &
    &                                   ptr_int, ptr_ext_data, p_metrics )

    ! INPUT PARAMETERS:
    TYPE(t_patch),TARGET,INTENT(IN) :: &  !< patch on which computation is performed
      &  ptr_patch
    TYPE(t_int_state),  INTENT(IN)  :: &
      &  ptr_int
    TYPE(t_nh_prog), INTENT(INOUT)  :: &  !< prognostic state vector
      &  ptr_nh_prog
    TYPE(t_nh_diag), INTENT(INOUT)  :: &  !< diagnostic state vector
      &  ptr_nh_diag
    TYPE(t_external_data), INTENT(INOUT) :: & !< external data
      &  ptr_ext_data
    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state

    TYPE(t_cartesian_coordinates) :: cc_cell, cc_center
  
    REAL(wp) :: z_dist, torus_len, torus_ht, rovcp, u_wind, v_wind, decay_width, &
                zvn1, zvn2, zscale_h

    REAL(wp), ALLOCATABLE :: z_sfc(:,:,:)

    INTEGER  :: jc,je,jk,jb,jt,i_startblk,i_startidx,i_endidx   !< loop indices
    INTEGER  :: nblks_c,npromz_c,nblks_e,npromz_e
    INTEGER  :: nlev, nlevp1                   !< number of full and half levels
    INTEGER  :: nlen, i_rcstartlev, jcn, jbn

  !-------------------------------------------------------------------------

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_int_c
    npromz_c = ptr_patch%npromz_int_c
    nblks_e  = ptr_patch%nblks_int_e
    npromz_e = ptr_patch%npromz_int_e

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    nlevp1 = ptr_patch%nlevp1

    ALLOCATE (z_sfc(nproma, 1, nblks_c))

    i_rcstartlev = 2
    i_startblk   = ptr_patch%edges%start_blk(i_rcstartlev,1)

    torus_len = ptr_patch%geometry_info%domain_length
    torus_ht  = ptr_patch%geometry_info%domain_height
    decay_width = min(torus_len,torus_ht)/10._wp

    zscale_h = rd*zt0/grav

    ! topography
    ptr_ext_data%atm%topography_c(:,:) = 0.0_wp

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = zp0

    ! init temperature
    ptr_nh_diag%temp(:,:,:)   = zt0

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen)
    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
 
        DO jk = 1, nlev

         ! compute full level pressure
         z_sfc(1:nlen,1,jb) = ptr_ext_data%atm%topography_c(1:nlen,jb)
         ptr_nh_diag%pres(1:nlen,jk,jb) = zp0                                     &
           &                  * exp(-(p_metrics%z_mc(1:nlen,jk,jb)-z_sfc(1:nlen,1,jb))/zscale_h)

         ! init virtual potential temperature
         ptr_nh_prog%theta_v(1:nlen,jk,jb) = zt0                                   &
           & * (p0ref/ptr_nh_diag%pres(1:nlen,jk,jb))**rd_o_cpd

         ! init density field rho
         ptr_nh_prog%rho(1:nlen,jk,jb) = ptr_nh_diag%pres(1:nlen,jk,jb)            &
           &                           / (rd * zt0)

         ! init rhotheta_v
         ptr_nh_prog%rhotheta_v(1:nlen,jk,jb) = ptr_nh_prog%rho(1:nlen,jk,jb)      &
           &                                  * ptr_nh_prog%theta_v(1:nlen,jk,jb)

         ! init exner pressure
         ptr_nh_prog%exner(1:nlen,jk,jb) = (ptr_nh_diag%pres(1:nlen,jk,jb)/p0ref)**rd_o_cpd

        END DO !jk
    ENDDO !jb
!$OMP END DO 
!$OMP END PARALLEL

 !First model level
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,je,i_startidx,i_endidx,jcn,jbn,cc_cell,z_dist, &
!$OMP            u_wind,v_wind,zvn1,zvn2) 
     DO jb = i_startblk, nblks_e

        CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e,  &
                           i_startidx, i_endidx, i_rcstartlev)

        DO je = i_startidx, i_endidx
 
          !First neighbour
          jcn  = ptr_patch%edges%cell_idx(je,jb,1)
          jbn  = ptr_patch%edges%cell_blk(je,jb,1)

          cc_center%x(:) = 0._wp
          cc_cell = ptr_patch%cells%cartesian_center(jcn,jbn) 
          z_dist  = plane_torus_distance(cc_center%x,cc_cell%x,ptr_patch%geometry_info)

          u_wind = u0 * (1._wp+exp(-(z_dist/decay_width)**2)) !* sin(waveno_x*cc_cell%x(1))**2)
          v_wind = v0 * (1._wp+exp(-(z_dist/decay_width)**2)) !* sin(waveno_y*cc_cell%x(2))**2)
                    
          zvn1   = u_wind * ptr_patch%edges%primal_normal_cell(je,jb,1)%v1 + &
                   v_wind * ptr_patch%edges%primal_normal_cell(je,jb,1)%v2

          !Second neighbour
          jcn  = ptr_patch%edges%cell_idx(je,jb,2)
          jbn  = ptr_patch%edges%cell_blk(je,jb,2)

          cc_center%x(:) = 0._wp
          cc_cell = ptr_patch%cells%cartesian_center(jcn,jbn) 
          z_dist  = plane_torus_distance(cc_center%x,cc_cell%x,ptr_patch%geometry_info)

          u_wind = u0 * (1._wp+exp(-(z_dist/decay_width)**2)) !* sin(waveno_x*cc_cell%x(1))**2)
          v_wind = v0 * (1._wp+exp(-(z_dist/decay_width)**2)) !* sin(waveno_y*cc_cell%x(2))**2)
                    
          zvn2   = u_wind * ptr_patch%edges%primal_normal_cell(je,jb,2)%v1 + &
                   v_wind * ptr_patch%edges%primal_normal_cell(je,jb,2)%v2

          ! compute normal wind component
          ptr_nh_prog%vn(je,nlev,jb) = ptr_int%c_lin_e(je,1,jb) * zvn1 + &
                                       ptr_int%c_lin_e(je,2,jb) * zvn2
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! copy to all vertical levels
    DO jk = nlev-1, 1, -1
      ptr_nh_prog%vn(:,jk,:) =  ptr_nh_prog%vn(:,nlev,:) !* real(max((jk-20)/10,0))
    END DO

    !
    ! initialize vertical velocity field
    !
    CALL init_w(ptr_patch, ptr_int, ptr_nh_prog%vn, p_metrics%z_ifc, ptr_nh_prog%w)
    CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)


  END SUBROUTINE init_nh_state_prog_advtest


!-------------------------------------------------------------------------
! 
  END MODULE mo_nh_torus_exp
