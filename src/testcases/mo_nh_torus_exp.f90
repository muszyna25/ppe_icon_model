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
  REAL(wp), PARAMETER :: zh0     = 1000._wp   !< height (m) above which temperature increases
  REAL(wp), PARAMETER :: dtdz    = 0.003_wp   !< lapse rate

  REAL(wp):: sst_cbl
  LOGICAL :: is_dry_cbl

  PUBLIC :: init_nh_state_cbl, sst_cbl, is_dry_cbl

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
  SUBROUTINE init_nh_state_cbl( ptr_patch, ptr_nh_prog, ptr_nh_diag,  &
    &                           ptr_int, ptr_ext_data, p_metrics )

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

    REAL(wp) :: rho_sfc
    INTEGER  :: jc,je,jk,jb,jt,i_startblk,i_startidx,i_endidx   !< loop indices
    INTEGER  :: nblks_c,npromz_c,nblks_e,npromz_e
    INTEGER  :: nlev, nlevp1                   !< number of full and half levels
    INTEGER  :: nlen, i_rcstartlev

  !-------------------------------------------------------------------------

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_int_c
    npromz_c = ptr_patch%npromz_int_c
    nblks_e  = ptr_patch%nblks_int_e
    npromz_e = ptr_patch%npromz_int_e

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    nlevp1 = ptr_patch%nlevp1

    i_rcstartlev = 2
    i_startblk   = ptr_patch%edges%start_blk(i_rcstartlev,1)

    ! topography
    ptr_ext_data%atm%topography_c(:,:) = 0.0_wp

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = zp0

    rho_sfc = zp0 / (rd * sst_cbl)

    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
 
        DO jk = 1, nlev

         ! init virtual potential temperature: keep it lower than sst
         ! to have upward heat flux for CBL
         ptr_nh_prog%theta_v(1:nlen,jk,jb) = sst_cbl - 2._wp + &
                     max(0._wp, (p_metrics%z_mc(1:nlen,jk,jb)-zh0)*dtdz)

         ! init density field rho
         ptr_nh_prog%rho(1:nlen,jk,jb) = rho_sfc

         ! compute full level pressure
         ptr_nh_diag%pres(1:nlen,jk,jb) = ptr_nh_prog%rho(1:nlen,jk,jb) * rd * &
                                          ptr_nh_prog%theta_v(1:nlen,jk,jb) 
         ! init rhotheta_v
         ptr_nh_prog%rhotheta_v(1:nlen,jk,jb) = ptr_nh_prog%rho(1:nlen,jk,jb)      &
           &                                  * ptr_nh_prog%theta_v(1:nlen,jk,jb)

         ! init exner pressure
         ptr_nh_prog%exner(1:nlen,jk,jb) = (ptr_nh_diag%pres(1:nlen,jk,jb)/p0ref)**rd_o_cpd

        END DO !jk
    ENDDO !jb

    !Zero mean wind to avoid Croiolis force and hence mean pressure gradient
    ptr_nh_prog%vn = 0._wp
    ptr_nh_prog%w = 0._wp

    IF(is_dry_cbl)THEN 
      ! Tracers set to 0 for DRY case
      ptr_nh_prog%tracer =  0._wp
    ELSE
      CALL finish('mo_nh_torus_exp: init_nh_state_cbl:',  &
        &      '   moist CBL not implemented yet- Stopping!')
    END IF 

  END SUBROUTINE init_nh_state_cbl


!-------------------------------------------------------------------------
! 
  END MODULE mo_nh_torus_exp
