!>
!!  Held-Suarez test for the NH-Core
!!
!!
!! @par Revision History
!! - first version by P. Ripodas , DWD, (2010-09)
!!
!! @par Literature
!! - Held, I. M. and Suarez, M. J. (1994): A Proposal for the Intercomparison
!!   of the Dynamical Cores of Atmospheric General Circulation Models
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
MODULE mo_nh_hs_test
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
USE mo_physical_constants,  ONLY: rd, rd_o_cpd, p0ref, grav
USE mo_model_domain,        ONLY: t_patch
USE mo_ext_data_types,      ONLY: t_external_data
USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
USE mo_parallel_config,     ONLY: nproma


IMPLICIT NONE

PRIVATE

REAL(wp), PARAMETER :: zp_hs      = 100000._wp            !< surface pressure
REAL(wp), PARAMETER :: zt_hs      = 300._wp               !< atmospheric temperature
!REAL(wp), PARAMETER :: zt_hs      = 250._wp               !< atmospheric temperature


PUBLIC :: init_nh_state_prog_held_suarez

!--------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!
  !>
  !!               Initialization of prognostic state vector.
  !!
  !!
  !! @par Revision History
  !!
  !!
  SUBROUTINE init_nh_state_prog_held_suarez( ptr_patch, ptr_nh_prog, ptr_nh_diag, &
    &                                   ptr_ext_data, p_metrics )

    TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
      &  ptr_nh_prog

    TYPE(t_nh_diag), INTENT(INOUT)      :: &  !< diagnostic state vector
      &  ptr_nh_diag

    TYPE(t_external_data), INTENT(INOUT):: &  !< external data
      &  ptr_ext_data

    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state


    INTEGER               :: nblks_c, npromz_c, nlen
    INTEGER               :: nlev                 !< number of full levels
    INTEGER               :: jk, jb  ! loop variables

    REAL(wp)              :: zscale_h             !< initialized variables
    REAL(wp), ALLOCATABLE :: z_sfc(:,:,:)



!--------------------------------------------------------------------
!
    zscale_h = rd*zt_hs/grav

    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c

    nlev = ptr_patch%nlev

    ALLOCATE (z_sfc(nproma, 1, ptr_patch%nblks_c))

    ! topography
    ptr_ext_data%atm%topography_c(:,:) = 0.0_wp

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = zp_hs

    ! init temperature
    ptr_nh_diag%temp(:,:,:)   = zt_hs


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
        ptr_nh_diag%pres(1:nlen,jk,jb) = zp_hs                                      &
          &                  * exp(-(p_metrics%z_mc(1:nlen,jk,jb)-z_sfc(1:nlen,1,jb))/zscale_h)

        ! init virtual potential temperature
        ptr_nh_prog%theta_v(1:nlen,jk,jb) = zt_hs                                   &
          & * (p0ref/ptr_nh_diag%pres(1:nlen,jk,jb))**rd_o_cpd

        ! init density field rho
        ptr_nh_prog%rho(1:nlen,jk,jb) = ptr_nh_diag%pres(1:nlen,jk,jb)            &
          &                           / (rd * zt_hs)

        ! init exner pressure
        ptr_nh_prog%exner(1:nlen,jk,jb) = (ptr_nh_diag%pres(1:nlen,jk,jb)/p0ref)**rd_o_cpd

      ENDDO !jk
    ENDDO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    !
    ! initialize horizontal velocity field
    !

    ptr_nh_prog%vn(:,:,:) =     0.0_wp

    !
    ! initialize vertical velocity field
    !

    ptr_nh_prog%w(:,:,:) =      0.0_wp


  END SUBROUTINE init_nh_state_prog_held_suarez


!--------------------------------------------------------------------
  END MODULE mo_nh_hs_test
