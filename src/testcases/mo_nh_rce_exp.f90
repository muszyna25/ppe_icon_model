!>
!!  Subroutine to initialize a Radiative Convective Equilibrium Exp 
!!
!!
!! @par Revision History
!! - first version by Levi Silvers , MPIM, (2013-4-24)
!! - revised for global RCE model: start from dry and isothermal atmosphere
!! @par Literature
!! -
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
MODULE mo_nh_rce_exp

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_physical_constants,  ONLY: rd, grav, p0ref,rd_o_cpd
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics, t_nh_ref
  USE mo_parallel_config,     ONLY: nproma
  USE mo_nh_testcases_nml,    ONLY: tpe_psfc, tpe_temp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_nh_state_rce_glb

  !DEFINED PARAMETERS (Stevens 2007 JAS):
!! for the initial state, it is unimportant what exact temperature is used
!! this temperature should be consistent with the temperature of the ocean in the RCE case
!! it can be set via the nh_testcase_nml namelist (see the use statement above
!!$  REAL(wp), PARAMETER :: zp0     = 100000._wp !< surface pressure
!!$  REAL(wp), PARAMETER :: zh0     = 0._wp      !< height (m) above which temperature increases
!!$  REAL(wp), PARAMETER :: dtdz    = 0.006_wp   !< lapse rate
!!$  REAL(wp), PARAMETER :: zt0     = 302.15_wp  ! Surface temperature (K)
!  REAL(wp), PARAMETER :: lambda  = 1500._wp   !moist height from Stevens(2007)

!--------------------------------------------------------------------
   CONTAINS
!-------------------------------------------------------------------------
  !>
  !! Initialization of prognostic state vector for the nh RCE test case with a constant T profile
  !!  
  !!
  SUBROUTINE init_nh_state_rce_glb( ptr_patch, ptr_nh_prog,  ptr_nh_ref, ptr_nh_diag,  &
  &                                 ptr_metrics)

    ! INPUT PARAMETERS:
    TYPE(t_patch),TARGET,  INTENT(IN)   :: &  !< patch on which computation is performed
      &  ptr_patch
    TYPE(t_nh_prog),       INTENT(INOUT):: &  !< prognostic state vector
      &  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT):: &  !< diagnostic state vector
      &  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   :: &
      &  ptr_metrics                          !< NH metrics state
    TYPE(t_nh_ref),        INTENT(INOUT):: &  !< reference state vector
      &  ptr_nh_ref

    INTEGER  :: jb,jk  ! loop indices
    INTEGER  :: nblks_c,npromz_c
    INTEGER  :: nlen,nlev

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c

    ! number of vertical levels
    nlev   = ptr_patch%nlev
  
    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = tpe_psfc
  
    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp

    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = 1, nlev
        ! start from an isothermal atmosphere with temperature tpe_temp
        ! set pressure and Exner function first
        ptr_nh_diag%pres(1:nlen,jk,jb)    = tpe_psfc * EXP(-grav/rd/tpe_temp*ptr_metrics%z_mc(1:nlen,jk,jb))
        ptr_nh_prog%rho(1:nlen,jk,jb)     = ptr_nh_diag%pres(1:nlen,jk,jb)/rd/tpe_temp
        ptr_nh_prog%exner(1:nlen,jk,jb)   = (ptr_nh_diag%pres(1:nlen,jk,jb)/p0ref)**rd_o_cpd
        ptr_nh_prog%theta_v(1:nlen,jk,jb) = tpe_temp/ptr_nh_prog%exner(1:nlen,jk,jb)
!!$!++jsr
!!$        write(0,*) 'ptr_nh_prog%exner(1,jk,jb)=',ptr_nh_prog%exner(1,jk,jb), &
!!$                   'ptr_nh_prog%rho(1,jk,jb)=',ptr_nh_prog%rho(1,jk,jb), 'tpe_psfc=',tpe_psfc      
!!$!--jsr
     END DO


!!$      DO jk = 1 , nlev
!!$        ptr_nh_prog%rho(1:nlen,jk,jb) = (ptr_nh_prog%exner(1:nlen,jk,jb)**cvd_o_rd)*p0ref/rd / &
!!$                                         ptr_nh_prog%theta_v(1:nlen,jk,jb)     
!!$!++jsr
!!$
!!$!--jsr
!!$      END DO 

  
    END DO ! jb

  !meridional and zonal wind
  ptr_nh_prog%vn = 0._wp
  ptr_nh_ref%vn_ref = ptr_nh_prog%vn

  !vertical wind
  ptr_nh_prog%w = 0._wp
  ptr_nh_ref%w_ref = ptr_nh_prog%w

  END SUBROUTINE init_nh_state_rce_glb

END MODULE mo_nh_rce_exp

