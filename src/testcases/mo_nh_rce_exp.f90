!>
!!  Subroutine to initialize a Radiative Convective Equilibrium Exp 
!!
!!
!! @par Revision History
!! - first version by Levi Silvers , MPIM, (2013-4-24)
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
MODULE mo_nh_rce_exp

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_physical_constants,  ONLY: rd, cpd, grav, p0ref,rd_o_cpd, cvd_o_rd
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics, t_nh_ref
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: iqv, iqc

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_nh_state_rce_glb

  !DEFINED PARAMETERS (Stevens 2007 JAS):
  REAL(wp), PARAMETER :: zp0     = 100000._wp !< surface pressure
  REAL(wp), PARAMETER :: zh0     = 0._wp      !< height (m) above which temperature increases
  REAL(wp), PARAMETER :: dtdz    = 0.006_wp   !< lapse rate
  REAL(wp), PARAMETER :: zt0     = 302.15_wp  ! Surface temperature (K)
!  REAL(wp), PARAMETER :: lambda  = 1500._wp   !moist height from Stevens(2007)

!--------------------------------------------------------------------
   CONTAINS
!-------------------------------------------------------------------------
  !>
  !! Initialization of prognostic state vector for the nh RCE test case 
  !!
  SUBROUTINE init_nh_state_rce_glb( ptr_patch, ptr_nh_prog,  ptr_nh_ref, ptr_nh_diag,  &
  &                           ptr_int, ptr_ext_data, ptr_metrics)

    ! INPUT PARAMETERS:
    TYPE(t_patch),TARGET,  INTENT(IN)   :: &  !< patch on which computation is performed
      &  ptr_patch
    TYPE(t_int_state),     INTENT(IN)   :: &
      &  ptr_int
    TYPE(t_nh_prog),       INTENT(INOUT):: &  !< prognostic state vector
      &  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT):: &  !< diagnostic state vector
      &  ptr_nh_diag
    TYPE(t_external_data), INTENT(INOUT):: &  !< external data
      &  ptr_ext_data
    TYPE(t_nh_metrics),    INTENT(IN)   :: &
      &  ptr_metrics                          !< NH metrics state
    TYPE(t_nh_ref),        INTENT(INOUT):: &  !< reference state vector
      &  ptr_nh_ref

    REAL(wp) :: rho_sfc, z_help(1:nproma)!, zvn1, zvn2, zu, zv
    INTEGER  :: jb,jk  ! loop indices
    INTEGER  :: nblks_c,npromz_c,nblks_e,npromz_e
    INTEGER  :: nlen,nlev,jg

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c
    nblks_e  = ptr_patch%nblks_e
    npromz_e = ptr_patch%npromz_e

    ! number of vertical levels
    nlev   = ptr_patch%nlev

    !patch id
    jg = ptr_patch%id

    ! use equation of state to set a reference density
    rho_sfc = zp0 / (rd * zt0 )
  
    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = zp0
  
    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp

    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF
! I think that as of now this is a dry configuration and that to introduce
! moisture I will need to initialize tracer(iqv) as something like that below

  !  ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = 0.8_wp * spec_humi(sat_pres_water(zt0),zp0) * &
  !             EXP(-ptr_metrics%z_mc(1:nlen,jk,jb)/lambda)

      DO jk = 1, nlev
        ! init potential temperature 
        ! z_mc is geometric height at full level center
        ! why isn't a reference pressure used here?  this looks like temp
        z_help(1:nlen) = zt0 + max(0._wp, (ptr_metrics%z_mc(1:nlen,jk,jb)-zh0)*dtdz)
    
        ! virtual potential temperature
        ptr_nh_prog%theta_v(1:nlen,jk,jb) = z_help(1:nlen) * ( 1._wp + &
            0.61_wp*ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) - ptr_nh_prog%tracer(1:nlen,jk,jb,iqc) ) 
      END DO 
 
      !Get hydrostatic pressure and exner at lowest level
      ptr_nh_diag%pres(1:nlen,nlev,jb) = zp0 - rho_sfc * ptr_metrics%geopot(1:nlen,nlev,jb)
      ptr_nh_prog%exner(1:nlen,nlev,jb) = (ptr_nh_diag%pres(1:nlen,nlev,jb)/p0ref)**rd_o_cpd 

      !Get exner at other levels
      DO jk = nlev-1, 1, -1
      ! average of virtual pot. temp.
         z_help(1:nlen) = 0.5_wp * ( ptr_nh_prog%theta_v(1:nlen,jk,jb) +  &
                                     ptr_nh_prog%theta_v(1:nlen,jk+1,jb) )
      ! exner : what is this i.c.? Look up the prog equation for Exner   
         ptr_nh_prog%exner(1:nlen,jk,jb) = ptr_nh_prog%exner(1:nlen,jk+1,jb) &
            &  -grav/cpd*ptr_metrics%ddqz_z_half(1:nlen,jk+1,jb)/z_help(1:nlen)
      END DO

      DO jk = 1 , nlev
         ! rhotheta has to have the same meaning as exner
         ptr_nh_prog%rhotheta_v(1:nlen,jk,jb) = &
              (ptr_nh_prog%exner(1:nlen,jk,jb)**cvd_o_rd)*p0ref/rd

         ptr_nh_prog%rho(1:nlen,jk,jb) = ptr_nh_prog%rhotheta_v(1:nlen,jk,jb) / &
                                         ptr_nh_prog%theta_v(1:nlen,jk,jb)     
      END DO 
  
    END DO ! jb

  !meridional and zonal wind
  ptr_nh_prog%vn = 0._wp
  ptr_nh_ref%vn_ref = ptr_nh_prog%vn

  !vertical wind
  ptr_nh_prog%w = 0._wp
  ptr_nh_ref%w_ref = ptr_nh_prog%w

  END SUBROUTINE init_nh_state_rce_glb

END MODULE mo_nh_rce_exp

