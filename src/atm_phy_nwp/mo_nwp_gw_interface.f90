!>
!! This module is the interface between nwp_nh_interface to the
!! gravity wave drag related parameterisations:
!! inwp_sso == 1 == COSMO subgrid scale orographic gravity wave drag
!! inwp_gwd == 1 == IFS non-orographic gravity wave drag
!!
!! @author    J. SCINOCCIA -- Original Fortran Code
!! @author    A. Orr  -- Rewritten in IFS format   E.C.M.W.F.     August 2008
!!
!! @par Revision History
!! Initial iplementation for ICON Kristina Froehlich, MPI-M, Hamburg (2011-05-19)
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

MODULE mo_nwp_gw_interface

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message

  USE mo_model_domain,         ONLY: t_patch

  USE mo_impl_constants,       ONLY: min_rlcell_int
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c

  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_nonhydro_types,       ONLY: t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_nwp_phy_state,        ONLY: phy_params
  USE mo_parallel_config,      ONLY: nproma
  USE mo_run_config,           ONLY: msg_level
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_sso_cosmo,            ONLY: sso
  USE mo_gwd_wms,              ONLY: gwdrag_wms

  IMPLICIT NONE

  PRIVATE



  PUBLIC  ::  nwp_gwdrag

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE nwp_gwdrag  (   tcall_sso_jg,              & !>input
                         &   lcall_sso_jg,              & !>input
                         &   tcall_gwd_jg,              & !>input
                         &   lcall_gwd_jg,              & !>input
                         &   p_patch,p_metrics,         & !>input
                         &   ext_data,                  & !>input
                         &   p_diag ,                   & !>inout
                         &   prm_diag,prm_nwp_tend      ) !>inout



    TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch         !<grid/patch info.
    TYPE(t_external_data),       INTENT(inout):: ext_data        !< external data, inout only for accomodating ext_data%atm%sso_gamma
    TYPE(t_nh_metrics)          ,INTENT(in)   :: p_metrics
    TYPE(t_nh_diag),      TARGET,INTENT(inout):: p_diag          !<the dyn diag vars
    TYPE(t_nwp_phy_diag),        INTENT(inout):: prm_diag        !<the atm phys vars
    TYPE(t_nwp_phy_tend), TARGET,INTENT(inout):: prm_nwp_tend    !< atm tend vars

    REAL(wp),  INTENT(in)   :: tcall_sso_jg    !< time interval for sso
    LOGICAL ,  INTENT(in)   :: lcall_sso_jg    !< .TRUE.: sso scheme is actually called
    REAL(wp),  INTENT(in)   :: tcall_gwd_jg    !< time interval for gwd
    LOGICAL ,  INTENT(in)   :: lcall_gwd_jg    !< .TRUE.: gwd scheme is actually called

    ! Local array bounds:

    INTEGER :: nlev, nlevp1            !< number of full and half levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !< blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    REAL(wp) ::            &           !< != zonal component of vertical momentum flux (Pa)
      &  z_fluxu(nproma,p_patch%nlevp1 , p_patch%nblks_c)
    REAL(wp) ::            &           !< != meridional component of vertical momentum flux (Pa)
      &  z_fluxv (nproma,p_patch%nlevp1, p_patch%nblks_c)
    REAL(wp) ::            &           !< total precipitation rate [kg/m2/s]
      &  ztot_prec_rate(nproma)

    INTEGER :: jk,jc,jb,jg             !<block indeces


    i_nchdom  = MAX(1,p_patch%n_childdom)

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! domain ID
    jg     = p_patch%id

    ! exclude boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ztot_prec_rate) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)


! Sub-grid Scale Orographic drag

      IF (lcall_sso_jg .AND. atm_phy_nwp_config(jg)%inwp_sso == 1) THEN

        ! ATTENTION: - geopot_agl is the full level geopotential height above ground 
        !              and not full level geopotential height above MSL.
        !            - geopot_agl_ifc is the half level geopotential height above ground 
        !              and not the half level geopotential height above MSL.
        !              I.e. geopot_agl_ifc(:,nlevp1,jb) is not the surface geopotential height,  
        !              but geopot_agl_ifc(:,nlevp1,jb) = 0
        !            Since both geopot_agl AND geopot_agl_ifc are defined above ground, the 
        !            resulting geopotential height which is computed inside of "SSO" is correct. 
        CALL sso(                                          &
          & ie        =nproma                           ,  & !> in:  actual array size
          & ke        =nlev                             ,  & !< in:  actual array size
          & ke1       =nlevp1                           ,  & !< in:  ke + 1
          & istart    =i_startidx                       ,  & !< in:  start index of calculation
          & iend      =i_endidx                         ,  & !< in:  end index of calculation
          & ppf       =p_diag%pres              (:,:,jb),  & !< in:  full level pressure
          & pph       =p_diag%pres_ifc          (:,:,jb),  & !< in:  half level pressure
          & pfif      =p_metrics%geopot_agl     (:,:,jb),  & !< in:  full level geopotential height
          & pt        =p_diag%temp              (:,:,jb),  & !< in:  temperature
          & pu        =p_diag%u                 (:,:,jb),  & !< in:  zonal wind component
          & pv        =p_diag%v                 (:,:,jb),  & !< in:  meridional wind component
          & pfis      =p_metrics%geopot_agl_ifc (:,nlevp1,jb),& !< in:surface geopotential height
          & psso_stdh =ext_data%atm%sso_stdh    (:,jb)  ,  & !< in:  standard deviation
          & psso_gamma=ext_data%atm%sso_gamma   (:,jb)  ,  & !< in:  anisotropy
          & psso_theta=ext_data%atm%sso_theta   (:,jb)  ,  & !< in:  angle
          & psso_sigma=ext_data%atm%sso_sigma   (:,jb)  ,  & !< in:  slope
          & pdt       = tcall_sso_jg                    ,  & !< in:  time step
          & ldebug    =.FALSE.                          ,  & !< in:  debug control switch
          & pdu_sso   =prm_nwp_tend%ddt_u_sso   (:,:,jb),  & !< out: u-tendency due to SSO
          & pdv_sso   =prm_nwp_tend%ddt_v_sso   (:,:,jb)   ) !< out: v-tendency due to SSO
 ! GZ: The computation of the frictional heating rate is now done in interface_nwp for
 ! SSO, GWD and Rayleigh friction together
 !         & pdt_sso   =prm_nwp_tend%ddt_temp_sso(:,:,jb)   ) !< out: temperature tendency
                                                             ! due to SSO

        ! Limit SSO wind tendencies. They can become numerically unstable in the upper stratosphere and mesosphere
        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            prm_nwp_tend%ddt_u_sso(jc,jk,jb) = MAX(-0.1_wp,prm_nwp_tend%ddt_u_sso(jc,jk,jb))
            prm_nwp_tend%ddt_u_sso(jc,jk,jb) = MIN( 0.1_wp,prm_nwp_tend%ddt_u_sso(jc,jk,jb))
            prm_nwp_tend%ddt_v_sso(jc,jk,jb) = MAX(-0.1_wp,prm_nwp_tend%ddt_v_sso(jc,jk,jb))
            prm_nwp_tend%ddt_v_sso(jc,jk,jb) = MIN( 0.1_wp,prm_nwp_tend%ddt_v_sso(jc,jk,jb))
          ENDDO
        ENDDO

      ENDIF

! Non-orgographic gravity wave drag

      IF (lcall_gwd_jg .AND. atm_phy_nwp_config(jg)%inwp_gwd == 1) THEN

        ! get total precipitation rate [kg/m2/s] ==> input for gwdrag_wms
        DO jc =  i_startidx, i_endidx
          ztot_prec_rate(jc) = prm_diag%rain_gsp_rate (jc,jb) &  ! rain_gsp
            &                + prm_diag%snow_gsp_rate (jc,jb) &  ! snow_gsp
            &                + prm_diag%rain_con_rate (jc,jb) &  ! rain_con
            &                + prm_diag%snow_con_rate (jc,jb)    ! snow_con
        ENDDO


        CALL gwdrag_wms(                                   &
           & kidia    = i_startidx                      ,  & 
           & kfdia    = i_endidx                        ,  &
           & klon     = nproma                          ,  &
           & klev     = nlev                            ,  & !< in:  actual array size
           & klaunch  = phy_params(jg)%klaunch          ,  & !< in:  launch level
           & ptstep   = tcall_gwd_jg                    ,  & !< in:  time step
           & ptm1     = p_diag%temp             (:,:,jb),  & !< in:  temperature
           & pum1     = p_diag%u                (:,:,jb),  & !< in:  zonal wind component
           & pvm1     = p_diag%v                (:,:,jb),  & !< in:  meridional wind component
           & papm1    = p_diag%pres             (:,:,jb),  & !< in:  full level pressure
           & paphm1   = p_diag%pres_ifc         (:,:,jb),  & !< in:  half level pressure
           & pgeo1    = p_metrics%geopot_agl    (:,:,jb),  & !< in:  full level geopotential
           & pgelat   = p_patch%cells%center    (:,jb)%lat,& !< in:  latitude (rad)
           & pprecip  = ztot_prec_rate          (:)     ,  & !< in:  total surface precipitation rate
           & ptenu    = prm_nwp_tend%ddt_u_gwd  (:,:,jb),  & !< out: u-tendency
           & ptenv    = prm_nwp_tend%ddt_v_gwd  (:,:,jb),  & !< out: v-tendency
           & pfluxu   = z_fluxu (:,:,jb)                ,  & !< out: zonal  GWD vertical mom flux
           & pfluxv   = z_fluxv (:,:,jb)   )                 !< out: merid. GWD vertical mom flux

        ! Limit also gwdrag wind tendencies. They can become numerically unstable in the upper mesosphere
        DO jk = 1, nlev
!DIR$ IVDEP
          DO jc = i_startidx, i_endidx
            prm_nwp_tend%ddt_u_gwd(jc,jk,jb) = MAX(-0.05_wp,prm_nwp_tend%ddt_u_gwd(jc,jk,jb))
            prm_nwp_tend%ddt_u_gwd(jc,jk,jb) = MIN( 0.05_wp,prm_nwp_tend%ddt_u_gwd(jc,jk,jb))
            prm_nwp_tend%ddt_v_gwd(jc,jk,jb) = MAX(-0.05_wp,prm_nwp_tend%ddt_v_gwd(jc,jk,jb))
            prm_nwp_tend%ddt_v_gwd(jc,jk,jb) = MIN( 0.05_wp,prm_nwp_tend%ddt_v_gwd(jc,jk,jb))
          ENDDO
        ENDDO

      ENDIF

    ENDDO ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE nwp_gwdrag

END MODULE mo_nwp_gw_interface

