!OPTION! -cont
!! this command should fix the problem of copying arrays in a subroutine call
!>
!! This module is the interface between nwp_nh_interface to the 
!! convection parameterisation(s):
!! inwp_conv == 1 == Tiedtke-Bechtold convection
!!
!! @author Kristina Froehlich, DWD, Offenbach (2010-01-25)
!!
!! @par Revision History
!! Initial Kristina Froehlich, DWD, Offenbach (2010-01-25)
!! Add gust   by Helmut Frank, DWD, Offenbach (2013-03-13)
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

MODULE mo_nwp_conv_interface

  USE mo_kind,                 ONLY: wp
  USE mo_parallel_config,      ONLY: nproma
  USE mo_model_domain,         ONLY: t_patch
  USE mo_impl_constants,       ONLY: min_rlcell_int, icosmo, igme, iedmf
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag,&
    &                                t_nh_metrics
  USE mo_nonhydrostatic_config,ONLY: kstart_moist
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_nwp_phy_state,        ONLY: phy_params
  USE mo_run_config,           ONLY: iqv, iqc, iqi, iqr, iqs, nqtendphy
  USE mo_physical_constants,   ONLY: grav, alf, cvd
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_cumaster,             ONLY: cumastrn
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_art_config,           ONLY: art_config
  USE mo_util_phys,            ONLY: nwp_con_gust
!!$  USE mo_cuparameters,         ONLY: lmfscv

  IMPLICIT NONE

  PRIVATE



  PUBLIC  ::  nwp_convection

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE nwp_convection ( tcall_conv_jg,             & !>input
    &                         p_patch,p_metrics,         & !>input
    &                         ext_data,                  & !>input
    &                         p_prog,                    & !>in
    &                         p_prog_rcf,                & !>inout
    &                         p_diag ,                   & !>inout
    &                         prm_diag,prm_nwp_tend      ) !>inout



    TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !!<grid/patch info.
    TYPE(t_external_data),       INTENT(in)   :: ext_data        !< external data
    TYPE(t_nh_metrics)          ,INTENT(in)   :: p_metrics
    TYPE(t_nh_prog),      TARGET,INTENT(in)   :: p_prog          !<the dyn prog vars
    TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog_rcf      !<call freq
    TYPE(t_nh_diag),      TARGET,INTENT(inout):: p_diag          !<the dyn diag vars
    TYPE(t_nwp_phy_diag),        INTENT(inout):: prm_diag        !<the atm phys vars
    TYPE(t_nwp_phy_tend),TARGET, INTENT(inout):: prm_nwp_tend    !< atm tend vars

    REAL(wp),                    INTENT(in)   :: tcall_conv_jg   !< time interval for 
                                                                 !< convection
    ! Local array bounds:

    INTEGER :: nlev, nlevp1            !< number of full and half levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !< blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    REAL(wp) :: z_omega_p(nproma,p_patch%nlev) !< vertical velocity in p-system
    REAL(wp) :: z_plitot (nproma,p_patch%nlev) !< cloud water + cloud ice
    REAL(wp) :: z_qhfl (nproma,p_patch%nlevp1) !< 3D moisture flux( convection)
    REAL(wp) :: z_shfl (nproma,p_patch%nlevp1) !< 3D sensible heat flux "-"
    REAL(wp) :: z_dtdqv  (nproma,p_patch%nlev) !< 3D moisture convergence
    REAL(wp) :: z_dtdt   (nproma,p_patch%nlev) !< temporal temperature tendency
    REAL(wp) :: z_dtdqv_sv(nproma,p_patch%nlev)!< save array for moisture convergence
    REAL(wp) :: z_dtdt_sv(nproma,p_patch%nlev) !< save array for temperature tendency

    ! Local scalars:

    INTEGER  :: jk,jc,jb,jg                !< block indeces
    INTEGER  :: zk850, zk950               !< level indices
    REAL(wp) :: u850, u950, v850, v950     !< zonal and meridional velocity at specific heights
    REAL(wp) :: ticeini, lfocvd, wfac
    INTEGER  :: iqrd, iqsd

    ! local variables related to the blocking
    i_nchdom  = MAX(1,p_patch%n_childdom)
    jg        = p_patch%id

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    !in order to account for mesh refinement
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    lfocvd  = alf/cvd
    ticeini = 258.15_wp

    ! IDs for optional arguments for detrainment of rain and snow
    IF (atm_phy_nwp_config(jg)%ldetrain_conv_prec) THEN
      iqrd = iqr
      iqsd = iqs
    ELSE
      iqrd = nqtendphy
      iqsd = nqtendphy
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_omega_p,z_plitot,z_qhfl,z_shfl,z_dtdqv,&
!$OMP            z_dtdt,z_dtdqv_sv,z_dtdt_sv,zk850,zk950,u850,u950,v850,v950,wfac), ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)


      !-------------------------------------------------------------------------
      !> Calculate vertical velocity in p-system
      !-------------------------------------------------------------------------
      !
      DO jk = kstart_moist(jg),nlev
        DO jc = i_startidx,i_endidx
          z_omega_p(jc,jk)= -0.5_wp*(p_prog%w(jc,jk,jb)+p_prog%w(jc,jk+1,jb)) &
                            *p_prog%rho(jc,jk,jb)*grav
        ENDDO
      ENDDO

      IF( atm_phy_nwp_config(jg)%inwp_convection == 1 ) THEN

        !>
        !! define convective-related fields
        !! NOTE: Heat fluxes are defined negative upwards in the convection scheme. 
        !!       In the turbulences scheme (1,2,3) they defined either positive or
        !!       negative upward. 
        !! Thus pass fluxes to cumastrn that are negative when upwards!!!

        SELECT CASE (atm_phy_nwp_config(jg)%inwp_turb)
        CASE (0)

          z_qhfl(i_startidx:i_endidx,nlevp1) = - 4.79846_wp*1.e-5_wp !> moisture flux kg/m2/s
          z_shfl(i_startidx:i_endidx,nlevp1) = - 17._wp              !! sens. heat fl W/m**2

        CASE (icosmo,igme,iedmf,10,11,12)

          ! In turb1,turb2 and turb3, the flux is positive downwards / negative upwards

          z_qhfl(i_startidx:i_endidx,nlevp1) = prm_diag%qhfl_s(i_startidx:i_endidx,jb)
          z_shfl(i_startidx:i_endidx,nlevp1) = prm_diag%shfl_s(i_startidx:i_endidx,jb)


        END SELECT

        z_omega_p (:,1:kstart_moist(jg)-1) = 0._wp
        z_dtdqv   (:,1:kstart_moist(jg)-1) = 0._wp
        z_dtdt    (:,1:kstart_moist(jg)-1) = 0._wp
        z_plitot  (:,1:kstart_moist(jg)-1) = 0._wp


        DO jk = kstart_moist(jg),nlev
          DO jc = i_startidx,i_endidx
            ! moisture convergence
            z_dtdqv(jc,jk) =                                                                 &
              p_diag%ddt_tracer_adv(jc,jk,jb,iqv) + prm_nwp_tend%ddt_tracer_turb(jc,jk,jb,iqv)
            z_dtdqv_sv(jc,jk) = z_dtdqv(jc,jk)

            ! temperature tendencies from other physical processes (used for mass flux closure)
            z_dtdt(jc,jk) =                              &
              &    prm_nwp_tend%ddt_temp_radsw(jc,jk,jb) &
              &  + prm_nwp_tend%ddt_temp_radlw(jc,jk,jb) &
              &  + prm_nwp_tend%ddt_temp_turb (jc,jk,jb) &
              &  + p_diag%ddt_temp_dyn(jc,jk,jb)
            z_dtdt_sv(jc,jk) = z_dtdt(jc,jk)

            ! cloud water + cloud ice for entrainment computation
            z_plitot(jc,jk) = p_prog_rcf%tracer(jc,jk,jb,iqc) &
                            + p_prog_rcf%tracer(jc,jk,jb,iqi)
          ENDDO
        ENDDO

        ! The following input fields must be reset to zero because the convective
        ! tendencies are added to them
        prm_nwp_tend%ddt_tracer_pconv(:,:,jb,2:)  = 0._wp
        prm_nwp_tend%ddt_u_pconv     (:,:,jb)     = 0._wp
        prm_nwp_tend%ddt_v_pconv     (:,:,jb)     = 0._wp

        prm_diag%rain_con_rate_3d(:,:,jb)         = 0._wp
        prm_diag%snow_con_rate_3d(:,:,jb)         = 0._wp

        !-------------------------------------------------------------------------
        !> Convection
        !-------------------------------------------------------------------------

        IF ( atm_phy_nwp_config(jg)%inwp_turb /= iedmf ) THEN ! DUALM is allow to turn of shallow convection
          DO jc = i_startidx,i_endidx
            prm_diag%ldshcv(jc,jb) = .TRUE.
          ENDDO
        ENDIF


        IF(art_config(jg)%nconv_tracer > 0) THEN
          CALL cumastrn &
&           (kidia  = i_startidx            , kfdia  = i_endidx               ,& !> IN
&            klon   = nproma ,     ktdia  = kstart_moist(jg)  , klev = nlev   ,& !! IN
&            ldland = ext_data%atm%llsm_atm_c(:,jb), ptsphy = tcall_conv_jg   ,& !! IN
&            ldlake = ext_data%atm%llake_c(:,jb), k950 = prm_diag%k950(:,jb)  ,& !! IN
&            phy_params = phy_params(jg), capdcfac=prm_diag%tropics_mask(:,jb),& !! IN
&            mtnmask=p_metrics%mask_mtnpoints(:,jb), pten=p_diag%temp(:,:,jb) ,& !! IN
&            pqen   = p_prog_rcf%tracer(:,:,jb,iqv)                           ,& !! IN
&            puen   = p_diag%u   (:,:,jb)   , pven   = p_diag%v( :,:,jb)      ,& !! IN
&            plitot = z_plitot              , pvervel= z_omega_p              ,& !! IN
&            pqhfl  = z_qhfl                , pahfs  = z_shfl                 ,& !! IN
&            pap    = p_diag%pres(:,:,jb)   , paph   = p_diag%pres_ifc(:,:,jb),& !! IN
&            pgeo   = p_metrics%geopot_agl (:,:,jb)                           ,& !! IN
&            pgeoh  = p_metrics%geopot_agl_ifc(:,:,jb)                        ,& !! IN
&            zdph   = p_diag%dpres_mc     (:,:,jb)                            ,& !! IN
&            zdgeoh = p_metrics%dgeopot_mc(:,:,jb)                            ,& !! IN
&            pcloudnum = prm_diag%cloud_num(:,jb)                             ,& !! IN
&            ptent  = z_dtdt                                                  ,& !! INOUT
&            ptenu  = prm_nwp_tend%ddt_u_pconv     (:,:,jb)                   ,& !! OUT
&            ptenv  = prm_nwp_tend%ddt_v_pconv     (:,:,jb)                   ,& !! OUT
&            ptenq  = z_dtdqv                                                 ,& !! INOUT
&            ptenl  = prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqc)               ,& !! OUT
&            pteni  = prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqi)               ,& !! OUT
&            ptenr  = prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqrd)              ,& !! OUT
&            ptens  = prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqsd)              ,& !! OUT
&            ldcum  = prm_diag%locum   (:,jb)                                 ,& !! OUT
&            ktype  = prm_diag%ktype   (:,jb)                                 ,& !! OUT
&            kcbot  = prm_diag%mbas_con(:,jb)                                 ,& !! OUT
&            kctop  = prm_diag%mtop_con(:,jb)                                 ,& !! OUT
&            LDSHCV = prm_diag%ldshcv  (:,jb)                                 ,& !! IN
&            pmfu   =      prm_diag%con_udd(:,:,jb,1)                         ,& !! OUT
&            pmfd   =      prm_diag%con_udd(:,:,jb,2)                         ,& !! OUT
&            pmfude_rate = prm_diag%con_udd(:,:,jb,3)                         ,& !! OUT
&            pmfdde_rate = prm_diag%con_udd(:,:,jb,4)                         ,& !! OUT
&            ptu    =      prm_diag%con_udd(:,:,jb,5)                         ,& !! OUT
&            pqu    =      prm_diag%con_udd(:,:,jb,6)                         ,& !! OUT
&            plu    =      prm_diag%con_udd(:,:,jb,7)                         ,& !! OUT
&            pmflxr =      prm_diag%rain_con_rate_3d(:,:,jb)                  ,& !! OUT
&            pmflxs =      prm_diag%snow_con_rate_3d(:,:,jb)                  ,& !! OUT
&            prain  =      prm_diag%rain_upd (:,jb)                           ,& !! OUT
&            pdtke_con =   prm_nwp_tend%ddt_tke_pconv(:,:,jb)                 ,& !! OUT
&            pcape  =      prm_diag%cape     (:,jb)                           ,& !! OUT
&            ktrac  = art_config(jg)%nconv_tracer                             ,& !! IN 
&            pcen   = p_prog_rcf%conv_tracer(jb,:)                            ,& !! IN 
&            ptenc  = prm_nwp_tend%conv_tracer_tend(jb,:) )                      !! OUT

        ELSE

          CALL cumastrn &
&           (kidia  = i_startidx            , kfdia  = i_endidx               ,& !> IN
&            klon   = nproma ,     ktdia  = kstart_moist(jg)  , klev = nlev   ,& !! IN
&            ldland = ext_data%atm%llsm_atm_c(:,jb), ptsphy = tcall_conv_jg   ,& !! IN
&            ldlake = ext_data%atm%llake_c(:,jb), k950 = prm_diag%k950(:,jb)  ,& !! IN
&            phy_params = phy_params(jg), capdcfac=prm_diag%tropics_mask(:,jb),& !! IN
&            mtnmask=p_metrics%mask_mtnpoints(:,jb), pten=p_diag%temp(:,:,jb) ,& !! IN
&            pqen   = p_prog_rcf%tracer(:,:,jb,iqv)                           ,& !! IN
&            puen   = p_diag%u   (:,:,jb)   , pven   = p_diag%v( :,:,jb)      ,& !! IN
&            plitot = z_plitot              , pvervel= z_omega_p              ,& !! IN
&            pqhfl  = z_qhfl                , pahfs  = z_shfl                 ,& !! IN
&            pap    = p_diag%pres(:,:,jb)   , paph   = p_diag%pres_ifc(:,:,jb),& !! IN
&            pgeo   = p_metrics%geopot_agl (:,:,jb)                           ,& !! IN
&            pgeoh  = p_metrics%geopot_agl_ifc(:,:,jb)                        ,& !! IN
&            zdph   = p_diag%dpres_mc     (:,:,jb)                            ,& !! IN
&            zdgeoh = p_metrics%dgeopot_mc(:,:,jb)                            ,& !! IN
&            pcloudnum = prm_diag%cloud_num(:,jb)                             ,& !! IN
&            ptent  = z_dtdt                                                  ,& !! INOUT
&            ptenu  = prm_nwp_tend%ddt_u_pconv     (:,:,jb)                   ,& !! OUT
&            ptenv  = prm_nwp_tend%ddt_v_pconv     (:,:,jb)                   ,& !! OUT
&            ptenq  = z_dtdqv                                                 ,& !! INOUT
&            ptenl  = prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqc)               ,& !! OUT
&            pteni  = prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqi)               ,& !! OUT
&            ptenr  = prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqrd)              ,& !! OUT
&            ptens  = prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqsd)              ,& !! OUT
&            ldcum  = prm_diag%locum   (:,jb)                                 ,& !! OUT
&            ktype  = prm_diag%ktype   (:,jb)                                 ,& !! OUT
&            kcbot  = prm_diag%mbas_con(:,jb)                                 ,& !! OUT
&            kctop  = prm_diag%mtop_con(:,jb)                                 ,& !! OUT
&            LDSHCV = prm_diag%ldshcv  (:,jb)                                 ,& !! IN
&            pmfu   =      prm_diag%con_udd(:,:,jb,1)                         ,& !! OUT
&            pmfd   =      prm_diag%con_udd(:,:,jb,2)                         ,& !! OUT
&            pmfude_rate = prm_diag%con_udd(:,:,jb,3)                         ,& !! OUT
&            pmfdde_rate = prm_diag%con_udd(:,:,jb,4)                         ,& !! OUT
&            ptu    =      prm_diag%con_udd(:,:,jb,5)                         ,& !! OUT
&            pqu    =      prm_diag%con_udd(:,:,jb,6)                         ,& !! OUT
&            plu    =      prm_diag%con_udd(:,:,jb,7)                         ,& !! OUT
&            pmflxr =      prm_diag%rain_con_rate_3d(:,:,jb)                  ,& !! OUT
&            pmflxs =      prm_diag%snow_con_rate_3d(:,:,jb)                  ,& !! OUT
&            prain  =      prm_diag%rain_upd (:,jb)                           ,& !! OUT
&            pdtke_con =   prm_nwp_tend%ddt_tke_pconv(:,:,jb)                 ,& !! OUT
&            pcape  =      prm_diag%cape     (:,jb)                           ,& !! OUT
&            ktrac  = 0                                                        ) !! IN 
        ENDIF


        ! Postprocessing on some fields


        prm_nwp_tend%ddt_temp_pconv  (i_startidx:i_endidx,kstart_moist(jg):,jb) =  &
          &    z_dtdt   (i_startidx:i_endidx,kstart_moist(jg):)                    &
          &  - z_dtdt_sv(i_startidx:i_endidx,kstart_moist(jg):)

        prm_nwp_tend%ddt_tracer_pconv(i_startidx:i_endidx,kstart_moist(jg):,jb,iqv) =  &
          &    z_dtdqv   (i_startidx:i_endidx,kstart_moist(jg):)                       &
          &  - z_dtdqv_sv(i_startidx:i_endidx,kstart_moist(jg):)

        ! Convert detrained cloud ice into cloud water if the temperature is only slightly below freezing
        ! and convective cloud top is not cold enough for substantial ice initiation
        DO jk = kstart_moist(jg),nlev
          DO jc = i_startidx,i_endidx
            IF (prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi) > 0._wp .AND. p_diag%temp(jc,jk,jb) > ticeini) THEN
              wfac = MAX(0._wp, MIN(1._wp,0.25_wp*(p_diag%temp(jc,jk,jb)-ticeini)) + &
                     0.25_wp*MIN(0._wp,p_diag%temp(jc,prm_diag%mtop_con(jc,jb),jb)-ticeini) )
              prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqc) = prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqc) + &
                wfac*prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi)
              prm_nwp_tend%ddt_temp_pconv(jc,jk,jb) = prm_nwp_tend%ddt_temp_pconv(jc,jk,jb) - &
                lfocvd*wfac*prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi)
              prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi) = (1._wp-wfac)*prm_nwp_tend%ddt_tracer_pconv(jc,jk,jb,iqi)
            ENDIF
          ENDDO
        ENDDO

        ! convective contribution to wind gust 
        ! (based on simple parameterization by Peter Bechthold)
        !
        DO jc=i_startidx,i_endidx
          IF ( prm_diag%ktype(jc,jb) == 1 )  THEN   ! penetrative convection
            zk850 = prm_diag%k850(jc,jb)
            zk950 = prm_diag%k950(jc,jb)
            ! We take the arithmetic mean of u(jc,zk850,jb) and u(jc,zk850-1,jb)
            ! as well as v(jc,zk850,jb) and v(jc,zk850-1,jb), since the levels 
            ! zk850 and zk950 are located just below the respective threshold heights.
            u850 = 0.5_wp * (p_diag%u(jc,zk850,jb) + p_diag%u(jc,zk850-1,jb))
            u950 = 0.5_wp * (p_diag%u(jc,zk950,jb) + p_diag%u(jc,zk950-1,jb))
            v850 = 0.5_wp * (p_diag%v(jc,zk850,jb) + p_diag%v(jc,zk850-1,jb))
            v950 = 0.5_wp * (p_diag%v(jc,zk950,jb) + p_diag%v(jc,zk950-1,jb))

            prm_diag%con_gust(jc,jb) = nwp_con_gust( u850, u950, v850, v950 )
          ELSE
            prm_diag%con_gust(jc,jb) = 0._wp
          ENDIF
        ENDDO  ! jc


      ENDIF !inwp_conv

    ENDDO  ! jb
!$OMP END DO
!$OMP END PARALLEL      

  END SUBROUTINE nwp_convection

END MODULE mo_nwp_conv_interface

