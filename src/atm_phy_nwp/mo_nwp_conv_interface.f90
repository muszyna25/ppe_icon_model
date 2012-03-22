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
!!
!! @par Copyright
!! 2002-2011 by DWD and MPI-M
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
MODULE mo_nwp_conv_interface

  USE mo_kind,                 ONLY: wp
  USE mo_parallel_config,      ONLY: nproma
  USE mo_model_domain,         ONLY: t_patch
  USE mo_impl_constants,       ONLY: min_rlcell_int
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag,&
    &                                t_nh_metrics
  USE mo_nonhydrostatic_config,ONLY: kstart_moist
  USE mo_nwp_phy_state,        ONLY: t_nwp_phy_diag, t_nwp_phy_tend, phy_params
  USE mo_run_config,           ONLY: iqv, iqc, iqi !, iqs
  USE mo_physical_constants,   ONLY: grav, alv
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_cumaster,             ONLY: cumastrn
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_icoham_sfc_indices,   ONLY: nsfc_type, iwtr, iice, ilnd

  IMPLICIT NONE

  PRIVATE

  PUBLIC  ::  nwp_convection

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

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

    INTEGER :: nblks_c, nblks_e        !< number of blocks for cells / edges
    INTEGER :: npromz_e, npromz_c      !< length of last block line
    INTEGER :: nlev, nlevp1            !< number of full and half levels
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !< blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    REAL(wp) :: z_omega_p(nproma,p_patch%nlev) !< vertical velocity in p-system
    REAL(wp) :: z_plitot (nproma,p_patch%nlev) !< cloud water + cloud ice
    REAL(wp) :: z_qhfl (nproma,p_patch%nlevp1) !< 3D moisture flux( convection)
    REAL(wp) :: z_shfl (nproma,p_patch%nlevp1) !< 3D sensible heat flux "-"
    REAL(wp) :: z_mflxr(nproma,p_patch%nlevp1) !< 3D conv. rain flux
    REAL(wp) :: z_mflxs(nproma,p_patch%nlevp1) !< 3D conv. snow flux
    REAL(wp) :: z_dtdqv  (nproma,p_patch%nlev) !< 3D moisture convergence
    REAL(wp) :: z_dtdt   (nproma,p_patch%nlev) !< temporal temperature tendency
    REAL(wp) :: z_dtdqv_sv(nproma,p_patch%nlev)!< save array for moisture convergence
    REAL(wp) :: z_dtdt_sv(nproma,p_patch%nlev) !< save array for temperature tendency

    ! Local scalars:

    INTEGER :: jk,jc,jb,jg                !<block indeces


    ! local variables related to the blocking
    nblks_c   = p_patch%nblks_int_c
    npromz_c  = p_patch%npromz_int_c
    nblks_e   = p_patch%nblks_int_e
    npromz_e  = p_patch%npromz_int_e
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

 

!$OMP PARALLEL

#ifdef __xlC__
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_omega_p,z_plitot,z_qhfl,z_shfl,z_mflxr,&
!$OMP            z_mflxs,z_dtdqv,z_dtdt,z_dtdqv_sv,z_dtdt_sv)
#else
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,z_omega_p,z_plitot,z_qhfl,z_shfl,z_mflxr,&
!$OMP            z_mflxs,z_dtdqv,z_dtdt,z_dtdqv_sv,z_dtdt_sv),SCHEDULE(guided)
#endif
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
        !! NOTA: heat fluxes are  defined positive when directed upwards
        !!       so POSITIVE during day over land
        !!      - in Call of convection code sign is reversed
        !! thus pass fluxes to cumastrn that are negative when upwards

        IF(atm_phy_nwp_config(jg)%inwp_turb == 0 ) THEN

          z_qhfl(i_startidx:i_endidx,nlevp1) = - 4.79846_wp*1.e-5_wp !> moisture flux W/m**2
          z_shfl(i_startidx:i_endidx,nlevp1) = - 17._wp              !! sens. heat fl W/m**2

        ELSEIF (atm_phy_nwp_config(jg)%inwp_turb == 1 ) THEN

          prm_diag%qhfl_s(i_startidx:i_endidx,jb) = prm_diag%lhfl_s(i_startidx:i_endidx,jb) & 
             &                                   / alv
          ! PR. In turb1, the flux is positive upwards

          z_qhfl(i_startidx:i_endidx,nlevp1) = - prm_diag%qhfl_s(i_startidx:i_endidx,jb)
          z_shfl(i_startidx:i_endidx,nlevp1) = - prm_diag%shfl_s(i_startidx:i_endidx,jb)

        ELSEIF (atm_phy_nwp_config(jg)%inwp_turb == 2 ) THEN

         ! PR. In turb2, the flux is negative upwards
          z_qhfl(i_startidx:i_endidx,nlevp1) = prm_diag%qhfl_s(i_startidx:i_endidx,jb)

          IF (nsfc_type == 1  ) THEN
            IF (ilnd <= nsfc_type ) THEN   ! sensible heat flux not implemented over land
              z_shfl(i_startidx:i_endidx,nlevp1) = - 17._wp  !! sens. heat fl W/m**2

            ELSE              ! These is the case for which sensible heat
                              ! is implementead in vdiff
              z_shfl(i_startidx:i_endidx,nlevp1) = prm_diag%shfl_s(i_startidx:i_endidx,jb)
            END IF
          ELSE
            z_shfl( i_startidx:i_endidx,nlevp1) = - 17._wp !! sens. heat fl W/m**2 not yet implemented
          ENDIF
        ENDIF

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
        prm_nwp_tend%ddt_temp_pconv  (:,:,jb)     = 0._wp
        prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqc) = 0._wp
        prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqi) = 0._wp
        prm_nwp_tend%ddt_u_pconv     (:,:,jb)     = 0._wp
        prm_nwp_tend%ddt_v_pconv     (:,:,jb)     = 0._wp

        z_mflxs(:,:)   = 0._wp
        z_mflxr(:,:)   = 0._wp

        !-------------------------------------------------------------------------
        !> Convection
        !-------------------------------------------------------------------------

!#ifdef __BOUNDCHECK
          CALL cumastrn &
&           (kidia  = i_startidx ,                kfdia  = i_endidx           ,& !> IN
&            klon   = nproma ,     ktdia  = kstart_moist(jg)  , klev = nlev   ,& !! IN
&            ldland = ext_data%atm%llsm_atm_c(:,jb), ptsphy = tcall_conv_jg   ,& !! IN
&            phy_params = phy_params(jg),                                      & !! IN
&            pten   = p_diag%temp      (:,:,jb)                              ,& !! IN
&            pqen   = p_prog_rcf%tracer(:,:,jb,iqv)                          ,& !! IN
&            puen   = p_diag%u         (:,:,jb), pven   = p_diag%v( :,:,jb) ,& !! IN
&            plitot = z_plitot                 , pvervel= z_omega_p         ,& !! IN
&            pqhfl  = z_qhfl                   , pahfs  = z_shfl            ,& !! IN
&            pap    = p_diag%pres      (:,:,jb), paph   = p_diag%pres_ifc(:,:,jb),& !! IN
&            pgeo   = p_metrics%geopot_agl (:,:,jb)                        ,& !! IN
&            pgeoh  = p_metrics%geopot_agl_ifc(:,:,jb)                        ,& !! IN
&            zdph   = p_diag%dpres_mc    (:,:,jb)                            ,& !! IN
&            zdgeoh = p_metrics%dgeopot_mc(:,:,jb)                            ,& !! IN
&            ptent  = z_dtdt                                                  ,& !! INOUT
&            ptenu  = prm_nwp_tend%ddt_u_pconv     (:,:,jb)                   ,& !! OUT
&            ptenv  = prm_nwp_tend%ddt_v_pconv     (:,:,jb)                   ,& !! OUT
&            ptenq  = z_dtdqv                                                 ,& !! INOUT
&            ptenl  = prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqc)               ,& !! OUT
&            pteni  = prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqi)               ,& !! OUT
!&            ptens  = prm_nwp_tend%ddt_tracer_pconv(:,:,jb,iqs)               ,& !! OUT
&            ldcum  = prm_diag%locum  (:,jb)                                  ,& !! OUT
&            ktype  = prm_diag%ktype   (:,jb)                                 ,& !! OUT
&            kcbot  = prm_diag%mbas_con(:,jb)                                 ,& !! OUT
&            kctop  = prm_diag%mtop_con(:,jb)                                 ,& !! OUT
&            pmfu   =      prm_diag%con_udd(:,:,jb,1)                         ,& !! OUT
&            pmfd   =      prm_diag%con_udd(:,:,jb,2)                         ,& !! OUT
&            pmfude_rate = prm_diag%con_udd(:,:,jb,3)                         ,& !! OUT
&            pmfdde_rate = prm_diag%con_udd(:,:,jb,4)                         ,& !! OUT
&            ptu    =      prm_diag%con_udd(:,:,jb,5)                         ,& !! OUT
&            pqu    =      prm_diag%con_udd(:,:,jb,6)                         ,& !! OUT
&            plu    =      prm_diag%con_udd(:,:,jb,7)                         ,& !! OUT
&            pmflxr =      z_mflxr                                            ,& !! OUT
&            pmflxs =      z_mflxs                                            ,& !! OUT
&            prain  =      prm_diag%rain_upd (:,jb)                           ,& !! OUT
&            pcape =       prm_diag%cape     (:,jb)                           ) !! OUT


        ! Postprocessing on some fields


!             --------------------------------------------
! KF calculate convective gusts according an advice from P. Bechtold
!               from winds in 850 hPa and 925 hPa
!          -  a quite simple aprroximation is used                  
!             --------------------------------------------

!        DO jc = i_startidx,i_endidx
!            IF(ktype(jc) == 1) THEN 
!              DO jk = kctop(jc ),nlev
!          if(( p_diag%pres (jc,jk,jb) =< 85000._wp)) THEN
!            u1= p_diag%u (jc,jk  ,jb)
!            u2= p_diag%u (jc,jk+1,jb)
!          ENDIF
!          if(( p_diag%pres (jc,jk,jb) =< 92500._wp)) THEN
!            u19= p_diag%u (jc,jk  ,jb)
!            u29= p_diag%u (jc,jk+1,jb)
!          ENDIF
!        ENDDO
!          U_850(jc)=0.5_ireals*(u1+u2)
!          U_925(jc)=0.5_ireals*(u19+u29)
!    
!       p_diag%con_gust(j1,jb) =mixfact*MAX(0._ireals,(U_850(jc)-U_925(jc)))
!        ENDIF
!     ENDDO


!          p_diag%extra(:,:,jb,1)= bkaba

          prm_nwp_tend%ddt_temp_pconv  (i_startidx:i_endidx,kstart_moist(jg):,jb) =  &
            &    z_dtdt   (i_startidx:i_endidx,kstart_moist(jg):)                    &
            &  - z_dtdt_sv(i_startidx:i_endidx,kstart_moist(jg):)

          prm_nwp_tend%ddt_tracer_pconv(i_startidx:i_endidx,kstart_moist(jg):,jb,iqv) =  &
            &    z_dtdqv   (i_startidx:i_endidx,kstart_moist(jg):)                       &
            &  - z_dtdqv_sv(i_startidx:i_endidx,kstart_moist(jg):)

          prm_diag%tracer_rate(i_startidx:i_endidx,jb,3) = z_mflxr(i_startidx:i_endidx,nlevp1)
          prm_diag%tracer_rate(i_startidx:i_endidx,jb,4) = z_mflxs(i_startidx:i_endidx,nlevp1)


        ENDIF !inwp_conv

      ENDDO
!$OMP END DO
!$OMP END PARALLEL      

  END SUBROUTINE nwp_convection

END MODULE mo_nwp_conv_interface

