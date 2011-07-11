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
  USE mo_exception,            ONLY: message, message_text, finish
  USE mo_mpi,                  ONLY: p_pe, p_nprocs
  USE mo_parallel_configuration, ONLY:  p_test_run, nproma

  USE mo_model_domain,         ONLY: t_patch
  USE mo_grf_interpolation,    ONLY: t_gridref_state
  USE mo_impl_constants,       ONLY: min_rlcell_int
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c

!  USE mo_ext_data,             ONLY: t_external_data
  USE mo_nonhydro_state,       ONLY: t_nh_prog, t_nh_diag,&
    &                                t_nh_metrics
  USE mo_nonhydrostatic_nml,   ONLY: kstart_moist
  USE mo_nwp_phy_state,        ONLY: t_nwp_phy_diag,prm_diag,&
    &                                t_nwp_phy_tend
  USE mo_run_nml,              ONLY: msg_level, ntracer, iqv, &
    &                                iqc, iqi, iqs
  USE mo_physical_constants,   ONLY:  vtmpc1, grav, alv
!  USE mo_atm_phy_nwp_nml,      ONLY: inwp_convection, inwp_turb
  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config
  USE mo_cumaster,             ONLY: cumastrn

  IMPLICIT NONE

  PRIVATE

  PUBLIC  ::  nwp_convection

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE nwp_convection  ( tcall_conv_jg,              & !>input
                            &   p_patch,p_metrics,         & !>input
!                            &   ext_data,                  & !>input
                            &   p_prog,                    & !>in
                            &   p_prog_rcf,                & !>inout
                            &   p_diag ,                   & !>inout
                            &   prm_diag,prm_nwp_tend      ) !>inout 



    TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !!<grid/patch info.
 !   TYPE(t_external_data),       INTENT(in)   :: ext_data        !< external data
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

    REAL(wp) ::            &           !< vertical velocity in p-system
      &  z_omega_p(nproma,p_patch%nlev  ,p_patch%nblks_c) 
    REAL(wp) ::            &           !< cloud water + cloud ice (set to 0)
      &  z_plitot (nproma,p_patch%nlev  ,p_patch%nblks_c) 
    REAL(wp) ::            &           !< 3D moisture flux( convection)
      &  z_qhfl   (nproma,p_patch%nlevp1,p_patch%nblks_c) 
    REAL(wp) ::            &           !< 3D sensible heat flux "-"
      &  z_shfl   (nproma,p_patch%nlevp1,p_patch%nblks_c) 
    REAL(wp) ::            &           !< 3D conv. rain flux
      &  z_mflxr  (nproma,p_patch%nlevp1,p_patch%nblks_c) 
    REAL(wp) ::            &           !< 3D conv. snow flux
      &  z_mflxs  (nproma,p_patch%nlevp1,p_patch%nblks_c) 
    REAL(wp) ::            &           !< 3D moisture convergence
      &  z_dtdqv  (nproma,p_patch%nlev  ,p_patch%nblks_c) 
    REAL(wp) ::            &           !< temporal temperature tendency
      &  z_dtdt  (nproma,p_patch%nlev  ,p_patch%nblks_c) 

    ! Local scalars:

    INTEGER :: jk,jc,jb,jg                !<block indeces

    !KF temporary field
    LOGICAL:: landseemask(nproma,p_patch%nblks_c)


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
!$OMP WORKSHARE
    landseemask(:,:)   = .FALSE.
    z_plitot (:,:,:)   = 0._wp
    z_mflxs  (:,:,:)   = 0._wp
    z_mflxr  (:,:,:)   = 0._wp
!$OMP END WORKSHARE

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx),SCHEDULE(guided)

      DO jb = i_startblk, i_endblk
        !
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
&                       i_startidx, i_endidx, rl_start, rl_end)


        !-------------------------------------------------------------------------
        !> Calculate vertical velocity in p-system
        !-------------------------------------------------------------------------
        !
        DO jk = 1,nlev
          DO jc = i_startidx,i_endidx
            z_omega_p(jc,jk,jb)= -p_prog%w(jc,jk,jb)*p_prog%rho(jc,jk,jb)*grav
          ENDDO
        ENDDO

        IF( atm_phy_nwp_config(jg)%inwp_convection == 1 ) THEN

        !>
        !! define convective-related fields
        !!KF testvalues taken from Peter Bechtold
        !! NOTA: heat fluxes are  defined positive when directed upwards
        !!       so POSITIVE during day over land
        !!      - in Call of convection code sign is reversed

          !KF preliminary setup for fluxes



!PR pass fluxes to cumastrn that are negative when upwards!!!!


        IF(atm_phy_nwp_config(jg)%inwp_turb == 0 ) THEN

        z_qhfl( i_startidx:i_endidx,nlevp1,jb) = - 4.79846_wp*1.e-5_wp !> moisture flux W/m**2
        z_shfl( i_startidx:i_endidx,nlevp1,jb) = -   17._wp              !! sens. heat fl W/m**2

        ELSEIF (atm_phy_nwp_config(jg)%inwp_turb == 1 ) THEN

        prm_diag%qhfl_s( i_startidx:i_endidx,jb) = prm_diag%lhfl_s( i_startidx:i_endidx,jb) & 
             &                                   / alv
        ! PR. In turb1, the flux is positive upwards

        z_qhfl( i_startidx:i_endidx,nlevp1,jb) = - prm_diag%qhfl_s( i_startidx:i_endidx,jb)
        z_shfl( i_startidx:i_endidx,nlevp1,jb) = - prm_diag%shfl_s( i_startidx:i_endidx,jb)
!
!KF      Fluxes are not yet ready for use

!        z_qhfl( i_startidx:i_endidx,nlevp1,jb) =  4.79846_wp*1.e-5_wp !> moisture flux W/m**2
!        z_shfl( i_startidx:i_endidx,nlevp1,jb) =    17._wp            !! sens. heat fl W/m**2

        ELSEIF (atm_phy_nwp_config(jg)%inwp_turb == 2 ) THEN

        ! PR. In turb2, the flux is negative upwards
        z_qhfl( i_startidx:i_endidx,nlevp1,jb) = prm_diag%qhfl_s (i_startidx:i_endidx,jb)
        z_shfl( i_startidx:i_endidx,nlevp1,jb) = -   17._wp            !! sens. heat fl W/m**2

        ENDIF

        z_dtdqv(i_startidx:i_endidx,:,jb) = 0.0_wp                                &
                           &   +p_diag%ddt_tracer_adv(i_startidx:i_endidx,:,jb,iqv)

        ! input from other physical processes on the convection
        z_dtdt(i_startidx:i_endidx,:,jb)= 0._wp                                                  &
          &                              + prm_nwp_tend%ddt_temp_radsw(i_startidx:i_endidx,:,jb) &
          &                              + prm_nwp_tend%ddt_temp_radlw(i_startidx:i_endidx,:,jb)
!&                                       + prm_nwp_tend%ddt_temp_turb(i_startidx:i_endidx,:,jb)

        !KF its a must to set them to zero!
        prm_nwp_tend%ddt_temp_pconv  (i_startidx:i_endidx,:,jb)   = 0._wp
        prm_nwp_tend%ddt_tracer_pconv(i_startidx:i_endidx,:,jb,:) = 0._wp
        prm_nwp_tend%ddt_u_pconv     (i_startidx:i_endidx,:,jb)   = 0._wp
        prm_nwp_tend%ddt_v_pconv     (i_startidx:i_endidx,:,jb)   = 0._wp


        !-------------------------------------------------------------------------
        !> Convection
        !-------------------------------------------------------------------------

!#ifdef __BOUNDCHECK
          CALL cumastrn &
&           (kidia  = i_startidx ,                kfdia  = i_endidx           ,& !> IN
&            klon   = nproma ,     ktdia  = kstart_moist(jg)  , klev = nlev   ,& !! IN
&            ldland =  landseemask        (:,jb), ptsphy = tcall_conv_jg      ,& !! IN
!&           ldland = ext_data%atm%lsm_atm_c(:,jb), ptsphy = tcall_phy_jg(itconv)  ,& !! IN
&            pten   = p_diag%temp      (:,:,jb)                              ,& !! IN
&            pqen   = p_prog_rcf%tracer(:,:,jb,iqv)                          ,& !! IN
&            puen   = p_diag%u         (:,:,jb), pven   = p_diag%v( :,:,jb) ,& !! IN
&            plitot = z_plitot          (:,:,jb), pvervel= z_omega_p (:,:,jb) ,& !! IN
&            pqhfl  = z_qhfl            (:,:,jb), pahfs  = z_shfl    (:,:,jb) ,& !! IN
&            pap    = p_diag%pres      (:,:,jb), paph   = p_diag%pres_ifc(:,:,jb),& !! IN
&            pgeo   = p_metrics%geopot_agl (:,:,jb)                        ,& !! IN
&            pgeoh  = p_metrics%geopot_agl_ifc(:,:,jb)                        ,& !! IN
&            zdph   = p_diag%dpres_mc    (:,:,jb)                            ,& !! IN
&            zdgeoh = p_metrics%dgeopot_mc(:,:,jb)                            ,& !! IN
&            ptent  = z_dtdt                       (:,:,jb)                   ,& !! INOUT
&            ptenu  = prm_nwp_tend%ddt_u_pconv     (:,:,jb)                   ,& !! OUT
&            ptenv  = prm_nwp_tend%ddt_v_pconv     (:,:,jb)                   ,& !! OUT
&            ptenq  = z_dtdqv                      (:,:,jb)                   ,& !! INOUT
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
&            pmflxr =      z_mflxr             (:,:,jb)                       ,& !! OUT
&            pmflxs =      z_mflxs             (:,:,jb)                       ,& !! OUT
&            prain  =      prm_diag%rain_upd (:,jb)                           ,& !! OUT
&            pcape =       prm_diag%cape     (:,jb)                           ) !! OUT

!#else
!          CALL cumastrn &
!&           (kidia  = i_startidx ,                kfdia  = i_endidx           ,& !> IN
!&            klon   = nproma ,     ktdia  = kstart_moist(jg)  , klev = nlev   ,& !! IN
!&            ldland =  landseemask        (1,jb), ptsphy = tcall_conv_jg      ,& !! IN
!!&           ldland = ext_data%atm%lsm_atm_c(1,jb), ptsphy = tcall_phy_jg(itconv)  ,& !! IN
!&            pten   = p_diag%temp      (1,1,jb)                              ,& !! IN
!&            pqen   = p_prog_rcf%tracer(1,1,jb,iqv)                          ,& !! IN
!&            puen   = p_diag%u         (1,1,jb), pven   = p_diag%v( 1,1,jb) ,& !! IN
!&            plitot = z_plitot          (1,1,jb), pvervel= z_omega_p (1,1,jb) ,& !! IN
!&            pqhfl  = z_qhfl            (1,1,jb), pahfs  = z_shfl    (1,1,jb) ,& !! IN
!&            pap    = p_diag%pres      (1,1,jb), paph   = p_diag%pres_ifc(1,1,jb),& !! IN
!&            pgeo   = p_metrics%geopot_agl(1,1,jb)                        ,& !! IN
!&            pgeoh  = p_metrics%geopot_agl_ifc(1,1,jb)                        ,& !! IN
!&            zdph   = p_diag%dpres_mc    (1,1,jb)                            ,& !! IN
!&            zdgeoh = p_metrics%dgeopot_mc(1,1,jb)                            ,& !! IN
!&            ptent  = z_dtdt                       (1,1,jb)                   ,& !! INOUT
!&            ptenu  = prm_nwp_tend%ddt_u_pconv     (1,1,jb)                   ,& !! OUT
!&            ptenv  = prm_nwp_tend%ddt_v_pconv     (1,1,jb)                   ,& !! OUT
!&            ptenq  = z_dtdqv                      (1,1,jb)                   ,& !! INOUT
!&            ptenl  = prm_nwp_tend%ddt_tracer_pconv(1,1,jb,iqc)               ,& !! OUT
!&            pteni  = prm_nwp_tend%ddt_tracer_pconv(1,1,jb,iqi)               ,& !! OUT
!!&            ptens  = prm_nwp_tend%ddt_tracer_pconv(1,1,jb,iqs)               ,& !! OUT
!&            ldcum  = prm_diag%locum  (1,jb)                                  ,& !! OUT
!&            ktype  = prm_diag%ktype   (1,jb)                                 ,& !! OUT
!&            kcbot  = prm_diag%mbas_con(1,jb)                                 ,& !! OUT
!&            kctop  = prm_diag%mtop_con(1,jb)                                 ,& !! OUT
!&            pmfu   =      prm_diag%con_udd(1,1,jb,1)                         ,& !! OUT
!&            pmfd   =      prm_diag%con_udd(1,1,jb,2)                         ,& !! OUT
!&            pmfude_rate = prm_diag%con_udd(1,1,jb,3)                         ,& !! OUT
!&            pmfdde_rate = prm_diag%con_udd(1,1,jb,4)                         ,& !! OUT
!&            ptu    =      prm_diag%con_udd(1,1,jb,5)                         ,& !! OUT
!&            pqu    =      prm_diag%con_udd(1,1,jb,6)                         ,& !! OUT
!&            plu    =      prm_diag%con_udd(1,1,jb,7)                         ,& !! OUT
!&            pmflxr =      z_mflxr             (1,1,jb)                       ,& !! OUT
!&            pmflxs =      z_mflxs             (1,1,jb)                       ,& !! OUT
!&            prain  =      prm_diag%rain_upd (1,jb)                           ,& !! OUT
!&            pcape =       prm_diag%cape     (1,jb)                           ) !! OUT

!#endif

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
          prm_nwp_tend%ddt_temp_pconv  (i_startidx:i_endidx,1:kstart_moist(jg)-1,jb) = 0._wp

          prm_nwp_tend%ddt_temp_pconv  (i_startidx:i_endidx,kstart_moist(jg):,jb) =   &
&           z_dtdt(i_startidx:i_endidx,kstart_moist(jg):,jb)                          &
&           -  prm_nwp_tend%ddt_temp_radsw (i_startidx:i_endidx,kstart_moist(jg):,jb) &
&           -  prm_nwp_tend%ddt_temp_radlw (i_startidx:i_endidx,kstart_moist(jg):,jb)

          prm_nwp_tend%ddt_tracer_pconv(i_startidx:i_endidx,kstart_moist(jg):,jb,iqv) = &
&           z_dtdqv(i_startidx:i_endidx,kstart_moist(jg):,jb)                           &
&           - p_diag%ddt_tracer_adv(i_startidx:i_endidx,kstart_moist(jg):,jb,iqv)

          prm_diag%tracer_rate(i_startidx:i_endidx,jb,3) = z_mflxr(i_startidx:i_endidx,nlevp1,jb)
          prm_diag%tracer_rate(i_startidx:i_endidx,jb,4) = z_mflxs(i_startidx:i_endidx,nlevp1,jb)


          ENDIF !inwp_conv

      ENDDO
!$OMP END DO
!$OMP END PARALLEL      

  END SUBROUTINE nwp_convection

END MODULE mo_nwp_conv_interface

