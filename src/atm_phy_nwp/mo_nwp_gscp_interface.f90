!OPTION! -cont
!! this command should fix the problem of copying arrays in a subroutine call
!>
!! This module is the interface between nwp_nh_interface to the 
!! microphysical parameterisations:
!! inwp_gscp == 1 == one_moment bulk microphysics by Doms and Schaettler(2004) 
!!                                                and Seifert and Beheng(2006)
!! inwp_gscp == 2 == ...
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nwp_gscp_interface

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message, message_text
  USE mo_parallel_config,      ONLY: nproma

  USE mo_model_domain,         ONLY: t_patch
  USE mo_impl_constants,       ONLY: min_rlcell_int
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c

  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nonhydrostatic_config,ONLY: kstart_moist
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_run_config,           ONLY: msg_level, iqv, iqc, iqi, iqr, iqs,  &
                                     iqni, iqni_nuc, iqg, iqh, iqnr, iqns,&
                                     iqng, iqnh    
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_gscp_cosmo,           ONLY: hydci_pp_old, kessler_pp
  USE gscp_hydci_pp,           ONLY: hydci_pp
  USE gscp_hydci_pp_ice,       ONLY: hydci_pp_ice
  USE mo_exception,            ONLY: finish
  USE mo_mcrph_sb,             ONLY: two_moment_mcrph
  USE mo_sync,                 ONLY: global_max


  IMPLICIT NONE

  PRIVATE

  PUBLIC  ::  nwp_microphysics

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE nwp_microphysics( tcall_gscp_jg,                & !>input
                            &   p_patch,p_metrics,            & !>input
                            &   p_prog,                       & !>inout
                            &   p_prog_rcf,                   & !>inout
                            &   p_diag ,                      & !>inout
                            &   prm_diag                      ) !>inout 



    TYPE(t_patch),        TARGET,INTENT(in)   :: p_patch        !!<grid/patch info.
    TYPE(t_nh_metrics)          ,INTENT(in)   :: p_metrics
    TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog          !<the dyn prog vars
    TYPE(t_nh_prog),      TARGET,INTENT(inout):: p_prog_rcf      !<call freq
    TYPE(t_nh_diag),      TARGET,INTENT(inout):: p_diag          !<the dyn diag vars
    TYPE(t_nwp_phy_diag),        INTENT(inout):: prm_diag        !<the atm phys vars

    REAL(wp),                    INTENT(in)   :: tcall_gscp_jg   !< time interval for 
                                                                 !< microphysics
    ! Local array bounds:

    INTEGER :: nlev, nlevp1            !< number of full levels !CK<
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !< blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index

    ! Local scalars:

    INTEGER :: jc,jb,jg               !<block indeces
    REAL(wp) :: qvmax, qcmax, qrmax, qimax, qsmax, qhmax, qgmax, tmax, wmax

    ! local variables
    !
    i_nchdom  = MAX(1,p_patch%n_childdom)

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1 !CK<
    
    ! domain ID
    jg = p_patch%id

    ! exclude boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    !Some run time diagnostics for two moment scheme. It would be better
    !to have similar diagnostics for other schemes as well
    IF (msg_level>15.AND.atm_phy_nwp_config(jg)%inwp_gscp==4) THEN

      !local max val 
      tmax  = MAXVAL(p_diag%temp(:,:,:))
      qvmax = MAXVAL(p_prog_rcf%tracer(:,:,:,iqv))
      qcmax = MAXVAL(p_prog_rcf%tracer(:,:,:,iqc))
      qrmax = MAXVAL(p_prog_rcf%tracer(:,:,:,iqr))
      qsmax = MAXVAL(p_prog_rcf%tracer(:,:,:,iqs))
      qimax = MAXVAL(p_prog_rcf%tracer(:,:,:,iqi))
      qhmax = MAXVAL(p_prog_rcf%tracer(:,:,:,iqh))
      qgmax = MAXVAL(p_prog_rcf%tracer(:,:,:,iqg))
      wmax  = MAXVAL(p_prog%w(:,:,:))

      !now calculate global max
      tmax  = global_max(tmax)
      qvmax = global_max(qvmax)
      qcmax = global_max(qcmax)
      qrmax = global_max(qrmax)
      qsmax = global_max(qsmax)
      qimax = global_max(qimax)
      qgmax = global_max(qgmax)
      qhmax = global_max(qhmax)
      wmax  = global_max(wmax)

      CALL message ("mo_gscp_interface:", "output max values before microphysics")
      WRITE(message_text,'(A10,9A11)')   '   ', 'w','qv','qc','qr','qi','qs','qg','qh','tmax'
      CALL message("",TRIM(message_text))
      WRITE(message_text,'(A10,9E11.3)') '   ', wmax,qvmax,qcmax,qrmax,qimax,qsmax,qgmax,qhmax,tmax
      CALL message("",TRIM(message_text))

    END IF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, rl_start, rl_end)

          !>  prognostic microphysics and precipitation scheme from COSMO
          !!  NOTE: since microphysics is a fast process it is
          !!        allowed to do a sequential updating!

        SELECT CASE (atm_phy_nwp_config(jg)%inwp_gscp)

        CASE(0)  ! no micro physics scheme
          
          WRITE(0,*) "                           "

        CASE(1)  ! COSMO-EU scheme (2-cat ice: cloud ice, snow)
!                ! version unified with COSMO scheme
                 ! unification version: COSMO_V4_23

          CALL hydci_pp (                                   &
            & nvec   =nproma                           ,    & !> in:  actual array size
            & ke     =nlev                             ,    & !< in:  actual array size
            & ivstart=i_startidx                       ,    & !< in:  start index of calculation
            & ivend  =i_endidx                         ,    & !< in:  end index of calculation
            & kstart =kstart_moist(jg)                 ,    & !< in:  vertical start index
            & zdt    =tcall_gscp_jg                    ,    & !< in:  timestep
            & qi0    =atm_phy_nwp_config(jg)%qi0       ,    & 
            & qc0    =atm_phy_nwp_config(jg)%qc0       ,    & 
            & dz     =p_metrics%ddqz_z_full(:,:,jb)    ,    & !< in:  vertical layer thickness
            & t      =p_diag%temp   (:,:,jb)           ,    & !< inout:  temp,tracer,...
            & p      =p_diag%pres   (:,:,jb)           ,    & !< in:  full level pres
            & rho    =p_prog%rho    (:,:,jb  )         ,    & !< in:  density
            & qv     =p_prog_rcf%tracer (:,:,jb,iqv)   ,    & !< inout:  spec. humidity
            & qc     =p_prog_rcf%tracer (:,:,jb,iqc)   ,    & !< inout:  cloud water
            & qi     =p_prog_rcf%tracer (:,:,jb,iqi)   ,    & !< inout:  cloud ice
            & qr     =p_prog_rcf%tracer (:,:,jb,iqr)   ,    & !< inout:  rain water
            & qs     =p_prog_rcf%tracer (:,:,jb,iqs)   ,    & !< inout:  snow
            & prr_gsp=prm_diag%rain_gsp_rate (:,jb)    ,    & !< out: precipitation rate of rain
            & prs_gsp=prm_diag%snow_gsp_rate (:,jb)    ,    & !< out: precipitation rate of snow
            & idbg=msg_level/2                         ,    &
            & l_cv=.TRUE. )
          

        CASE(2)  ! COSMO-DE (3-cat ice: snow, cloud ice, graupel)

          CALL finish('mo_nwp_gscp_interface', 'Graupel scheme not implemented.')


        CASE(3)  ! improved ice nucleation scheme by C. Koehler based on hydci_pp

          CALL hydci_pp_ice (                                &
            & nvec   =nproma                            ,    & !> in:  actual array size
            & ke     =nlev                              ,    & !< in:  actual array size
            & ke1    =nlevp1                            ,    & !< in:  half model levels (start index is 1)
            & ivstart=i_startidx                        ,    & !< in:  start index of calculation
            & ivend  =i_endidx                          ,    & !< in:  end index of calculation
            & kstart =kstart_moist(jg)                  ,    & !< in:  vertical start index
            & zdt    =tcall_gscp_jg                     ,    & !< in:  timestep
            & qi0    =atm_phy_nwp_config(jg)%qi0        ,    & 
            & qc0    =atm_phy_nwp_config(jg)%qc0        ,    & 
            & dz     =p_metrics%ddqz_z_full(:,:,jb)     ,    & !< in:  vertical layer thickness
            & t      =p_diag%temp   (:,:,jb)            ,    & !< in:  temp,tracer,...
            & p      =p_diag%pres   (:,:,jb)            ,    & !< in:  full level pres
            & w      =p_prog%w(:,:,jb)                  ,    & !< in:  vertical wind speed, half levs (m/s)
            & tke    =p_prog_rcf%tke(:,:,jb)            ,    & !< in:  turbulent kinetik energy
            & rho    =p_prog%rho    (:,:,jb  )          ,    & !< in:  density
            & qv     =p_prog_rcf%tracer (:,:,jb,iqv)    ,    & !< in:  spec. humidity
            & qc     =p_prog_rcf%tracer (:,:,jb,iqc)    ,    & !< in:  cloud water
            & qi     =p_prog_rcf%tracer (:,:,jb,iqi)    ,    & !< in:  cloud ice
            & qni    =p_prog_rcf%tracer (:,:,jb,iqni)   ,    & !< in:  cloud ice number     ( 1/kg)
            & qni_nuc=p_prog_rcf%tracer (:,:,jb,iqni_nuc),   & !< in:  activated ice nuclei ( 1/kg)            
            & qr     =p_prog_rcf%tracer (:,:,jb,iqr)    ,    & !< in:  rain water
            & qs     =p_prog_rcf%tracer (:,:,jb,iqs)    ,    & !< in:  snow
            & prr_gsp=prm_diag%rain_gsp_rate (:,jb)     ,    & !< out: precipitation rate of rain
            & prs_gsp=prm_diag%snow_gsp_rate (:,jb)     ,    & !< out: precipitation rate of snow
            & idbg=msg_level/2                          ,    &
            & l_cv=.TRUE. )

        CASE(4)  ! two-moment scheme 

          CALL two_moment_mcrph(                       &
                       isize  = nproma,                &!in: array size
                       ke     = nlev,                  &!in: end level/array size
                       is     = i_startidx,            &!in: start index
                       ie     = i_endidx,              &!in: end index
                       ks     = kstart_moist(jg),      &!in: start level
                       dt     = tcall_gscp_jg ,        &!in: time step
                       dz     = p_metrics%ddqz_z_full(:,:,jb),  &!in: vertical layer thickness
                       rho    = p_prog%rho(:,:,jb  )       ,    &!in:  density
                       pres   = p_diag%pres(:,:,jb  )      ,    &!in:  pressure
                       qv     = p_prog_rcf%tracer (:,:,jb,iqv), &!inout:sp humidity
                       qc     = p_prog_rcf%tracer (:,:,jb,iqc), &!inout:cloud water
                       qr     = p_prog_rcf%tracer (:,:,jb,iqr), &!inout:rain
                       qnr    = p_prog_rcf%tracer (:,:,jb,iqnr),&!inout:rain droplet number 
                       qi     = p_prog_rcf%tracer (:,:,jb,iqi), &!inout: ice
                       qni    = p_prog_rcf%tracer (:,:,jb,iqni),&!inout: cloud ice number
                       qs     = p_prog_rcf%tracer (:,:,jb,iqs), &!inout: snow 
                       qns    = p_prog_rcf%tracer (:,:,jb,iqns),&!inout: snow number
                       qg     = p_prog_rcf%tracer (:,:,jb,iqg), &!inout: graupel 
                       qng    = p_prog_rcf%tracer (:,:,jb,iqng),&!inout: graupel number
                       qh     = p_prog_rcf%tracer (:,:,jb,iqh), &!inout: hail 
                       qnh    = p_prog_rcf%tracer (:,:,jb,iqnh),&!inout: hail number
                       tk     = p_diag%temp(:,:,jb),            &!inout: temp 
                       w      = p_prog%w(:,:,jb),               &!inout: w
                       prec_r = prm_diag%rain_gsp_rate (:,jb),  &!inout precp rate rain
                       prec_i = prm_diag%ice_gsp_rate (:,jb),   &!inout precp rate ice
                       prec_s = prm_diag%snow_gsp_rate (:,jb),  &!inout precp rate snow
                       prec_g = prm_diag%graupel_gsp_rate (:,jb),&!inout precp rate graupel
                       prec_h = prm_diag%hail_gsp_rate (:,jb),   &!inout precp rate hail
                       msg_level = msg_level                ,    &
                       l_cv=.TRUE.          )    

        CASE(9)  ! Kessler scheme (warm rain scheme)

          CALL kessler_pp (                                  &
            & ie     =nproma                            ,    & ! in:  actual array size
            & ke     =nlev                              ,    & ! in:  actual array size
            & istart =i_startidx                        ,    & ! in:  start index of calculation
            & iend   =i_endidx                          ,    & ! in:  end index of calculation
            & kstart =kstart_moist(jg)                  ,    & ! in:  vertical start index
            & zdt    =tcall_gscp_jg                     ,    & ! in:  timestep
            & qc0    = atm_phy_nwp_config(jg)%qc0       ,    & 
            & dz     =p_metrics%ddqz_z_full(:,:,jb)     ,    & ! in:  vertical layer thickness
            & t      =p_diag%temp   (:,:,jb)            ,    & ! in:  temp,tracer,...
            & p      =p_diag%pres   (:,:,jb)            ,    & ! in:  full level pres
            & rho    =p_prog%rho    (:,:,jb  )          ,    & ! in:  density
            & qv     =p_prog_rcf%tracer (:,:,jb,iqv)    ,    & ! in:  spec. humidity
            & qc     =p_prog_rcf%tracer (:,:,jb,iqc)    ,    & ! in:  cloud water
            & qr     =p_prog_rcf%tracer (:,:,jb,iqr)    ,    & ! in:  rain water
            & prr_gsp=prm_diag%rain_gsp_rate (:,jb)     ,    & ! out: precipitation rate of rain
            & idbg   =msg_level/2                       ,    &
            & l_cv    =.TRUE. )


        CASE(10)  ! COSMO-EU scheme (2-cat ice: cloud ice, snow)
                  ! old version from COSMO V4.14

          CALL hydci_pp_old (                               &
            & ie     =nproma                           ,    & !> in:  actual array size
            & ke     =nlev                             ,    & !< in:  actual array size
            & istart =i_startidx                       ,    & !< in:  start index of calculation
            & iend   =i_endidx                         ,    & !< in:  end index of calculation
            & kstart =kstart_moist(jg)                 ,    & !< in:  vertical start index
            & zdt    =tcall_gscp_jg                    ,    & !< in:  timestep
            & qi0    =atm_phy_nwp_config(jg)%qi0       ,    & 
            & qc0    =atm_phy_nwp_config(jg)%qc0       ,    & 
            & dz     =p_metrics%ddqz_z_full(:,:,jb)    ,    & !< in:  vertical layer thickness
            & t      =p_diag%temp   (:,:,jb)           ,    & !< in:  temp,tracer,...
            & p      =p_diag%pres   (:,:,jb)           ,    & !< in:  full level pres
            & rho    =p_prog%rho    (:,:,jb  )         ,    & !< in:  density
            & qv     =p_prog_rcf%tracer (:,:,jb,iqv)   ,    & !< in:  spec. humidity
            & qc     =p_prog_rcf%tracer (:,:,jb,iqc)   ,    & !< in:  cloud water
            & qi     =p_prog_rcf%tracer (:,:,jb,iqi)   ,    & !< in:  cloud ice
            & qr     =p_prog_rcf%tracer (:,:,jb,iqr)   ,    & !< in:  rain water
            & qs     =p_prog_rcf%tracer (:,:,jb,iqs)   ,    & !< in:  snow
            & prr_gsp=prm_diag%rain_gsp_rate (:,jb)    ,    & !< out: precipitation rate of rain
            & prs_gsp=prm_diag%snow_gsp_rate (:,jb)    ,    & !< out: precipitation rate of snow
            & idbg=msg_level/2                         ,    &
            & l_cv=.TRUE. )





        CASE DEFAULT

          CALL finish('mo_nwp_gscp_interface', 'Unknown cloud physics scheme [1-5].')

        END SELECT



        !-------------------------------------------------------------------------
        !>
        !! Calculate surface precipitation
        !!
        !-------------------------------------------------------------------------
      
        IF(atm_phy_nwp_config(jg)%inwp_gscp==4)THEN

!DIR$ IVDEP
          DO jc =  i_startidx, i_endidx

            prm_diag%rain_gsp(jc,jb) = prm_diag%rain_gsp(jc,jb)        &
   &                                 + tcall_gscp_jg                   &
   &                                 * prm_diag%rain_gsp_rate (jc,jb)

            prm_diag%snow_gsp(jc,jb) = prm_diag%snow_gsp(jc,jb)        &
   &                                 + tcall_gscp_jg                   &
   &                                 * prm_diag%snow_gsp_rate (jc,jb)

            prm_diag%ice_gsp(jc,jb) = prm_diag%ice_gsp(jc,jb)          &
   &                                 + tcall_gscp_jg                   &
   &                                 * prm_diag%ice_gsp_rate (jc,jb)

            prm_diag%hail_gsp(jc,jb) = prm_diag%hail_gsp(jc,jb)        &
   &                                 + tcall_gscp_jg                   &
   &                                 * prm_diag%hail_gsp_rate (jc,jb)

            prm_diag%graupel_gsp(jc,jb) = prm_diag%graupel_gsp(jc,jb)  &
   &                                 + tcall_gscp_jg                   &
   &                                 * prm_diag%graupel_gsp_rate (jc,jb)

            prm_diag%tot_prec(jc,jb) = prm_diag%tot_prec(jc,jb)         &
   &                                 +  tcall_gscp_jg                   &
   &                                 * ( prm_diag%rain_gsp_rate (jc,jb) &
   &                                 +   prm_diag%snow_gsp_rate (jc,jb) &
   &                                 +   prm_diag%ice_gsp_rate (jc,jb)  &
   &                                 +   prm_diag%hail_gsp_rate (jc,jb) &
   &                                 +   prm_diag%graupel_gsp_rate (jc,jb) )

          ENDDO
 
        ELSE !For all other schemes

!DIR$ IVDEP
          DO jc =  i_startidx, i_endidx

            prm_diag%rain_gsp(jc,jb) = prm_diag%rain_gsp(jc,jb)        &
   &                                 + tcall_gscp_jg                   &
   &                                 * prm_diag%rain_gsp_rate (jc,jb)

            prm_diag%snow_gsp(jc,jb) = prm_diag%snow_gsp(jc,jb)        &
   &                                 + tcall_gscp_jg                   &
   &                                 * prm_diag%snow_gsp_rate (jc,jb)

            prm_diag%tot_prec(jc,jb) = prm_diag%tot_prec(jc,jb)         &
   &                                 +  tcall_gscp_jg                   &
   &                                 * ( prm_diag%rain_gsp_rate (jc,jb) &
   &                                 +   prm_diag%snow_gsp_rate (jc,jb) )
          ENDDO

        END IF 

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
 
    !Some run time diagnostics for two moment scheme
    IF (msg_level>15.AND.atm_phy_nwp_config(jg)%inwp_gscp==4) THEN

      !local max val 
      tmax  = MAXVAL(p_diag%temp(:,:,:))
      qvmax = MAXVAL(p_prog_rcf%tracer(:,:,:,iqv))
      qcmax = MAXVAL(p_prog_rcf%tracer(:,:,:,iqc))
      qrmax = MAXVAL(p_prog_rcf%tracer(:,:,:,iqr))
      qsmax = MAXVAL(p_prog_rcf%tracer(:,:,:,iqs))
      qimax = MAXVAL(p_prog_rcf%tracer(:,:,:,iqi))
      qhmax = MAXVAL(p_prog_rcf%tracer(:,:,:,iqh))
      qgmax = MAXVAL(p_prog_rcf%tracer(:,:,:,iqg))
      wmax  = MAXVAL(p_prog%w(:,:,:))

      !now calculate global max
      tmax  = global_max(tmax)
      qvmax = global_max(qvmax)
      qcmax = global_max(qcmax)
      qrmax = global_max(qrmax)
      qsmax = global_max(qsmax)
      qimax = global_max(qimax)
      qgmax = global_max(qgmax)
      qhmax = global_max(qhmax)
      wmax  = global_max(wmax)

      CALL message ("mo_gscp_interface:", "output max values after microphysics")
      WRITE(message_text,'(A10,9A11)')   '   ', 'w','qv','qc','qr','qi','qs','qg','qh','tmax'
      CALL message("",TRIM(message_text))
      WRITE(message_text,'(A10,9E11.3)') '   ', wmax,qvmax,qcmax,qrmax,qimax,qsmax,qgmax,qhmax,tmax
      CALL message("",TRIM(message_text))

    END IF

     
  END SUBROUTINE nwp_microphysics

END MODULE mo_nwp_gscp_interface

