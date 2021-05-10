!+ Module for latent heat nudging 
!-------------------------------------------------------------------------------
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_latent_heat_nudging

!-------------------------------------------------------------------------------
!>
!! Description:
!!   The module "lheat_nudge" performs the latent heat nudging (lhn).
!!   The lhn adds temperature increments to the prognostic variable t
!!   so that the total temperature increase due to latent heat release
!!   in the current timestep corresponds to the amount of analyzed
!!   (or observed) precipitation.
!!   The temperature increments added due to lhn are derived from the 
!!   model heating rate profiles (large scale condensation and convective 
!!   heating) scaled by the ratio of analyzed to modelled precipitation
!!   (total precipitation: rain and snow from large scale and 
!!   convective processes). The analyzed precipitation is based on radar
!!   data merged with the model (total) precipitation fields.
!!    
!!   The module contains as an organizational unit the subroutine
!!   "organize_lhn" which is called from the module organize_assimilation_config(jg).
!!   Further module procedures (subroutines) called by organize_lhn:
!!   -> lhn_obs_prep : reading+preparing the precip radar data
!!      |
!!      |--> lhn_obs_open  : open the radar data file and read general header
!!      |--> lhn_obs_read  : read a record (header + data) from radar data file
!!      |--> distribute_field : distribute field to all PE's
!!
!!   -> lhn_skill_scores   : verification of precipitation model against radar
!!
!!   -> lhn_t_inc  : derivation of temperature increments by scaling of
!!      |           model latent heating profiles 
!!      |--> assimilation_config(jg)%lhn_artif     : apply artificial profile
!!      |--> assimilation_config(jg)%lhn_filt      : vertical filtering of local ttend_lhn profile
!!      |--> assimilation_config(jg)%lhn_limit     : limiting of the ttend_lhn
!!      |--> assimilation_config(jg)%lhn_relax     : horizontal filtering of ttend_lhn
!!           |--> hor_filt
!!                |--> init_horizontal_filtering
!!                |--> horizontal_filtering
!!                |--> exchange_boundaries
!!   -> lhn_q_inc  : adjust humidity (i.e. qv) to new temperature (t+ttend_lhn)
!!
!!   Note: The names of input/output variables/arrays defined only once 
!!         in the module declaration section but used and "filled" by the 
!!         different subroutines are documented in the description parts
!!         of each procedure for clarity.
!!
!===============================================================================

! Modules used:

!USE mo_datetime,                ONLY: t_datetime,print_datetime
USE mtime,                      ONLY: datetime, newDatetime, timedelta, &
                                      PROLEPTIC_GREGORIAN, setCalendar,                       &
                                      newTimedelta, &
                                      OPERATOR(-), OPERATOR (<), OPERATOR(+), OPERATOR(==), OPERATOR(*), &
                                      assignment(=), OPERATOR (>=), OPERATOR (>), datetimeToString
USE mtime,                       ONLY: timedeltaToString, MAX_TIMEDELTA_STR_LEN, &
  &                                    getPTStringFromSeconds, deallocateTimedelta, deallocateDatetime
!USE mo_mtime_extensions,        ONLY: get_datetime_string

USE mo_kind,               ONLY: wp, i4, i8

USE mo_parallel_config,    ONLY: nproma

USE mo_exception,          ONLY: message, message_text, finish, print_value, open_log, close_log

USE mo_physical_constants, ONLY: r_v   => rv    , & !> gas constant for water vapour
                               rvd_m_o => vtmpc1, & !! rv/rd-1._wp
                                 o_m_rdv        , & !! 1 - r_d/r_v
                                 rdv            , & !! r_d / r_v
                                 cvd            , & !!
                                 cl    => clw   , & !! specific heat of water
                                 lwd   => alv   , & !! latent heat of vapourization
                                 b3    => tmelt , & !!
                                 tmelt

USE mo_convect_tables,     ONLY: b1    => c1es  , & !! constants for computing the sat. vapour
                                 b2w   => c3les , & !! pressure over water (l) and ice (i)
                                 b2i   => c3ies , & !!               -- " --
                                 b4w   => c4les , & !!               -- " --
                                 b4i   => c4ies , & !!               -- " --
                                 b234w => c5les     !!               -- " --

USE mo_math_constants    , ONLY: pi

USE mo_assimilation_config ,ONLY: assimilation_config

USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config

USE mo_loopindices,             ONLY: get_indices_c
USE mo_nonhydro_types,          ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
USE mo_nwp_phy_types,           ONLY: t_nwp_phy_diag, t_nwp_phy_tend
USE mo_impl_constants_grf,      ONLY: grf_bdywidth_c
USE mo_impl_constants,          ONLY: min_rlcell_int, SUCCESS
USE mo_nonhydrostatic_config,   ONLY: kstart_moist
USE mo_model_domain,            ONLY: t_patch
USE mo_radar_data_types,        ONLY: t_radar_fields, t_lhn_diag
USE mo_time_config,             ONLY: time_config
USE mo_mpi,                     ONLY: my_process_is_stdio, p_io
USE mo_io_units,                ONLY: find_next_free_unit
USE mo_run_config,              ONLY: msg_level, iqv, iqc, iqi
USE mo_math_laplace,            ONLY: nabla2_scalar
USE mo_sync,                    ONLY: SYNC_C, sync_patch_array_mult,global_sum
USE mo_intp_data_strc,          ONLY: t_int_state

!===============================================================================

IMPLICIT NONE

PUBLIC :: organize_lhn 

!===============================================================================

! Local scalars:
!---------------

  INTEGER (KIND=i4) ::  &
    nulhn                ! unit of lhn output file

  REAL  (KIND=wp)              ::           &
    zdt                ,& ! valid time step for integration
    rnlhn              ,& ! ration of lhn time step
    sec_per_hr=3600.   ,& ! seconds per hour
    sec_per_hr_inv        ! inverse of seconds per hour

  CHARACTER (LEN=20)    ::  yroutine    ! name of the subroutine
  CHARACTER (LEN=80)    ::  yerrmsg     ! text message for model abort
  CHARACTER (LEN=12)    ::  yulhn       ! name of lhn output file


! Local arrays  
!--------------

  INTEGER,PARAMETER :: ndiag_max=20
  INTEGER (KIND=i4) :: diag_sum(ndiag_max), g_diag_sum(ndiag_max)


!-------------------------------------------------------------------------------
! End of declarations.    Public and private subroutines :
!-------------------------------------------------------------------------------


!===============================================================================

CONTAINS

!===============================================================================
!+ Module procedure "organize_lhn" in "lheat_nudge" to organize the
!  different program steps/subroutines for the latent heat nudging
!-------------------------------------------------------------------------------

SUBROUTINE organize_lhn ( &
                            & dt_loc, p_sim_time,             & !>in
                            & pt_patch, p_metrics,              & !>in
                            & pt_int_state,                     & !>in
                            & pt_prog_rcf,                      & !>inout
                            & pt_diag ,                         & !>inout
                            & prm_diag,                         & !>in
                            & lhn_fields,                       & !>inout
                            & radar_data,                       &
                            & prm_nwp_tend,                     &
                            & datetime_current,                 &
                            & ltlhn, ltlhnverif, lvalid_data    ) !>in / out

!-------------------------------------------------------------------------------
!
! Description:
! This subroutine is the main routine of LHN and is called by lmorg. It contains
! the organization of the LHN:
!
! 1. Read and prepare radar precipitation data
! 2. Get total model latent heating profiles: add convective latent heating 
!    contributions to large scale latent heating (in terms of temperature 
!    tendency as K/s)
! 3. Determine total model precipitation rate, analyze precipitation, i.e. 
!    combine model and observation values
! 4. Determine the latent heat nudging temperature increment by
!    scaling the model profiles; apply artificial profile if requested
! 5. Apply the latent heat nudging temperature increment
! 6. Apply a corresponding humidity increment to maintain or produce
!    (nearly) saturation ("humidity enhancement") or reduce qv in case
!    of negative temperature increments (leave rel. hum. unchanged)
! 7. Output of values and profiles at diagnostic grid points
!    each timestep
!
!-------------------------------------------------------------------------------

  TYPE(t_patch),   TARGET, INTENT(inout) :: pt_patch     !<grid/patch info.
  TYPE(t_nh_metrics),      INTENT(in)    :: p_metrics
  TYPE(t_nh_diag), TARGET, INTENT(inout) :: pt_diag     !<the diagnostic variables
  TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog_rcf !<the prognostic variables (with
                                                          !< red. calling frequency for tracers!
  TYPE(t_nwp_phy_diag),    INTENT(in)    :: prm_diag
  TYPE(t_radar_fields),    INTENT(in)    :: radar_data
  TYPE(t_lhn_diag),        INTENT(inout)    :: lhn_fields
  TYPE(t_int_state)   ,    INTENT(IN)    :: pt_int_state

  TYPE(t_nwp_phy_tend),TARGET,INTENT(in) :: prm_nwp_tend
  TYPE (datetime),INTENT(in) :: datetime_current

  REAL(KIND=wp), INTENT(IN) :: &
   dt_loc, p_sim_time

  LOGICAL, INTENT(IN)  :: &
    ltlhn,ltlhnverif                 ! true if latent heat nudging is active

  LOGICAL, INTENT(OUT) :: lvalid_data

! Local parameters, scalars, arrays:
!-------------------------------------------------------------------------------
! Local scalars:
!---------------

  LOGICAL                   ::           &
    lopen_log  
 
  INTEGER ::  kqrs(nproma)       ! upper layer with qrs_flux > 0.0

  INTEGER :: jg          ! domain ID
  INTEGER :: jb,jc,jk,i_rlstart, i_rlend, ndiag, iter, nlev
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !< slices

  REAL (KIND=wp) :: zdcoeff

  REAL (KIND=wp)               ::           &
    zdt_1                ! inverse of the timestep for physics ( = 1/dt )


  REAL (KIND=wp), DIMENSION(nproma)  ::  qrsgmax,qrsgthres,ttmin,ttmax,vcoordsum,qrsflux_int

! Local arrays:

  REAL (KIND=wp) ::       &
    zprmod       (nproma,pt_patch%nblks_c)  ,&
    zprmod_ref   (nproma,pt_patch%nblks_c)  ,&
    zprrad       (nproma,pt_patch%nblks_c)  ,&
    zprmod_ref_f (nproma,pt_patch%nblks_c)  ,&
    zprrad_f     (nproma,pt_patch%nblks_c)

  REAL  (KIND=wp) ::           &
    pr_obs(nproma,pt_patch%nblks_c)     ,& ! observed (radar) precipitation rate         (kg/m2*s)
    pr_mod(nproma,pt_patch%nblks_c)     ,& ! total model precipitation rate              (kg/m2*s)
    pr_ref(nproma,pt_patch%nblks_c)     ,& ! total reference precipitation rate              (kg/m2*s)
    pr_ana(nproma,pt_patch%nblks_c)     ,& ! analyzed precipitation rate                 (kg/m2*s)
    pr_mod_nofilt(nproma,pt_patch%nblks_c),& !
    pr_obs_nofilt(nproma,pt_patch%nblks_c),& !
    z_pr_mod(nproma,1,pt_patch%nblks_c),& !
    z_pr_obs(nproma,1,pt_patch%nblks_c),& !
    z_nabla2_prmod(nproma,pt_patch%nlev,pt_patch%nblks_c),& !
    z_nabla2_probs(nproma,pt_patch%nlev,pt_patch%nblks_c),& !
    z_nabla2_ttlh(nproma,pt_patch%nlev,pt_patch%nblks_c),& !
    wobs_space(nproma,pt_patch%nblks_c) ,& ! weights (spatial) for the precip obs          ( 1 )
    wobs_time(nproma,pt_patch%nblks_c)  ,& ! weights (temporal) for the precip obs         ( 1 )
    lhn_diag(nproma,pt_patch%nlev-15:pt_patch%nlev,pt_patch%nblks_c)   ,& ! array for test output of diverse 2D fields
    tt_lheat(nproma,pt_patch%nlev,pt_patch%nblks_c)   ,& ! tt_lheat
    qrsflux(nproma,pt_patch%nlev,pt_patch%nblks_c)    ,& ! qrsflux
    scale_diag(nproma,pt_patch%nblks_c)   ,& ! global distribution of scale_fac
    treat_diag(nproma,pt_patch%nblks_c)  ! ,& ! diagnose of treatment
!    windcor_diag(nproma,pt_patch%nblks_c) ,& ! weight with respect to the mean wind

  INTEGER (KIND=i4)  ::        &
    diag_out(pt_patch%nblks_c,ndiag_max) ! array for exchange between PE's (used by global_values)

  LOGICAL  :: &
    scale_fac_index(nproma,pt_patch%nblks_c)

  LOGICAL :: ltoold, ltoyoung
  INTEGER :: izlocstat        ! error status on allocation of fields

!- End of header
!===============================================================================
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Section 0 : Preliminaries : 
!             Check if lhn should be executed at the current timestep
!             Allocate space for fields; determine dt
!-------------------------------------------------------------------------------
  yroutine='organize_lhn'

  izlocstat  = 0   ! initialization
  jg         = pt_patch%id
  nlev       = pt_patch%nlev
  zdt        = dt_loc

  zdt_1 = 1.0_wp/zdt
  sec_per_hr_inv = 1.0_wp/sec_per_hr

  IF (msg_level > 12) THEN
     WRITE(message_text,'(a,f10.2,3i10)' ) 'intent(in) parameter: ', p_sim_time, nproma, pt_patch%nlev,jg
     CALL message(yroutine,message_text)
  ENDIF

  rnlhn      = (p_sim_time)/REAL(assimilation_config(jg)%nlhn_end)

  ! settings to exclude boundary interpolation zone of nested domains
  i_rlstart = grf_bdywidth_c+1
  i_rlend   = min_rlcell_int

  i_startblk = pt_patch%cells%start_block(i_rlstart)
  i_endblk   = pt_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb)
  DO jb = 1, pt_patch%nblks_c
    pr_obs(:,jb) = -0.1_wp
    pr_ana(:,jb) = 0.0_wp
    wobs_space(:,jb) = -1.0_wp
    wobs_time(:,jb) = -1.0_wp
    lhn_diag(:,:,jb) = -99.0_wp
    scale_diag(:,jb) = 0.0_wp
    treat_diag(:,jb) = 0.0_wp
    diag_out(jb,:) = 0
    scale_fac_index(:,jb) =.FALSE.
    lhn_fields%ttend_lhn(:,:,jb) = 0.0_wp
    lhn_fields%qvtend_lhn(:,:,jb) = 0.0_wp
  END DO
!$OMP END DO

!$OMP DO PRIVATE(jb)
  DO jb = 1, i_startblk    ! initialization along nest boundaries
    pr_mod(:,jb) = 0.0_wp
    pr_ref(:,jb) = 0.0_wp
    tt_lheat(:,:,jb) = 0._wp
  END DO
!$OMP END DO

  ! diagnose freezing level - this is otherwise done only at output times
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, i_rlstart, i_rlend)

    prm_diag%hzerocl(i_startidx:i_endidx,jb) = p_metrics%z_ifc(i_startidx:i_endidx,nlev+1,jb)

    DO jk = kstart_moist(jg)+1, nlev
      DO jc = i_startidx, i_endidx 
        IF ( prm_diag%hzerocl(jc,jb) /= p_metrics%z_ifc(jc,nlev+1,jb)) THEN ! freezing level found
          CYCLE
        ELSE IF (pt_diag%temp(jc,jk-1,jb) < tmelt .AND. pt_diag%temp(jc,jk,jb) >= tmelt) THEN
          prm_diag%hzerocl(jc,jb) = p_metrics%z_mc(jc,jk-1,jb) -            &
         &      ( p_metrics%z_mc(jc,jk-1,jb) - p_metrics%z_mc(jc,jk,jb) )*  &
         &      (    pt_diag%temp(jc,jk-1,jb) - tmelt ) /                   &
         &      (    pt_diag%temp(jc,jk-1,jb) - pt_diag%temp(jc,jk,jb) )
        END IF
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  IF (my_process_is_stdio() .AND. (assimilation_config(jg)%lhn_diag) ) THEN
     INQUIRE (file=yulhn,OPENED=lopen_log)
     IF (.NOT. lopen_log ) THEN
       CALL open_lhn_log()
       WRITE (nulhn,'(a,f10.2,3i10)' ) 'LHN : intent(in) parameter: ', p_sim_time, nproma, pt_patch%nlev,jg
       WRITE (nulhn,'(a,2f10.2)' ) 'LHN : relevant time step/time now : ',zdt,p_sim_time*sec_per_hr_inv
       WRITE(nulhn, *)' parameters set for LHN :'
       WRITE(nulhn, *)' Climatological Profile enable : assimilation_config(jg)%lhn_artif = ',assimilation_config(jg)%lhn_artif
       WRITE(nulhn, *)' Vertical Filtering of increments  : assimilation_config(jg)%lhn_filt   = ',assimilation_config(jg)%lhn_filt
       WRITE(nulhn, *)' Horizontal Filtering of increments : assimilation_config(jg)%lhn_relax  = ', &
                        assimilation_config(jg)%lhn_relax,assimilation_config(jg)%nlhn_relax
       WRITE(nulhn, *)' Absolute limit of incs.  :  assimilation_config(jg)%lhn_limit = ',assimilation_config(jg)%lhn_limit,&
                        assimilation_config(jg)%abs_lhn_lim,' (K/second)'
       WRITE(nulhn, *)' Absolute limit of incs.  :  assimilation_config(jg)%lhn_limitp = ',assimilation_config(jg)%lhn_limitp,&
                        assimilation_config(jg)%abs_lhn_lim,' (K/second)'
       WRITE(nulhn, *)' Humidity enhancement :     assimilation_config(jg)%lhn_hum_adj = ',assimilation_config(jg)%lhn_hum_adj
       WRITE(nulhn, *)' Diagnostic output :        assimilation_config(jg)%lhn_diag    = ',assimilation_config(jg)%lhn_diag
     ENDIF  
  ENDIF  



!-------------------------------------------------------------------------------
! Section 1 : Read and prepare radar precipitation data 
!             - Determine if new data have to be read
!             - Read new data, project the observation data onto the model grid
!             - Distribute gridded observations to the PE's (if parallel)
!             - Determine spatial and temporal weights for observations 
!               (includes interpolation in time between consecutive obs)
!-------------------------------------------------------------------------------

     CALL lhn_obs_prep (pt_patch,radar_data,prm_diag,lhn_fields,pr_obs(:,:),wobs_space(:,:),wobs_time(:,:), &
       &                ltoyoung,ltoold,datetime_current)

     lvalid_data = .NOT. ltoold
     IF (ltoold.OR.ltoyoung) RETURN


!-------------------------------------------------------------------------------
! Section 2 : Determine total model precipitation rate
!             Analyze precipitation, i.e. combine model and observation values
!-------------------------------------------------------------------------------

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
          DO jb = i_startblk, i_endblk

            CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, i_rlstart, i_rlend)

            SELECT CASE (atm_phy_nwp_config(jg)%inwp_gscp)
            CASE(1,3)

              DO jc = i_startidx, i_endidx
                pr_mod(jc,jb) = prm_diag%rain_gsp_rate(jc,jb) + prm_diag%snow_gsp_rate(jc,jb)
              ENDDO

            CASE(2)

              DO jc = i_startidx, i_endidx
                pr_mod(jc,jb) = prm_diag%rain_gsp_rate(jc,jb) + prm_diag%snow_gsp_rate(jc,jb) + prm_diag%graupel_gsp_rate(jc,jb)
              ENDDO

            CASE(4,5) ! 2-mom schemes

              DO jc = i_startidx, i_endidx
                pr_mod(jc,jb) = prm_diag%rain_gsp_rate(jc,jb) + prm_diag%snow_gsp_rate(jc,jb) + &
                                prm_diag%graupel_gsp_rate(jc,jb) + prm_diag%hail_gsp_rate(jc,jb)
              ENDDO

            CASE(9)

              DO jc = i_startidx, i_endidx
                pr_mod(jc,jb) = prm_diag%rain_gsp_rate(jc,jb)
              ENDDO

            END SELECT

            DO jk = kstart_moist(jg), nlev
              DO jc = i_startidx, i_endidx
                qrsflux(jc,jk,jb) = prm_diag%qrs_flux(jc,jk,jb)
              ENDDO
            ENDDO

            IF (atm_phy_nwp_config(jg)%inwp_convection > 0) THEN

              DO jc = i_startidx, i_endidx
                pr_mod(jc,jb) = pr_mod(jc,jb) + prm_diag%rain_con_rate_3d (jc,nlev,jb) + prm_diag%snow_con_rate_3d (jc,nlev,jb)
              ENDDO

              DO jk = kstart_moist(jg), nlev
                DO jc = i_startidx, i_endidx
                  qrsflux(jc,jk,jb) = qrsflux(jc,jk,jb) &
                       + prm_diag%rain_con_rate_3d(jc,jk,jb) + prm_diag%snow_con_rate_3d(jc,jk,jb)
                ENDDO
              ENDDO
            ENDIF

            pr_ref(:,jb) = pr_mod(:,jb)
            zprmod(:,jb) = pr_mod(:,jb)

!-------------------------------------------------------------------------------
! Section 3: Get total model latent heating profiles
!            (in terms of temperature tendency as K/s)
!-------------------------------------------------------------------------------

            IF (atm_phy_nwp_config(jg)%inwp_convection > 0) THEN
              DO jk = kstart_moist(jg), nlev
                DO jc = i_startidx, i_endidx
                  tt_lheat(jc,jk,jb) = prm_diag%tt_lheat(jc,jk,jb)*zdt_1 + prm_nwp_tend%ddt_temp_pconv(jc,jk,jb)
                ENDDO
              ENDDO
            ELSE
              DO jk = kstart_moist(jg), nlev
                DO jc = i_startidx, i_endidx
                  tt_lheat(jc,jk,jb) = prm_diag%tt_lheat(jc,jk,jb)*zdt_1
                ENDDO
              ENDDO
            ENDIF

          ENDDO
!$OMP END DO 
!$OMP END PARALLEL


! ------------------------------------------------------------------------------
! Section 4: get reference precipition for comparison of radar and model
! ------------------------------------------------------------------------------

!   take the vertical integral of the precipitation flux as reference.
!   It is computed in src_gscp.hydci_pp or src_gscp.hydci_pp_gr

!$OMP PARALLEL

   IF (assimilation_config(jg)%lhn_qrs) THEN

!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,qrsflux_int,qrsgmax,qrsgthres,vcoordsum,kqrs) ICON_OMP_GUIDED_SCHEDULE
     DO jb=i_startblk,i_endblk

       CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, i_rlstart, i_rlend)

       qrsgmax(i_startidx:i_endidx) = 0._wp
       DO jk = kstart_moist(jg), nlev
         DO jc=i_startidx,i_endidx
           qrsgmax(jc) = MAX(qrsgmax(jc),qrsflux(jc,jk,jb))
         ENDDO
       ENDDO
       DO jc=i_startidx,i_endidx
         qrsgthres(jc) = MAX(assimilation_config(jg)%thres_lhn,assimilation_config(jg)%rqrsgmax*qrsgmax(jc))
         qrsflux_int(jc) = 0.0_wp
         vcoordsum(jc)=0.0_wp
         kqrs(jc)=nlev+1
       ENDDO

       DO jk=kstart_moist(jg),nlev
         DO jc=i_startidx,i_endidx
           IF (qrsgmax(jc) >= qrsgthres(jc)) THEN
             IF (qrsflux(jc,jk,jb) >= qrsgthres(jc) .AND. kqrs(jc)==nlev+1) then
               kqrs(jc)=jk
             ENDIF
           ENDIF
         ENDDO
       ENDDO
       DO jk=MINVAL(kqrs(i_startidx:i_endidx)),nlev
         DO jc=i_startidx,i_endidx
           IF (jk >= kqrs(jc)) THEN
             qrsflux_int(jc) = qrsflux_int(jc) + qrsflux(jc,jk,jb)  &
                              * (p_metrics%z_ifc(jc,jk,jb)-p_metrics%z_ifc(jc,jk+1,jb))
             vcoordsum(jc)=vcoordsum(jc)+(p_metrics%z_ifc(jc,jk,jb)-p_metrics%z_ifc(jc,jk+1,jb))
           ENDIF
         ENDDO
       ENDDO
       DO jc=i_startidx,i_endidx
         IF (vcoordsum(jc) /= 0.0_wp) qrsflux_int(jc) = qrsflux_int(jc) / vcoordsum(jc)
         pr_ref(jc,jb) = qrsflux_int(jc)
       ENDDO

     ENDDO
!$OMP END DO 
   ENDIF

!$OMP DO PRIVATE(jb)
   DO jb=i_startblk,i_endblk
     pr_obs_nofilt(:,jb) = pr_obs(:,jb)
     pr_mod_nofilt(:,jb) = pr_ref(:,jb)
   ENDDO
!$OMP END DO

!$OMP END PARALLEL


   IF (assimilation_config(jg)%lhn_relax) THEN

      z_pr_obs(:,1,:) = pr_obs(:,:)
      z_pr_mod(:,1,:) = pr_ref(:,:)

      CALL sync_patch_array_mult(SYNC_C, pt_patch, 3, z_pr_mod, z_pr_obs, tt_lheat)

      zdcoeff = 0.05_wp ! diffusion coefficient for nabla2 diffusion
  
      DO iter = 1, assimilation_config(jg)%nlhn_relax ! perform niter iterations
                         ! note: a variable number of iterations (with an exit condition) potentially 
                         ! causes trouble with MPI reproducibility
  
        CALL nabla2_scalar(z_pr_mod, pt_patch, pt_int_state, z_nabla2_prmod, 1, 1, grf_bdywidth_c+1, min_rlcell_int)
        CALL nabla2_scalar(z_pr_obs, pt_patch, pt_int_state, z_nabla2_probs, 1, 1, grf_bdywidth_c+1, min_rlcell_int)
        CALL nabla2_scalar(tt_lheat, pt_patch, pt_int_state, z_nabla2_ttlh,        &
                           kstart_moist(jg), nlev, grf_bdywidth_c+1, min_rlcell_int)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk

          CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
                             i_startidx, i_endidx, i_rlstart, i_rlend)

          DO jc = i_startidx, i_endidx
            pr_ref(jc,jb) = MAX(0.0_wp, z_pr_mod(jc,1,jb) + zdcoeff *           &
                            pt_patch%cells%area(jc,jb) * z_nabla2_prmod(jc,1,jb))
   
            pr_obs(jc,jb) = z_pr_obs(jc,1,jb) + zdcoeff *           &
                            pt_patch%cells%area(jc,jb) * z_nabla2_probs(jc,1,jb)
          ENDDO

          DO jk = kstart_moist(jg),pt_patch%nlev
            DO jc = i_startidx, i_endidx
              tt_lheat(jc,jk,jb) = tt_lheat(jc,jk,jb) + zdcoeff *                    &
                                  pt_patch%cells%area(jc,jb) * z_nabla2_ttlh(jc,jk,jb)
            ENDDO
          ENDDO

        ENDDO
!$OMP END DO 
!$OMP END PARALLEL
   
        z_pr_mod(:,1,:) = pr_ref(:,:)
        z_pr_obs(:,1,:) = pr_obs(:,:)

        CALL sync_patch_array_mult(SYNC_C, pt_patch, 3, z_pr_mod, z_pr_obs, tt_lheat)

        pr_ref(:,:) = z_pr_mod(:,1,:)
        pr_obs(:,:) = z_pr_obs(:,1,:)
      ENDDO

! Clipping negative values and reset to original radar domain
      WHERE (pr_obs        < 0.0_wp ) pr_obs =  0.0_wp
      WHERE (pr_obs_nofilt < 0.0_wp ) pr_obs = -1.0_wp

   ENDIF
  
 IF (ltlhnverif) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb)
  DO jb = i_startblk, i_endblk
    zprmod_ref  (:,jb) = pr_mod_nofilt(:,jb)
    zprrad      (:,jb) = pr_obs_nofilt(:,jb)
    zprmod_ref_f(:,jb) = pr_ref       (:,jb)
    zprrad_f    (:,jb) = pr_obs       (:,jb)
  END DO
!$OMP END DO
!$OMP END PARALLEL

   CALL lhn_verification ('SW',pt_patch,radar_data,lhn_fields,p_sim_time,wobs_space,&
                          zprmod,zprmod_ref,zprrad,zprmod_ref_f,zprrad_f)
 ENDIF


!-------------------------------------------------------------------------------
! Section 7 : Determine the latent heat nudging temperature increment by
!             scaling the model profiles; apply artificial profile if requested
!             (This step needs communication between nodes to exchange the
!             heating profiles that are found outside the PE's subdomain)
!-------------------------------------------------------------------------------

   IF (izlocstat /= 0) THEN
     yerrmsg =' ERROR  *** allocation of space for lhn - fields failed'
     CALL finish (yroutine,yerrmsg)
   ENDIF
   scale_fac_index = .FALSE.



   IF (assimilation_config(jg)%lhn_diag) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
     DO jb=i_startblk,i_endblk
       CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, i_rlstart, i_rlend)
       CALL lhn_t_inc (i_startidx,i_endidx,jg,pt_patch%nlev,p_metrics%z_ifc(:,:,jb), &
                       tt_lheat(:,:,jb),wobs_time(:,jb), wobs_space(:,jb), pr_obs(:,jb), pr_ref(:,jb),pr_ana(:,jb),     &
                       lhn_fields%ttend_lhn(:,:,jb), treat_diag(:,jb),scale_diag(:,jb),scale_fac_index(:,jb),           &
!                       prm_nwp_tend%ddt_temp_pconv(:,:,jb),diag_out(jb,:))
                       pt_diag%u(:,:,jb),pt_diag%v(:,:,jb),prm_diag%k850(:,jb),prm_diag%k950(:,jb),prm_diag%k700(:,jb), &
                       diag_out(jb,:))
     ENDDO
!$OMP END DO 
!$OMP END PARALLEL
    ELSE
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
     DO jb=i_startblk,i_endblk
       CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
          &                i_startidx, i_endidx, i_rlstart, i_rlend)
       CALL lhn_t_inc (i_startidx,i_endidx,jg,pt_patch%nlev,p_metrics%z_ifc(:,:,jb), &
                       tt_lheat(:,:,jb),wobs_time(:,jb), wobs_space(:,jb), pr_obs(:,jb), pr_ref(:,jb),pr_ana(:,jb),     &
                       lhn_fields%ttend_lhn(:,:,jb), treat_diag(:,jb),scale_diag(:,jb),scale_fac_index(:,jb),           &
!                       prm_nwp_tend%ddt_temp_pconv(:,:,jb))
                       pt_diag%u(:,:,jb),pt_diag%v(:,:,jb),prm_diag%k850(:,jb),prm_diag%k950(:,jb),prm_diag%k700(:,jb)  &
                       )
     ENDDO
!$OMP END DO 
!$OMP END PARALLEL
    ENDIF

    IF (assimilation_config(jg)%lhn_diag) THEN
       diag_sum = 0
       g_diag_sum = 0
       DO ndiag=1,15
         diag_sum(ndiag) = SUM(diag_out(:,ndiag))
       ENDDO
       g_diag_sum(1:15) = global_sum(diag_sum(1:15),opt_iroot=p_io) 

      IF (my_process_is_stdio()) THEN
        WRITE(nulhn, *)
        WRITE(nulhn, *)' Diagnostics of LHN - nudging scheme, subroutine lhn_t_inc'

        WRITE(nulhn, *)'Diagnostics of LHN, lhn_t_inc, timestep : ',p_sim_time*zdt_1
        WRITE(nulhn, '(A,L3,f6.2)' )' Latent Heat Nudging active          : ',ltlhn, REAL(diag_out(i_endblk,16))/100.
        WRITE(nulhn, *)' n of points with increments         : ',g_diag_sum(1)
        WRITE(nulhn, *)' n of points with local profiles     : ',g_diag_sum(2)
        WRITE(nulhn, *)' n of points with upscaling          : ',g_diag_sum(3)
        WRITE(nulhn, *)' n of points with limited upscaling  : ',g_diag_sum(4)
        WRITE(nulhn, *)' n of points with downscaling        : ',g_diag_sum(5)
        WRITE(nulhn, *)' n of points with limited downscaling: ',g_diag_sum(6)
        WRITE(nulhn, *)' n of points with artif. prof         : ',g_diag_sum(7)
        WRITE(nulhn, *)' points with imposed positive limit  : ',g_diag_sum(8)
        WRITE(nulhn, *)' points with imposed negative limit  : ',g_diag_sum(9)
        WRITE(nulhn, *)' points with wind weighting < 1 and > 0 : ',g_diag_sum(10)
        WRITE(nulhn, *)' points with wind weighting equal 0     : ',g_diag_sum(11)
        WRITE(nulhn, *)' points with applied in_cloud treatment : ',g_diag_sum(15)
        WRITE(nulhn, *)
        WRITE(nulhn, *)' Vert. Filtering : n points eliminate oscillations: ',g_diag_sum(12)
        WRITE(nulhn, *)' Vert. Filtering : n points eliminate isolate peaks: ',g_diag_sum(13)
        WRITE(nulhn, *)' Vert. Filtering : n points smoothed : ',g_diag_sum(14)
      ENDIF
    ENDIF

!-------------------------------------------------------------------------------
! Section 9 : Apply a corresponding humidity increment to maintain or produce
!             (nearly) saturation ("humidity enhancement") or reduce qv in case
!             of negative temperature increments (leave rel. hum. unchanged)
!             or
!             Apply saturation adjustment to balance the thermodynamic fields
!             after the LHN
!-------------------------------------------------------------------------------


!$OMP PARALLEL

  IF (assimilation_config(jg)%lhn_hum_adj) THEN

!$OMP DO PRIVATE(jb,jc,i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
     DO jb=i_startblk,i_endblk
       CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
         &                i_startidx, i_endidx, i_rlstart, i_rlend)

        CALL lhn_q_inc( &
             i_startidx,i_endidx,jg,zdt,nlev, &
             pt_diag%temp(:,:,jb),lhn_fields%ttend_lhn(:,:,jb), &
             pt_diag%pres(:,:,jb), &
             pt_prog_rcf%tracer (:,:,jb,iqv), &
             pt_prog_rcf%tracer (:,:,jb,iqc), &
             pt_prog_rcf%tracer (:,:,jb,iqi), &
             lhn_fields%qvtend_lhn(:,:,jb), &
             scale_fac_index(:,jb))
     ENDDO
!$OMP END DO 

  ENDIF

!$OMP END PARALLEL

!-------------------------------------------------------------------------------
! Section 10 : Diagnostic procedure...
! a) integrate observed precipitation rates over one hour
! b) control output for some variables via lhn_diag
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Section 10a : integrate observed precipitation rates over one hour
!-------------------------------------------------------------------------------
! for verification integrate observed precipitation rates over one hour:

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,ttmin,ttmax)
  DO jb=i_startblk,i_endblk

    CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
      &                i_startidx, i_endidx, i_rlstart, i_rlend)

   DO jc = i_startidx, i_endidx
     IF (pr_obs(jc,jb) > 0.0_wp)     lhn_fields%pr_obs_sum(jc,jb) = lhn_fields%pr_obs_sum(jc,jb)       &
                                                   + pr_obs(jc,jb) * zdt
     IF (pr_mod(jc,jb) > 0.0_wp)     lhn_fields%pr_mod_sum(jc,jb) = lhn_fields%pr_mod_sum(jc,jb)       &
                                                   + pr_mod(jc,jb) * zdt
     IF (pr_ref(jc,jb) > 0.0_wp)     lhn_fields%pr_ref_sum(jc,jb) = lhn_fields%pr_ref_sum(jc,jb)       &
                                                   + pr_ref(jc,jb) * zdt

!-------------------------------------------------------------------------------
! Section 10b : control output for some variables via lhn_diag
!-------------------------------------------------------------------------------
! control output of pr_mod, pr_obs, pr_ana in lhn_diag (upper 3 levels)
! ouput is in kg/(s*m**2)

     lhn_diag(jc,nlev,jb) = pr_obs(jc,jb)

     lhn_diag(jc,nlev-1,jb) = pr_mod(jc,jb)                 ! ive: 2
     lhn_diag(jc,nlev-2,jb) = pr_ref(jc,jb)                 ! ive: 3
     lhn_diag(jc,nlev-3,jb) = pr_ana(jc,jb)                 ! ive: 4
     lhn_diag(jc,nlev-4,jb) = wobs_space(jc,jb)             ! ive: 5
     lhn_diag(jc,nlev-5,jb) = wobs_time(jc,jb)              ! ive: 6
     IF (pr_obs(jc,jb) > assimilation_config(jg)%thres_lhn .AND. pr_mod(jc,jb) > assimilation_config(jg)%thres_lhn) &
       lhn_diag(jc,nlev-6,jb) = 1.0_wp          ! ive: 7
     IF (pr_obs(jc,jb) > assimilation_config(jg)%thres_lhn .AND. pr_mod(jc,jb) <= assimilation_config(jg)%thres_lhn) &
       lhn_diag(jc,nlev-6,jb) = 2.0_wp          ! ive: 7
     IF (pr_obs(jc,jb) <= assimilation_config(jg)%thres_lhn .AND. pr_mod(jc,jb) > assimilation_config(jg)%thres_lhn) &
       lhn_diag(jc,nlev-6,jb) = 3.0_wp          ! ive: 7
     IF (pr_obs(jc,jb) > pr_mod(jc,jb)) THEN
       lhn_diag(jc,nlev-7,jb) = 1.0_wp          ! ive: 8
     ELSE
       lhn_diag(jc,nlev-7,jb) = 0.0_wp          ! ive: 8
     ENDIF
     IF (pr_obs(jc,jb) < pr_mod(jc,jb)) THEN
       lhn_diag(jc,nlev-8,jb) = -1.0_wp         ! ive: 9
     ELSE
       lhn_diag(jc,nlev-8,jb) = 0.0_wp          ! ive: 9
     ENDIF
     lhn_diag(jc,nlev-10,jb) = lhn_fields%pr_obs_sum(jc,jb)            ! ive: 11
     lhn_diag(jc,nlev-11,jb) = lhn_fields%pr_mod_sum(jc,jb)            ! ive: 12
     lhn_diag(jc,nlev-12,jb) = lhn_fields%pr_ref_sum(jc,jb)            ! ive: 13
     lhn_diag(jc,nlev-13,jb) = scale_diag(jc,jb)       ! ive: 14
     lhn_diag(jc,nlev-14,jb) = treat_diag(jc,jb)            ! ive: 15
     lhn_diag(jc,nlev-15,jb) = lhn_fields%brightband(jc,jb) ! ive: 16
   ENDDO

   ttmin(i_startidx:i_endidx) = 0._wp
   ttmax(i_startidx:i_endidx) = 0._wp
   DO jk = kstart_moist(jg), nlev
     DO jc = i_startidx, i_endidx
       ttmin(jc) = MIN(ttmin(jc),lhn_fields%ttend_lhn(jc,jk,jb))
       ttmax(jc) = MAX(ttmax(jc),lhn_fields%ttend_lhn(jc,jk,jb))

       prm_diag%tt_lheat(jc,jk,jb) = tt_lheat(jc,jk,jb)
       prm_diag%ttend_lhn(jc,jk,jb) = lhn_fields%ttend_lhn(jc,jk,jb)
       prm_diag%qvtend_lhn(jc,jk,jb) = lhn_fields%qvtend_lhn(jc,jk,jb)
     ENDDO
   ENDDO
   DO jc = i_startidx, i_endidx
     lhn_diag(jc,nlev-9,jb) = MERGE(0._wp,1._wp,ttmin(jc)==0._wp .AND. ttmax(jc)==0._wp)
   ENDDO
   prm_diag%lhn_diag(i_startidx:i_endidx,nlev-15:nlev,jb) = lhn_diag(i_startidx:i_endidx,nlev-15:nlev,jb)

  ENDDO

!$OMP END DO 
!$OMP END PARALLEL

   IF (datetime_current%time%minute  == 0) THEN
      IF (ltlhnverif) THEN
         CALL lhn_verification ('HR',pt_patch,radar_data,lhn_fields,p_sim_time/3600.,wobs_space,&
                                lhn_fields%pr_mod_sum, lhn_fields%pr_ref_sum,lhn_fields%pr_obs_sum)
      ENDIF
      lhn_fields%pr_obs_sum(:,:)  = 0.0_wp
      lhn_fields%pr_mod_sum(:,:)  = 0.0_wp
      lhn_fields%pr_ref_sum(:,:)  = 0.0_wp
   ENDIF

!-------------------------------------------------------------------------------
! Section 11 : Deallocate fields needed only during the lhn - step
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! End of module procedure organize_lhn
!-------------------------------------------------------------------------------

END SUBROUTINE organize_lhn

!===============================================================================
!+ Module procedure in "lheat_nudge" preparing radar precip data for lhn
!-------------------------------------------------------------------------------

SUBROUTINE lhn_obs_prep (pt_patch,radar_data,prm_diag,lhn_fields,pr_obs, &
                      &  wobs_space, wobs_time, ltoyoung,ltoold,datetime_current)

!-------------------------------------------------------------------------------
!
! Description:
!   This procedure of the module "lheat_nudge" is called from organize_lhn.
!   It reads in observations (rain rates in mm/h) from the radar data input
!   file and prepares them for further use:
!   The observations are interpolation in time
!   between consecutive observation times is done and spatial and temporal
!   weighting factors are assigned to the resulting observations at each
!   grid point.
!
!   Namelist parameters used: 
!            lhn_black, blacklist_file, radar_in
!   Input arrays : none. (Use of general information on model grid)
!   Output arrays: pr_obs,wobs_space,wobs_time
!
! Method:
!    1: Interpolate observation data in time and assign temporal weight
!    2: Determine the spatial normalized weight
!
!-------------------------------------------------------------------------------

! Subroutine / Function arguments: None
!-------------------------------------------------------------------------------
! Scalar arguments, intent(inout) :
!-------------------------------
!  INTEGER   (KIND=i4), INTENT(IN)     ::       &
!    i_startidx, i_endidx,jg
  TYPE(t_patch),   TARGET, INTENT(in)    :: pt_patch     !<grid/patch info.
  TYPE(t_radar_fields),    INTENT(in)    :: radar_data
  TYPE(t_nwp_phy_diag),    INTENT(in)    :: prm_diag
  TYPE(t_lhn_diag),        INTENT(inout) :: lhn_fields

  REAL (KIND=wp), DIMENSION(nproma,pt_patch%nblks_c),INTENT(OUT)    :: &
    wobs_time, wobs_space, pr_obs

  LOGICAL, INTENT(OUT) :: ltoyoung, ltoold
  TYPE (datetime),INTENT(in) :: datetime_current


! Local parameters, scalars, arrays :
!-------------------------------------------------------------------------------
! Local scalars:
!---------------
! 1. Variables for organizing the code

  INTEGER (KIND=i4), SAVE             ::       &
     i              ,&
     iread          ,& ! number of data fields in time cache
     weight_index_0 ,& !
     weight_index_p1,& !
     weight_index_p2,& !
     weight_index_p3,& !
     weight_index_m1, weight_index_m1lim, & !
     weight_index_m2, weight_index_m2lim   !

  INTEGER :: jb,jc,i_rlstart, i_rlend
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !< slices

  REAL (KIND=wp)                   ::       &
    pr_time_limit  

  INTEGER (KIND=i4)             ::       &
    icenter         ,&
! concerning the time interpolation of DX data
    num_t_obs (nproma,pt_patch%nblks_c,0:4)
! concerning the spatial interpolation of radar data

  TYPE(datetime),  POINTER             ::  &
     tnow           ,&
     center_time    ,& ! time of actual radar observation
     next_time_1    ,& ! time of next observation
     next_time_2    ,& ! time of next observation
     next_time_3    ,& ! time of next observation
     prev_time_1    ,& ! time of previos observation
     prev_time_2       ! time of previos observation

  TYPE(timedelta), POINTER             ::  &
     time_delta     ,& !
     inc_time_p, inc_time_m, &
     titer

  REAL (KIND=wp),ALLOCATABLE :: td_in_min(:)

  INTEGER,PARAMETER :: ndiag=9
!  INTEGER :: n,diag_out(ndiag),g_diag_sum(ndiag)
  INTEGER :: diag_out(ndiag),g_diag_sum(ndiag)

  LOGICAL :: lp1, lp2, lp3, lm1, lm2

  CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN):: lhn_dt_obs_PTstr
  CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN):: zdt_PTstr

  INTEGER :: jg    ! domain ID

  REAL  (KIND=wp) ::           &
    obs_sum(nproma,pt_patch%nblks_c) ,&      ! observed (radar) precipitation rate         (kg/m2*s)
    obs_ratio(nproma,pt_patch%nblks_c) ,&
    obs_sum_g

  INTEGER         ::           &
    obs_cnt(nproma,pt_patch%nblks_c)      ! total model precipitation rate              (kg/m2*s)

  INTEGER :: ns, nsums, nsume, nsum_g

! Local arrays:
!--------------
!- End of header
!===============================================================================
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Section 0: Initialization of some variables
!            At first call: Get start date and time of this model run
!                           Open data file and read general header information
!            At each call : Allocate space for temporal fields
!
!-------------------------------------------------------------------------------
  yroutine='lhn_obs_prep'

  jg = pt_patch%id

  tnow     => newDatetime(datetime_current)
  IF (msg_level > 12) &
   CALL print_value("mod_time (mmddhhmm)",tnow%date%month*1000000+tnow%date%day*10000+tnow%time%hour*100+tnow%time%minute)

  center_time => newDatetime(time_config%tc_current_date)
  next_time_1 => newDatetime(time_config%tc_current_date) ! create datetime pointer
  next_time_2 => newDatetime(time_config%tc_current_date) !0,0,0,0,0,0,0)
  next_time_3 => newDatetime(time_config%tc_current_date) !0,0,0,0,0,0,0)
  prev_time_1 => newDatetime(time_config%tc_current_date) !0,0,0,0,0,0,0)
  prev_time_2 => newDatetime(time_config%tc_current_date) !0,0,0,0,0,0,0)


  CALL getPTStringFromSeconds(INT(assimilation_config(jg)%lhn_dt_obs*60._wp,i8),lhn_dt_obs_PTstr)
  inc_time_p  => newTimedelta(lhn_dt_obs_PTstr)
  !
  ! negative counterpart
  CALL getPTStringFromSeconds(INT(assimilation_config(jg)%lhn_dt_obs*60._wp*(-1._wp),i8),lhn_dt_obs_PTstr)
  inc_time_m  => newTimedelta(lhn_dt_obs_PTstr)
  !
  CALL getPTStringFromSeconds(INT(zdt,i8),zdt_PTstr)
  titer       => newTimedelta(zdt_PTstr)
  !
  time_delta  => newTimedelta('PT0S')


  IF (msg_level > 12) &
   CALL print_value("mod_time_2 (mmddhhmm)",tnow%date%month*1000000+tnow%date%day*10000+tnow%time%hour*100+tnow%time%minute)

  iread=assimilation_config(jg)%nobs_times

  IF (iread <= 1) THEN
    CALL message (yroutine,'too few radar data available')
    RETURN 
  ENDIF

  ALLOCATE (td_in_min(iread))
  td_in_min = -9999.9_wp

  DO i=1,iread
     time_delta = tnow - radar_data%radar_td%obs_date(i) 
     td_in_min(i)=time_delta%second/60.+time_delta%minute+time_delta%hour*60.+time_delta%day*1440.
     if ( radar_data%radar_td%obs_date(i) > tnow) td_in_min(i) = -1.0_wp * td_in_min(i) 
     ! Note: tnow - obs_date gives always positive values in mtime !!!
     ! To be consistent with the linear interpolation below, younger dates are set to negative !
  ENDDO

  ltoold=.false.
  ltoyoung=.false.
  IF (ALL(td_in_min > 0.5*assimilation_config(jg)%lhn_dt_obs) ) THEN
     IF (msg_level > 12) THEN
       CALL message (yroutine,'obsvervations are too old')
       CALL print_value('Max time_diff',MAXVAL(td_in_min))
       CALL print_value('Min time_diff',MINVAL(td_in_min))
     ENDIF
     ltoold=.true.
     RETURN
  ELSE IF (ALL(td_in_min < -0.5*assimilation_config(jg)%lhn_dt_obs)) THEN
     IF (msg_level > 12) THEN
       CALL message (yroutine,'obsvervations are too young')
       CALL print_value('Max time_diff',MAXVAL(td_in_min))
       CALL print_value('Min time_diff',MINVAL(td_in_min))
     ENDIF
     ltoyoung=.true.
     RETURN

  ELSE
     icenter=MINLOC(ABS(td_in_min),1)
     center_time=radar_data%radar_td%obs_date(icenter)

     IF (ABS(td_in_min(icenter)) > 0.5*assimilation_config(jg)%lhn_dt_obs) THEN
       ltoyoung=.true.
       RETURN
     ENDIF
     next_time_1=center_time+inc_time_p
     next_time_2=center_time+2*inc_time_p
     next_time_3=center_time+3*inc_time_p
     prev_time_1=center_time+inc_time_m
     prev_time_2=center_time+2*inc_time_m
     weight_index_0=icenter
     IF (msg_level > 12) CALL print_value ('o_ct',radar_data%radar_td%obs_date(icenter)%time%minute)
     IF (msg_level > 12) CALL print_value ('td_ct',td_in_min(icenter))
     IF (msg_level > 12) CALL print_value ('icenter',icenter)
  ENDIF

  lp1=.FALSE.
  lp2=.FALSE.
  lp3=.FALSE.
  lm1=.FALSE.
  lm2=.FALSE.

  DO i=1,iread
    IF (.NOT.lp1 .AND. radar_data%radar_td%obs_date(i) == next_time_1) THEN 
      weight_index_p1=i
      IF (msg_level > 12) CALL print_value ('o_p1',radar_data%radar_td%obs_date(i)%time%minute)
      IF (msg_level > 12) CALL print_value ('o_p1',weight_index_p1)
      IF (msg_level > 12) CALL print_value ('td_p1',td_in_min(weight_index_p1))
      lp1=.true.
    ELSE IF (.NOT.lp2 .AND. radar_data%radar_td%obs_date(i) == next_time_2) THEN 
      weight_index_p2=i
      IF (msg_level > 12) CALL print_value ('o_p2',radar_data%radar_td%obs_date(i)%time%minute)
      IF (msg_level > 12) CALL print_value ('o_p2',weight_index_p2)
      IF (msg_level > 12) CALL print_value ('td_p2',td_in_min(weight_index_p2))
      lp2=.true.
    ELSE IF (.NOT.lp3 .AND. radar_data%radar_td%obs_date(i) == next_time_3) THEN 
      weight_index_p3=i
      IF (msg_level > 12) CALL print_value ('o_p3',radar_data%radar_td%obs_date(i)%time%minute)
      IF (msg_level > 12) CALL print_value ('o_p3',weight_index_p3)
      IF (msg_level > 12) CALL print_value ('td_p3',td_in_min(weight_index_p3))
      lp3=.true.
    ELSE IF (.NOT.lm1 .AND. radar_data%radar_td%obs_date(i) == prev_time_1) THEN 
      weight_index_m1=i
      IF (msg_level > 12) CALL print_value ('o_m1',radar_data%radar_td%obs_date(i)%time%minute)
      IF (msg_level > 12) CALL print_value ('o_m1',weight_index_m1)
      IF (msg_level > 12) CALL print_value ('td_m1',td_in_min(weight_index_m1))
      lm1=.true.
    ELSE IF (.NOT.lm2 .AND. radar_data%radar_td%obs_date(i) == prev_time_2) THEN
      weight_index_m2=i
      IF (msg_level > 12) CALL print_value ('o_m2',radar_data%radar_td%obs_date(i)%time%minute)
      IF (msg_level > 12) CALL print_value ('o_m2',weight_index_m2)
      IF (msg_level > 12) CALL print_value ('td_m2',td_in_min(weight_index_m2))
      lm2=.true.
    ENDIF
    IF (lp1 .AND. lp2 .AND. lp3 .AND. lm1 .AND. lm2) EXIT
  ENDDO

  weight_index_m1lim = MAX(1,weight_index_m1) ! to avoid errors with array bound checking
  weight_index_m2lim = MAX(1,weight_index_m2) ! to avoid errors with array bound checking

  IF (assimilation_config(jg)%lhn_bright .AND. &
   & ( (MOD(datetime_current%time%minute,5)  == 0 .AND. datetime_current%time%second < 10 ) &
   &   .OR. MAXVAL(lhn_fields%brightband(:,:)) < 0)) &
   & THEN
   ! brightband detection is called at every observation time, ie. every 5 minutes
   ! the hourly precipitation sum is calculated from the hour around icenter time [icenter - 5;icenter + 6]
      ! exclude boundary interpolation zone of nested domains
      i_rlstart = grf_bdywidth_c+1
      i_rlend   = min_rlcell_int

      i_startblk = pt_patch%cells%start_block(i_rlstart)
      i_endblk   = pt_patch%cells%end_block(i_rlend)


     nsums=MAX(1,icenter-5)
     nsume=MIN(iread,nsums+12)

     obs_sum   (:,:) =  0.0_wp
     obs_ratio (:,:) =  0.0_wp
     obs_cnt   (:,:) =  0_i4
     obs_sum_g       =  0.0_wp
     nsum_g          =  0_i4

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,ns) ICON_OMP_GUIDED_SCHEDULE
     DO jb=i_startblk,i_endblk
       CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
         &                i_startidx, i_endidx, i_rlstart, i_rlend)
        DO ns=nsums,nsume
          DO jc=i_startidx,i_endidx
            IF (NINT(radar_data%radar_ct%blacklist(jc,jb)) /= 0_i4 ) CYCLE
            IF (radar_data%radar_td%obs(jc,jb,ns) >= 0.0_wp) THEN
              obs_sum (jc,jb) = obs_sum (jc,jb) + radar_data%radar_td%obs(jc,jb,ns)
              obs_cnt (jc,jb) = obs_cnt (jc,jb) + 1_i4
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO 
!$OMP END PARALLEL

     DO jb=i_startblk,i_endblk
       CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
         &                i_startidx, i_endidx, i_rlstart, i_rlend)
        DO jc=i_startidx,i_endidx
            IF (obs_cnt (jc,jb) <  1_i4 ) CYCLE
            obs_sum (jc,jb) = obs_sum (jc,jb) / REAL(obs_cnt (jc,jb),wp)
            IF ( obs_sum (jc,jb) > 2.5_wp ) THEN
               obs_sum_g = obs_sum_g + obs_sum (jc,jb)
               nsum_g    = nsum_g + 1
            ENDIF
        ENDDO
      ENDDO

      obs_sum_g = global_sum(obs_sum_g)
      nsum_g = global_sum(nsum_g)
      IF (nsum_g < 1 ) THEN
         obs_ratio(:,:) = 0.0
      ELSE
         obs_sum_g = obs_sum_g / nsum_g
         IF ( obs_sum_g > 0._wp ) obs_ratio (:,:) = obs_sum (:,:) / obs_sum_g
      ENDIF
      lhn_fields%brightband(:,:) = 0.0_wp
      CALL detect_bright_band (pt_patch,radar_data,prm_diag,lhn_fields,obs_ratio)
 
  ENDIF

! reset counters
  num_t_obs       = 0

!  print *, MAXVAL(obs)


  pr_time_limit = 0.0_wp
! If the data is in high frequency take into account observations that
! are within the interval [-2,3]*lhn_dt_obs fore the time interpolation
! of obs and wobs_space
  ! exclude boundary interpolation zone of nested domains
  i_rlstart = grf_bdywidth_c+1
  i_rlend   = min_rlcell_int

  i_startblk = pt_patch%cells%start_block(i_rlstart)
  i_endblk   = pt_patch%cells%end_block(i_rlend)

  IF (td_in_min(weight_index_0) >= 0) THEN
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_GUIDED_SCHEDULE
    DO jb=i_startblk,i_endblk
      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
      &                i_startidx, i_endidx, i_rlstart, i_rlend)
      DO jc=i_startidx,i_endidx
        IF (NINT(radar_data%radar_ct%blacklist(jc,jb)) /= 1_i4 .AND. NINT(lhn_fields%brightband(jc,jb)) /= 1_i4) THEN
          IF ((radar_data%radar_td%obs(jc,jb,weight_index_0) >= pr_time_limit) .AND. &
              (lp1 .AND. radar_data%radar_td%obs(jc,jb,weight_index_p1) >= pr_time_limit)) THEN
            ! observation is valid between t>=0 and t<=+lhn_dt_obs
            pr_obs(jc,jb)    = radar_data%radar_td%obs(jc,jb,weight_index_0)                           &
              + (radar_data%radar_td%obs(jc,jb,weight_index_p1)-radar_data%radar_td%obs(jc,jb,weight_index_0)) &
              * (abs(td_in_min(weight_index_0)))/ assimilation_config(jg)%lhn_dt_obs
            pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
            wobs_time(jc,jb) = 1.0_wp
            num_t_obs (jc,jb,1) = 1
            IF (assimilation_config(jg)%lhn_spqual) THEN
              wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_0)                                   &
                + (radar_data%radar_td%spqual(jc,jb,weight_index_p1)-radar_data%radar_td%spqual(jc,jb,weight_index_0)) &
                * (abs(td_in_min(weight_index_0)))/ assimilation_config(jg)%lhn_dt_obs
            ELSE
              wobs_space(jc,jb) = 1.0_wp
            ENDIF
          ELSEIF ((radar_data%radar_td%obs(jc,jb,weight_index_0) >= pr_time_limit) .AND. &
                  (lm1 .AND. radar_data%radar_td%obs(jc,jb,weight_index_m1lim) >= pr_time_limit)) THEN
            ! observation is valid between t>=0 and t<=+lhn_dt_obs
            pr_obs(jc,jb)    = radar_data%radar_td%obs(jc,jb,weight_index_0)                           &
              + (radar_data%radar_td%obs(jc,jb,weight_index_0)-radar_data%radar_td%obs(jc,jb,weight_index_m1lim)) &
              * (abs(td_in_min(weight_index_0)))/ assimilation_config(jg)%lhn_dt_obs
            pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
            wobs_time(jc,jb) = 1.0_wp
            num_t_obs (jc,jb,1) = 1
            IF (assimilation_config(jg)%lhn_spqual) THEN
              wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_0)                                   &
                + (radar_data%radar_td%spqual(jc,jb,weight_index_0)-radar_data%radar_td%spqual(jc,jb,weight_index_m1lim)) &
                * (abs(td_in_min(weight_index_0)))/ assimilation_config(jg)%lhn_dt_obs
            ELSE
              wobs_space(jc,jb) = 1.0_wp
            ENDIF
          ELSEIF ((lm1 .AND. radar_data%radar_td%obs(jc,jb,weight_index_m1lim) >= pr_time_limit) .AND. &
                  (lp1 .AND. radar_data%radar_td%obs(jc,jb,weight_index_p1) >= pr_time_limit)) THEN
            ! observation is valid between t>=-lhn_dt_obs and t<=+lhn_dt_obs
            pr_obs(jc,jb)     = radar_data%radar_td%obs(jc,jb,weight_index_m1lim)                               &
              + (radar_data%radar_td%obs(jc,jb,weight_index_p1)-radar_data%radar_td%obs(jc,jb,weight_index_m1lim))  &
              * (abs(td_in_min(weight_index_m1lim)))/ (2.0_wp*assimilation_config(jg)%lhn_dt_obs)
            pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
            wobs_time(jc,jb) = 0.75_wp   
            num_t_obs (jc,jb,2) = 1
            IF (assimilation_config(jg)%lhn_spqual) THEN
              wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_m1lim)                               &
                + (radar_data%radar_td%spqual(jc,jb,weight_index_p1)-radar_data%radar_td%spqual(jc,jb,weight_index_m1lim))&
                * (abs(td_in_min(weight_index_m1lim)))/ (2.0_wp*assimilation_config(jg)%lhn_dt_obs)
            ELSE
              wobs_space(jc,jb) = 1.0_wp
            ENDIF
          ELSEIF ((radar_data%radar_td%obs(jc,jb,weight_index_0) >= pr_time_limit) .AND. &
                  (lp2 .AND. radar_data%radar_td%obs(jc,jb,weight_index_p2)>= pr_time_limit)) THEN
            ! observation is valid between t>=0 and t<=+2assimilation_config(jg)%lhn_dt_obs
            pr_obs(jc,jb) = radar_data%radar_td%obs(jc,jb,weight_index_0)                                   &
              + (radar_data%radar_td%obs(jc,jb,weight_index_p2)-radar_data%radar_td%obs(jc,jb,weight_index_0))      &
              * (abs(td_in_min(weight_index_0)))/ (2.0_wp*assimilation_config(jg)%lhn_dt_obs)
            pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
            wobs_time(jc,jb) = 0.75_wp
            num_t_obs (jc,jb,2)= 1
            IF (assimilation_config(jg)%lhn_spqual) THEN
              wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_0)                                   &
                + (radar_data%radar_td%spqual(jc,jb,weight_index_p2)-radar_data%radar_td%spqual(jc,jb,weight_index_0)) &
                * (abs(td_in_min(weight_index_0)))/ (2.0_wp*assimilation_config(jg)%lhn_dt_obs)
            ELSE
              wobs_space(jc,jb) = 1.0_wp
            ENDIF
          ELSEIF((lm1 .AND. radar_data%radar_td%obs(jc,jb,weight_index_m1lim) >= pr_time_limit) .AND. &
                 (lp2 .AND. radar_data%radar_td%obs(jc,jb,weight_index_p2) >= pr_time_limit)) THEN
            ! observation is valid between t>=-lhn_dt_obs and t<=+2lhn_dt_obs
            pr_obs(jc,jb)    = radar_data%radar_td%obs(jc,jb,weight_index_m1lim)                                  &
              + (radar_data%radar_td%obs(jc,jb,weight_index_p2)-radar_data%radar_td%obs(jc,jb,weight_index_m1lim))     &
              * (abs(td_in_min(weight_index_m1lim)))/ (3.0_wp*assimilation_config(jg)%lhn_dt_obs)
            pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
            wobs_time(jc,jb) = 0.5_wp
            num_t_obs (jc,jb,3) = 1
            IF (assimilation_config(jg)%lhn_spqual) THEN
              wobs_space(jc,jb)= radar_data%radar_td%spqual(jc,jb,weight_index_m1lim)                                  &
                + (radar_data%radar_td%spqual(jc,jb,weight_index_p2)-radar_data%radar_td%spqual(jc,jb,weight_index_m1lim)) &
                * (abs(td_in_min(weight_index_m1lim)))/ (3.0_wp*assimilation_config(jg)%lhn_dt_obs)
            ELSE
              wobs_space(jc,jb) = 1.0_wp
            ENDIF
          ELSEIF((lm2 .AND. radar_data%radar_td%obs(jc,jb,weight_index_m2lim) >= pr_time_limit) .AND. &
                 (lp1 .AND. radar_data%radar_td%obs(jc,jb,weight_index_p1) >= pr_time_limit)) THEN
            ! observation is valid between t>=-2lhn_dt_obs and t<=+lhn_dt_obs
            pr_obs(jc,jb) = radar_data%radar_td%obs(jc,jb,weight_index_m2lim)                                  &
              + (radar_data%radar_td%obs(jc,jb,weight_index_p1)-radar_data%radar_td%obs(jc,jb,weight_index_m2lim))     &
              * (abs(td_in_min(weight_index_m2lim)))/ (3.0_wp*assimilation_config(jg)%lhn_dt_obs)
            pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
            wobs_time(jc,jb) = 0.5_wp
            num_t_obs (jc,jb,3) = 1
            IF (assimilation_config(jg)%lhn_spqual) THEN
              wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_m2lim)                                  &
                + (radar_data%radar_td%spqual(jc,jb,weight_index_p1)-radar_data%radar_td%spqual(jc,jb,weight_index_m2lim))&
                * (abs(td_in_min(weight_index_m2lim)))/(3.0_wp*assimilation_config(jg)%lhn_dt_obs)
            ELSE
              wobs_space(jc,jb) = 1.0_wp
            ENDIF
          ELSEIF((radar_data%radar_td%obs(jc,jb,weight_index_0) >= pr_time_limit) .AND. &
                 (lp3 .AND. radar_data%radar_td%obs(jc,jb,weight_index_p3)>= pr_time_limit)) THEN
            ! observation is valid between t>=0 and t<=+3lhn_dt_obs
            pr_obs(jc,jb) = radar_data%radar_td%obs(jc,jb,weight_index_0)                                  &
              + (radar_data%radar_td%obs(jc,jb,weight_index_p3)-radar_data%radar_td%obs(jc,jb,weight_index_0))     &
              * (abs(td_in_min(weight_index_0)))/ (3.0_wp*assimilation_config(jg)%lhn_dt_obs)
            pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
            wobs_time(jc,jb) = 0.5_wp
            num_t_obs (jc,jb,3) = 1
            IF (assimilation_config(jg)%lhn_spqual) THEN
              wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_0)                                  &
                + (radar_data%radar_td%spqual(jc,jb,weight_index_p3)-radar_data%radar_td%spqual(jc,jb,weight_index_0)) &
                * (abs(td_in_min(weight_index_0)))/ (3.0_wp*assimilation_config(jg)%lhn_dt_obs)
            ELSE
              wobs_space(jc,jb) = 1.0_wp
            ENDIF
         ELSEIF((lm2 .AND. radar_data%radar_td%obs(jc,jb,weight_index_m2lim) >= pr_time_limit) .AND. &
                (lp2 .AND. radar_data%radar_td%obs(jc,jb,weight_index_p2) >= pr_time_limit)) THEN
           ! observation is valid between t>=-2lhn_dt_obs and t<=+2lhn_dt_obs
           pr_obs(jc,jb) = radar_data%radar_td%obs(jc,jb,weight_index_m2lim)                                  &
             + (radar_data%radar_td%obs(jc,jb,weight_index_p2)-radar_data%radar_td%obs(jc,jb,weight_index_m2lim))     &
             * (abs(td_in_min(weight_index_m2lim)))/ (4.0_wp*assimilation_config(jg)%lhn_dt_obs)
           pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
           wobs_time(jc,jb) = 0.25_wp
           num_t_obs (jc,jb,4) = 1
           IF (assimilation_config(jg)%lhn_spqual) THEN
             wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_m2lim)                                  &
               + (radar_data%radar_td%spqual(jc,jb,weight_index_p2)-radar_data%radar_td%spqual(jc,jb,weight_index_m2lim)) &
               * (abs(td_in_min(weight_index_m2lim)))/ (4.0_wp*assimilation_config(jg)%lhn_dt_obs)
           ELSE
             wobs_space(jc,jb) = 1.0_wp
           ENDIF
        ELSEIF((lm1 .AND. radar_data%radar_td%obs(jc,jb,weight_index_m1lim) >= pr_time_limit) .AND. &
               (lp3 .AND. radar_data%radar_td%obs(jc,jb,weight_index_p3) >= pr_time_limit)) THEN
           ! observation is valid between t>=-lhn_dt_obs and t<=+3lhn_dt_obs
           pr_obs(jc,jb) = radar_data%radar_td%obs(jc,jb,weight_index_m1lim)                                  &
             + (radar_data%radar_td%obs(jc,jb,weight_index_p3)-radar_data%radar_td%obs(jc,jb,weight_index_m1lim))     &
             * (abs(td_in_min(weight_index_m1lim)))/ (4.0_wp*assimilation_config(jg)%lhn_dt_obs)
           pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
           wobs_time(jc,jb) = 0.25_wp
           num_t_obs (jc,jb,4) = 1
           IF (assimilation_config(jg)%lhn_spqual) THEN
             wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_m1lim)                                  &
               + (radar_data%radar_td%spqual(jc,jb,weight_index_p3)-radar_data%radar_td%spqual(jc,jb,weight_index_m1lim)) &
               * (abs(td_in_min(weight_index_m1lim)))/ (4.0_wp*assimilation_config(jg)%lhn_dt_obs)
           ELSE
              wobs_space(jc,jb) = 1.0_wp
           ENDIF
         ELSE
           ! observation is not valid
           pr_obs(jc,jb) = -0.1_wp
           wobs_space(jc,jb) = 0.0_wp
           wobs_time(jc,jb) = 0.0_wp
           num_t_obs (jc,jb,0) = 1
         ENDIF
        ENDIF
       ENDDO
     ENDDO
  
!$OMP END DO 
!$OMP END PARALLEL

   ELSE  ! td_min < 0 !!!

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx) ICON_OMP_GUIDED_SCHEDULE
     DO jb=i_startblk,i_endblk
       CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
       &                i_startidx, i_endidx, i_rlstart, i_rlend)
       DO jc=i_startidx,i_endidx

        IF (NINT(radar_data%radar_ct%blacklist(jc,jb)) /= 1_i4 .AND. NINT(lhn_fields%brightband(jc,jb)) /= 1_i4) THEN
         IF ((radar_data%radar_td%obs(jc,jb,weight_index_0) >= pr_time_limit) .AND. &
             (lm1 .AND. radar_data%radar_td%obs(jc,jb,weight_index_m1lim) >= pr_time_limit)) THEN
           ! observation is valid between t>=0 and t<=+lhn_dt_obs
           pr_obs(jc,jb)    = radar_data%radar_td%obs(jc,jb,weight_index_m1lim)                           &
             + (radar_data%radar_td%obs(jc,jb,weight_index_0)-radar_data%radar_td%obs(jc,jb,weight_index_m1lim)) &
             * (ABS(td_in_min(weight_index_m1lim)))/ assimilation_config(jg)%lhn_dt_obs
           pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
           wobs_time(jc,jb) = 1.0_wp
           num_t_obs (jc,jb,1) = 1
           IF (assimilation_config(jg)%lhn_spqual) THEN
             wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_m1lim)                                   &
               + (radar_data%radar_td%spqual(jc,jb,weight_index_0)-radar_data%radar_td%spqual(jc,jb,weight_index_m1lim)) &
               * (ABS(td_in_min(weight_index_m1lim)))/ assimilation_config(jg)%lhn_dt_obs
           ELSE
             wobs_space(jc,jb) = 1.0_wp
           ENDIF
         ELSEIF ((radar_data%radar_td%obs(jc,jb,weight_index_0) >= pr_time_limit) .AND. &
                 (lp1 .AND. radar_data%radar_td%obs(jc,jb,weight_index_p1) >= pr_time_limit)) THEN
           ! observation is valid between t>=0 and t<=+lhn_dt_obs
           pr_obs(jc,jb)    = radar_data%radar_td%obs(jc,jb,weight_index_0)                           &
             - (radar_data%radar_td%obs(jc,jb,weight_index_p1)-radar_data%radar_td%obs(jc,jb,weight_index_0)) &
             * (ABS(td_in_min(weight_index_0)))/ assimilation_config(jg)%lhn_dt_obs
           pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
           wobs_time(jc,jb) = 1.0_wp
           num_t_obs (jc,jb,1) = 1
           IF (assimilation_config(jg)%lhn_spqual) THEN
             wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_0)                                   &
               - (radar_data%radar_td%spqual(jc,jb,weight_index_p1)-radar_data%radar_td%spqual(jc,jb,weight_index_0)) &
               * (ABS(td_in_min(weight_index_0)))/ assimilation_config(jg)%lhn_dt_obs
           ELSE
             wobs_space(jc,jb) = 1.0_wp
           ENDIF
         ELSEIF((lm1 .AND. radar_data%radar_td%obs(jc,jb,weight_index_m1lim) >= pr_time_limit) .AND. &
                (lp1 .AND. radar_data%radar_td%obs(jc,jb,weight_index_p1) >= pr_time_limit)) THEN
           ! observation is valid between t>=-lhn_dt_obs and t<=+lhn_dt_obs
           pr_obs(jc,jb)     = radar_data%radar_td%obs(jc,jb,weight_index_m1lim)                               &
             + (radar_data%radar_td%obs(jc,jb,weight_index_p1)-radar_data%radar_td%obs(jc,jb,weight_index_m1lim))  &
             * (ABS(td_in_min(weight_index_m1lim)))/ (2.0_wp*assimilation_config(jg)%lhn_dt_obs)
           pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
           wobs_time(jc,jb) = 0.75_wp   
           num_t_obs (jc,jb,2) = 1
           IF (assimilation_config(jg)%lhn_spqual) THEN
             wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_m1lim)                               &
               + (radar_data%radar_td%spqual(jc,jb,weight_index_p1)-radar_data%radar_td%spqual(jc,jb,weight_index_m1lim))&
               * (ABS(td_in_min(weight_index_m1lim)))/ (2.0_wp*assimilation_config(jg)%lhn_dt_obs)
           ELSE
             wobs_space(jc,jb) = 1.0_wp
           ENDIF
         ELSEIF((radar_data%radar_td%obs(jc,jb,weight_index_0) >= pr_time_limit) .AND. &
                (lm2 .AND. radar_data%radar_td%obs(jc,jb,weight_index_m2lim)>= pr_time_limit)) THEN
           ! observation is valid between t>=0 and t<=+2assimilation_config(jg)%lhn_dt_obs
           pr_obs(jc,jb) = radar_data%radar_td%obs(jc,jb,weight_index_m2lim)                                   &
             + (radar_data%radar_td%obs(jc,jb,weight_index_0)-radar_data%radar_td%obs(jc,jb,weight_index_m2lim))      &
             * (ABS(td_in_min(weight_index_m2lim)))/ (2.0_wp*assimilation_config(jg)%lhn_dt_obs)
           pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
           wobs_time(jc,jb) = 0.75_wp
           num_t_obs (jc,jb,2)= 1
           IF (assimilation_config(jg)%lhn_spqual) THEN
             wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_m2lim)                                   &
               + (radar_data%radar_td%spqual(jc,jb,weight_index_0)-radar_data%radar_td%spqual(jc,jb,weight_index_m2lim)) &
               * (ABS(td_in_min(weight_index_m2lim)))/ (2.0_wp*assimilation_config(jg)%lhn_dt_obs)
           ELSE
             wobs_space(jc,jb) = 1.0_wp
           ENDIF
         ELSEIF((radar_data%radar_td%obs(jc,jb,weight_index_0) >= pr_time_limit) .AND. &
                (lp2 .AND.radar_data%radar_td%obs(jc,jb,weight_index_p2)>= pr_time_limit)) THEN
           ! observation is valid between t>=0 and t<=+2assimilation_config(jg)%lhn_dt_obs
           pr_obs(jc,jb) = radar_data%radar_td%obs(jc,jb,weight_index_0)                                   &
             - (radar_data%radar_td%obs(jc,jb,weight_index_p2)-radar_data%radar_td%obs(jc,jb,weight_index_0))      &
             * (ABS(td_in_min(weight_index_0)))/ (2.0_wp*assimilation_config(jg)%lhn_dt_obs)
           pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
           wobs_time(jc,jb) = 0.75_wp
           num_t_obs (jc,jb,2)= 1
           IF (assimilation_config(jg)%lhn_spqual) THEN
             wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_0)                                   &
               - (radar_data%radar_td%spqual(jc,jb,weight_index_p2)-radar_data%radar_td%spqual(jc,jb,weight_index_0)) &
               * (ABS(td_in_min(weight_index_0)))/ (2.0_wp*assimilation_config(jg)%lhn_dt_obs)
           ELSE
             wobs_space(jc,jb) = 1.0_wp
           ENDIF
         ELSEIF((lm1.AND.radar_data%radar_td%obs(jc,jb,weight_index_m1lim) >= pr_time_limit) .AND. &
                (lp2.AND.radar_data%radar_td%obs(jc,jb,weight_index_p2) >= pr_time_limit)) THEN
           ! observation is valid between t>=-lhn_dt_obs and t<=+2lhn_dt_obs
           pr_obs(jc,jb)    = radar_data%radar_td%obs(jc,jb,weight_index_m1lim)                                  &
             + (radar_data%radar_td%obs(jc,jb,weight_index_p2)-radar_data%radar_td%obs(jc,jb,weight_index_m1lim))     &
             * (ABS(td_in_min(weight_index_m1lim)))/ (3.0_wp*assimilation_config(jg)%lhn_dt_obs)
           pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
           wobs_time(jc,jb) = 0.5_wp
           num_t_obs (jc,jb,3) = 1
           IF (assimilation_config(jg)%lhn_spqual) THEN
             wobs_space(jc,jb)= radar_data%radar_td%spqual(jc,jb,weight_index_m1lim)                                  &
               + (radar_data%radar_td%spqual(jc,jb,weight_index_p2)-radar_data%radar_td%spqual(jc,jb,weight_index_m1lim)) &
               * (ABS(td_in_min(weight_index_m1lim)))/ (3.0_wp*assimilation_config(jg)%lhn_dt_obs)
           ELSE
             wobs_space(jc,jb) = 1.0_wp
           ENDIF
         ELSEIF((lm2 .AND. radar_data%radar_td%obs(jc,jb,weight_index_m2lim) >= pr_time_limit) .AND. &
                (lp1 .AND. radar_data%radar_td%obs(jc,jb,weight_index_p1) >= pr_time_limit)) THEN
           ! observation is valid between t>=-2lhn_dt_obs and t<=+lhn_dt_obs
           pr_obs(jc,jb) = radar_data%radar_td%obs(jc,jb,weight_index_m2lim)                                  &
             + (radar_data%radar_td%obs(jc,jb,weight_index_p1)-radar_data%radar_td%obs(jc,jb,weight_index_m2lim))     &
             * (ABS(td_in_min(weight_index_m2lim)))/ (3.0_wp*assimilation_config(jg)%lhn_dt_obs)
           pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
           wobs_time(jc,jb) = 0.5_wp
           num_t_obs (jc,jb,3) = 1
           IF (assimilation_config(jg)%lhn_spqual) THEN
             wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_m2lim)                                  &
             + (radar_data%radar_td%spqual(jc,jb,weight_index_p1)-radar_data%radar_td%spqual(jc,jb,weight_index_m2lim))&
             * (ABS(td_in_min(weight_index_m2lim)))/(3.0_wp*assimilation_config(jg)%lhn_dt_obs)
           ELSE
             wobs_space(jc,jb) = 1.0_wp
           ENDIF
         ELSEIF((radar_data%radar_td%obs(jc,jb,weight_index_0) >= pr_time_limit) .AND. &
                (lp3 .AND. radar_data%radar_td%obs(jc,jb,weight_index_p3)>= pr_time_limit)) THEN
           ! observation is valid between t>=0 and t<=+3lhn_dt_obs
           pr_obs(jc,jb) = radar_data%radar_td%obs(jc,jb,weight_index_0)                                  &
             - (radar_data%radar_td%obs(jc,jb,weight_index_p3)-radar_data%radar_td%obs(jc,jb,weight_index_0))     &
             * (ABS(td_in_min(weight_index_0)))/ (3.0_wp*assimilation_config(jg)%lhn_dt_obs)
           pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
           wobs_time(jc,jb) = 0.5_wp
           num_t_obs (jc,jb,3) = 1
           IF (assimilation_config(jg)%lhn_spqual) THEN
             wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_0)                                  &
               - (radar_data%radar_td%spqual(jc,jb,weight_index_p3)-radar_data%radar_td%spqual(jc,jb,weight_index_0)) &
               * (ABS(td_in_min(weight_index_0)))/ (3.0_wp*assimilation_config(jg)%lhn_dt_obs)
           ELSE
             wobs_space(jc,jb) = 1.0_wp
           ENDIF
         ELSEIF((lm1 .AND. radar_data%radar_td%obs(jc,jb,weight_index_m1lim) >= pr_time_limit) .AND. &
                (lp3 .AND. radar_data%radar_td%obs(jc,jb,weight_index_p3) >= pr_time_limit)) THEN
           ! observation is valid between t>=-lhn_dt_obs and t<=+3lhn_dt_obs
           pr_obs(jc,jb) = radar_data%radar_td%obs(jc,jb,weight_index_m1lim)                                  &
             + (radar_data%radar_td%obs(jc,jb,weight_index_p3)-radar_data%radar_td%obs(jc,jb,weight_index_m1lim))     &
             * (ABS(td_in_min(weight_index_m1lim)))/ (4.0_wp*assimilation_config(jg)%lhn_dt_obs)
           pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
           wobs_time(jc,jb) = 0.25_wp
           num_t_obs (jc,jb,4) = 1
           IF (assimilation_config(jg)%lhn_spqual) THEN
             wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_m1lim)                                  &
               + (radar_data%radar_td%spqual(jc,jb,weight_index_p3)-radar_data%radar_td%spqual(jc,jb,weight_index_m1lim)) &
               * (ABS(td_in_min(weight_index_m1lim)))/ (4.0_wp*assimilation_config(jg)%lhn_dt_obs)
           ELSE
             wobs_space(jc,jb) = 1.0_wp
           ENDIF
         ELSEIF((lm2 .AND. radar_data%radar_td%obs(jc,jb,weight_index_m2lim) >= pr_time_limit) .AND. &
                (lp2 .AND. radar_data%radar_td%obs(jc,jb,weight_index_p2) >= pr_time_limit)) THEN
           ! observation is valid between t>=-2lhn_dt_obs and t<=+2lhn_dt_obs
           pr_obs(jc,jb) = radar_data%radar_td%obs(jc,jb,weight_index_m2lim)                                  &
             + (radar_data%radar_td%obs(jc,jb,weight_index_p2)-radar_data%radar_td%obs(jc,jb,weight_index_m2lim))     &
             * (ABS(td_in_min(weight_index_m2lim)))/ (4.0_wp*assimilation_config(jg)%lhn_dt_obs)
           pr_obs(jc,jb) = pr_obs(jc,jb)*sec_per_hr_inv
           wobs_time(jc,jb) = 0.25_wp
           num_t_obs (jc,jb,4) = 1
           IF (assimilation_config(jg)%lhn_spqual) THEN
             wobs_space(jc,jb) = radar_data%radar_td%spqual(jc,jb,weight_index_m2lim)                                  &
               + (radar_data%radar_td%spqual(jc,jb,weight_index_p2)-radar_data%radar_td%spqual(jc,jb,weight_index_m2lim)) &
               * (ABS(td_in_min(weight_index_m2lim)))/ (4.0_wp*assimilation_config(jg)%lhn_dt_obs)
           ELSE
             wobs_space(jc,jb) = 1.0_wp
           ENDIF
         ELSE
           ! observation is not valid
           pr_obs(jc,jb) = -0.1_wp
           wobs_space(jc,jb) = 0.0_wp
           wobs_time(jc,jb) = 0.0_wp
           num_t_obs (jc,jb,0) = 1
         ENDIF
        ENDIF
       ENDDO
     ENDDO
!$OMP END DO 
!$OMP END PARALLEL
  
   ENDIF

  DEALLOCATE (td_in_min)

! if no spatial quality function is used, set wobs_space constant to one
!IF (.NOT.assimilation_config(jg)%lhn_spqual) wobs_space(:,:) = 1.0_wp

! determine statistics about spatial weights

  IF ( assimilation_config(jg)%lhn_diag ) THEN
    diag_out(:) = 0
    diag_out(1) = count( num_t_obs(:,:,1) == 1) !num1delta_t_obs
    diag_out(2) = count( num_t_obs(:,:,2) == 1) !num2delta_t_obs
    diag_out(3) = count( num_t_obs(:,:,3) == 1) !num3delta_t_obs
    diag_out(4) = count( num_t_obs(:,:,4) == 1) !num4delta_t_obs
    diag_out(5) = count( num_t_obs(:,:,0) == 1)         !numnone
    diag_out(6) = count( radar_data%radar_ct%blacklist(:,:) > 0.0_wp ) !blacklist
    diag_out(7) = count( wobs_space(:,:) == 1.0_wp)    !numfull
    diag_out(8) = count( (wobs_space(:,:) <  1.0_wp) &  !numred
                         .and. (wobs_space(:,:) >  0.0_wp) )
    diag_out(9) = count( wobs_space(:,:) == 0.0_wp)    !numzero

    g_diag_sum = 0
    g_diag_sum(1:ndiag) = global_sum( diag_out(1:ndiag),opt_iroot=p_io )

    IF (my_process_is_stdio()) THEN
     WRITE(nulhn, *)
     WRITE(nulhn, *)' Diagnostics of RADAR obs time interpolation, lhn_dt_obs = ',assimilation_config(jg)%lhn_dt_obs
     WRITE(nulhn, *)' number of treated obs points in time weighting : ',  &
          g_diag_sum( 1)+g_diag_sum( 2)+g_diag_sum( 3)+g_diag_sum( 4)+g_diag_sum( 5)
     WRITE(nulhn, *)' n of points with obs-dist 1 delta_t  : num1delta_t_obs = ',g_diag_sum( 1)
     WRITE(nulhn, *)' n of points with obs-dist 2 delta_t  : num2delta_t_obs = ',g_diag_sum( 2)
     WRITE(nulhn, *)' n of points with obs-dist 3 delta_t  : num3delta_t_obs = ',g_diag_sum( 3)
     WRITE(nulhn, *)' n of points with obs-dist 4 delta_t  : num4delta_t_obs = ',g_diag_sum( 4)
     WRITE(nulhn, *)' n of points without obs              : numnone = ',g_diag_sum( 5)
     WRITE(nulhn, *)' n of points which are blacklisted    : numblack = ',g_diag_sum( 6)
     WRITE(nulhn, *)
     WRITE(nulhn, *)' Diagnostics of RADAR obs space weighting'
     WRITE(nulhn, *)' n of points with full    spatial weight : numfull = ',g_diag_sum( 7)
     WRITE(nulhn, *)' n of points with reduced spatial weight : numred  = ',g_diag_sum( 8)
     WRITE(nulhn, *)' n of points with zero    spatial weight : numzero = ',g_diag_sum( 9)
     WRITE(nulhn, *)
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Section 5: Deallocate space of temporal arrays
!-------------------------------------------------------------------------------

  CALL deallocateTimedelta(inc_time_p)
  CALL deallocateTimedelta(inc_time_m)
  CALL deallocateTimedelta(titer)
  CALL deallocateTimedelta(time_delta)

  CALL deallocateDatetime(tnow)
  CALL deallocateDatetime(center_time)
  CALL deallocateDatetime(next_time_1)
  CALL deallocateDatetime(next_time_2)
  CALL deallocateDatetime(next_time_3)
  CALL deallocateDatetime(prev_time_1)
  CALL deallocateDatetime(prev_time_2)

!-------------------------------------------------------------------------------
! End of subroutine 
!-------------------------------------------------------------------------------

END SUBROUTINE lhn_obs_prep

SUBROUTINE detect_bright_band(pt_patch,radar_data,prm_diag,lhn_fields,sumrad)
!-------------------------------------------------------------------------------
!
! Description:
!   This procedure in the module "lheat_nudge" is called by "lhn_obs_prep" and
!   detects all grid points which are possibly influenced by bright band effects
!  
!-------------------------------------------------------------------------------

  TYPE(t_patch),   TARGET, INTENT(in)    :: pt_patch     !<grid/patch info.
  TYPE(t_radar_fields),    INTENT(in)    :: radar_data
  TYPE(t_nwp_phy_diag),    INTENT(in)    :: prm_diag
  TYPE(t_lhn_diag),        INTENT(inout)    :: lhn_fields

  REAL (KIND=wp), DIMENSION(nproma,pt_patch%nblks_c),INTENT(IN)    :: &
     sumrad

  INTEGER (KIND=i4) :: &
    jc,jb,nh,jg,nbright,nbrightg

  INTEGER :: i_rlstart, i_rlend
  INTEGER :: i_startblk, i_endblk    !> blocks
  INTEGER :: i_startidx, i_endidx    !< slices

!- End of header
!===============================================================================
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
  yroutine = 'detect_bright_band'

     !-------------------------------------------------------------------------------
     ! Evaluate radar height with resprect to temperature
     !-------------------------------------------------------------------------------

   ! exclude boundary interpolation zone of nested domains
   jg = pt_patch%id
   i_rlstart = grf_bdywidth_c+1
   i_rlend   = min_rlcell_int

   i_startblk = pt_patch%cells%start_block(i_rlstart)
   i_endblk   = pt_patch%cells%end_block(i_rlend)

   nbright    = 0
   nbrightg    = 0


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,nh) ICON_OMP_GUIDED_SCHEDULE
   DO jb=i_startblk,i_endblk
       CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
         &                i_startidx, i_endidx, i_rlstart, i_rlend)
       DO jc=i_startidx,i_endidx
           IF(NINT(radar_data%radar_ct%blacklist(jc,jb)) /= 1_i4 .AND. &
             MAXVAL(radar_data%radar_td%radheight(jc,jb,:)) > 0.0_wp ) THEN
              DO nh=1,assimilation_config(jg)%nradar
                 IF  ( radar_data%radar_td%radheight(jc,jb,nh) > 0.0_wp                   &
               .AND.  (prm_diag%hzerocl(jc,jb)-radar_data%radar_td%radheight(jc,jb,nh)) >= -100_wp &
               .AND.  (prm_diag%hzerocl(jc,jb)-radar_data%radar_td%radheight(jc,jb,nh)) <= 1000_wp &
               .AND.  sumrad(jc,jb) > 1.0_wp) THEN
!             ) THEN
                       lhn_fields%brightband(jc,jb)=1.0_wp
                       EXIT
                 ENDIF
              ENDDO
           ENDIF
       ENDDO
   ENDDO
!$OMP END DO 
!$OMP END PARALLEL
  IF ( assimilation_config(jg)%lhn_diag ) THEN
    nbright=COUNT(lhn_fields%brightband > 0.0)
    nbrightg=global_sum(nbright,opt_iroot=p_io)
    IF (my_process_is_stdio()) &
     WRITE(nulhn, *)' n of points which are possibly brightband    : numbright = ',nbrightg
  ENDIF
   
!-------------------------------------------------------------------------------
! End of subroutine
!-------------------------------------------------------------------------------

END SUBROUTINE detect_bright_band


!===============================================================================
!+ Module procedure in "lheat_nudge" determining T - increments due to LHN
!-------------------------------------------------------------------------------
 
SUBROUTINE lhn_t_inc (i_startidx, i_endidx,jg,ke,zlev,tt_lheat,wobs_time, wobs_space, &
                      pr_obs, pr_mod,pr_ana,ttend_lhn,treat_diag,scale_diag, &
                      scale_fac_index,u,v,k850,k950,k700,diag_out)

!-------------------------------------------------------------------------------
!
! Description:
!
!   Namelist parameters used: assimilation_config(jg)%lhn_artif,
!                             assimilation_config(jg)%fac_lhn_up,assimilation_config(jg)%fac_lhn_down
!                             assimilation_config(jg)%lhn_relax,nassimilation_config(jg)%lhn_relax
!                             assimilation_config(jg)%lhn_limit,assimilation_config(jg)%abs_lhn_lim
!                             assimilation_config(jg)%lhn_filt,assimilation_config(jg)%lhn_diag
!                             assimilation_config(jg)%lhn_incloud
!                             
!   Input arrays : tt_lheat, 
!                  pr_ana,pr_mod
!   Output arrays: ttend_lhn
!
! Method:
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Subroutine arguments :
!-------------------------------------------------------------------------------
! Scalar arguments, intent(in) :
!-------------------------------
  INTEGER   (KIND=i4), INTENT(IN)     ::       &
    i_startidx, i_endidx
  INTEGER   (KIND=i4), INTENT(IN)     ::       &  ! domain ID
    jg
  INTEGER   (KIND=i4), INTENT(IN)     ::       &
    ke             ! number of grid points to be treated by lhn

  REAL (KIND=wp), DIMENSION(:),INTENT(IN)    :: &
    wobs_time, wobs_space, pr_obs, pr_mod
  REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) ::   &   ! dim (ie,ke)
    tt_lheat,zlev
  REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT) ::   &   ! dim (ie,ke)
    ttend_lhn
  REAL(KIND=wp), DIMENSION(:), INTENT(OUT) ::   &   ! dim (ie,ke)
    treat_diag, pr_ana, scale_diag

  LOGICAL, INTENT(INOUT) :: &
    scale_fac_index(:)

  REAL(KIND=wp), DIMENSION(:,:), INTENT(IN) ::  u,v

  INTEGER (KIND=i4), DIMENSION(:),  INTENT(IN) :: k850, k950, k700

  INTEGER (KIND=i4), OPTIONAL, INTENT(OUT) :: &
    diag_out(:)

!-------------------------------------------------------------------------------
! Local parameters, scalars, arrays :
!-------------------------------------------------------------------------------
! Local scalars:
!---------------

  INTEGER   (KIND=i4)                 ::       &
    i,ip,k            ,& ! indeces for position in local subdomain
    n_local           ,& ! diagnostic : number of points with local profiles
    n_artif            ,& ! diagnostic : number of points with artificial profile
    n_up_lim,n_up     ,& ! diagnostic : number of points with limited upsclaing
    n_down_lim,n_down ,& ! diagnostic : number of points with limited downscaling
    n_ex_lim_p        ,& ! diagnostic : number of pts with inc. above limit
    n_ex_lim_n           ! diagnostic : number of pts with inc. below (neg) limit
    
  REAL (KIND=wp)                         ::       &
    epsilon=1.0E-35_wp ,& ! small number : add to avoid division by zero
    eps=0.2/3600.   ,& ! limit used for profile filterering (flags/
                              ! eliminates values below 0.2 deg heating/h)
    pr_quot         ,& ! ratio of analyzed (observed) to model precipitation
    pr_artif         ,& ! precip corresponding to artificial heating profile
    abs_lim_neg     ,& ! negative absolut limit for increments (= -assimilation_config(jg)%abs_lhn_lim)
    abs_lim_pos     ,& ! positive absolut limit for increments (= assimilation_config(jg)%abs_lhn_lim)
    prmax           ,&
    prmax_th        

! Local arrays:
!--------------

  LOGICAL                             ::       &
    lelim         ,& !  for elimination of isolated peaks in profile filter
    lsmooth          !  for smoothing in profile filter

  INTEGER   (KIND=i4)                              ::       &
    nelimosc, nelimiso, nsmooth ,& ! diagnostic variables for filter_prof
    nelimosc_proc, nelimiso_proc, nsmooth_proc, &
    n_windcor, n_windcor0, n_incloud, ntreat, treat_list(nproma)

  REAL (KIND=wp)                         ::       &
    scale_fac                 ,& ! scaling factor determined for each profile
    tt_artif(nproma,ke)       ,& ! artificial heating profile
    abs_lim_prof(nproma,ke)   ,& ! vertical profile for limitation the temperature increment
    w, ntcoeff, fac


  REAL (KIND=wp)    ::       &
    umean, vmean, zvb, zvb_llim, zvb_ulim, w950, w850, w700, wind_corr(nproma), &
    rfade

  REAL (KIND=wp)    ::       &
    tt_artif_max, zlev_artif_max, std_zlev, dz, tt_max_2_pr_artif

!- End of header
!===============================================================================
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
  yroutine='lhn_t_inc'

!-------------------------------------------------------------------------------
! Section 0 : Preliminaries : 
!             - determine (approx.) max number of profiles from neighbour nodes
!             Initialize fields
!-------------------------------------------------------------------------------

! calculate the height of the levels (middle of layers), location, 
! where lh-values are defined

!Values of climatological profile, used originally in COSMO
!  pr_artif = 0.0008333        ! climatological precip of 3 mm/h
!  tt_max_2_pr_artif = 0.00015        ! correspondend latent heat release (maximum value within column)

! corresponding values for ILAM, taken global mean values into account:
  pr_artif = 0.0003   ! 0.00021
  tt_max_2_pr_artif = 0.0005 ! 0.00037

! It might be possible to deirve the values as average over entire domain in the first call of LHN routine
!  pr_artif = global_sum_array (pr_mod)
!  do k=1,ke
!    tt_lheat_sum(k)=global_sum_array (tt_lheat(:,k))
!  enddo
!  tt_max_2_pr_artif = maxval(tt_lheat_sum(:))

  tt_artif_max = assimilation_config(jg)%tt_artif_max 
!20190524
  pr_artif = pr_artif * ABS(tt_artif_max)/tt_max_2_pr_artif 

!  pr_artif = tt_artif_max ! ILAM scheint in etwa ein Verhaltnis von 1:1 zu erzeugen. Sieht man im Mittelwert, wie auch in Extremfallen!

  zlev_artif_max = assimilation_config(jg)%zlev_artif_max
  std_zlev = zlev_artif_max/assimilation_config(jg)%std_artif_max


! set temperature increments to be determined to zero
  ttend_lhn = 0.0_wp

! flags for profile filter
  lelim=.TRUE.
  lsmooth=.TRUE.

! initializations of counters
  n_local    = 0
  n_artif     = 0
  n_up_lim   = 0
  n_up       = 0
  n_down_lim = 0
  n_down     = 0
  n_windcor  = 0
  n_windcor0 = 0
  nelimosc_proc = 0_i4
  nelimiso_proc = 0_i4
  nsmooth_proc  = 0_i4
  n_ex_lim_p = 0_i4
  n_ex_lim_n = 0_i4
  n_incloud  = 0_i4


! other initializations (assume that profiles are at local nproma,node first)
  scale_fac = 1._wp
  pr_quot = 1._wp

  abs_lim_pos = assimilation_config(jg)%abs_lhn_lim*assimilation_config(jg)%lhn_coef
  abs_lim_neg = -1. * abs_lim_pos

  IF (assimilation_config(jg)%lhn_wweight) THEN
    w950 = 0.5_wp
    w850 = 0.25_wp
    w700 = 0.25_wp
    zvb_llim = 20.0_wp
    zvb_ulim = 30.0_wp
  ENDIF


!-------------------------------------------------------------------------------
! Section 2 : Select the latent heating profile to scale and set scaling factor:
!             (only points with pr_ana != pr_mod are treated)
!             - determine if the local or a nearby heating profile is taken
!             - retain the indices of the selected profile and assign them to a
!               special list if it belongs to a subdomain held at another node
!             - set the scaling factor (depending on the ratio of analyzed to
!               model precipitation and respecting the scaling limits)
!-------------------------------------------------------------------------------

  ntreat = 0_i4
  
  rfade=assimilation_config(jg)%start_fadeout
  IF ( rnlhn >= rfade .AND. rfade < 1.0_wp ) then
     ntcoeff = ( 1.0_wp + rfade/( 1.0_wp - rfade ) ) * ( 1.0_wp - rnlhn )
     ntcoeff  = MAX(0.0_wp,ntcoeff)
  ELSE
     ntcoeff = 1.0_wp
  ENDIF

  DO ip=i_startidx,i_endidx

! Calculate pr_ana when wobs_time*wobs_spce > 0
     w = wobs_time(ip) * wobs_space(ip)
     IF ( w == 0.0_wp                                             &
       .OR. (pr_obs(ip)  < 0.0 )                                  &
       .OR. ((pr_obs(ip) < assimilation_config(jg)%thres_lhn)     &
       .AND. (pr_mod(ip) < assimilation_config(jg)%thres_lhn) )   &
        ) THEN
        pr_ana(ip) = pr_mod(ip)
        treat_diag(ip)=-1_wp
        CYCLE
     ELSE
        IF ( w < 0.0_wp .OR. w > 1.0_wp ) then
#ifndef __SX__
           WRITE(message_text,'(a,f8.4,i6)') 'lhn_pr_ana w unvalid : ',w,ip
           CALL message ('',message_text,level=2)
#endif
           CYCLE
        ELSE
           pr_ana(ip) = w * pr_obs(ip) + (1.0_wp-w) * pr_mod(ip)
           ntreat = ntreat +1
           treat_list(ntreat) = ip
           treat_diag(ip)=0_wp
        ENDIF
     ENDIF
  ENDDO

!NEC$ ivdep
  DO i = 1, ntreat
     ip = treat_list(i)
     pr_quot       = (pr_ana(ip)+epsilon)/(pr_mod(ip)+epsilon)

     IF(assimilation_config(jg)%lhn_logscale) THEN
       IF (pr_quot > 0.) THEN
          pr_quot = 1._wp + LOG(pr_quot)
       ELSE
!          PRINT *, 'pr_quot = 0 ', ip, pr_ana(ip), pr_mod(ip)
          pr_quot = 0.
       ENDIF
     ENDIF

! local model precip is within range [ assimilation_config(jg)%fac_lhn_down*pr_ana , assimilation_config(jg)%fac_lhn_up*pr_ana ]
! -> scaling of a local profile (upscaling or downscaling)
     IF ( pr_quot >= assimilation_config(jg)%fac_lhn_down .AND. pr_quot <= assimilation_config(jg)%fac_lhn_up ) THEN
        scale_fac = pr_quot
        n_local = n_local + 1
        treat_diag(ip)=1_wp
        IF ( pr_quot < 1 ) n_down = n_down + 1
        IF ( pr_quot > 1 ) n_up   = n_up   + 1

! local model precip is too large -> limited downscaling of local profile
     ELSEIF ( pr_quot < assimilation_config(jg)%fac_lhn_down ) THEN
        scale_fac = assimilation_config(jg)%fac_lhn_down 
        n_down_lim = n_down_lim + 1
        treat_diag(ip)=2_wp

     ELSEIF ( pr_quot > assimilation_config(jg)%fac_lhn_artif .AND. assimilation_config(jg)%lhn_artif) THEN

       IF (assimilation_config(jg)%lhn_logscale) THEN
          fac = exp(assimilation_config(jg)%fac_lhn_artif -1)
       ELSE
          fac = assimilation_config(jg)%fac_lhn_artif
       ENDIF
        prmax    = pr_mod(ip) * fac
        prmax_th = assimilation_config(jg)%thres_lhn * fac

        prmax_th = MIN(prmax_th,pr_ana(ip))
        prmax    = MAX(prmax,prmax_th)

        scale_fac = assimilation_config(jg)%fac_lhn_artif_tune * (pr_ana(ip))/pr_artif 

        n_artif = n_artif + 1
        treat_diag(ip)=4_wp

! local model precip is too small -> limited 
! upscaling of local profile
     ELSE
        scale_fac = assimilation_config(jg)%fac_lhn_up
        n_up_lim = n_up_lim + 1
        treat_diag(ip)=3_wp
     ENDIF
    
     IF (pr_quot > assimilation_config(jg)%fac_lhn_up) THEN
        scale_fac_index(ip)=.TRUE.
     ELSE
        scale_fac_index(ip)=.FALSE.
     ENDIF
     scale_diag(ip)=scale_fac

     IF (assimilation_config(jg)%lhn_artif_only) treat_diag(ip) = 4_wp

  ENDDO

!-------------------------------------------------------------------------------
! Section 4 : Scale the heating profiles with the predetermined factors
!             Insert the artificial profile if no suitable profile was found
!             -> Get the temperature increment due _only_ to lhn - correction
!             tt_lheat (= dT due to saturation adjustment has already been added to T
!             (src_leapfrog) so only the part due to LHN has to be added now:
!             ttend_lhn = tt_lheat * (scale_fac - 1)
!             Filter the profiles to exclude spurious peaks/noise
!-------------------------------------------------------------------------------

! filter and scale the local profiles and nearby profiles from this node
! or artificial profiles

  tt_artif(:,ke-4:ke) = 0.0
  abs_lim_prof(:,ke-4:ke) = 0.0
  DO k = 1, ke-5 ! do not touch the lowest model levels
!NEC$ ivdep
    DO i = 1, ntreat
      ip = treat_list(i)
      dz = zlev(ip,k) - zlev_artif_max
      tt_artif(ip,k) = tt_artif_max * exp(-0.5*((dz/std_zlev)**2))
      IF (tt_artif(ip,k) < 1.e-7) tt_artif(ip,k) = 0.0
      abs_lim_prof(ip,k) = abs_lim_pos * exp(-0.5*((dz/(std_zlev))**2))
      IF (abs_lim_prof(ip,k) < 1.e-7) abs_lim_prof(ip,k) = 0.0
    ENDDO
  ENDDO

!      determine the temperature correction due to lhn

  DO k = 1, ke
!NEC$ ivdep
    DO i = 1, ntreat
      ip = treat_list(i)
      IF (treat_diag(ip) > 0 .AND. treat_diag(ip) < 4) THEN

        ttend_lhn(ip,k) = (scale_diag(ip)-1.) * tt_lheat(ip,k)
        ttend_lhn(ip,k) = assimilation_config(jg)%lhn_coef * ttend_lhn(ip,k) * ntcoeff

      ELSE IF (treat_diag(ip) == 4) THEN

        ttend_lhn(ip,k) = ABS(scale_diag(ip) * tt_artif(ip,k)) * assimilation_config(jg)%lhn_coef * ntcoeff

        IF (ttend_lhn(ip,k) < 0._wp) ttend_lhn(ip,k) = 0.0_wp

      ENDIF
    ENDDO
  ENDDO
!-------------------------------------------------------------------------------
! Section 5 : Vertical filtering of increments (if assimilation_config(jg)%lhn_filt)
!-------------------------------------------------------------------------------

  IF (assimilation_config(jg)%lhn_filt) THEN

     CALL filter_prof (ttend_lhn,ntreat,treat_list,1,ke,      &
                       eps,lelim,lsmooth,nelimosc,nelimiso,nsmooth)
     nelimosc_proc = nelimosc_proc + nelimosc
     nelimiso_proc = nelimiso_proc + nelimiso
     nsmooth_proc  = nsmooth_proc  + nsmooth

  ENDIF


  DO k = 1, ke
!NEC$ ivdep
    DO i = 1, ntreat
      ip = treat_list(i)

!-------------------------------------------------------------------------------
! Section 6 : Set all increments to zero where tt_lheat in negative (i.e. layers without cloud)
!-------------------------------------------------------------------------------

      IF (assimilation_config(jg)%lhn_incloud .AND. .NOT. assimilation_config(jg)%lhn_artif_only) THEN
        IF (tt_lheat(ip,k) < 0.0_wp) THEN
          ttend_lhn(ip,k) = 0.0_wp
          n_incloud = n_incloud + 1
        ENDIF
      ENDIF

!-------------------------------------------------------------------------------
! Section 7 : Impose absolute limit on increments if requested (if assimilation_config(jg)%lhn_limit)
!-------------------------------------------------------------------------------

      IF (assimilation_config(jg)%lhn_limitp) THEN
        IF (ttend_lhn(ip,k) > abs_lim_prof(ip,k)) THEN
           ttend_lhn(ip,k) = abs_lim_prof(ip,k)
           n_ex_lim_p = n_ex_lim_p + 1
        ELSEIF (ttend_lhn(ip,k) < -1.*abs_lim_prof(ip,k)) THEN
           ttend_lhn(ip,k) = -1.*abs_lim_prof(ip,k)
           n_ex_lim_n = n_ex_lim_n + 1
        ENDIF

      ELSE IF (assimilation_config(jg)%lhn_limit) THEN

        IF (ttend_lhn(ip,k) > abs_lim_pos) THEN
           ttend_lhn(ip,k) = abs_lim_pos
           n_ex_lim_p = n_ex_lim_p + 1
        ELSEIF (ttend_lhn(ip,k) < abs_lim_neg) THEN
           ttend_lhn(ip,k) = abs_lim_neg
           n_ex_lim_n = n_ex_lim_n + 1
        ENDIF

      ENDIF
    ENDDO
  ENDDO
!-------------------------------------------------------------------------------
! Section 8 : Weighting of the temperature increment with respect to
!             the mean horizontal wind within the column
!             The mean horizontal wind is defined as linear combination of
!             0.5*wind(950 hPa)+0.25*(wind(850 hPa)+wind(700 hPa))
!
!-------------------------------------------------------------------------------

  IF (assimilation_config(jg)%lhn_wweight) THEN
!NEC$ ivdep
    DO i = 1, ntreat
      ip = treat_list(i)
      umean = w950 * u(ip,k950(ip))       &
            + w850 * u(ip,k850(ip))       &
            + w700 * u(ip,k700(ip))
      vmean = w950 * v(ip,k950(ip))       &
            + w850 * v(ip,k850(ip))       &
            + w700 * v(ip,k700(ip))
      zvb   = SQRT(umean * umean + vmean * vmean)

      IF (zvb <= zvb_llim) THEN
         wind_corr(ip)=1.0_wp
      ELSE IF (zvb <= zvb_ulim) THEN
         wind_corr(ip)=1.0_wp - (1.0_wp/(zvb_ulim-zvb_llim))   &
                                * (zvb - zvb_llim)
         n_windcor=n_windcor+1
      ELSE
         wind_corr(ip)=0.0_wp
         n_windcor0=n_windcor0+1
      ENDIF

    ENDDO

    DO k = 1, ke
!NEC$ ivdep
      DO i = 1, ntreat
        ip = treat_list(i)
        ttend_lhn(ip,k) = ttend_lhn(ip,k) * wind_corr(ip)
      ENDDO
    ENDDO

  ENDIF



!-------------------------------------------------------------------------------
! Section 9 : Diagnostic output on lhn - increments 
!             Summing information at PE0 for printout
!-------------------------------------------------------------------------------

   IF (PRESENT(diag_out)) THEN
     ! get summed diagnostics at PE0 from all PE's using collect_values 
     ! (collect_values needs real vector as input)
        diag_out = 0
        diag_out( 1) = ntreat
        diag_out( 2) = n_local
        diag_out( 3) = n_up + n_up_lim
        diag_out( 4) = n_up_lim
        diag_out( 5) = n_down + n_down_lim
        diag_out( 6) = n_down_lim
        diag_out( 7) = n_artif
        diag_out( 8) = n_ex_lim_p
        diag_out( 9) = n_ex_lim_n
        diag_out(10) = n_windcor
        diag_out(11) = n_windcor0
        diag_out(12)= nelimosc_proc
        diag_out(13)= nelimiso_proc
        diag_out(14)= nsmooth_proc
        diag_out(15)= n_incloud
        diag_out(16)= INT(assimilation_config(jg)%lhn_coef * ntcoeff * 100.)
   ENDIF
       
!-------------------------------------------------------------------------------
! End of subroutine 
!-------------------------------------------------------------------------------

END SUBROUTINE lhn_t_inc

!-------------------------------------------------------------------------------
 
!===============================================================================
!+ Module procedure in "lheat_nudge" adjusting humidity to LHN T - increments
!-------------------------------------------------------------------------------
 
SUBROUTINE lhn_q_inc(i_startidx,i_endidx,jg,zdt,ke,t,ttend_lhn,p,qv,qc,qi, &
                     qvtend_lhn,scale_fac_index,diag_out)

!-------------------------------------------------------------------------------
!
! Description:
!   Subroutine which computes q-increments for the previously added temperature
!   increments resulting from LHN.
!
! Method:
!   In areas of negative T-increments, qv is adjusted so that the relative
!   humidity from before is retained unaltered.
!   In areas of positive  T-increments, qv is raised, so that saturation is
!   reached. No adjustment to qc is made, since LHN assumes that the
!   condensed water vapour has precipitated.
!
!   relhum == const ---> relhum(told) = relhum(tnew)
!
!              rdv * (relhum(told)*esat(tnew))
!       qv  = ---------------------------------------
!              p - o_m_rdv * (relhum(told)*esat(tnew))
!
!
!   Used fields : ttend_lhn, t, qv, p0, pp
!   Modified fields : qv(.,.,.)
!
!-------------------------------------------------------------------------------
 
! Scalar arguments, intent(in) :
!-------------------------------
  INTEGER   (KIND=i4), INTENT(IN)     ::       &
     i_startidx,i_endidx,ke
  INTEGER   (KIND=i4), INTENT(IN)     ::       &  ! domain ID
     jg
  REAL (KIND=wp), INTENT(IN)  ::  &
     zdt,                             &
     t(:,:),                          &
     p(:,:),                          &
     qc(:,:),                         &
     qi(:,:),                         &
     qv(:,:)!,                         &
!     ttend_lhn(:,:)

  REAL (KIND=wp), INTENT(INOUT)  :: &
     qvtend_lhn(:,:),                         &
     ttend_lhn(:,:)

  LOGICAL, INTENT(IN)                :: &
     scale_fac_index(:)

  INTEGER (KIND=i4), OPTIONAL, INTENT(OUT) :: &
    diag_out(:)


! Local parameters, scalars, arrays :
!------------------------------------

  REAL (KIND=wp):: qv_new(SIZE(qv,1),SIZE(qv,2))


! Local parameters:
  REAL    (KIND=wp   ), SAVE ::  &
   epsy = 1.E-8_wp      ,& ! very small value > 0.
   delt_minn= -3.E-6_wp ,& ! minimal T-change before applying T-adjustment
   delt_minp= 3.E-6_wp  ,& ! minimal T-change before applying T-adjustment
   f_raise  = 1._wp     ,& ! relative humidity in positive adjustment areas
   tau_nudge = 1._wp/1800._wp     ,& ! time weight for nudging of the humidity
                               ! increment (increment spread over 30. min=1800 sec.)
   fac_q_max= 2._wp        ! maximal factor allowed in change of qv

! Local scalars:
  REAL    (KIND=wp   ) ::  &
   f_esat  ,& ! Name of satement function (saturation water vapour pressure)
   f_qv    ,& ! Name of satement function (specific humidity from p,e)
   f_e     ,& ! Name of satement function (water vapour pressure from q,p)
   zt      ,& ! Dummy argument for statement functions
   ze      ,& ! ...
   zqv     ,& ! ...
   zp      ,& ! ...
   esat    ,& ! saturation water vapour pressure
   relhum     ! relative humidity

  INTEGER (KIND=i4)  ::  &
   jc,k,        & ! loop indices
   nred,ninc   ,& ! number of points where humidity qv is reduced/increased
   ninc2          ! number of points where humidity qv is reduced/increased

!- End of header
!===============================================================================
!-------------------------------------------------------------------------------
! Begin of subroutine
!-------------------------------------------------------------------------------
! STATEMENT FUNCTIONS (identical to those used in meteo_utilities (satad) )
  f_esat(zt)     = b1*EXP( b2w*(zt-b3)/(zt-b4w) )      ! Magnus formula
  f_qv(ze,zp)    = rdv*ze/( zp - o_m_rdv*ze )          ! spec. hum. qv from e, p
  f_e (zqv,zp)   = MAX( epsy , zqv ) * zp / (rdv + zqv*o_m_rdv)  ! e from  qv, p

!-------------------------------------------------------------------------------
! 1. Increase / decrease specific humidity
!-------------------------------------------------------------------------------


  nred = 0
  ninc = 0
  ninc2 = 0

  zt = tau_nudge * zdt

  DO   k=1,ke
    DO jc=i_startidx,i_endidx

     IF ( ttend_lhn(jc,k) < delt_minn .OR. ttend_lhn(jc,k) > delt_minp) THEN
       zp = p(jc,k)

!       !saturation pressure before temperature increment
!       esat = f_esat ( t(jc,k) - ttend_lhn(jc,k) * zdt) ! for ICON temperature increment is still not applied yet
       esat = f_esat ( t(jc,k) )

!       ! relhum before temperature increment
       relhum = f_e ( qv(jc,k) , zp ) / esat
       relhum = MIN ( f_raise,relhum)

!       ! specific humidity after temperature increment so that relhum is unchanged
       qv_new(jc,k) = f_qv ( relhum * f_esat(t(jc,k)+ttend_lhn(jc,k) * zdt) , zp )

       IF ( ttend_lhn(jc,k) > delt_minp ) THEN
        
         ninc = ninc + 1
!ks: if criteria changed
         IF ( (scale_fac_index(jc)) .AND. (qc(jc,k)+qi(jc,k) <= epsy) ) THEN
!! add an additional increment to qv at gridpoints where the precipitation rate should
!! be increased and f has not reached 100% so far!
!! Attention: do not add these increments at points with a positive temperature increment
!! generally, because positive temperature increments can also occur at gridpoints where
!! the precipitation rate should be decreased (pos. temp. inc. below clouds, higher 
!! evaporation)

          ninc2= ninc2 + 1
          zqv = f_qv ( f_raise * f_esat (t(jc,k)+ttend_lhn(jc,k) * zdt), zp )
          zqv = MIN ( zqv , fac_q_max * qv_new(jc,k))
          zqv = zqv - qv_new(jc,k)
           !! All grid points with qc+qi>0 will not be adjusted!
          qv_new(jc,k) = qv_new(jc,k) + zqv*zt

         ENDIF
! Diagnostics Output

       ELSE IF ( ttend_lhn(jc,k) < delt_minn ) THEN
         nred = nred + 1
       ENDIF
     ELSE
       qv_new(jc,k) = qv(jc,k)
     ENDIF

     qvtend_lhn(jc,k) = (qv_new(jc,k) - qv(jc,k))/zdt

     IF ( assimilation_config(jg)%lhn_limit) qvtend_lhn(jc,k) = MIN ( qvtend_lhn(jc,k), 0.1*qv(jc,k)/zdt)

     IF ( assimilation_config(jg)%lhn_no_ttend ) ttend_lhn(jc,k) = 0.0 
     

    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! 2. Diagostic output
!-------------------------------------------------------------------------------
 IF (assimilation_config(jg)%lhn_diag) THEN
  IF (PRESENT(diag_out)) THEN
     diag_out    = 0
     diag_out( 1)= ninc
     diag_out( 2)= nred
     diag_out( 3)= ninc2
  ENDIF
 ENDIF


!-------------------------------------------------------------------------------
! End of subroutine 
!-------------------------------------------------------------------------------

END SUBROUTINE lhn_q_inc

!===============================================================================
!+ Saturation Adjustment after changing the temperature without adjusting the moisture
!-------------------------------------------------------------------------------

!===============================================================================
!+ Filtering of vertical (heating) profiles, elimination of isolated peaks
!-------------------------------------------------------------------------------

SUBROUTINE filter_prof (prof_filt,ntreat,treat_list,kup,klow,eps,lelim,lsmooth, &
                        nelimosc,nelimiso,nsmooth)

!-------------------------------------------------------------------------------
! Description:
!  This subroutines filters a vertical profile (e.g. heating profile to be used
!  in latent heat nudging to elimiate computational noise).
!
! Method: 
!  lelim : eliminate isolated peaks of small vertical extent
!    a) value on one level below and above is below specified eps
!    b) value two levels above is below eps and less than one level below the
!       profile value was idntified as very small or isolated peak
!
!  lsmooth : apply a simple one-dimensional shapiro-filter with S=1/2
!            to levels where the value is above eps
!
!-------------------------------------------------------------------------------
 
! Subroutine arguments, intent=in and intent=inout
  LOGICAL, INTENT(IN)       ::    &
    lelim             ,& ! flag 0/1 for elimination of isolated peaks
    lsmooth              ! flag 0/1 for smoothing

  INTEGER (KIND=i4), INTENT(IN)       ::    &
    ntreat, treat_list(:), & ! number of horizontal grid points and index list
    kup               ,& ! array dimensions of profile prof
    klow                 ! array dimensions of profile prof

  REAL (KIND=wp),    INTENT(INOUT)       ::      &
    prof_filt(nproma,kup:klow)    ,& ! on output: filtered profile array
    eps                  ! limits above which values are modified

! Local subroutine variables and arrays
  LOGICAL                                   ::    &
    lflag(nproma,kup:klow)          ! array to flag levels to be set to zero
  INTEGER   (KIND=i4)                ::    &
    i,ip,k,nheat(nproma), &           ! loop index, acceptable level counter
    nelimosc,nelimiso,nsmooth ! diag. output
  REAL (KIND=wp), DIMENSION(nproma, kup:klow) :: proffilt      ! profile array

!- End of header
!-------------------------------------------------------------------------------
! Begin Subroutine filter_prof
!-------------------------------------------------------------------------------

   nelimosc = 0
   nelimiso = 0
   nsmooth  = 0

! eliminate isolated peaks
   IF (lelim) THEN

!NEC$ ivdep
     DO i = 1, ntreat
       ip = treat_list(i)
!    eliminate oscillations around zero between adjacent levels
       IF ( ABS(prof_filt(ip,klow)) < eps ) THEN
         proffilt(ip,klow) = 0.
       ELSE
         proffilt(ip,klow) = prof_filt(ip,klow)
       ENDIF
       IF ( ABS(prof_filt(ip,kup)) < eps ) THEN
         proffilt(ip,kup) = 0.
       ELSE
         proffilt(ip,kup) = prof_filt(ip,kup)
       ENDIF
     ENDDO
     DO k=klow-1,kup+1,-1
!NEC$ ivdep
       DO i = 1, ntreat
         ip = treat_list(i)
         IF ( (prof_filt(ip,k-1)*prof_filt(ip,k)) <= 0. .AND.   &
           (prof_filt(ip,k)*prof_filt(ip,k+1)) <= 0.         ) THEN
            proffilt(ip,k) = 0._wp
            nelimosc = nelimosc + 1
         ELSE
            proffilt(ip,k) = prof_filt(ip,k)
         ENDIF
       ENDDO
     ENDDO


!    eliminate isolated peaks of small vertical extent
!    a) value on one level below and above is below specified eps
!    b) value two levels above is below eps and less than one level below the
!       profile value was idntified as very small or isolated peak
!    nheat : counter of lower levels with accepted heating rate (profile value)
!          - is reset to zero when an isolated peak is diagnosed or the heating
!            at the level considered is below specified eps
     nheat(:) = 0
     lflag(:,:) = .FALSE.

     DO k=klow-1,kup+2,-1
!NEC$ ivdep
       DO i = 1, ntreat
         ip = treat_list(i)
         nheat(ip)=nheat(ip)+1
         IF ( ABS(proffilt(ip,k-1))  <= eps .AND. ABS(proffilt(ip,k+1)) <= eps ) THEN
           lflag(ip,k) = .TRUE.
           nheat(ip) = 0
           nelimiso = nelimiso + 1
         ENDIF
         IF ( nheat(ip) < 2 .AND. ABS(proffilt(ip,k-2)) <= eps ) THEN
           lflag(ip,k-1) = .TRUE.
           lflag(ip,k) = .TRUE.
           nheat(ip) = 0
           nelimiso = nelimiso + 1
         ENDIF
       ENDDO
     ENDDO
     DO k=klow,kup,-1
!NEC$ ivdep
       DO i = 1, ntreat
         ip = treat_list(i)
         IF (lflag(ip,k)) proffilt(ip,k) = 0.
         prof_filt(ip,k) = proffilt(ip,k)
       ENDDO
     ENDDO
   ELSE
     proffilt(:,:)=prof_filt(:,:)
   ENDIF

! smooth profile
   IF (lsmooth) THEN

     DO k=klow-1,kup+1,-1
!NEC$ ivdep
       DO i = 1, ntreat
         ip = treat_list(i)
         IF ( ABS(proffilt(ip,k)) >= eps ) THEN
            prof_filt(ip,k) = 0.5  * proffilt(ip,k) + 0.25 * (proffilt(ip,k+1)+proffilt(ip,k-1))
            nsmooth = nsmooth + 1
         ENDIF
       ENDDO
     ENDDO
!NEC$ ivdep
     DO i = 1, ntreat
       ip = treat_list(i)
       IF ( ABS(proffilt(ip,klow)) >= eps) THEN
         prof_filt(ip,klow)=0.66 * proffilt(ip,klow) + 0.33 * proffilt(ip,klow-1)
         nsmooth = nsmooth + 1
       ENDIF
     ENDDO
!NEC$ ivdep
     DO i = 1, ntreat
       ip = treat_list(i)
       IF ( ABS(proffilt(ip,kup)) >= eps) THEN
         prof_filt(ip,kup)=0.66 * proffilt(ip,kup) + 0.33 * proffilt(ip,kup+1)
         nsmooth = nsmooth + 1
       ENDIF
     ENDDO
   ENDIF

!-------------------------------------------------------------------------------

END SUBROUTINE filter_prof

!===============================================================================

SUBROUTINE lhn_verification (ytime,pt_patch,radar_data,lhn_fields,nsteps,wobs_space,zprmod,zprmod_ref,zprrad,zprmod_ref_f,zprrad_f)

!-------------------------------------------------------------------------------
!
! Description:
! This subroutine calculates different parameter for verification of model precipitaion
! against radar observations. The values are stored in the LHN log file.
!-------------------------------------------------------------------------------

! Subroutine / Function arguments
!-------------------------------------------------------------------------------

 CHARACTER (LEN=*), INTENT(IN)  ::       &
   ytime

 TYPE(t_patch),   TARGET, INTENT(in)    :: pt_patch     !<grid/patch info.
 TYPE(t_radar_fields),    INTENT(in)    :: radar_data
 TYPE(t_lhn_diag),        INTENT(in)    :: lhn_fields

 REAL (KIND=wp), INTENT (IN) :: nsteps

 REAL (KIND=wp), INTENT(IN)  ::       &
   zprmod(nproma,pt_patch%nblks_c),                    &
   zprmod_ref(nproma,pt_patch%nblks_c),                &
   zprrad(nproma,pt_patch%nblks_c)

 REAL (KIND=wp), DIMENSION(nproma,pt_patch%nblks_c),INTENT(IN)    :: &
   wobs_space

 REAL (KIND=wp), INTENT(IN), OPTIONAL  ::       &
   zprmod_ref_f(nproma,pt_patch%nblks_c),              &
   zprrad_f(nproma,pt_patch%nblks_c)

! Local scalars:
! -------------

 REAL (KIND=wp)         ::       &
   zflar,                            & !sum of diagnostic/prognostic precipitation
   timefac,                          &
   zprmod_s,                         &
   zprmod_ref_s,                     &
   zprrad_s,                         &
   zprmod_ref_f_s,                   &
   zprrad_f_s

 REAL (KIND=wp)         ::       &
   realbuf  (7),      & ! for communication
   realbuf_g(7)         ! for communication


 INTEGER (KIND=i4) :: &
   jb,jc !,i_ver(nproma,pt_patch%nblks_c)

 INTEGER :: i_rlstart, i_rlend
 INTEGER :: i_startblk, i_endblk    !> blocks
 INTEGER :: i_startidx, i_endidx    !< slices

 INTEGER (KIND=i4)         ::           &
   zpranz, zprcount

 INTEGER, PARAMETER :: nthre=6

 INTEGER (KIND=i4) ::  &
   ass,bss,css,dss,zss       ! table of contengency

 INTEGER (KIND=i4) ::  &
   i,j,ii,jj,  &
   itab(7,7),i1,i2,ith, &
   histmod(7),histobs(7),anzobs,anzmod

 REAL (KIND=wp) ::  &
   rass,rbss,rcss,rdss    ! table of contengency as real

 REAL (KIND=wp)               ::           &
   hr, far, fr, fbi, ts, rets, rhss, ets, pod, tss, hss ! skill scores

 REAL (KIND=wp)               ::           &
   thr(6), thr_o(6)

   IF ( ytime == 'SW') THEN
      timefac=3600.0_wp
   ELSE
      timefac=1.0_wp
   ENDIF


   zprmod_s       = 0.0_wp
   zprmod_ref_s   = 0.0_wp
   zprmod_ref_f_s = 0.0_wp
   zprrad_s       = 0.0_wp
   zprrad_f_s     = 0.0_wp
   zprcount       = 0_i4
   zpranz         = 0_i4

   ! exclude boundary interpolation zone of nested domains
   i_rlstart = grf_bdywidth_c+1
   i_rlend   = min_rlcell_int

   i_startblk = pt_patch%cells%start_block(i_rlstart)
   i_endblk   = pt_patch%cells%end_block(i_rlend)


!!$OMP PARALLEL
!!$OMP DO PRIVATE(jc,i_startidx,i_endidx,zprmod_s,zprmod_ref_s,zprrad_s,zprcount,zprmod_ref_f_s,zprrad_f_s) ICON_OMP_GUIDED_SCHEDULE
   DO jb=i_startblk,i_endblk
       CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
         &                i_startidx, i_endidx, i_rlstart, i_rlend)
       DO jc=i_startidx,i_endidx
          IF (wobs_space(jc,jb) > 0.75_wp                         .AND.  &
              NINT(radar_data%radar_ct%blacklist(jc,jb)) /= 1_i4  .AND.  &
              NINT(lhn_fields%brightband(jc,jb)) /= 1_i4 .AND.  &
              zprrad(jc,jb) >= 0.0_wp) THEN
              zprmod_s        = zprmod_s        + zprmod(jc,jb)
              zprmod_ref_s    = zprmod_ref_s    + zprmod_ref(jc,jb)
              zprrad_s        = zprrad_s        + zprrad(jc,jb)
              IF (PRESENT (zprmod_ref_f)) zprmod_ref_f_s  = zprmod_ref_f_s  + zprmod_ref_f(jc,jb)
              IF (PRESENT (zprrad_f))     zprrad_f_s      = zprrad_f_s      + MAX(0._wp,zprrad_f(jc,jb))
              zprcount        = zprcount        + 1_i4
          ENDIF
       ENDDO
   ENDDO
!!$OMP END DO 
!!$OMP END PARALLEL

   zpranz = global_sum(zprcount,opt_iroot=p_io)
   zprmod_s = global_sum(zprmod_s,opt_iroot=p_io)
   zprmod_ref_s = global_sum(zprmod_ref_s,opt_iroot=p_io)
   zprrad_s = global_sum(zprrad_s,opt_iroot=p_io)
   IF (PRESENT (zprmod_ref_f)) zprmod_ref_f_s = global_sum(zprmod_ref_f_s,opt_iroot=p_io)
   IF (PRESENT (zprrad_f)) zprrad_f_s = global_sum(zprrad_f_s,opt_iroot=p_io)

   IF (my_process_is_stdio()) THEN

     IF (zpranz > 0) THEN
       zflar          = 1.0_wp / REAL (zpranz)
       zprmod_s       = zprmod_s       * zflar * timefac
       zprmod_ref_s   = zprmod_ref_s   * zflar * timefac
       zprrad_s       = zprrad_s       * zflar * timefac
       IF (PRESENT (zprmod_ref_f)) zprmod_ref_f_s = zprmod_ref_f_s * zflar * timefac
       IF (PRESENT (zprrad_f))     zprrad_f_s     = zprrad_f_s     * zflar * timefac

     ENDIF

      WRITE(nulhn, *)'Verification:'
      IF (ytime == "HR") THEN
        WRITE(nulhn, '(a,a3,f6.1,3f8.4)')'Modell (mod,ref,filt)',ytime,nsteps,zprmod_s,zprmod_ref_s,zprmod_ref_f_s
        WRITE(nulhn, '(a,a3,f6.1,2f9.4)')'Radar      (obs,filt)',ytime,nsteps,zprrad_s,zprrad_f_s
        WRITE(nulhn, '(a,a3,f6.1,f9.4)')'Statist (bias)',ytime,nsteps,zprmod_s-zprrad_s
      ELSE
        WRITE(nulhn, '(a,a3,f10.0,3f8.4)')'Modell (mod,ref,filt)',ytime,nsteps,zprmod_s,zprmod_ref_s,zprmod_ref_f_s
        WRITE(nulhn, '(a,a3,f10.0,2f9.4)')'Radar      (obs,filt)',ytime,nsteps,zprrad_s,zprrad_f_s
        WRITE(nulhn, '(a,a3,f10.0,f9.4)')'Statist (bias)',ytime,nsteps,zprmod_s-zprrad_s
      ENDIF


   ENDIF



   hr=0._wp
   fr=0._wp
   far=0._wp
   fbi=0._wp
   ts=0._wp
   ets=0._wp
   hss=0._wp
   tss=0._wp
   pod=0._wp
   rhss=0.0_wp
   rets=0.0_wp
   thr_o = (/ 0.1, 0.2, 0.5, 1.0, 2.0, 5.0 /)
   thr = thr_o / timefac
   histobs(:) = 0_i4
   histmod(:) = 0_i4
   anzobs=0_i4
   anzmod=0_i4

! contingence table:
!             |   Observed    |
!       -----------------------
!       |     |  yes  |   no  |
! -----------------------------
! Mod   | yes |  ass  |  bss  |
!       -----------------------
! elled | no  |  css  |  dss  |
! -----------------------------

   DO jj=1,nthre+1
    DO ii=1,nthre+1
     itab(ii,jj)=0
    ENDDO
   ENDDO

   i_startblk = pt_patch%cells%start_block(i_rlstart)
   i_endblk   = pt_patch%cells%end_block(i_rlend)


!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,i1,i2,histobs,histmod,ith,itab,anzobs,anzmod) ICON_OMP_GUIDED_SCHEDULE
   DO jb=i_startblk,i_endblk
       CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
         &                i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc=i_startidx,i_endidx
       IF (wobs_space(jc,jb) > 0.75_wp .AND. &
           NINT(radar_data%radar_ct%blacklist(jc,jb)) /= 1_i4 .AND. &
           NINT(lhn_fields%brightband(jc,jb)) /= 1_i4) THEN 
           i1=1
           i2=1
           IF (zprrad(jc,jb) > 0.0_wp) THEN
              histobs(1)=histobs(1)+1_i4
           ENDIF
           IF (zprmod(jc,jb) > 0.0_wp) THEN
              histmod(1)=histmod(1)+1_i4
           ENDIF
           DO ith=1,nthre
            IF (zprrad(jc,jb).GE.thr(ith)) THEN
             i1=i1+1
             histobs(ith+1)=histobs(ith+1)+1_i4
            ENDIF
            IF (zprmod(jc,jb).GE.thr(ith)) THEN
             i2=i2+1
             histmod(ith+1)=histmod(ith+1)+1_i4
            ENDIF
           ENDDO
           itab(i1,i2)=itab(i1,i2)+1
           anzobs=anzobs+1_i4
           anzmod=anzmod+1_i4
       ENDIF
     ENDDO
   ENDDO
!!$OMP END DO 
!!$OMP END PARALLEL

   DO ith=1,nthre

     ass=0_i4
     bss=0_i4
     css=0_i4
     dss=0_i4
     zss=0_i4
!  Observation no / Forecast no
     DO j=1,ith
      DO i=1,ith
       dss=dss+itab(i,j)
      ENDDO
     ENDDO
!  Observation no / Forecast yes
     DO j=ith+1,nthre+1
      DO i=1,ith
       bss=bss+itab(i,j)
      ENDDO
     ENDDO
!  Observation yes / Forecast no
     DO j=1,ith
      DO i=ith+1,nthre+1
       css=css+itab(i,j)
      ENDDO
     ENDDO
!  Observation yes / Forecast yes
     DO  j=ith+1,nthre+1
      DO  i=ith+1,nthre+1
       ass=ass+itab(i,j)
      ENDDO
     ENDDO


  zss=ass+bss+css+dss

! calculate skill scores

   realbuf = 0
   realbuf( 1)= REAL(ass)
   realbuf( 2)= REAL(bss)
   realbuf( 3)= REAL(css)
   realbuf( 4)= REAL(dss)
   realbuf( 5)= REAL(zss)
   realbuf( 6)= REAL(histmod(ith))
   realbuf( 7)= REAL(histobs(ith))

   realbuf_g(1:7) = global_sum(realbuf(1:7),opt_iroot=p_io)

  IF(my_process_is_stdio()) THEN

   ass    = INT(realbuf_g( 1))
   bss    = INT(realbuf_g( 2))
   css    = INT(realbuf_g( 3))
   dss    = INT(realbuf_g( 4))
   zss    = INT(realbuf_g( 5))
   histmod(ith) = INT(realbuf_g( 6))
   histobs(ith) = INT(realbuf_g( 7))
   rass=REAL(ass)
   rbss=REAL(bss)
   rcss=REAL(css)
   rdss=REAL(dss)

   ! new nomenclature according to U. Damrath in "Die neue Modellkette des DWD II"
   ! hitrate or percent correct score
   IF ((rass+rbss+rcss+rdss) > 0._wp) &
    hr  = 100.0_wp * (rass+rdss)/(rass+rbss+rcss+rdss)
   ! Probability of detection
   IF ((rass+rcss) > 0._wp) &
    pod  = 100.0_wp * (rass) / (rass+rcss)
   ! false alarm rate
   IF ((rbss+rdss) > 0._wp) &
    fr = 100.0_wp * (rbss) / (rdss+rbss)
   ! false alarm ratio
   IF ((rass+rbss) > 0._wp) &
    far = 100.0_wp * (rbss) / (rass+rbss)
   ! frequency bias
   IF ((rass+rcss) > 0._wp) &
    fbi = (rass+rbss) / (rass+rcss)
   ! threat score
   IF ((rbss+rcss+rass) > 0._wp) &
    ts  = 100.0_wp * (rass) / (rbss+rcss+rass)
   ! equitable threat score
   ! randomly correct forecasted wet points
   IF ((rass+rbss+rcss+rdss) > 0._wp) &
    rets = ((rass+rbss) * (rass+rcss)) / ((rass+rbss+rcss+rdss))
   IF ((rass+rbss+rcss-rets) > 0._wp) &
    ets = 100.0_wp * (rass-rets) / (rass+rbss+rcss-rets)
   ! True skill statistics or Hanssen-Kuipers discriminant or Kuipers score
   IF ((rass+rcss) > 0._wp .AND. (rbss+rdss) > 0._wp) &
    tss = ( (rass)/(rass+rcss) + (rdss)/(rbss+rdss) - 1._wp) * 100_wp
   ! Heidke skill score
   IF ((rass+rbss+rcss+rdss) > 0._wp) &
    rhss = (((rass+rbss)*(rass+rcss)+(rcss+rdss)*(rbss+rdss))/(rass+rbss+rcss+rdss))
   IF ((rass+rbss+rcss+rdss-rhss) > 0._wp) &
    hss = 100.0_wp * (rass+rdss-rhss)/(rass+rbss+rcss+rdss-rhss)

    IF (ytime == "HR") THEN
      WRITE(nulhn,'(a25,a3,f6.1,5i7,f5.2)')'skill scores (a,b,c,d):',ytime,nsteps,ass,bss,css,dss,zss,thr_o(ith)
      WRITE(nulhn,'(a14,a3,f6.1,9f8.2)')'skill scores:',ytime,nsteps,hr,pod,far,fr,fbi,ts,ets,hss,tss
    ELSE
      WRITE(nulhn,'(a25,a3,f10.0,5i7,f5.2)')'skill scores (a,b,c,d):',ytime,nsteps,ass,bss,css,dss,zss,thr_o(ith)
      WRITE(nulhn,'(a14,a3,f10.0,9f8.2)')'skill scores:',ytime,nsteps,hr,pod,far,fr,fbi,ts,ets,hss,tss
    ENDIF

  ENDIF

 ENDDO ! loop over thresholds

 realbuf = 0
 realbuf( 1)= REAL(histmod(7))
 realbuf( 2)= REAL(histobs(7))
 realbuf( 3)= REAL(anzobs)
 realbuf( 4)= REAL(anzmod)

 realbuf_g(1:4) = global_sum(realbuf(1:4),opt_iroot=p_io)

 IF(my_process_is_stdio()) THEN

   histmod(7) = INT(realbuf_g( 1))
   histobs(7) = INT(realbuf_g( 2))
   anzobs = INT(realbuf_g( 3))
   anzmod = INT(realbuf_g( 4))

   DO ith=1,6
      histobs(ith)=histobs(ith)-histobs(ith+1)
      histmod(ith)=histmod(ith)-histmod(ith+1)
   ENDDO
   IF (ytime == "HR") THEN
     WRITE(nulhn,'(a17,a3,f6.1,8i12)')'Histogramm model:',ytime,nsteps,anzmod,(histmod(i),i=1,7)
     WRITE(nulhn,'(a17,a3,f6.1,8i12)')'Histogramm radar:',ytime,nsteps,anzobs,(histobs(i),i=1,7)
   ELSE
     WRITE(nulhn,'(a17,a3,f10.0,8i12)')'Histogramm model:',ytime,nsteps,anzmod,(histmod(i),i=1,7)
     WRITE(nulhn,'(a17,a3,f10.0,8i12)')'Histogramm radar:',ytime,nsteps,anzobs,(histobs(i),i=1,7)
   ENDIF
 ENDIF

END SUBROUTINE lhn_verification

  !-------------------------------------------------------------------------
  SUBROUTINE open_lhn_log( )

!    CHARACTER (len=MAX_CHAR_LENGTH) :: file_ti    ! file name
    INTEGER :: istatus
    LOGICAL :: lopened

    yulhn = 'lhn.log'
    INQUIRE (FILE=TRIM(yulhn),OPENED=lopened)
    IF (lopened) THEN
      CALL message('LHN','lhn.log already open!')
      RETURN
    ENDIF
      
    nulhn = find_next_free_unit(100,1000)
    OPEN(UNIT=nulhn,FILE=TRIM(yulhn),ACTION="write", FORM='FORMATTED',IOSTAT=istatus)
    IF (istatus/=SUCCESS) THEN
      CALL finish('LHN','lhn.log already open')
    ELSE
      CALL message('LHN','opened lhn.log')
    ENDIF

  END SUBROUTINE open_lhn_log

!===============================================================================
! End of module
!-------------------------------------------------------------------------------

END MODULE mo_latent_heat_nudging
