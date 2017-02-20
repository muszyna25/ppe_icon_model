!>
!! Contains subroutines to calculate ozone climatology depending on day of the year.
!! Taken from DWD's GME.
!!
!! Also routines overtaken from ECHAM providing vertical interpolation between
!! obs-data and model levels (another time interpolation will be added later)
!!
!! @author Thorsten Reinhardt, AGeoBw, Offenbach
!! @author  M.A. Giorgetta, MPI, May  2000
!!
!! @par Revision History
!! Initial Release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-01-31)
!! ECHAM routine implemented by Kristina Froehlich, MPI-M (2011-08-17)
!! MACC data added by Martin Koehler, DWD (2015-05-12)
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

MODULE mo_o3_util

  USE mo_kind,                 ONLY: i8
  USE mo_exception,            ONLY: finish
  USE mo_parallel_config,      ONLY: nproma
  USE mo_impl_constants,       ONLY: min_rlcell_int, max_dom
  USE mo_kind,                 ONLY: wp
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_math_constants,       ONLY: pi,deg2rad,rad2deg
  USE mo_model_domain,         ONLY: t_patch
  USE mo_nh_vert_interp,       ONLY: prepare_lin_intp, lin_intp
  USE mo_nonhydro_types,       ONLY: t_nh_diag
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_o3_gems_data,         ONLY: rghg7
  USE mo_o3_macc_data,         ONLY: rghg7_macc
  USE mo_radiation_config,     ONLY: irad_o3
  USE mo_physical_constants,   ONLY: amd,amo3,rd,grav
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config, ltuning_ozone, icpl_o3_tp
  USE mo_time_config,          ONLY: time_config
  USE mo_impl_constants,       ONLY: io3_art
  USE mtime,                   ONLY: datetime, newDatetime, timedelta, newTimedelta, &
       &                             getPTStringFromMS, OPERATOR(+),                 &
       &                             NO_OF_MS_IN_A_MINUTE, NO_OF_MS_IN_A_HOUR,       &
       &                             NO_OF_MS_IN_A_SECOND,                           &       
       &                             getDayOfYearFromDatetime, MAX_TIMEDELTA_STR_LEN,&
       &                             deallocateTimedelta, deallocateDatetime
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights, &
       &                               calculate_time_interpolation_weights
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: calc_o3_clim, o3_timeint, o3_pl2ml, o3_zl2ml, calc_o3_gems


CONTAINS

  !=======================================================================

  SUBROUTINE o3_timeint( kproma,kbdim,nlev_pres,         & ! IN
                       & ext_o3, current_date,           & ! IN kproma,nlev_p,jb,nmonth
                       & o3_time_int                      )! OUT kproma,nlev_p

 ! In prior version, just a certain month was selected. This is not used anymore
 ! (revision13218) an the routine is now doing true time interpolation
 ! The time interpolation is done to the radiation time step, monthly ozone
 ! values provided that are given at the middle of each month (monthly averages)

    INTEGER, INTENT(in)         :: kproma ! 
    INTEGER, INTENT(in)         :: kbdim  ! 
    INTEGER, INTENT(in)         :: nlev_pres   ! number of o3 data levels
    TYPE(datetime), POINTER, INTENT(in) :: current_date
    REAL(wp), INTENT(in) , DIMENSION(kbdim,nlev_pres,0:13) :: ext_o3
    REAL(wp), INTENT(out), DIMENSION(kbdim,nlev_pres)      :: o3_time_int
    TYPE(t_time_interpolation_weights) :: tiw
    
    tiw = calculate_time_interpolation_weights(current_date)
    
    o3_time_int(1:kproma,:)=tiw%weight1*ext_o3(1:kproma,:,tiw%month1_index)+ &
                            tiw%weight2*ext_o3(1:kproma,:,tiw%month2_index)

  END SUBROUTINE o3_timeint

  SUBROUTINE o3_pl2ml ( kproma,kbdim,nlev_pres,klev,&
    &                   pfoz,phoz,ppf,pph,   &
    &                   o3_time_int, o3_clim)

    !- Description: o3 pressure levels to o3 sigma hybrid levels
    !
    !  The time interpolated profile of ozone is interpolated to the
    !  model full levels and integrated again fromp=0 to p=ps.
    !  Finally the ozone profile on the model levels is normalized such
    !  that the integrated amount on the model grid is identical to
    !  that on the grid of the climatology.
    !
    !- Author:
    !
    !  M.A. Giorgetta, MPI, May 2000
    !  K.   Froehlich, MPI, August 2011, adjusted to ICON

    ! INPUT
    ! -----

    INTEGER, INTENT(in)                             :: kproma ! 
    INTEGER, INTENT(in)                             :: kbdim  ! first dimension of 2-d arrays
    INTEGER, INTENT(in)                             :: klev   ! number of levels
    INTEGER, INTENT(in)                             :: nlev_pres   ! number of o3 data levels
    REAL(wp),INTENT(in) ,DIMENSION(nlev_pres)       :: pfoz  ! full level pressure of o3 data
    REAL(wp),INTENT(in) ,DIMENSION(nlev_pres+1)     :: phoz  ! half level pressure of o3 data
    REAL(wp),INTENT(in) ,DIMENSION(kbdim,klev)      :: ppf  ! full level pressure 
    REAL(wp), INTENT(in) ,DIMENSION(kbdim,klev+1)    :: pph  ! half level pressure
    REAL(wp),INTENT(in) ,DIMENSION(kbdim,nlev_pres) :: o3_time_int !zozonec_x
    REAL(wp),INTENT(OUT),DIMENSION(kbdim,klev)      :: o3_clim ! ozone in g/g


    ! LOCAL
    ! -----

    ! pressure integrated ozone at half levels of ozone grid,
    ! and integral at surface
    REAL(wp), DIMENSION(kproma)                 :: zozintc

    ! time interpolated ozone at full levels of model grid,
    ! pressure integrated ozone at half levels of model,
    ! and integral at surface
    REAL(wp), DIMENSION(kproma,klev)           :: zozonem 
    REAL(wp), DIMENSION(kproma)                :: zozintm

    REAL(wp) :: zdp1,zdp2 ! pressure weights in linear interpolation
    INTEGER  :: jl,jk,jkk ! loop indices
    INTEGER,DIMENSION(kproma)  :: jk1,jkn   ! first and last model level in interpolation

    INTEGER,DIMENSION(kproma)            :: kwork
    LOGICAL,DIMENSION(kproma)            :: kk_flag


    ! interpolate ozone profile to model grid
    ! ---------------------------------------
    ! set ozone concentration at levels above the uppermost level of
    ! the ozone climatology to the value in the uppermost level of 
    ! the ozone climatology

    jk1(:)     = 1
    kk_flag(:) = .TRUE.
    DO jk = 1,klev
       DO jl=1,kproma
          IF (ppf(jl,jk)<=pfoz(1) .AND. kk_flag(jl)) THEN 
             zozonem(jl,jk)= o3_time_int(jl,1)
             jk1(jl)=jk+1
          ELSE
             kk_flag(jl) = .FALSE.
          END IF
       END DO
    END DO

    ! set ozone concentration at levels below the lowermost level of
    ! the ozone climatology to the value in the lowermost level of 
    ! the ozone climatology
    jkn(:)=klev
    kk_flag(:) = .TRUE.
    DO jk = klev,1,-1
       DO jl=1,kproma
          IF (ppf(jl,jk)>=pfoz(nlev_pres).AND. kk_flag(jl)) THEN
             zozonem(jl,jk)=o3_time_int(jl,nlev_pres)
             jkn(jl)=jk-1
          ELSE
             kk_flag(jl) = .FALSE.
          END IF
       END DO
    ENDDO

    DO jk=1,klev
       kk_flag(:) = .TRUE.
       kwork(:)   = 1
       DO jkk = 1,nlev_pres
          DO jl=1,kproma
             IF(jk >= jk1(jl) .AND. jk <= jkn(jl))  THEN
                IF (ppf(jl,jk) <= pfoz(jkk) .AND. jkk >= kwork(jl) &
                                            .AND. kk_flag(jl)) THEN
                   kwork(jl)   = jkk
                   kk_flag(jl) = .FALSE.
                END IF
             END IF
          END DO
       END DO

       DO jl=1,kproma
          IF(jk >= jk1(jl) .AND. jk <= jkn(jl))  THEN
                jkk = kwork(jl)
                ! model level is in interval ]pfoz(jkk-1),pfoz(jkk)]
                ! -> make interpolation
                zdp1=pfoz(jkk)-ppf(jl,jk)
                zdp2=ppf(jl,jk)-pfoz(jkk-1)
                zozonem(jl,jk)=(zdp1*o3_time_int(jl,jkk-1) &
                          &   + zdp2*o3_time_int(jl,jkk))  &
                          &   /  (zdp1+zdp2)
          END IF
       END DO
    END DO

       ! integrate ozone profile on grid of climatology
       ! from top to surface
       ! ----------------------------------------------
    zozintc=0._wp
    kk_flag(:) = .TRUE.
    jk1(:)     = 2
    DO jk=2,nlev_pres+1
          ! integrate layers of climatology above surface
       DO jl=1,kproma
          IF (phoz(jk)<=pph(jl,klev) .AND. kk_flag(jl) ) THEN
              zozintc(jl)=zozintc(jl)+ &
               &  o3_time_int(jl,jk-1 )*(phoz(jk)-phoz(jk-1))
              jk1(jl) = jk+1
          ELSE
              kk_flag(jl) = .FALSE.
          END IF
       END DO
    END DO
       ! integrate layer of climatology that is intersected
       ! by the surface from upper boundary to surface
    DO jl=1,kproma
       zozintc(jl)=zozintc(jl)+ &
            &  o3_time_int(jl,jk1(jl)-1 )*(pph(jl,klev)-phoz(jk1(jl)-1))
    END DO

       ! integrate ozone profile on grid of model
       ! from top to surface
       ! ----------------------------------------
    zozintm=0._wp
    DO jk=2,klev+1
       DO jl=1,kproma
         zozintm(jl)=zozintm(jl) + zozonem(jl,jk-1)*(pph(jl,jk)-pph(jl,jk-1))
       END DO
    END DO

       ! normalize interpolated ozone profile such that the
       ! ozone integral computed on the model grid is equal
       ! to that integrated on the grid of the climatology
       ! --------------------------------------------------
    DO jk=1,klev
       DO jl=1,kproma
         o3_clim(jl,jk)=zozonem(jl,jk)/zozintm(jl) * zozintc(jl)
       END DO
    END DO

  END SUBROUTINE o3_pl2ml

  SUBROUTINE o3_zl2ml(nblks,npromz, nlev_o3,nlev, &
    &                 zf_o3, z_mc, ape_o3,  model_o3)

    !< Description: convert o3 on z-level to model height levels for the
    !! NH-model 
    !! The verical interpolation routines from mo_nh_vert_intp are used

    !!  @author Kristina Froehlich, MPI-M (2011-09-06)

    ! Input dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlev_o3    ! Number of input levels
    INTEGER , INTENT(IN) :: nlev       ! Number of output levels

    ! Input fields
    REAL(wp), DIMENSION(nproma,nlev_o3,nblks),INTENT(IN) :: zf_o3 ! height coordinate field of input data (m)
    REAL(wp), DIMENSION(nproma,nlev   ,nblks),INTENT(IN) :: z_mc  ! height coordinate field of model data (m)

    ! Input/Output field
    REAL(wp), DIMENSION(nproma,nlev_o3,nblks),INTENT(IN) :: ape_o3   ! external APE o3 data
    REAL(wp), DIMENSION(nproma,nlev   ,nblks),INTENT(OUT):: model_o3 ! converted o3 field for the model

    ! help field, since ape_o3 must remain as INTENT IN
    REAL(wp), DIMENSION(nproma,nlev_o3,nblks)            :: aux_o3   ! external APE o3 data

    ! weighting factors and indizes
    REAL(wp), DIMENSION(nproma,nlev,nblks) :: wfac ! weighting factor of upper level
    REAL(wp), DIMENSION(nproma,     nblks) :: wfacpbl1, wfacpbl2

    INTEGER,  DIMENSION(nproma,nlev,nblks) :: idx0    ! index of upper level
    INTEGER, DIMENSION(nproma,      nblks) :: kpbl1, kpbl2, bot_idx

    aux_o3(:,:,:)= ape_o3(:,:,:) 

    CALL prepare_lin_intp(z3d_in=zf_o3,  z3d_out=z_mc,               &
          &                  nblks= nblks,npromz= npromz,            &
          &                  nlevs_in=nlev_o3,      nlevs_out=nlev,  &
          &                  wfac=wfac, idx0=idx0  ,bot_idx=bot_idx  )

    CALL lin_intp        (f3d_in=aux_o3, f3d_out=model_o3,                  &
          &                  nblks= nblks,npromz=npromz,                    &
          &                  nlevs_in=nlev_o3,      nlevs_out=nlev,         &
          &                  wfac=wfac, idx0=idx0  ,bot_idx=bot_idx        ,&
          &                  kpbl1=kpbl1, wfacpbl1=wfacpbl1,                &
          &                  kpbl2=kpbl2, wfacpbl2=wfacpbl2                ,& 
          &                  l_loglin=.false., l_pd_limit=.true.,           &
          &                  l_extrapol=.false.                              )


  END SUBROUTINE o3_zl2ml


 !=======================================================================

  SUBROUTINE calc_o3_clim(kbdim,jg,p_inc_rad,z_sim_time,pt_patch,zvio3,zhmo3)

    INTEGER, INTENT(in)   :: &
      & kbdim, jg

    REAL(wp), INTENT(in)   :: &
      & p_inc_rad, & !radiation time step in seconds
      & z_sim_time

    TYPE(t_patch),      INTENT(in)    :: pt_patch    ! Patch

    REAL(wp),      INTENT(inout)  :: &
      & zvio3(kbdim,pt_patch%nblks_c) , &
      & zhmo3(kbdim,pt_patch%nblks_c)

    REAL(wp) ::                     &
      & z_sim_time_rad,  &
      & zstunde,                   &
      & ztwo  , ztho

    ! output from o3_par_t5:
    REAL(wp) :: zvio3_c(21),zvio3_s(15)
    REAL(wp) :: zhmo3_c(21),zhmo3_s(15)

    ! output from legtri_vec
    REAL(wp) :: zalp   (kbdim,21)

    REAL(wp) :: zvio3_f(kbdim,11)
    REAL(wp) :: zhmo3_f(kbdim,11)

    REAL(wp) :: zsin1(pt_patch%nblks_c),zcos1(pt_patch%nblks_c)   ! sine and cosine of longitude
    REAL(wp) :: zsin2(pt_patch%nblks_c),zcos2(pt_patch%nblks_c)   ! dito for 2*longitude
    REAL(wp) :: zsin3(pt_patch%nblks_c),zcos3(pt_patch%nblks_c)   ! dito ...
    REAL(wp) :: zsin4(pt_patch%nblks_c),zcos4(pt_patch%nblks_c)   !
    REAL(wp) :: zsin5(pt_patch%nblks_c),zcos5(pt_patch%nblks_c)   !

    
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
!!$    INTEGER :: jk
    
    INTEGER :: &
      & jj, itaja, jb

    INTEGER , SAVE :: itaja_o3_previous(max_dom) = 0
    
    INTEGER :: jmm,mmm,mnc,mns,jnn,jc

    TYPE(datetime), POINTER :: current => NULL()
    TYPE(timedelta), POINTER :: td => NULL()
    CHARACTER(len=MAX_TIMEDELTA_STR_LEN) :: td_string 
    
    i_nchdom  = MAX(1,pt_patch%n_childdom)
    
    z_sim_time_rad = z_sim_time + 0.5_wp*p_inc_rad

    current => newDatetime(time_config%tc_exp_startdate)
    CALL getPTStringFromMS(INT(1000.0_wp*z_sim_time_rad,i8), td_string)
    td => newTimedelta(td_string)
    current = time_config%tc_exp_startdate + td
    jj = INT(current%date%year)
    itaja = getDayOfYearFromDateTime(current)
    zstunde = current%time%hour+( &
         &    REAL(current%time%minute*NO_OF_MS_IN_A_MINUTE &
         &        +current%time%second*NO_OF_MS_IN_A_SECOND &
         &        +current%time%ms,wp)/REAL(NO_OF_MS_IN_A_HOUR,wp))
    CALL deallocateDatetime(current)
    CALL deallocateTimedelta(td)

    !decide whether new ozone calculation is necessary
    IF ( itaja == itaja_o3_previous(jg) ) RETURN
    itaja_o3_previous(jg) = itaja
    
    ztwo    = 0.681_wp + 0.2422_wp*REAL(jj-1949,wp)-REAL((jj-1949)/4,wp)
    ztho    = 2._wp*pi*( REAL(itaja, wp) -1.0_wp + ztwo )/365.2422_wp


!
!     ozone parameters
!     ----------------

!     a) time of the year dependency
    CALL o3_par_t5 (ztho,zvio3_c,zvio3_s,zhmo3_c,zhmo3_s)

    !     b) inverse Legendre transform

    ! nest boudaries have to be included for reduced-grid option
    rl_start = 1
    rl_end   = min_rlcell_int

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)
      
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, zalp, jmm, jc, zvio3_f ,zhmo3_f, mmm, mnc,&
!$OMP  mns, jnn )  ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
&                         i_startidx, i_endidx, rl_start, rl_end)
     
      CALL legtri_vec( kbdim,i_startidx,i_endidx,pt_patch%cells%center(1:kbdim,jb)%lat, 6, zalp )
   
        DO jmm=1,11    
         DO jc=i_startidx,i_endidx
           zvio3_f(jc,jmm)=0._wp
           zhmo3_f(jc,jmm)=0._wp
         END DO
        END DO    ! loop over first horizontal index
 
        mmm=0 
        mnc=0 
        mns=0 
        DO jmm=1,6 
          mmm=mmm+1 
          DO jnn=jmm,6 
            mnc=mnc+1   
            DO jc=i_startidx,i_endidx
              zvio3_f(jc,mmm)=zvio3_f(jc,mmm)+zalp(jc,mnc)*zvio3_c(mnc)  
              zhmo3_f(jc,mmm)=zhmo3_f(jc,mmm)+zalp(jc,mnc)*zhmo3_c(mnc)  
            END DO    ! loop over first horizontal index
          END DO                                           
          IF(jmm.NE.1) THEN                             
            mmm=mmm+1                                    
            DO jnn=jmm,6                              
              mns=mns+1                                
              DO jc=i_startidx,i_endidx
                zvio3_f(jc,mmm)=zvio3_f(jc,mmm)+zalp(jc,mns+6)*zvio3_s(mns)
                zhmo3_f(jc,mmm)=zhmo3_f(jc,mmm)+zalp(jc,mns+6)*zhmo3_s(mns)
              END DO    ! loop over first horizontal index
           END DO                                           
          END IF                                            
        END DO

!     c) inverse Fourier transform

        DO jc=i_startidx,i_endidx
          zcos1(jb)=COS (pt_patch%cells%center(jc,jb)%lon)  
          zsin1(jb)=SIN (pt_patch%cells%center(jc,jb)%lon) 
          zcos2(jb)=zcos1(jb)*zcos1(jb)-zsin1(jb)*zsin1(jb)   
          zsin2(jb)=zsin1(jb)*zcos1(jb)+zcos1(jb)*zsin1(jb)   
          zcos3(jb)=zcos2(jb)*zcos1(jb)-zsin2(jb)*zsin1(jb)  
          zsin3(jb)=zsin2(jb)*zcos1(jb)+zcos2(jb)*zsin1(jb)  
          zcos4(jb)=zcos3(jb)*zcos1(jb)-zsin3(jb)*zsin1(jb)  
          zsin4(jb)=zsin3(jb)*zcos1(jb)+zcos3(jb)*zsin1(jb)  
          zcos5(jb)=zcos4(jb)*zcos1(jb)-zsin4(jb)*zsin1(jb) 
          zsin5(jb)=zsin4(jb)*zcos1(jb)+zcos4(jb)*zsin1(jb) 

!       vertically integrated ozone amount (Pa O3)
          zvio3(jc,jb) = zvio3_f(jc, 1) + 2._wp *                       &
                       ( zvio3_f(jc, 2)*zcos1(jb) + zvio3_f(jc, 3)*zsin1(jb)+   &
                         zvio3_f(jc, 4)*zcos2(jb) + zvio3_f(jc, 5)*zsin2(jb)+   &
                         zvio3_f(jc, 6)*zcos3(jb) + zvio3_f(jc, 7)*zsin3(jb)+   &
                         zvio3_f(jc, 8)*zcos4(jb) + zvio3_f(jc, 9)*zsin4(jb)+   &
                         zvio3_f(jc,10)*zcos5(jb) + zvio3_f(jc,11)*zsin5(jb) )

!       pressure height of ozone maximum position (Pa)
          zhmo3(jc,jb) = zhmo3_f(jc, 1) + 2._wp *                       &
                       ( zhmo3_f(jc, 2)*zcos1(jb) + zhmo3_f(jc, 3)*zsin1(jb)+   &
                         zhmo3_f(jc, 4)*zcos2(jb) + zhmo3_f(jc, 5)*zsin2(jb)+   &
                         zhmo3_f(jc, 6)*zcos3(jb) + zhmo3_f(jc, 7)*zsin3(jb)+   &
                         zhmo3_f(jc, 8)*zcos4(jb) + zhmo3_f(jc, 9)*zsin4(jb)+   &
                         zhmo3_f(jc,10)*zcos5(jb) + zhmo3_f(jc,11)*zsin5(jb) )

        END DO    ! loop over first horizontal index

!!$        ! 3-dimensional O3
!!$        ! top level
!!$        DO jc = i_startidx,i_endidx
!!$          zptop32  (jc,jb) = (SQRT(z_pres_ifc(jc,1,jb)))**3
!!$          zo3_hm   (jc,jb) = (SQRT(zhmo3(jc,jb)))**3
!!$          zo3_top  (jc,jb) = prm_diag%vio3(jc,jb)*zptop32(jc,jb)/(zptop32(jc,jb)+zo3_hm(jc,jb))
!!$        ENDDO
!!$        ! loop over layers
!!$        DO jk = 1,nlev   
!!$          DO jc = i_startidx,i_endidx
!!$            zpbot32  (jc,jb) = (SQRT(z_pres_ifc(jc,jk+1,jb)))**3
!!$            zo3_bot  (jc,jb) = prm_diag%vio3(jc,jb)* zpbot32(jc,jb)    &
!!$              /( zpbot32(jc,jb) + zo3_hm(jc,jb))        
!!$            !O3 content
!!$            zduo3(jc,jk,jb)= zo3_bot (jc,jb)-zo3_top (jc,jb)
!!$            z_o3(jc,jk,jb) = (amo3/amd) * &
!!$              & ((z_pres_ifc(jc,jk+1,jb)-z_pres_ifc(jc,jk,jb))/zduo3(jc,jk,jb))
!!$            ! store previous bottom values in arrays for top of next layer
!!$            zo3_top (jc,jb) = zo3_bot (jc,jb)
!!$          ENDDO
!!$        ENDDO

    ENDDO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL 

  END SUBROUTINE calc_o3_clim
  
  
  !! Calculates the parameters of a T5 spectral distribution of o3 depending on the time of year.
  !! It is a copy of SUBROUTINE ozone (in ozone.f90) of DWD's GME (contentually identically, only 
  !! formally slightly adapted).
  !!
  SUBROUTINE o3_par_t5 (pytime,pvio3_c,pvio3_s,phmo3_c,phmo3_s)

!=======================================================================
!
! Current Code Owner: DWD, B. Ritter
!    phone: +49-69-8062-2703, fax: +49-69-8062-3721
!    email: bodo.ritter@dwd.de
!
!=======================================================================
!
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        2000/01/28 B. Ritter   based on J.F.Geleyn
!  Initial Release
! 1.14       2002/01/16 Helmut P. Frank
!  Introduce KIND-notation for all reals with real2kind.pl
! V2_24        2010/04/28 Michael Gertz
!  Adaptions to SVN
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!     The routine needs to be called only once per radiation time step,
!     but is called from the radiation interface routine *parrad* each
!     time this interface routine is called.
!     The routine computes seperate parameters for the vertically
!     integrated ozone amount (pvio3_c_) and for the height of the ozone
!     maximum (phmo3_) from the actual time of the year.

!     Input:

      REAL(wp), INTENT(in) :: pytime      ! time of year   (sr)

!     Output:

!     arrays for T5-distributions of ozone integral and height of maximum
!     (vi: vertical integral, hm: height of maximum; c:cosine, s:sine)

      REAL(wp), INTENT(out) :: pvio3_c(21),pvio3_s(15)
      REAL(wp), INTENT(out) :: phmo3_c(21),phmo3_s(15)

!     Local variables and arrays:

!     coefficients for the time fourier development
!     (.vio3../.hmo3.. stand for vertical integral and height of maximum;
!      .....cn/.....sn stand for the (n=0,4) terms of the fourier
!      development)

      REAL(wp) :: zvio3c0 (21),zvio3s0(15)
      REAL(wp) :: zvio3c1 (21),zvio3s1(15)
      REAL(wp) :: zvio3c2 (21),zvio3s2(15)
      REAL(wp) :: zvio3c3 (21),zvio3s3(15)
      REAL(wp) :: zvio3c4 (21),zvio3s4(15)
      REAL(wp) :: zhmo3c0 (21),zhmo3s0(15)
      REAL(wp) :: zhmo3c1 (21),zhmo3s1(15)
      REAL(wp) :: zhmo3c2 (21),zhmo3s2(15)
      REAL(wp) :: zhmo3c3 (21),zhmo3s3(15)
      REAL(wp) :: zhmo3c4 (21),zhmo3s4(15)

      REAL(wp) :: zc1yt          ! cosine of time of year
      REAL(wp) :: zs1yt          ! sine   of time of year
      REAL(wp) :: zc2yt          ! 2*cosine of time of year
      REAL(wp) :: zs2yt          ! 2*sine   of time of year

      INTEGER jmn         ! loop index
!
!

!     data statements:

      DATA zvio3c0/                                                                            &
     &+.6012E-01_wp,+.1887E-02_wp,+.7410E-02_wp,+.9950E-03_wp,-.1426E-02_wp,-.2072E-03_wp,     &
     &                  -.4954E-03_wp,+.7955E-05_wp,-.3701E-03_wp,+.4116E-04_wp,-.4163E-04_wp, &
     &                                -.2933E-03_wp,+.2154E-04_wp,-.2849E-03_wp,-.1604E-03_wp, &
     &                                              -.1054E-03_wp,+.4974E-03_wp,+.1047E-03_wp, &
     &                                                            +.8323E-04_wp,+.2874E-03_wp, &
     &                                                                          +.1333E-03_wp/
      DATA zvio3s0/                                                                            &
     &                  +.4210E-03_wp,-.9591E-03_wp,+.2811E-03_wp,-.2257E-03_wp,-.1713E-03_wp, &
     &                                -.3538E-03_wp,+.1095E-03_wp,-.4390E-03_wp,-.5605E-05_wp, &
     &                                              +.1478E-03_wp,+.2849E-03_wp,+.3430E-03_wp, &
     &                                                            +.8248E-04_wp,+.1442E-03_wp, &
     &                                                                            -.1375E-04_wp/
      DATA zhmo3c0/                                                                            &
     &    +.3166E+04_wp,+.8663E+02_wp,+.9401E+03_wp,+.1999E+02_wp,-.3530E+03_wp,-.3311E+02_wp, &
     &                  -.4903E+02_wp,-.4015E+00_wp,-.1333E+02_wp,+.5675E+01_wp,+.7221E+01_wp, &
     &                                -.3001E+02_wp,+.7570E+01_wp,-.1142E+02_wp,-.1365E+02_wp, &
     &                                              -.1502E+02_wp,+.4911E+02_wp,+.1425E+02_wp, &
     &                                                            +.8983E+01_wp,+.3064E+02_wp, &
     &                                                                            +.1693E+02_wp/
      DATA zhmo3s0/                                                                            &
     &                  +.4231E+02_wp,-.7391E+02_wp,+.1273E+02_wp,+.2086E+02_wp,-.1597E+02_wp, &
     &                                -.3591E+02_wp,+.1059E+02_wp,-.2779E+02_wp,-.6923E+01_wp, &
     &                                              +.1397E+02_wp,+.2387E+02_wp,+.2883E+02_wp, &
     &                                                            +.8626E+01_wp,+.1607E+02_wp, &
     &                                                                            -.2676E+01_wp/
      DATA zvio3c1/                                                                            &
     &    +.7090E-04_wp,+.4930E-05_wp,+.6829E-03_wp,+.1897E-03_wp,+.7226E-04_wp,-.2807E-03_wp, &
     &                  +.4970E-04_wp,-.1753E-03_wp,-.7843E-04_wp,-.1649E-03_wp,-.1037E-03_wp, &
     &                                -.4830E-04_wp,-.6304E-04_wp,-.1100E-03_wp,-.7952E-04_wp, &
     &                                              +.1326E-04_wp,+.2599E-04_wp,+.9926E-05_wp, &
     &                                                            -.9247E-05_wp,-.3521E-05_wp, &
     &                                                                            -.1780E-04_wp/
      DATA zvio3s1/                                                      &
     &                  +.6333E-04_wp,+.1145E-03_wp,+.1192E-03_wp,+.4934E-04_wp,+.2699E-04_wp, &
     &                                +.3684E-04_wp,-.2395E-05_wp,+.2045E-04_wp,-.8684E-04_wp, &
     &                                              +.5301E-04_wp,-.4176E-05_wp,+.4103E-04_wp, &
     &                                                            +.2783E-04_wp,+.1754E-04_wp, &
     &                                                                            +.1116E-04_wp/
      DATA zhmo3c1/                                                                            &
     &-.3450E+02_wp,+.2148E+03_wp,+.3376E+02_wp,+.6535E+02_wp,-.1564E+02_wp,-.4273E+02_wp,     &
     &                  +.9553E+01_wp,-.4647E+01_wp,-.6129E+01_wp,-.6727E+01_wp,-.6761E+01_wp, &
     &                                -.2467E+01_wp,-.2181E+01_wp,-.5361E+01_wp,-.2395E+01_wp, &
     &                                              +.5952E+00_wp,+.2106E+01_wp,-.1367E+01_wp, &
     &                                                            -.2349E+01_wp,+.3532E+00_wp, &
     &                                                                            -.3169E+01_wp/
      DATA zhmo3s1/                                                                            &
     &                  +.3977E+01_wp,+.5032E+01_wp,+.6226E+01_wp,-.3625E+00_wp,-.1373E+01_wp, &
     &                                +.4600E+01_wp,+.4312E+01_wp,+.2882E+01_wp,-.6351E+01_wp, &
     &                                              +.5731E+01_wp,-.2574E+01_wp,+.3235E+00_wp, &
     &                                                            +.2806E+01_wp,+.8133E+00_wp, &
     &                                                                            +.2032E+01_wp/
      DATA zvio3c2/                                                                            &
     &    +.8571E-03_wp,+.3086E-02_wp,+.9287E-03_wp,+.2787E-03_wp,+.1826E-03_wp,-.1006E-03_wp, &
     &                  +.1092E-03_wp,-.1266E-03_wp,+.5372E-04_wp,-.1188E-03_wp,-.3285E-04_wp, &
     &                                -.1783E-04_wp,-.3018E-05_wp,-.8709E-04_wp,-.8707E-04_wp, &
     &                                              +.8633E-04_wp,+.3530E-04_wp,+.4863E-04_wp, &
     &                                                            +.3917E-05_wp,-.3252E-04_wp, &
     &                                                                            -.1936E-06_wp/
      DATA zvio3s2/                                                                            &
     &                  -.8822E-04_wp,+.1341E-03_wp,+.3095E-04_wp,+.8230E-04_wp,+.2735E-04_wp, &
     &                                +.1714E-04_wp,-.9406E-04_wp,+.1912E-04_wp,-.5402E-04_wp, &
     &                                              +.3571E-04_wp,+.3897E-04_wp,+.4487E-04_wp, &
     &                                                            +.3079E-04_wp,+.3196E-04_wp, &
     &                                                                            -.2391E-05_wp/
      DATA zhmo3c2/                                                                            &
     &    +.5216E+02_wp,+.1613E+03_wp,+.3284E+02_wp,-.7670E+02_wp,-.9548E+01_wp,+.1608E+02_wp, &
     &                  +.1023E+02_wp,-.1090E+02_wp,+.2748E+01_wp,-.3846E+01_wp,-.4135E+01_wp, &
     &                                +.1255E+01_wp,-.3301E-01_wp,-.5273E+01_wp,-.7247E+01_wp, &
     &                                              +.1387E+02_wp,+.4184E+01_wp,+.6495E+01_wp, &
     &                                                            +.2944E+01_wp,-.1947E+01_wp, &
     &                                                                            +.1132E+01_wp/
      DATA zhmo3s2/                                                                            &
     &                  -.1968E+02_wp,+.1192E+02_wp,-.1194E+01_wp,+.1084E+01_wp,+.2946E+01_wp, &
     &                                +.2630E+01_wp,-.1256E+02_wp,+.1395E+01_wp,-.2222E+01_wp, &
     &                                              +.4864E+01_wp,+.6450E+01_wp,+.5568E+01_wp, &
     &                                                            +.5292E+01_wp,+.4876E+01_wp, &
     &                                                                            -.7579E+00_wp/
      DATA zvio3c3/                                                                            &
     &    -.2759E-03_wp,-.2781E-03_wp,-.1087E-03_wp,-.1633E-03_wp,-.3627E-04_wp,-.4242E-04_wp, &
     &                  +.6045E-05_wp,-.1703E-04_wp,+.4562E-04_wp,-.1009E-04_wp,+.2663E-04_wp, &
     &                                -.1786E-04_wp,+.1550E-04_wp,-.9135E-06_wp,+.2372E-04_wp, &
     &                                              +.1100E-05_wp,+.2299E-04_wp,+.4659E-05_wp, &
     &                                                            +.2423E-05_wp,+.7321E-05_wp, &
     &                                                                            +.8852E-05_wp/
      DATA zvio3s3/                                                                            &
     &                  -.3678E-04_wp,-.2219E-04_wp,-.3911E-04_wp,-.4398E-04_wp,-.1142E-04_wp, &
     &                                -.9121E-05_wp,-.2011E-04_wp,+.4711E-06_wp,-.3775E-05_wp, &
     &                                              +.3866E-05_wp,+.2400E-04_wp,+.2043E-04_wp, &
     &                                                            -.1824E-05_wp,-.5550E-05_wp, &
     &                                                                            +.2506E-05_wp/
      DATA zhmo3c3/                                                                            &
     &    -.1534E+03_wp,-.2095E+02_wp,-.1006E+03_wp,-.7385E+01_wp,+.5203E+01_wp,+.9434E+00_wp, &
     &                  -.3814E+00_wp,-.3175E+01_wp,+.3366E+01_wp,+.3378E+00_wp,+.2740E+00_wp, &
     &                                -.2669E+01_wp,+.8452E+00_wp,+.3498E+00_wp,+.2192E+01_wp, &
     &                                              -.4024E+00_wp,+.1544E+01_wp,-.4588E+00_wp, &
     &                                                            +.6998E+00_wp,+.6263E+00_wp, &
     &                                                                            +.1228E+01_wp/
      DATA zhmo3s3/                                                                            &
     &                  -.3588E+01_wp,+.2076E+00_wp,-.2088E+01_wp,-.4159E+01_wp,+.2244E+00_wp, &
     &                                -.7751E+00_wp,-.2749E+01_wp,+.7234E+00_wp,+.4390E+00_wp, &
     &                                              -.1646E+00_wp,+.1700E+01_wp,+.1046E+01_wp, &
     &                                                            -.7856E+00_wp,-.1644E+01_wp, &
     &                                                                            +.2648E+00_wp/
      DATA zvio3c4/                                                      &
     &    -.1460E-03_wp,+.3422E-03_wp,-.3529E-04_wp,+.1791E-03_wp,-.1917E-03_wp,-.2558E-04_wp, &
     &                  +.6547E-04_wp,+.6401E-04_wp,+.4823E-04_wp,+.7084E-05_wp,+.2895E-04_wp, &
     &                                -.1561E-04_wp,+.8179E-06_wp,+.1028E-04_wp,-.7667E-05_wp, &
     &                                              -.4347E-05_wp,+.7293E-05_wp,-.5735E-05_wp, &
     &                                                            +.7838E-05_wp,-.2933E-05_wp, &
     &                                                                            +.3686E-05_wp/
      DATA zvio3s4/                                                                            &
     &                  -.4560E-05_wp,-.5292E-04_wp,-.1252E-04_wp,+.1850E-04_wp,-.2273E-04_wp, &
     &                                +.6552E-05_wp,+.1422E-04_wp,-.6545E-05_wp,+.7998E-06_wp, &
     &                                              +.2845E-04_wp,+.2497E-04_wp,+.2844E-04_wp, &
     &                                                            +.3855E-06_wp,-.1487E-04_wp, &
     &                                                                            +.1954E-05_wp/
      DATA zhmo3c4/                                                                            &
     &    +.9260E+01_wp,-.9055E+01_wp,+.5460E+01_wp,-.7603E+01_wp,-.3329E+02_wp,-.1048E+02_wp, &
     &                  +.9328E+01_wp,+.4597E+01_wp,+.3827E+01_wp,-.3201E+01_wp,+.1708E+01_wp, &
     &                                -.1548E+01_wp,-.5323E+00_wp,+.3039E+01_wp,+.5740E+00_wp, &
     &                                              +.1353E+00_wp,-.2354E+01_wp,+.2818E+00_wp, &
     &                                                            +.1113E+01_wp,-.1891E+01_wp, &
     &                                                                            -.3074E+00_wp/
      DATA zhmo3s4/                                                                            &
     &                  -.2446E+01_wp,+.4199E+01_wp,-.2571E+01_wp,+.8194E+01_wp,+.4206E+00_wp, &
     &                                +.3856E+01_wp,+.1159E+01_wp,+.2547E+01_wp,-.1314E+01_wp, &
     &                                              +.2331E+01_wp,+.1144E+01_wp,-.4408E+00_wp, &
     &                                                            -.6797E+00_wp,-.2598E+01_wp, &
     &                                                                            +.8953E+00_wp/
!
!     sine and cosine of time of year
      zc1yt = COS(pytime)
      zs1yt = SIN(pytime)
      zc2yt = zc1yt**2-zs1yt**2
      zs2yt = 2._wp*zs1yt*zc1yt

!     second order Fourier development
!     ================================

!     cosine coefficients
      DO jmn=1,21
        pvio3_c(jmn)=zvio3c0(jmn) + 2._wp*                   &
          &            (zvio3c1(jmn)*zc1yt+zvio3c2(jmn)*zs1yt   &
          &            +zvio3c3(jmn)*zc2yt+zvio3c4(jmn)*zs2yt)
        phmo3_c(jmn)=zhmo3c0(jmn) + 2._wp*                   &
          &            (zhmo3c1(jmn)*zc1yt+zhmo3c2(jmn)*zs1yt   &
          &            +zhmo3c3(jmn)*zc2yt+zhmo3c4(jmn)*zs2yt)
      END DO

!     sine coefficients
      DO jmn=1,15
        pvio3_s(jmn)=zvio3s0(jmn) + 2._wp*                   &
          &            (zvio3s1(jmn)*zc1yt+zvio3s2(jmn)*zs1yt   &
          &            +zvio3s3(jmn)*zc2yt+zvio3s4(jmn)*zs2yt)
        phmo3_s(jmn)=zhmo3s0(jmn) + 2._wp*                   &
          &            (zhmo3s1(jmn)*zc1yt+zhmo3s2(jmn)*zs1yt   &
          &            +zhmo3s3(jmn)*zc2yt+zhmo3s4(jmn)*zs2yt)
      END DO

  END SUBROUTINE o3_par_t5


  !! This a copy from legtri_vec.f90 of DWD's GME (content identical, only formally
  !! slightly adapted).
  !!
  SUBROUTINE legtri_vec (kbdim, ki1sc, ki1ec, plat, kcp,  palp )
!
!**** *LEGTRI* - *LEGENDRE FUNKTION FUER EINE DREIECK-ABSCHNEIDUNG
!  This is a vectorizable version of the subroutine legtri.
!  The directives are for NEC SX vector computers.
!
!=======================================================================
!
! Current Code Owner: DWD, B. Ritter
!    phone: +49-69-8062-2703, fax: +49-69-8062-3721
!    email: bodo.ritter@dwd.de
!
!=======================================================================
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V2_16        2008/06/23 J.-O. Beismann, D. Liermann, H. Frank
!  Initial release
! V2_24        2010/04/28 Michael Gertz
!  Adaptions to SVN
!
!     PURPOSE.
!     --------
!
!     THIS ROUTINE COMPUTES THE VALUES *PALP* FOR THE ARGUMENT
!     *PSIN* OF THE NORMALISED *LEGENDRE ASSOCIATED FUNCTIONS IN THE
!     ORDER ((JN1=JM1,KCP),JM1=1,KCP) FOR JN=JN1-1 AND JM=JM1-1 .
!
!     plat is the latitude in radians.
!
!**   INTERFACE.
!     ----------
!
!          THERE ARE THREE ARGUMENTS: *plat* IS THE LATITUDE.
!                                     *KCP*  IS ONE PLUS THE LIMIT WAVE NUMBER.
!                                     *PALP* IS THE ARRAY OF THE RESULTS.
!                                     ki1sc,,...,ki2ec are dimensional limits of plat, palp.
!
!     METHOD.
!     -------
!
!          SIMPLE RECURENCE FORMULA.
!

   
      INTEGER,  INTENT(IN)  :: kbdim, kcp, ki1sc, ki1ec
      REAL(wp), INTENT(IN)  :: plat( kbdim )

      REAL(wp), INTENT(OUT) :: palp( kbdim, kcp*(kcp+1)/2 )

      REAL(wp)              ::  zsin( kbdim ), &
                                zcos( kbdim ), &
                                zf1m( kbdim ), &
                                zf2m( kbdim )

      REAL(wp)              :: ZE1, ZE2, Z2M, ZM, ZN, ZN2, ZRE1, zsqrt3, zrs2m

      INTEGER ::  IC, ICP, JJ, JM, JM1, JM2, JN
      INTEGER ::  j1

!**********************************************************************
!
!*         1.     PRELIMINARY SETTING.
!                 ----------- --------
!
      ICP= KCP
      IC = ICP-1
      JJ = 2
      zsqrt3 = SQRT(3.0_wp)
      zf1m(ki1sc:ki1ec) = zsqrt3
!
!     ------------------------------------------------------------------
!
!*         2.     COMPUTATIONS.
!                 -------------
!
      palp(ki1sc:ki1ec,1) = 1._wp
      zsin(ki1sc:ki1ec)   = SIN( plat(ki1sc:ki1ec) )
      zcos(ki1sc:ki1ec)   = SQRT ( MAX ( 0._wp, 1._wp-zsin(ki1sc:ki1ec)**2 ) )
      palp(ki1sc:ki1ec,2) = zsqrt3*zsin(ki1sc:ki1ec)

      DO JM1=1,ICP
        JM  = JM1-1
        ZM  = REAL(JM,wp)
        Z2M = ZM+ZM
        ZRE1= SQRT(Z2M+3._wp)
        ZE1 = 1._wp/ZRE1
        IF(JM.EQ.0) GO TO 201

        zrs2m = 1._wp/SQRT(Z2M)
!CDIR COLLAPSE
        DO j1=ki1sc,ki1ec
          !  The following line produces results identical to legtri
          !           ZF2M(j1,j2) = ZF1M(j1,j2)*ZCOS(j1,j2)/SQRT(Z2M)
          !  However, the next line is approx. 30 % faster on ibm power 5
          ZF2M(j1) = ZF1M(j1)*ZCOS(j1)*zrs2m
          ZF1M(j1) = ZF2M(j1)*ZRE1
        ENDDO
        JJ=JJ+1
!CDIR COLLAPSE
        DO j1=ki1sc,ki1ec
          PALP(j1,JJ) = ZF2M(j1)
        ENDDO
        IF(JM.EQ.IC) CYCLE

        JJ = JJ+1
!CDIR COLLAPSE
        DO j1=ki1sc,ki1ec
          PALP(j1,JJ) = ZF1M(j1)*ZSIN(j1)
        ENDDO
        IF(JM1.EQ.IC) CYCLE

  201   CONTINUE
        JM2 = JM+2

        DO JN=JM2,IC
          ZN = REAL(JN,wp)
          ZN2= ZN**2
          ZE2= SQRT((4._wp*ZN2-1._wp)/(ZN2-ZM**2))
          JJ = JJ+1
!CDIR COLLAPSE
            DO j1=ki1sc,ki1ec
              PALP(j1,JJ)=ZE2*( ZSIN(j1)*PALP(j1,JJ-1)   &
                                          -ZE1*PALP(j1,JJ-2) )
            ENDDO
          ZE1 = 1._wp/ZE2
        END DO
      END DO
!

  END SUBROUTINE legtri_vec

  !>
  !! Calculates GEMS ozone climatology.
  !! Taken and adapted from ECMWF's IFS (37r2).
  !!
  !! @par Revision History
  !! Initial Release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-10-18)
  !!
  SUBROUTINE calc_o3_gems(pt_patch,mtime_datetime,p_diag,prm_diag,ext_data)

    TYPE(t_patch),           INTENT(in) :: pt_patch    ! Patch
    TYPE(datetime), POINTER, INTENT(in) :: mtime_datetime
    TYPE(t_nh_diag),         INTENT(in) :: p_diag  !!the diagostic variables
    TYPE(t_nwp_phy_diag),    INTENT(in):: prm_diag

    TYPE(t_external_data),   INTENT(inout) :: ext_data  !!the external data state

    CHARACTER(len=*), PARAMETER :: routine =  'calc_o3_gems'

    INTEGER , PARAMETER :: ilat=64,nlev_gems=91

    ! Taken from su_ghgclim.F90 of ECMWF's IFS (37r2).
    REAL(wp), PARAMETER     :: zrefp(nlev_gems) = (/ &
      & 0.10000000E+01_wp, 0.29900000E+01_wp, 0.56834998E+01_wp, 0.10147500E+02_wp, &
      & 0.17160500E+02_wp, 0.27682501E+02_wp, 0.42848999E+02_wp, 0.63956501E+02_wp, &
      & 0.92440987E+02_wp, 0.12985049E+03_wp, 0.17781149E+03_wp, 0.23799651E+03_wp, &
      & 0.31209000E+03_wp, 0.40175449E+03_wp, 0.50860199E+03_wp, 0.63416602E+03_wp, &
      & 0.77987903E+03_wp, 0.94705548E+03_wp, 0.11368750E+04_wp, 0.13503740E+04_wp, &
      & 0.15884360E+04_wp, 0.18517920E+04_wp, 0.21410149E+04_wp, 0.24565259E+04_wp, &
      & 0.27986001E+04_wp, 0.31673640E+04_wp, 0.35628110E+04_wp, 0.39848059E+04_wp, &
      & 0.44330962E+04_wp, 0.49073169E+04_wp, 0.54070068E+04_wp, 0.59314971E+04_wp, &
      & 0.64797832E+04_wp, 0.70505981E+04_wp, 0.76428970E+04_wp, 0.82572246E+04_wp, &
      & 0.88957646E+04_wp, 0.95614326E+04_wp, 0.10257570E+05_wp, 0.10988080E+05_wp, &
      & 0.11757620E+05_wp, 0.12571580E+05_wp, 0.13435160E+05_wp, 0.14352070E+05_wp, &
      & 0.15325200E+05_wp, 0.16357520E+05_wp, 0.17452211E+05_wp, 0.18612539E+05_wp, &
      & 0.19841900E+05_wp, 0.21144000E+05_wp, 0.22522650E+05_wp, 0.23981760E+05_wp, &
      & 0.25525529E+05_wp, 0.27158301E+05_wp, 0.28884641E+05_wp, 0.30709301E+05_wp, &
      & 0.32637240E+05_wp, 0.34673641E+05_wp, 0.36823859E+05_wp, 0.39093570E+05_wp, &
      & 0.41487949E+05_wp, 0.44010520E+05_wp, 0.46651148E+05_wp, 0.49386160E+05_wp, &
      & 0.52190051E+05_wp, 0.55035488E+05_wp, 0.57894629E+05_wp, 0.60745699E+05_wp, &
      & 0.63574230E+05_wp, 0.66367703E+05_wp, 0.69114523E+05_wp, 0.71801031E+05_wp, &
      & 0.74412922E+05_wp, 0.76938641E+05_wp, 0.79368109E+05_wp, 0.81689688E+05_wp, &
      & 0.83891563E+05_wp, 0.85965023E+05_wp, 0.87903594E+05_wp, 0.89700203E+05_wp, &
      & 0.91347422E+05_wp, 0.92841422E+05_wp, 0.94181656E+05_wp, 0.95367977E+05_wp, &
      & 0.96399797E+05_wp, 0.97280203E+05_wp, 0.98034609E+05_wp, 0.98677672E+05_wp, &
      & 0.99198242E+05_wp, 0.99594977E+05_wp, 0.99881500E+05_wp /)

    ! Taken from su_ghgclim.F90 of ECMWF's IFS (37r2).
    REAL(wp), PARAMETER     :: zmday(12) = (/ &
      &  31._wp,    59.25_wp,  90.25_wp, 120.25_wp, 151.25_wp, 181.25_wp, &
      & 212.25_wp, 243.25_wp, 273.25_wp, 304.25_wp, 334.25_wp, 365.25_wp /)

    ! Taken from su_ghgclim.F90 of ECMWF's IFS (37r2).
    REAL(wp), PARAMETER     :: zytime(12) = (/ &
      &  22320._wp,  64980._wp, 107640._wp, 151560._wp, 195480._wp, 239400._wp, &
      & 283320._wp, 327960._wp, 371880._wp, 415800._wp, 459720._wp, 503640._wp /)    

    ! local fields
    INTEGER  :: idx0(nproma,0:pt_patch%nlev,pt_patch%nblks_c)
    REAL(wp) :: zlat(0:ilat+1)
    REAL(wp) :: zozn(0:ilat+1,1:nlev_gems)
    REAL(wp) :: zpresh(0:nlev_gems)
    REAL(wp) :: rclpr(0:nlev_gems)
    REAL(wp) :: zo3(nproma,1:nlev_gems,pt_patch%nblks_c)
    REAL(wp) :: zviozo(nproma,0:pt_patch%nlev)
    REAL(wp) :: zozovi(nproma,0:nlev_gems,pt_patch%nblks_c)
    REAL(wp) :: deltaz(nproma,pt_patch%nlev),dtdz(nproma,pt_patch%nlev),o3_clim(pt_patch%nlev)
    LOGICAL  :: l_found(nproma)

    ! local scalars
    INTEGER  :: jk,jkk,jk1,jl,jc,jb !loop indices
    INTEGER  :: idy,im,imn,im1,im2,jk_start,i_startidx,i_endidx,i_nchdom,i_startblk,i_endblk
    INTEGER  :: rl_start,rl_end,k375,k100,ktp
    REAL(wp) :: ztimi,zxtime,zjl,zlatint,zint,zadd_o3,tuneo3_1(nlev_gems),tuneo3_2(nlev_gems),&
                o3_macc1,o3_macc2,o3_gems1,o3_gems2
    REAL(wp) :: dzsum,dtdzavg,tpshp,wfac,wfac_lat(ilat),wfac_p(nlev_gems),wfac_tr(ilat),wfac_p_tr(nlev_gems)
    LOGICAL  :: lfound_all


    !Time index. Taken from su_ghgclim.F90 of ECMWF's IFS (37r2).

    IDY = mtime_datetime%date%day - 1 !NDD(KINDAT)-1
    IMN = mtime_datetime%date%month ! NMM(KINDAT)
    IF (IMN == 1) THEN
      ZXTIME=REAL(IDY*1440 + mtime_datetime%time%hour*60 + mtime_datetime%time%minute, wp) !KMINUT      
    ELSEIF (IMN == 2) THEN
      IF(IDY == 28) IDY=IDY-1
      ! A DAY IN FEB. IS 28.25*24*60/28=1452.8571min LONG.
      ZXTIME=44640._wp+REAL(IDY,wp)*1452.8571_wp + REAL(mtime_datetime%time%hour*60 + mtime_datetime%time%minute, wp) !KMINUT
    ELSE
      ZXTIME=(ZMDAY(IMN-1)+REAL(IDY,KIND(ZXTIME)))*1440._wp + REAL(mtime_datetime%time%hour*60 + mtime_datetime%time%minute, wp) !KMINUT
    ENDIF
    ! 525960=MINUTES IN A SIDERAL YEAR (365.25d)
    ZXTIME=MOD(ZXTIME,525960._wp)

    IM1=0
    IM2=0
    IF (ZXTIME <= ZYTIME(1)) THEN
      IM1=12
      IM2=1
      ZTIMI=(ZYTIME(1)-ZXTIME)/44640._wp
    ELSEIF(ZXTIME > ZYTIME(12)) THEN
      IM1=12
      IM2=1
      ZTIMI=(548280._wp-ZXTIME)/44640._wp
      ! 548280.=(365.25d + 15.5d)*24*60
    ELSE
      DO IM=1,11
        IF (ZXTIME > ZYTIME(IM) .AND. ZXTIME <= ZYTIME(IM+1)) THEN
          IM1=IM
          IM2=IM+1
          ZTIMI=(ZXTIME-ZYTIME(IM2))/(ZYTIME(IM1)-ZYTIME(IM2))
        ENDIF
      ENDDO
      IF ( IM1 == 0 .OR. IM2 == 0 ) THEN
        CALL finish(TRIM(routine),'Problem with time interpolation in suecozc!')
      ENDIF
    ENDIF

    ! Pressure levels of climatological ozone fields. From su_ghgclim.F90 of ECMWF's IFS (37r2).

    ZPRESH(0)=0.0_wp
    RCLPR(0) =0.0_wp
    DO JK=1,NLEV_GEMS-1
      ZPRESH(JK)=(ZREFP(JK)+ZREFP(JK+1))*0.5_wp
      RCLPR(JK) =ZPRESH(JK)
    ENDDO
    ZPRESH(NLEV_GEMS)=110000._wp
    RCLPR(nlev_gems) =ZPRESH(nlev_gems)

    ! Preparations for latitude interpolations

    zlatint=180._wp/REAL(ilat,wp)

    DO jl=0,ilat+1
      zlat(jl)=(-90._wp+0.5_wp*zlatint+(jl-1)*zlatint)*deg2rad
    ENDDO

    ! volume mixing ratio to ozone pressure thickness

    SELECT CASE (irad_o3)
    CASE (7)
      DO jk=1,nlev_gems
        DO jl=1,ilat
          zozn(JL,JK) = amo3/amd * (RGHG7(JL,JK,IM2)&
            & +ZTIMI*(RGHG7(JL,JK,IM1)-RGHG7(JL,JK,IM2)))
          zozn(JL,JK) = zozn(JL,JK) * (ZPRESH(JK)-ZPRESH(JK-1))
        ENDDO
      ENDDO
    CASE (9)
      DO jk=1,nlev_gems
        DO jl=1,ilat
          zozn(JL,JK) = amo3/amd * (RGHG7_MACC(JL,JK,IM2)&
            & +ZTIMI*(RGHG7_MACC(JL,JK,IM1)-RGHG7_MACC(JL,JK,IM2)))
          zozn(JL,JK) = zozn(JL,JK) * (ZPRESH(JK)-ZPRESH(JK-1))
        ENDDO
      ENDDO
    CASE (79) ! blending between GEMS and MACC

      ! Latitude-dependent weight for using MACC in the Antarctic region
      DO jl = 1, ilat
        IF (zlat(jl)*rad2deg > -45._wp) THEN
          wfac_lat(jl) = 0._wp
        ELSE IF (zlat(jl)*rad2deg > -60._wp) THEN
          wfac_lat(jl) = -(45._wp+zlat(jl)*rad2deg)/15._wp
        ELSE
          wfac_lat(jl) = 1._wp
        ENDIF
      ENDDO

      ! Pressure-dependent weight for using MACC in the upper stratosphere and mesosphere
      DO jk = 1, nlev_gems
        IF (zrefp(jk) > 500._wp) THEN
          wfac_p(jk) = 0._wp
        ELSE IF (zrefp(jk) > 100._wp) THEN
          wfac_p(jk) = 1._wp - (zrefp(jk)-100._wp)/400._wp
        ELSE
          wfac_p(jk) = 1._wp
        ENDIF
      ENDDO

      ! Latitude mask field for tropics (used for ozone enhancement in January and February)
      DO jl = 1, ilat
        IF (ABS(zlat(jl))*rad2deg > 30._wp) THEN
          wfac_tr(jl) = 0._wp
        ELSE IF (ABS(zlat(jl))*rad2deg > 20._wp) THEN
          wfac_tr(jl) = (30._wp-ABS(zlat(jl))*rad2deg)/10._wp
        ELSE
          wfac_tr(jl) = 1._wp
        ENDIF
      ENDDO

      ! Pressure mask field for tropics (used for ozone enhancement in January and February)
      DO jk = 1, nlev_gems
        IF (zrefp(jk) >= 7000._wp .AND. zrefp(jk) <= 10000._wp) THEN
          wfac_p_tr(jk) = 1._wp
        ELSE IF (zrefp(jk) < 7000._wp .AND. zrefp(jk) >= 5000._wp) THEN
          wfac_p_tr(jk) = (zrefp(jk)-5000._wp)/2000._wp
        ELSE IF (zrefp(jk) > 10000._wp .AND. zrefp(jk) < 15000._wp) THEN
          wfac_p_tr(jk) = (15000._wp-zrefp(jk))/5000._wp
        ELSE
          wfac_p_tr(jk) = 0._wp
        ENDIF
      ENDDO

      ! Profile functions for accelerated ozone hole filling in November
      ! (Accomplished by taking a weighted average between November and December climatologies)
      DO jk = 1, nlev_gems
        IF (zrefp(jk) >= 2000._wp .AND. zrefp(jk) <= 10000._wp) THEN
          tuneo3_1(jk) = 1._wp
        ELSE IF (zrefp(jk) < 2000._wp .AND. zrefp(jk) >= 1000._wp) THEN
          tuneo3_1(jk) = (zrefp(jk)-1000._wp)/1000._wp
        ELSE IF (zrefp(jk) > 10000._wp .AND. zrefp(jk) < 15000._wp) THEN
          tuneo3_1(jk) = (15000._wp-zrefp(jk))/5000._wp
        ELSE
          tuneo3_1(jk) = 0._wp
        ENDIF
      ENDDO
      DO jk = 1, nlev_gems
        IF (zrefp(jk) <= 2500._wp) THEN
          tuneo3_2(jk) = 1._wp
        ELSE IF (zrefp(jk) <= 12500._wp) THEN
          tuneo3_2(jk) = (12500._wp-zrefp(jk))/10000._wp
        ELSE
          tuneo3_2(jk) = 0._wp
        ENDIF
      ENDDO

      DO jk=1,nlev_gems
        DO jl=1,ilat
          wfac = MAX(wfac_lat(jl),wfac_p(jk))
          o3_macc1 = RGHG7_MACC(JL,JK,IM1) + MERGE(wfac_lat(jl)*tuneo3_1(jk)*tuneo3_2(jk)*&
                     MAX(0._wp,RGHG7_MACC(JL,JK,12)-RGHG7_MACC(JL,JK,IM1)), 0._wp, im1==11)
          o3_macc2 = RGHG7_MACC(JL,JK,IM2) + MERGE(wfac_lat(jl)*tuneo3_1(jk)*tuneo3_2(jk)*&
                     MAX(0._wp,RGHG7_MACC(JL,JK,12)-RGHG7_MACC(JL,JK,IM2)), 0._wp, im2==11)

          o3_gems1 = RGHG7(JL,JK,IM1) + MERGE(wfac_tr(jl)*wfac_p_tr(jk)*&
                     MAX(0._wp,RGHG7(JL,JK,12)-RGHG7(JL,JK,IM1)), 0._wp, im1==1 .OR. im1==2)
          o3_gems2 = RGHG7(JL,JK,IM2) + MERGE(wfac_tr(jl)*wfac_p_tr(jk)*&
                     MAX(0._wp,RGHG7(JL,JK,12)-RGHG7(JL,JK,IM2)), 0._wp, im1==1 .OR. im1==2)

          zozn(JL,JK) = amo3/amd * ( wfac * (o3_macc2+ZTIMI*(o3_macc1-o3_macc2)) + &
                              (1._wp-wfac)* (o3_gems2+ZTIMI*(o3_gems1-o3_gems2)) )
          zozn(JL,JK) = zozn(JL,JK) * (ZPRESH(JK)-ZPRESH(JK-1))
        ENDDO
      ENDDO

    CASE (io3_art)
      DO jk=1,nlev_gems
        DO jl=1,ilat
          zozn(JL,JK) = amo3/amd * (RGHG7(JL,JK,IM2)&
            & +ZTIMI*(RGHG7(JL,JK,IM1)-RGHG7(JL,JK,IM2)))
          zozn(JL,JK) = zozn(JL,JK) * (ZPRESH(JK)-ZPRESH(JK-1))
        ENDDO
      ENDDO
    END SELECT 

    DO jk=1,nlev_gems
      zozn(0,JK)      = zozn(1,jk)
      zozn(ilat+1,jk) = zozn(ilat,jk)
    ENDDO
    
    ! nest boudaries have to be included for reduced-grid option
    rl_start = 1
    rl_end   = min_rlcell_int

    i_nchdom  = MAX(1,pt_patch%n_childdom)

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,jkk,jk1,i_startidx,i_endidx,zjl,jk_start,l_found,lfound_all,&
!$OMP zint,zviozo,zadd_o3,deltaz,dtdz,dzsum,dtdzavg,ktp,tpshp,wfac,k375,k100,o3_clim) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(pt_patch, jb, i_startblk, i_endblk, &
        &                 i_startidx, i_endidx, rl_start, rl_end)

      ! Latitude interpolation
      
      DO jkk=1,nlev_gems

        DO jc=i_startidx,i_endidx

          zjl=1._wp+(rad2deg*pt_patch%cells%center(jc,jb)%lat+90._wp-.5_wp*zlatint)/zlatint

          !just select nearest value
          !zo3(jc,jkk,jb) = zozn(NINT(zjl),jkk)

          !linear interpolation
          zo3(jc,jkk,jb) = zozn(INT(zjl),jkk) &
            & + (zozn(INT(zjl)+1,jkk)-zozn(INT(zjl),jkk)) &
            &  /  (zlat(INT(zjl)) - zlat(INT(zjl)+1)) &
            &  * (zlat(INT(zjl)) - pt_patch%cells%center(jc,jb)%lat)

        ENDDO !jc

      ENDDO !jkk

      ! ACCUMULATE FROM TOP TO BOTTOM THE LATITUDE INTERPOLATED FIELDS
      ! From radghg.F90 of ECMWF's IFS.

      zozovi(i_startidx:i_endidx,0,jb) = 0._wp
      DO jkk=1,nlev_gems
        zozovi(i_startidx:i_endidx,jkk,jb) = zozovi(i_startidx:i_endidx,jkk-1,jb) &
          &                                 + zo3(i_startidx:i_endidx,jkk,jb)
      ENDDO

      ! REDISTRIBUTE THE VERTIC. INTEGR. CLIM. O3 ON THE MODEL GRID
      ! Adapted from radghg.F90 of ECMWF's IFS.

      jk = 0
      zviozo(i_startidx:i_endidx,0) = 0._wp

      jk_start = 0

      DO jk = 0,pt_patch%nlev
        l_found(:) = .FALSE.
        lfound_all = .FALSE.
        DO jkk = jk_start,nlev_gems-1
          DO jc = i_startidx,i_endidx
            IF( p_diag%pres_ifc(jc,jk+1,jb) >= RCLPR(jkk)  &
              & .AND. p_diag%pres_ifc(jc,jk+1,jb) < RCLPR(jkk+1)) THEN
              ZINT=(p_diag%pres_ifc(jc,jk+1,jb)-RCLPR(jkk))/(RCLPR(jkk+1)-RCLPR(jkk)) 
              ZVIOZO(jc,JK) = ZOZOVI(jc,jkk,jb) + ZINT * (ZOZOVI(jc,jkk+1,jb)-ZOZOVI(jc,jkk,jb))
              l_found(jc) = .TRUE.
              idx0(jc,jk,jb) = jkk
            ELSEIF ( p_diag%pres_ifc(jc,jk+1,jb) > RCLPR(nlev_gems) ) THEN
              l_found(jc) = .TRUE.
              idx0(jc,jk,jb) = nlev_gems
            ENDIF
          ENDDO !jc
          IF (ALL(l_found(i_startidx:i_endidx))) THEN
            lfound_all = .TRUE.
            EXIT
          ENDIF
        ENDDO !jkk
        IF (lfound_all) THEN
          jk_start = MIN(MINVAL(idx0(i_startidx:i_endidx,jk,jb)),nlev_gems-1)
        ENDIF
      ENDDO !jk

      ! COMPUTE THE MASS MIXING RATIO

      DO jk = 1,pt_patch%nlev
!DIR$ IVDEP
        DO jc = i_startidx,i_endidx
          ext_data%atm%o3(jc,jk,jb)=(ZVIOZO(jc,jk)-ZVIOZO(jc,jk-1)) / p_diag%dpres_mc(jc,jk,jb)

! Tuning to increase stratospheric ozone in order to reduce temperature biases;
! the tuning factors are computed in atm_phy_nwp_config
          IF ( ltuning_ozone ) THEN
            zadd_o3 = MIN(atm_phy_nwp_config(pt_patch%id)%ozone_maxinc,                   &
              ext_data%atm%o3(jc,jk,jb) * atm_phy_nwp_config(pt_patch%id)%fac_ozone(jk) * &
              atm_phy_nwp_config(pt_patch%id)%shapefunc_ozone(jc,jb) )
            ext_data%atm%o3(jc,jk,jb) = ext_data%atm%o3(jc,jk,jb) + zadd_o3
          ENDIF

        ENDDO
      ENDDO

      ! Ozone adaptation around the extratropical tropopause
      ! The procedure starts with diagnosing the thermal tropopause including its sharpness (dTdz_up - dTdz_down),
      ! the latter is used to compute a weighting factor. Afterwards, climatological O3 mixing ratios are modified
      ! in order to get a sharp jump at the tropopause, using the climatological values at 100 hPa (375 hPa) as
      ! proxies for the lower stratospheric (tropospheric) ozone values
      IF (icpl_o3_tp == 1) THEN
        DO jk = 1,pt_patch%nlev-1
          DO jc = i_startidx,i_endidx
            deltaz(jc,jk) = -rd/grav*(p_diag%temp(jc,jk,jb)+p_diag%temp(jc,jk+1,jb))*                       &
              (p_diag%pres(jc,jk,jb)-p_diag%pres(jc,jk+1,jb))/(p_diag%pres(jc,jk,jb)+p_diag%pres(jc,jk+1,jb))
            dtdz(jc,jk) = (p_diag%temp(jc,jk,jb)-p_diag%temp(jc,jk+1,jb))/deltaz(jc,jk)
          ENDDO
        ENDDO

        DO jc = i_startidx,i_endidx
          l_found(jc) = .FALSE.
          IF (ABS(pt_patch%cells%center(jc,jb)%lat)*rad2deg > 30._wp) THEN
            ! Determine thermal tropopause according to WMO definition
            DO jk = prm_diag%k850(jc,jb), 3, -1
              IF (p_diag%pres(jc,jk,jb) < 37500._wp .AND. p_diag%pres(jc,jk,jb) > 10000._wp) THEN
                IF (dtdz(jc,jk) < -2.e-3_wp .AND. dtdz(jc,jk-1) > -2.e-3_wp) THEN
                  l_found(jc) = .TRUE.
                  DO jk1 = jk-2, 2, -1
                    dzsum = SUM(deltaz(jc,jk1:jk-1))
                    IF (dzsum > 2.e3_wp) EXIT
                    dtdzavg = (p_diag%temp(jc,jk1,jb)-p_diag%temp(jc,jk,jb))/dzsum
                    IF (dtdzavg < -2.e-3_wp) l_found(jc) = .FALSE.
                  ENDDO
                ENDIF
              ENDIF
              IF (l_found(jc)) THEN
                ! weighting factor; the latitude-dependent component increases from 0 to 1 between 30 and 40 deg lat,
                ! the component depending on the tropopause sharpness increases from 0 to 1 between 2.5 K/km and 7.5 K/km
                ktp = jk
                tpshp = 0.5_wp*(dtdz(jc,jk-1)+dtdz(jc,jk-2)) - 0.5_wp*(dtdz(jc,jk)+dtdz(jc,jk+1))
                wfac = MIN(1._wp,0.1_wp*(ABS(pt_patch%cells%center(jc,jb)%lat)*rad2deg-30._wp)) * &
                  MIN(1._wp,200._wp*MAX(0._wp,tpshp-2.5e-3_wp))
                EXIT
              ENDIF
            ENDDO
          ENDIF
          IF (l_found(jc)) THEN
            ! Determine level indices right below 100 and 375 hPa
            k375 = prm_diag%k850(jc,jb)  ! to be safe over very high mountains
            DO jk = prm_diag%k850(jc,jb), 3, -1
              IF (p_diag%pres(jc,jk,jb) > 37500._wp .AND. p_diag%pres(jc,jk-1,jb) <= 37500._wp) k375 = jk
              IF (p_diag%pres(jc,jk,jb) > 10000._wp .AND. p_diag%pres(jc,jk-1,jb) <= 10000._wp) k100 = jk
              IF (p_diag%pres(jc,jk,jb) < 10000._wp) EXIT
            ENDDO
            o3_clim(k100:k375) = ext_data%atm%o3(jc,k100:k375,jb)
            jkk = k100+1
            DO jk = k100, k375
              ! Modify ozone profiles; the climatological profile is shifted down by at most 125 hPa
              IF (jk < ktp) THEN ! levels above the tropopause
                IF (p_diag%pres(jc,jk,jb) - p_diag%pres(jc,k100,jb) < 12500._wp) THEN
                  ext_data%atm%o3(jc,jk,jb) = (1._wp-wfac)*ext_data%atm%o3(jc,jk,jb) + wfac*o3_clim(k100)
                ELSE
                  DO jk1 = jkk, k375
                    IF (p_diag%pres(jc,jk,jb) - p_diag%pres(jc,jk1-1,jb) >= 12500._wp .AND. &
                        p_diag%pres(jc,jk,jb) - p_diag%pres(jc,jk1,jb)   < 12500._wp ) THEN
                      ext_data%atm%o3(jc,jk,jb) = (1._wp-wfac)*ext_data%atm%o3(jc,jk,jb) + wfac*o3_clim(jk1)
                      jkk = jk1
                      EXIT
                    ENDIF
                  ENDDO
                ENDIF
              ELSE IF (jk > ktp) THEN ! levels below the tropopause
                ext_data%atm%o3(jc,jk,jb) = (1._wp-wfac)*ext_data%atm%o3(jc,jk,jb) + wfac*o3_clim(k375)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF

    ENDDO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL     

  END SUBROUTINE calc_o3_gems
  
END MODULE mo_o3_util

