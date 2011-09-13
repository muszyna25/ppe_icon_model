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
MODULE mo_o3_util

  USE mo_get_utc_date_tr,      ONLY: get_utc_date_tr
  USE mo_parallel_config,      ONLY: nproma
  USE mo_impl_constants,       ONLY: min_rlcell
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_kind,                 ONLY: wp
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_math_constants,       ONLY: pi
  USE mo_model_domain,         ONLY: t_patch
  USE mo_nh_vert_interp,       ONLY: prepare_lin_intp, lin_intp

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: calc_o3_clim, o3_timeint, o3_pl2ml, o3_zl2ml

!  CHARACTER(len=*), PARAMETER :: version = '$Id$'

CONTAINS

  !=======================================================================

SUBROUTINE o3_timeint( kproma,kbdim, nlev_pres,nmonths,&  ! IN
                     & selmon,                          & ! IN
                     & ext_O3 ,                         & ! IN kproma,nlev_p,jb,nmonth
                     & o3_time_int                       )! OUT kproma,nlev_p

  ! NOTE the very first approach ist no time avareaging but selection of 
  ! a special month (september) for aqua planet simulations

    INTEGER, INTENT(in)         :: kproma ! 
    INTEGER, INTENT(in)         :: kbdim  ! 
    INTEGER, INTENT(in)         :: nmonths     ! number of months in o3
    INTEGER, INTENT(in),OPTIONAL:: selmon      ! selected month for experiment
    INTEGER, INTENT(in)         :: nlev_pres   ! number of o3 data levels

    REAL(wp), INTENT(in) , DIMENSION(kbdim,nlev_pres,nmonths):: ext_o3
    REAL(wp), INTENT(out), DIMENSION(kbdim,nlev_pres)        :: o3_time_int

    IF(PRESENT(selmon))THEN
      o3_time_int(1:kproma,1:nlev_pres) =ext_o3(1:kproma,1:nlev_pres,selmon)
    ELSE
      o3_time_int(1:kproma,1:nlev_pres) =ext_o3(1:kproma,1:nlev_pres,9) ! September
    ENDIF

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
    REAL(wp), DIMENSION(nproma,nlev   ,nblks),INTENT(OUT):: model_o3 !converted o3 field for the model

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

  SUBROUTINE calc_o3_clim(kbdim,p_inc_rad,z_sim_time,pt_patch,zvio3,zhmo3)

    INTEGER, INTENT(in)   :: &
      & kbdim

    REAL(wp), INTENT(in)   :: &
      & p_inc_rad, & !radiation time step in seconds
      & z_sim_time

    TYPE(t_patch),      INTENT(in)    :: pt_patch    ! Patch

    REAL(wp),      INTENT(inout)  :: &
      & zvio3(kbdim,pt_patch%nblks_c) , &
      & zhmo3(kbdim,pt_patch%nblks_c)

    REAL(wp) ::                     &
      & z_sim_time_rad,  &
      & zstunde,                   & ! output from routine get_utc_date_tr
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
      & jj, itaja, jb !& ! output from routine get_utc_date_tr

    INTEGER , SAVE :: itaja_o3_previous = 0
    
    INTEGER :: jmm,mmm,mnc,mns,jnn,jc

    i_nchdom  = MAX(1,pt_patch%n_childdom)
    
    z_sim_time_rad = z_sim_time + 0.5_wp*p_inc_rad

    CALL get_utc_date_tr (                 &
      &   p_sim_time     = z_sim_time_rad, &
      &   itype_calendar = 0,              &
      &   iyear          = jj,             &
      &   nactday        = itaja,          &
      &   acthour        = zstunde )

    !decide whether new ozone calculation is necessary
    IF ( itaja == itaja_o3_previous ) RETURN
    itaja_o3_previous = itaja

    
    ztwo    = 0.681_wp + 0.2422_wp*REAL(jj-1949,wp)-REAL((jj-1949)/4,wp)
    ztho    = 2._wp*pi*( REAL(itaja, wp) -1.0_wp + ztwo )/365.2422_wp


!
!     ozone parameters
!     ----------------

!     a) time of the year dependency
    CALL o3_par_t5 (ztho,zvio3_c,zvio3_s,zhmo3_c,zhmo3_s)

    !     b) inverse Legendre transform

      !in order to account for mesh refinement
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell

    i_startblk = pt_patch%cells%start_blk(rl_start,1)
    i_endblk   = pt_patch%cells%end_blk(rl_end,i_nchdom)
      
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, zalp, jmm, jc, zvio3_f ,zhmo3_f, mmm, mnc, mns, jnn ) 
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
!$OMP END DO
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


  !! This a copy from legtri_vec.f90 of DWD's GME (contentually identically, only formally
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

END MODULE mo_o3_util

