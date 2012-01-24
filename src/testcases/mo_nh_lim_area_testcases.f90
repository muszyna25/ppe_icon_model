!>
!!  Subroutine to initialized several test cases 
!!   for the NH-Core in limited area mode
!!
!!
!! @par Revision History
!! - first version by P. Ripodas , DWD, (2012-01)
!!
!! @par Literature
!! -M. L. Weisman and J. B. Klemp, 1982
!!  The Dependence of Numerically Simulated Convective Storms on 
!!  Vertical Wind Shear and Buoyancy.
!!  Monthly Weather Review, 110,504-520
!! - COSMO documentation
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
  MODULE mo_nh_lim_area_testcases
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2008
!
!-------------------------------------------------------------------------
!

   USE mo_kind,                ONLY: wp
   USE mo_physical_constants,  ONLY: rd_o_cpd, p0ref, grav, tmelt,  &
                                   & cvd_o_rd, re, omega, cpd ,     &
                                     vtmpc1 , rdv,  rd            
                                     
   USE mo_math_constants,      ONLY: pi, deg2rad, rad2deg
   USE mo_model_domain,        ONLY: t_patch
   USE mo_nonhydro_state,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
   USE mo_run_config,          ONLY: iqv,iqc, ntracer
   USE mo_impl_constants,      ONLY: inwp, MAX_CHAR_LENGTH, min_rlcell_int
   USE mo_parallel_config,     ONLY: nproma, p_test_run
   USE mo_satad,               ONLY:  sat_pres_water, &  !! saturation vapor pressure w.r.t. water
            &                         sat_pres_ice,   &  !! saturation vapor pressure w.r.t. ice
            &                         spec_humi          !! Specific humidity
   USE mo_exception,           ONLY: message, finish, message_text
   USE mo_advection_config,    ONLY: advection_config
   USE mo_interpolation,       ONLY: t_int_state, cells2edges_scalar, edges2cells_scalar, &
                                   & rbf_vec_interpol_cell
   USE mo_loopindices,         ONLY: get_indices_e
   USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
   USE mo_extpar_config,        ONLY: itopo
   USE mo_sync,                 ONLY: global_sum_array, sync_patch_array,  sync_patch_array_mult, &
                                   &  SYNC_C, SYNC_E
   USE mo_nh_init_utils,        ONLY: init_w, hydro_adjust, convert_thdvars, virtual_temp
   USE mo_vertical_coord_table, ONLY: vct_a

   IMPLICIT NONE

   PUBLIC  :: init_nh_atmo_ana_nconstlayers, init_nh_anaprof_uv

   PRIVATE

   REAL(wp), PARAMETER :: cpd_o_rd = 1._wp / rd_o_cpd
   REAL(wp), PARAMETER :: grav_o_cpd = grav / cpd
   REAL(wp), PARAMETER :: grav_o_rd = grav / rd


! !DEFINED PARAMETERS for the piecewise const. Brunt-Vaisala-freq (N) 
!  layers atmosphere:
   INTEGER, PARAMETER, PUBLIC  :: max_nlayers_nconst = 10
! !DEFINED PARAMETERS for the analytical wind profiles 
!  layers atmosphere:
   INTEGER, PARAMETER, PUBLIC  :: max_nlayers_linwind = 10

! !DEFINED namelist variables for the piecewise const. Brunt-Vaisala-freq (N) 
!  layers atmosphere:
   INTEGER, PUBLIC      :: nlayers_nconst      ! Number of desired layers with a constant N
   REAL(wp), PUBLIC     :: p_base_nconst       ! pressure at h_nconst(1)
   REAL(wp) , PUBLIC    :: theta0_base_nconst  ! potential temperature at h_nconst(1)
   REAL(wp), PUBLIC     :: h_nconst(max_nlayers_nconst) ! base height for each layer in m
   REAL(wp), PUBLIC     :: N_nconst(max_nlayers_nconst) ! N in 1/s for each layer
   REAL(wp), PUBLIC     :: rh_nconst(max_nlayers_nconst) ! base relative humidity for each layer
   REAL(wp), PUBLIC     :: rhgr_nconst(max_nlayers_nconst) ! gradien tof relative humidity for 
                                                           ! each layer in 1/m, positive for 
                                                           ! decreasing rel hum with height
! !DEFINED namelist variables for the analytical wind profile
  INTEGER, PUBLIC       :: itype_anaprof_uv    ! kind of wind profile
  ! For  itype_anaprof_uv == 1 (arbitrary number of constant gradient U(z) layers ):
  INTEGER, PUBLIC       :: nlayers_linwind     ! Number of layers
  REAL(wp), PUBLIC      :: h_linwind(max_nlayers_linwind) ! base height for each layer in m
  REAL(wp), PUBLIC      :: u_linwind(max_nlayers_linwind) ! U at base layer in m/s
  REAL(wp), PUBLIC      :: ugr_linwind(max_nlayers_linwind) ! gradient of U for each layer,
                                                            ! positive for increasing windspeed
                                                            ! with height, in 1/s
  ! For  itype_anaprof_uv == 2/3, constant U/V 
  REAl(wp),PUBLIC       :: vel_const                        ! 
     

!--------------------------------------------------------------------

  CONTAINS
!-------------------------------------------------------------------------
!
  !>
  !! Initialization of prognostic state for an 
  !! atmosphere with layers exhibiting constant 
  !! Brunt-Vaisala-Frequency N 
  !! (N takes water vapor relhum into account).
  !!
  !! It is assumed that N refers to a moist subsaturated atmosphere 
  !!
  !! @par Revision History
  !!
  !!
  SUBROUTINE init_nh_atmo_ana_nconstlayers( ptr_patch, ptr_nh_prog, ptr_nh_diag, &
    &                                p_metrics, p_int, l_hydro_adjust )


    TYPE(t_patch), TARGET, INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
      &  ptr_nh_prog

    TYPE(t_nh_diag), INTENT(INOUT)      :: &  !< diagnostic state vector
      &  ptr_nh_diag


    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state
    TYPE(t_int_state), INTENT(IN)       :: p_int
    LOGICAL, INTENT(IN)                 :: l_hydro_adjust !if .TRUE. hydrostatically balanced 
                                                         ! initial condition

    INTEGER        ::  jc, jb, jk,    &
                       nlen, nblks_e,  nblks_c, npromz_c, &
                       jl, jN
    INTEGER        :: nlev        !< number of full levels

    REAL(wp), ALLOCATABLE :: theta(:,:,:), theta_v(:,:,:),exner(:,:,:)
    REAL(wp), ALLOCATABLE :: relhum(:,:,:)
    REAL(wp), ALLOCATABLE :: z_qv(:,:,:)
    INTEGER,  ALLOCATABLE :: jNlayer(:,:,:)

    REAL(wp), DIMENSION(nlayers_nconst) :: thetab, exnerb, tempb, qvb, &
                             theta_vb, rhb, presb

  
    REAL(wp)  :: z_h, e_aux, temp_aux, pres_aux,  tempv_aux

!--------------------------------------------------------------------
!


    ! number of vertical levels
    nlev   = ptr_patch%nlev
    ALLOCATE (theta(nproma,nlev,ptr_patch%nblks_c), &
              theta_v(nproma,nlev,ptr_patch%nblks_c), &
              exner(nproma,nlev,ptr_patch%nblks_c), &
              relhum(nproma,nlev,ptr_patch%nblks_c) )
    ALLOCATE ( z_qv(nproma,nlev,ptr_patch%nblks_c) )
    ALLOCATE ( jNlayer(nproma,nlev,ptr_patch%nblks_c) )
   nblks_c   = ptr_patch%nblks_int_c
   npromz_c  = ptr_patch%npromz_int_c
   nblks_e   = ptr_patch%nblks_int_e
   !npromz_e  = ptr_patch%npromz_int_e

! set the values at h_nconst(1)
    thetab(1) = theta0_base_nconst
    exnerb(1) = (p_base_nconst/p0ref)**rd_o_cpd
    rhb(1)    = rh_nconst(1)
    presb(1)  = p_base_nconst
    tempb(1)  =  thetab(1) * exnerb(1)
    qvb(1)      = qv_rhtp(rhb(1), tempb(1), presb(1) )
    theta_vb(1) = thetab(1)*(1._wp+vtmpc1*qvb(1)) 
   
! set the values at the base of each layer
! for the base at one layer, integrate the hydrostatic equation in the layer below
  WRITE(*,*) nlayers_nconst
  DO  jl = 2, nlayers_nconst
    thetab(jl) = thetab(jl-1)*           &
               & EXP( N_nconst(jl-1)**2*(h_nconst(jl)-h_nconst(jl-1))/grav )
    !here I use rh of the layer below
    rhb(jl)    = rh_nconst(jl-1)-rhgr_nconst(jl-1)*(h_nconst(jl)-h_nconst(jl-1))
    rhb(jl)    = MAX(0._wp,MIN(1._wp,rhb(jl)))

!first consider constant qv (= qvb) in the layer jl-1

    exnerb(jl) = exner_nconst(exnerb(jl-1), N_nconst(jl-1),thetab(jl-1), & 
                           &  thetab(jl), qvb(jl-1), qvb(jl-1))
    presb(jl)  = p0ref*(exnerb(jl)**cpd_o_rd)
    tempb(jl)  = thetab(jl)*exnerb(jl)
    qvb(jl)    = qv_rhtp(rhb(jl), tempb(jl), presb(jl) )
! recalculate exner with the fisrt estimation of qvb, then recalculate qvb
    exnerb(jl) = exner_nconst(exnerb(jl-1), N_nconst(jl-1),thetab(jl-1), & 
                           &  thetab(jl), qvb(jl-1), qvb(jl))
    presb(jl)  = p0ref*(exnerb(jl)**cpd_o_rd)
    tempb(jl)  = thetab(jl)*exnerb(jl)
    qvb(jl)    = qv_rhtp(rhb(jl), tempb(jl), presb(jl) )
! last estimation of exner
    exnerb(jl) = exner_nconst(exnerb(jl-1), N_nconst(jl-1),thetab(jl-1), & 
                           &  thetab(jl), qvb(jl-1), qvb(jl))
    presb(jl)  = p0ref*(exnerb(jl)**cpd_o_rd)
    tempb(jl)  = thetab(jl)*exnerb(jl)
! now consider the relative humidity of the base of the current layer
    rhb(jl)    = rh_nconst(jl)
    qvb(jl)    = qv_rhtp(rhb(jl), tempb(jl), presb(jl) )
    theta_vb(jl)= thetab(jl)*(1._wp+vtmpc1*qvb(jl))

  END DO

  DO  jl = 1, nlayers_nconst
    WRITE(*,*) jl, presb(jl), tempb(jl), rhb(jl), qvb(jl), exnerb(jl), thetab(jl), theta_vb(jl)
  END DO
  
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jk,jc,nlen,z_h,jN,tempv_aux, pres_aux,temp_aux)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = nlev, 1, -1
         DO jc = 1, nlen

            z_h = p_metrics%z_mc(jc,jk,jb)
             !set the layer corresponding to this point
             IF (z_h <  h_nconst(1) ) THEN
               jNlayer(jc,jk,jb)=0
             ELSEIF (z_h >=  h_nconst(nlayers_nconst) ) THEN 
              jNlayer(jc,jk,jb)= nlayers_nconst
             ELSE
              DO jl=1,nlayers_nconst-1
               IF (z_h >= h_nconst(jl) .AND. z_h < h_nconst(jl+1)  ) THEN
                jNlayer(jc,jk,jb)=jl 
                !EXIT
               END IF          
              END DO
             END IF
             jN=jNlayer(jc,jk,jb)
             !WRITE(*,*) jc, jk, jb, jN, z_h, h_nconst(jN)
             IF (jN < 0 .OR. jN > nlayers_nconst) THEN
                CALL finish ('corresponding layer has not been found')
             END IF

             IF (jN == 0) THEN
              !consider a isothermal atmosphere below h_nconst(1) with qv=qvb(1)
              z_qv(jc,jk,jb) = qvb(1)
              tempv_aux = tempb(1)*(1._wp+vtmpc1*qvb(1))
              pres_aux  = presb(1) * &
                & EXP(-grav_o_rd*(z_h-h_nconst(1))/tempv_aux)
              exner(jc,jk,jb)= (pres_aux/p0ref)**rd_o_cpd  
              theta(jc,jk,jb)= tempb(1)/exner(jc,jk,jb)
              theta_v(jc,jk,jb)=theta(jc,jk,jb)*(1._wp+vtmpc1*qvb(1))
             ELSE
              theta(jc,jk,jb)=thetab(jN)*          &
               & EXP( N_nconst(jN)**2*(z_h-h_nconst(jN))/grav )
              relhum(jc,jk,jb)=rh_nconst(jN)-rhgr_nconst(jN)*(z_h-h_nconst(jN))
              relhum(jc,jk,jb)=MAX(0._wp,MIN(1._wp,relhum(jc,jk,jb)))
              !WRITE(*,*) theta(jc,jk,jb), relhum(jc,jk,jb), jN, jNlayer(jc,jk+1,jb)
              IF (jk < nlev .AND. jN == jNlayer(jc,jk+1,jb) ) THEN
               ! in this case we integrate starting in the level bellow, instead of 
               !  starting in the base of the layer

               ! first consider qv constant (= z_qv(jc, jk+1,jb))
               !WRITE(*,*) "from the previous level"
               exner(jc,jk,jb)= exner_nconst(exner(jc,jk+1,jb), N_nconst(jN),  &
                           &  theta(jc,jk+1,jb), theta(jc,jk,jb),              &
                           &  z_qv(jc, jk+1,jb), z_qv(jc, jk+1,jb))
               pres_aux  = p0ref*(exner(jc,jk,jb)**cpd_o_rd)
               temp_aux  = theta(jc,jk,jb)*exner(jc,jk,jb)
               ! 1st estimation of  qv, then re-estimate exner
               z_qv(jc,jk,jb)= qv_rhtp(relhum(jc,jk,jb), temp_aux, pres_aux )
               exner(jc,jk,jb)= exner_nconst(exner(jc,jk+1,jb), N_nconst(jN),  &
                           &  theta(jc,jk+1,jb), theta(jc,jk,jb),              &
                           &  z_qv(jc, jk+1,jb), z_qv(jc, jk,jb))
               pres_aux  = p0ref*(exner(jc,jk,jb)**cpd_o_rd)
               temp_aux  = theta(jc,jk,jb)*exner(jc,jk,jb)
               ! re-estimate qv, then re-estimate exner
               z_qv(jc,jk,jb)= qv_rhtp(relhum(jc,jk,jb), temp_aux, pres_aux )
               exner(jc,jk,jb)= exner_nconst(exner(jc,jk+1,jb), N_nconst(jN),  &
                           &  theta(jc,jk+1,jb), theta(jc,jk,jb),              &
                           &  z_qv(jc, jk+1,jb), z_qv(jc, jk,jb))
               temp_aux  = theta(jc,jk,jb)*exner(jc,jk,jb)
               ! re-estimate qv
               z_qv(jc,jk,jb)= qv_rhtp(relhum(jc,jk,jb), temp_aux, pres_aux )
               theta_v(jc,jk,jb)=theta(jc,jk,jb)*(1._wp+vtmpc1*z_qv(jc,jk,jb))

               
              ELSE
               !WRITE(*,*) "from the layer base"
               ! here we integrate from the base of the layer
               !
               ! here I could re-estimate the values at the base of the layer
               ! integrating from (jc,jk+1,jb) to the base of the layer 
               ! in the case jN>1
                IF (jN>1  .AND. exnerb(jN) > exner(jc,jk+1,jb)) THEN
                  ! check that exner(jN)<exner(jc,jk+1,jb)
                  CALL finish ('base layer has larger pressure than the level     &
                           &   below, exnerb(jN) was not well approximated,      &
                           &   it should be recalculated')
                             
                END IF
               !
               ! first consider qv constant (= qvb(jN))
               exner(jc,jk,jb)= exner_nconst(exnerb(jN), N_nconst(jN),thetab(jN), & 
                           &  theta(jc,jk,jb), qvb(jN), qvb(jN))
               pres_aux  = p0ref*(exner(jc,jk,jb)**cpd_o_rd)
               temp_aux  = theta(jc,jk,jb)*exner(jc,jk,jb)
               ! 1st estimation of  qv, then re-estimate exner
               z_qv(jc,jk,jb)= qv_rhtp(relhum(jc,jk,jb), temp_aux, pres_aux )
               exner(jc,jk,jb)= exner_nconst(exnerb(jN), N_nconst(jN),thetab(jN), & 
                           &  theta(jc,jk,jb), qvb(jN),z_qv(jc,jk,jb) )
               pres_aux  = p0ref*(exner(jc,jk,jb)**cpd_o_rd)
               temp_aux  = theta(jc,jk,jb)*exner(jc,jk,jb)
               ! re-estimate qv, then re-estimate exner
               z_qv(jc,jk,jb)= qv_rhtp(relhum(jc,jk,jb), temp_aux, pres_aux )
               exner(jc,jk,jb)= exner_nconst(exnerb(jN), N_nconst(jN),thetab(jN), & 
                           &  theta(jc,jk,jb), qvb(jN),z_qv(jc,jk,jb) )
               temp_aux  = theta(jc,jk,jb)*exner(jc,jk,jb)
               ! re-estimate qv
               z_qv(jc,jk,jb)= qv_rhtp(relhum(jc,jk,jb), temp_aux, pres_aux )
               theta_v(jc,jk,jb)=theta(jc,jk,jb)*(1._wp+vtmpc1*z_qv(jc,jk,jb))
               
              END IF
             END IF
             
             !WRITE(*,*) pres_aux, temp_aux, z_qv(jc,jk,jb), exner(jc,jk,jb), theta_v(jc,jk,jb)

         END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL


   ! Copy to prognostic model fields
!$OMP PARALLEL 
!$OMP DO PRIVATE(jb,jk,jc,nlen)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = nlev, 1, -1
        DO jc = 1, nlen  
          ptr_nh_prog%theta_v(jc,jk,jb)    = theta_v(jc,jk,jb)
          ptr_nh_prog%exner(jc,jk,jb)      = exner(jc,jk,jb)
          ptr_nh_prog%tracer(jc,jk,jb,iqv) = z_qv(jc,jk,jb)

          ptr_nh_prog%rho(jc,jk,jb)  = ptr_nh_prog%exner(jc,jk,jb)**cvd_o_rd*p0ref &
                                       /rd/ptr_nh_prog%theta_v(jc,jk,jb)
          ptr_nh_prog%rhotheta_v(jc,jk,jb) = ptr_nh_prog%rho(jc,jk,jb) *          &
                                             ptr_nh_prog%theta_v(jc,jk,jb)
        ENDDO !jc
      ENDDO !jk     
    ENDDO !jb
!$OMP END DO
!$OMP END PARALLEL

  CALL diagnose_pres_temp (p_metrics, ptr_nh_prog,ptr_nh_prog, ptr_nh_diag,     &
                             ptr_patch, opt_calc_pres=.TRUE., opt_calc_temp=.TRUE.)


  IF (l_hydro_adjust) THEN

   CALL hydro_adjust ( ptr_patch, p_metrics, ptr_nh_prog%rho,     &
                     & ptr_nh_prog%exner, ptr_nh_prog%theta_v,    &
                     & ptr_nh_prog%rhotheta_v  )

  END IF   

    DEALLOCATE(theta, theta_v, relhum)
    DEALLOCATE(exner, z_qv)
    DEALLOCATE(jNlayer)
  END SUBROUTINE init_nh_atmo_ana_nconstlayers

!-------------------------------------------------------------------------
!
  !>
  !! Initialization of the wind profile for the 
  !!  limited area testcases 
  !! 
  !!
  !! @par Revision History
  !!
  !!
  SUBROUTINE init_nh_anaprof_uv( ptr_patch, vn, w,    &
    &                                p_metrics, p_int )
    TYPE(t_patch), TARGET, INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    REAL(wp), INTENT(INOUT) :: vn(:,:,:)    ! edge-normal wind component (m/s)

    REAL(wp), INTENT(INOUT) :: w(:,:,:)  ! vertical wind component (m/s)
    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state
    TYPE(t_int_state), INTENT(IN)       :: p_int

    REAL(wp), ALLOCATABLE :: z_me(:,:,:)

    REAL(wp)              :: z_u, z_v
    INTEGER               :: je,jb,jk,jl,jN
    INTEGER               :: i_startblk, i_startidx, i_endidx, nblks_e 
    INTEGER               :: nlev

!--------------------------------------------------------------------
!
    ! number of vertical levels
    nlev   = ptr_patch%nlev

!    nblks_c   = ptr_patch%nblks_int_c
!    npromz_c  = ptr_patch%npromz_int_c
    nblks_e   = ptr_patch%nblks_int_e

    SELECT CASE (itype_anaprof_uv)

      CASE(1)
      ! arbitrary number of constant gradient U(z) layers 
      ! horizontal normal components of the velocity
      ! initialize horizontal velocities

       ALLOCATE(z_me(nproma,nlev,ptr_patch%nblks_e))
       z_me(:,:,:) = 0.0_wp

       ! Compute geometric height at edge points
       CALL cells2edges_scalar(p_metrics%z_mc, ptr_patch, &
               p_int%c_lin_e, z_me)

       CALL sync_patch_array(SYNC_E,ptr_patch,z_me)
       i_startblk = ptr_patch%edges%start_blk(2,1)

       ! horizontal normal components of the velocity
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je,z_u,jl,jN)
       DO jb = i_startblk, nblks_e

        CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx

             !set the layer corresponding to this point
             IF (z_me(je,jk,jb) <  h_linwind(1) ) THEN
               jN=0
             ELSEIF (z_me(je,jk,jb) >=  h_linwind(nlayers_linwind) ) THEN 
              jN= nlayers_linwind
             ELSE
              DO jl=1,nlayers_linwind-1
               IF (z_me(je,jk,jb) >= h_linwind(jl) .AND. z_me(je,jk,jb) < h_linwind(jl+1)) THEN
                jN=jl 
                !EXIT
               END IF          
              END DO
             END IF
             IF (jN <= 0 .OR. jN > nlayers_linwind) THEN
                CALL finish ('corresponding layer has not been found')
             END IF
            
            z_u = u_linwind(jN) +  ugr_linwind(jN)*    &
                                &  (z_me(je,jk,jb)-h_linwind(jN))   !v component is zero
            vn(je,jk,jb) = &
             z_u * ptr_patch%edges%primal_normal(je,jb)%v1
          ENDDO !je
        ENDDO !jk
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL
       DEALLOCATE(z_me)

      CASE(2)
      ! constant zonal wind
       z_u = vel_const
       i_startblk = ptr_patch%edges%start_blk(2,1)

       ! horizontal normal components of the velocity
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je)
       DO jb = i_startblk, nblks_e

        CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx

           vn(je,jk,jb) =  &
             z_u * ptr_patch%edges%primal_normal(je,jb)%v1
          ENDDO !je
        ENDDO !jk
     ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

      CASE(3)
      ! constant meridional wind
       z_v = vel_const
       i_startblk = ptr_patch%edges%start_blk(2,1)

       ! horizontal normal components of the velocity
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,je)
       DO jb = i_startblk, nblks_e

        CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

        DO jk = 1, nlev
          DO je = i_startidx, i_endidx

           vn(je,jk,jb) =  &
             z_v * ptr_patch%edges%primal_normal(je,jb)%v2
          ENDDO !je
        ENDDO !jk
     ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

    END SELECT

 ! initialize vertical velocity
   CALL init_w(ptr_patch, p_int, vn, p_metrics%z_ifc, w)
   CALL sync_patch_array(SYNC_C, ptr_patch, w)

  END SUBROUTINE init_nh_anaprof_uv


!--------------------------------------------------------------------

! Calculate the exner function integrating the hydrostatic equation.
!  It considers we are within a layer of constant Brunt Vaisala frequency.
!  It considers that qv is constant within the two levels of the integration, 
!  in practice it uses the mean value of qv
  REAL(wp)  FUNCTION exner_nconst(exnerb,N,thetab, theta, qvb, qv)
 
   REAL(wp), INTENT (IN):: exnerb, N, thetab, theta, qvb, qv !exnerb, thetab and qvb are 
                                                             ! the values at one level
                                                             ! N is the Brunt-Vaisala-freq
                                                             ! theta and qv are the values at 
                                                             ! the level for which we want to calculate
                                                             ! the exner function
   REAL(wp)             :: qv_mean, factor1, factor2

   qv_mean = 0.5_wp * (qvb+qv)
   factor1  = (1._wp + vtmpc1*qv_mean )
   factor2  = (grav/N)**2/(cpd*factor1)
   exner_nconst = exnerb+ factor2*(1._wp/theta - 1._wp/thetab)

  END FUNCTION exner_nconst
!--------------------------------------------------------------------
!
! Calculate specific humidity from relative humidity, temp and pres
  REAL(wp)  FUNCTION qv_rhtp(rh, temp, pres)
   REAL(wp), INTENT (IN)   :: rh, temp, pres
   REAL (wp)               :: e_aux

      IF (temp > tmelt) THEN
        e_aux     = rh*sat_pres_water(temp)
      ELSE
        e_aux     = rh*sat_pres_ice(temp)
      ENDIF
      qv_rhtp     = MAX(0._wp,MIN(1._wp,spec_humi(e_aux,pres)))
   END FUNCTION qv_rhtp

!--------------------------------------------------------------------
  END MODULE mo_nh_lim_area_testcases
