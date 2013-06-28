!>
!!  Subroutine to initialize the CBL case for HDCP2
!!
!!
!! @par Revision History
!! - first version by Anurag Dipankar , MPIM, (2012-12-12)
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
MODULE mo_nh_torus_exp
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2008
!
!-------------------------------------------------------------------------
!
!
!

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, min_rledge, min_rlcell, SUCCESS
  USE mo_io_units,            ONLY: find_next_free_unit
  USE mo_physical_constants,  ONLY: rd, rv, cpd, p0ref, cvd_o_rd, rd_o_cpd, grav, alv
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_math_constants,      ONLY: pi, pi2
  USE mo_math_utilities,      ONLY: gnomonic_proj, t_geographical_coordinates, &
     &                              t_cartesian_coordinates, gc2cc, az_eqdist_proj
  USE mo_loopindices,         ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics, t_nh_ref
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_math_utilities,      ONLY: plane_torus_distance, arc_length
  USE mo_grid_config,         ONLY: grid_sphere_radius, is_plane_torus
  USE mo_sync,                ONLY: global_sum_array, sync_patch_array, SYNC_C, SYNC_E
  USE mo_nh_init_utils,       ONLY: init_w
  USE mo_run_config,          ONLY: iqv, iqc
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_e, grf_bdywidth_c
  USE mo_satad,               ONLY: spec_humi, sat_pres_water
  USE mo_les_config,          ONLY: les_config
  USE mo_nh_diagnose_pres_temp,ONLY: diagnose_pres_temp
  USE mo_vert_utilities,      ONLY: vert_intp_linear_1d

  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: init_nh_state_cbl, cbl_stevens_fluxes, init_nh_state_rico, u_cbl, v_cbl, th_cbl
  PUBLIC :: init_nh_state_gate

  !Linear profiles of variables for LES testcases
  REAL(wp) :: u_cbl(2)   !u_cbl(1) = constant, u_cbl(2) = gradient
  REAL(wp) :: v_cbl(2)   !v_cbl(1) = constant, v_cbl(2) = gradient
  REAL(wp) :: th_cbl(2)  !th_cbl(1) = constant,th_cbl(2) = gradient

  !DEFINED PARAMETERS (Stevens 2007 JAS) for init_nh_state_cbl:
  REAL(wp), PARAMETER :: zp0     = p0ref      !< surface pressure
  REAL(wp), PARAMETER :: zh0     = 0._wp      !< height (m) above which temperature increases
  REAL(wp), PARAMETER :: lambda  = 1500._wp   !moist height from Stevens(2007)
  REAL(wp), PARAMETER :: dtdz_st = 0.04_wp    !< theta lapse rate in stratosphere (T>0!)
  REAL(wp), PARAMETER :: z_tropo = 11000._wp  !height tropopause
  REAL(wp), PARAMETER :: rh_sfc  = 0.71_wp    !RH at surface [1], default 0.8

  !DEFINED PARAMETERS (RICO case):
  REAL(wp), PARAMETER :: zpsfc   = 101540._wp   !< surface pressure
  REAL(wp), PARAMETER :: zh1     = 740._wp      !< height (m) above which temperature increases
  REAL(wp), PARAMETER :: ztsfc   = 297.9_wp
  REAL(wp), PARAMETER :: zh2     = 3260._wp     !< moist height for RICO
  REAL(wp), PARAMETER :: zh3     = 17500._wp    !< height for extrapolated RICO profiles
  REAL(wp), PARAMETER :: zh4     = 30000._wp    !< height for extrapolated RICO profiles
  REAL(wp), PARAMETER :: zh5     = 60000._wp    !< height for extrapolated RICO profiles

  !DEFINED PARAMETERS (GATE case):
  REAL(wp), PARAMETER :: zg_tropo = 17000._wp  !height tropopause
  REAL(wp), PARAMETER :: zp_sfc   = 101200._wp !< surface pressure
  REAL(wp), PARAMETER :: zt_sfc   = 299.88_wp

!--------------------------------------------------------------------

   CONTAINS
!-------------------------------------------------------------------------
!
!
  !>
  !! Initialization of prognostic state vector for the nh CBL test case 
  !!  without moisture
  !!
  !! @par Revision History
  !!
  !!
  SUBROUTINE init_nh_state_cbl( ptr_patch, ptr_nh_prog,  ptr_nh_ref, ptr_nh_diag,  &
    &                           ptr_int, ptr_metrics)

    ! INPUT PARAMETERS:
    TYPE(t_patch),TARGET,  INTENT(IN)   :: &  !< patch on which computation is performed
      &  ptr_patch
    TYPE(t_int_state),     INTENT(IN)   :: &
      &  ptr_int
    TYPE(t_nh_prog),       INTENT(INOUT):: &  !< prognostic state vector
      &  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT):: &  !< diagnostic state vector
      &  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   :: &
      &  ptr_metrics                          !< NH metrics state
    TYPE(t_nh_ref),        INTENT(INOUT):: &  !< reference state vector
      &  ptr_nh_ref

    REAL(wp) :: rho_sfc, z_help(1:nproma), zvn1, zvn2, zu, zv, zt00, zh00
    INTEGER  :: jc,jk,jb,i_startblk,i_startidx,i_endidx   !< loop indices
    INTEGER  :: nblks_c,npromz_c,nblks_e,npromz_e
    INTEGER  :: nlev, nlevp1                  !< number of full and half levels
    INTEGER  :: nlen, jcn, jbn, jg, ntropo

  !-------------------------------------------------------------------------

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_int_c
    npromz_c = ptr_patch%npromz_int_c
    nblks_e  = ptr_patch%nblks_int_e
    npromz_e = ptr_patch%npromz_int_e

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    nlevp1 = ptr_patch%nlevp1

    !patch id
    jg = ptr_patch%id

    !Set some reference density
    IF(les_config(jg)%isrfc_type==1)THEN
      rho_sfc = zp0 / (rd * les_config(jg)%sst)
    ELSE
      rho_sfc = zp0 / (rd * th_cbl(1) )
    END IF

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

      !Tracers
      IF(.NOT.les_config(jg)%is_dry_cbl)THEN
        DO jk = 1, nlev
          ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = rh_sfc * spec_humi(sat_pres_water(th_cbl(1)),zp0) * &
                    EXP(-ptr_metrics%z_mc(1:nlen,jk,jb)/lambda)
        END DO
      END IF

      ntropo = 0
      DO jk = nlev, 1, -1
         ! init potential temperature
         z_help(1:nlen) = th_cbl(1) + max(0._wp, (ptr_metrics%z_mc(1:nlen,jk,jb)-zh0)*th_cbl(2))

         ! constant temperature above tropopause
         if ((ptr_metrics%z_mc(1,jk,jb) > z_tropo) .and. (ntropo == 0)) then
            ntropo = 1
            zt00   = z_help(1)
            zh00   = ptr_metrics%z_mc(1,jk,jb)
         endif
         if (ptr_metrics%z_mc(1,jk,jb) > z_tropo) then
            z_help(1:nlen) = zt00 + (ptr_metrics%z_mc(1:nlen,jk,jb)-zh00) * dtdz_st
         endif

         ! virtual potential temperature
         ptr_nh_prog%theta_v(1:nlen,jk,jb) = z_help(1:nlen) * ( 1._wp + &
           0.61_wp*ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) - ptr_nh_prog%tracer(1:nlen,jk,jb,iqc) ) 
      END DO

      !Get hydrostatic pressure and exner at lowest level
      ptr_nh_diag%pres(1:nlen,nlev,jb)  = zp0 - rho_sfc * ptr_metrics%geopot(1:nlen,nlev,jb)
      ptr_nh_prog%exner(1:nlen,nlev,jb) = (ptr_nh_diag%pres(1:nlen,nlev,jb)/p0ref)**rd_o_cpd 

      !Get exner at other levels
      DO jk = nlev-1, 1, -1
         z_help(1:nlen) = 0.5_wp * ( ptr_nh_prog%theta_v(1:nlen,jk,jb) +  &
                                     ptr_nh_prog%theta_v(1:nlen,jk+1,jb) )
   
         ptr_nh_prog%exner(1:nlen,jk,jb) = ptr_nh_prog%exner(1:nlen,jk+1,jb) &
            &  -grav/cpd*ptr_metrics%ddqz_z_half(1:nlen,jk+1,jb)/z_help(1:nlen)
      END DO

      IF ( jb == 1 ) THEN
        DO jk = 1,nlev
          write(*,*) 'CBL case setup: level, p, T, theta,v, qv: ', jk,   &
            & ptr_nh_prog%exner(1,jk,jb)**(cpd/rd) * p0ref,              &
            & ptr_nh_prog%exner(1,jk,jb) * ptr_nh_prog%theta_v(1,jk,jb), &
            & ptr_nh_prog%theta_v(1,jk,jb),                              &
            & ptr_nh_prog%tracer(1,jk,jb,iqv)
        ENDDO
      ENDIF

      DO jk = 1 , nlev
         ! rhotheta has to have the same meaning as exner
         ptr_nh_prog%rhotheta_v(1:nlen,jk,jb) = &
              (ptr_nh_prog%exner(1:nlen,jk,jb)**cvd_o_rd)*p0ref/rd

         ptr_nh_prog%rho(1:nlen,jk,jb) = ptr_nh_prog%rhotheta_v(1:nlen,jk,jb) / &
                                         ptr_nh_prog%theta_v(1:nlen,jk,jb)     
      END DO !jk

    ENDDO !jb

!--------------------------------------------------------------------------------
    !Mean wind 
!--------------------------------------------------------------------------------
    i_startblk = ptr_patch%edges%start_blk(2,1)
    DO jb = i_startblk , nblks_e
     CALL get_indices_e(ptr_patch, jb, i_startblk, nblks_e, i_startidx, i_endidx, 2)
     DO jk = 1 , nlev 
      DO jc = i_startidx, i_endidx

        !Torus geometry is flat so zu is only function of height which is same for all cells
        !But it is kept varyign with jc,jb to introduce topography lateron
        jcn  =   ptr_patch%edges%cell_idx(jc,jb,1)
        jbn  =   ptr_patch%edges%cell_blk(jc,jb,1)
        zu   =   u_cbl(1) + u_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
        zv   =   v_cbl(1) + v_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)

        zvn1 =  zu * ptr_patch%edges%primal_normal_cell(jc,jb,1)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(jc,jb,1)%v2      
 
        jcn  =   ptr_patch%edges%cell_idx(jc,jb,2)
        jbn  =   ptr_patch%edges%cell_blk(jc,jb,2)
        zu   =   u_cbl(1) + u_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
        zv   =   v_cbl(1) + v_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
      
        zvn2 =  zu * ptr_patch%edges%primal_normal_cell(jc,jb,2)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(jc,jb,2)%v2      

        ptr_nh_prog%vn(jc,jk,jb) = ptr_int%c_lin_e(jc,1,jb)*zvn1 + &
                                   ptr_int%c_lin_e(jc,2,jb)*zvn2

        ptr_nh_ref%vn_ref(jc,jk,jb) = ptr_nh_prog%vn(jc,jk,jb)
      END DO
     END DO
    END DO     
    
    !W wind and reference
    CALL init_w(ptr_patch, ptr_int, ptr_nh_prog%vn, ptr_metrics%z_ifc, ptr_nh_prog%w)
    CALL sync_patch_array(SYNC_C, ptr_patch, ptr_nh_prog%w)
    ptr_nh_ref%w_ref = ptr_nh_prog%w


  END SUBROUTINE init_nh_state_cbl

  !>
  !! Initialization of prognostic state vector for the nh GATE test case 
  !!  with moisture
  !!
  !! @par Revision History
  !!
  !!
  SUBROUTINE init_nh_state_gate( ptr_patch, ptr_nh_prog,  ptr_nh_ref, ptr_nh_diag,  &
    &                           ptr_int, ptr_metrics)

    ! INPUT PARAMETERS:
    TYPE(t_patch),TARGET,  INTENT(IN)   :: &  !< patch on which computation is performed
      &  ptr_patch
    TYPE(t_int_state),     INTENT(IN)   :: &
      &  ptr_int
    TYPE(t_nh_prog),       INTENT(INOUT):: &  !< prognostic state vector
      &  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT):: &  !< diagnostic state vector
      &  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   :: &
      &  ptr_metrics                          !< NH metrics state
    TYPE(t_nh_ref),        INTENT(INOUT):: &  !< reference state vector
      &  ptr_nh_ref

    REAL(wp) :: rho_sfc, z_help(1:nproma), zvn1, zvn2, zu, zv, zt00, zh00
    INTEGER  :: je,jk,jb,i_startblk,i_startidx,i_endidx   !< loop indices
    INTEGER  :: nblks_c,npromz_c,nblks_e,npromz_e
    INTEGER  :: nlev, nlevp1                        !< number of full and half levels
    INTEGER  :: nlen, i_rcstartlev, jcn, jbn, ist, jg, ntropo, iunit, nk

    REAL(wp), ALLOCATABLE, DIMENSION(:) :: zz, z_tp, z_qv, z_u, z_v

    !Constant in horizontal
    REAL(wp), ALLOCATABLE ::  &
      tp_gate(:),             & !Potential temperature [K]
      u_gate (:),             & !Zonal wind speed      [m/s]
      v_gate (:),             & !Meridional            [m/s]
      qv_gate(:)                !specific humidity     [g/kg]
    CHARACTER(len=max_char_length), PARAMETER :: routine = 'mo_nh_torus_exp:init_nh_state_gate'
  !-------------------------------------------------------------------------
    
    !Open formatted file to read ls forcing data
    iunit = find_next_free_unit(10,20)
    OPEN (unit=iunit,file='ini_profils.dat',access='SEQUENTIAL', &
          form='FORMATTED', action='READ', status='OLD', IOSTAT=ist)

    IF(ist/=success)THEN
      CALL finish (TRIM(routine), 'open init_profils.dat failed')
    ENDIF  

    !Read the input file till end. The order of file assumed is:
    !Z(m) - w_ls - u_geo - v_geo - ddt_temp_hadv_ls - ddt_temp_rad_ls - ddt_qv_hadv_ls    
    
    !Read first line
    READ(iunit,*,IOSTAT=ist)
    IF(ist/=success)CALL finish(TRIM(routine), 'problem reading first line in ls_forcing.dat')

    READ(iunit,*,IOSTAT=ist)nk
    IF(ist/=success)CALL finish(TRIM(routine), 'problem reading second line in ls_forcing.dat')
   
    ALLOCATE( zz(nk), z_tp(nk), z_qv(nk), z_u(nk), z_v(nk))

    DO jk = nk , 1, -1
      READ(iunit,*,IOSTAT=ist)zz(jk),z_tp(jk),z_qv(jk),z_u(jk),z_v(jk)
      IF(ist/=success)THEN
          CALL finish (TRIM(routine), 'something wrong in ini_profils.dat')
      END IF
    END DO

    !Check if the file is written in descending order
    IF(zz(1) < zz(nk))CALL finish (TRIM(routine), 'Writing initial profils data in descending order!')

    CLOSE(iunit)

    nlev   = ptr_patch%nlev
    ALLOCATE( tp_gate(nlev), u_gate(nlev), v_gate(nlev), qv_gate(nlev) )

    !Now perform interpolation to grid levels assuming:
    !a) linear interpolation
    !b) Beyond the last Z level the values are linearly extrapolated 
    !c) Assuming model grid is flat-NOT on sphere

    CALL vert_intp_linear_1d(zz,z_tp,ptr_metrics%z_mc(2,:,2),tp_gate)
    CALL vert_intp_linear_1d(zz,z_u,ptr_metrics%z_mc(2,:,2),u_gate)
    CALL vert_intp_linear_1d(zz,z_v,ptr_metrics%z_mc(2,:,2),v_gate)
    CALL vert_intp_linear_1d(zz,z_qv,ptr_metrics%z_mc(2,:,2),qv_gate)

    DEALLOCATE( zz, z_tp, z_qv, z_u, z_v)

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_int_c
    npromz_c = ptr_patch%npromz_int_c
    nblks_e  = ptr_patch%nblks_int_e
    npromz_e = ptr_patch%npromz_int_e

    ! number of vertical levels
    nlevp1 = ptr_patch%nlevp1

    !patch id
    jg = ptr_patch%id

    !
    i_rcstartlev = 2
    i_startblk   = ptr_patch%cells%start_blk(i_rcstartlev,1)

    !Set some reference density    
    rho_sfc = zp_sfc / (rd * les_config(jg)%sst)

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = zp_sfc

    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp

    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = 1, nlev
        ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = qv_gate(jk)/1000._wp
      END DO

      ntropo = 0
      DO jk = nlev, 1, -1
         z_help(1:nlen) = tp_gate(jk)
         ! constant temperature above tropopause
         if ((ptr_metrics%z_mc(1,jk,jb) > zg_tropo) .and. (ntropo == 0)) then
            ntropo = 1
            zt00   = z_help(1)
            zh00   = ptr_metrics%z_mc(1,jk,jb)
         endif
         if (ptr_metrics%z_mc(1,jk,jb) > zg_tropo) then
            z_help(1:nlen) = zt00 + (ptr_metrics%z_mc(1:nlen,jk,jb)-zh00) * dtdz_st
         endif
         ptr_nh_prog%theta_v(1:nlen,jk,jb) = z_help(1:nlen) * ( 1._wp + &
           0.61_wp*ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) - ptr_nh_prog%tracer(1:nlen,jk,jb,iqc) ) 
      END DO


      !Get hydrostatic pressure and exner at lowest level
      ptr_nh_diag%pres(1:nlen,nlev,jb) = zp_sfc - rho_sfc * ptr_metrics%geopot(1:nlen,nlev,jb)
      ptr_nh_prog%exner(1:nlen,nlev,jb) = (ptr_nh_diag%pres(1:nlen,nlev,jb)/p0ref)**rd_o_cpd 

      !Get exner at other levels
      DO jk = nlev-1, 1, -1
         z_help(1:nlen) = 0.5_wp * ( ptr_nh_prog%theta_v(1:nlen,jk,jb) +  &
                                     ptr_nh_prog%theta_v(1:nlen,jk+1,jb) )
   
         ptr_nh_prog%exner(1:nlen,jk,jb) = ptr_nh_prog%exner(1:nlen,jk+1,jb) &
            &  -grav/cpd*ptr_metrics%ddqz_z_half(1:nlen,jk+1,jb)/z_help(1:nlen)
      END DO

      DO jk = 1 , nlev
         ! rhotheta has to have the same meaning as exner
         ptr_nh_prog%rhotheta_v(1:nlen,jk,jb) = &
              (ptr_nh_prog%exner(1:nlen,jk,jb)**cvd_o_rd)*p0ref/rd

         ptr_nh_prog%rho(1:nlen,jk,jb) = ptr_nh_prog%rhotheta_v(1:nlen,jk,jb) / &
                                         ptr_nh_prog%theta_v(1:nlen,jk,jb)     
      END DO !jk


      CALL diagnose_pres_temp (ptr_metrics, ptr_nh_prog,ptr_nh_prog, ptr_nh_diag,     &
                             ptr_patch, opt_calc_pres=.TRUE., opt_calc_temp=.TRUE.)
    ! Print case set up
    IF ( jb == 1 ) THEN
      DO jk = 1,nlev
       write(*,*) 'GATE case setup: level, p, T_diag, T, theta,v, qv: ', jk,&
          &ptr_nh_prog%exner(1,jk,jb)**(cpd/rd) * p0ref,                   &
          & ptr_nh_diag%temp(1,jk,jb),                                      &
          & ptr_nh_prog%exner(1,jk,jb) * ptr_nh_prog%theta_v(1,jk,jb),      &
          & ptr_nh_prog%theta_v(1,jk,jb),                                   &
          & ptr_nh_prog%tracer(1,jk,jb,iqv)
      ENDDO
      write(*,*) 'GATE case setup surface: Ps: ',     ptr_nh_diag%pres_sfc(1,jb)
    ENDIF


    ENDDO !jb

    !Mean wind 
    DO jb = 1 , nblks_e
     CALL get_indices_e( ptr_patch, jb, 1, nblks_e, i_startidx, i_endidx, grf_bdywidth_e+1)
     DO jk = 1 , nlev 
      DO je = i_startidx, i_endidx

        !Torus geometry is flat so zu is only function of height which is same for all cells
        !But it is kept varyign with jc,jb to introduce topography lateron
        jcn  =   ptr_patch%edges%cell_idx(je,jb,1)
        jbn  =   ptr_patch%edges%cell_blk(je,jb,1)
        zu   =   u_gate(jk)
        zv   =   v_gate(jk)

        zvn1 =  zu * ptr_patch%edges%primal_normal_cell(je,jb,1)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(je,jb,1)%v2      
 
        jcn  =   ptr_patch%edges%cell_idx(je,jb,2)
        jbn  =   ptr_patch%edges%cell_blk(je,jb,2)
        zu   =   u_gate(jk)
        zv   =   v_gate(jk)
      
        zvn2 =  zu * ptr_patch%edges%primal_normal_cell(je,jb,2)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(je,jb,2)%v2      

        ptr_nh_prog%vn(je,jk,jb) = ptr_int%c_lin_e(je,1,jb)*zvn1 + &
                                   ptr_int%c_lin_e(je,2,jb)*zvn2

        ptr_nh_ref%vn_ref(je,jk,jb) = ptr_nh_prog%vn(je,jk,jb)
      END DO
     END DO
    END DO     
    
    !W wind and reference
    ptr_nh_prog%w  = 0._wp; ptr_nh_ref%w_ref  = ptr_nh_prog%w

    DEALLOCATE( tp_gate, u_gate, v_gate, qv_gate )

  END SUBROUTINE init_nh_state_gate

  
  !>
  !! Initialization of prognostic state vector for the nh RICO test case 
  !!  with moisture
  !!
  !! @par Revision History
  !!
  !!
  SUBROUTINE init_nh_state_rico( ptr_patch, ptr_nh_prog,  ptr_nh_ref, ptr_nh_diag,  &
    &                           ptr_int, ptr_metrics)

    ! INPUT PARAMETERS:
    TYPE(t_patch),TARGET,  INTENT(IN)   :: &  !< patch on which computation is performed
      &  ptr_patch
    TYPE(t_int_state),     INTENT(IN)   :: &
      &  ptr_int
    TYPE(t_nh_prog),       INTENT(INOUT):: &  !< prognostic state vector
      &  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT):: &  !< diagnostic state vector
      &  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   :: &
      &  ptr_metrics                          !< NH metrics state
    TYPE(t_nh_ref),        INTENT(INOUT):: &  !< reference state vector
      &  ptr_nh_ref

    REAL(wp) :: rho_sfc, z_help(1:nproma), zvn1, zvn2, zu, zv
    INTEGER  :: je,jk,jb,i_startblk,i_startidx,i_endidx   !< loop indices
    INTEGER  :: nblks_c,npromz_c,nblks_e,npromz_e
    INTEGER  :: nlev, nlevp1                        !< number of full and half levels
    INTEGER  :: nlen, i_rcstartlev, jcn, jbn, jg

  !-------------------------------------------------------------------------

    ! values for the blocking
    nblks_c  = ptr_patch%nblks_int_c
    npromz_c = ptr_patch%npromz_int_c
    nblks_e  = ptr_patch%nblks_int_e
    npromz_e = ptr_patch%npromz_int_e

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    nlevp1 = ptr_patch%nlevp1

    !patch id
    jg = ptr_patch%id

    !
    i_rcstartlev = 2
    i_startblk   = ptr_patch%cells%start_blk(i_rcstartlev,1)

    !Set some reference density    
    rho_sfc = zpsfc / (rd * les_config(jg)%sst)

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = zpsfc

    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp

    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      !Tracers
        DO jk = 1, nlev
          ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = min((0.016_wp - 0.0022_wp * ptr_metrics%z_mc(1:nlen,jk,jb)/zh1), &
                                                (0.0138_wp - 0.0114_wp * (ptr_metrics%z_mc(1:nlen,jk,jb)-zh1)/(zh2 - zh1)))
          ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = max(ptr_nh_prog%tracer(1:nlen,jk,jb,iqv),             &
                                                (0.0024_wp - 0.0006_wp * (ptr_metrics%z_mc(1:nlen,jk,jb) - zh2)/(4000._wp - zh2)))
          ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) = max(ptr_nh_prog%tracer(1:nlen,jk,jb,iqv), 3e-6_wp)                            
        END DO

      DO jk = 1, nlev
       ! init potential temperature
       z_help(1:nlen) = ztsfc + max(0._wp, (ptr_metrics%z_mc(1:nlen,jk,jb)-zh1)*19.1_wp/(4000._wp-zh1))
       z_help(1:nlen) = max(z_help(1:nlen), ztsfc +(ptr_metrics%z_mc(1:nlen,jk,jb)-zh3)*502._wp/(zh4-zh3))
       z_help(1:nlen) = max(z_help(1:nlen), ztsfc +(ptr_metrics%z_mc(1:nlen,jk,jb)-zh4)*2600._wp/(zh5-zh4))

       ! virtual potential temperature
       ptr_nh_prog%theta_v(1:nlen,jk,jb) = z_help(1:nlen) * ( 1._wp + &
           0.61_wp*ptr_nh_prog%tracer(1:nlen,jk,jb,iqv) - ptr_nh_prog%tracer(1:nlen,jk,jb,iqc) ) 
      END DO

      !Get hydrostatic pressure and exner at lowest level
      ptr_nh_diag%pres(1:nlen,nlev,jb) = zpsfc - rho_sfc * ptr_metrics%geopot(1:nlen,nlev,jb)
      ptr_nh_prog%exner(1:nlen,nlev,jb) = (ptr_nh_diag%pres(1:nlen,nlev,jb)/p0ref)**rd_o_cpd 

      !Get exner at other levels
      DO jk = nlev-1, 1, -1
         z_help(1:nlen) = 0.5_wp * ( ptr_nh_prog%theta_v(1:nlen,jk,jb) +  &
                                     ptr_nh_prog%theta_v(1:nlen,jk+1,jb) )
   
         ptr_nh_prog%exner(1:nlen,jk,jb) = ptr_nh_prog%exner(1:nlen,jk+1,jb) &
            &  -grav/cpd*ptr_metrics%ddqz_z_half(1:nlen,jk+1,jb)/z_help(1:nlen)
      END DO

      DO jk = 1 , nlev
         ! rhotheta has to have the same meaning as exner
         ptr_nh_prog%rhotheta_v(1:nlen,jk,jb) = &
              (ptr_nh_prog%exner(1:nlen,jk,jb)**cvd_o_rd)*p0ref/rd

         ptr_nh_prog%rho(1:nlen,jk,jb) = ptr_nh_prog%rhotheta_v(1:nlen,jk,jb) / &
                                         ptr_nh_prog%theta_v(1:nlen,jk,jb)     
      END DO !jk

    ENDDO !jb


    !Mean wind 
    DO jb = 1 , nblks_e
     CALL get_indices_e( ptr_patch, jb, 1, nblks_e, i_startidx, i_endidx, grf_bdywidth_e+1)
     DO jk = 1 , nlev 
      DO je = i_startidx, i_endidx

        !Torus geometry is flat so zu is only function of height which is same for all cells
        !But it is kept varyign with jc,jb to introduce topography lateron
        jcn  =   ptr_patch%edges%cell_idx(je,jb,1)
        jbn  =   ptr_patch%edges%cell_blk(je,jb,1)
        zu   =   u_cbl(1) + u_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
        zu   =   min(zu, -2._wp)
        zv   =   v_cbl(1) + v_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)

        zvn1 =  zu * ptr_patch%edges%primal_normal_cell(je,jb,1)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(je,jb,1)%v2      
 
        jcn  =   ptr_patch%edges%cell_idx(je,jb,2)
        jbn  =   ptr_patch%edges%cell_blk(je,jb,2)
        zu   =   u_cbl(1) + u_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
        zu   =   min(zu, -2._wp)
        zv   =   v_cbl(1) + v_cbl(2) * ptr_metrics%z_mc(jcn,jk,jbn)
      
        zvn2 =  zu * ptr_patch%edges%primal_normal_cell(je,jb,2)%v1 + &
                zv * ptr_patch%edges%primal_normal_cell(je,jb,2)%v2      

        ptr_nh_prog%vn(je,jk,jb) = ptr_int%c_lin_e(je,1,jb)*zvn1 + &
                                   ptr_int%c_lin_e(je,2,jb)*zvn2

        ptr_nh_ref%vn_ref(je,jk,jb) = ptr_nh_prog%vn(je,jk,jb)
      END DO
     END DO
    END DO     
    
    !W wind and reference
    ptr_nh_prog%w  = 0._wp; ptr_nh_ref%w_ref  = ptr_nh_prog%w
    
  END SUBROUTINE init_nh_state_rico
  
!-------------------------------------------------------------------------


  SUBROUTINE cbl_stevens_fluxes( t_l, qv_l, p_l, rho_l, tsk, shfx, lhfx )

  !-------------------------------------------------------------------------
  ! Calculate sensible and latent heat fluxes from buoyancy flux for Stevens 
  ! (2007) case. (code from Wayne Angevine 2013) 
  !
  ! Variable explanations
  ! tsk  = skin temperature
  ! zp0  = surface pressure (Pa)
  ! qv1d = water vapor mixing ratio in lowest level
  ! th1d = potential temperature in lowest level
  ! kqfx = exchange coefficient for water vapor (used in qfx calculation below)
  !
  ! Constants for saturation calculation:
  ! svp1  = 0.6112
  ! svp2  = 17.67
  ! svp3  = 29.65
  !-------------------------------------------------------------------------

    ! INPUT PARAMETERS:
    REAL(wp), INTENT(IN) ::    t_l       , & ! temperature at lowest level [K]
     &                         qv_l      , & ! moisture at lowest level    [kg/kg]
     &                         p_l       , & ! pressure at lowest level    [Pa]
     &                         rho_l         ! rho at lowest level         [kg/m3]

    ! INPUT/OUTPUT PARAMETERS:
    REAL(wp), INTENT(INOUT) :: tsk           ! skin temperature [K]

    ! OUTPUT PARAMETERS:
    REAL(wp), INTENT(OUT) ::   shfx      , & ! sensible heat flux [W/m2]
      &                        lhfx          ! latent heat flux   [W/m2]

    ! LOCAL VARIABLES:
    REAL(wp) :: Beta0, Vs, mav, qsfc, qsfc_air, kqfx, khfx, th_l

  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
    Beta0 = 7e-4_wp       ! Fixed specified buoyancy flux (25 W/m^2)
   !Beta0 = 4.2e-4_wp     ! Fixed specified buoyancy flux (15 W/m^2)
   !Beta0 = 11.2e-4_wp    ! Fixed specified buoyancy flux (40 W/m^2)
    Vs    = 0.01_wp       ! Specified velocity scale
    mav   = 0.90_wp       ! Surface moisture availability

  ! Calculate saturation mixing ratio at SST (from previous timestep)
  ! The rest of the buoyancy flux goes to sensible heat
  ! Calculate surface saturated q and q in air at surface 
   !e1=svp1*exp(svp2*(tsk-tmelt)/(tsk-svp3))                       
   !qsfc=rd/rv*e1/((zp0/1000.)-e1)
    qsfc     = spec_humi(sat_pres_water(tsk),zp0)
    qsfc_air = qsfc * mav

    th_l =  t_l * (p0ref/p_l)**rd_o_cpd 

    ! Calculate hfx,qfx, and SST to keep buoyancy flux constant
    ! Could calculate moisture flux first, but should be OK either way
    kqfx = Vs * (qsfc_air - qv_l)   ! exchange coefficient for water vapor
    khfx = Beta0 * th_l/grav - 0.608_wp * th_l * kqfx
    tsk  = khfx / Vs + th_l

    ! Convert units
    shfx = rho_l * cpd * khfx
    lhfx = rho_l * alv * kqfx

  END SUBROUTINE cbl_stevens_fluxes


!-------------------------------------------------------------------------

  END MODULE mo_nh_torus_exp
