!>
!! mo_ls_forcing
!!
!! This module initializes and applies the large-scale forcing for idealized simulations
!! This module assumes that the model grid is FLAT
!! 2013-JUNE-04: AT THIS STAGE LS FORCING WILL WORK IN RESTART MODE ONLY IF ITS CALLED EVERY
!! DYNAMIC TIMESTEP. TO MAKE IT WORK "SMOOTHLY" ADD_VAR HAS TO WORK ON 1D VARS
!!
!! @author Anurag Dipankar, MPI-M
!!
!! @par Revision History
!! Initial release by Anurag Dipankar, MPI-M (2013-May-30)
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

MODULE mo_ls_forcing

  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: filename_max, find_next_free_unit
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_impl_constants,      ONLY: success, max_char_length
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c, get_indices_v
  USE mo_les_utilities
  USE mo_ls_forcing_nml
  USE mo_physical_constants,  ONLY: rd, cpd, alv, cvd
  USE mo_fortran_tools,       ONLY: init

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_ls_forcing, apply_ls_forcing

  !Anurag Dipankar, MPIM (2013-May): variables for large-scale (LS)
  !forcing to be used in periodic domain or other relveant cases.

  !Dependent on z and time only
  REAL(wp), SAVE, ALLOCATABLE ::  &
    w_ls(:,:),            & !subsidence          [m/s]
    u_geo  (:,:),         & !u-geostrophic wind  [m/s]
    v_geo  (:,:),         & !v-geostrophic wind  [m/s]
    ddt_temp_hadv_ls(:,:),& !LS horizontal advective tendency for temp [k/s]
    ddt_temp_rad(:,:),    & !radiative tendency for temp [k/s]
    ddt_qv_hadv_ls(:,:),  & !LS horizontal advective tendency for qv [1/s]
    u_nudg(:,:),          & !LS nudging u    [m/s]
    v_nudg (:,:),         & !LS nudging v    [m/s]
    temp_nudg(:,:),       & !LS nudging temp [k]
    qv_nudg(:,:)            !LS nudging qv   [kg/kg]

    ! To make sure that LS forcing is not applied before it is read in
    REAL(wp), SAVE :: dt_forcing = 0._wp
    REAL(wp), SAVE :: dt_relax   = 0._wp
    REAL(wp), SAVE :: dt_nudging = 0._wp

  CONTAINS


  !>
  !! init_ls_forcing
  !!------------------------------------------------------------------------
  !! Initialize large-scale forcings from input file
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-May-30)
  SUBROUTINE init_ls_forcing(p_metrics)

    TYPE(t_nh_metrics),INTENT(in),   TARGET :: p_metrics

    CHARACTER(len=max_char_length), PARAMETER :: routine = 'mo_ls_forcing:init_ls_forcing'

    REAL(wp), ALLOCATABLE, DIMENSION(:) :: zz, zw, zu, zv, z_dt_temp_adv, &
                                           z_dt_temp_rad, z_dt_qv_adv
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: z_temp, z_qv 
    REAL(wp)  :: end_time
    INTEGER   :: iunit, ist, nk, nt, jk, nlev, i, nskip, n

    nlev   = SIZE(p_metrics%z_mc,2)

    !------------------------------------------------------------------------------
    ! READ LS FORCING FILE
    !------------------------------------------------------------------------------

    IF(is_ls_forcing)THEN

      !Open formatted file to read ls forcing data
      iunit = find_next_free_unit(10,20)
      OPEN (unit=iunit,file='ls_forcing.dat',access='SEQUENTIAL', &
            form='FORMATTED', action='READ', status='OLD', IOSTAT=ist)

      IF(ist/=success)THEN
        CALL finish (TRIM(routine), 'open ls_forcing.dat failed')
      ENDIF

      !Read the input file till end. The order of file assumed is:
      !Z(m) - u_geo(m/s) - v_geo(m/s) - w_ls(m/s) - dt_temp_rad(K/s) - ddt_qv_hadv_ls(1/s) - ddt_temp_hadv_ls(K/s)

      !Skip the first line and read next 2 lines with information about vertical and time levels
      READ(iunit,*,IOSTAT=ist)                       !skip
      READ(iunit,*,IOSTAT=ist)nk,dt_forcing,end_time !vertical levels, forcing interval, end of forcing data
      IF(ist/=success) &
        CALL finish (TRIM(routine), &
          'Must provide vertical level, forcing interval, end of forcing time info in forcing file')

      nt = INT(end_time/dt_forcing)+1

      IF(nt>1)THEN
        READ(iunit,*,IOSTAT=ist)nskip !lines to skip between successive time levels
        IF(ist/=success) &
    	    CALL finish (TRIM(routine), &
        	  'Time levels > 1 so must provide number lines to skip between successive time levels in forcing file')
      END IF

      IF(nt<=1)dt_forcing=-999._wp

      ALLOCATE( zz(nk), zw(nk), zu(nk), zv(nk), z_dt_temp_adv(nk), &
        	z_dt_temp_rad(nk), z_dt_qv_adv(nk) )

      ALLOCATE( w_ls(nlev,nt), u_geo(nlev,nt), v_geo(nlev,nt), ddt_temp_hadv_ls(nlev,nt), &
    		ddt_temp_rad(nlev,nt), ddt_qv_hadv_ls(nlev,nt) )

      DO n = 1,nt

        DO jk = nk , 1, -1
          READ(iunit,*,IOSTAT=ist)zz(jk),zu(jk),zv(jk),zw(jk),z_dt_temp_rad(jk), &
        	  		  z_dt_qv_adv(jk),z_dt_temp_adv(jk)
          IF(ist/=success)THEN
            CALL finish (TRIM(routine), 'something wrong in forcing.dat')
          END IF
        END DO

        !Skip lines
        IF(nt>1)THEN
          DO i = 1 , nskip
            READ(iunit,*,IOSTAT=ist)
          END DO
        END IF

        !Check if the file is written in descending order
        IF(zz(1) < zz(nk))CALL finish (TRIM(routine), 'Write LS forcing data in descending order!')

        !Now perform interpolation to grid levels assuming:
        !a) linear interpolation
        !b) Beyond the last Z level the values are linearly extrapolated
        !c) Assuming model grid is flat-NOT on sphere

        CALL vert_intp_linear_1d(zz,zw,p_metrics%z_mc(1,:,1),w_ls(:,n))
        CALL vert_intp_linear_1d(zz,zu,p_metrics%z_mc(1,:,1),u_geo(:,n))
        CALL vert_intp_linear_1d(zz,zv,p_metrics%z_mc(1,:,1),v_geo(:,n))
        CALL vert_intp_linear_1d(zz,z_dt_temp_adv,p_metrics%z_mc(1,:,1),ddt_temp_hadv_ls(:,n))
        CALL vert_intp_linear_1d(zz,z_dt_temp_rad,p_metrics%z_mc(1,:,1),ddt_temp_rad(:,n))
        CALL vert_intp_linear_1d(zz,z_dt_qv_adv,p_metrics%z_mc(1,:,1),ddt_qv_hadv_ls(:,n))

      END DO !n

      CLOSE(iunit)

      DEALLOCATE( zz, zw, zu, zv, z_dt_temp_adv, z_dt_temp_rad, z_dt_qv_adv )

    END IF !is_ls_forcing

    !------------------------------------------------------------------------------
    ! READ NUDGING FILE
    !------------------------------------------------------------------------------


    IF(is_nudging) THEN

      !Open formatted file to read nudging data
      iunit = find_next_free_unit(10,20)
      OPEN (unit=iunit,file='nudging.dat',access='SEQUENTIAL', &
            form='FORMATTED', action='READ', status='OLD', IOSTAT=ist)

      IF(ist/=success)THEN
        CALL finish (TRIM(routine), 'open nudging.dat failed')
      ENDIF

       !Read the input file till end. The order of file assumed is:
       !Z(m) - tau(s) - u(m/s) - v(m/s) - temp(K) - qv(kg/kg) 

       !Skip the first line and read next 2 lines with information about vertical and time levels
       READ(iunit,*,IOSTAT=ist)                         !skip
       READ(iunit,*,IOSTAT=ist)nk,dt_nudging,end_time   !vertical levels,time levels and interval
       IF(ist/=success) &
          CALL finish (TRIM(routine), &
          'Must provide vertical level, time level and time interval info in nudging file')

       IF(nt>1)THEN
         READ(iunit,*,IOSTAT=ist)nskip              !lines to skip between successive time levels
         IF(ist/=success) &
           CALL finish (TRIM(routine), &
           'Time levels > 1 so must provide number lines to skip between successive time levels in nudging file')
       END IF
       IF(nt<=1)dt_nudging=-999._wp

       ALLOCATE( zz(nk), zu(nk), zv(nk), z_temp(nk), z_qv(nk) )

       ALLOCATE( u_nudg(nlev,nt), v_nudg(nlev,nt), temp_nudg(nlev,nt), qv_nudg(nlev,nt))

       DO n = 1 , nt

         DO jk = nk , 1, -1
           READ(iunit,*,IOSTAT=ist)zz(jk),dt_relax,zu(jk),zv(jk),z_temp(jk),z_qv(jk)
           IF(ist/=success)THEN
               CALL finish (TRIM(routine), 'something wrong in nudging.dat')
           END IF
         END DO

         !Skip lines
         IF(nt>1)THEN
          DO i = 1 , nskip
            READ(iunit,*,IOSTAT=ist)
          END DO
         END IF

         !Check if the file is written in descending order
         IF(zz(1) < zz(nk))CALL finish (TRIM(routine), 'Write nuding data in descending order!')

         !Now perform interpolation to grid levels assuming:
         !a) linear interpolation
         !b) Beyond the last Z level the values are linearly extrapolated
         !c) Assuming model grid is flat-NOT on sphere

         CALL vert_intp_linear_1d(zz,zu,p_metrics%z_mc(1,:,1),u_nudg(:,n))
         CALL vert_intp_linear_1d(zz,zv,p_metrics%z_mc(1,:,1),v_nudg(:,n))
         CALL vert_intp_linear_1d(zz,z_temp,p_metrics%z_mc(1,:,1),temp_nudg(:,n))
         CALL vert_intp_linear_1d(zz,z_qv,p_metrics%z_mc(1,:,1),qv_nudg(:,n))

       END DO

       CLOSE(iunit)

       DEALLOCATE( zz, zu, zv, z_temp, z_qv )

     END IF !is_nudging

     WRITE(message_text,*)dt_forcing
     CALL message('Time varying LS forcing read in:',message_text)

  END SUBROUTINE init_ls_forcing

  !>
  !! apply_ls_forcing
  !!------------------------------------------------------------------------
  !! Apply large-scale forcing: called in the end from mo_nh_interface_nwp
  !! It uses the most updated u,v to calculate the large-scale subsidence induced
  !! advective tendencies. All other tendencies don't need any computation.
  !! All tendencies are then accumulated in the slow physics tendency terms
  !!
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-May-30)
  SUBROUTINE apply_ls_forcing(p_patch, p_metrics, curr_sim_time,    & !in
                              p_prog, p_diag, qv, rl_start, rl_end, & !in
                              ddt_u_ls, ddt_v_ls,                   & !out
                              ddt_temp_ls, ddt_qv_ls,               & !out
                              ddt_temp_subs_ls, ddt_qv_subs_ls,     & !output
                              ddt_temp_adv_ls, ddt_qv_adv_ls,       & !output
                              ddt_temp_nud_ls, ddt_qv_nud_ls, wsub)   !output

    TYPE(t_patch),INTENT(in),        TARGET :: p_patch
    TYPE(t_nh_metrics),INTENT(in),   TARGET :: p_metrics
    TYPE(t_nh_prog),INTENT(in),      TARGET :: p_prog
    TYPE(t_nh_diag),INTENT(in),      TARGET :: p_diag
    INTEGER,        INTENT(in)              :: rl_start
    INTEGER,        INTENT(in)              :: rl_end
    REAL(wp),       INTENT(in)              :: curr_sim_time
    REAL(wp),       INTENT(in)              :: qv(:,:,:)
    REAL(wp),       INTENT(out)             :: ddt_u_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_v_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_temp_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_qv_ls(:)

    !Individual LS tendencies for profile output
    REAL(wp),       INTENT(out)             :: ddt_temp_subs_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_qv_subs_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_temp_adv_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_qv_adv_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_temp_nud_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_qv_nud_ls(:)
    REAL(wp),       INTENT(out)             :: wsub(:)

    CHARACTER(len=max_char_length), PARAMETER :: routine = 'mo_ls_forcing:apply_ls_forcing'
    REAL(wp) :: varin(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: inv_no_gb_cells, inv_dt_relax
    REAL(wp), DIMENSION(p_patch%nlev)  :: u_gb, v_gb, temp_gb, qv_gb, w_gb
    REAL(wp), DIMENSION(p_patch%nlev+1):: u_gb_hl, v_gb_hl, temp_gb_hl, qv_gb_hl
    REAL(wp), DIMENSION(p_patch%nlev)  :: inv_dz, exner_gb
    REAL(wp), DIMENSION(p_patch%nlev)  :: z_ddt_t_rad, z_ugeo, z_vgeo

    REAL(wp), DIMENSION(p_patch%nlev)  :: nudg_tend_u,nudg_tend_v,nudg_tend_qv,nudg_tend_temp !,z_wsub

    REAL(wp) :: int_weight
    INTEGER  :: i_nchdom, i_startblk, i_endblk, jk, nlev, jb, jc
    INTEGER  :: n_curr, n_next, i_startidx, i_endidx

    i_nchdom  = MAX(1,p_patch%n_childdom)
    nlev      = p_patch%nlev
    inv_no_gb_cells = 1._wp / REAL(p_patch%n_patch_cells_g,wp)

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    !0) Initialize all passed ddt's to 0
!$OMP PARALLEL
    CALL init(ddt_u_ls); CALL init(ddt_v_ls);  CALL init(ddt_temp_ls)
    call init(ddt_qv_ls)
!$OMP END PARALLEL

    !use theta instead of temperature for subsidence (Anurag, Christopher)
!$OMP PARALLEL
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
       DO jk = 1 , nlev
         DO jc = i_startidx, i_endidx
           varin(jc,jk,jb)  = p_diag%temp(jc,jk,jb)/p_prog%exner(jc,jk,jb)
         END DO
       END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    CALL global_hor_mean(p_patch, p_diag%u, u_gb, inv_no_gb_cells, i_nchdom)
    CALL global_hor_mean(p_patch, p_diag%v, v_gb, inv_no_gb_cells, i_nchdom)
    CALL global_hor_mean(p_patch, varin   , temp_gb, inv_no_gb_cells, i_nchdom)
    CALL global_hor_mean(p_patch, qv, qv_gb, inv_no_gb_cells, i_nchdom)
    CALL global_hor_mean(p_patch, p_prog%exner(:,:,:), exner_gb, inv_no_gb_cells, i_nchdom)

    IF(is_subsidence_moment.OR.is_subsidence_heat) &
      inv_dz(:) = 1._wp / p_metrics%ddqz_z_full(1,:,1)

    !Find where in the array of LS forcing and nudging current time stands
    !and do linear interpolation in time
    IF(dt_forcing > 0._wp)THEN
      n_curr = FLOOR(curr_sim_time/dt_forcing)+1
      n_next = n_curr+1
      int_weight = curr_sim_time/dt_forcing-n_curr+1

      ! Christopher Moseley: temporary catch
      IF (int_weight.LT.0 .OR.int_weight.GT.1) THEN
        WRITE(message_text,*)n_curr,n_next,int_weight
        CALL message('INTERPOLATION ERROR in LS forcing:',message_text)
      END IF

      wsub            = w_ls(:,n_curr)*(1.-int_weight)+w_ls(:,n_next)*int_weight
      ddt_qv_adv_ls   = ddt_qv_hadv_ls(:,n_curr)*(1.-int_weight)+ddt_qv_hadv_ls(:,n_next)*int_weight
      ddt_temp_adv_ls = ddt_temp_hadv_ls(:,n_curr)*(1.-int_weight)+ddt_temp_hadv_ls(:,n_next)*int_weight
      z_ddt_t_rad     = ddt_temp_rad(:,n_curr)*(1.-int_weight)+ddt_temp_rad(:,n_next)*int_weight
      z_ugeo          = u_geo(:,n_curr)*(1.-int_weight)+u_geo(:,n_next)*int_weight
      z_vgeo          = v_geo(:,n_curr)*(1.-int_weight)+v_geo(:,n_next)*int_weight
    ELSE
      wsub            = w_ls(:,1)
      ddt_qv_adv_ls   = ddt_qv_hadv_ls(:,1)
      ddt_temp_adv_ls = ddt_temp_hadv_ls(:,1)
      z_ddt_t_rad     = ddt_temp_rad(:,1)
      z_ugeo          = u_geo(:,1)
      z_vgeo          = v_geo(:,1)
    END IF

    !1a) Horizontal interpolation and their advective tendency - momentum

    IF(is_subsidence_moment)THEN
      CALL vert_intp_linear_1d(p_metrics%z_mc(1,:,1),u_gb,p_metrics%z_ifc(1,:,1),u_gb_hl)
      CALL vert_intp_linear_1d(p_metrics%z_mc(1,:,1),v_gb,p_metrics%z_ifc(1,:,1),v_gb_hl)

      ddt_u_ls =  ddt_u_ls - wsub*vertical_derivative(u_gb_hl,inv_dz)
      ddt_v_ls =  ddt_v_ls - wsub*vertical_derivative(v_gb_hl,inv_dz)
    END IF

    !1b) Horizontal interpolation and their advective tendency - temperature and water vapor

    IF(is_subsidence_heat)THEN
      CALL vert_intp_linear_1d(p_metrics%z_mc(1,:,1),temp_gb,p_metrics%z_ifc(1,:,1),temp_gb_hl)
      CALL vert_intp_linear_1d(p_metrics%z_mc(1,:,1),qv_gb,p_metrics%z_ifc(1,:,1),qv_gb_hl)

      ddt_temp_subs_ls = -wsub*vertical_derivative(temp_gb_hl,inv_dz)  ! vertical advective tendency for theta
      ddt_qv_subs_ls   = -wsub*vertical_derivative(qv_gb_hl  ,inv_dz)  ! vertical advective tendenciy for qv

      ddt_temp_ls = ddt_temp_ls + ddt_temp_subs_ls
      ddt_qv_ls   = ddt_qv_ls   + ddt_qv_subs_ls
    END IF

    !2)Horizontal advective forcing: at present only for tracers and temperature
    IF(is_advection)THEN
      ddt_qv_ls   = ddt_qv_ls + ddt_qv_adv_ls
      ddt_temp_ls = ddt_temp_ls + ddt_temp_adv_ls
    END IF

    !3)Coriolis and geostrophic wind
    IF(is_geowind)THEN
      !Remember model grid is flat
      ddt_u_ls = ddt_u_ls - p_patch%cells%f_c(1,1) * z_vgeo
      ddt_v_ls = ddt_v_ls + p_patch%cells%f_c(1,1) * z_ugeo
    END IF

    !4)Radiative forcing
    IF(is_rad_forcing)THEN
      ddt_temp_ls = ddt_temp_ls + z_ddt_t_rad
    END IF

    !5)Nudging
    IF(is_nudging) THEN

     !WRITE(message_text,*)dt_nudging
      !CALL message('dt_nudging =',message_text)

      IF (dt_nudging > 0._wp) THEN
        n_curr = FLOOR(curr_sim_time/dt_nudging)+1
        n_next = n_curr+1
        int_weight = curr_sim_time/dt_nudging-n_curr+1

        inv_dt_relax    = 1._wp / dt_relax

        ! Interpolation error catch
        IF (int_weight.LT.0 .OR. int_weight.GT.1) THEN
          WRITE(message_text,*)n_curr,n_next,int_weight
          CALL message('INTERPOLATION ERROR in nudging:',message_text)
        END IF

        ddt_u_ls    =  ddt_u_ls  - ( u_gb  - &
                       u_nudg(:,n_curr)*(1.-int_weight)-u_nudg(:,n_next)*int_weight)*inv_dt_relax
        ddt_v_ls    =  ddt_v_ls  - ( v_gb  - &
                       v_nudg(:,n_curr)*(1.-int_weight)-v_nudg(:,n_next)*int_weight)*inv_dt_relax
        ddt_qv_nud_ls   =  - (qv_gb - qv_nudg(:,n_curr)*(1.-int_weight)-qv_nudg(:,n_next)*int_weight)*inv_dt_relax
        ddt_temp_nud_ls =  - (temp_gb - temp_nudg(:,n_curr)*(1.-int_weight)-temp_nudg(:,n_next)*int_weight)*inv_dt_relax
      ELSE
        ddt_u_ls    =  ddt_u_ls  - ( u_gb  - u_nudg(:,1) ) * inv_dt_relax
        ddt_v_ls    =  ddt_v_ls  - ( v_gb  - v_nudg(:,1) ) * inv_dt_relax
        ddt_qv_nud_ls   =  - (qv_gb - qv_nudg(:,1)) * inv_dt_relax
        ddt_temp_nud_ls =  - (temp_gb - temp_nudg(:,1)) * inv_dt_relax
      END IF

      ddt_qv_ls   =  ddt_qv_ls   + ddt_qv_nud_ls
      ddt_temp_ls =  ddt_temp_ls + ddt_temp_nud_ls

    END IF

    !Convert theta tendency to temp at once
    ddt_temp_ls = exner_gb * ddt_temp_ls

  END SUBROUTINE apply_ls_forcing

!-------------------------------------------------------------------------------

END MODULE mo_ls_forcing



