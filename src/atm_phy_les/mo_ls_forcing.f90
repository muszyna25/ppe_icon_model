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

  !Constant in time and horizontal space
  REAL(wp), SAVE, ALLOCATABLE ::  &
    w_ls(:),            & !subsidence          [m/s]
    u_geo  (:),         & !u-geostrophic wind  [m/s]
    v_geo  (:),         & !v-geostrophic wind  [m/s]
    ddt_temp_hadv_ls(:),& !LS horizontal advective tendency for temp [k/s]
    ddt_temp_rad_ls(:), & !LS (not large-scale actually!) radiative tendency for temp [k/s]
    ddt_qv_hadv_ls(:)     !LS horizontal advective tendency for qv [1/s]

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
    INTEGER   :: iunit, ist, nk, jk, nlev

    !Open formatted file to read ls forcing data
    iunit = find_next_free_unit(10,20)
    OPEN (unit=iunit,file='ls_forcing.dat',access='SEQUENTIAL', &
          form='FORMATTED', action='READ', status='OLD', IOSTAT=ist)

    IF(ist/=success)THEN
      CALL finish (routine, 'open ls_forcing.dat failed')
    ENDIF

    !Read the input file till end. The order of file assumed is:
    !Z(m) - w_ls - u_geo - v_geo - ddt_temp_hadv_ls - ddt_temp_rad_ls - ddt_qv_hadv_ls

    !Read first line
    READ(iunit,*,IOSTAT=ist)
    IF(ist/=success)CALL finish(routine, 'problem reading first line in ls_forcing.dat')

    READ(iunit,*,IOSTAT=ist)nk
    IF(ist/=success)CALL finish(routine, 'problem reading second line in ls_forcing.dat')

    ALLOCATE( zz(nk), zw(nk), zu(nk), zv(nk), z_dt_temp_adv(nk), &
              z_dt_temp_rad(nk), z_dt_qv_adv(nk) )

    DO jk = nk , 1, -1
      READ(iunit,*,IOSTAT=ist)zz(jk),zw(jk),zu(jk),zv(jk),z_dt_temp_adv(jk), &
                              z_dt_temp_rad(jk),z_dt_qv_adv(jk)
      IF(ist/=success)THEN
          CALL finish (routine, 'something wrong in forcing.dat')
      END IF
    END DO

    !Check if the file is written in descending order
    IF(zz(1) < zz(nk))CALL finish (routine, 'Write LS forcing data in descending order!')

    CLOSE(iunit)

    nlev = SIZE(p_metrics%z_mc,2)
    ALLOCATE( w_ls(nlev), u_geo(nlev), v_geo(nlev), ddt_temp_hadv_ls(nlev), &
              ddt_temp_rad_ls(nlev), ddt_qv_hadv_ls(nlev) )

    !Now perform interpolation to grid levels assuming:
    !a) linear interpolation
    !b) Beyond the last Z level the values are linearly extrapolated
    !c) Assuming model grid is flat-NOT on sphere

    CALL vert_intp_linear_1d(zz,zw,p_metrics%z_mc(1,:,1),w_ls)
    CALL vert_intp_linear_1d(zz,zu,p_metrics%z_mc(1,:,1),u_geo)
    CALL vert_intp_linear_1d(zz,zv,p_metrics%z_mc(1,:,1),v_geo)
    CALL vert_intp_linear_1d(zz,z_dt_temp_adv,p_metrics%z_mc(1,:,1),ddt_temp_hadv_ls)
    CALL vert_intp_linear_1d(zz,z_dt_temp_rad,p_metrics%z_mc(1,:,1),ddt_temp_rad_ls)
    CALL vert_intp_linear_1d(zz,z_dt_qv_adv,p_metrics%z_mc(1,:,1),ddt_qv_hadv_ls)

    DEALLOCATE( zz, zw, zu, zv, z_dt_temp_adv, z_dt_temp_rad, z_dt_qv_adv )


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
  SUBROUTINE apply_ls_forcing(p_patch, p_metrics, p_prog, p_diag,      & !in
                              qv, rl_start, rl_end,                    & !in
                              ddt_u_ls, ddt_v_ls, ddt_temp_ls,         & !out
                              ddt_qv_ls                      )           !out

    TYPE(t_patch),INTENT(in),        TARGET :: p_patch
    TYPE(t_nh_metrics),INTENT(in),   TARGET :: p_metrics
    TYPE(t_nh_prog),INTENT(in),      TARGET :: p_prog
    TYPE(t_nh_diag),INTENT(in),      TARGET :: p_diag
    INTEGER,        INTENT(in)              :: rl_start
    INTEGER,        INTENT(in)              :: rl_end
    REAL(wp),       INTENT(in)              :: qv(:,:,:)
    REAL(wp),       INTENT(out)             :: ddt_u_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_v_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_temp_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_qv_ls(:)

    CHARACTER(len=max_char_length), PARAMETER :: routine = 'mo_ls_forcing:apply_ls_forcing'
    REAL(wp) :: varin(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp), DIMENSION(p_patch%nlev+1) :: u_gb_hl, v_gb_hl, th_gb_hl, qv_gb_hl
    REAL(wp), DIMENSION(p_patch%nlev)   :: inv_dz, exner_gb, inv_exner_gb, var_gb
    REAL(wp) :: inv_no_gb_cells
    INTEGER  :: i_nchdom, i_startblk, i_endblk, jk, nlev, nlevp1

    !0) Initialize all passed ddt's to 0
!$OMP PARALLEL
    CALL init(ddt_u_ls); CALL init(ddt_v_ls);  CALL init(ddt_temp_ls)
    call init(ddt_qv_ls)
!$OMP END PARALLEL

    i_nchdom  = MAX(1,p_patch%n_childdom)
    nlev      = p_patch%nlev
    nlevp1    = p_patch%nlev+1
    inv_no_gb_cells = 1._wp / REAL(p_patch%n_patch_cells_g,wp)

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    !Some precalculations
    IF(is_subsidence_heat .OR. is_advection .OR. is_rad_forcing) &
      CALL global_hor_mean(p_patch, p_prog%exner(:,:,:), exner_gb, inv_no_gb_cells, i_nchdom)

    inv_exner_gb = 1./exner_gb

    IF(is_subsidence_moment.OR.is_subsidence_heat) &
      inv_dz(:) = 1._wp / p_metrics%ddqz_z_full(1,:,1)

    !1a) Horizontal mean of variables and their vertical advective tendency - momentum

    IF(is_subsidence_moment)THEN
      CALL global_hor_mean(p_patch, p_diag%u, var_gb, inv_no_gb_cells, i_nchdom)
      CALL vert_intp_linear_1d(p_metrics%z_mc(1,:,1),var_gb,p_metrics%z_ifc(1,:,1),u_gb_hl)

      CALL global_hor_mean(p_patch, p_diag%v, var_gb, inv_no_gb_cells, i_nchdom)
      CALL vert_intp_linear_1d(p_metrics%z_mc(1,:,1),var_gb,p_metrics%z_ifc(1,:,1),v_gb_hl)

      ddt_u_ls =  ddt_u_ls - w_ls*vertical_derivative(u_gb_hl,inv_dz)
      ddt_v_ls =  ddt_v_ls - w_ls*vertical_derivative(v_gb_hl,inv_dz)
    END IF

    !1b) Horizontal mean of variables and their vertical advective tendency

    IF(is_subsidence_heat)THEN
      !use theta instead of temperature for subsidence
!$OMP PARALLEL WORKSHARE
        varin(:,:,i_startblk:i_endblk) = p_diag%temp(:,:,i_startblk:i_endblk) / &
                                         p_prog%exner(:,:,i_startblk:i_endblk)
!$OMP END PARALLEL WORKSHARE
      CALL global_hor_mean(p_patch, varin, var_gb, inv_no_gb_cells, i_nchdom)
      CALL vert_intp_linear_1d(p_metrics%z_mc(1,:,1),var_gb,p_metrics%z_ifc(1,:,1),th_gb_hl)

      CALL global_hor_mean(p_patch, qv, var_gb, inv_no_gb_cells, i_nchdom)
      CALL vert_intp_linear_1d(p_metrics%z_mc(1,:,1),var_gb,p_metrics%z_ifc(1,:,1),qv_gb_hl)

      ddt_temp_ls = ddt_temp_ls - w_ls*vertical_derivative(th_gb_hl,inv_dz)
      ddt_qv_ls   = ddt_qv_ls   - w_ls*vertical_derivative(qv_gb_hl,inv_dz)
    END IF

    !2)Horizontal advective forcing
    IF(is_advection)THEN

      ddt_qv_ls   = ddt_qv_ls + ddt_qv_hadv_ls

      IF(is_theta)THEN
        ddt_temp_ls = ddt_temp_ls + ddt_temp_hadv_ls
      ELSE
        ddt_temp_ls = ddt_temp_ls + inv_exner_gb*ddt_temp_hadv_ls
      END IF

    END IF

    !3)Coriolis and geostrophic wind
    IF(is_geowind)THEN
      !Remember model grid is flat
      ddt_u_ls = ddt_u_ls - p_patch%cells%f_c(1,1) * v_geo
      ddt_v_ls = ddt_v_ls + p_patch%cells%f_c(1,1) * u_geo
    END IF

    !4)Radiative forcing
    IF(is_rad_forcing)THEN
      IF(is_theta)THEN
        ddt_temp_ls = ddt_temp_ls + ddt_temp_rad_ls
      ELSE
        ddt_temp_ls = ddt_temp_ls + ddt_temp_rad_ls * inv_exner_gb
      END IF
    END IF

    !Convert theta tendency to temp at once
    ddt_temp_ls = exner_gb * ddt_temp_ls

  END SUBROUTINE apply_ls_forcing

!-------------------------------------------------------------------------------

END MODULE mo_ls_forcing



