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
!! @par Copyright
!! 2002-2013 by DWD and MPI-M
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

MODULE mo_ls_forcing

  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: filename_max, find_next_free_unit
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_impl_constants,      ONLY: success, max_char_length, min_rlcell_int
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_model_domain,        ONLY: t_patch
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_parallel_config,     ONLY: nproma
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c, get_indices_v
  USE mo_vert_utilities,      ONLY: vert_intp_full2half_cell_3d, vert_intp_linear_1d
  USE mo_ls_forcing_nml
  USE mo_physical_constants,  ONLY: rd
  USE mo_sync,                ONLY: global_sum_array, omp_global_sum_array
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'

  PUBLIC :: init_ls_forcing, apply_ls_forcing

  !Anurag Dipankar, MPIM (2013-May): variables for large-scale (LS) 
  !forcing to be used in periodic domain or other relveant cases.
    
  !Constant in time and horizontal space
  REAL(wp), ALLOCATABLE ::  &
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
      CALL finish (TRIM(routine), 'open ls_forcing.dat failed')
    ENDIF  

    !Read the input file till end. The order of file assumed is:
    !Z(m) - w_ls - u_geo - v_geo - ddt_temp_hadv_ls - ddt_temp_rad_ls - ddt_qv_hadv_ls    
    
    !Read first line
    READ(iunit,*,IOSTAT=ist)
    IF(ist/=success)CALL finish(TRIM(routine), 'problem reading first line in ls_forcing.dat')

    READ(iunit,*,IOSTAT=ist)nk
    IF(ist/=success)CALL finish(TRIM(routine), 'problem reading second line in ls_forcing.dat')
   
    ALLOCATE( zz(nk), zw(nk), zu(nk), zv(nk), z_dt_temp_adv(nk), &
              z_dt_temp_rad(nk), z_dt_qv_adv(nk) )

    DO jk = nk , 1, -1
      READ(iunit,*,IOSTAT=ist)zz(jk),zw(jk),zu(jk),zv(jk),z_dt_temp_adv(jk), &
                              z_dt_temp_rad(jk),z_dt_qv_adv(jk)
      IF(ist/=success)THEN
          CALL finish (TRIM(routine), 'something wrong in forcing.dat')
      END IF
    END DO
    
    !Check if the file is written in descending order
    IF(zz(1) < zz(nk))CALL finish (TRIM(routine), 'Write LS forcing data in descending order!')

    CLOSE(iunit)

    nlev = SIZE(p_metrics%z_mc,2)
    ALLOCATE( w_ls(nlev), u_geo(nlev), v_geo(nlev), ddt_temp_hadv_ls(nlev), &
              ddt_temp_rad_ls(nlev), ddt_qv_hadv_ls(nlev) )

    !Now perform interpolation to grid levels assuming:
    !a) linear interpolation
    !b) Beyond the last Z level the values are linearly extrapolated 
    !c) Assuming model grid is flat-NOT on sphere

    CALL vert_intp_linear_1d(zz,zw,p_metrics%z_mc(2,:,2),w_ls)
    CALL vert_intp_linear_1d(zz,zu,p_metrics%z_mc(2,:,2),u_geo)
    CALL vert_intp_linear_1d(zz,zv,p_metrics%z_mc(2,:,2),v_geo)
    CALL vert_intp_linear_1d(zz,z_dt_temp_adv,p_metrics%z_mc(2,:,2),ddt_temp_hadv_ls)
    CALL vert_intp_linear_1d(zz,z_dt_temp_rad,p_metrics%z_mc(2,:,2),ddt_temp_rad_ls)
    CALL vert_intp_linear_1d(zz,z_dt_qv_adv,p_metrics%z_mc(2,:,2),ddt_qv_hadv_ls)

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
                              qv, ql, qv_sfc, t_sfc, rl_start, rl_end, & !in
                              ddt_u_ls, ddt_v_ls, ddt_temp_ls,         & !out
                              ddt_qv_ls, ddt_ql_ls, ddt_qi_ls)           !out

    TYPE(t_patch),INTENT(in),        TARGET :: p_patch
    TYPE(t_nh_metrics),INTENT(in),   TARGET :: p_metrics 
    TYPE(t_nh_prog),INTENT(in),      TARGET :: p_prog
    TYPE(t_nh_diag),INTENT(in),      TARGET :: p_diag
    INTEGER,        INTENT(in)              :: rl_start
    INTEGER,        INTENT(in)              :: rl_end
    REAL(wp),       INTENT(in)              :: qv(:,:,:)
    REAL(wp),       INTENT(in)              :: ql(:,:,:)
    REAL(wp),       INTENT(in)              :: t_sfc(:,:)
    REAL(wp),       INTENT(in)              :: qv_sfc(:,:)
    REAL(wp),       INTENT(out)             :: ddt_u_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_v_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_temp_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_qv_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_ql_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_qi_ls(:)
    
    CHARACTER(len=max_char_length), PARAMETER :: routine = 'mo_ls_forcing:apply_ls_forcing'
    REAL(wp) :: varin(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: varout(nproma,p_patch%nlev+1,p_patch%nblks_c)
    REAL(wp) :: rhos(nproma,p_patch%nblks_c)
    REAL(wp) :: inv_no_gb_cells
    REAL(wp), DIMENSION(p_patch%nlev+1) :: u_gb, v_gb, temp_gb, qv_gb, ql_gb
    REAL(wp), DIMENSION(p_patch%nlev)   :: inv_rho_gb, rho_gb, inv_dz, exner_gb

    INTEGER  :: i_nchdom, i_startblk, i_endblk, jk, nlev, nlevp1

    !0) Initialize all passed ddt's to 0 
!$OMP PARALLEL WORKSHARE
    ddt_u_ls  = 0._wp; ddt_v_ls = 0._wp;  ddt_temp_ls = 0._wp
    ddt_qv_ls = 0._wp; ddt_ql_ls = 0._wp; ddt_qi_ls  = 0._wp
!$OMP END PARALLEL WORKSHARE
   
    i_nchdom  = MAX(1,p_patch%n_childdom)
    nlev      = p_patch%nlev
    nlevp1    = p_patch%nlev+1
    inv_no_gb_cells = 1._wp / REAL(p_patch%n_patch_cells_g,wp)

    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)

    !Some precalculations   
    IF(is_subsidence_heat .OR. is_advection .OR. is_rad_forcing) &
      CALL global_hor_mean(p_patch, p_prog%exner(:,:,:), exner_gb, inv_no_gb_cells, i_nchdom)

    IF(is_subsidence_moment.OR.is_subsidence_heat)THEN
!$OMP PARALLEL WORKSHARE
      !Get surface density
      rhos(:,i_startblk:i_endblk) = p_diag%pres_sfc(:,i_startblk:i_endblk) / (rd * &
                      t_sfc(:,i_startblk:i_endblk)*(1._wp+0.61_wp*qv_sfc(:,i_startblk:i_endblk)) )  
!$OMP END PARALLEL WORKSHARE

      !rho_gb(nlev)
      CALL global_hor_mean(p_patch, p_prog%rho(:,:,:), rho_gb, inv_no_gb_cells, i_nchdom)
      inv_rho_gb  = 1._wp / rho_gb

      !Vertical advective forcing
      inv_dz(:) = 1._wp / p_metrics%ddqz_z_full(2,:,2)        
    END IF

    !1a) Horizontal mean of variables and their vertical advective tendency - momentum

    IF(is_subsidence_moment)THEN
!$OMP PARALLEL WORKSHARE
      !rho*u
      varin(:,:,i_startblk:i_endblk) = p_diag%u(:,:,i_startblk:i_endblk)*p_prog%rho(:,:,i_startblk:i_endblk)
!$OMP END PARALLEL WORKSHARE
      CALL vert_intp_full2half_cell_3d(p_patch, p_metrics, varin, varout, rl_start, rl_end)
      CALL global_hor_mean(p_patch, varout, u_gb, inv_no_gb_cells, i_nchdom)

!$OMP PARALLEL WORKSHARE
      varin(:,:,i_startblk:i_endblk) = p_diag%v(:,:,i_startblk:i_endblk)*p_prog%rho(:,:,i_startblk:i_endblk)
!$OMP END PARALLEL WORKSHARE
      CALL vert_intp_full2half_cell_3d(p_patch, p_metrics, varin, varout, rl_start, rl_end)
      CALL global_hor_mean(p_patch, varout, v_gb, inv_no_gb_cells, i_nchdom)
        
      ddt_u_ls    =  ddt_u_ls - w_ls*vertical_derivative(u_gb,inv_dz)*inv_rho_gb
      ddt_v_ls    =  ddt_v_ls - w_ls*vertical_derivative(v_gb,inv_dz)*inv_rho_gb
    END IF

    !1b) Horizontal mean of variables and their vertical advective tendency 

    IF(is_subsidence_heat)THEN
      !use theta instead of temperature for subsidence
!$OMP PARALLEL WORKSHARE
      varin(:,:,i_startblk:i_endblk) = p_diag%temp(:,:,i_startblk:i_endblk) * &
                         p_prog%rho(:,:,i_startblk:i_endblk) / p_prog%exner(:,:,i_startblk:i_endblk)
!$OMP END PARALLEL WORKSHARE
      CALL vert_intp_full2half_cell_3d(p_patch, p_metrics, varin, varout, rl_start, rl_end)
!$OMP PARALLEL WORKSHARE
      !assume exner at surface ==1
      varout(:,nlevp1,i_startblk:i_endblk) = t_sfc(:,i_startblk:i_endblk)*rhos(:,i_startblk:i_endblk)
!$OMP END PARALLEL WORKSHARE
      CALL global_hor_mean(p_patch, varout, temp_gb, inv_no_gb_cells, i_nchdom)

!$OMP PARALLEL WORKSHARE
      varin(:,:,i_startblk:i_endblk) = qv(:,:,i_startblk:i_endblk)*p_prog%rho(:,:,i_startblk:i_endblk)
!$OMP END PARALLEL WORKSHARE
      CALL vert_intp_full2half_cell_3d(p_patch, p_metrics, varin, varout, rl_start, rl_end)
!$OMP PARALLEL WORKSHARE
      varout(:,nlevp1,i_startblk:i_endblk) = qv_sfc(:,i_startblk:i_endblk)*rhos(:,i_startblk:i_endblk)
!$OMP END PARALLEL WORKSHARE
      CALL global_hor_mean(p_patch, varout, qv_gb, inv_no_gb_cells, i_nchdom)

!$OMP PARALLEL WORKSHARE
      varin(:,:,i_startblk:i_endblk) = ql(:,:,i_startblk:i_endblk)*p_prog%rho(:,:,i_startblk:i_endblk)
!$OMP END PARALLEL WORKSHARE
      CALL vert_intp_full2half_cell_3d(p_patch, p_metrics, varin, varout, rl_start, rl_end)
!$OMP PARALLEL WORKSHARE
      varout(:,nlevp1,i_startblk:i_endblk) = 0._wp
!$OMP END PARALLEL WORKSHARE
      CALL global_hor_mean(p_patch, varout, ql_gb, inv_no_gb_cells, i_nchdom)

      ddt_temp_ls(1:nlev) =  ddt_temp_ls(1:nlev) - w_ls(1:nlev)*vertical_derivative(temp_gb,inv_dz) * &
                                        inv_rho_gb(1:nlev) / exner_gb(1:nlev)

      ddt_qv_ls   =  ddt_qv_ls - w_ls*vertical_derivative(qv_gb,inv_dz)*inv_rho_gb
      ddt_ql_ls   =  ddt_ql_ls - w_ls*vertical_derivative(ql_gb,inv_dz)*inv_rho_gb
    END IF

    !2)Horizontal advective forcing: at present only for tracers
    !  Advective tendencies for ql is always assumed 0
    IF(is_advection)THEN
      IF(is_theta)THEN
        ddt_temp_ls = ddt_temp_ls + ddt_temp_hadv_ls * exner_gb      
      ELSE
        ddt_temp_ls = ddt_temp_ls + ddt_temp_hadv_ls 
      END IF
      ddt_qv_ls   = ddt_qv_ls   + ddt_qv_hadv_ls
    END IF

    !3)Coriolis and geostrophic wind
    IF(is_geowind)THEN
      !Remember model grid is flat   
      ddt_u_ls = ddt_u_ls - p_patch%cells%f_c(2,2) * v_geo
      ddt_v_ls = ddt_v_ls + p_patch%cells%f_c(2,2) * u_geo
    END IF

    !4)Radiative forcing
    IF(is_rad_forcing)THEN
      IF(is_theta)THEN
        ddt_temp_ls = ddt_temp_ls + ddt_temp_rad_ls * exner_gb
      ELSE
        ddt_temp_ls = ddt_temp_ls + ddt_temp_rad_ls 
      END IF
    END IF

  END SUBROUTINE apply_ls_forcing


  !>
  !! global_hor_mean: only called for interior points
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-May-30)
  SUBROUTINE global_hor_mean(p_patch, var, varout, inv_no_cells, nchdom)

    TYPE(t_patch),     INTENT(in), TARGET :: p_patch
    REAL(wp), INTENT(in)                  :: var(:,:,:), inv_no_cells
    INTEGER,  INTENT(in)                  :: nchdom
    REAL(wp), INTENT(out)                 :: varout(:)                     

    REAL(wp) :: var_aux(SIZE(var,1),SIZE(var,2),SIZE(var,3))
    INTEGER  :: i_startblk, i_endblk, rl_start
    INTEGER  :: i_endidx, i_startidx
    INTEGER  :: jk, jc, jb, nz

    !Put all fields to 0
    var_aux(:,:,:) = 0._wp

    rl_start   = grf_bdywidth_c+1
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_endblk   = p_patch%cells%end_blk(min_rlcell_int,nchdom)
    nz         = SIZE(var,2)

   !Now put values in interior nodes
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
       CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, min_rlcell_int)
       DO jk = 1 , nz
         DO jc = i_startidx , i_endidx
             var_aux(jc,jk,jb) = var(jc,jk,jb)
         END DO
       END DO
    END DO 
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   DO jk = 1 , nz
    varout(jk) = global_sum_array(var_aux(:,jk,:)) * inv_no_cells
   END DO


  END SUBROUTINE global_hor_mean

  !>
  !! vertical_derivative
  !!------------------------------------------------------------------------
  !! @par Revision History
  !! Initial release by Anurag Dipankar, MPI-M (2013-May-30)
  FUNCTION vertical_derivative (var, inv_dz) RESULT(dvardz)

    REAL(wp), INTENT(in) :: var(:), inv_dz(:)
                     
    REAL(wp) :: dvardz(SIZE(inv_dz))                     
    INTEGER  :: jk

!$OMP PARALLEL
!$OMP DO PRIVATE(jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jk = 1 , SIZE(inv_dz)
      dvardz(jk) = ( var(jk) - var(jk+1) ) * inv_dz(jk)
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END FUNCTION vertical_derivative

!-------------------------------------------------------------------------------
     
END MODULE mo_ls_forcing



