!>
!!  This module contains (not yet) utility programs for boundary interpolation and feedback
!!  of diagnostic variables and the upscaling and downscaling routines needed
!!  for the reduced physics grid
!!
!! @par Revision History
!!  Developed and tested by Guenther Zaengl, DWD (2010-02-10)
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!!

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_phys_nest_utilities
!
!
USE mo_kind,                ONLY: wp
USE mo_exception,           ONLY: message_text, message
USE mo_model_domain,        ONLY: t_patch, t_grid_cells, p_patch_local_parent, p_patch
USE mo_intp_data_strc,      ONLY: t_int_state, p_int_state_local_parent
USE mo_grf_intp_data_strc,  ONLY: t_gridref_state, t_gridref_single_state, &
                                  p_grf_state, p_grf_state_local_parent
USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag
USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_lnd_diag, t_wtr_prog
USE mo_ext_data_types,      ONLY: t_external_data
USE mo_nwp_lnd_state,       ONLY: p_lnd_state
USE mo_grf_bdyintp,         ONLY: interpol_scal_grf
USE mo_grf_nudgintp,        ONLY: interpol_scal_nudging
USE mo_parallel_config,     ONLY: nproma, p_test_run
USE mo_dynamics_config,     ONLY: nnow_rcf
USE mo_run_config,          ONLY: msg_level, iqv, iqc, iqi
USE mo_grid_config,         ONLY: l_limited_area
USE mo_nwp_phy_state,       ONLY: prm_diag
USE mo_nonhydro_state,      ONLY: p_nh_state
USE mo_impl_constants,      ONLY: min_rlcell, min_rlcell_int, nexlevs_rrg_vnest, dzsoil
USE mo_physical_constants,  ONLY: rd, grav, stbo, vtmpc1, tmelt
USE mo_satad,               ONLY: qsat_rho
USE mo_loopindices,         ONLY: get_indices_c
USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_ovlparea_start_c, grf_fbk_start_c
USE mo_vertical_coord_table,ONLY: vct_a
USE mo_communication,       ONLY: exchange_data, exchange_data_mult
USE mo_sync,                ONLY: SYNC_C, sync_patch_array, sync_patch_array_mult
USE mo_lnd_nwp_config,      ONLY: nlev_soil, nlev_snow, lmulti_snow, lseaice, llake, &
                                  frlake_thrhld, frsea_thrhld, isub_lake, ntiles_total
USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
USE mo_phyparam_soil,       ONLY: cadp
USE mo_mpi,                 ONLY: my_process_is_mpi_seq
USE mo_seaice_nwp,          ONLY: frsi_min
USE mo_flake,               ONLY: flake_coldinit
USE mo_data_flake,          ONLY: tpl_T_r, C_T_min, rflk_depth_bs_ref
USE mo_fortran_tools,       ONLY: init, copy

IMPLICIT NONE

PRIVATE


PUBLIC :: upscale_rad_input        ! for RRTM
PUBLIC :: upscale_rad_input_rg     ! for Ritter-Geleyn
PUBLIC :: downscale_rad_output     ! for RRTM
PUBLIC :: downscale_rad_output_rg  ! for Ritter-Geleyn
PUBLIC :: interpol_phys_grf
PUBLIC :: feedback_phys_diag
PUBLIC :: interpol_rrg_grf
PUBLIC :: copy_rrg_ubc

CONTAINS

!>
!! This routine averages the input fields for RRTM radiation to the next coarser grid level.
!!
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2010-12-01
!!
SUBROUTINE upscale_rad_input(jg, jgp, nlev_rg, fr_land, fr_glac, emis_rad, &
  cosmu0, albvisdir, albnirdir, albvisdif, albnirdif, albdif,              &
  tsfc, ktype, pres_ifc, pres, temp, acdnc, tot_cld, clc, q_o3,            &
  aeq1, aeq2, aeq3, aeq4, aeq5,                                            &
  rg_fr_land, rg_fr_glac, rg_emis_rad,                                     &
  rg_cosmu0, rg_albvisdir, rg_albnirdir, rg_albvisdif, rg_albnirdif,       &
  rg_albdif, rg_tsfc, rg_rtype, rg_pres_ifc, rg_pres, rg_temp, rg_acdnc,   &
  rg_tot_cld, rg_clc, rg_q_o3, rg_aeq1, rg_aeq2, rg_aeq3, rg_aeq4, rg_aeq5,&
  z_pres_ifc, z_tot_cld, buffer_rrg                                        )

  ! Input grid parameters
  INTEGER, INTENT(IN)  :: jg, jgp  ! domain IDs of main and reduced grids
  INTEGER, INTENT(IN)  :: nlev_rg  ! number of model levels on reduced grid

  ! Input fields (on full grid)
  REAL(wp), INTENT(IN) ::                                                                 &
    fr_land(:,:), fr_glac(:,:), emis_rad(:,:),                                            &
    cosmu0(:,:), albvisdir(:,:), albnirdir(:,:), albvisdif(:,:), albnirdif(:,:),          &
    albdif(:,:), tsfc(:,:), pres_ifc(:,:,:), pres(:,:,:), temp(:,:,:), acdnc(:,:,:),      &
    tot_cld(:,:,:,:), clc(:,:,:), q_o3(:,:,:), aeq1(:,:,:), aeq2(:,:,:), aeq3(:,:,:),     &
    aeq4(:,:,:), aeq5(:,:,:)

  INTEGER, INTENT(IN) :: ktype(:,:)

  ! Corresponding output fields (on reduced grid)
  REAL(wp), TARGET, INTENT(OUT) ::                                           &
    rg_fr_land(:,:),rg_fr_glac(:,:), rg_emis_rad(:,:),                       &
    rg_cosmu0(:,:), rg_albvisdir(:,:), rg_albnirdir(:,:), rg_albvisdif(:,:), &
    rg_albnirdif(:,:), rg_albdif(:,:), rg_tsfc(:,:), rg_rtype(:,:),          &
    rg_pres_ifc(:,:,:), rg_pres(:,:,:), rg_temp(:,:,:), rg_acdnc(:,:,:),     &
    rg_tot_cld(:,:,:,:), rg_clc(:,:,:), rg_q_o3(:,:,:), rg_aeq1(:,:,:),      &
    rg_aeq2(:,:,:), rg_aeq3(:,:,:), rg_aeq4(:,:,:),rg_aeq5(:,:,:),           &
    ! these have the same function as the intermediate storage fields below but are passed to the calling routine
    z_pres_ifc(:,:,:), z_tot_cld(:,:,:,:)


  ! Intermediate storage fields needed in the case of MPI parallelization
  REAL(wp), ALLOCATABLE, TARGET ::                                             &
    z_fr_land(:,:),z_fr_glac(:,:), z_emis_rad(:,:),                            &
    z_cosmu0(:,:), z_albvisdir(:,:), z_albnirdir(:,:), z_albvisdif(:,:),       &
    z_albnirdif(:,:), z_albdif(:,:), z_tsfc(:,:), z_rtype(:,:),                &
    z_pres(:,:,:), z_temp(:,:,:), z_acdnc(:,:,:), z_clc(:,:,:), z_q_o3(:,:,:), &
    z_aeq1(:,:,:), z_aeq2(:,:,:), z_aeq3(:,:,:), z_aeq4(:,:,:), z_aeq5(:,:,:), &
    z_aux3d(:,:,:), zrg_aux3d(:,:,:)


  ! Pointers to output fields (no MPI) or intermediate fields (MPI)
  REAL(wp), POINTER ::                                                   &
    p_fr_land(:,:),p_fr_glac(:,:), p_emis_rad(:,:),                      &
    p_cosmu0(:,:), p_albvisdir(:,:), p_albnirdir(:,:), p_albvisdif(:,:), &
    p_albnirdif(:,:), p_albdif(:,:), p_tsfc(:,:), p_rtype(:,:),          &
    p_pres_ifc(:,:,:), p_pres(:,:,:), p_temp(:,:,:), p_acdnc(:,:,:),     &
    p_tot_cld(:,:,:,:), p_clc(:,:,:), p_q_o3(:,:,:), p_aeq1(:,:,:),      &
    p_aeq2(:,:,:), p_aeq3(:,:,:), p_aeq4(:,:,:),p_aeq5(:,:,:)

  ! Buffer for auxiliary temperature and pressure levels above the vertical nest interface
  REAL(wp), INTENT(IN) :: buffer_rrg(:,:,:)

  ! Pointers to types needed to minimize code duplication for MPI/no-MPI cases
  TYPE(t_grid_cells), POINTER     :: p_gcp
  TYPE(t_gridref_state), POINTER  :: p_grf
  TYPE(t_patch),      POINTER     :: p_pp

  ! Indices
  INTEGER :: jb, jc, jk, jk1, i_chidx, i_nchdom, i_startrow, &
             i_startblk, i_endblk, i_startidx, i_endidx, nblks_c_lp

  INTEGER :: nlev, nlevp1      !< number of full and half levels
  INTEGER :: nshift, nlevp1_rg, nst, nstart, nend
  REAL(wp) :: exdist_h, exdist_f

  INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
  REAL(wp), POINTER :: p_fbkwgt(:,:,:)

  !-----------------------------------------------------------------------

  IF (msg_level >= 10) THEN
    WRITE(message_text,'(a,i2,a,i2)') 'Upscaling of radiation input fields',&
      jg,' =>',jgp
    CALL message('upscale_rad_input',message_text)
  ENDIF

  ! Number of levels of the full grid
  nlev   = p_patch(jg)%nlev
  nlevp1 = p_patch(jg)%nlevp1

  ! layer shift w.r.t. global grid (> 0 in case of vertical nesting)
  nst = p_patch(jg)%nshift_total

  ! nlev_rg carries the number of model levels of the reduced grid,
  ! which may be larger than nlev
  nlevp1_rg = nlev_rg + 1
  nshift = nlev_rg - nlev ! resulting shift parameter

  ! Parameters used in case of latm_above_top = .TRUE.:
  !
  ! start and end levels for parent to local parent copying in case of vertical nesting
  nstart = p_patch(jgp)%nlev - nlev_rg + 1
  nend   = p_patch(jgp)%nlev - nlev_rg + nshift
  !
  ! extrapolation distances for passive layer above model top (m) if there is no vertical nesting
  exdist_h = 1.5_wp*(vct_a(nst+1)-vct_a(nst+2))
  exdist_f = vct_a(nst+1)+0.5_wp*exdist_h - 0.5_wp*(vct_a(nst+1)+vct_a(nst+2))

  p_grf => p_grf_state_local_parent(jg)
  p_gcp => p_patch_local_parent(jg)%cells
  p_pp  => p_patch_local_parent(jg)

  i_chidx  = p_patch(jg)%parent_child_index
  i_nchdom = MAX(1,p_patch(jgp)%n_childdom)

  ! Set pointers to index and coefficient fields for cell-based variables
  iidx => p_gcp%child_idx
  iblk => p_gcp%child_blk

  p_fbkwgt => p_grf%fbk_wgt_bln

  ! Allocation of local storage fields at local parent level in MPI-case
  IF (jgp == 0 .AND. .NOT. l_limited_area) THEN
    nblks_c_lp = p_gcp%end_blk(min_rlcell,i_chidx)

    ALLOCATE(z_fr_land(nproma,nblks_c_lp), z_fr_glac(nproma,nblks_c_lp),                &
             z_emis_rad(nproma,nblks_c_lp),                                             &
             z_cosmu0(nproma,nblks_c_lp), z_albvisdir(nproma,nblks_c_lp),               &
             z_albnirdir(nproma,nblks_c_lp), z_albvisdif(nproma,nblks_c_lp),            &
             z_albnirdif(nproma,nblks_c_lp), z_albdif(nproma,nblks_c_lp),               &
             z_tsfc(nproma,nblks_c_lp), z_rtype(nproma,nblks_c_lp),                     &
             z_pres(nproma,nlev_rg,nblks_c_lp), z_temp(nproma,nlev_rg,nblks_c_lp),      &
             z_acdnc(nproma,nlev_rg,nblks_c_lp), z_clc(nproma,nlev_rg,nblks_c_lp),      &
             z_q_o3(nproma,nlev_rg,nblks_c_lp), z_aeq1(nproma,nlev_rg,nblks_c_lp),      &
             z_aeq2(nproma,nlev_rg,nblks_c_lp), z_aeq3(nproma,nlev_rg,nblks_c_lp),      &
             z_aeq4(nproma,nlev_rg,nblks_c_lp), z_aeq5(nproma,nlev_rg,nblks_c_lp),      &
             z_aux3d(nproma,11,nblks_c_lp), zrg_aux3d(nproma,11,p_patch(jgp)%nblks_c) )

    ! Set pointers to either the parent-level variables (non-MPI case) or to the
    ! intermediate storage fields (MPI case)
    p_fr_land    => z_fr_land
    p_fr_glac    => z_fr_glac
    p_emis_rad   => z_emis_rad
    p_cosmu0     => z_cosmu0
    p_albvisdir  => z_albvisdir
    p_albnirdir  => z_albnirdir
    p_albvisdif  => z_albvisdif
    p_albnirdif  => z_albnirdif
    p_albdif     => z_albdif
    p_tsfc       => z_tsfc
    p_rtype      => z_rtype
    p_pres_ifc   => z_pres_ifc
    p_pres       => z_pres
    p_temp       => z_temp
    p_acdnc      => z_acdnc
    p_tot_cld    => z_tot_cld
    p_clc        => z_clc
    p_q_o3       => z_q_o3
    p_aeq1       => z_aeq1
    p_aeq2       => z_aeq2
    p_aeq3       => z_aeq3
    p_aeq4       => z_aeq4
    p_aeq5       => z_aeq5
  ELSE
    p_fr_land    => rg_fr_land
    p_fr_glac    => rg_fr_glac
    p_emis_rad   => rg_emis_rad
    p_cosmu0     => rg_cosmu0
    p_albvisdir  => rg_albvisdir
    p_albnirdir  => rg_albnirdir
    p_albvisdif  => rg_albvisdif
    p_albnirdif  => rg_albnirdif
    p_albdif     => rg_albdif
    p_tsfc       => rg_tsfc
    p_rtype      => rg_rtype
    p_pres_ifc   => rg_pres_ifc
    p_pres       => rg_pres
    p_temp       => rg_temp
    p_acdnc      => rg_acdnc
    p_tot_cld    => rg_tot_cld
    p_clc        => rg_clc
    p_q_o3       => rg_q_o3
    p_aeq1       => rg_aeq1
    p_aeq2       => rg_aeq2
    p_aeq3       => rg_aeq3
    p_aeq4       => rg_aeq4
    p_aeq5       => rg_aeq5
  ENDIF

  IF (p_test_run) THEN
    p_fr_land    = 0._wp
    p_fr_glac    = 0._wp
    p_emis_rad   = 0._wp
    p_cosmu0     = 0._wp
    p_albvisdir  = 0._wp
    p_albnirdir  = 0._wp
    p_albvisdif  = 0._wp
    p_albnirdif  = 0._wp
    p_albdif     = 0._wp
    p_tsfc       = 0._wp
    p_rtype      = 0._wp
    p_pres_ifc   = 0._wp
    p_pres       = 0._wp
    p_temp       = 0._wp
    p_acdnc      = 0._wp
    p_tot_cld    = 0._wp
    p_clc        = 0._wp
    p_q_o3       = 0._wp
    p_aeq1       = 0._wp
    p_aeq2       = 0._wp
    p_aeq3       = 0._wp
    p_aeq4       = 0._wp
    p_aeq5       = 0._wp
  ENDIF


  ! Now average input fields to parent grid cells

  ! Start/End block in the parent domain
  IF (jg == 1 .AND. l_limited_area) THEN
    i_startrow = grf_fbk_start_c
  ELSE
    i_startrow = grf_ovlparea_start_c
  ENDIF
  i_startblk = p_gcp%start_blk(i_startrow,i_chidx)
  i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jk1) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,                           &
                       i_startidx, i_endidx, i_startrow, min_rlcell_int)
!DIR$ IVDEP
    DO jc = i_startidx, i_endidx

      p_fr_land(jc,jb) =                                         &
        fr_land(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        fr_land(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        fr_land(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        fr_land(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_fr_glac(jc,jb) =                                         &
        fr_glac(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        fr_glac(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        fr_glac(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        fr_glac(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_emis_rad(jc,jb) =                                         &
        emis_rad(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        emis_rad(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        emis_rad(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        emis_rad(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_cosmu0(jc,jb) =                                         &
        cosmu0(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        cosmu0(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        cosmu0(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        cosmu0(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      ! Limiting positive values of the cosine of the zenith angle to a minimum
      ! of 0.05 is taken over from SR pre_radiation_nwp_steps. This avoids numerical
      ! trouble along the day-night boundary near the model top
      IF (p_cosmu0(jc,jb) > 1.e-5_wp) p_cosmu0(jc,jb) = MAX(0.05_wp,p_cosmu0(jc,jb))

      p_albvisdir(jc,jb) =                                         &
        albvisdir(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        albvisdir(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        albvisdir(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        albvisdir(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_albnirdir(jc,jb) =                                         &
        albnirdir(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        albnirdir(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        albnirdir(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        albnirdir(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_albvisdif(jc,jb) =                                         &
        albvisdif(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        albvisdif(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        albvisdif(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        albvisdif(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_albnirdif(jc,jb) =                                         &
        albnirdif(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        albnirdif(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        albnirdif(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        albnirdif(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_albdif(jc,jb) =                                         &
        albdif(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        albdif(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        albdif(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        albdif(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_tsfc(jc,jb) = SQRT(SQRT(                                 &
        tsfc(iidx(jc,jb,1),iblk(jc,jb,1))**4*p_fbkwgt(jc,jb,1) + &
        tsfc(iidx(jc,jb,2),iblk(jc,jb,2))**4*p_fbkwgt(jc,jb,2) + &
        tsfc(iidx(jc,jb,3),iblk(jc,jb,3))**4*p_fbkwgt(jc,jb,3) + &
        tsfc(iidx(jc,jb,4),iblk(jc,jb,4))**4*p_fbkwgt(jc,jb,4) ) )

      p_rtype(jc,jb) =                                        &
        REAL(ktype(iidx(jc,jb,1),iblk(jc,jb,1)),wp)*p_fbkwgt(jc,jb,1) + &
        REAL(ktype(iidx(jc,jb,2),iblk(jc,jb,2)),wp)*p_fbkwgt(jc,jb,2) + &
        REAL(ktype(iidx(jc,jb,3),iblk(jc,jb,3)),wp)*p_fbkwgt(jc,jb,3) + &
        REAL(ktype(iidx(jc,jb,4),iblk(jc,jb,4)),wp)*p_fbkwgt(jc,jb,4)

      p_pres_ifc(jc,nlevp1_rg,jb) =                                      &
        pres_ifc(iidx(jc,jb,1),nlevp1,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        pres_ifc(iidx(jc,jb,2),nlevp1,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        pres_ifc(iidx(jc,jb,3),nlevp1,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        pres_ifc(iidx(jc,jb,4),nlevp1,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

    ENDDO

    IF (jgp == 0 .AND. .NOT. l_limited_area) THEN ! combine 2D fields in a 3D field to speed up MPI communication
      DO jc = i_startidx, i_endidx
        z_aux3d(jc, 1,jb) = p_cosmu0(jc,jb)
        z_aux3d(jc, 2,jb) = p_albvisdir(jc,jb)
        z_aux3d(jc, 3,jb) = p_albnirdir(jc,jb)
        z_aux3d(jc, 4,jb) = p_albvisdif(jc,jb)
        z_aux3d(jc, 5,jb) = p_albnirdif(jc,jb)
        z_aux3d(jc, 6,jb) = p_albdif(jc,jb)
        z_aux3d(jc, 7,jb) = p_tsfc(jc,jb)
        z_aux3d(jc, 8,jb) = p_fr_land(jc,jb)
        z_aux3d(jc, 9,jb) = p_fr_glac(jc,jb)
        z_aux3d(jc,10,jb) = p_emis_rad(jc,jb)
        z_aux3d(jc,11,jb) = p_rtype(jc,jb)
      ENDDO
    ENDIF

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
!DIR$ IVDEP
      DO jk = 1, nlev
        jk1 = jk + nshift
#else
    DO jk = 1, nlev
      jk1 = jk + nshift
      DO jc = i_startidx, i_endidx
#endif
        p_pres_ifc(jc,jk1,jb) =                                        &
          pres_ifc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          pres_ifc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          pres_ifc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          pres_ifc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_pres(jc,jk1,jb) =                                        &
          pres(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          pres(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          pres(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          pres(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_temp(jc,jk1,jb) =                                        &
          temp(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          temp(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          temp(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          temp(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_q_o3(jc,jk1,jb) =                                        &
          q_o3(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          q_o3(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          q_o3(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          q_o3(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_aeq1(jc,jk1,jb) =                                        &
          aeq1(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          aeq1(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          aeq1(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          aeq1(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_aeq2(jc,jk1,jb) =                                        &
          aeq2(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          aeq2(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          aeq2(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          aeq2(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_aeq3(jc,jk1,jb) =                                        &
          aeq3(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          aeq3(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          aeq3(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          aeq3(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_aeq4(jc,jk1,jb) =                                        &
          aeq4(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          aeq4(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          aeq4(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          aeq4(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_aeq5(jc,jk1,jb) =                                        &
          aeq5(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          aeq5(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          aeq5(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          aeq5(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      ENDDO
    ENDDO

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
!DIR$ IVDEP
      DO jk = 1, nlev
        jk1 = jk + nshift
#else
    DO jk = 1, nlev
      jk1 = jk + nshift
      DO jc = i_startidx, i_endidx
#endif

        p_acdnc(jc,jk1,jb) =                                        &
          acdnc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          acdnc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          acdnc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          acdnc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_clc(jc,jk1,jb) =                                        &
          clc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          clc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          clc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          clc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

!CDIR EXPAND=3
        p_tot_cld(jc,jk1,jb,1:3) =                                        &
          tot_cld(iidx(jc,jb,1),jk,iblk(jc,jb,1),1:3)*p_fbkwgt(jc,jb,1) + &
          tot_cld(iidx(jc,jb,2),jk,iblk(jc,jb,2),1:3)*p_fbkwgt(jc,jb,2) + &
          tot_cld(iidx(jc,jb,3),jk,iblk(jc,jb,3),1:3)*p_fbkwgt(jc,jb,3) + &
          tot_cld(iidx(jc,jb,4),jk,iblk(jc,jb,4),1:3)*p_fbkwgt(jc,jb,4)

      ENDDO
    ENDDO

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = 1, nlev
        jk1 = jk + nshift
#else
    DO jk = 1, nlev
      jk1 = jk + nshift
      DO jc = i_startidx, i_endidx
#endif

        ! enhance averaged QC and QI in order to be more consistent with cloud cover scheme
        IF (p_clc(jc,jk1,jb) > 0._wp .AND. p_clc(jc,jk1,jb) < 0.95_wp) THEN
!CDIR EXPAND=2
          p_tot_cld(jc,jk1,jb,2:3) =  0.5_wp*(p_tot_cld(jc,jk1,jb,2:3) + SQRT( &
            tot_cld(iidx(jc,jb,1),jk,iblk(jc,jb,1),2:3)**2*p_fbkwgt(jc,jb,1) + &
            tot_cld(iidx(jc,jb,2),jk,iblk(jc,jb,2),2:3)**2*p_fbkwgt(jc,jb,2) + &
            tot_cld(iidx(jc,jb,3),jk,iblk(jc,jb,3),2:3)**2*p_fbkwgt(jc,jb,3) + &
            tot_cld(iidx(jc,jb,4),jk,iblk(jc,jb,4),2:3)**2*p_fbkwgt(jc,jb,4)  ))
        ENDIF

      ENDDO
    ENDDO

    IF (nshift > 0) THEN ! set values for extra layer(s) above the top of the computational model grid
      !
      ! assume zero-gradient condition for aerosols, ozone and clouds (the latter are zero anyway in practice)
      jk1 = nshift + 1
      DO jk = 1, nshift
        DO jc = i_startidx, i_endidx
          p_aeq1(jc,jk,jb) = p_aeq1(jc,jk1,jb)
          p_aeq2(jc,jk,jb) = p_aeq2(jc,jk1,jb)
          p_aeq3(jc,jk,jb) = p_aeq3(jc,jk1,jb)
          p_aeq4(jc,jk,jb) = p_aeq4(jc,jk1,jb)
          p_aeq5(jc,jk,jb) = p_aeq5(jc,jk1,jb)
          p_acdnc(jc,jk,jb) = p_acdnc(jc,jk1,jb)
          p_tot_cld(jc,jk,jb,1:3) = p_tot_cld(jc,jk1,jb,1:3)
          p_clc (jc,jk,jb) = p_clc(jc,jk1,jb)
          p_q_o3(jc,jk,jb) = p_q_o3(jc,jk1,jb)
        ENDDO
      ENDDO

      IF (jgp == 0 .OR. p_patch(jg)%nshift == 0) THEN ! settings for passive extra layer above model top for global grid (nshift=1 in this case)
        DO jc = i_startidx, i_endidx
          ! Temperature is extrapolated linearly assuming a vertical temperature gradient of -5.0 K/km
          p_temp(jc,1,jb) = p_temp(jc,2,jb) - 5.0e-3_wp*exdist_f
          p_pres_ifc(jc,1,jb) = p_pres_ifc(jc,2,jb)*EXP(-grav*exdist_h/(rd*p_temp(jc,1,jb)))
          p_pres(jc,1,jb) = SQRT(p_pres_ifc(jc,1,jb)*p_pres_ifc(jc,2,jb))
        ENDDO
      ELSE ! get information from buffer field
        DO jk = 1, nshift
          DO jc = i_startidx, i_endidx
            rg_pres_ifc(jc,jk,jb)   = buffer_rrg(jc,jk,jb)
            rg_pres    (jc,jk,jb)   = buffer_rrg(jc,nshift+jk,jb)
            rg_temp    (jc,jk,jb)   = buffer_rrg(jc,2*nshift+jk,jb)
          ENDDO
        ENDDO
      ENDIF

    ENDIF

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  IF (jgp == 0 .AND. .NOT. l_limited_area) THEN
    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 6, 5*nlev_rg+12, &
                            RECV1=rg_pres_ifc, SEND1=z_pres_ifc,             &
                            RECV2=rg_pres,     SEND2=z_pres,                 &
                            RECV3=rg_temp,     SEND3=z_temp,                 &
                            RECV4=rg_acdnc,    SEND4=z_acdnc,                &
                            RECV5=zrg_aux3d,   SEND5=z_aux3d,                &
                            RECV6=rg_q_o3,     SEND6=z_q_o3                  )

    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 9, 9*nlev_rg, &
                            RECV1=rg_aeq1,     SEND1=z_aeq1,              &
                            RECV2=rg_aeq2,     SEND2=z_aeq2,              &
                            RECV3=rg_aeq3,     SEND3=z_aeq3,              &
                            RECV4=rg_aeq4,     SEND4=z_aeq4,              &
                            RECV5=rg_aeq5,     SEND5=z_aeq5,              &
                            RECV6=rg_clc,      SEND6=z_clc,               &
                            RECV4D=rg_tot_cld, SEND4D=z_tot_cld           )


    i_startblk = p_patch(jgp)%cells%start_blk(1,1)
    i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch(jgp), jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rlcell)

      DO jc = i_startidx, i_endidx
        rg_cosmu0(jc,jb)    = zrg_aux3d(jc,1,jb)
        rg_albvisdir(jc,jb) = zrg_aux3d(jc,2,jb)
        rg_albnirdir(jc,jb) = zrg_aux3d(jc,3,jb)
        rg_albvisdif(jc,jb) = zrg_aux3d(jc,4,jb)
        rg_albnirdif(jc,jb) = zrg_aux3d(jc,5,jb)
        rg_albdif(jc,jb)    = zrg_aux3d(jc,6,jb)
        rg_tsfc(jc,jb)      = zrg_aux3d(jc,7,jb)
        rg_fr_land(jc,jb)   = zrg_aux3d(jc,8,jb)
        rg_fr_glac(jc,jb)   = zrg_aux3d(jc,9,jb)
        rg_emis_rad(jc,jb)  = zrg_aux3d(jc,10,jb)
        rg_rtype(jc,jb)     = zrg_aux3d(jc,11,jb)
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    DEALLOCATE(z_fr_land, z_fr_glac, z_emis_rad, z_cosmu0, z_albvisdir, z_albnirdir, &
      & z_albvisdif, z_albnirdif, z_albdif, z_tsfc, z_rtype, z_pres, z_temp, z_acdnc,&
      & z_clc, z_q_o3, z_aeq1, z_aeq2, z_aeq3, z_aeq4, z_aeq5, z_aux3d, zrg_aux3d )

  ENDIF

END SUBROUTINE upscale_rad_input

!>
!! This routine interpolates the output fields of RRTM radiation from the reduced
!! grid to the full grid.
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2010-12-03
!!
SUBROUTINE downscale_rad_output(jg, jgp, nlev_rg,                           &
  rg_aclcov, rg_lwflxall, rg_trsolall, rg_trsol_clr_sfc, rg_lwflx_up_sfc,   &
  rg_trsol_up_toa, rg_trsol_up_sfc, rg_trsol_par_sfc, rg_trsol_dn_sfc_diff, &
  tsfc_rg, albdif_rg, emis_rad_rg, cosmu0_rg, tot_cld_rg, z_tot_cld,        &
  pres_ifc_rg, z_pres_ifc, tsfc, albdif, aclcov, lwflxall, trsolall,        &
  lwflx_up_sfc, trsol_up_toa, trsol_up_sfc, trsol_par_sfc,                  &
  trsol_dn_sfc_diff, trsol_clr_sfc)


  ! Input grid parameters
  INTEGER, INTENT(IN)  :: jg, jgp  ! domain IDs of main and reduced grids
  INTEGER, INTENT(IN)  :: nlev_rg  ! number of model levels on reduced grid

  ! Input fields (on reduced grid) to be downscaled to full grid
  REAL(wp), TARGET, INTENT(IN) ::                                                              &
    rg_aclcov(:,:), rg_lwflxall(:,:,:), rg_trsolall(:,:,:), rg_lwflx_up_sfc(:,:),              &
    rg_trsol_up_toa(:,:), rg_trsol_up_sfc(:,:), rg_trsol_par_sfc(:,:), rg_trsol_dn_sfc_diff(:,:)

  ! Auxiliary input fields on reduced grid needed for downscaling corrections
  REAL(wp), INTENT(IN), TARGET ::                                                &
    tsfc_rg(:,:), albdif_rg(:,:), emis_rad_rg(:,:), cosmu0_rg(:,:),              &
    tot_cld_rg(:,:,:,:), pres_ifc_rg(:,:,:), z_pres_ifc(:,:,:), z_tot_cld(:,:,:,:)

  ! Auxiliary input fields on full grid needed for downscaling corrections
  REAL(wp), INTENT(IN) :: tsfc(:,:), albdif(:,:), rg_trsol_clr_sfc(:,:)

  ! Downscaled output fields (on full grid)
  REAL(wp), INTENT(INOUT) :: aclcov(:,:), lwflxall(:,:,:), trsolall(:,:,:), lwflx_up_sfc(:,:),         &
    trsol_up_toa(:,:), trsol_up_sfc(:,:), trsol_par_sfc(:,:), trsol_dn_sfc_diff(:,:), trsol_clr_sfc(:,:)

  ! Intermediate storage fields needed in the case of MPI parallelization
  REAL(wp), ALLOCATABLE, TARGET ::  z_lwflxall(:,:,:), z_trsolall(:,:,:)

  ! Storage fields needed to downscale transmissitivity differences for solar radiation
  REAL(wp), ALLOCATABLE, TARGET ::             &
  zrg_trdiffsolall(:,:,:), z_trdiffsolall(:,:,:)

  ! Pointers to output fields (no MPI) or intermediate fields (MPI)
  REAL(wp), POINTER ::                                                        &
    p_lwflxall(:,:,:), p_trsolall(:,:,:), p_pres_ifc(:,:,:), p_tot_cld(:,:,:,:)

  ! Additional storage fields to map 2D array(s) to 3D array
  REAL(wp), ALLOCATABLE :: zpg_aux3d(:,:,:), zrg_aux3d(:,:,:), z_aux3d(:,:,:)

  ! Auxiliary fields on full grid for back-interpolated values
  REAL(wp), DIMENSION(nproma,p_patch(jg)%nblks_c) :: tsfc_backintp, alb_backintp

  !> pressure scale for longwave downscaling correction
  REAL(wp), PARAMETER :: pscal = 1._wp/4000._wp

  ! More local variables
  REAL(wp) :: dpresg, pfaclw, intqctot

  REAL(wp), DIMENSION(nproma) :: tqv, dlwem_o_dtg, swfac1, swfac2, lwfac1, lwfac2, logtqv, &
                                 dtrans_o_dalb_clrsfc

  REAL(wp), DIMENSION(nproma,p_patch(jg)%nlevp1) :: intclw, intcli, &
    dtrans_o_dalb_all, dlwflxall_o_dtg, pfacswa

  ! Pointers to types needed to minimize code duplication for MPI/no-MPI cases
  TYPE(t_grid_cells), POINTER     :: p_gcp
  TYPE(t_int_state),  POINTER     :: p_int
  TYPE(t_gridref_state), POINTER  :: p_grf
  TYPE(t_patch),      POINTER     :: p_pp

  ! Indices
  INTEGER :: i_chidx, i_nchdom, nblks_c_lp, i_startrow
  INTEGER :: jb, jk, jk1, jc, i_startblk, i_endblk, i_startidx, i_endidx
  INTEGER :: jc1, jc2, jc3, jc4, jb1, jb2, jb3, jb4

  INTEGER :: nlev, nlevp1      !< number of full and half levels
  INTEGER :: nshift, nlevp1_rg, nlev_tot
  INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk

  INTEGER :: n2dvars, iclcov, itsfc, ialb, iemis, icosmu0, ilwsfc, itrutoa, itrusfc, itrdiff, itrclrsfc, itrparsfc

  LOGICAL :: l_limit(3)
!-----------------------------------------------------------------------

  IF (msg_level >= 10) THEN
    WRITE(message_text,'(a,i2,a,i2)') 'Downscaling of radiation output fields',&
      jgp,' =>',jg
    CALL message('downscale_rad_output',message_text)
  ENDIF

  ! For the time being, the radiation grid is assumed to have the same levels
  ! as the full grid
  nlev   = p_patch(jg)%nlev
  nlevp1 = p_patch(jg)%nlevp1

  ! nlev_rg carries the number of model levels of the reduced grid,
  ! which may be larger than nlev
  nlevp1_rg = nlev_rg + 1
  nshift = nlev_rg - nlev ! resulting shift parameter
  p_grf => p_grf_state_local_parent(jg)
  p_int => p_int_state_local_parent(jg)
  p_gcp => p_patch_local_parent(jg)%cells
  p_pp  => p_patch_local_parent(jg)

  i_nchdom = MAX(1,p_patch(jg)%n_childdom)
  i_chidx  = p_patch(jg)%parent_child_index

  nblks_c_lp = p_pp%nblks_c

  ! pointers to child index/block
  iidx => p_gcp%child_idx
  iblk => p_gcp%child_blk

  ! named constants for accessing 2D variables contained in zrg_aux3d
  n2dvars = nshift+11
  iclcov  = nshift+1
  itsfc   = nshift+2
  ialb    = nshift+3
  iemis   = nshift+4
  icosmu0 = nshift+5
  ilwsfc  = nshift+6
  itrutoa = nshift+7
  itrusfc = nshift+8
  itrdiff = nshift+9
  itrclrsfc = nshift+10
  itrparsfc = nshift+11


  ! Allocation of local storage fields at local parent level in MPI-case
  IF (jgp == 0 .AND. .NOT. l_limited_area) THEN

    ALLOCATE( z_lwflxall(nproma,nlevp1_rg,nblks_c_lp),        z_trsolall(nproma,nlevp1_rg,nblks_c_lp),          &
              zpg_aux3d(nproma,n2dvars,p_patch(jgp)%nblks_c),                                                   &
              zrg_aux3d(nproma,n2dvars,nblks_c_lp),           z_aux3d(nproma,11,p_patch(jg)%nblks_c),           &
              zrg_trdiffsolall(nproma,nlevp1_rg,nblks_c_lp),  z_trdiffsolall(nproma,nlevp1,p_patch(jg)%nblks_c) )


    ! Perform communication from parent to local parent grid in the MPI case,
    ! and set pointers such that further processing is the same for MPI / non-MPI cases

!$OMP PARALLEL
    CALL init(zpg_aux3d(:,1:nshift,:))
    CALL copy(rg_aclcov(:,:), zpg_aux3d(:,iclcov,:))
    CALL copy(tsfc_rg(:,:), zpg_aux3d(:,itsfc,:))
    CALL copy(albdif_rg(:,:), zpg_aux3d(:,ialb,:))
    CALL copy(emis_rad_rg(:,:), zpg_aux3d(:,iemis,:))
    CALL copy(cosmu0_rg(:,:), zpg_aux3d(:,icosmu0,:))
    CALL copy(rg_lwflx_up_sfc(:,:), zpg_aux3d(:,ilwsfc,:))
    CALL copy(rg_trsol_up_toa(:,:), zpg_aux3d(:,itrutoa,:))
    CALL copy(rg_trsol_up_sfc(:,:), zpg_aux3d(:,itrusfc,:))
    CALL copy(rg_trsol_dn_sfc_diff(:,:), zpg_aux3d(:,itrdiff,:))
    CALL copy(rg_trsol_clr_sfc(:,:), zpg_aux3d(:,itrclrsfc,:))
    CALL copy(rg_trsol_par_sfc(:,:), zpg_aux3d(:,itrparsfc,:))
!$OMP END PARALLEL

    nlev_tot = 2*nlevp1_rg + n2dvars

    CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, 3, nlev_tot, &
                            RECV1=z_lwflxall, SEND1=rg_lwflxall,     &
                            RECV2=z_trsolall, SEND2=rg_trsolall,     &
                            RECV3=zrg_aux3d , SEND3=zpg_aux3d        )

    p_lwflxall   => z_lwflxall
    p_trsolall   => z_trsolall
    p_pres_ifc   => z_pres_ifc
    p_tot_cld    => z_tot_cld
  ELSE
    ALLOCATE(zrg_aux3d(nproma,n2dvars,nblks_c_lp),          z_aux3d(nproma,11,p_patch(jg)%nblks_c),          &
             zrg_trdiffsolall(nproma,nlevp1_rg,nblks_c_lp), z_trdiffsolall(nproma,nlevp1,p_patch(jg)%nblks_c))

    p_lwflxall   => rg_lwflxall
    p_trsolall   => rg_trsolall
    p_pres_ifc   => pres_ifc_rg
    p_tot_cld    => tot_cld_rg

!$OMP PARALLEL
    CALL init(zrg_aux3d(:,1:nshift,:))
    CALL copy(rg_aclcov(:,:), zrg_aux3d(:,iclcov,:))
    CALL copy(tsfc_rg(:,:), zrg_aux3d(:,itsfc,:))
    CALL copy(albdif_rg(:,:), zrg_aux3d(:,ialb,:))
    CALL copy(emis_rad_rg(:,:), zrg_aux3d(:,iemis,:))
    CALL copy(cosmu0_rg(:,:), zrg_aux3d(:,icosmu0,:))
    CALL copy(rg_lwflx_up_sfc(:,:), zrg_aux3d(:,ilwsfc,:))
    CALL copy(rg_trsol_up_toa(:,:), zrg_aux3d(:,itrutoa,:))
    CALL copy(rg_trsol_up_sfc(:,:), zrg_aux3d(:,itrusfc,:))
    CALL copy(rg_trsol_dn_sfc_diff(:,:), zrg_aux3d(:,itrdiff,:))
    CALL copy(rg_trsol_clr_sfc(:,:), zrg_aux3d(:,itrclrsfc,:))
    CALL copy(rg_trsol_par_sfc(:,:), zrg_aux3d(:,itrparsfc,:))
!$OMP END PARALLEL
  ENDIF

  ! Compute transmissivity differences before downscaling
  ! Note: transmissivities must monotonically decrease from top to bottom in order to obtain
  ! non-negative solar heating rates. The limiter in interpol_scal_nudging can ensure this
  ! only when passing transmissivity differences, rather than raw transmissivities, into the
  ! routine

  ! Start/End block in the local parent domain - same extent to nest boundary region as used for radiation calculation
  IF (jg == 1 .AND. l_limited_area) THEN
    i_startrow = grf_fbk_start_c
  ELSE
    i_startrow = grf_ovlparea_start_c
  ENDIF

!$OMP PARALLEL PRIVATE(i_startblk, i_endblk)

  i_startblk = p_gcp%start_blk(i_startrow,i_chidx)
  i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,                           &
                       i_startidx, i_endidx, i_startrow, min_rlcell_int)

    DO jk = 1, nshift+1
      DO jc = i_startidx, i_endidx
        zrg_trdiffsolall(jc,jk,jb) = p_trsolall(jc,jk,jb)
      ENDDO
    ENDDO

    DO jk = nshift+2, nlevp1_rg
      DO jc = i_startidx, i_endidx
        zrg_trdiffsolall(jc,jk,jb) = p_trsolall(jc,jk-1,jb) - p_trsolall(jc,jk,jb)
      ENDDO
    ENDDO

  ENDDO
!$OMP END DO

  ! For the limited-area mode, possibly undefined grid points entering into the subsequent interpolation
  ! need to be set to zero
  IF (jg == 1 .AND. l_limited_area) THEN
    i_startblk = p_gcp%start_block(grf_ovlparea_start_c)
    i_endblk   = p_gcp%end_block(grf_ovlparea_start_c)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
                         grf_ovlparea_start_c, grf_ovlparea_start_c)

      DO jc = i_startidx, i_endidx
        zrg_trdiffsolall(jc,:,jb) = 0._wp
        p_lwflxall(jc,:,jb)       = 0._wp
        zrg_aux3d(jc,:,jb)        = 0._wp
      ENDDO

    ENDDO
!$OMP END DO
  ENDIF

!$OMP END PARALLEL

  ! Interpolate reduced-grid fields to full grid

  ! interpol_vec/scal_nudging needs all fields with full boundaries, so the arrays
  ! have to be sync'd before calling these routines.
  ! Please note that we cannot use sync_patch_array here (comparing parallel/non parallel results)
  ! since the arrays don't start with lower bound 1 in the non parallel case!

  nlev_tot = 2*nlevp1_rg + n2dvars

  IF (.NOT. my_process_is_mpi_seq()) THEN
    CALL exchange_data_mult(p_pp%comm_pat_c, 3, nlev_tot, recv1=p_lwflxall, &
      &                     recv2=zrg_trdiffsolall, recv3=zrg_aux3d   )
  END IF

  IF (p_test_run) THEN
    z_trdiffsolall = 0._wp
    lwflxall       = 0._wp
    z_aux3d        = 0._wp
  ENDIF


  l_limit(1) = .TRUE.      ! limit transmissivity differences to non-negative values
  l_limit(2:3) = .FALSE.

  CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, nshift, 3, 1, &
    &                         f3din1=zrg_trdiffsolall, f3dout1=z_trdiffsolall,          &
    &                         f3din2=p_lwflxall,       f3dout2=lwflxall,                &
    &                         f3din3=zrg_aux3d,        f3dout3=z_aux3d,                 &
    &                         llimit_nneg=l_limit,     overshoot_fac=1.0_wp)

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

  ! For the limited-area mode, near-boundary grid points filled with partly invalid data need to be reset
  IF (jg == 1 .AND. l_limited_area) THEN
    i_startblk = p_gcp%start_block(grf_fbk_start_c)
    i_endblk   = p_gcp%end_block(grf_fbk_start_c)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
                         grf_fbk_start_c, grf_fbk_start_c)


      DO jk = nshift+1, n2dvars
        DO jc = i_startidx, i_endidx
          z_aux3d(iidx(jc,jb,1),jk-nshift,iblk(jc,jb,1)) = zrg_aux3d(jc,jk,jb)
          z_aux3d(iidx(jc,jb,2),jk-nshift,iblk(jc,jb,2)) = zrg_aux3d(jc,jk,jb)
          z_aux3d(iidx(jc,jb,3),jk-nshift,iblk(jc,jb,3)) = zrg_aux3d(jc,jk,jb)
          z_aux3d(iidx(jc,jb,4),jk-nshift,iblk(jc,jb,4)) = zrg_aux3d(jc,jk,jb)
        ENDDO
      ENDDO

      DO jk = nshift+1, nlevp1_rg
        DO jc = i_startidx, i_endidx
          z_trdiffsolall(iidx(jc,jb,1),jk-nshift,iblk(jc,jb,1)) = zrg_trdiffsolall(jc,jk,jb)
          z_trdiffsolall(iidx(jc,jb,2),jk-nshift,iblk(jc,jb,2)) = zrg_trdiffsolall(jc,jk,jb)
          z_trdiffsolall(iidx(jc,jb,3),jk-nshift,iblk(jc,jb,3)) = zrg_trdiffsolall(jc,jk,jb)
          z_trdiffsolall(iidx(jc,jb,4),jk-nshift,iblk(jc,jb,4)) = zrg_trdiffsolall(jc,jk,jb)

          lwflxall(iidx(jc,jb,1),jk-nshift,iblk(jc,jb,1)) = p_lwflxall(jc,jk,jb)
          lwflxall(iidx(jc,jb,2),jk-nshift,iblk(jc,jb,2)) = p_lwflxall(jc,jk,jb)
          lwflxall(iidx(jc,jb,3),jk-nshift,iblk(jc,jb,3)) = p_lwflxall(jc,jk,jb)
          lwflxall(iidx(jc,jb,4),jk-nshift,iblk(jc,jb,4)) = p_lwflxall(jc,jk,jb)
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO
  ENDIF

  CALL copy(z_aux3d(:,1,:), aclcov(:,:))
  CALL copy(z_aux3d(:,2,:), tsfc_backintp(:,:))
  CALL copy(z_aux3d(:,3,:), alb_backintp(:,:))
!$OMP BARRIER

  ! Reconstruct solar transmissivities from interpolated transmissivity differences

  ! Start/End block in the child domain
  i_startblk = p_patch(jg)%cells%start_blk(grf_bdywidth_c+1,1)
  i_endblk   = p_patch(jg)%cells%end_blk(min_rlcell_int,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell_int)

    DO jc = i_startidx, i_endidx
      trsolall(jc,1,jb) = z_trdiffsolall(jc,1,jb)
    ENDDO

    DO jk = 2, nlevp1
      DO jc = i_startidx, i_endidx
        trsolall(jc,jk,jb) = trsolall(jc,jk-1,jb) - z_trdiffsolall(jc,jk,jb)
      ENDDO
    ENDDO
  ENDDO
!$OMP END DO


  ! Apply empirical downscaling corrections depending on variations
  ! in ground temperature and albedo

  ! Start/End block in the local parent domain
  i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
  i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jk1,tqv,intclw,intcli,dpresg,pfacswa, &
!$OMP            dlwem_o_dtg,swfac1,swfac2,lwfac1,lwfac2,dtrans_o_dalb_all,         &
!$OMP            pfaclw,intqctot,dlwflxall_o_dtg,jc1,jc2,jc3,jc4,jb1,jb2,jb3,       &
!$OMP            jb4,logtqv,dtrans_o_dalb_clrsfc) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,                      &
                       i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int)

    tqv(:)            = 0._wp
    intclw(:,nlevp1)  = 0._wp
    intcli(:,nlevp1)  = 0._wp
    pfacswa(:,nlevp1) = 1._wp
    DO jk = nlev,1,-1
      jk1 = jk + nshift
      DO jc = i_startidx, i_endidx
        dpresg        = (p_pres_ifc(jc,jk1+1,jb) - p_pres_ifc(jc,jk1,jb))/grav
        tqv(jc)       = tqv(jc)+p_tot_cld(jc,jk1,jb,iqv)*dpresg
        intclw(jc,jk) = intclw(jc,jk+1)+p_tot_cld(jc,jk1,jb,iqc)*dpresg
        intcli(jc,jk) = intcli(jc,jk+1)+p_tot_cld(jc,jk1,jb,iqi)*dpresg
        pfacswa(jc,jk)= 1._wp-0.16_wp*(p_pres_ifc(jc,nlevp1_rg,jb)-p_pres_ifc(jc,jk1,jb))/&
                        p_pres_ifc(jc,nlevp1_rg,jb)
      ENDDO
    ENDDO

    DO jc = i_startidx, i_endidx
      logtqv(jc) = LOG(MAX(1._wp,tqv(jc)))
      dlwem_o_dtg(jc) = zrg_aux3d(jc,iemis,jb)*4._wp*stbo*zrg_aux3d(jc,itsfc,jb)**3
      swfac1(jc) = EXP(0.36_wp*LOG( MAX(1.e-3_wp,p_trsolall(jc,nlevp1_rg,jb))/    &
                   MAX(1.e-3_wp,zrg_aux3d(jc,itrclrsfc,jb)) ))
      swfac2(jc) = EXP(0.1_wp*LOG( MAX(0.25_wp,3._wp*zrg_aux3d(jc,icosmu0,jb)) ))
      lwfac2(jc) = 0.92_wp*EXP(-0.07_wp*logtqv(jc))
    ENDDO
    DO jc = i_startidx, i_endidx
      lwfac1(jc) = MERGE(1.677_wp, 0.4388_wp, tqv(jc) > 15._wp) &
        * EXP(MERGE(-0.72_wp, -0.225_wp, tqv(jc) > 15._wp)*logtqv(jc))
    ENDDO

    DO jk = 1,nlevp1
      jk1 = jk + nshift
      DO jc = i_startidx, i_endidx
        dtrans_o_dalb_all(jc,jk) = - p_trsolall(jc,nlevp1_rg,jb)*pfacswa(jc,jk)*swfac1(jc)/ &
                                   ( (1._wp-zrg_aux3d(jc,ialb,jb))*swfac2(jc) )

        pfaclw = lwfac1(jc)+(lwfac2(jc)-lwfac1(jc))*EXP(-SQRT((p_pres_ifc(jc,nlevp1_rg,jb)- &
                 p_pres_ifc(jc,jk1,jb))*pscal))
        intqctot = MIN(0.30119_wp,MAX(1.008e-3_wp,intclw(jc,jk)+0.2_wp*intcli(jc,jk)))

        dlwflxall_o_dtg(jc,jk) = -dlwem_o_dtg(jc)*pfaclw*(1._wp-(6.9_wp+LOG(intqctot))/5.7_wp)
      ENDDO
    ENDDO

    DO jc = i_startidx, i_endidx
      dtrans_o_dalb_clrsfc(jc) = - zrg_aux3d(jc,itrclrsfc,jb)/((1._wp-zrg_aux3d(jc,ialb,jb))*swfac2(jc))
    ENDDO

    ! Now apply the corrections
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      jc1 = iidx(jc,jb,1)
      jc2 = iidx(jc,jb,2)
      jc3 = iidx(jc,jb,3)
      jc4 = iidx(jc,jb,4)
      jb1 = iblk(jc,jb,1)
      jb2 = iblk(jc,jb,2)
      jb3 = iblk(jc,jb,3)
      jb4 = iblk(jc,jb,4)
      DO jk = 1,nlevp1
#else
!CDIR NOLOOPCHG
    DO jk = 1,nlevp1
!CDIR NODEP,VOVERTAKE,VOB
      DO jc = i_startidx, i_endidx
        jc1 = iidx(jc,jb,1)
        jc2 = iidx(jc,jb,2)
        jc3 = iidx(jc,jb,3)
        jc4 = iidx(jc,jb,4)
        jb1 = iblk(jc,jb,1)
        jb2 = iblk(jc,jb,2)
        jb3 = iblk(jc,jb,3)
        jb4 = iblk(jc,jb,4)
#endif

        trsolall(jc1,jk,jb1) = MAX(trsolall(jc1,jk,jb1) + dtrans_o_dalb_all(jc,jk)* &
          ( albdif(jc1,jb1) - alb_backintp(jc1,jb1) ), 0._wp)
        trsolall(jc2,jk,jb2) = MAX(trsolall(jc2,jk,jb2) + dtrans_o_dalb_all(jc,jk)* &
          ( albdif(jc2,jb2) - alb_backintp(jc2,jb2) ), 0._wp)
        trsolall(jc3,jk,jb3) = MAX(trsolall(jc3,jk,jb3) + dtrans_o_dalb_all(jc,jk)* &
          ( albdif(jc3,jb3) - alb_backintp(jc3,jb3) ), 0._wp)
        trsolall(jc4,jk,jb4) = MAX(trsolall(jc4,jk,jb4) + dtrans_o_dalb_all(jc,jk)* &
          ( albdif(jc4,jb4) - alb_backintp(jc4,jb4) ), 0._wp)

        lwflxall(jc1,jk,jb1) = lwflxall(jc1,jk,jb1) + dlwflxall_o_dtg(jc,jk)* &
          ( tsfc(jc1,jb1) - tsfc_backintp(jc1,jb1) )
        lwflxall(jc2,jk,jb2) = lwflxall(jc2,jk,jb2) + dlwflxall_o_dtg(jc,jk)* &
          ( tsfc(jc2,jb2) - tsfc_backintp(jc2,jb2) )
        lwflxall(jc3,jk,jb3) = lwflxall(jc3,jk,jb3) + dlwflxall_o_dtg(jc,jk)* &
          ( tsfc(jc3,jb3) - tsfc_backintp(jc3,jb3) )
        lwflxall(jc4,jk,jb4) = lwflxall(jc4,jk,jb4) + dlwflxall_o_dtg(jc,jk)* &
          ( tsfc(jc4,jb4) - tsfc_backintp(jc4,jb4) )

      ENDDO
    ENDDO

!DIR$ IVDEP
    DO jc = i_startidx, i_endidx
      jc1 = iidx(jc,jb,1)
      jc2 = iidx(jc,jb,2)
      jc3 = iidx(jc,jb,3)
      jc4 = iidx(jc,jb,4)
      jb1 = iblk(jc,jb,1)
      jb2 = iblk(jc,jb,2)
      jb3 = iblk(jc,jb,3)
      jb4 = iblk(jc,jb,4)

      lwflx_up_sfc(jc1,jb1) = z_aux3d(jc1,6,jb1) + dlwem_o_dtg(jc)* (tsfc(jc1,jb1) - tsfc_backintp(jc1,jb1))
      lwflx_up_sfc(jc2,jb2) = z_aux3d(jc2,6,jb2) + dlwem_o_dtg(jc)* (tsfc(jc2,jb2) - tsfc_backintp(jc2,jb2))
      lwflx_up_sfc(jc3,jb3) = z_aux3d(jc3,6,jb3) + dlwem_o_dtg(jc)* (tsfc(jc3,jb3) - tsfc_backintp(jc3,jb3))
      lwflx_up_sfc(jc4,jb4) = z_aux3d(jc4,6,jb4) + dlwem_o_dtg(jc)* (tsfc(jc4,jb4) - tsfc_backintp(jc4,jb4))

      ! The downscaling corrections for the upward fluxes assume that the downward radiation is unaffected by
      ! the surface albedo, so that the correction on upward radiation is the same as on net radiation
      ! except for the sign because the convention is that upward fluxes are upward-positive
      trsol_up_toa(jc1,jb1) = MAX(z_aux3d(jc1,7,jb1) - dtrans_o_dalb_all(jc,1)* &
          ( albdif(jc1,jb1) - alb_backintp(jc1,jb1) ), 0._wp)
      trsol_up_toa(jc2,jb2) = MAX(z_aux3d(jc2,7,jb2) - dtrans_o_dalb_all(jc,1)* &
          ( albdif(jc2,jb2) - alb_backintp(jc2,jb2) ), 0._wp)
      trsol_up_toa(jc3,jb3) = MAX(z_aux3d(jc3,7,jb3) - dtrans_o_dalb_all(jc,1)* &
          ( albdif(jc3,jb3) - alb_backintp(jc3,jb3) ), 0._wp)
      trsol_up_toa(jc4,jb4) = MAX(z_aux3d(jc4,7,jb4) - dtrans_o_dalb_all(jc,1)* &
          ( albdif(jc4,jb4) - alb_backintp(jc4,jb4) ), 0._wp)

      trsol_up_sfc(jc1,jb1) = MAX(z_aux3d(jc1,8,jb1) - dtrans_o_dalb_all(jc,nlevp1)* &
          ( albdif(jc1,jb1) - alb_backintp(jc1,jb1) ), 0._wp)
      trsol_up_sfc(jc2,jb2) = MAX(z_aux3d(jc2,8,jb2) - dtrans_o_dalb_all(jc,nlevp1)* &
          ( albdif(jc2,jb2) - alb_backintp(jc2,jb2) ), 0._wp)
      trsol_up_sfc(jc3,jb3) = MAX(z_aux3d(jc3,8,jb3) - dtrans_o_dalb_all(jc,nlevp1)* &
          ( albdif(jc3,jb3) - alb_backintp(jc3,jb3) ), 0._wp)
      trsol_up_sfc(jc4,jb4) = MAX(z_aux3d(jc4,8,jb4) - dtrans_o_dalb_all(jc,nlevp1)* &
          ( albdif(jc4,jb4) - alb_backintp(jc4,jb4) ), 0._wp)

      ! Downscaling of clear-sky net transmissivities at the surface
      ! They are used in radheat for the tile corrections
      trsol_clr_sfc(jc1,jb1) = MAX(z_aux3d(jc1,10,jb1) + dtrans_o_dalb_clrsfc(jc)* &
          ( albdif(jc1,jb1) - alb_backintp(jc1,jb1) ), 0._wp)
      trsol_clr_sfc(jc2,jb2) = MAX(z_aux3d(jc2,10,jb2) + dtrans_o_dalb_clrsfc(jc)* &
          ( albdif(jc2,jb2) - alb_backintp(jc2,jb2) ), 0._wp)
      trsol_clr_sfc(jc3,jb3) = MAX(z_aux3d(jc3,10,jb3) + dtrans_o_dalb_clrsfc(jc)* &
          ( albdif(jc3,jb3) - alb_backintp(jc3,jb3) ), 0._wp)
      trsol_clr_sfc(jc4,jb4) = MAX(z_aux3d(jc4,10,jb4) + dtrans_o_dalb_clrsfc(jc)* &
          ( albdif(jc4,jb4) - alb_backintp(jc4,jb4) ), 0._wp)

      ! Note: the downward diffuse radiation must not undergo a correction based on the surface albedo!
      ! Only upward and net radiation are affected by the surface properties
      trsol_dn_sfc_diff(jc1,jb1) = MAX(z_aux3d(jc1,9,jb1), 0._wp)
      trsol_dn_sfc_diff(jc2,jb2) = MAX(z_aux3d(jc2,9,jb2), 0._wp)
      trsol_dn_sfc_diff(jc3,jb3) = MAX(z_aux3d(jc3,9,jb3), 0._wp)
      trsol_dn_sfc_diff(jc4,jb4) = MAX(z_aux3d(jc4,9,jb4), 0._wp)

      ! Ensure that the diffuse downward radiation does not exceed the total radiation
      trsol_dn_sfc_diff(jc1,jb1) = MIN(trsol_dn_sfc_diff(jc1,jb1),trsolall(jc1,nlevp1,jb1)/(1._wp-albdif(jc1,jb1)))
      trsol_dn_sfc_diff(jc2,jb2) = MIN(trsol_dn_sfc_diff(jc2,jb2),trsolall(jc2,nlevp1,jb2)/(1._wp-albdif(jc2,jb2)))
      trsol_dn_sfc_diff(jc3,jb3) = MIN(trsol_dn_sfc_diff(jc3,jb3),trsolall(jc3,nlevp1,jb3)/(1._wp-albdif(jc3,jb3)))
      trsol_dn_sfc_diff(jc4,jb4) = MIN(trsol_dn_sfc_diff(jc4,jb4),trsolall(jc4,nlevp1,jb4)/(1._wp-albdif(jc4,jb4)))

      ! No albedo correction for transmissivity of photosynthetically active radiation
      trsol_par_sfc(jc1,jb1) = MAX(z_aux3d(jc1,11,jb1), 0._wp)
      trsol_par_sfc(jc2,jb2) = MAX(z_aux3d(jc2,11,jb2), 0._wp)
      trsol_par_sfc(jc3,jb3) = MAX(z_aux3d(jc3,11,jb3), 0._wp)
      trsol_par_sfc(jc4,jb4) = MAX(z_aux3d(jc4,11,jb4), 0._wp)
    ENDDO

  ENDDO
!$OMP END DO

  ! Finally ensure again that transmissivities increase monotonically with height,
  ! which may be broken by the downscaling correction applied above

  ! Start/End block in the child domain
  i_startblk = p_patch(jg)%cells%start_blk(grf_bdywidth_c+1,1)
  i_endblk   = p_patch(jg)%cells%end_blk(min_rlcell_int,i_nchdom)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell_int)

    DO jk = nlev, 1, -1
      DO jc = i_startidx, i_endidx
        trsolall(jc,jk,jb) = MIN(1._wp,MAX(trsolall(jc,jk,jb),trsolall(jc,jk+1,jb)))
      ENDDO
    ENDDO

  ENDDO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

  IF (jgp == 0 .AND. .NOT. l_limited_area) THEN
    DEALLOCATE(z_lwflxall, z_trsolall, zpg_aux3d)
  ENDIF

  DEALLOCATE(zrg_aux3d, z_aux3d, zrg_trdiffsolall, z_trdiffsolall)

END SUBROUTINE downscale_rad_output

!>
!! This routine averages the input fields for Ritter-Geleyn radiation to the next coarser grid
!! level. It is an adaptation of upscale_rad_input to the Ritter-Geleyn radiation.
!!
!! @par Revision History
!! Thorsten Reinhardt, AGeoBw, 2010-12-06
!!
SUBROUTINE upscale_rad_input_rg(jg, jgp, nlev_rg, nlevp1_rg,         &
  &  cosmu0, albdif, alb_ther, temp_ifc, dpres_mc,                   &
  &  tot_cld, clc, sqv, duco2, duo3,                                 &
  &  aeq1, aeq2, aeq3, aeq4, aeq5, pres_sfc,pres_ifc,                &
  &  rg_cosmu0, rg_albdif, rg_alb_ther, rg_temp_ifc, rg_dpres_mc,    &
  &  rg_tot_cld, rg_clc, rg_sqv, rg_duco2, rg_duo3,                  &
  &  rg_aeq1, rg_aeq2, rg_aeq3, rg_aeq4, rg_aeq5, rg_pres_sfc )


  ! Input grid parameters
  INTEGER, INTENT(IN)  :: jg, jgp  ! domain IDs of main and reduced grids
  INTEGER, INTENT(IN)  :: nlev_rg, nlevp1_rg  ! number of model levels on reduced grid

  ! Input fields (on full grid)
  REAL(wp), INTENT(IN) ::                                                             &
    & cosmu0(:,:), albdif(:,:), alb_ther(:,:), temp_ifc(:,:,:),                       &
    & dpres_mc(:,:,:), tot_cld(:,:,:,:), clc(:,:,:), sqv(:,:,:), duco2(:,:,:),        &
    & duo3(:,:,:), aeq1(:,:,:),aeq2(:,:,:),aeq3(:,:,:),aeq4(:,:,:),aeq5(:,:,:),       &
    & pres_sfc(:,:)

  ! Input field (on full grid) without corresponding output fields on reduced grid
  REAL(wp), INTENT(IN) ::  &
    & pres_ifc(:,:,:)

  ! Corresponding output fields (on reduced grid)
  REAL(wp), TARGET, INTENT(OUT) ::                                                            &
    & rg_cosmu0(:,:), rg_albdif(:,:), rg_alb_ther(:,:), rg_temp_ifc(:,:,:),                   &
    & rg_dpres_mc(:,:,:), rg_tot_cld(:,:,:,:), rg_clc(:,:,:), rg_sqv(:,:,:), rg_duco2(:,:,:), &
    & rg_duo3(:,:,:), rg_aeq1(:,:,:), rg_aeq2(:,:,:), rg_aeq3(:,:,:), rg_aeq4(:,:,:),         &
    & rg_aeq5(:,:,:), rg_pres_sfc(:,:)

  ! Intermediate storage fields needed in the case of MPI parallelization
  REAL(wp), ALLOCATABLE, TARGET ::                                                        &
    & z_cosmu0(:,:), z_albdif(:,:), z_alb_ther(:,:), z_temp_ifc(:,:,:),                   &
    & z_dpres_mc(:,:,:), z_tot_cld(:,:,:,:), z_clc(:,:,:), z_sqv(:,:,:), z_duco2(:,:,:),  &
    & z_duo3(:,:,:), z_aeq1(:,:,:), z_aeq2(:,:,:), z_aeq3(:,:,:), z_aeq4(:,:,:),          &
    & z_aeq5(:,:,:), z_pres_sfc(:,:), z_aux3d(:,:,:),  zrg_aux3d(:,:,:)


  ! Pointers to output fields (no MPI) or intermediate fields (MPI)
  REAL(wp), POINTER ::                                                                    &
    & p_cosmu0(:,:), p_albdif(:,:), p_alb_ther(:,:), p_temp_ifc(:,:,:),                   &
    & p_dpres_mc(:,:,:), p_tot_cld(:,:,:,:), p_clc(:,:,:), p_sqv(:,:,:), p_duco2(:,:,:),  &
    & p_duo3(:,:,:), p_aeq1(:,:,:), p_aeq2(:,:,:), p_aeq3(:,:,:), p_aeq4(:,:,:),          &
    & p_aeq5(:,:,:), p_pres_sfc(:,:)


  ! Pointers to types needed to minimize code duplication for MPI/no-MPI cases
  TYPE(t_grid_cells), POINTER     :: p_gcp
  TYPE(t_gridref_state), POINTER  :: p_grf
  TYPE(t_patch),      POINTER     :: p_pp

  ! Indices
  INTEGER :: jb, jc, jk, jk1, i_chidx, i_nchdom, &
             i_startblk, i_endblk, i_startidx, i_endidx, nblks_c_lp

  INTEGER :: nlev, nlevp1      !< number of full and half levels
  INTEGER :: nshift
  REAL(wp) :: z_help_pres_ratio, z_rho_1

  INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
  REAL(wp), POINTER :: p_fbkwgt(:,:,:)


  !-----------------------------------------------------------------------

  IF (msg_level >= 10) THEN
    WRITE(message_text,'(a,i2,a,i2)') 'Upscaling of radiation input fields',&
      jg,' =>',jgp
    CALL message('upscale_rad_input',message_text)
  ENDIF

  ! For the time being, the radiation grid is assumed to have the same levels
  ! as the full grid
  ! Number of levels of the full grid
  nlev   = p_patch(jg)%nlev
  nlevp1 = p_patch(jg)%nlevp1

  ! nlev_rg carries the number of model levels of the reduced grid,
  ! which may be larger than nlev
  nshift = nlev_rg - nlev ! resulting shift parameter

  p_grf => p_grf_state_local_parent(jg)
  p_gcp => p_patch_local_parent(jg)%cells
  p_pp  => p_patch_local_parent(jg)

  i_chidx  = p_patch(jg)%parent_child_index
  i_nchdom = MAX(1,p_patch(jgp)%n_childdom)

  ! Set pointers to index and coefficient fields for cell-based variables
  iidx => p_gcp%child_idx
  iblk => p_gcp%child_blk

  p_fbkwgt => p_grf%fbk_wgt_bln

  ! Allocation of local storage fields at local parent level in MPI-case
  IF (jgp == 0) THEN
    nblks_c_lp = p_gcp%end_blk(min_rlcell,i_chidx)

    ALLOCATE(z_cosmu0(nproma,nblks_c_lp), z_albdif(nproma,nblks_c_lp),                     &
             z_alb_ther(nproma,nblks_c_lp), z_temp_ifc(nproma,nlevp1_rg,nblks_c_lp),       &
             z_dpres_mc(nproma,nlev_rg,nblks_c_lp),z_tot_cld(nproma,nlev_rg,nblks_c_lp,3), &
             z_clc(nproma,nlev_rg,nblks_c_lp),                                             &
             z_sqv(nproma,nlev_rg,nblks_c_lp),z_duco2(nproma,nlev_rg,nblks_c_lp),          &
             z_duo3(nproma,nlev_rg,nblks_c_lp), z_aeq1(nproma,nlev_rg,nblks_c_lp),         &
             z_aeq2(nproma,nlev_rg,nblks_c_lp), z_aeq3(nproma,nlev_rg,nblks_c_lp),         &
             z_aeq4(nproma,nlev_rg,nblks_c_lp), z_aeq5(nproma,nlev_rg,nblks_c_lp),         &
             z_pres_sfc(nproma,nblks_c_lp), z_aux3d(nproma,6,nblks_c_lp),                  &
             zrg_aux3d(nproma,4,p_patch(jgp)%nblks_c)                                      )

    ! Set pointers to either the parent-level variables (non-MPI case) or to the
    ! intermediate storage fields (MPI case)
    p_cosmu0     => z_cosmu0
    p_albdif     => z_albdif
    p_alb_ther   => z_alb_ther
    p_temp_ifc   => z_temp_ifc
    p_dpres_mc   => z_dpres_mc
    p_tot_cld    => z_tot_cld
    p_clc        => z_clc
    p_sqv        => z_sqv
    p_duco2      => z_duco2
    p_duo3       => z_duo3
    p_aeq1       => z_aeq1
    p_aeq2       => z_aeq2
    p_aeq3       => z_aeq3
    p_aeq4       => z_aeq4
    p_aeq5       => z_aeq5
    p_pres_sfc   => z_pres_sfc
  ELSE
    p_cosmu0     => rg_cosmu0
    p_albdif     => rg_albdif
    p_alb_ther   => rg_alb_ther
    p_temp_ifc   => rg_temp_ifc
    p_dpres_mc   => rg_dpres_mc
    p_tot_cld    => rg_tot_cld
    p_clc        => rg_clc
    p_sqv        => rg_sqv
    p_duco2      => rg_duco2
    p_duo3       => rg_duo3
    p_aeq1       => rg_aeq1
    p_aeq2       => rg_aeq2
    p_aeq3       => rg_aeq3
    p_aeq4       => rg_aeq4
    p_aeq5       => rg_aeq5
    p_pres_sfc   => rg_pres_sfc
  ENDIF

  IF (p_test_run) THEN
    p_cosmu0     = 0._wp
    p_albdif     = 0._wp
    p_alb_ther   = 0._wp
    p_temp_ifc   = 0._wp
    p_dpres_mc   = 0._wp
    p_tot_cld    = 0._wp
    p_clc        = 0._wp
    p_sqv        = 0._wp
    p_duco2      = 0._wp
    p_duo3       = 0._wp
    p_aeq1       = 0._wp
    p_aeq2       = 0._wp
    p_aeq3       = 0._wp
    p_aeq4       = 0._wp
    p_aeq5       = 0._wp
    p_pres_sfc   = 0._wp
  ENDIF

  ! Now average input fields to parent grid cells

  ! Start/End block in the parent domain
  i_startblk = p_gcp%start_blk(grf_ovlparea_start_c,i_chidx)
  i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jk1,z_help_pres_ratio,&
!$OMP z_rho_1) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,                           &
                       i_startidx, i_endidx, grf_ovlparea_start_c, min_rlcell_int)

    DO jc = i_startidx, i_endidx

      p_cosmu0(jc,jb) =                                         &
        cosmu0(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        cosmu0(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        cosmu0(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        cosmu0(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_albdif(jc,jb) =                                            &
        albdif(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        albdif(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        albdif(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        albdif(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_alb_ther(jc,jb) =                                         &
        alb_ther(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        alb_ther(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        alb_ther(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        alb_ther(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_pres_sfc(jc,jb) =                                         &
        pres_sfc(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        pres_sfc(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        pres_sfc(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        pres_sfc(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_temp_ifc(jc,nlevp1_rg,jb) = SQRT(SQRT(                              &
        temp_ifc(iidx(jc,jb,1),nlevp1,iblk(jc,jb,1))**4*p_fbkwgt(jc,jb,1) + &
        temp_ifc(iidx(jc,jb,2),nlevp1,iblk(jc,jb,2))**4*p_fbkwgt(jc,jb,2) + &
        temp_ifc(iidx(jc,jb,3),nlevp1,iblk(jc,jb,3))**4*p_fbkwgt(jc,jb,3) + &
        temp_ifc(iidx(jc,jb,4),nlevp1,iblk(jc,jb,4))**4*p_fbkwgt(jc,jb,4) ) )

    ENDDO

    ! combine 2D fields in a 3D field to speed up MPI communication
    IF (jgp == 0) THEN
      DO jc = i_startidx, i_endidx
        z_aux3d(jc,1,jb) = p_cosmu0(jc,jb)
        z_aux3d(jc,2,jb) = p_albdif(jc,jb)
        z_aux3d(jc,3,jb) = p_alb_ther(jc,jb)
        z_aux3d(jc,4,jb) = p_pres_sfc(jc,jb)
      ENDDO
    ENDIF

#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      DO jk = 1, nlev
        jk1 = jk + nshift
#else
    DO jk = 1, nlev
      jk1 = jk + nshift
      DO jc = i_startidx, i_endidx
#endif
        p_temp_ifc(jc,jk1,jb) =                                        &
          temp_ifc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          temp_ifc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          temp_ifc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          temp_ifc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_dpres_mc(jc,jk1,jb) =                                        &
          dpres_mc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          dpres_mc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          dpres_mc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          dpres_mc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_sqv(jc,jk1,jb) =                                        &
          sqv(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          sqv(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          sqv(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          sqv(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_duco2(jc,jk1,jb) =                                        &
          duco2(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          duco2(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          duco2(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          duco2(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_duo3(jc,jk1,jb) =                                        &
          duo3(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          duo3(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          duo3(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          duo3(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_aeq1(jc,jk1,jb) =                                        &
          aeq1(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          aeq1(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          aeq1(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          aeq1(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_aeq2(jc,jk1,jb) =                                        &
          aeq2(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          aeq2(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          aeq2(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          aeq2(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_aeq3(jc,jk1,jb) =                                        &
          aeq3(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          aeq3(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          aeq3(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          aeq3(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_aeq4(jc,jk1,jb) =                                        &
          aeq4(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          aeq4(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          aeq4(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          aeq4(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_aeq5(jc,jk1,jb) =                                        &
          aeq5(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          aeq5(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          aeq5(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          aeq5(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

!CDIR EXPAND=3
        p_tot_cld(jc,jk1,jb,1:3) =                                        &
          tot_cld(iidx(jc,jb,1),jk,iblk(jc,jb,1),1:3)*p_fbkwgt(jc,jb,1) + &
          tot_cld(iidx(jc,jb,2),jk,iblk(jc,jb,2),1:3)*p_fbkwgt(jc,jb,2) + &
          tot_cld(iidx(jc,jb,3),jk,iblk(jc,jb,3),1:3)*p_fbkwgt(jc,jb,3) + &
          tot_cld(iidx(jc,jb,4),jk,iblk(jc,jb,4),1:3)*p_fbkwgt(jc,jb,4)

        p_clc(jc,jk1,jb) =                                        &
          clc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          clc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          clc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          clc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      ENDDO
    ENDDO

    IF (nshift == 1) THEN ! set values for passive top layer if present

      DO jc = i_startidx, i_endidx

        p_dpres_mc(jc,1,jb) =                                            &
          &  pres_ifc(iidx(jc,jb,1),1,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          &  pres_ifc(iidx(jc,jb,2),1,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          &  pres_ifc(iidx(jc,jb,3),1,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          &  pres_ifc(iidx(jc,jb,4),1,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        z_help_pres_ratio = p_dpres_mc(jc,1,jb) / p_dpres_mc(jc,2,jb)

        p_temp_ifc(jc,1,jb) = p_temp_ifc(jc,2,jb)  &
          & + ( p_temp_ifc(jc,2,jb)-p_temp_ifc(jc,3,jb) ) * z_help_pres_ratio

        ! For o3, co2, aerosols and cloud fields, a no-gradient condition is assumed

        p_duo3(jc,1,jb)  =  p_duo3(jc,2,jb)  * z_help_pres_ratio
        p_duco2(jc,1,jb) =  p_duco2(jc,2,jb) * z_help_pres_ratio

        z_rho_1          = ( 0.5_wp*p_dpres_mc(jc,1,jb) ) / ( rd*p_temp_ifc(jc,1,jb) &
          &                *(1.0_wp + vtmpc1*p_tot_cld(jc,2,jb,1) - p_tot_cld(jc,2,jb,2)) )

        p_sqv(jc,1,jb)   =  qsat_rho(p_temp_ifc(jc,1,jb),z_rho_1)
        p_aeq1(jc,1,jb)  =  p_aeq1(jc,2,jb)
        p_aeq2(jc,1,jb)  =  p_aeq2(jc,2,jb)
        p_aeq3(jc,1,jb)  =  p_aeq3(jc,2,jb)
        p_aeq4(jc,1,jb)  =  p_aeq4(jc,2,jb)
        p_aeq5(jc,1,jb)  =  p_aeq5(jc,2,jb)

!CDIR EXPAND=3
        p_tot_cld(jc,1,jb,1:3) = p_tot_cld(jc,2,jb,1:3)

        p_clc(jc,1,jb)   = p_clc(jc,2,jb)
      ENDDO
    ENDIF

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  IF (jgp == 0) THEN

    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 6, 5*nlev_rg+5, &
                            RECV1=rg_temp_ifc, SEND1=z_temp_ifc,            &
                            RECV2=rg_dpres_mc, SEND2=z_dpres_mc,            &
                            RECV3=rg_sqv,      SEND3=z_sqv,                 &
                            RECV4=rg_duco2,    SEND4=z_duco2,               &
                            RECV5=rg_duo3,     SEND5=z_duo3,                &
                            RECV6=zrg_aux3d,   SEND6=z_aux3d              )

    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 6, 6*nlev_rg, &
                            RECV1=rg_aeq1,     SEND1=z_aeq1,              &
                            RECV2=rg_aeq2,     SEND2=z_aeq2,              &
                            RECV3=rg_aeq3,     SEND3=z_aeq3,              &
                            RECV4=rg_aeq4,     SEND4=z_aeq4,              &
                            RECV5=rg_aeq5,     SEND5=z_aeq5,              &
                            RECV6=rg_clc,      SEND6=z_clc                )


    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 3, 3*nlev_rg, &
                            RECV4D=rg_tot_cld, SEND4D=z_tot_cld           )


    i_startblk = p_patch(jgp)%cells%start_blk(1,1)
    i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell,i_nchdom)

    ! OpenMP section commented because the DO loop does almost no work (overhead larger than benefit)
!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch(jgp), jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rlcell)

      DO jc = i_startidx, i_endidx
        rg_cosmu0(jc,jb)    = zrg_aux3d(jc,1,jb)
        rg_albdif(jc,jb)    = zrg_aux3d(jc,2,jb)
        rg_alb_ther(jc,jb)  = zrg_aux3d(jc,3,jb)
        rg_pres_sfc(jc,jb)  = zrg_aux3d(jc,4,jb)
      ENDDO

    ENDDO
!!$OMP END DO NOWAIT
!!$OMP END PARALLEL

    DEALLOCATE(z_cosmu0, z_albdif, z_alb_ther, z_temp_ifc, z_dpres_mc,z_tot_cld, &
      &  z_clc, z_sqv, z_duco2, z_duo3, z_aeq1, z_aeq2, z_aeq3, z_aeq4, z_aeq5,  &
      &  z_pres_sfc, z_aux3d,zrg_aux3d )

  ENDIF

END SUBROUTINE upscale_rad_input_rg

SUBROUTINE downscale_rad_output_rg( jg, jgp, nlev_rg,                &
  &  rg_lwflxall, rg_trsolall, tsfc_rg,                              &
  &  albeff_rg, albefffac_rg,flsp_rg, flsd_rg,                       &
  &  alb_ther_rg, cosmu0_rg, tot_cld_rg,                             &
  &  dpres_mc_rg, pres_sfc_rg,tsfc, albdif, zsct, lwflxall, trsolall )


  ! Input grid parameters
  INTEGER, INTENT(IN)  :: jg, jgp  ! domain IDs of main and reduced grids
  INTEGER, INTENT(IN)  :: nlev_rg  ! number of model levels on reduced grid

  ! Other input fields (on reduced grid)
  REAL(wp), TARGET, INTENT(IN) ::                                               &
    rg_lwflxall(:,:,:), rg_trsolall(:,:,:)

  ! Auxiliary input fields on reduced grid needed for downscaling corrections
  REAL(wp), INTENT(IN), TARGET ::                                      &
    & tsfc_rg(:,:), albeff_rg(:,:), albefffac_rg(:,:),                 &
    & flsp_rg(:,:),flsd_rg(:,:),                                       &
    & alb_ther_rg(:,:), cosmu0_rg(:,:),                                &
    & tot_cld_rg(:,:,:,:), dpres_mc_rg(:,:,:), pres_sfc_rg(:,:)

  ! Auxiliary input fields on full grid needed for downscaling corrections
  REAL(wp), INTENT(IN) :: tsfc(:,:), albdif(:,:)

  REAL(wp), INTENT(IN) :: zsct

  ! Downscaled output fields (on full grid)
  REAL(wp), INTENT(OUT) ::                                                &
    lwflxall(:,:,:), trsolall(:,:,:)

  ! Intermediate storage fields needed in the case of MPI parallelization
  REAL(wp), ALLOCATABLE, TARGET ::                                          &
    z_lwflxall(:,:,:), z_trsolall(:,:,:), z_dpres_mc(:,:,:), z_tot_cld(:,:,:,:)


  ! Pointers to output fields (no MPI) or intermediate fields (MPI)
  REAL(wp), POINTER ::                                                             &
    & p_lwflxall(:,:,:), p_trsolall(:,:,:), p_dpres_mc(:,:,:), p_tot_cld(:,:,:,:)

  ! Additional storage fields to map 2D array(s) to 3D array
  REAL(wp), ALLOCATABLE :: zpg_aux3d(:,:,:), z_aux3d(:,:,:)

  REAL(wp), ALLOCATABLE :: pres_ifc(:,:,:)

  ! Additional storage fields to map 2D array(s) to 3D array
  REAL(wp), ALLOCATABLE, TARGET :: zrg_aux3d(:,:,:)

  ! Auxiliary fields on full grid for back-interpolated values
  REAL(wp), DIMENSION(nproma,p_patch(jg)%nblks_c) :: tsfc_backintp, alb_backintp, albfac_backintp

  !> pressure scale for longwave downscaling correction
  REAL(wp), PARAMETER :: pscal = 1._wp/4000._wp

  ! More local variables
  REAL(wp) :: dpresg, pfaclw, intqctot, dlwflxclr_o_dtg

  REAL(wp), DIMENSION(nproma) :: tqv, dlwem_o_dtg, swfac2, lwfac1, lwfac2

  REAL(wp), DIMENSION(nproma,p_patch(jg)%nlevp1) :: intclw, intcli,     &
    dtrans_o_dalb_all, dlwflxall_o_dtg


  ! Pointers to types needed to minimize code duplication for MPI/no-MPI cases
  TYPE(t_grid_cells), POINTER     :: p_gcp
  TYPE(t_int_state),  POINTER     :: p_int
  TYPE(t_gridref_state), POINTER  :: p_grf
  TYPE(t_patch),      POINTER     :: p_pp

  ! Indices
  INTEGER :: i_chidx, nblks_c_lp
  INTEGER :: jb, jk, jk1, jc, i_startblk, i_endblk, i_startidx, i_endidx
  INTEGER :: jc1, jc2, jc3, jc4, jb1, jb2, jb3, jb4


  INTEGER :: nlev, nlevp1      !< number of full and half levels
  INTEGER :: nshift, nlevp1_rg, nlev_tot
  INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk

  INTEGER :: n2dvars, ipsfc, itsfc, ialb, ialbfac, ialbther, icosmu0, iflsp, iflsd

  LOGICAL :: l_limit(3)
  REAL(wp) :: rlimval(3)
!-----------------------------------------------------------------------

  IF (msg_level >= 10) THEN
    WRITE(message_text,'(a,i2,a,i2)') 'Downscaling of radiation output fields',&
      jgp,' =>',jg
    CALL message('downscale_rad_output',message_text)
  ENDIF

  ! For the time being, the radiation grid is assumed to have the same levels
  ! as the full grid
  nlev   = p_patch(jg)%nlev
  nlevp1 = p_patch(jg)%nlevp1

  ! nlev_rg carries the number of model levels of the reduced grid,
  ! which may be larger than nlev
  nlevp1_rg = nlev_rg + 1
  nshift = nlev_rg - nlev ! resulting shift parameter

  p_grf => p_grf_state_local_parent(jg)
  p_int => p_int_state_local_parent(jg)
  p_gcp => p_patch_local_parent(jg)%cells
  p_pp  => p_patch_local_parent(jg)

  i_chidx  = p_patch(jg)%parent_child_index

  nblks_c_lp = p_pp%nblks_c

  ! pointers to child index/block
  iidx => p_gcp%child_idx
  iblk => p_gcp%child_blk

  ! named constants for accessing 2D variables contained in zrg_aux3d
  n2dvars  = nshift+8
  ipsfc    = nshift+1
  itsfc    = nshift+2
  ialb     = nshift+3
  ialbfac  = nshift+4
  iflsp    = nshift+5
  iflsd    = nshift+6
  ialbther = nshift+7
  icosmu0  = nshift+8

  ! Allocation of local storage fields at local parent level in MPI-case
  IF (jgp == 0) THEN

    ALLOCATE(z_lwflxall(nproma,nlevp1_rg,nblks_c_lp), z_trsolall(nproma,nlevp1_rg,nblks_c_lp) , &
      &      z_dpres_mc(nproma,nlev_rg,nblks_c_lp), z_tot_cld(nproma,nlev_rg,nblks_c_lp,3),     &
      &      zpg_aux3d(nproma,n2dvars,p_patch(jgp)%nblks_c),                                    &
      &      zrg_aux3d(nproma,n2dvars,nblks_c_lp),   z_aux3d(nproma,8,p_patch(jg)%nblks_c),     &
      &      pres_ifc(nproma,nlevp1_rg,nblks_c_lp) )

    ! Perform communication from parent to local parent grid in the MPI case,
    ! and set pointers such that further processing is the same for MPI / non-MPI cases

    IF (nshift > 0) zpg_aux3d(:,1:nshift,:) = 0._wp
    zpg_aux3d(:,ipsfc,:)  = pres_sfc_rg(:,:)
    zpg_aux3d(:,itsfc,:)  = tsfc_rg(:,:)
    zpg_aux3d(:,ialb,:)   = albeff_rg(:,:)
    zpg_aux3d(:,ialbfac,:)= albefffac_rg(:,:)
    zpg_aux3d(:,iflsp,:)  = flsp_rg(:,:)
    zpg_aux3d(:,iflsd,:)  = flsd_rg(:,:)
    zpg_aux3d(:,ialbther,:)  = alb_ther_rg(:,:)
    zpg_aux3d(:,icosmu0,:)= cosmu0_rg(:,:)

    nlev_tot = 2*nlevp1_rg + nlev_rg + n2dvars


    CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, 4, nlev_tot, &
                            RECV1=z_lwflxall, SEND1=rg_lwflxall,     &
                            RECV2=z_trsolall, SEND2=rg_trsolall,     &
                            RECV3=z_dpres_mc, SEND3=dpres_mc_rg,     &
                            RECV4=zrg_aux3d , SEND4=zpg_aux3d        )

    CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, 3, 3*nlev_rg, &
                            RECV4D=z_tot_cld, SEND4D=tot_cld_rg       )



    p_lwflxall   => z_lwflxall
    p_trsolall   => z_trsolall
    p_dpres_mc   => z_dpres_mc
    p_tot_cld    => z_tot_cld

  ELSE

    ALLOCATE(zrg_aux3d(nproma,n2dvars,nblks_c_lp),z_aux3d(nproma,8,p_patch(jg)%nblks_c), &
      &      pres_ifc(nproma,nlevp1_rg,nblks_c_lp) )

    p_lwflxall   => rg_lwflxall
    p_trsolall   => rg_trsolall
    p_dpres_mc   => dpres_mc_rg
    p_tot_cld    => tot_cld_rg

    IF (nshift > 0) zrg_aux3d(:,1:nshift,:) = 0._wp
    zrg_aux3d(:,ipsfc,:)  = pres_sfc_rg(:,:)
    zrg_aux3d(:,itsfc,:)  = tsfc_rg(:,:)
    zrg_aux3d(:,ialb,:)   = albeff_rg(:,:)
    zrg_aux3d(:,ialbfac,:)= albefffac_rg(:,:)
    zrg_aux3d(:,iflsp,:)   = flsp_rg(:,:)
    zrg_aux3d(:,iflsd,:)   = flsd_rg(:,:)
    zrg_aux3d(:,ialbther,:)  = alb_ther_rg(:,:)
    zrg_aux3d(:,icosmu0,:)= cosmu0_rg(:,:)

  ENDIF

  ! Interpolate reduced-grid fields to full grid

  ! interpol_vec/scal_nudging needs all fields with full boundaries, so the arrays
  ! have to be sync'd before calling these routines.
  ! Please note that we cannot use sync_patch_array here (comparing parallel/non parallel results)
  ! since the arrays don't start with lower bound 1 in the non paralellel case!

  nlev_tot = 2*nlevp1_rg + n2dvars

  IF (.NOT. my_process_is_mpi_seq()) THEN
    CALL exchange_data_mult(p_pp%comm_pat_c, 3 ,nlev_tot, recv1=p_lwflxall, recv2=p_trsolall, &
      &                     recv3=zrg_aux3d )
  END IF

  IF (p_test_run) THEN
    trsolall = 0._wp
    lwflxall = 0._wp
    z_aux3d  = 0._wp
  ENDIF

  l_limit(1)   = .TRUE.      ! limit transmissivity to positive values
  l_limit(2:3) = .FALSE.
  rlimval(:)   = 2.94e-37_wp ! seems to be the lower threshold for SW transmissivity

  CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, nshift, 3, 1, &
    &                         f3din1=p_trsolall, f3dout1=trsolall,                      &
    &                         f3din2=p_lwflxall, f3dout2=lwflxall,                      &
    &                         f3din3=zrg_aux3d,  f3dout3=z_aux3d,                       &
                              llimit_nneg=l_limit, rlimval=rlimval, overshoot_fac=1.1_wp)

  tsfc_backintp(:,:)   = z_aux3d(:,itsfc-nshift,:)
  alb_backintp(:,:)    = z_aux3d(:,ialb-nshift,:)
  albfac_backintp(:,:) = z_aux3d(:,ialbfac-nshift,:)

  ! Finally apply empirical downscaling corrections depending on variations
  ! in ground temperature and albedo

  ! Start/End block in the parent domain
  i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
  i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jk1,tqv,intclw,intcli,dpresg,                   &
!$OMP            dlwem_o_dtg,swfac2,lwfac1,lwfac2,dtrans_o_dalb_all,                          &
!$OMP            pfaclw,intqctot,dlwflxclr_o_dtg,dlwflxall_o_dtg,jc1,jc2,jc3,jc4,jb1,jb2,jb3, &
!$OMP jb4) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,                      &
                       i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int)

    tqv(:)            = 0._wp
    intclw(:,nlevp1)  = 0._wp
    intcli(:,nlevp1)  = 0._wp
    pres_ifc(:,1,jb)     = 0._wp
    pres_ifc(:,nlevp1_rg,jb) = zrg_aux3d(:,ipsfc,jb)

    DO jk = nlev,1,-1
      jk1 = jk + nshift
      DO jc = i_startidx, i_endidx
        pres_ifc(jc,jk1,jb) = pres_ifc(jc,jk1+1,jb) - p_dpres_mc(jc,jk1,jb)

        dpresg        = (pres_ifc(jc,jk1+1,jb) - pres_ifc(jc,jk1,jb))/grav
        tqv(jc)       = tqv(jc)+p_tot_cld(jc,jk1,jb,iqv)*dpresg
        intclw(jc,jk) = intclw(jc,jk+1)+p_tot_cld(jc,jk1,jb,iqc)*dpresg
        intcli(jc,jk) = intcli(jc,jk+1)+p_tot_cld(jc,jk1,jb,iqi)*dpresg
      ENDDO
    ENDDO

    DO jc = i_startidx, i_endidx
      dlwem_o_dtg(jc) = ( 1._wp-zrg_aux3d(jc,ialbther,jb) )*4._wp*stbo*zrg_aux3d(jc,itsfc,jb)**3
      swfac2(jc) =  MAX(0.25_wp,3._wp*zrg_aux3d(jc,icosmu0,jb))**0.0677988_wp
      IF (tqv(jc) > 15._wp) THEN
        lwfac1(jc) = 1.677_wp*MAX(1._wp,tqv(jc))**(-0.72_wp)
      ELSE
        lwfac1(jc) = 0.4388_wp*MAX(1._wp,tqv(jc))**(-0.225_wp)
      ENDIF
      lwfac2(jc) = 0.92_wp*MAX(1._wp,tqv(jc))**(-0.07_wp)
    ENDDO

    DO jk = 1,nlevp1
      jk1 = jk + nshift
      DO jc = i_startidx, i_endidx

        dtrans_o_dalb_all(jc,jk) = (-zrg_aux3d(jc,iflsp,jb)-zrg_aux3d(jc,iflsd,jb)) &
          &                      / ( zsct * zrg_aux3d(jc,icosmu0,jb) * swfac2(jc) ) &
          &                      * ( 0.482994_wp + 0.436556_wp                      &
          &                      * (zrg_aux3d(jc,iflsp,jb)/zrg_aux3d(jc,iflsd,jb))**0.0953119_wp )

        pfaclw = lwfac1(jc)+(lwfac2(jc)-lwfac1(jc))*EXP(-SQRT((pres_ifc(jc,nlevp1_rg,jb)- &
                 pres_ifc(jc,jk1,jb))*pscal ))
        intqctot = MIN(0.30119_wp,MAX(1.008e-3_wp,intclw(jc,jk)+0.2_wp*intcli(jc,jk)))

        dlwflxclr_o_dtg        = -dlwem_o_dtg(jc)*pfaclw
        dlwflxall_o_dtg(jc,jk) = dlwflxclr_o_dtg*(1._wp-(6.9_wp+LOG(intqctot))/5.7_wp)
      ENDDO
    ENDDO

    ! Now apply the corrections
#ifdef __LOOP_EXCHANGE
    DO jc = i_startidx, i_endidx
      jc1 = iidx(jc,jb,1)
      jc2 = iidx(jc,jb,2)
      jc3 = iidx(jc,jb,3)
      jc4 = iidx(jc,jb,4)
      jb1 = iblk(jc,jb,1)
      jb2 = iblk(jc,jb,2)
      jb3 = iblk(jc,jb,3)
      jb4 = iblk(jc,jb,4)
      DO jk = 1,nlevp1
#else
!CDIR NOLOOPCHG
    DO jk = 1,nlevp1
!CDIR NODEP,VOVERTAKE,VOB
      DO jc = i_startidx, i_endidx
        jc1 = iidx(jc,jb,1)
        jc2 = iidx(jc,jb,2)
        jc3 = iidx(jc,jb,3)
        jc4 = iidx(jc,jb,4)
        jb1 = iblk(jc,jb,1)
        jb2 = iblk(jc,jb,2)
        jb3 = iblk(jc,jb,3)
        jb4 = iblk(jc,jb,4)
#endif

        trsolall(jc1,jk,jb1) = MAX(trsolall(jc1,jk,jb1) + dtrans_o_dalb_all(jc,jk)* &
          ( albfac_backintp(jc1,jb1)*albdif(jc1,jb1) - alb_backintp(jc1,jb1) ), rlimval(1))
        trsolall(jc2,jk,jb2) = MAX(trsolall(jc2,jk,jb2) + dtrans_o_dalb_all(jc,jk)* &
          ( albfac_backintp(jc2,jb2)*albdif(jc2,jb2) - alb_backintp(jc2,jb2) ), rlimval(1))
        trsolall(jc3,jk,jb3) = MAX(trsolall(jc3,jk,jb3) + dtrans_o_dalb_all(jc,jk)* &
          ( albfac_backintp(jc3,jb3)*albdif(jc3,jb3) - alb_backintp(jc3,jb3) ), rlimval(1))
        trsolall(jc4,jk,jb4) = MAX(trsolall(jc4,jk,jb4) + dtrans_o_dalb_all(jc,jk)* &
          ( albfac_backintp(jc4,jb4)*albdif(jc4,jb4) - alb_backintp(jc4,jb4) ), rlimval(1))

        lwflxall(jc1,jk,jb1) = lwflxall(jc1,jk,jb1) + dlwflxall_o_dtg(jc,jk)* &
          ( tsfc(jc1,jb1) - tsfc_backintp(jc1,jb1) )
        lwflxall(jc2,jk,jb2) = lwflxall(jc2,jk,jb2) + dlwflxall_o_dtg(jc,jk)* &
          ( tsfc(jc2,jb2) - tsfc_backintp(jc2,jb2) )
        lwflxall(jc3,jk,jb3) = lwflxall(jc3,jk,jb3) + dlwflxall_o_dtg(jc,jk)* &
          ( tsfc(jc3,jb3) - tsfc_backintp(jc3,jb3) )
        lwflxall(jc4,jk,jb4) = lwflxall(jc4,jk,jb4) + dlwflxall_o_dtg(jc,jk)* &
          ( tsfc(jc4,jb4) - tsfc_backintp(jc4,jb4) )

      ENDDO
    ENDDO

  ENDDO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  IF (jgp == 0) THEN
    DEALLOCATE( z_lwflxall, z_trsolall, z_dpres_mc, z_tot_cld, zpg_aux3d )
  ENDIF

  DEALLOCATE( zrg_aux3d,z_aux3d, pres_ifc )

END SUBROUTINE downscale_rad_output_rg



!>
!! This routine optimizes the boundary interpolation of diagnostic physics fields for output
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2010-12-03
!!
SUBROUTINE interpol_phys_grf (ext_data, jg, jgc, jn)

  USE mo_nwp_phy_state,      ONLY: prm_diag
  USE mo_nonhydro_state,     ONLY: p_nh_state

  ! Input:
  TYPE(t_external_data), INTENT(in) :: ext_data(:)
  INTEGER              , INTENT(in) :: jg,jgc,jn

  ! Pointers
  TYPE(t_patch),                POINTER :: ptr_pp
  TYPE(t_patch),                POINTER :: ptr_pc
  TYPE(t_gridref_single_state), POINTER :: ptr_grf
  TYPE(t_lnd_diag),             POINTER :: ptr_ldiagp ! parent level land diag state
  TYPE(t_lnd_diag),             POINTER :: ptr_ldiagc ! child level land diag state
  TYPE(t_lnd_prog),             POINTER :: ptr_lprogp ! parent level land prog state
  TYPE(t_lnd_prog),             POINTER :: ptr_lprogc ! child level land prog state
  TYPE(t_wtr_prog),             POINTER :: ptr_wprogp ! parent level water prog state
  TYPE(t_wtr_prog),             POINTER :: ptr_wprogc ! child level water prog state

  ! Local fields
  INTEGER, PARAMETER  :: nfields_p1=66   ! Number of positive-definite 2D physics fields for which boundary interpolation is needed
  INTEGER, PARAMETER  :: nfields_p2=18   ! Number of remaining 2D physics fields for which boundary interpolation is needed
  INTEGER, PARAMETER  :: nfields_l2=18   ! Number of 2D land state fields

  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, jb, jc, jk, jt, nlev_c, ic, i_count, indlist(nproma)
  INTEGER :: styp                        ! soiltype at child level

  LOGICAL :: lsfc_interp

  ! Temporary storage to do boundary interpolation for all 2D fields in one step
  REAL(wp) :: z_aux3dp1_p(nproma,nfields_p1,p_patch(jg)%nblks_c),       &  ! 2D pos. def. physics diag fields
    &         z_aux3dp1_c(nproma,nfields_p1,p_patch(jgc)%nblks_c),      &
    &         z_aux3dp2_p(nproma,nfields_p2,p_patch(jg)%nblks_c),       &  ! 2D physics diag fields
    &         z_aux3dp2_c(nproma,nfields_p2,p_patch(jgc)%nblks_c),      &
    &         z_aux3dl2_p(nproma,nfields_l2,p_patch(jg)%nblks_c),       &  ! 2D land state fields
    &         z_aux3dl2_c(nproma,nfields_l2,p_patch(jgc)%nblks_c),      &
    &         z_aux3dso_p(nproma,3*nlev_soil,p_patch(jg)%nblks_c),      &  ! 3D land state fields for soil
    &         z_aux3dso_c(nproma,3*nlev_soil,p_patch(jgc)%nblks_c),     &
    &         z_aux3dsn_p(nproma,5*nlev_snow,p_patch(jg)%nblks_c),      &  ! 3D land state fields for multi-layer snow
    &         z_aux3dsn_c(nproma,5*nlev_snow,p_patch(jgc)%nblks_c)         ! (used if lmulti_snow = true))

  ! set pointers
  ptr_pp  => p_patch(jg)
  ptr_pc  => p_patch(jgc)
  ptr_grf => p_grf_state(jg)%p_dom(jn)

  nlev_c = ptr_pc%nlev

  IF (atm_phy_nwp_config(jg)%inwp_surface == 1) THEN
    lsfc_interp = .TRUE.
  ELSE
    lsfc_interp = .FALSE.
  ENDIF

  IF (lsfc_interp) THEN
    ptr_ldiagp => p_lnd_state(jg)%diag_lnd
    ptr_ldiagc => p_lnd_state(jgc)%diag_lnd
    ptr_lprogp => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
    ptr_lprogc => p_lnd_state(jgc)%prog_lnd(nnow_rcf(jgc))
    ptr_wprogp => p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))
    ptr_wprogc => p_lnd_state(jgc)%prog_wtr(nnow_rcf(jgc))
  ENDIF

  IF (p_test_run) THEN
     z_aux3dp1_p(:,:,:) = 0._wp
     z_aux3dp2_p(:,:,:) = 0._wp
     z_aux3dl2_p(:,:,:) = 0._wp
     z_aux3dso_p(:,:,:) = 0._wp
     z_aux3dsn_p(:,:,:) = 0._wp
  ENDIF

  i_startblk = ptr_pp%cells%start_block(grf_bdywidth_c+1)
  i_endblk   = ptr_pp%cells%end_block(min_rlcell_int)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,ic,i_count) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell_int)

    DO jc = i_startidx, i_endidx
      z_aux3dp1_p(jc,1,jb) = prm_diag(jg)%tot_prec(jc,jb)
      z_aux3dp1_p(jc,2,jb) = prm_diag(jg)%rain_gsp(jc,jb)
      z_aux3dp1_p(jc,3,jb) = prm_diag(jg)%snow_gsp(jc,jb)
      z_aux3dp1_p(jc,4,jb) = prm_diag(jg)%rain_con(jc,jb)
      z_aux3dp1_p(jc,5,jb) = prm_diag(jg)%snow_con(jc,jb)
      z_aux3dp1_p(jc,6,jb) = prm_diag(jg)%rain_gsp_rate(jc,jb)
      z_aux3dp1_p(jc,7,jb) = prm_diag(jg)%snow_gsp_rate(jc,jb)
      z_aux3dp1_p(jc,8,jb) = prm_diag(jg)%rain_con_rate(jc,jb)
      z_aux3dp1_p(jc,9,jb) = prm_diag(jg)%snow_con_rate(jc,jb)
      z_aux3dp1_p(jc,10,jb) = prm_diag(jg)%gz0(jc,jb)
      z_aux3dp1_p(jc,11,jb) = prm_diag(jg)%tcm(jc,jb)
      z_aux3dp1_p(jc,12,jb) = prm_diag(jg)%tch(jc,jb)
      z_aux3dp1_p(jc,13,jb) = prm_diag(jg)%tfm(jc,jb)
      z_aux3dp1_p(jc,14,jb) = prm_diag(jg)%tfh(jc,jb)
      z_aux3dp1_p(jc,15,jb) = prm_diag(jg)%tfv(jc,jb)
      z_aux3dp1_p(jc,16,jb) = prm_diag(jg)%t_2m(jc,jb)
      z_aux3dp1_p(jc,17,jb) = prm_diag(jg)%qv_2m(jc,jb)
      z_aux3dp1_p(jc,18,jb) = prm_diag(jg)%td_2m(jc,jb)
      z_aux3dp1_p(jc,19,jb) = prm_diag(jg)%rh_2m(jc,jb)
      z_aux3dp1_p(jc,20,jb) = prm_diag(jg)%gust10(jc,jb)
      z_aux3dp1_p(jc,21,jb) = prm_diag(jg)%sp_10m(jc,jb)
      z_aux3dp1_p(jc,22,jb) = prm_diag(jg)%swflxsfc(jc,jb)
      z_aux3dp1_p(jc,23,jb) = prm_diag(jg)%swflx_dn_sfc_diff(jc,jb)
      z_aux3dp1_p(jc,24,jb) = prm_diag(jg)%lwflx_up_sfc(jc,jb)
      z_aux3dp1_p(jc,25,jb) = prm_diag(jg)%swflxtoa(jc,jb)
      z_aux3dp1_p(jc,26,jb) = prm_diag(jg)%flxdwswtoa(jc,jb)
      z_aux3dp1_p(jc,27,jb) = prm_diag(jg)%swflxsfc_a(jc,jb)
      z_aux3dp1_p(jc,28,jb) = prm_diag(jg)%asodifd_s(jc,jb)
      z_aux3dp1_p(jc,29,jb) = prm_diag(jg)%asodifu_s(jc,jb)
      z_aux3dp1_p(jc,30,jb) = prm_diag(jg)%athu_s(jc,jb)
      z_aux3dp1_p(jc,31,jb) = prm_diag(jg)%athd_s(jc,jb)
      z_aux3dp1_p(jc,32,jb) = prm_diag(jg)%swflxtoa_a(jc,jb)
      z_aux3dp1_p(jc,33,jb) = prm_diag(jg)%asod_t(jc,jb)
      z_aux3dp1_p(jc,34,jb) = prm_diag(jg)%asou_t(jc,jb)
      z_aux3dp1_p(jc,35,jb) = prm_diag(jg)%asodird_s(jc,jb)
      z_aux3dp1_p(jc,36,jb) = prm_diag(jg)%htop_con(jc,jb) - prm_diag(jg)%hbas_con(jc,jb)
      z_aux3dp1_p(jc,37,jb) = prm_diag(jg)%htop_dc(jc,jb)
      z_aux3dp1_p(jc,38,jb) = prm_diag(jg)%snowlmt(jc,jb) + 999._wp
      z_aux3dp1_p(jc,39,jb) = prm_diag(jg)%hzerocl(jc,jb) + 999._wp
      z_aux3dp1_p(jc,40,jb) = prm_diag(jg)%clcl(jc,jb)
      z_aux3dp1_p(jc,41,jb) = prm_diag(jg)%clcm(jc,jb)
      z_aux3dp1_p(jc,42,jb) = prm_diag(jg)%clch(jc,jb)
      z_aux3dp1_p(jc,43,jb) = prm_diag(jg)%clct(jc,jb)
      z_aux3dp1_p(jc,44,jb) = prm_diag(jg)%cape(jc,jb)
      z_aux3dp1_p(jc,45,jb) = prm_diag(jg)%swflx_up_toa(jc,jb)
      z_aux3dp1_p(jc,46,jb) = prm_diag(jg)%swflx_up_sfc(jc,jb)
      z_aux3dp1_p(jc,47,jb) = prm_diag(jg)%tmax_2m(jc,jb)
      z_aux3dp1_p(jc,48,jb) = prm_diag(jg)%tmin_2m(jc,jb)
      z_aux3dp1_p(jc,49:51,jb) = prm_diag(jg)%tot_cld_vi(jc,jb,1:3)
      z_aux3dp1_p(jc,52:56,jb) = p_nh_state(jg)%diag%tracer_vi(jc,jb,1:5)
      z_aux3dp1_p(jc,57,jb) = prm_diag(jg)%clct_mod(jc,jb)

      IF (atm_phy_nwp_config(jg)%inwp_gscp == 2) THEN
        z_aux3dp1_p(jc,58,jb) = prm_diag(jg)%graupel_gsp(jc,jb)
        z_aux3dp1_p(jc,59,jb) = prm_diag(jg)%graupel_gsp_rate(jc,jb)
      ELSE
        z_aux3dp1_p(jc,58,jb) = 0._wp
        z_aux3dp1_p(jc,59,jb) = 0._wp
      ENDIF
      z_aux3dp1_p(jc,60,jb) = prm_diag(jg)%tvm(jc,jb)
      z_aux3dp1_p(jc,61,jb) = prm_diag(jg)%tvh(jc,jb)
      z_aux3dp1_p(jc,62,jb) = prm_diag(jg)%tkr(jc,jb)
      z_aux3dp1_p(jc,63,jb) = prm_diag(jg)%cldepth(jc,jb)
      z_aux3dp1_p(jc,64,jb) = prm_diag(jg)%t_2m_land(jc,jb)
      z_aux3dp1_p(jc,65,jb) = prm_diag(jg)%td_2m_land(jc,jb)
      z_aux3dp1_p(jc,66,jb) = prm_diag(jg)%rh_2m_land(jc,jb)

      z_aux3dp2_p(jc,1,jb) = prm_diag(jg)%u_10m(jc,jb)
      z_aux3dp2_p(jc,2,jb) = prm_diag(jg)%v_10m(jc,jb)
      z_aux3dp2_p(jc,3,jb) = prm_diag(jg)%lhfl_s(jc,jb)
      z_aux3dp2_p(jc,4,jb) = prm_diag(jg)%lhfl_bs(jc,jb)
      z_aux3dp2_p(jc,5,jb) = prm_diag(jg)%shfl_s(jc,jb)
      z_aux3dp2_p(jc,6,jb) = prm_diag(jg)%qhfl_s(jc,jb)
      z_aux3dp2_p(jc,7,jb) = prm_diag(jg)%umfl_s(jc,jb)
      z_aux3dp2_p(jc,8,jb) = prm_diag(jg)%vmfl_s(jc,jb)
      z_aux3dp2_p(jc,9,jb) = prm_diag(jg)%alhfl_s(jc,jb)
      z_aux3dp2_p(jc,10,jb) = prm_diag(jg)%alhfl_bs(jc,jb)
      z_aux3dp2_p(jc,11,jb) = prm_diag(jg)%ashfl_s(jc,jb)
      z_aux3dp2_p(jc,12,jb) = prm_diag(jg)%aqhfl_s(jc,jb)
      z_aux3dp2_p(jc,13,jb) = prm_diag(jg)%aumfl_s(jc,jb)
      z_aux3dp2_p(jc,14,jb) = prm_diag(jg)%avmfl_s(jc,jb)
      z_aux3dp2_p(jc,15,jb) = prm_diag(jg)%lwflxsfc(jc,jb)
      z_aux3dp2_p(jc,16,jb) = prm_diag(jg)%lwflxsfc_a(jc,jb)
      z_aux3dp2_p(jc,17,jb) = prm_diag(jg)%lwflxtoa_a(jc,jb)
      z_aux3dp2_p(jc,18,jb) = prm_diag(jg)%hbas_con(jc,jb)
    ENDDO

    IF (lsfc_interp) THEN
      DO jc = i_startidx, i_endidx
        z_aux3dl2_p(jc,1,jb) = ptr_ldiagp%t_snow(jc,jb)
        z_aux3dl2_p(jc,2,jb) = ptr_ldiagp%t_s(jc,jb)
        z_aux3dl2_p(jc,3,jb) = ptr_ldiagp%w_snow(jc,jb)
        z_aux3dl2_p(jc,4,jb) = ptr_ldiagp%rho_snow(jc,jb)
        z_aux3dl2_p(jc,5,jb) = ptr_ldiagp%w_i(jc,jb)
        z_aux3dl2_p(jc,6,jb) = ptr_ldiagp%h_snow(jc,jb)
        z_aux3dl2_p(jc,7,jb) = ptr_ldiagp%freshsnow(jc,jb)
        z_aux3dl2_p(jc,8,jb) = ptr_ldiagp%snowfrac(jc,jb)
        z_aux3dl2_p(jc,9,jb) = ptr_ldiagp%runoff_s(jc,jb)
        z_aux3dl2_p(jc,10,jb) = ptr_ldiagp%runoff_g(jc,jb)
        z_aux3dl2_p(jc,11,jb) = ptr_ldiagp%t_so(jc,nlev_soil+1,jb)
        IF (lmulti_snow) THEN
          z_aux3dl2_p(jc,12,jb) = ptr_ldiagp%t_snow_mult(jc,nlev_snow+1,jb)
        ELSE
          z_aux3dl2_p(jc,12,jb) = 0._wp
        ENDIF
      ENDDO

      IF (lseaice) THEN
        DO jc = i_startidx, i_endidx
          IF (ptr_wprogp%t_ice(jc,jb) > 10._wp) THEN
            z_aux3dl2_p(jc,13,jb) = ptr_wprogp%t_ice(jc,jb)
          ELSE
            z_aux3dl2_p(jc,13,jb) = ptr_lprogp%t_g(jc,jb)
          ENDIF
          z_aux3dl2_p(jc,14,jb) = ptr_wprogp%h_ice(jc,jb)
          z_aux3dl2_p(jc,15,jb) = ptr_wprogp%t_snow_si(jc,jb)
          z_aux3dl2_p(jc,16,jb) = ptr_wprogp%h_snow_si(jc,jb)
          z_aux3dl2_p(jc,17,jb) = ptr_ldiagp%fr_seaice(jc,jb)
        ENDDO
      ELSE
        z_aux3dl2_p(:,13:17,jb) = 0._wp
      ENDIF

      IF (llake) THEN
        DO jc = i_startidx, i_endidx
          z_aux3dl2_p(jc,18,jb) = ptr_lprogp%t_g(jc,jb)
        ENDDO
        i_count = ext_data(jg)%atm%fp_count(jb)
        DO ic = 1, i_count
          jc = ext_data(jg)%atm%idx_lst_fp(ic,jb)
          z_aux3dl2_p(jc,18,jb) = ptr_lprogp%t_g_t(jc,jb,isub_lake)
        ENDDO
      ELSE
        z_aux3dl2_p(:,18,jb) = 0._wp
      ENDIF

      DO jk = 1, nlev_soil
        DO jc = i_startidx, i_endidx
          z_aux3dso_p(jc,3*(jk-1)+1,jb) = ptr_ldiagp%t_so(jc,jk,jb)
          z_aux3dso_p(jc,3*(jk-1)+2,jb) = ptr_ldiagp%w_so(jc,jk,jb)
          z_aux3dso_p(jc,3*(jk-1)+3,jb) = ptr_ldiagp%w_so_ice(jc,jk,jb)
        ENDDO
      ENDDO

      IF (lmulti_snow) THEN
        DO jk = 1, nlev_snow
          DO jc = i_startidx, i_endidx
            z_aux3dsn_p(jc,5*(jk-1)+1,jb) = ptr_ldiagp%t_snow_mult(jc,jk,jb)
            z_aux3dsn_p(jc,5*(jk-1)+2,jb) = ptr_ldiagp%rho_snow_mult(jc,jk,jb)
            z_aux3dsn_p(jc,5*(jk-1)+3,jb) = ptr_ldiagp%wliq_snow(jc,jk,jb)
            z_aux3dsn_p(jc,5*(jk-1)+4,jb) = ptr_ldiagp%wtot_snow(jc,jk,jb)
            z_aux3dsn_p(jc,5*(jk-1)+5,jb) = ptr_ldiagp%dzh_snow(jc,jk,jb)
          ENDDO
        ENDDO
      ENDIF

    ENDIF

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! Halo update is needed before interpolation
    IF (lsfc_interp .AND. lmulti_snow) THEN

      CALL sync_patch_array_mult(SYNC_C,ptr_pp,5,z_aux3dp1_p,z_aux3dp2_p,z_aux3dl2_p,z_aux3dso_p,z_aux3dsn_p)
      CALL interpol_scal_grf (ptr_pp, ptr_pc, ptr_grf, 5, z_aux3dp1_p, z_aux3dp1_c, z_aux3dp2_p, z_aux3dp2_c,&
        z_aux3dl2_p, z_aux3dl2_c, z_aux3dso_p, z_aux3dso_c, z_aux3dsn_p, z_aux3dsn_c, &
        llimit_nneg=(/.TRUE.,.FALSE.,.TRUE.,.TRUE.,.TRUE./), lnoshift=.TRUE.)

    ELSE IF (lsfc_interp) THEN

      CALL sync_patch_array_mult(SYNC_C,ptr_pp,4,z_aux3dp1_p,z_aux3dp2_p,z_aux3dl2_p,z_aux3dso_p)
      CALL interpol_scal_grf (ptr_pp, ptr_pc, ptr_grf, 4, z_aux3dp1_p, z_aux3dp1_c, z_aux3dp2_p, z_aux3dp2_c,&
        z_aux3dl2_p, z_aux3dl2_c, z_aux3dso_p, z_aux3dso_c,                       &
        llimit_nneg=(/.TRUE.,.FALSE.,.TRUE.,.TRUE./), lnoshift=.TRUE.)

    ELSE
      CALL sync_patch_array_mult(SYNC_C,ptr_pp,2,z_aux3dp1_p,z_aux3dp2_p)
      CALL interpol_scal_grf (ptr_pp, ptr_pc, ptr_grf, 2, z_aux3dp1_p, z_aux3dp1_c, z_aux3dp2_p, z_aux3dp2_c, &
        llimit_nneg=(/.TRUE.,.FALSE./), lnoshift=.TRUE.)

    ENDIF

    CALL sync_patch_array_mult(SYNC_C,ptr_pp,3,prm_diag(jg)%tkvm,prm_diag(jg)%tkvh,prm_diag(jg)%rcld)
    CALL interpol_scal_grf (ptr_pp, ptr_pc, ptr_grf, 3, prm_diag(jg)%tkvm, prm_diag(jgc)%tkvm, &
      prm_diag(jg)%tkvh, prm_diag(jgc)%tkvh, prm_diag(jg)%rcld, prm_diag(jgc)%rcld,            &
      llimit_nneg=(/.TRUE.,.TRUE.,.TRUE./))

  i_startblk = ptr_pc%cells%start_blk(1,1)
  i_endblk   = ptr_pc%cells%end_blk(grf_bdywidth_c,1)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jt,styp,ic,i_count,indlist) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_pc, jb, i_startblk, i_endblk,        &
                       i_startidx, i_endidx, 1, grf_bdywidth_c)

    DO jc = i_startidx, i_endidx          ! to avoid undershoots when taking time differences:
      prm_diag(jgc)%tot_prec(jc,jb)       = MAX(z_aux3dp1_c(jc,1,jb),prm_diag(jgc)%tot_prec(jc,jb))
      prm_diag(jgc)%rain_gsp(jc,jb)       = MAX(z_aux3dp1_c(jc,2,jb),prm_diag(jgc)%rain_gsp(jc,jb))
      prm_diag(jgc)%snow_gsp(jc,jb)       = MAX(z_aux3dp1_c(jc,3,jb),prm_diag(jgc)%snow_gsp(jc,jb))
      prm_diag(jgc)%rain_con(jc,jb)       = MAX(z_aux3dp1_c(jc,4,jb),prm_diag(jgc)%rain_con(jc,jb))
      prm_diag(jgc)%snow_con(jc,jb)       = MAX(z_aux3dp1_c(jc,5,jb),prm_diag(jgc)%snow_con(jc,jb))
      prm_diag(jgc)%rain_gsp_rate(jc,jb)  = z_aux3dp1_c(jc,6,jb)
      prm_diag(jgc)%snow_gsp_rate(jc,jb)  = z_aux3dp1_c(jc,7,jb)
      prm_diag(jgc)%rain_con_rate(jc,jb)  = z_aux3dp1_c(jc,8,jb)
      prm_diag(jgc)%snow_con_rate(jc,jb)  = z_aux3dp1_c(jc,9,jb)
      prm_diag(jgc)%gz0(jc,jb)            = z_aux3dp1_c(jc,10,jb)
      prm_diag(jgc)%tcm(jc,jb)            = z_aux3dp1_c(jc,11,jb)
      prm_diag(jgc)%tch(jc,jb)            = z_aux3dp1_c(jc,12,jb)
      prm_diag(jgc)%tfm(jc,jb)            = z_aux3dp1_c(jc,13,jb)
      prm_diag(jgc)%tfh(jc,jb)            = z_aux3dp1_c(jc,14,jb)
      prm_diag(jgc)%tfv(jc,jb)            = z_aux3dp1_c(jc,15,jb)
      prm_diag(jgc)%t_2m(jc,jb)           = z_aux3dp1_c(jc,16,jb)
      prm_diag(jgc)%qv_2m(jc,jb)          = z_aux3dp1_c(jc,17,jb)
      prm_diag(jgc)%td_2m(jc,jb)          = MIN(prm_diag(jgc)%t_2m(jc,jb),z_aux3dp1_c(jc,18,jb))
      prm_diag(jgc)%rh_2m(jc,jb)          = MIN(100._wp,z_aux3dp1_c(jc,19,jb)) ! unit is %
      prm_diag(jgc)%gust10(jc,jb)         = z_aux3dp1_c(jc,20,jb)
      prm_diag(jgc)%sp_10m(jc,jb)         = z_aux3dp1_c(jc,21,jb)
      prm_diag(jgc)%swflxsfc(jc,jb)       = z_aux3dp1_c(jc,22,jb)
      prm_diag(jgc)%swflx_dn_sfc_diff(jc,jb) = z_aux3dp1_c(jc,23,jb)
      prm_diag(jgc)%lwflx_up_sfc(jc,jb)   = z_aux3dp1_c(jc,24,jb)
      prm_diag(jgc)%swflxtoa(jc,jb)       = z_aux3dp1_c(jc,25,jb)
      prm_diag(jgc)%flxdwswtoa(jc,jb)     = z_aux3dp1_c(jc,26,jb)
      prm_diag(jgc)%swflxsfc_a(jc,jb)     = z_aux3dp1_c(jc,27,jb)
      prm_diag(jgc)%asodifd_s(jc,jb)      = z_aux3dp1_c(jc,28,jb)
      prm_diag(jgc)%asodifu_s(jc,jb)      = z_aux3dp1_c(jc,29,jb)
      prm_diag(jgc)%athu_s(jc,jb)         = z_aux3dp1_c(jc,30,jb)
      prm_diag(jgc)%athd_s(jc,jb)         = z_aux3dp1_c(jc,31,jb)
      prm_diag(jgc)%swflxtoa_a(jc,jb)     = z_aux3dp1_c(jc,32,jb)
      prm_diag(jgc)%asod_t(jc,jb)         = z_aux3dp1_c(jc,33,jb)
      prm_diag(jgc)%asou_t(jc,jb)         = z_aux3dp1_c(jc,34,jb)
      prm_diag(jgc)%asodird_s(jc,jb)      = z_aux3dp1_c(jc,35,jb)
      prm_diag(jgc)%htop_dc(jc,jb)        = z_aux3dp1_c(jc,37,jb)
      prm_diag(jgc)%snowlmt(jc,jb)        = z_aux3dp1_c(jc,38,jb) - 999._wp
      IF (prm_diag(jgc)%snowlmt(jc,jb) < p_nh_state(jgc)%metrics%z_ifc(jc,nlev_c+1,jb)) &
        prm_diag(jgc)%snowlmt(jc,jb) = -999._wp
      prm_diag(jgc)%hzerocl(jc,jb)        = z_aux3dp1_c(jc,39,jb) - 999._wp
      IF (prm_diag(jgc)%hzerocl(jc,jb) < p_nh_state(jgc)%metrics%z_ifc(jc,nlev_c+1,jb)) &
        prm_diag(jgc)%hzerocl(jc,jb) = -999._wp
      prm_diag(jgc)%clcl(jc,jb)           = MIN(1._wp,z_aux3dp1_c(jc,40,jb))
      prm_diag(jgc)%clcm(jc,jb)           = MIN(1._wp,z_aux3dp1_c(jc,41,jb))
      prm_diag(jgc)%clch(jc,jb)           = MIN(1._wp,z_aux3dp1_c(jc,42,jb))
      prm_diag(jgc)%clct(jc,jb)           = MIN(1._wp,z_aux3dp1_c(jc,43,jb))
      prm_diag(jgc)%cape(jc,jb)           = z_aux3dp1_c(jc,44,jb)
      prm_diag(jgc)%swflx_up_toa(jc,jb)   = z_aux3dp1_c(jc,45,jb)
      prm_diag(jgc)%swflx_up_sfc(jc,jb)   = z_aux3dp1_c(jc,46,jb)
      prm_diag(jgc)%tmax_2m(jc,jb)        = MAX(z_aux3dp1_c(jc,47,jb),prm_diag(jgc)%t_2m(jc,jb))
      prm_diag(jgc)%tmin_2m(jc,jb)        = MIN(z_aux3dp1_c(jc,48,jb),prm_diag(jgc)%t_2m(jc,jb),prm_diag(jgc)%tmax_2m(jc,jb))
      prm_diag(jgc)%tot_cld_vi(jc,jb,1:3) = z_aux3dp1_c(jc,49:51,jb)
      p_nh_state(jgc)%diag%tracer_vi(jc,jb,1:5) = z_aux3dp1_c(jc,52:56,jb)
      prm_diag(jgc)%clct_mod(jc,jb)       = MIN(1._wp,z_aux3dp1_c(jc,57,jb))

      IF (atm_phy_nwp_config(jgc)%inwp_gscp == 2) THEN
        prm_diag(jgc)%graupel_gsp(jc,jb)      = MAX(z_aux3dp1_c(jc,58,jb),prm_diag(jgc)%graupel_gsp(jc,jb))
        prm_diag(jgc)%graupel_gsp_rate(jc,jb) = z_aux3dp1_c(jc,59,jb)
      ENDIF
      prm_diag(jgc)%tvm(jc,jb)            = z_aux3dp1_c(jc,60,jb)
      prm_diag(jgc)%tvh(jc,jb)            = z_aux3dp1_c(jc,61,jb)
      prm_diag(jgc)%tkr(jc,jb)            = z_aux3dp1_c(jc,62,jb)
      prm_diag(jgc)%cldepth(jc,jb)        = MIN(1._wp,z_aux3dp1_c(jc,63,jb))
      prm_diag(jgc)%t_2m_land(jc,jb)      = z_aux3dp1_c(jc,64,jb)
      prm_diag(jgc)%td_2m_land(jc,jb)     = z_aux3dp1_c(jc,65,jb)
      prm_diag(jgc)%rh_2m_land(jc,jb)     = z_aux3dp1_c(jc,66,jb)

      prm_diag(jgc)%u_10m(jc,jb)          = z_aux3dp2_c(jc,1,jb)
      prm_diag(jgc)%v_10m(jc,jb)          = z_aux3dp2_c(jc,2,jb)
      prm_diag(jgc)%lhfl_s(jc,jb)         = z_aux3dp2_c(jc,3,jb)
      prm_diag(jgc)%lhfl_bs(jc,jb)        = z_aux3dp2_c(jc,4,jb)
      prm_diag(jgc)%shfl_s(jc,jb)         = z_aux3dp2_c(jc,5,jb)
      prm_diag(jgc)%qhfl_s(jc,jb)         = z_aux3dp2_c(jc,6,jb)
      prm_diag(jgc)%umfl_s(jc,jb)         = z_aux3dp2_c(jc,7,jb)
      prm_diag(jgc)%vmfl_s(jc,jb)         = z_aux3dp2_c(jc,8,jb)
      prm_diag(jgc)%alhfl_s(jc,jb)        = z_aux3dp2_c(jc,9,jb)
      prm_diag(jgc)%alhfl_bs(jc,jb)       = z_aux3dp2_c(jc,10,jb)
      prm_diag(jgc)%ashfl_s(jc,jb)        = z_aux3dp2_c(jc,11,jb)
      prm_diag(jgc)%aqhfl_s(jc,jb)        = z_aux3dp2_c(jc,12,jb)
      prm_diag(jgc)%aumfl_s(jc,jb)        = z_aux3dp2_c(jc,13,jb)
      prm_diag(jgc)%avmfl_s(jc,jb)        = z_aux3dp2_c(jc,14,jb)
      prm_diag(jgc)%lwflxsfc(jc,jb)       = z_aux3dp2_c(jc,15,jb)
      prm_diag(jgc)%lwflxsfc_a(jc,jb)     = z_aux3dp2_c(jc,16,jb)
      prm_diag(jgc)%lwflxtoa_a(jc,jb)     = z_aux3dp2_c(jc,17,jb)

      ! Special treatment for convection base and top height (no convection => - 500 m)
      prm_diag(jgc)%hbas_con(jc,jb)       = z_aux3dp2_c(jc,18,jb)
      IF (prm_diag(jgc)%hbas_con(jc,jb) < p_nh_state(jgc)%metrics%z_ifc(jc,nlev_c+1,jb)) THEN
        prm_diag(jgc)%hbas_con(jc,jb) = -500._wp
        prm_diag(jgc)%htop_con(jc,jb) = -500._wp
      ELSE
        prm_diag(jgc)%htop_con(jc,jb) = z_aux3dp1_c(jc,36,jb) + prm_diag(jgc)%hbas_con(jc,jb)
      ENDIF
    ENDDO

    IF (lsfc_interp) THEN
      DO jc = i_startidx, i_endidx
        ptr_ldiagc%t_snow(jc,jb)           = z_aux3dl2_c(jc,1,jb)
        ptr_ldiagc%t_s(jc,jb)              = z_aux3dl2_c(jc,2,jb)
        ptr_ldiagc%w_snow(jc,jb)           = z_aux3dl2_c(jc,3,jb)
        ptr_ldiagc%rho_snow(jc,jb)         = z_aux3dl2_c(jc,4,jb)
        ptr_ldiagc%w_i(jc,jb)              = z_aux3dl2_c(jc,5,jb)
        ptr_ldiagc%h_snow(jc,jb)           = z_aux3dl2_c(jc,6,jb)
        ptr_ldiagc%freshsnow(jc,jb)        = MIN(1._wp,z_aux3dl2_c(jc,7,jb))
        ptr_ldiagc%snowfrac(jc,jb)         = MIN(1._wp,z_aux3dl2_c(jc,8,jb))
        ptr_ldiagc%snowfrac_lc(jc,jb)      = MIN(1._wp,z_aux3dl2_c(jc,8,jb))
        ptr_ldiagc%runoff_s(jc,jb)         = z_aux3dl2_c(jc,9,jb)
        ptr_ldiagc%runoff_g(jc,jb)         = z_aux3dl2_c(jc,10,jb)
        ptr_ldiagc%t_so(jc,nlev_soil+1,jb) = z_aux3dl2_c(jc,11,jb)
        IF (lmulti_snow) &
          ptr_ldiagc%t_snow_mult(jc,nlev_snow+1,jb) = z_aux3dl2_c(jc,12,jb)
      ENDDO

      IF (lseaice) THEN
        DO jc = i_startidx, i_endidx
          ptr_wprogc%t_ice(jc,jb)     = MIN(tmelt,z_aux3dl2_c(jc,13,jb))
          ptr_wprogc%h_ice(jc,jb)     = MAX(0._wp,z_aux3dl2_c(jc,14,jb))
          ptr_wprogc%t_snow_si(jc,jb) = MIN(tmelt,z_aux3dl2_c(jc,15,jb))
          ptr_wprogc%h_snow_si(jc,jb) = MAX(0._wp,z_aux3dl2_c(jc,16,jb))
          ptr_ldiagc%fr_seaice(jc,jb) = MAX(0._wp,MIN(1._wp,z_aux3dl2_c(jc,17,jb)))
          IF (ptr_ldiagc%fr_seaice(jc,jb) < frsi_min )         ptr_ldiagc%fr_seaice(jc,jb) = 0._wp
          IF (ptr_ldiagc%fr_seaice(jc,jb) > (1._wp-frsi_min) ) ptr_ldiagc%fr_seaice(jc,jb) = 1._wp
          IF (ext_data(jgc)%atm%fr_land(jc,jb) >= 1._wp-MAX(frlake_thrhld,frsea_thrhld)) THEN ! pure land point
            ptr_wprogc%h_ice(jc,jb) = 0._wp
            ptr_wprogc%h_snow_si(jc,jb) = 0._wp
            ptr_ldiagc%fr_seaice(jc,jb) = 0._wp
          ENDIF
        ENDDO
      ENDIF

      IF (llake) THEN

        ! preset lake variables on nest boundary points with default values
        DO jc = i_startidx, i_endidx
          ptr_wprogc%t_snow_lk(jc,jb) = tmelt
          ptr_wprogc%h_snow_lk(jc,jb) = 0._wp
          ptr_wprogc%t_mnw_lk (jc,jb) = tmelt
          ptr_wprogc%t_wml_lk (jc,jb) = tmelt
          ptr_wprogc%t_bot_lk (jc,jb) = tmelt
          ptr_wprogc%c_t_lk   (jc,jb) = C_T_min
          ptr_wprogc%h_ml_lk  (jc,jb) = 0._wp
          ptr_wprogc%t_b1_lk  (jc,jb) = tpl_T_r
          ptr_wprogc%h_b1_lk  (jc,jb) = rflk_depth_bs_ref
        ENDDO

        ! ensure that only nest boundary points are processed
        i_count = 0
        DO ic = 1, ext_data(jgc)%atm%fp_count(jb)
          jc = ext_data(jgc)%atm%idx_lst_fp(ic,jb)
          IF (jc >= i_startidx .AND. jc <= i_endidx) THEN
            i_count = i_count + 1
            indlist(i_count) = jc
          ENDIF
        ENDDO

        CALL flake_coldinit(                                    &
             nflkgb      = i_count,                             &
             idx_lst_fp  = indlist,                             &
             depth_lk    = ext_data(jgc)%atm%depth_lk  (:,jb),  &
             tskin       = z_aux3dl2_c(:,18,jb),                &  ! estimate for lake sfc temp
             t_snow_lk_p = ptr_wprogc%t_snow_lk(:,jb),          &
             h_snow_lk_p = ptr_wprogc%h_snow_lk(:,jb),          &
             t_ice_p     = ptr_wprogc%t_ice    (:,jb),          &
             h_ice_p     = ptr_wprogc%h_ice    (:,jb),          &
             t_mnw_lk_p  = ptr_wprogc%t_mnw_lk (:,jb),          &
             t_wml_lk_p  = ptr_wprogc%t_wml_lk (:,jb),          &
             t_bot_lk_p  = ptr_wprogc%t_bot_lk (:,jb),          &
             c_t_lk_p    = ptr_wprogc%c_t_lk   (:,jb),          &
             h_ml_lk_p   = ptr_wprogc%h_ml_lk  (:,jb),          &
             t_b1_lk_p   = ptr_wprogc%t_b1_lk  (:,jb),          &
             h_b1_lk_p   = ptr_wprogc%h_b1_lk  (:,jb),          &
             t_g_lk_p    = ptr_lprogc%t_g_t    (:,jb,isub_lake) )

      ENDIF



      DO jk = 1, nlev_soil
        DO jc = i_startidx, i_endidx
          ptr_ldiagc%t_so(jc,jk,jb)     = z_aux3dso_c(jc,3*(jk-1)+1,jb)
          ptr_ldiagc%w_so(jc,jk,jb)     = z_aux3dso_c(jc,3*(jk-1)+2,jb)
          !
          ! Make sure that aggregated w_so is always larger than air dryness point
          ! at points where the soiltype allows infiltration of water.
          ! Same ad hoc fix as in mo_nwp_sfc_utils:aggregate_landvars
          ! w_so_ice is neglected
          styp = ext_data(jgc)%atm%soiltyp(jc,jb)
          IF ( (styp>=3) .AND. (styp<=8)) THEN   ! 3:sand; 8:peat
            ptr_ldiagc%w_so(jc,jk,jb) = MAX(ptr_ldiagc%w_so(jc,jk,jb),dzsoil(jk)*cadp(styp))
          ENDIF
          !
          ptr_ldiagc%w_so_ice(jc,jk,jb) = z_aux3dso_c(jc,3*(jk-1)+3,jb)
        ENDDO
      ENDDO

      ! Copy interpolated values to tile-based variables; this is actually needed in order
      ! to avoid loss of grib encoding accuracy for t_so_t
      DO jt = 1, ntiles_total
        DO jk = 1, nlev_soil
          DO jc = i_startidx, i_endidx
            ptr_lprogc%t_so_t(jc,jk,jb,jt)     = ptr_ldiagc%t_so(jc,jk,jb)
            ptr_lprogc%w_so_t(jc,jk,jb,jt)     = ptr_ldiagc%w_so(jc,jk,jb)
            ptr_lprogc%w_so_ice_t(jc,jk,jb,jt) = ptr_ldiagc%w_so_ice(jc,jk,jb)
          ENDDO
        ENDDO
      ENDDO

      IF (lmulti_snow) THEN
        DO jk = 1, nlev_snow
          DO jc = i_startidx, i_endidx
            ptr_ldiagc%t_snow_mult(jc,jk,jb)   = z_aux3dsn_c(jc,5*(jk-1)+1,jb)
            ptr_ldiagc%rho_snow_mult(jc,jk,jb) = z_aux3dsn_c(jc,5*(jk-1)+2,jb)
            ptr_ldiagc%wliq_snow(jc,jk,jb)     = z_aux3dsn_c(jc,5*(jk-1)+3,jb)
            ptr_ldiagc%wtot_snow(jc,jk,jb)     = z_aux3dsn_c(jc,5*(jk-1)+4,jb)
            ptr_ldiagc%dzh_snow(jc,jk,jb)      = z_aux3dsn_c(jc,5*(jk-1)+5,jb)
          ENDDO
        ENDDO
      ENDIF

    ENDIF



  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE interpol_phys_grf


!>
!! This routine performs lateral boundary interpolation of non-advected physics variables
!! entering into the computation of radiation. Calling this routine is required when
!! radiation is computed on a reduced grid
!!
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2011-09-19
!!
SUBROUTINE interpol_rrg_grf (jg, jgc, jn, ntl_rcf)

  ! Input grid parameters
  INTEGER, INTENT(in) :: jg, jgc, jn, ntl_rcf

  ! Pointers
  TYPE(t_patch),                POINTER :: ptr_pp
  TYPE(t_patch),                POINTER :: ptr_pc
  TYPE(t_gridref_single_state), POINTER :: ptr_grf
  TYPE(t_nwp_phy_diag),         POINTER :: prm_diagp
  TYPE(t_nwp_phy_diag),         POINTER :: prm_diagc
  TYPE(t_lnd_prog),             POINTER :: ptr_lprogp
  TYPE(t_lnd_prog),             POINTER :: ptr_lprogc_t1
  TYPE(t_lnd_prog),             POINTER :: ptr_lprogc_t2


  ! Local fields
  INTEGER, PARAMETER  :: nfields=4    ! Number of 2D fields for which boundary interpolation is needed
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, jb, jc

  ! Temporary storage to do boundary interpolation for all 2D fields in one step
  REAL(wp) :: z_aux3d_p(nproma,nfields,p_patch(jg)%nblks_c), &
              z_aux3d_c(nproma,nfields,p_patch(jgc)%nblks_c)

  ! set pointers
  ptr_pp        => p_patch(jg)
  ptr_pc        => p_patch(jgc)
  ptr_grf       => p_grf_state(jg)%p_dom(jn)
  prm_diagp     => prm_diag(jg)
  prm_diagc     => prm_diag(jgc)
  ptr_lprogp    => p_lnd_state(jg)%prog_lnd(ntl_rcf)
  ! We just need to set both time levels of the child land state,
  ! without having to care which one is now and new
  ptr_lprogc_t1 => p_lnd_state(jgc)%prog_lnd(1)
  ptr_lprogc_t2 => p_lnd_state(jgc)%prog_lnd(2)

  IF (p_test_run) THEN
     z_aux3d_p(:,:,:) = 0._wp
  ENDIF

  i_startblk = ptr_pp%cells%start_blk(1,1)
  i_endblk   = ptr_pp%nblks_c

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, 1)

    DO jc = i_startidx, i_endidx

      z_aux3d_p(jc,1,jb) = ptr_lprogp%t_g(jc,jb)
      z_aux3d_p(jc,2,jb) = prm_diagp%albdif(jc,jb)
      z_aux3d_p(jc,3,jb) = prm_diagp%albvisdif(jc,jb)
      z_aux3d_p(jc,4,jb) = prm_diagp%albnirdif(jc,jb)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ! Halo update is needed before interpolation
    CALL sync_patch_array(SYNC_C,ptr_pp,z_aux3d_p)

    CALL interpol_scal_grf (ptr_pp, ptr_pc, ptr_grf, 1, z_aux3d_p, z_aux3d_c, llimit_nneg=(/.TRUE./),&
      &                     lnoshift=.TRUE.)


  i_startblk = ptr_pc%cells%start_blk(1,1)
  i_endblk   = ptr_pc%cells%end_blk(grf_bdywidth_c,1)

  ! Note: prognostic land fields are set on both time levels to safely avoid
  ! errors when radiation calls for parent and child grids are not properly
  ! synchronized
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_pc, jb, i_startblk, i_endblk,        &
                       i_startidx, i_endidx, 1, grf_bdywidth_c)

    DO jc = i_startidx, i_endidx

      ptr_lprogc_t1%t_g(jc,jb)     = z_aux3d_c(jc,1,jb)
      ptr_lprogc_t2%t_g(jc,jb)     = z_aux3d_c(jc,1,jb)
      prm_diagc%albdif(jc,jb)      = z_aux3d_c(jc,2,jb)
      prm_diagc%albvisdif(jc,jb)   = z_aux3d_c(jc,3,jb)
      prm_diagc%albnirdif(jc,jb)   = z_aux3d_c(jc,4,jb)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE interpol_rrg_grf

!>
!! This routine copies additional model levels to the local parent grid if vertical nesting
!! is combined with a reduced radiation grid and the option latm_above_top = .TRUE.
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2015-01-26
!!
SUBROUTINE copy_rrg_ubc (jg, jgc)

  ! Input grid parameters
  INTEGER, INTENT(in) :: jg, jgc

  ! Local fields

  INTEGER :: jks, jke, nshift

  jks = MAX(0, p_patch(jgc)%nshift - nexlevs_rrg_vnest) + 1
  jke = p_patch(jgc)%nshift
  nshift = MIN(nexlevs_rrg_vnest, p_patch(jgc)%nshift)

  IF (nshift > 0) THEN
    CALL exchange_data_mult(p_patch_local_parent(jgc)%comm_pat_glb_to_loc_c, 3, 3*nshift,                        &
       RECV1=prm_diag(jgc)%buffer_rrg(:,         1:  nshift,:), SEND1=p_nh_state(jg)%diag%pres_ifc(:,jks:jke,:), &
       RECV2=prm_diag(jgc)%buffer_rrg(:,  nshift+1:2*nshift,:), SEND2=p_nh_state(jg)%diag%pres(:,jks:jke,:),     &
       RECV3=prm_diag(jgc)%buffer_rrg(:,2*nshift+1:3*nshift,:), SEND3=p_nh_state(jg)%diag%temp(:,jks:jke,:)      )
  ENDIF

END SUBROUTINE copy_rrg_ubc



!>
!! This routine performs the feedback of diagnostic physics fields for output
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2010-12-03
!!
SUBROUTINE feedback_phys_diag(jg, jgp)

  INTEGER, INTENT(IN) :: jg   ! child grid level
  INTEGER, INTENT(IN) :: jgp  ! parent grid level


  ! Pointers to types needed to minimize code duplication for MPI/no-MPI cases
  TYPE(t_grid_cells), POINTER     :: p_gcp
  TYPE(t_gridref_state), POINTER  :: p_grf
  TYPE(t_patch),      POINTER     :: p_pp

  ! Indices
  INTEGER :: jb, jc, i_chidx, i_nchdom, &
             i_startblk, i_endblk, i_startidx, i_endidx, nblks_c_lp

  INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
  REAL(wp), POINTER :: p_fbkwgt(:,:,:), p_aux3d(:,:,:)
  REAL(wp), ALLOCATABLE, TARGET :: z_aux3d_lp(:,:,:), z_aux3d_par(:,:,:)

  !-----------------------------------------------------------------------

  IF (msg_level >= 10) THEN
    WRITE(message_text,'(a,i2,a,i2)') 'Feedback of diagnostic physics fields',&
      p_patch(jg)%id,' =>',p_patch(jgp)%id
    CALL message('feedback_phys_diag',message_text)
  ENDIF

  p_grf => p_grf_state_local_parent(jg)
  p_gcp => p_patch_local_parent(jg)%cells
  p_pp  => p_patch_local_parent(jg)

  i_chidx  = p_patch(jg)%parent_child_index
  i_nchdom = MAX(1,p_patch(jgp)%n_childdom)

  ! Set pointers to index and coefficient fields for cell-based variables
  iidx => p_gcp%child_idx
  iblk => p_gcp%child_blk

  p_fbkwgt => p_grf%fbk_wgt_aw

  ! Allocation of local storage fields
  nblks_c_lp = p_gcp%end_blk(min_rlcell,i_chidx)
  ALLOCATE(z_aux3d_lp(nproma,5,nblks_c_lp), z_aux3d_par(nproma,5,p_patch(jgp)%nblks_c))

  p_aux3d => z_aux3d_lp

  ! Now average input fields to parent grid cells

  ! Start/End block in the parent domain
  i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
  i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,                      &
                       i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int)

    DO jc = i_startidx, i_endidx

      p_aux3d(jc,1,jb) =                                         &
        prm_diag(jg)%tot_prec(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        prm_diag(jg)%tot_prec(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        prm_diag(jg)%tot_prec(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        prm_diag(jg)%tot_prec(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_aux3d(jc,2,jb) =                                         &
        prm_diag(jg)%rain_con(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        prm_diag(jg)%rain_con(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        prm_diag(jg)%rain_con(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        prm_diag(jg)%rain_con(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_aux3d(jc,3,jb) =                                         &
        prm_diag(jg)%snow_con(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        prm_diag(jg)%snow_con(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        prm_diag(jg)%snow_con(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        prm_diag(jg)%snow_con(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_aux3d(jc,4,jb) =                                         &
        prm_diag(jg)%rain_gsp(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        prm_diag(jg)%rain_gsp(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        prm_diag(jg)%rain_gsp(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        prm_diag(jg)%rain_gsp(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_aux3d(jc,5,jb) =                                         &
        prm_diag(jg)%snow_gsp(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        prm_diag(jg)%snow_gsp(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        prm_diag(jg)%snow_gsp(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        prm_diag(jg)%snow_gsp(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

    ENDDO

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  CALL exchange_data(p_pp%comm_pat_loc_to_glb_c_fbk, RECV=z_aux3d_par, SEND=z_aux3d_lp)
  p_aux3d => z_aux3d_par

  i_startblk = p_patch(jgp)%cells%start_blk(1,1)
  i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell_int, i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch(jgp), jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, 1, min_rlcell_int)

    DO jc = i_startidx, i_endidx

      IF (p_grf_state(jgp)%mask_ovlp_c(jc,jb,i_chidx)) THEN
        prm_diag(jgp)%tot_prec(jc,jb)      = p_aux3d(jc,1,jb)
        prm_diag(jgp)%rain_con(jc,jb)      = p_aux3d(jc,2,jb)
        prm_diag(jgp)%snow_con(jc,jb)      = p_aux3d(jc,3,jb)
        prm_diag(jgp)%rain_gsp(jc,jb)      = p_aux3d(jc,4,jb)
        prm_diag(jgp)%snow_gsp(jc,jb)      = p_aux3d(jc,5,jb)
      END IF

    ENDDO

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  DEALLOCATE(z_aux3d_lp, z_aux3d_par)

END SUBROUTINE feedback_phys_diag


END MODULE mo_phys_nest_utilities



