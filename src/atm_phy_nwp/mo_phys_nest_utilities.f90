!>
!!  This module contains (not yet) utility programs for boundary interpolation and feedback
!!  of diagnostic variables and the upscaling and downscaling routines needed 
!!  for the reduced physics grid
!!
!! @par Revision History
!!  Developed and tested by Guenther Zaengl, DWD (2010-02-10)
!!
!! @par Copyright
!! 2002-2006 by DWD and MPI-M
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_phys_nest_utilities
!
!
USE mo_kind,                ONLY: wp
USE mo_exception,           ONLY: message_text, message
USE mo_model_domain,        ONLY: t_patch, t_grid_cells, p_patch_local_parent, p_patch
USE mo_grid_config,         ONLY: n_dom, n_dom_start
USE mo_intp_data_strc,      ONLY: t_int_state, p_int_state, p_int_state_local_parent
USE mo_grf_intp_data_strc,  ONLY: t_gridref_state, t_gridref_single_state, &
                                  p_grf_state, p_grf_state_local_parent
USE mo_nwp_phy_state,       ONLY: t_nwp_phy_diag
USE mo_nwp_lnd_types,       ONLY: t_lnd_prog, t_lnd_diag
USE mo_nwp_lnd_state,       ONLY: p_lnd_state
USE mo_grf_bdyintp,         ONLY: interpol_scal_grf
USE mo_grf_nudgintp,        ONLY: interpol_scal_nudging
USE mo_parallel_config,     ONLY: nproma, p_test_run
USE mo_run_config,          ONLY: msg_level, iqv, iqc, iqi
USE mo_nwp_phy_state,       ONLY: prm_diag
USE mo_impl_constants,      ONLY: min_rlcell, min_rlcell_int
USE mo_physical_constants,  ONLY: rd, grav, stbo, vtmpc1
USE mo_satad,               ONLY: qsat_rho
USE mo_loopindices,         ONLY: get_indices_c
USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_ovlparea_start_c, grf_fbk_start_c
USE mo_vertical_coord_table,ONLY: vct_a
USE mo_mpi,                 ONLY: my_process_is_mpi_seq
USE mo_communication,       ONLY: exchange_data, exchange_data_mult
USE mo_sync,                ONLY: SYNC_C, sync_patch_array, sync_patch_array_mult
USE mo_lnd_nwp_config,      ONLY: nlev_soil, nlev_snow, lmulti_snow

IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

PUBLIC :: upscale_rad_input, downscale_rad_output, interpol_phys_grf, feedback_phys_diag, &
  &       upscale_rad_input_rg, downscale_rad_output_rg, interpol_rrg_grf

CONTAINS

!>
!! This routine averages the input fields for RRTM radiation to the next coarser grid level.
!!
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2010-12-01
!!
SUBROUTINE upscale_rad_input(jg, jgp, nlev_rg, fr_land, fr_glac, emis_rad, &
  cosmu0, albvisdir, albnirdir, albvisdif, albnirdif,                      &
  tsfc, pres_ifc, pres, temp, acdnc, tot_cld, q_o3,                        &
  aeq1, aeq2, aeq3, aeq4, aeq5,                                            &
  rg_fr_land, rg_fr_glac, rg_emis_rad,                                     &
  rg_cosmu0, rg_albvisdir, rg_albnirdir, rg_albvisdif, rg_albnirdif,       &
  rg_tsfc, rg_pres_ifc, rg_pres, rg_temp, rg_acdnc, rg_tot_cld, rg_q_o3,   &
  rg_aeq1, rg_aeq2, rg_aeq3, rg_aeq4, rg_aeq5 )

  ! Input grid parameters
  INTEGER, INTENT(IN)  :: jg, jgp  ! domain IDs of main and reduced grids
  INTEGER, INTENT(IN)  :: nlev_rg  ! number of model levels on reduced grid

  ! Input fields (on full grid)
  REAL(wp), INTENT(IN) ::                                                                 &
    fr_land(:,:), fr_glac(:,:), emis_rad(:,:),                                            &
    cosmu0(:,:), albvisdir(:,:), albnirdir(:,:), albvisdif(:,:), albnirdif(:,:),          &
    tsfc(:,:), pres_ifc(:,:,:), pres(:,:,:), temp(:,:,:), acdnc(:,:,:), tot_cld(:,:,:,:), &
    q_o3(:,:,:), aeq1(:,:,:), aeq2(:,:,:), aeq3(:,:,:), aeq4(:,:,:), aeq5(:,:,:)

  ! Corresponding output fields (on reduced grid)
  REAL(wp), TARGET, INTENT(OUT) ::                                           &
    rg_fr_land(:,:),rg_fr_glac(:,:), rg_emis_rad(:,:),                       &
    rg_cosmu0(:,:), rg_albvisdir(:,:), rg_albnirdir(:,:), rg_albvisdif(:,:), &
    rg_albnirdif(:,:), rg_tsfc(:,:), rg_pres_ifc(:,:,:), rg_pres(:,:,:),     &
    rg_temp(:,:,:), rg_acdnc(:,:,:), rg_tot_cld(:,:,:,:), rg_q_o3(:,:,:),    &
    rg_aeq1(:,:,:),rg_aeq2(:,:,:),rg_aeq3(:,:,:),rg_aeq4(:,:,:),rg_aeq5(:,:,:)

  ! Intermediate storage fields needed in the case of MPI parallelization
  REAL(wp), ALLOCATABLE, TARGET ::                                         &
    z_fr_land(:,:),z_fr_glac(:,:), z_emis_rad(:,:),                        &
    z_cosmu0(:,:), z_albvisdir(:,:), z_albnirdir(:,:), z_albvisdif(:,:),   &
    z_albnirdif(:,:), z_tsfc(:,:), z_pres_ifc(:,:,:), z_pres(:,:,:),       &
    z_temp(:,:,:), z_acdnc(:,:,:), z_tot_cld(:,:,:,:), z_q_o3(:,:,:),      &
    z_aeq1(:,:,:),z_aeq2(:,:,:),z_aeq3(:,:,:),z_aeq4(:,:,:),z_aeq5(:,:,:), &
    z_aux3d(:,:,:), zrg_aux3d(:,:,:)

  ! Pointers to output fields (no MPI) or intermediate fields (MPI)
  REAL(wp), POINTER ::                                                   &
    p_fr_land(:,:),p_fr_glac(:,:), p_emis_rad(:,:),                      &
    p_cosmu0(:,:), p_albvisdir(:,:), p_albnirdir(:,:), p_albvisdif(:,:), &
    p_albnirdif(:,:), p_tsfc(:,:), p_pres_ifc(:,:,:), p_pres(:,:,:),     &
    p_temp(:,:,:), p_acdnc(:,:,:), p_tot_cld(:,:,:,:), p_q_o3(:,:,:),    &
    p_aeq1(:,:,:),p_aeq2(:,:,:),p_aeq3(:,:,:),p_aeq4(:,:,:),p_aeq5(:,:,:)


  ! Pointers to types needed to minimize code duplication for MPI/no-MPI cases
  TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
  TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
  TYPE(t_patch),      POINTER     :: p_pp => NULL()

  ! Indices
  INTEGER :: jb, jc, jk, jk1, i_chidx, i_nchdom, &
             i_startblk, i_endblk, i_startidx, i_endidx, nblks_c_lp

  INTEGER :: nlev, nlevp1      !< number of full and half levels
  INTEGER :: nshift, nlevp1_rg, nst
  REAL(wp) :: scalfac, exdist, rdelta_z

  INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
  LOGICAL :: l_parallel
  REAL(wp), POINTER :: p_fbkwgt(:,:,:)

  !-----------------------------------------------------------------------

  IF (msg_level >= 10) THEN
    WRITE(message_text,'(a,i2,a,i2)') 'Upscaling of radiation input fields',&
      jg,' =>',jgp
    CALL message('upscale_rad_input',message_text)
  ENDIF

  IF (my_process_is_mpi_seq()) THEN
    l_parallel = .FALSE.
  ELSE
    l_parallel = .TRUE.
  ENDIF

  ! Number of levels of the full grid
  nlev   = p_patch(jg)%nlev
  nlevp1 = p_patch(jg)%nlevp1

  ! nlev_rg carries the number of model levels of the reduced grid,
  ! which may be larger than nlev
  nlevp1_rg = nlev_rg + 1
  nshift = nlev_rg - nlev ! resulting shift parameter

  IF (l_parallel) THEN
    p_grf => p_grf_state_local_parent(jg)
    p_gcp => p_patch_local_parent(jg)%cells
    p_pp  => p_patch_local_parent(jg)
  ELSE
    p_grf => p_grf_state(jgp)
    p_gcp => p_patch(jgp)%cells
    p_pp  => p_patch(jgp)
  ENDIF

  i_chidx  = p_patch(jg)%parent_child_index
  i_nchdom = MAX(1,p_patch(jgp)%n_childdom)

  ! Set pointers to index and coefficient fields for cell-based variables
  iidx => p_gcp%child_idx
  iblk => p_gcp%child_blk

  p_fbkwgt => p_grf%fbk_wgt_c

  ! layer shift w.r.t. global grid (> 0 in case of vertical nesting)
  nst = p_patch(jg)%nshift_total
  ! inverse height difference between layer 1 and 2
  rdelta_z = 1._wp/(0.5_wp*(vct_a(nst+1)-vct_a(nst+3))) ! note: vct refers to half levels

  ! scale factor for extrapolation to the barycenter of the passive top layer
  ! (gives multiplied with T an estimate for the height distance between the model top
  !  and the mass center of the passive layer)
  scalfac = LOG(2._wp)*rd/grav

  ! height distance between uppermost full and half levels
  exdist = 0.5_wp*(vct_a(nst+1)-vct_a(nst+2))

  ! Allocation of local storage fields at local parent level in MPI-case
  IF (l_parallel .AND. jgp == 0) THEN
    nblks_c_lp = p_gcp%end_blk(min_rlcell,i_chidx)

    ALLOCATE(z_fr_land(nproma,nblks_c_lp), z_fr_glac(nproma,nblks_c_lp),                &
             z_emis_rad(nproma,nblks_c_lp),                                             &
             z_cosmu0(nproma,nblks_c_lp), z_albvisdir(nproma,nblks_c_lp),               &
             z_albnirdir(nproma,nblks_c_lp), z_albvisdif(nproma,nblks_c_lp),            &
             z_albnirdif(nproma,nblks_c_lp), z_tsfc(nproma,nblks_c_lp),                 &
             z_pres_ifc(nproma,nlevp1_rg,nblks_c_lp), z_pres(nproma,nlev_rg,nblks_c_lp),&
             z_temp(nproma,nlev_rg,nblks_c_lp), z_acdnc(nproma,nlev_rg,nblks_c_lp),     &
             z_tot_cld(nproma,nlev_rg,nblks_c_lp,4), z_q_o3(nproma,nlev_rg,nblks_c_lp), &
             z_aeq1(nproma,nlev_rg,nblks_c_lp), z_aeq2(nproma,nlev_rg,nblks_c_lp),      &
             z_aeq3(nproma,nlev_rg,nblks_c_lp), z_aeq4(nproma,nlev_rg,nblks_c_lp),      &
             z_aeq5(nproma,nlev_rg,nblks_c_lp),                                         &
             z_aux3d(nproma,9,nblks_c_lp), zrg_aux3d(nproma,9,p_patch(jgp)%nblks_c) )

  ENDIF

  ! Set pointers to either the parent-level variables (non-MPI case) or to the
  ! intermediate storage fields (MPI case)
  IF (l_parallel .AND. jgp == 0) THEN
    p_fr_land    => z_fr_land
    p_fr_glac    => z_fr_glac
    p_emis_rad   => z_emis_rad
    p_cosmu0     => z_cosmu0
    p_albvisdir  => z_albvisdir
    p_albnirdir  => z_albnirdir
    p_albvisdif  => z_albvisdif
    p_albnirdif  => z_albnirdif
    p_tsfc       => z_tsfc
    p_pres_ifc   => z_pres_ifc
    p_pres       => z_pres
    p_temp       => z_temp
    p_acdnc      => z_acdnc
    p_tot_cld    => z_tot_cld
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
    p_tsfc       => rg_tsfc
    p_pres_ifc   => rg_pres_ifc
    p_pres       => rg_pres
    p_temp       => rg_temp
    p_acdnc      => rg_acdnc
    p_tot_cld    => rg_tot_cld
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
    p_tsfc       = 0._wp
    p_pres_ifc   = 0._wp
    p_pres       = 0._wp
    p_temp       = 0._wp
    p_acdnc      = 0._wp
    p_tot_cld    = 0._wp
    p_q_o3       = 0._wp
    p_aeq1       = 0._wp
    p_aeq2       = 0._wp
    p_aeq3       = 0._wp
    p_aeq4       = 0._wp
    p_aeq5       = 0._wp
  ENDIF


  ! Now average input fields to parent grid cells

  ! Start/End block in the parent domain
  i_startblk = p_gcp%start_blk(grf_ovlparea_start_c,i_chidx)
  i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jk1) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,                                &
                       i_startidx, i_endidx, grf_ovlparea_start_c, min_rlcell_int, i_chidx)

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

      p_tsfc(jc,jb) =                                         &
        tsfc(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        tsfc(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        tsfc(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        tsfc(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_pres_ifc(jc,nlevp1_rg,jb) =                                      &
        pres_ifc(iidx(jc,jb,1),nlevp1,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        pres_ifc(iidx(jc,jb,2),nlevp1,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        pres_ifc(iidx(jc,jb,3),nlevp1,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        pres_ifc(iidx(jc,jb,4),nlevp1,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

    ENDDO

    IF (l_parallel .AND. jgp == 0) THEN ! combine 2D fields in a 3D field to speed up MPI communication
      DO jc = i_startidx, i_endidx
        z_aux3d(jc,1,jb) = p_cosmu0(jc,jb)
        z_aux3d(jc,2,jb) = p_albvisdir(jc,jb)
        z_aux3d(jc,3,jb) = p_albnirdir(jc,jb)
        z_aux3d(jc,4,jb) = p_albvisdif(jc,jb)
        z_aux3d(jc,5,jb) = p_albnirdif(jc,jb)
        z_aux3d(jc,6,jb) = p_tsfc(jc,jb)
        z_aux3d(jc,7,jb) = p_fr_land(jc,jb)
        z_aux3d(jc,8,jb) = p_fr_glac(jc,jb)
        z_aux3d(jc,9,jb) = p_emis_rad(jc,jb)        
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
        
        p_acdnc(jc,jk1,jb) =                                        &
          acdnc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          acdnc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          acdnc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          acdnc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

!CDIR EXPAND=4
        p_tot_cld(jc,jk1,jb,1:4) =                                        &
          tot_cld(iidx(jc,jb,1),jk,iblk(jc,jb,1),1:4)*p_fbkwgt(jc,jb,1) + &
          tot_cld(iidx(jc,jb,2),jk,iblk(jc,jb,2),1:4)*p_fbkwgt(jc,jb,2) + &
          tot_cld(iidx(jc,jb,3),jk,iblk(jc,jb,3),1:4)*p_fbkwgt(jc,jb,3) + &
          tot_cld(iidx(jc,jb,4),jk,iblk(jc,jb,4),1:4)*p_fbkwgt(jc,jb,4)

      ENDDO
    ENDDO

    IF (nshift == 1) THEN ! set values for passive top layer if present
      DO jc = i_startidx, i_endidx
        p_pres_ifc(jc,1,jb) = 0._wp ! TOA
        p_pres(jc,1,jb) = 0.5_wp*p_pres_ifc(jc,2,jb)
        ! Temperature is linearly extrapolated to the barycenter of the passive layer
        p_temp(jc,1,jb) = p_temp(jc,2,jb) + (scalfac*p_temp(jc,2,jb)+exdist)* &
          MIN(4.e-3_wp,MAX(-6.5e-3_wp,rdelta_z*(p_temp(jc,2,jb)-p_temp(jc,3,jb))))
        ! For ozone, aerosols and cloud fields, a no-gradient condition is assumed
        p_q_o3(jc,1,jb) = p_q_o3(jc,2,jb)
        p_aeq1(jc,1,jb) = p_aeq1(jc,2,jb)
        p_aeq2(jc,1,jb) = p_aeq2(jc,2,jb)
        p_aeq3(jc,1,jb) = p_aeq3(jc,2,jb)
        p_aeq4(jc,1,jb) = p_aeq4(jc,2,jb)
        p_aeq5(jc,1,jb) = p_aeq5(jc,2,jb)
        p_acdnc(jc,1,jb) = p_acdnc(jc,2,jb)
!CDIR EXPAND=4
        p_tot_cld(jc,1,jb,1:4) = p_tot_cld(jc,2,jb,1:4)
      ENDDO
    ENDIF

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  IF (l_parallel .AND. jgp == 0) THEN

    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 6, 5*nlev_rg+10, &
                            RECV1=rg_pres_ifc, SEND1=z_pres_ifc,             &
                            RECV2=rg_pres,     SEND2=z_pres,                 &
                            RECV3=rg_temp,     SEND3=z_temp,                 &
                            RECV4=rg_acdnc,    SEND4=z_acdnc,                &
                            RECV5=zrg_aux3d,   SEND5=z_aux3d,                &
                            RECV6=rg_q_o3,     SEND6=z_q_o3                  )

    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 5, 5*nlev_rg, &
                            RECV1=rg_aeq1,     SEND1=z_aeq1,              &
                            RECV2=rg_aeq2,     SEND2=z_aeq2,              &
                            RECV3=rg_aeq3,     SEND3=z_aeq3,              &
                            RECV4=rg_aeq4,     SEND4=z_aeq4,              &
                            RECV5=rg_aeq5,     SEND5=z_aeq5               )
    
    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 4, 4*nlev_rg, &
                            RECV4D=rg_tot_cld, SEND4D=z_tot_cld           )


    i_startblk = p_patch(jgp)%cells%start_blk(1,1)
    i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell,i_nchdom)

    ! OpenMP section commented because the DO loop does almost no work (overhead larger than benefit)
!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch(jgp), jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rlcell, i_nchdom)

      DO jc = i_startidx, i_endidx
        rg_cosmu0(jc,jb)    = zrg_aux3d(jc,1,jb)
        rg_albvisdir(jc,jb) = zrg_aux3d(jc,2,jb)
        rg_albnirdir(jc,jb) = zrg_aux3d(jc,3,jb)
        rg_albvisdif(jc,jb) = zrg_aux3d(jc,4,jb)
        rg_albnirdif(jc,jb) = zrg_aux3d(jc,5,jb)
        rg_tsfc(jc,jb)      = zrg_aux3d(jc,6,jb)
        rg_fr_land(jc,jb)   = zrg_aux3d(jc,7,jb)
        rg_fr_glac(jc,jb)   = zrg_aux3d(jc,8,jb)
        rg_emis_rad(jc,jb)  = zrg_aux3d(jc,9,jb)
      ENDDO

    ENDDO
!!$OMP END DO NOWAIT
!!$OMP END PARALLEL

    DEALLOCATE(z_fr_land, z_fr_glac, z_emis_rad, z_cosmu0, z_albvisdir, z_albnirdir, z_albvisdif,&
      & z_albnirdif, z_tsfc, z_pres_ifc, z_pres, z_temp, z_acdnc, z_tot_cld, z_q_o3,   &
      & z_aeq1, z_aeq2, z_aeq3, z_aeq4, z_aeq5, z_aux3d, zrg_aux3d )

  ENDIF

END SUBROUTINE upscale_rad_input

!>
!! This routine interpolates the output fields of RRTM radiation from the reduced 
!! grid to the full grid.
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2010-12-03
!!
SUBROUTINE downscale_rad_output(jg, jgp, nlev_rg,                         &
  rg_aclcov, rg_lwflxclr, rg_lwflxall, rg_trsolclr, rg_trsolall,          &
  tsfc_rg, albvisdif_rg, emis_rad_rg, cosmu0_rg, tot_cld_rg, pres_ifc_rg, &
  tsfc, albvisdif, aclcov, lwflxclr, lwflxall, trsolclr, trsolall         )


  ! Input grid parameters
  INTEGER, INTENT(IN)  :: jg, jgp  ! domain IDs of main and reduced grids
  INTEGER, INTENT(IN)  :: nlev_rg  ! number of model levels on reduced grid

  ! Input fields (on reduced grid) to be downscaled to full grid
  REAL(wp), TARGET, INTENT(IN) ::                                               &
    rg_aclcov(:,:), rg_lwflxclr(:,:,:), rg_lwflxall(:,:,:), rg_trsolclr(:,:,:), &
    rg_trsolall(:,:,:)

  ! Auxiliary input fields on reduced grid needed for downscaling corrections
  REAL(wp), INTENT(IN), TARGET ::                                      &
    tsfc_rg(:,:), albvisdif_rg(:,:), emis_rad_rg(:,:), cosmu0_rg(:,:), &
    tot_cld_rg(:,:,:,:), pres_ifc_rg(:,:,:)

  ! Auxiliary input fields on full grid needed for downscaling corrections
  REAL(wp), INTENT(IN) :: tsfc(:,:), albvisdif(:,:)

  ! Downscaled output fields (on full grid)
  REAL(wp), INTENT(OUT) ::                                                        &
    aclcov(:,:), lwflxclr(:,:,:), lwflxall(:,:,:), trsolclr(:,:,:), trsolall(:,:,:)

  ! Intermediate storage fields needed in the case of MPI parallelization
  REAL(wp), ALLOCATABLE, TARGET ::                                              &
    z_lwflxclr(:,:,:), z_lwflxall(:,:,:), z_trsolclr(:,:,:), z_trsolall(:,:,:), &
    z_pres_ifc(:,:,:), z_tot_cld(:,:,:,:)

  ! Pointers to output fields (no MPI) or intermediate fields (MPI)
  REAL(wp), POINTER ::                                                          &
    p_lwflxclr(:,:,:), p_lwflxall(:,:,:), p_trsolclr(:,:,:), p_trsolall(:,:,:), &
    p_pres_ifc(:,:,:), p_tot_cld(:,:,:,:)

  ! Additional storage fields to map 2D array(s) to 3D array
  REAL(wp), ALLOCATABLE :: zpg_aux3d(:,:,:), zrg_aux3d(:,:,:), z_aux3d(:,:,:)

  ! Auxiliary fields on full grid for back-interpolated values
  REAL(wp), DIMENSION(nproma,p_patch(jg)%nblks_c) :: tsfc_backintp, alb_backintp

  ! More local variables
  REAL(wp) :: pscal, dpresg, pfaclw, intqctot

  REAL(wp), DIMENSION(nproma) :: tqv, dlwem_o_dtg, swfac1, swfac2, lwfac1, lwfac2

  REAL(wp), DIMENSION(nproma,p_patch(jg)%nlevp1) :: intclw, intcli, dtrans_o_dalb_clr, &
    dtrans_o_dalb_all, dlwflxclr_o_dtg, dlwflxall_o_dtg, pfacswc, pfacswa

  ! Pointers to types needed to minimize code duplication for MPI/no-MPI cases
  TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
  TYPE(t_int_state),  POINTER     :: p_int => NULL()
  TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
  TYPE(t_patch),      POINTER     :: p_pp => NULL()

  ! Indices
  INTEGER :: i_chidx, i_nchdom, nblks_c_lp
  INTEGER :: jb, jk, jk1, jc, i_startblk, i_endblk, i_startidx, i_endidx
  INTEGER :: jc1, jc2, jc3, jc4, jb1, jb2, jb3, jb4

  INTEGER :: nlev, nlevp1      !< number of full and half levels
  INTEGER :: nshift, nlevp1_rg, nlev_tot
  INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk

  INTEGER :: n2dvars, iclcov, itsfc, ialb, iemis, icosmu0

  LOGICAL :: l_parallel, l_limit(5)
  REAL(wp) :: rlimval(5)
!-----------------------------------------------------------------------

  IF (msg_level >= 10) THEN
    WRITE(message_text,'(a,i2,a,i2)') 'Downscaling of radiation output fields',&
      jgp,' =>',jg
    CALL message('downscale_rad_output',message_text)
  ENDIF

  IF (my_process_is_mpi_seq()) THEN
    l_parallel = .FALSE.
  ELSE
    l_parallel = .TRUE.
  ENDIF

  ! For the time being, the radiation grid is assumed to have the same levels
  ! as the full grid
  nlev   = p_patch(jg)%nlev
  nlevp1 = p_patch(jg)%nlevp1

  ! nlev_rg carries the number of model levels of the reduced grid,
  ! which may be larger than nlev
  nlevp1_rg = nlev_rg + 1
  nshift = nlev_rg - nlev ! resulting shift parameter

  IF (l_parallel) THEN
    p_grf => p_grf_state_local_parent(jg)
    p_int => p_int_state_local_parent(jg)
    p_gcp => p_patch_local_parent(jg)%cells
    p_pp  => p_patch_local_parent(jg)
  ELSE
    p_grf => p_grf_state(jgp)
    p_int => p_int_state(jgp)
    p_gcp => p_patch(jgp)%cells
    p_pp  => p_patch(jgp)
  ENDIF

  i_nchdom = MAX(1,p_patch(jg)%n_childdom)
  i_chidx  = p_patch(jg)%parent_child_index

  nblks_c_lp = p_pp%nblks_c

  ! pointers to child index/block
  iidx => p_gcp%child_idx
  iblk => p_gcp%child_blk

  ! named constants for accessing 2D variables contained in zrg_aux3d
  n2dvars = nshift+5
  iclcov  = nshift+1
  itsfc   = nshift+2
  ialb    = nshift+3
  iemis   = nshift+4
  icosmu0 = nshift+5


  ! Allocation of local storage fields at local parent level in MPI-case
  IF (l_parallel .AND. jgp == 0) THEN

    ALLOCATE( z_lwflxclr(nproma,nlevp1_rg,nblks_c_lp), z_lwflxall(nproma,nlevp1_rg,nblks_c_lp), &
              z_trsolclr(nproma,nlevp1_rg,nblks_c_lp), z_trsolall(nproma,nlevp1_rg,nblks_c_lp), &
              z_pres_ifc(nproma,nlevp1_rg,nblks_c_lp), z_tot_cld(nproma,nlev_rg,nblks_c_lp,4),  &
              zpg_aux3d(nproma,n2dvars,p_patch(jgp)%nblks_c),                                   &
              zrg_aux3d(nproma,n2dvars,nblks_c_lp),   z_aux3d(nproma,5,p_patch(jg)%nblks_c)     )

  ELSE

    ALLOCATE(zrg_aux3d(nproma,n2dvars,nblks_c_lp),z_aux3d(nproma,5,p_patch(jg)%nblks_c))

  ENDIF

  ! Perform communication from parent to local parent grid in the MPI case,
  ! and set pointers such that further processing is the same for MPI / non-MPI cases
  IF (l_parallel .AND. jgp == 0) THEN

    IF (nshift > 0) zpg_aux3d(:,1:nshift,:) = 0._wp
    zpg_aux3d(:,iclcov,:) = rg_aclcov(:,:)
    zpg_aux3d(:,itsfc,:)  = tsfc_rg(:,:)
    zpg_aux3d(:,ialb,:)   = albvisdif_rg(:,:)
    zpg_aux3d(:,iemis,:)  = emis_rad_rg(:,:)
    zpg_aux3d(:,icosmu0,:)= cosmu0_rg(:,:)

    nlev_tot = 5*nlevp1_rg + n2dvars

    CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, 6, nlev_tot, &
                            RECV1=z_lwflxclr, SEND1=rg_lwflxclr,     &
                            RECV2=z_lwflxall, SEND2=rg_lwflxall,     &
                            RECV3=z_trsolclr, SEND3=rg_trsolclr,     &
                            RECV4=z_trsolall, SEND4=rg_trsolall,     &
                            RECV5=z_pres_ifc, SEND5=pres_ifc_rg,     &
                            RECV6=zrg_aux3d , SEND6=zpg_aux3d        )

    CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, 4, 4*nlev_rg, &
                            RECV4D=z_tot_cld, SEND4D=tot_cld_rg       )

    p_lwflxclr   => z_lwflxclr
    p_lwflxall   => z_lwflxall
    p_trsolclr   => z_trsolclr
    p_trsolall   => z_trsolall
    p_pres_ifc   => z_pres_ifc
    p_tot_cld    => z_tot_cld
  ELSE
    p_lwflxclr   => rg_lwflxclr
    p_lwflxall   => rg_lwflxall
    p_trsolclr   => rg_trsolclr
    p_trsolall   => rg_trsolall
    p_pres_ifc   => pres_ifc_rg
    p_tot_cld    => tot_cld_rg

    IF (nshift > 0) zrg_aux3d(:,1:nshift,:) = 0._wp
    zrg_aux3d(:,iclcov,:) = rg_aclcov(:,:)
    zrg_aux3d(:,itsfc,:)  = tsfc_rg(:,:)
    zrg_aux3d(:,ialb,:)   = albvisdif_rg(:,:)
    zrg_aux3d(:,iemis,:)  = emis_rad_rg(:,:)
    zrg_aux3d(:,icosmu0,:)= cosmu0_rg(:,:)

  ENDIF


  ! Interpolate reduced-grid fields to full grid

  ! interpol_vec/scal_nudging needs all fields with full boundaries, so the arrays
  ! have to be sync'd before calling these routines.
  ! Please note that we cannot use sync_patch_array here (comparing parallel/non parallel results)
  ! since the arrays don't start with lower bound 1 in the non paralellel case!

  IF  (l_parallel) THEN

    nlev_tot = 4*nlevp1_rg + n2dvars

    CALL exchange_data_mult(p_pp%comm_pat_c, 5, nlev_tot, recv1=p_lwflxclr,       &
                            recv2=p_lwflxall, recv3=p_trsolclr, recv4=p_trsolall, &
                            recv5=zrg_aux3d                                       )

  ENDIF

  IF (p_test_run) THEN
    trsolall = 0._wp
    trsolclr = 0._wp
    lwflxall = 0._wp
    lwflxclr = 0._wp
    z_aux3d  = 0._wp
  ENDIF

  l_limit(1:2) = .TRUE.      ! limit transmissivity to positive values
  l_limit(3:5) = .FALSE.
  rlimval(:)   = 2.94e-37_wp ! seems to be the lower threshold for SW transmissivity
                             ! in the RRTM scheme

  CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, nshift, 5, 1, &
    &                         f3din1=p_trsolall, f3dout1=trsolall,                      &
    &                         f3din2=p_trsolclr, f3dout2=trsolclr,                      &
    &                         f3din3=p_lwflxall, f3dout3=lwflxall,                      &
    &                         f3din4=p_lwflxclr, f3dout4=lwflxclr,                      &
    &                         f3din5=zrg_aux3d,  f3dout5=z_aux3d,                       &
    !                              llimit_nneg=l_limit, rlimval=rlimval, overshoot_fac=1.1_wp)
    &                         llimit_nneg=l_limit, rlimval=rlimval, overshoot_fac=1.1_wp)

  aclcov(:,:)        = z_aux3d(:,1,:)
  tsfc_backintp(:,:) = z_aux3d(:,2,:)
  alb_backintp(:,:)  = z_aux3d(:,3,:)

  ! Finally apply empirical downscaling corrections depending on variations
  ! in ground temperature and albedo

  ! Start/End block in the parent domain
  i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
  i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

  pscal = 1._wp/4000._wp ! pressure scale for longwave correction

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jk1,tqv,intclw,intcli,dpresg,pfacswc,pfacswa,     &
!$OMP            dlwem_o_dtg,swfac1,swfac2,lwfac1,lwfac2,dtrans_o_dalb_clr,dtrans_o_dalb_all,   &
!$OMP            pfaclw,intqctot,dlwflxclr_o_dtg,dlwflxall_o_dtg,jc1,jc2,jc3,jc4,jb1,jb2,jb3,&
!$OMP  jb4) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,                               &
                       i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int, i_chidx)

    tqv(:)               = 0._wp
    intclw(:,nlevp1)  = 0._wp
    intcli(:,nlevp1)  = 0._wp
    pfacswc(:,nlevp1) = 1._wp
    pfacswa(:,nlevp1) = 1._wp
    DO jk = nlev,1,-1
      jk1 = jk + nshift
      DO jc = i_startidx, i_endidx
        dpresg        = (p_pres_ifc(jc,jk1+1,jb) - p_pres_ifc(jc,jk1,jb))/grav
        tqv(jc)       = tqv(jc)+p_tot_cld(jc,jk1,jb,iqv)*dpresg
        intclw(jc,jk) = intclw(jc,jk+1)+p_tot_cld(jc,jk1,jb,iqc)*dpresg
        intcli(jc,jk) = intcli(jc,jk+1)+p_tot_cld(jc,jk1,jb,iqi)*dpresg
        pfacswc(jc,jk)= 1._wp-0.08_wp*(p_pres_ifc(jc,nlevp1_rg,jb)-p_pres_ifc(jc,jk1,jb))/&
                        p_pres_ifc(jc,nlevp1_rg,jb)
        pfacswa(jc,jk)= 1._wp-0.16_wp*(p_pres_ifc(jc,nlevp1_rg,jb)-p_pres_ifc(jc,jk1,jb))/&
                        p_pres_ifc(jc,nlevp1_rg,jb)
      ENDDO
    ENDDO

    DO jc = i_startidx, i_endidx
      dlwem_o_dtg(jc) = zrg_aux3d(jc,iemis,jb)*4._wp*stbo*zrg_aux3d(jc,itsfc,jb)**3
      swfac1(jc) = (MAX(1.e-3_wp,p_trsolall(jc,nlevp1_rg,jb))/           &
                    MAX(1.e-3_wp,p_trsolclr(jc,nlevp1_rg,jb)) )**0.36_wp
      swfac2(jc) =  MAX(0.25_wp,3._wp*zrg_aux3d(jc,icosmu0,jb))**0.1_wp
      IF (tqv(jc) > 15._wp) then
        lwfac1(jc) = 1.677_wp*MAX(1._wp,tqv(jc))**(-0.72_wp)
      ELSE
        lwfac1(jc) = 0.4388_wp*MAX(1._wp,tqv(jc))**(-0.225_wp)
      ENDIF
      lwfac2(jc) = 0.92_wp*MAX(1._wp,tqv(jc))**(-0.07_wp)
    ENDDO

    DO jk = 1,nlevp1
      jk1 = jk + nshift
      DO jc = i_startidx, i_endidx
        dtrans_o_dalb_clr(jc,jk) = - p_trsolclr(jc,nlevp1_rg,jb)*pfacswc(jc,jk) / &
                                   ( (1._wp-zrg_aux3d(jc,ialb,jb))*swfac2(jc) )
        dtrans_o_dalb_all(jc,jk) = - p_trsolall(jc,nlevp1_rg,jb)*pfacswa(jc,jk)*swfac1(jc)/ &
                                   ( (1._wp-zrg_aux3d(jc,ialb,jb))*swfac2(jc) )

        pfaclw = lwfac1(jc)+(lwfac2(jc)-lwfac1(jc))*EXP(-SQRT((p_pres_ifc(jc,nlevp1_rg,jb)- &
                 p_pres_ifc(jc,jk1,jb))*pscal))
        intqctot = MIN(0.30119_wp,MAX(1.008e-3_wp,intclw(jc,jk)+0.2_wp*intcli(jc,jk)))

        dlwflxclr_o_dtg(jc,jk) = -dlwem_o_dtg(jc)*pfaclw
        dlwflxall_o_dtg(jc,jk) = dlwflxclr_o_dtg(jc,jk)*(1._wp-(6.9_wp+LOG(intqctot))/5.7_wp)
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

        trsolclr(jc1,jk,jb1) = MAX(trsolclr(jc1,jk,jb1) + dtrans_o_dalb_clr(jc,jk)* &
          ( albvisdif(jc1,jb1) - alb_backintp(jc1,jb1) ), rlimval(2))
        trsolclr(jc2,jk,jb2) = MAX(trsolclr(jc2,jk,jb2) + dtrans_o_dalb_clr(jc,jk)* &
          ( albvisdif(jc2,jb2) - alb_backintp(jc2,jb2) ), rlimval(2))
        trsolclr(jc3,jk,jb3) = MAX(trsolclr(jc3,jk,jb3) + dtrans_o_dalb_clr(jc,jk)* &
          ( albvisdif(jc3,jb3) - alb_backintp(jc3,jb3) ), rlimval(2))
        trsolclr(jc4,jk,jb4) = MAX(trsolclr(jc4,jk,jb4) + dtrans_o_dalb_clr(jc,jk)* &
          ( albvisdif(jc4,jb4) - alb_backintp(jc4,jb4) ), rlimval(2))

        trsolall(jc1,jk,jb1) = MAX(trsolall(jc1,jk,jb1) + dtrans_o_dalb_all(jc,jk)* &
          ( albvisdif(jc1,jb1) - alb_backintp(jc1,jb1) ), rlimval(1))
        trsolall(jc2,jk,jb2) = MAX(trsolall(jc2,jk,jb2) + dtrans_o_dalb_all(jc,jk)* &
          ( albvisdif(jc2,jb2) - alb_backintp(jc2,jb2) ), rlimval(1))
        trsolall(jc3,jk,jb3) = MAX(trsolall(jc3,jk,jb3) + dtrans_o_dalb_all(jc,jk)* &
          ( albvisdif(jc3,jb3) - alb_backintp(jc3,jb3) ), rlimval(1))
        trsolall(jc4,jk,jb4) = MAX(trsolall(jc4,jk,jb4) + dtrans_o_dalb_all(jc,jk)* &
          ( albvisdif(jc4,jb4) - alb_backintp(jc4,jb4) ), rlimval(1))


        lwflxclr(jc1,jk,jb1) = lwflxclr(jc1,jk,jb1) + dlwflxclr_o_dtg(jc,jk)* &
          ( tsfc(jc1,jb1) - tsfc_backintp(jc1,jb1) )
        lwflxclr(jc2,jk,jb2) = lwflxclr(jc2,jk,jb2) + dlwflxclr_o_dtg(jc,jk)* &
          ( tsfc(jc2,jb2) - tsfc_backintp(jc2,jb2) )
        lwflxclr(jc3,jk,jb3) = lwflxclr(jc3,jk,jb3) + dlwflxclr_o_dtg(jc,jk)* &
          ( tsfc(jc3,jb3) - tsfc_backintp(jc3,jb3) )
        lwflxclr(jc4,jk,jb4) = lwflxclr(jc4,jk,jb4) + dlwflxclr_o_dtg(jc,jk)* &
          ( tsfc(jc4,jb4) - tsfc_backintp(jc4,jb4) )

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

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  IF (l_parallel .AND. jgp == 0) THEN
    DEALLOCATE(z_lwflxclr, z_lwflxall, z_trsolclr, z_trsolall, z_pres_ifc, z_tot_cld, zpg_aux3d)
  ENDIF

  DEALLOCATE(zrg_aux3d,z_aux3d)

END SUBROUTINE downscale_rad_output

!>
!! This routine averages the input fields for Ritter-Geleyn radiation to the next coarser grid
!! level. It is an adaptation of upscale_rad_input to the Ritter-Geleyn radiation.  
!!
!! @par Revision History
!! Thorsten Reinhardt, AGeoBw, 2010-12-06
!!
SUBROUTINE upscale_rad_input_rg(jg, jgp, nlev_rg, nlevp1_rg,         &
  &  cosmu0, albvisdir, alb_ther, temp_ifc, dpres_mc,                &
  &  tot_cld, sqv, duco2, duo3,                                      &
  &  aeq1, aeq2, aeq3, aeq4, aeq5, pres_sfc,pres_ifc,                &
  &  rg_cosmu0, rg_albvisdir, rg_alb_ther, rg_temp_ifc, rg_dpres_mc, &
  &  rg_tot_cld, rg_sqv, rg_duco2, rg_duo3,                          &
  &  rg_aeq1, rg_aeq2, rg_aeq3, rg_aeq4, rg_aeq5, rg_pres_sfc )


  ! Input grid parameters
  INTEGER, INTENT(IN)  :: jg, jgp  ! domain IDs of main and reduced grids
  INTEGER, INTENT(IN)  :: nlev_rg, nlevp1_rg  ! number of model levels on reduced grid
  
  ! Input fields (on full grid)
  REAL(wp), INTENT(IN) ::                                                             &
    & cosmu0(:,:), albvisdir(:,:), alb_ther(:,:), temp_ifc(:,:,:),                    &
    & dpres_mc(:,:,:), tot_cld(:,:,:,:), sqv(:,:,:), duco2(:,:,:), duo3(:,:,:),       &
    & aeq1(:,:,:),aeq2(:,:,:),aeq3(:,:,:),aeq4(:,:,:),aeq5(:,:,:),                    &
    & pres_sfc(:,:)

  ! Input field (on full grid) without corresponding output fields on reduced grid
  REAL(wp), INTENT(IN) ::  &
    & pres_ifc(:,:,:)
  
  ! Corresponding output fields (on reduced grid)
  REAL(wp), TARGET, INTENT(OUT) ::                                                            &
    & rg_cosmu0(:,:), rg_albvisdir(:,:), rg_alb_ther(:,:), rg_temp_ifc(:,:,:),                &
    & rg_dpres_mc(:,:,:), rg_tot_cld(:,:,:,:), rg_sqv(:,:,:), rg_duco2(:,:,:), rg_duo3(:,:,:),&
    & rg_aeq1(:,:,:),rg_aeq2(:,:,:),rg_aeq3(:,:,:),rg_aeq4(:,:,:),rg_aeq5(:,:,:),             &
    & rg_pres_sfc(:,:)

  ! Intermediate storage fields needed in the case of MPI parallelization
  REAL(wp), ALLOCATABLE, TARGET ::                                                        &
    & z_cosmu0(:,:), z_albvisdir(:,:), z_alb_ther(:,:), z_temp_ifc(:,:,:),                &
    & z_dpres_mc(:,:,:), z_tot_cld(:,:,:,:), z_sqv(:,:,:), z_duco2(:,:,:), z_duo3(:,:,:), &
    & z_aeq1(:,:,:),z_aeq2(:,:,:),z_aeq3(:,:,:),z_aeq4(:,:,:),z_aeq5(:,:,:),              &
    & z_pres_sfc(:,:), z_aux3d(:,:,:),  zrg_aux3d(:,:,:)

  
  ! Pointers to output fields (no MPI) or intermediate fields (MPI)
  REAL(wp), POINTER ::                                                                    &
    & p_cosmu0(:,:), p_albvisdir(:,:), p_alb_ther(:,:), p_temp_ifc(:,:,:),                &
    & p_dpres_mc(:,:,:), p_tot_cld(:,:,:,:), p_sqv(:,:,:), p_duco2(:,:,:), p_duo3(:,:,:), &
    & p_aeq1(:,:,:),p_aeq2(:,:,:),p_aeq3(:,:,:),p_aeq4(:,:,:),p_aeq5(:,:,:),              &
    & p_pres_sfc(:,:)


  ! Pointers to types needed to minimize code duplication for MPI/no-MPI cases
  TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
  TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
  TYPE(t_patch),      POINTER     :: p_pp  => NULL()

  ! Indices
  INTEGER :: jb, jc, jk, jk1, i_chidx, i_nchdom, &
             i_startblk, i_endblk, i_startidx, i_endidx, nblks_c_lp

  INTEGER :: nlev, nlevp1      !< number of full and half levels
  INTEGER :: nshift, nst
  REAL(wp) :: z_help_pres_ratio, z_rho_1

  INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
  LOGICAL :: l_parallel
  REAL(wp), POINTER :: p_fbkwgt(:,:,:)

  
  !-----------------------------------------------------------------------

  IF (msg_level >= 10) THEN
    WRITE(message_text,'(a,i2,a,i2)') 'Upscaling of radiation input fields',&
      jg,' =>',jgp
    CALL message('upscale_rad_input',message_text)
  ENDIF

  IF (my_process_is_mpi_seq()) THEN
    l_parallel = .FALSE.
  ELSE
    l_parallel = .TRUE.
  ENDIF

  ! For the time being, the radiation grid is assumed to have the same levels
  ! as the full grid
  ! Number of levels of the full grid
  nlev   = p_patch(jg)%nlev
  nlevp1 = p_patch(jg)%nlevp1

  ! nlev_rg carries the number of model levels of the reduced grid,
  ! which may be larger than nlev
  nshift = nlev_rg - nlev ! resulting shift parameter

  IF (l_parallel) THEN
    p_grf => p_grf_state_local_parent(jg)
    p_gcp => p_patch_local_parent(jg)%cells
    p_pp  => p_patch_local_parent(jg)
  ELSE
    p_grf => p_grf_state(jgp)
    p_gcp => p_patch(jgp)%cells
    p_pp  => p_patch(jgp)
  ENDIF

  i_chidx  = p_patch(jg)%parent_child_index
  i_nchdom = MAX(1,p_patch(jgp)%n_childdom)

  ! Set pointers to index and coefficient fields for cell-based variables
  iidx => p_gcp%child_idx
  iblk => p_gcp%child_blk

  p_fbkwgt => p_grf%fbk_wgt_c

  ! layer shift w.r.t. global grid (> 0 in case of vertical nesting)
  nst = p_patch(jg)%nshift_total
  
  ! Allocation of local storage fields at local parent level in MPI-case
  IF (l_parallel .AND. jgp == 0) THEN
    nblks_c_lp = p_gcp%end_blk(min_rlcell,i_chidx)

    ALLOCATE(z_cosmu0(nproma,nblks_c_lp), z_albvisdir(nproma,nblks_c_lp),                  &
             z_alb_ther(nproma,nblks_c_lp), z_temp_ifc(nproma,nlevp1_rg,nblks_c_lp),       &
             z_dpres_mc(nproma,nlev_rg,nblks_c_lp),z_tot_cld(nproma,nlev_rg,nblks_c_lp,4), &
             z_sqv(nproma,nlev_rg,nblks_c_lp),z_duco2(nproma,nlev_rg,nblks_c_lp),          &
             z_duo3(nproma,nlev_rg,nblks_c_lp), z_aeq1(nproma,nlev_rg,nblks_c_lp),         &
             z_aeq2(nproma,nlev_rg,nblks_c_lp), z_aeq3(nproma,nlev_rg,nblks_c_lp),         &
             z_aeq4(nproma,nlev_rg,nblks_c_lp), z_aeq5(nproma,nlev_rg,nblks_c_lp),         &
             z_pres_sfc(nproma,nblks_c_lp), z_aux3d(nproma,6,nblks_c_lp),                  &
             zrg_aux3d(nproma,4,p_patch(jgp)%nblks_c)                                      )

  ENDIF
  
  ! Set pointers to either the parent-level variables (non-MPI case) or to the
  ! intermediate storage fields (MPI case)
  IF (l_parallel .AND. jgp == 0) THEN
    p_cosmu0     => z_cosmu0
    p_albvisdir  => z_albvisdir
    p_alb_ther   => z_alb_ther
    p_temp_ifc   => z_temp_ifc
    p_dpres_mc   => z_dpres_mc
    p_tot_cld    => z_tot_cld
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
    p_albvisdir  => rg_albvisdir
    p_alb_ther   => rg_alb_ther
    p_temp_ifc   => rg_temp_ifc
    p_dpres_mc   => rg_dpres_mc
    p_tot_cld    => rg_tot_cld
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
    p_albvisdir  = 0._wp
    p_alb_ther   = 0._wp
    p_temp_ifc   = 0._wp
    p_dpres_mc   = 0._wp
    p_tot_cld    = 0._wp
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

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,                                &
                       i_startidx, i_endidx, grf_ovlparea_start_c, min_rlcell_int, i_chidx)

    DO jc = i_startidx, i_endidx

      p_cosmu0(jc,jb) =                                         &
        cosmu0(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        cosmu0(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        cosmu0(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        cosmu0(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      p_albvisdir(jc,jb) =                                         &
        albvisdir(iidx(jc,jb,1),iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        albvisdir(iidx(jc,jb,2),iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        albvisdir(iidx(jc,jb,3),iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        albvisdir(iidx(jc,jb,4),iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

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
      
      p_temp_ifc(jc,nlevp1_rg,jb) =                                      &
        temp_ifc(iidx(jc,jb,1),nlevp1,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        temp_ifc(iidx(jc,jb,2),nlevp1,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        temp_ifc(iidx(jc,jb,3),nlevp1,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        temp_ifc(iidx(jc,jb,4),nlevp1,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

    ENDDO

    ! combine 2D fields in a 3D field to speed up MPI communication
    IF (l_parallel .AND. jgp == 0) THEN
      DO jc = i_startidx, i_endidx
        z_aux3d(jc,1,jb) = p_cosmu0(jc,jb)
        z_aux3d(jc,2,jb) = p_albvisdir(jc,jb)
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

!CDIR EXPAND=4
        p_tot_cld(jc,jk1,jb,1:4) =                                        &
          tot_cld(iidx(jc,jb,1),jk,iblk(jc,jb,1),1:4)*p_fbkwgt(jc,jb,1) + &
          tot_cld(iidx(jc,jb,2),jk,iblk(jc,jb,2),1:4)*p_fbkwgt(jc,jb,2) + &
          tot_cld(iidx(jc,jb,3),jk,iblk(jc,jb,3),1:4)*p_fbkwgt(jc,jb,3) + &
          tot_cld(iidx(jc,jb,4),jk,iblk(jc,jb,4),1:4)*p_fbkwgt(jc,jb,4)

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

!CDIR EXPAND=4
        p_tot_cld(jc,1,jb,1:4) = p_tot_cld(jc,2,jb,1:4)
                
      ENDDO
    ENDIF
    
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  IF (l_parallel .AND. jgp == 0) THEN

    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 6, 5*nlev_rg+5, &
                            RECV1=rg_temp_ifc, SEND1=z_temp_ifc,            &
                            RECV2=rg_dpres_mc, SEND2=z_dpres_mc,            &
                            RECV3=rg_sqv,      SEND3=z_sqv,                 &
                            RECV4=rg_duco2,    SEND4=z_duco2,               &
                            RECV5=rg_duo3,     SEND5=z_duo3,                &
                            RECV6=zrg_aux3d,   SEND6=z_aux3d              )
    
    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 5, 5*nlev_rg, &
                            RECV1=rg_aeq1,     SEND1=z_aeq1,              &
                            RECV2=rg_aeq2,     SEND2=z_aeq2,              &
                            RECV3=rg_aeq3,     SEND3=z_aeq3,              &
                            RECV4=rg_aeq4,     SEND4=z_aeq4,              &
                            RECV5=rg_aeq5,     SEND5=z_aeq5               )
    

    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 4, 4*nlev_rg, &
                            RECV4D=rg_tot_cld, SEND4D=z_tot_cld           )


    i_startblk = p_patch(jgp)%cells%start_blk(1,1)
    i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell,i_nchdom)

    ! OpenMP section commented because the DO loop does almost no work (overhead larger than benefit)
!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch(jgp), jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rlcell, i_nchdom)

      DO jc = i_startidx, i_endidx
        rg_cosmu0(jc,jb)    = zrg_aux3d(jc,1,jb)
        rg_albvisdir(jc,jb) = zrg_aux3d(jc,2,jb)
        rg_alb_ther(jc,jb)  = zrg_aux3d(jc,3,jb)
        rg_pres_sfc(jc,jb)  = zrg_aux3d(jc,4,jb)
      ENDDO

    ENDDO
!!$OMP END DO NOWAIT
!!$OMP END PARALLEL

    DEALLOCATE(z_cosmu0, z_albvisdir, z_alb_ther, z_temp_ifc, z_dpres_mc,z_tot_cld, &
      &  z_sqv, z_duco2, z_duo3, z_aeq1, z_aeq2, z_aeq3, z_aeq4, z_aeq5,            &
      &  z_pres_sfc, z_aux3d,zrg_aux3d )

  ENDIF

END SUBROUTINE upscale_rad_input_rg

SUBROUTINE downscale_rad_output_rg( jg, jgp, nlev_rg,                   &
  &  rg_lwflxall, rg_trsolall, tsfc_rg,                                 &
  &  albeff_rg, albefffac_rg,flsp_rg, flsd_rg,                          &
  &  alb_ther_rg, cosmu0_rg, tot_cld_rg,                                &
  &  dpres_mc_rg, pres_sfc_rg,tsfc, albvisdif, zsct, lwflxall, trsolall )


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
  REAL(wp), INTENT(IN) :: tsfc(:,:), albvisdif(:,:)

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

  ! More local variables
  REAL(wp) :: pscal, dpresg, pfaclw, intqctot, dlwflxclr_o_dtg
  
  REAL(wp), DIMENSION(nproma) :: tqv, dlwem_o_dtg, swfac2, lwfac1, lwfac2

  REAL(wp), DIMENSION(nproma,p_patch(jg)%nlevp1) :: intclw, intcli,     &
    dtrans_o_dalb_all, dlwflxall_o_dtg, pfacswc, pfacswa

  
  ! Pointers to types needed to minimize code duplication for MPI/no-MPI cases
  TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
  TYPE(t_int_state),  POINTER     :: p_int => NULL()
  TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
  TYPE(t_patch),      POINTER     :: p_pp => NULL()

  ! Indices
  INTEGER :: i_chidx, i_nchdom, nblks_c_lp
  INTEGER :: jb, jk, jk1, jc, i_startblk, i_endblk, i_startidx, i_endidx
  INTEGER :: jc1, jc2, jc3, jc4, jb1, jb2, jb3, jb4
  

  INTEGER :: nlev, nlevp1      !< number of full and half levels
  INTEGER :: nshift, nlevp1_rg, nlev_tot
  INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk

  INTEGER :: n2dvars, ipsfc, itsfc, ialb, ialbfac, ialbther, icosmu0, iflsp, iflsd
  
  LOGICAL :: l_parallel, l_limit(3)
  REAL(wp) :: rlimval(3)
!-----------------------------------------------------------------------

  IF (msg_level >= 10) THEN
    WRITE(message_text,'(a,i2,a,i2)') 'Downscaling of radiation output fields',&
      jgp,' =>',jg
    CALL message('downscale_rad_output',message_text)
  ENDIF

  IF (my_process_is_mpi_seq()) THEN
    l_parallel = .FALSE.
  ELSE
    l_parallel = .TRUE.
  ENDIF


  ! For the time being, the radiation grid is assumed to have the same levels
  ! as the full grid
  nlev   = p_patch(jg)%nlev
  nlevp1 = p_patch(jg)%nlevp1

  ! nlev_rg carries the number of model levels of the reduced grid,
  ! which may be larger than nlev
  nlevp1_rg = nlev_rg + 1
  nshift = nlev_rg - nlev ! resulting shift parameter

  IF (l_parallel) THEN
    p_grf => p_grf_state_local_parent(jg)
    p_int => p_int_state_local_parent(jg)
    p_gcp => p_patch_local_parent(jg)%cells
    p_pp  => p_patch_local_parent(jg)
  ELSE
    p_grf => p_grf_state(jgp)
    p_int => p_int_state(jgp)
    p_gcp => p_patch(jgp)%cells
    p_pp  => p_patch(jgp)
  ENDIF

  i_nchdom = MAX(1,p_patch(jg)%n_childdom)
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
  IF (l_parallel .AND. jgp == 0) THEN

    ALLOCATE(z_lwflxall(nproma,nlevp1_rg,nblks_c_lp), z_trsolall(nproma,nlevp1_rg,nblks_c_lp) , &
      &      z_dpres_mc(nproma,nlev_rg,nblks_c_lp), z_tot_cld(nproma,nlev_rg,nblks_c_lp,4),     &
      &      zpg_aux3d(nproma,n2dvars,p_patch(jgp)%nblks_c),                                    &
      &      zrg_aux3d(nproma,n2dvars,nblks_c_lp),   z_aux3d(nproma,8,p_patch(jg)%nblks_c),     &
      &      pres_ifc(nproma,nlevp1_rg,nblks_c_lp) )

  ELSE

    ALLOCATE(zrg_aux3d(nproma,n2dvars,nblks_c_lp),z_aux3d(nproma,8,p_patch(jg)%nblks_c), &
      &      pres_ifc(nproma,nlevp1_rg,nblks_c_lp) )
    
  ENDIF

  ! Perform communication from parent to local parent grid in the MPI case,
  ! and set pointers such that further processing is the same for MPI / non-MPI cases
  IF (l_parallel .AND. jgp == 0) THEN

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

    CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, 4, 4*nlev_rg, &
                            RECV4D=z_tot_cld, SEND4D=tot_cld_rg       )
    
    

    p_lwflxall   => z_lwflxall
    p_trsolall   => z_trsolall
    p_dpres_mc   => z_dpres_mc
    p_tot_cld    => z_tot_cld
     
 
  ELSE

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

  IF  (l_parallel) THEN

    nlev_tot = 2*nlevp1_rg + n2dvars

    CALL exchange_data_mult(p_pp%comm_pat_c, 3 ,nlev_tot, recv1=p_lwflxall, recv2=p_trsolall, &
      & recv3=zrg_aux3d )

  ENDIF
  
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

  pscal = 1._wp/4000._wp ! pressure scale for longwave correction


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,jk1,tqv,intclw,intcli,dpresg,pfacswc,pfacswa,     &
!$OMP            dlwem_o_dtg,swfac2,lwfac1,lwfac2,dtrans_o_dalb_all,            &
!$OMP            pfaclw,intqctot,dlwflxclr_o_dtg,dlwflxall_o_dtg,jc1,jc2,jc3,jc4,jb1,jb2,jb3, &
!$OMP jb4) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,                               &
                       i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int, i_chidx)

    tqv(:)            = 0._wp
    intclw(:,nlevp1)  = 0._wp
    intcli(:,nlevp1)  = 0._wp
    pfacswc(:,nlevp1) = 1._wp
    pfacswa(:,nlevp1) = 1._wp
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
        pfacswc(jc,jk)= 1._wp-0.08_wp*(pres_ifc(jc,nlevp1_rg,jb)-pres_ifc(jc,jk1,jb))/&
                        pres_ifc(jc,nlevp1_rg,jb)
        pfacswa(jc,jk)= 1._wp-0.16_wp*(pres_ifc(jc,nlevp1_rg,jb)-pres_ifc(jc,jk1,jb))/&
                        pres_ifc(jc,nlevp1_rg,jb)
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
          ( albfac_backintp(jc1,jb1)*albvisdif(jc1,jb1) - alb_backintp(jc1,jb1) ), rlimval(1))
        trsolall(jc2,jk,jb2) = MAX(trsolall(jc2,jk,jb2) + dtrans_o_dalb_all(jc,jk)* &
          ( albfac_backintp(jc2,jb2)*albvisdif(jc2,jb2) - alb_backintp(jc2,jb2) ), rlimval(1))
        trsolall(jc3,jk,jb3) = MAX(trsolall(jc3,jk,jb3) + dtrans_o_dalb_all(jc,jk)* &
          ( albfac_backintp(jc3,jb3)*albvisdif(jc3,jb3) - alb_backintp(jc3,jb3) ), rlimval(1))
        trsolall(jc4,jk,jb4) = MAX(trsolall(jc4,jk,jb4) + dtrans_o_dalb_all(jc,jk)* &
          ( albfac_backintp(jc4,jb4)*albvisdif(jc4,jb4) - alb_backintp(jc4,jb4) ), rlimval(1))

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
  
  IF (l_parallel .AND. jgp == 0) THEN
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
SUBROUTINE interpol_phys_grf (jg,jgc,jn,lsfc_interp)

  USE mo_nwp_phy_state,      ONLY: prm_diag

  ! Input:
  INTEGER, INTENT(in) :: jg,jgc,jn
  LOGICAL, INTENT(in) :: lsfc_interp

  ! Pointers
  TYPE(t_patch),                POINTER :: ptr_pp
  TYPE(t_patch),                POINTER :: ptr_pc
  TYPE(t_gridref_single_state), POINTER :: ptr_grf
  TYPE(t_int_state),            POINTER :: ptr_int
  TYPE(t_lnd_diag),             POINTER :: ptr_ldiagp ! parent level land diag state
  TYPE(t_lnd_diag),             POINTER :: ptr_ldiagc ! child level land diag state

  ! Local fields
  INTEGER, PARAMETER  :: nfields_p=9    ! Number of 2D phyiscs fields for which boundary interpolation is needed
  INTEGER, PARAMETER  :: nfields_l2=12   ! Number of 2D land state fields

  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, jb, jc, jk


  ! Temporary storage to do boundary interpolation for all 2D fields in one step
  REAL(wp) :: z_aux3dp_p(nproma,nfields_p,p_patch(jg)%nblks_c),         &  ! 2D physics diag fields
              z_aux3dp_c(nproma,nfields_p,p_patch(jgc)%nblks_c),        &
              z_aux3dl2_p(nproma,nfields_l2,p_patch(jg)%nblks_c),       &  ! 2D land state fields
              z_aux3dl2_c(nproma,nfields_l2,p_patch(jgc)%nblks_c),      &
              z_aux3dso_p(nproma,3*(nlev_soil+1),p_patch(jg)%nblks_c),  &  ! 3D land state fields for soil
              z_aux3dso_c(nproma,3*(nlev_soil+1),p_patch(jgc)%nblks_c), &
              z_aux3dsn_p(nproma,5*nlev_snow,p_patch(jg)%nblks_c),      &  ! 3D land state fields for multi-layer snow
              z_aux3dsn_c(nproma,5*nlev_snow,p_patch(jgc)%nblks_c)         ! (used if lmulti_snow = ture))

  ! set pointers
  ptr_pp  => p_patch(jg)
  ptr_pc  => p_patch(jgc)
  ptr_grf => p_grf_state(jg)%p_dom(jn)
  ptr_int => p_int_state(jg)

  IF (lsfc_interp) THEN
    ptr_ldiagp => p_lnd_state(jg)%diag_lnd
    ptr_ldiagc => p_lnd_state(jgc)%diag_lnd
  ENDIF

  IF (p_test_run) THEN
     z_aux3dp_p(:,:,:) = 0._wp
     z_aux3dl2_p(:,:,:) = 0._wp
     z_aux3dso_p(:,:,:) = 0._wp
     z_aux3dsn_p(:,:,:) = 0._wp
  ENDIF

  i_startblk = ptr_pp%cells%start_blk(1,1)
  i_endblk   = ptr_pp%nblks_c


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, 1)

    DO jc = i_startidx, i_endidx
      z_aux3dp_p(jc,1,jb) = prm_diag(jg)%tot_prec(jc,jb)
      z_aux3dp_p(jc,2,jb) = prm_diag(jg)%rain_gsp(jc,jb)
      z_aux3dp_p(jc,3,jb) = prm_diag(jg)%snow_gsp(jc,jb)
      z_aux3dp_p(jc,4,jb) = prm_diag(jg)%rain_con(jc,jb)
      z_aux3dp_p(jc,5,jb) = prm_diag(jg)%snow_con(jc,jb)
      z_aux3dp_p(jc,6,jb) = prm_diag(jg)%tracer_rate(jc,jb,1)
      z_aux3dp_p(jc,7,jb) = prm_diag(jg)%tracer_rate(jc,jb,2)
      z_aux3dp_p(jc,8,jb) = prm_diag(jg)%tracer_rate(jc,jb,3)
      z_aux3dp_p(jc,9,jb) = prm_diag(jg)%tracer_rate(jc,jb,4)
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
        z_aux3dl2_p(jc,11,jb) = ptr_ldiagp%t_so(jc,nlev_soil+2,jb)
        IF (lmulti_snow) &
          z_aux3dl2_p(jc,12,jb) = ptr_ldiagp%t_snow_mult(jc,nlev_snow+1,jb)
      ENDDO

      DO jk = 1, nlev_soil+1
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

      CALL sync_patch_array_mult(SYNC_C,ptr_pp,4,z_aux3dp_p,z_aux3dl2_p,z_aux3dso_p,z_aux3dsn_p)
      CALL interpol_scal_grf (ptr_pp, ptr_pc, ptr_int, ptr_grf, jn, 4, z_aux3dp_p, z_aux3dp_c, &
        z_aux3dl2_p, z_aux3dl2_c, z_aux3dso_p, z_aux3dso_c, z_aux3dsn_p, z_aux3dsn_c,          &
        llimit_nneg=(/.TRUE.,.TRUE.,.TRUE.,.TRUE./), lnoshift=.TRUE.)

    ELSE IF (lsfc_interp) THEN

      CALL sync_patch_array_mult(SYNC_C,ptr_pp,3,z_aux3dp_p,z_aux3dl2_p,z_aux3dso_p)
      CALL interpol_scal_grf (ptr_pp, ptr_pc, ptr_int, ptr_grf, jn, 3, z_aux3dp_p, z_aux3dp_c, &
        z_aux3dl2_p, z_aux3dl2_c, z_aux3dso_p, z_aux3dso_c,                                    &
        llimit_nneg=(/.TRUE.,.TRUE.,.TRUE./), lnoshift=.TRUE.)

    ELSE

      CALL sync_patch_array(SYNC_C,ptr_pp,z_aux3dp_p)
      CALL interpol_scal_grf (ptr_pp, ptr_pc, ptr_int, ptr_grf, jn, 1, &
        z_aux3dp_p, z_aux3dp_c, llimit_nneg=(/.TRUE./), lnoshift=.TRUE.)

    ENDIF

  i_startblk = ptr_pc%cells%start_blk(1,1)
  i_endblk   = ptr_pc%cells%end_blk(grf_bdywidth_c,1)


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_pc, jb, i_startblk, i_endblk,        &
                       i_startidx, i_endidx, 1, grf_bdywidth_c)

    DO jc = i_startidx, i_endidx
      prm_diag(jgc)%tot_prec(jc,jb)       = z_aux3dp_c(jc,1,jb)
      prm_diag(jgc)%rain_gsp(jc,jb)       = z_aux3dp_c(jc,2,jb)
      prm_diag(jgc)%snow_gsp(jc,jb)       = z_aux3dp_c(jc,3,jb)
      prm_diag(jgc)%rain_con(jc,jb)       = z_aux3dp_c(jc,4,jb)
      prm_diag(jgc)%snow_con(jc,jb)       = z_aux3dp_c(jc,5,jb)
      prm_diag(jgc)%tracer_rate(jc,jb,1)  = z_aux3dp_c(jc,6,jb)
      prm_diag(jgc)%tracer_rate(jc,jb,2)  = z_aux3dp_c(jc,7,jb)
      prm_diag(jgc)%tracer_rate(jc,jb,3)  = z_aux3dp_c(jc,8,jb)
      prm_diag(jgc)%tracer_rate(jc,jb,4)  = z_aux3dp_c(jc,9,jb)
    ENDDO



    IF (lsfc_interp) THEN
      DO jc = i_startidx, i_endidx
        ptr_ldiagc%t_snow(jc,jb)           = z_aux3dl2_c(jc,1,jb)
        ptr_ldiagc%t_s(jc,jb)              = z_aux3dl2_c(jc,2,jb) 
        ptr_ldiagc%w_snow(jc,jb)           = z_aux3dl2_c(jc,3,jb)
        ptr_ldiagc%rho_snow(jc,jb)         = z_aux3dl2_c(jc,4,jb) 
        ptr_ldiagc%w_i(jc,jb)              = z_aux3dl2_c(jc,5,jb) 
        ptr_ldiagc%h_snow(jc,jb)           = z_aux3dl2_c(jc,6,jb)
        ptr_ldiagc%freshsnow(jc,jb)        = z_aux3dl2_c(jc,7,jb)
        ptr_ldiagc%snowfrac(jc,jb)         = z_aux3dl2_c(jc,8,jb) 
        ptr_ldiagc%runoff_s(jc,jb)         = z_aux3dl2_c(jc,9,jb)
        ptr_ldiagc%runoff_g(jc,jb)         = z_aux3dl2_c(jc,10,jb) 
        ptr_ldiagc%t_so(jc,nlev_soil+2,jb) = z_aux3dl2_c(jc,11,jb)
        IF (lmulti_snow) &
          ptr_ldiagc%t_snow_mult(jc,nlev_snow+1,jb) = z_aux3dl2_c(jc,12,jb) 
      ENDDO

      DO jk = 1, nlev_soil+1
        DO jc = i_startidx, i_endidx
          ptr_ldiagc%t_so(jc,jk,jb)     = z_aux3dso_c(jc,3*(jk-1)+1,jb) 
          ptr_ldiagc%w_so(jc,jk,jb)     = z_aux3dso_c(jc,3*(jk-1)+2,jb) 
          ptr_ldiagc%w_so_ice(jc,jk,jb) = z_aux3dso_c(jc,3*(jk-1)+3,jb) 
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
  TYPE(t_int_state),            POINTER :: ptr_int
  TYPE(t_nwp_phy_diag),         POINTER :: prm_diagp
  TYPE(t_nwp_phy_diag),         POINTER :: prm_diagc
  TYPE(t_lnd_prog),             POINTER :: ptr_lprogp
  TYPE(t_lnd_prog),             POINTER :: ptr_lprogc_t1
  TYPE(t_lnd_prog),             POINTER :: ptr_lprogc_t2


  ! Local fields
  INTEGER, PARAMETER  :: nfields=2    ! Number of 2D fields for which boundary interpolation is needed
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, jb, jc

  ! Temporary storage to do boundary interpolation for all 2D fields in one step
  REAL(wp) :: z_aux3d_p(nproma,nfields,p_patch(jg)%nblks_c), &
              z_aux3d_c(nproma,nfields,p_patch(jgc)%nblks_c)

  ! set pointers
  ptr_pp        => p_patch(jg)
  ptr_pc        => p_patch(jgc)
  ptr_grf       => p_grf_state(jg)%p_dom(jn)
  ptr_int       => p_int_state(jg)
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

  ! OpenMP section commented because the DO loop does almost no work (overhead larger than benefit)
!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, 1)

    DO jc = i_startidx, i_endidx

      z_aux3d_p(jc,1,jb) = ptr_lprogp%t_g(jc,jb)
      z_aux3d_p(jc,2,jb) = prm_diagp%albvisdif(jc,jb)
    ENDDO
  ENDDO
!!$OMP END DO NOWAIT
!!$OMP END PARALLEL

    ! Halo update is needed before interpolation
    CALL sync_patch_array(SYNC_C,ptr_pp,z_aux3d_p)

    CALL interpol_scal_grf (ptr_pp, ptr_pc, ptr_int, ptr_grf, jn, 1,     &
      &                     z_aux3d_p, z_aux3d_c, llimit_nneg=(/.TRUE./),&
      &                     lnoshift=.TRUE.)


  i_startblk = ptr_pc%cells%start_blk(1,1)
  i_endblk   = ptr_pc%cells%end_blk(grf_bdywidth_c,1)

  ! Note: prognostic land fields are set on both time levels to safely avoid
  ! errors when radiation calls for parent and child grids are not properly 
  ! synchronized
  ! OpenMP section commented because the DO loop does almost no work (overhead larger than benefit)
!!$OMP PARALLEL
!!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_pc, jb, i_startblk, i_endblk,        &
                       i_startidx, i_endidx, 1, grf_bdywidth_c)

    DO jc = i_startidx, i_endidx

      ptr_lprogc_t1%t_g(jc,jb)      = z_aux3d_c(jc,1,jb)
      ptr_lprogc_t2%t_g(jc,jb)      = z_aux3d_c(jc,1,jb)
      prm_diagc%albvisdif(jc,jb)    = z_aux3d_c(jc,2,jb)

    ENDDO
  ENDDO
!!$OMP END DO NOWAIT
!!$OMP END PARALLEL

END SUBROUTINE interpol_rrg_grf



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
  TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
  TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
  TYPE(t_patch),      POINTER     :: p_pp => NULL()

  ! Indices
  INTEGER :: jb, jc, i_chidx, i_nchdom, &
             i_startblk, i_endblk, i_startidx, i_endidx, nblks_c_lp

  INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
  LOGICAL :: l_parallel
  REAL(wp), POINTER :: p_fbkwgt(:,:,:), p_aux3d(:,:,:)
  REAL(wp), ALLOCATABLE, TARGET :: z_aux3d_lp(:,:,:), z_aux3d_par(:,:,:)

  !-----------------------------------------------------------------------

  IF (msg_level >= 10) THEN
    WRITE(message_text,'(a,i2,a,i2)') 'Feedback of diagnostic physics fields',&
      p_patch(jg)%id,' =>',p_patch(jgp)%id
    CALL message('feedback_phys_diag',message_text)
  ENDIF

  IF (my_process_is_mpi_seq()) THEN
    l_parallel = .FALSE.
  ELSE
    l_parallel = .TRUE.
  ENDIF

  IF (l_parallel) THEN
    p_grf => p_grf_state_local_parent(jg)
    p_gcp => p_patch_local_parent(jg)%cells
    p_pp  => p_patch_local_parent(jg)
  ELSE
    p_grf => p_grf_state(jgp)
    p_gcp => p_patch(jgp)%cells
    p_pp  => p_patch(jgp)
  ENDIF

  i_chidx  = p_patch(jg)%parent_child_index
  i_nchdom = MAX(1,p_patch(jgp)%n_childdom)

  ! Set pointers to index and coefficient fields for cell-based variables
  iidx => p_gcp%child_idx
  iblk => p_gcp%child_blk

  p_fbkwgt => p_grf%fbk_wgt_c

  ! Allocation of local storage fields 
  nblks_c_lp = p_gcp%end_blk(min_rlcell,i_chidx)
  ALLOCATE(z_aux3d_lp(nproma,7,nblks_c_lp), z_aux3d_par(nproma,7,p_patch(jgp)%nblks_c))
  IF (l_parallel) THEN
    p_aux3d => z_aux3d_lp
  ELSE
    p_aux3d => z_aux3d_par
  ENDIF

  ! Now average input fields to parent grid cells

  ! Start/End block in the parent domain
  i_startblk = p_gcp%start_blk(grf_fbk_start_c,i_chidx)
  i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_pp, jb, i_startblk, i_endblk,                           &
                       i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int, i_chidx)

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

      p_aux3d(jc,6,jb) =                                         &
        prm_diag(jg)%tracer_rate(iidx(jc,jb,1),iblk(jc,jb,1),1)*p_fbkwgt(jc,jb,1) + &
        prm_diag(jg)%tracer_rate(iidx(jc,jb,2),iblk(jc,jb,2),1)*p_fbkwgt(jc,jb,2) + &
        prm_diag(jg)%tracer_rate(iidx(jc,jb,3),iblk(jc,jb,3),1)*p_fbkwgt(jc,jb,3) + &
        prm_diag(jg)%tracer_rate(iidx(jc,jb,4),iblk(jc,jb,4),1)*p_fbkwgt(jc,jb,4)

      p_aux3d(jc,7,jb) =                                         &
        prm_diag(jg)%tracer_rate(iidx(jc,jb,1),iblk(jc,jb,1),2)*p_fbkwgt(jc,jb,1) + &
        prm_diag(jg)%tracer_rate(iidx(jc,jb,2),iblk(jc,jb,2),2)*p_fbkwgt(jc,jb,2) + &
        prm_diag(jg)%tracer_rate(iidx(jc,jb,3),iblk(jc,jb,3),2)*p_fbkwgt(jc,jb,3) + &
        prm_diag(jg)%tracer_rate(iidx(jc,jb,4),iblk(jc,jb,4),2)*p_fbkwgt(jc,jb,4)

    ENDDO

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  IF (l_parallel) THEN

    CALL exchange_data(p_pp%comm_pat_loc_to_glb_c_fbk, RECV=z_aux3d_par, SEND=z_aux3d_lp)
    p_aux3d => z_aux3d_par

  ENDIF

  i_startblk = p_patch(jgp)%cells%start_blk(grf_fbk_start_c,i_chidx)
  i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell_int,i_chidx)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(p_patch(jgp), jb, i_startblk, i_endblk, &
                       i_startidx, i_endidx, grf_fbk_start_c, min_rlcell_int, i_chidx)

    DO jc = i_startidx, i_endidx
      prm_diag(jgp)%tot_prec(jc,jb)      = p_aux3d(jc,1,jb)
      prm_diag(jgp)%rain_con(jc,jb)      = p_aux3d(jc,2,jb)
      prm_diag(jgp)%snow_con(jc,jb)      = p_aux3d(jc,3,jb)
      prm_diag(jgp)%rain_gsp(jc,jb)      = p_aux3d(jc,4,jb)
      prm_diag(jgp)%snow_gsp(jc,jb)      = p_aux3d(jc,5,jb)
      prm_diag(jgp)%tracer_rate(jc,jb,1) = p_aux3d(jc,6,jb)
      prm_diag(jgp)%tracer_rate(jc,jb,2) = p_aux3d(jc,7,jb)
    ENDDO

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  DEALLOCATE(z_aux3d_lp, z_aux3d_par)

END SUBROUTINE feedback_phys_diag


END MODULE mo_phys_nest_utilities

