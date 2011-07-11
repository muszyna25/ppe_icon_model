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
MODULE mo_phys_nest_utilities
!
!
USE mo_kind,                ONLY: wp
USE mo_exception,           ONLY: message_text, message
USE mo_model_domain,        ONLY: t_patch, t_grid_cells
USE mo_model_domain_import, ONLY: n_dom, n_dom_start
USE mo_interpolation,       ONLY: t_int_state
USE mo_grf_interpolation,   ONLY: t_gridref_state, t_gridref_single_state
USE mo_grf_bdyintp,         ONLY: interpol_scal_grf
USE mo_grf_nudgintp,        ONLY: interpol_scal_nudging
USE mo_parallel_configuration,  ONLY: nproma, p_test_run
USE mo_run_nml,             ONLY: msg_level
USE mo_nwp_phy_state,       ONLY: prm_diag
USE mo_impl_constants,      ONLY: min_rlcell, min_rlcell_int
USE mo_loopindices,         ONLY: get_indices_c
USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c, grf_ovlparea_start_c, grf_fbk_start_c
USE mo_mpi,                 ONLY: my_process_is_mpi_seq
USE mo_communication,       ONLY: exchange_data, exchange_data_mult
USE mo_sync,                ONLY: SYNC_C, sync_patch_array
USE mo_subdivision,         ONLY: p_patch_local_parent, p_int_state_local_parent, &
                                  p_grf_state_local_parent

IMPLICIT NONE

PRIVATE

CHARACTER(len=*), PARAMETER :: version = '$Id$'

PUBLIC :: upscale_rad_input, downscale_rad_output, interpol_phys_grf, feedback_phys_diag, &
  &       upscale_rad_input_rg, downscale_rad_output_rg

CONTAINS

!>
!! This routine averages the input fields for RRTM radiation to the next coarser grid level.
!!
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2010-12-01
!!
SUBROUTINE upscale_rad_input(p_patch, p_par_patch, p_par_grf,        &
  fr_land, fr_glac,                                                  &
  cosmu0, albvisdir, albnirdir, albvisdif, albnirdif,                &
  tsfc, pres_ifc, pres, temp, acdnc, tot_cld, q_o3,                  &
  rg_fr_land, rg_fr_glac,                                            &
  rg_cosmu0, rg_albvisdir, rg_albnirdir, rg_albvisdif, rg_albnirdif, &
  rg_tsfc, rg_pres_ifc, rg_pres, rg_temp, rg_acdnc, rg_tot_cld, rg_q_o3 )

  ! Input types
  TYPE(t_patch),       TARGET, INTENT(IN)    ::  p_patch
  TYPE(t_patch),       TARGET, INTENT(IN)    ::  p_par_patch
  TYPE(t_gridref_state), TARGET, INTENT(IN)  ::  p_par_grf

  ! Other input fields (on full grid)
  REAL(wp), INTENT(IN) ::                                                             &
    fr_land(:,:), fr_glac(:,:),                                                       &
    cosmu0(:,:), albvisdir(:,:), albnirdir(:,:), albvisdif(:,:), albnirdif(:,:),      &
    tsfc(:,:), pres_ifc(:,:,:), pres(:,:,:), temp(:,:,:), acdnc(:,:,:), tot_cld(:,:,:,:), &
    q_o3(:,:,:)

  ! Corresponding output fields (on reduced grid)
  REAL(wp), TARGET, INTENT(OUT) ::                                           &
    rg_fr_land(:,:),rg_fr_glac(:,:),                                         &
    rg_cosmu0(:,:), rg_albvisdir(:,:), rg_albnirdir(:,:), rg_albvisdif(:,:), &
    rg_albnirdif(:,:), rg_tsfc(:,:), rg_pres_ifc(:,:,:), rg_pres(:,:,:),     &
    rg_temp(:,:,:), rg_acdnc(:,:,:), rg_tot_cld(:,:,:,:), rg_q_o3(:,:,:)

  ! Intermediate storage fields needed in the case of MPI parallelization
  REAL(wp), ALLOCATABLE, TARGET ::                                       &
    z_fr_land(:,:),z_fr_glac(:,:),                                            &
    z_cosmu0(:,:), z_albvisdir(:,:), z_albnirdir(:,:), z_albvisdif(:,:), &
    z_albnirdif(:,:), z_tsfc(:,:), z_pres_ifc(:,:,:), z_pres(:,:,:),     &
    z_temp(:,:,:), z_acdnc(:,:,:), z_tot_cld(:,:,:,:), z_q_o3(:,:,:),    &
    z_aux3d(:,:,:), zrg_aux3d(:,:,:)

  ! Pointers to output fields (no MPI) or intermediate fields (MPI)
  REAL(wp), POINTER ::                                                   &
    p_fr_land(:,:),p_fr_glac(:,:),                                       &
    p_cosmu0(:,:), p_albvisdir(:,:), p_albnirdir(:,:), p_albvisdif(:,:), &
    p_albnirdif(:,:), p_tsfc(:,:), p_pres_ifc(:,:,:), p_pres(:,:,:),     &
    p_temp(:,:,:), p_acdnc(:,:,:), p_tot_cld(:,:,:,:), p_q_o3(:,:,:)


  ! Pointers to types needed to minimize code duplication for MPI/no-MPI cases
  TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
  TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
  TYPE(t_patch),      POINTER     :: p_pp => NULL()

  ! Indices
  INTEGER :: jb, jc, jk, jg, jgp, i_chidx, i_nchdom, &
             i_startblk, i_endblk, i_startidx, i_endidx, nblks_c_lp

  INTEGER :: nlev, nlevp1      !< number of full and half levels

  INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
  LOGICAL :: l_parallel
  REAL(wp), POINTER :: p_fbkwgt(:,:,:)

  !-----------------------------------------------------------------------

  IF (msg_level >= 10) THEN
    WRITE(message_text,'(a,i2,a,i2)') 'Upscaling of radiation input fields',&
      p_patch%id,' =>',p_par_patch%id
    CALL message('upscale_rad_input',message_text)
  ENDIF

  IF (my_process_is_mpi_seq()) THEN
    l_parallel = .FALSE.
  ELSE
    l_parallel = .TRUE.
  ENDIF

  jgp = p_par_patch%id

  ! For the time being, the radiation grid is assumed to have the same levels
  ! as the full grid
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1

  IF (l_parallel) THEN
    jg = p_patch%id
    p_grf => p_grf_state_local_parent(jg)
    p_gcp => p_patch_local_parent(jg)%cells
    p_pp  => p_patch_local_parent(jg)
  ELSE
    p_grf => p_par_grf
    p_gcp => p_par_patch%cells
    p_pp  => p_par_patch
  ENDIF

  i_chidx  = p_patch%parent_child_index
  i_nchdom = MAX(1,p_par_patch%n_childdom)

  ! Set pointers to index and coefficient fields for cell-based variables
  iidx => p_gcp%child_idx
  iblk => p_gcp%child_blk

  p_fbkwgt => p_grf%fbk_wgt_c

  ! Allocation of local storage fields at local parent level in MPI-case
  IF (l_parallel .AND. jgp == 0) THEN
    nblks_c_lp = p_gcp%end_blk(min_rlcell,i_chidx)

    ALLOCATE(z_fr_land(nproma,nblks_c_lp), z_fr_glac(nproma,nblks_c_lp),           &
             z_cosmu0(nproma,nblks_c_lp), z_albvisdir(nproma,nblks_c_lp),          &
             z_albnirdir(nproma,nblks_c_lp), z_albvisdif(nproma,nblks_c_lp),       &
             z_albnirdif(nproma,nblks_c_lp), z_tsfc(nproma,nblks_c_lp),            &
             z_pres_ifc(nproma,nlevp1,nblks_c_lp), z_pres(nproma,nlev,nblks_c_lp), &
             z_temp(nproma,nlev,nblks_c_lp), z_acdnc(nproma,nlev,nblks_c_lp),      &
             z_tot_cld(nproma,nlev,nblks_c_lp,4), z_q_o3(nproma,nlev,nblks_c_lp),  &
             z_aux3d(nproma,8,nblks_c_lp), zrg_aux3d(nproma,8,p_par_patch%nblks_c) )

  ENDIF

  ! Set pointers to either the parent-level variables (non-MPI case) or to the
  ! intermediate storage fields (MPI case)
  IF (l_parallel .AND. jgp == 0) THEN
    p_fr_land    => z_fr_land
    p_fr_glac    => z_fr_glac
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
  ELSE
    p_fr_land    => rg_fr_land
    p_fr_glac    => rg_fr_glac
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
  ENDIF

  IF (p_test_run) THEN
    p_fr_land    = 0._wp
    p_fr_glac    = 0._wp
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
  ENDIF


  ! Now average input fields to parent grid cells

  ! Start/End block in the parent domain
  i_startblk = p_gcp%start_blk(grf_ovlparea_start_c,i_chidx)
  i_endblk   = p_gcp%end_blk(min_rlcell_int,i_chidx)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
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

      p_pres_ifc(jc,nlevp1,jb) =                                           &
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
      ENDDO
    ENDIF

    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx

        p_pres_ifc(jc,jk,jb) =                                         &
          pres_ifc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          pres_ifc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          pres_ifc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          pres_ifc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_pres(jc,jk,jb) =                                         &
          pres(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          pres(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          pres(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          pres(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_temp(jc,jk,jb) =                                         &
          temp(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          temp(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          temp(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          temp(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_q_o3(jc,jk,jb) =                                         &
          q_o3(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          q_o3(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          q_o3(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          q_o3(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_acdnc(jc,jk,jb) =                                         &
          acdnc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          acdnc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          acdnc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          acdnc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

      ENDDO
    ENDDO

    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx

!CDIR EXPAND=4
        p_tot_cld(jc,jk,jb,1:4) =                                         &
          tot_cld(iidx(jc,jb,1),jk,iblk(jc,jb,1),1:4)*p_fbkwgt(jc,jb,1) + &
          tot_cld(iidx(jc,jb,2),jk,iblk(jc,jb,2),1:4)*p_fbkwgt(jc,jb,2) + &
          tot_cld(iidx(jc,jb,3),jk,iblk(jc,jb,3),1:4)*p_fbkwgt(jc,jb,3) + &
          tot_cld(iidx(jc,jb,4),jk,iblk(jc,jb,4),1:4)*p_fbkwgt(jc,jb,4)

      ENDDO
    ENDDO

  ENDDO
!$OMP END DO
!$OMP END PARALLEL


  IF (l_parallel .AND. jgp == 0) THEN

    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 6, 6*nlev+9,  &
                            RECV1=rg_pres_ifc, SEND1=z_pres_ifc,          &
                            RECV2=rg_pres,     SEND2=z_pres,              &
                            RECV3=rg_temp,     SEND3=z_temp,              &
                            RECV4=rg_acdnc,    SEND4=z_acdnc,             &
                            RECV5=zrg_aux3d,   SEND5=z_aux3d,             &
                            RECV6=rg_q_o3,     SEND6=z_q_o3               )

    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 4, 4*nlev,    &
                            RECV4D=rg_tot_cld, SEND4D=z_tot_cld           )


    i_startblk = p_par_patch%cells%start_blk(1,1)
    i_endblk   = p_par_patch%cells%end_blk(min_rlcell,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_par_patch, jb, i_startblk, i_endblk, &
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
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    DEALLOCATE(z_fr_land, z_fr_glac, z_cosmu0, z_albvisdir, z_albnirdir, z_albvisdif,  &
      & z_albnirdif, z_tsfc, z_pres_ifc, z_pres, z_temp, z_acdnc, z_tot_cld, z_q_o3,   &
      & z_aux3d, zrg_aux3d )


  ENDIF

END SUBROUTINE upscale_rad_input

!>
!! This routine interpolates the output fields of RRTM radiation from the reduced 
!! grid to the full grid.
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2010-12-03
!!
SUBROUTINE downscale_rad_output(p_patch, p_par_patch, p_par_int, p_par_grf,  &
  rg_aclcov, rg_lwflxclr, rg_lwflxall, rg_trsolclr, rg_trsolall,             &
  aclcov, lwflxclr, lwflxall, trsolclr, trsolall                             )


  ! Input types
  TYPE(t_patch),       TARGET, INTENT(IN)    ::  p_patch
  TYPE(t_patch),       TARGET, INTENT(IN)    ::  p_par_patch
  TYPE(t_int_state),   TARGET, INTENT(IN)    ::  p_par_int
  TYPE(t_gridref_state), TARGET, INTENT(IN)  ::  p_par_grf

  ! Other input fields (on reduced grid)
  REAL(wp), TARGET, INTENT(IN) ::                                               &
    rg_aclcov(:,:), rg_lwflxclr(:,:,:), rg_lwflxall(:,:,:), rg_trsolclr(:,:,:), &
    rg_trsolall(:,:,:)

  ! Corresponding output fields (on full grid)
  REAL(wp), INTENT(OUT) ::                                                &
    aclcov(:,:), lwflxclr(:,:,:), lwflxall(:,:,:), trsolclr(:,:,:), trsolall(:,:,:)

  ! Intermediate storage fields needed in the case of MPI parallelization
  REAL(wp), ALLOCATABLE, TARGET ::                                          &
    z_aclcov(:,:), z_lwflxclr(:,:,:), z_lwflxall(:,:,:), z_trsolclr(:,:,:), &
    z_trsolall(:,:,:)

  ! Pointers to output fields (no MPI) or intermediate fields (MPI)
  REAL(wp), POINTER ::                                                      &
    p_aclcov(:,:), p_lwflxclr(:,:,:), p_lwflxall(:,:,:), p_trsolclr(:,:,:), &
    p_trsolall(:,:,:)

  ! Additional storage fields to map 2D array(s) to 3D array
  REAL(wp), ALLOCATABLE :: zrg_aux3d(:,:,:), z_aux3d(:,:,:)

  ! Pointers to types needed to minimize code duplication for MPI/no-MPI cases
  TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
  TYPE(t_int_state),  POINTER     :: p_int => NULL()
  TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
  TYPE(t_patch),      POINTER     :: p_pp => NULL()

  ! Indices
  INTEGER :: jg, jgp, i_chidx, i_nchdom, nblks_c_lp,jk

  INTEGER :: nlev, nlevp1      !< number of full and half levels

  LOGICAL :: l_parallel, l_limit(5)
  REAL(wp) :: rlimval(5)
!-----------------------------------------------------------------------

  IF (msg_level >= 10) THEN
    WRITE(message_text,'(a,i2,a,i2)') 'Downscaling of radiation output fields',&
      p_par_patch%id,' =>',p_patch%id
    CALL message('downscale_rad_output',message_text)
  ENDIF

  IF (my_process_is_mpi_seq()) THEN
    l_parallel = .FALSE.
  ELSE
    l_parallel = .TRUE.
  ENDIF

  jgp = p_par_patch%id

  ! For the time being, the radiation grid is assumed to have the same levels
  ! as the full grid
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1

  IF (l_parallel) THEN
    jg = p_patch%id
    p_grf => p_grf_state_local_parent(jg)
    p_int => p_int_state_local_parent(jg)
    p_gcp => p_patch_local_parent(jg)%cells
    p_pp  => p_patch_local_parent(jg)
  ELSE
    p_grf => p_par_grf
    p_int => p_par_int
    p_gcp => p_par_patch%cells
    p_pp  => p_par_patch
  ENDIF

  i_nchdom = MAX(1,p_patch%n_childdom)
  i_chidx  = p_patch%parent_child_index

  nblks_c_lp = p_pp%nblks_c

  ! Allocation of local storage fields at local parent level in MPI-case
  IF (l_parallel .AND. jgp == 0) THEN

    ALLOCATE(z_aclcov(nproma,nblks_c_lp),          z_lwflxclr(nproma,nlevp1,nblks_c_lp), &
             z_lwflxall(nproma,nlevp1,nblks_c_lp), z_trsolclr(nproma,nlevp1,nblks_c_lp), &
             z_trsolall(nproma,nlevp1,nblks_c_lp)                                        )

  ENDIF

  ALLOCATE(zrg_aux3d(nproma,1,nblks_c_lp),z_aux3d(nproma,1,p_patch%nblks_c))

  ! Perform communication from parent to local parent grid in the MPI case,
  ! and set pointers such that further processing is the same for MPI / non-MPI cases
  IF (l_parallel .AND. jgp == 0) THEN

    CALL exchange_data(p_pp%comm_pat_glb_to_loc_c, RECV=z_aclcov, SEND=rg_aclcov  )

    CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, 4, 4*nlev+4, &
                            RECV1=z_lwflxclr, SEND1=rg_lwflxclr,     &
                            RECV2=z_lwflxall, SEND2=rg_lwflxall,     &
                            RECV3=z_trsolclr, SEND3=rg_trsolclr,     &
                            RECV4=z_trsolall, SEND4=rg_trsolall      )

     p_aclcov     => z_aclcov
     p_lwflxclr   => z_lwflxclr
     p_lwflxall   => z_lwflxall
     p_trsolclr   => z_trsolclr
     p_trsolall   => z_trsolall
   ELSE
     p_aclcov     => rg_aclcov
     p_lwflxclr   => rg_lwflxclr
     p_lwflxall   => rg_lwflxall
     p_trsolclr   => rg_trsolclr
     p_trsolall   => rg_trsolall
   ENDIF

   zrg_aux3d(:,1,:) = p_aclcov(:,:)

  ! Interpolate reduced-grid fields to full grid

  ! interpol_vec/scal_nudging needs all fields with full boundaries, so the arrays
  ! have to be sync'd before calling these routines.
  ! Please note that we cannot use sync_patch_array here (comparing parallel/non parallel results)
  ! since the arrays don't start with lower bound 1 in the non paralellel case!

  IF  (l_parallel) THEN

    CALL exchange_data_mult(p_pp%comm_pat_c, 5, 4*nlev+5, recv1=p_lwflxclr,      &
                            recv2=p_lwflxall, recv3=p_trsolclr, recv4=p_trsolall,&
                            recv5=zrg_aux3d                                      )

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

  CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, 0, 5, 1, &
    &                         f3din1=p_trsolall, f3dout1=trsolall,                 &
    &                         f3din2=p_trsolclr, f3dout2=trsolclr,                 &
    &                         f3din3=p_lwflxall, f3dout3=lwflxall,                 &
    &                         f3din4=p_lwflxclr, f3dout4=lwflxclr,                 &
    &                         f3din5=zrg_aux3d,  f3dout5=z_aux3d,                  &
                              llimit_nneg=l_limit, rlimval=rlimval                 )

   aclcov(:,:) = z_aux3d(:,1,:)

  IF (msg_level >= 14) THEN
    DO jk = 1, nlevp1
      WRITE(message_text,'(a,i3,2E10.3)') 'max/min SW rg transmissivity  = ',&
           & jk, MAXVAL (p_trsolall(:,jk,:)), MINVAL(p_trsolall(:,jk,:))
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i3,2E10.3)') 'max/min SW fg transmissivity  = ',&
           & jk, MAXVAL (trsolall(:,jk,:)), MINVAL(trsolall(:,jk,:))
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i3,2E10.3)') 'max/min LW rg net flux  = ',&
           & jk, MAXVAL (p_lwflxall(:,jk,:)), MINVAL(p_lwflxall(:,jk,:))
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i3,2E10.3)') 'max/min LW fg net flux  = ',&
           & jk, MAXVAL (lwflxall(:,jk,:)), MINVAL(lwflxall(:,jk,:))
      CALL message('', TRIM(message_text))
    ENDDO
  ENDIF

  IF (l_parallel .AND. jgp == 0) THEN
    DEALLOCATE(z_aclcov, z_lwflxclr, z_lwflxall, z_trsolclr, z_trsolall)
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
SUBROUTINE upscale_rad_input_rg(p_patch, p_par_patch, p_par_grf,        &
  &  cosmu0, albvisdir, alb_ther, temp_ifc, dpres_mc,                   &
  &  tot_cld, sqv, duco2, duo3,                                         &
  &  aeq1, aeq2, aeq3, aeq4, aeq5, pres_sfc,                            &
  &  rg_cosmu0, rg_albvisdir, rg_alb_ther, rg_temp_ifc, rg_dpres_mc,    &
  &  rg_tot_cld, rg_sqv, rg_duco2, rg_duo3,                             &
  &  rg_aeq1, rg_aeq2, rg_aeq3, rg_aeq4, rg_aeq5, rg_pres_sfc )


  ! Input types
  TYPE(t_patch),       TARGET, INTENT(IN)    ::  p_patch
  TYPE(t_patch),       TARGET, INTENT(IN)    ::  p_par_patch
  TYPE(t_gridref_state), TARGET, INTENT(IN)  ::  p_par_grf

  ! Other input fields (on full grid)
  REAL(wp), INTENT(IN) ::                                                             &
    & cosmu0(:,:), albvisdir(:,:), alb_ther(:,:), temp_ifc(:,:,:), &
    & dpres_mc(:,:,:), tot_cld(:,:,:,:), sqv(:,:,:), duco2(:,:,:), duo3(:,:,:),&
    & aeq1(:,:,:),aeq2(:,:,:),aeq3(:,:,:),aeq4(:,:,:),aeq5(:,:,:), &
    & pres_sfc(:,:)

  ! Corresponding output fields (on reduced grid)
  REAL(wp), TARGET, INTENT(OUT) ::                                           &
    & rg_cosmu0(:,:), rg_albvisdir(:,:), rg_alb_ther(:,:), rg_temp_ifc(:,:,:), &
    & rg_dpres_mc(:,:,:), rg_tot_cld(:,:,:,:), rg_sqv(:,:,:), rg_duco2(:,:,:), rg_duo3(:,:,:),&
    & rg_aeq1(:,:,:),rg_aeq2(:,:,:),rg_aeq3(:,:,:),rg_aeq4(:,:,:),rg_aeq5(:,:,:), &
    & rg_pres_sfc(:,:)

  ! Intermediate storage fields needed in the case of MPI parallelization
  REAL(wp), ALLOCATABLE, TARGET ::                                       &
    & z_cosmu0(:,:), z_albvisdir(:,:), z_alb_ther(:,:), z_temp_ifc(:,:,:), &
    & z_dpres_mc(:,:,:), z_tot_cld(:,:,:,:), z_sqv(:,:,:), z_duco2(:,:,:), z_duo3(:,:,:),&
    & z_aeq1(:,:,:),z_aeq2(:,:,:),z_aeq3(:,:,:),z_aeq4(:,:,:),z_aeq5(:,:,:), &
    & z_pres_sfc(:,:), z_aux3d(:,:,:),  zrg_aux3d(:,:,:)

  ! Pointers to output fields (no MPI) or intermediate fields (MPI)
  REAL(wp), POINTER ::                                                   &
    & p_cosmu0(:,:), p_albvisdir(:,:), p_alb_ther(:,:), p_temp_ifc(:,:,:), &
    & p_dpres_mc(:,:,:), p_tot_cld(:,:,:,:), p_sqv(:,:,:), p_duco2(:,:,:), p_duo3(:,:,:),&
    & p_aeq1(:,:,:),p_aeq2(:,:,:),p_aeq3(:,:,:),p_aeq4(:,:,:),p_aeq5(:,:,:), &
    & p_pres_sfc(:,:)


  ! Pointers to types needed to minimize code duplication for MPI/no-MPI cases
  TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
  TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
  TYPE(t_patch),      POINTER     :: p_pp => NULL()

  ! Indices
  INTEGER :: jb, jc, jk, jg, jgp, i_chidx, i_nchdom, &
             i_startblk, i_endblk, i_startidx, i_endidx, nblks_c_lp

  INTEGER :: nlev, nlevp1      !< number of full and half levels

  INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk
  LOGICAL :: l_parallel
  REAL(wp), POINTER :: p_fbkwgt(:,:,:)

  !-----------------------------------------------------------------------

  IF (msg_level >= 10) THEN
    WRITE(message_text,'(a,i2,a,i2)') 'Upscaling of radiation input fields',&
      p_patch%id,' =>',p_par_patch%id
    CALL message('upscale_rad_input',message_text)
  ENDIF

  IF (my_process_is_mpi_seq()) THEN
    l_parallel = .FALSE.
  ELSE
    l_parallel = .TRUE.
  ENDIF

  jgp = p_par_patch%id

  ! For the time being, the radiation grid is assumed to have the same levels
  ! as the full grid
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1

  IF (l_parallel) THEN
    jg = p_patch%id
    p_grf => p_grf_state_local_parent(jg)
    p_gcp => p_patch_local_parent(jg)%cells
    p_pp  => p_patch_local_parent(jg)
  ELSE
    p_grf => p_par_grf
    p_gcp => p_par_patch%cells
    p_pp  => p_par_patch
  ENDIF

  i_chidx  = p_patch%parent_child_index
  i_nchdom = MAX(1,p_par_patch%n_childdom)

  ! Set pointers to index and coefficient fields for cell-based variables
  iidx => p_gcp%child_idx
  iblk => p_gcp%child_blk

  p_fbkwgt => p_grf%fbk_wgt_c

  ! Allocation of local storage fields at local parent level in MPI-case
  IF (l_parallel .AND. jgp == 0) THEN
    nblks_c_lp = p_gcp%end_blk(min_rlcell,i_chidx)

    ALLOCATE(z_cosmu0(nproma,nblks_c_lp), z_albvisdir(nproma,nblks_c_lp),    &
             z_alb_ther(nproma,nblks_c_lp), z_temp_ifc(nproma,nlevp1,nblks_c_lp),    &
             z_dpres_mc(nproma,nlev,nblks_c_lp),z_tot_cld(nproma,nlev,nblks_c_lp,4), &
             z_sqv(nproma,nlev,nblks_c_lp),z_duco2(nproma,nlev,nblks_c_lp),  &
             z_duo3(nproma,nlev,nblks_c_lp), z_aeq1(nproma,nlev,nblks_c_lp), &
             z_aeq2(nproma,nlev,nblks_c_lp), z_aeq3(nproma,nlev,nblks_c_lp), &
             z_aeq4(nproma,nlev,nblks_c_lp), z_aeq5(nproma,nlev,nblks_c_lp), &
             z_pres_sfc(nproma,nblks_c_lp), z_aux3d(nproma,6,nblks_c_lp),    &
             zrg_aux3d(nproma,4,p_par_patch%nblks_c)                               )

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
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
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
      
      p_temp_ifc(jc,nlevp1,jb) =                                         &
        temp_ifc(iidx(jc,jb,1),nlevp1,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
        temp_ifc(iidx(jc,jb,2),nlevp1,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
        temp_ifc(iidx(jc,jb,3),nlevp1,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
        temp_ifc(iidx(jc,jb,4),nlevp1,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

    ENDDO

    IF (l_parallel .AND. jgp == 0) THEN ! combine 2D fields in a 3D field to speed up MPI communication
      DO jc = i_startidx, i_endidx
        z_aux3d(jc,1,jb) = p_cosmu0(jc,jb)
        z_aux3d(jc,2,jb) = p_albvisdir(jc,jb)
        z_aux3d(jc,3,jb) = p_alb_ther(jc,jb)
        z_aux3d(jc,4,jb) = p_pres_sfc(jc,jb)
      ENDDO
    ENDIF

    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx

        p_temp_ifc(jc,jk,jb) =                                         &
          temp_ifc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          temp_ifc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          temp_ifc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          temp_ifc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_dpres_mc(jc,jk,jb) =                                         &
          dpres_mc(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          dpres_mc(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          dpres_mc(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          dpres_mc(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_sqv(jc,jk,jb) =                                         &
          sqv(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          sqv(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          sqv(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          sqv(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_duco2(jc,jk,jb) =                                         &
          duco2(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          duco2(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          duco2(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          duco2(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_duo3(jc,jk,jb) =                                         &
          duo3(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          duo3(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          duo3(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          duo3(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_aeq1(jc,jk,jb) =                                         &
          aeq1(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          aeq1(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          aeq1(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          aeq1(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_aeq2(jc,jk,jb) =                                         &
          aeq2(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          aeq2(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          aeq2(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          aeq2(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)
 
        p_aeq3(jc,jk,jb) =                                         &
          aeq3(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          aeq3(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          aeq3(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          aeq3(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)

        p_aeq4(jc,jk,jb) =                                         &
          aeq4(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          aeq4(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          aeq4(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          aeq4(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)
        
        p_aeq5(jc,jk,jb) =                                         &
          aeq5(iidx(jc,jb,1),jk,iblk(jc,jb,1))*p_fbkwgt(jc,jb,1) + &
          aeq5(iidx(jc,jb,2),jk,iblk(jc,jb,2))*p_fbkwgt(jc,jb,2) + &
          aeq5(iidx(jc,jb,3),jk,iblk(jc,jb,3))*p_fbkwgt(jc,jb,3) + &
          aeq5(iidx(jc,jb,4),jk,iblk(jc,jb,4))*p_fbkwgt(jc,jb,4)
 
        
      ENDDO
    ENDDO

    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx

!CDIR EXPAND=4
        p_tot_cld(jc,jk,jb,1:4) =                                         &
          tot_cld(iidx(jc,jb,1),jk,iblk(jc,jb,1),1:4)*p_fbkwgt(jc,jb,1) + &
          tot_cld(iidx(jc,jb,2),jk,iblk(jc,jb,2),1:4)*p_fbkwgt(jc,jb,2) + &
          tot_cld(iidx(jc,jb,3),jk,iblk(jc,jb,3),1:4)*p_fbkwgt(jc,jb,3) + &
          tot_cld(iidx(jc,jb,4),jk,iblk(jc,jb,4),1:4)*p_fbkwgt(jc,jb,4)

      ENDDO
    ENDDO

  ENDDO
!$OMP END DO
!$OMP END PARALLEL


  IF (l_parallel .AND. jgp == 0) THEN

    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 6, 5*nlev+5,  &
                            RECV1=rg_temp_ifc, SEND1=z_temp_ifc,          &
                            RECV2=rg_dpres_mc, SEND2=z_dpres_mc,          &
                            RECV3=rg_sqv,      SEND3=z_sqv,               &
                            RECV4=rg_duco2,    SEND4=z_duco2,             &
                            RECV5=rg_duo3,     SEND5=z_duo3,              &
                            RECV6=zrg_aux3d,   SEND6=z_aux3d              )
    
    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 5, 5*nlev  ,  &
                            RECV1=rg_aeq1,     SEND1=z_aeq1,              &
                            RECV2=rg_aeq2,     SEND2=z_aeq2,              &
                            RECV3=rg_aeq3,     SEND3=z_aeq3,              &
                            RECV4=rg_aeq4,     SEND4=z_aeq4,              &
                            RECV5=rg_aeq5,     SEND5=z_aeq5               )
    

    CALL exchange_data_mult(p_pp%comm_pat_loc_to_glb_c_fbk, 4, 4*nlev,    &
                            RECV4D=rg_tot_cld, SEND4D=z_tot_cld           )


    i_startblk = p_par_patch%cells%start_blk(1,1)
    i_endblk   = p_par_patch%cells%end_blk(min_rlcell,i_nchdom)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_par_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, 1, min_rlcell, i_nchdom)

      DO jc = i_startidx, i_endidx
        rg_cosmu0(jc,jb)    = zrg_aux3d(jc,1,jb)
        rg_albvisdir(jc,jb) = zrg_aux3d(jc,2,jb)
        rg_alb_ther(jc,jb)  = zrg_aux3d(jc,3,jb)
        rg_pres_sfc(jc,jb)  = zrg_aux3d(jc,4,jb)
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    DEALLOCATE(z_cosmu0, z_albvisdir, z_alb_ther, z_temp_ifc, z_dpres_mc,z_tot_cld,&
      &  z_sqv, z_duco2, z_duo3, z_aeq1, z_aeq2, z_aeq3, z_aeq4, z_aeq5,     &
      &  z_pres_sfc, z_aux3d,zrg_aux3d )

  ENDIF

END SUBROUTINE upscale_rad_input_rg

SUBROUTINE downscale_rad_output_rg(p_patch, p_par_patch, p_par_int, p_par_grf,  &
  rg_lwflxall, rg_trsolall, lwflxall, trsolall                             )


  ! Input types
  TYPE(t_patch),       TARGET, INTENT(IN)    ::  p_patch
  TYPE(t_patch),       TARGET, INTENT(IN)    ::  p_par_patch
  TYPE(t_int_state),   TARGET, INTENT(IN)    ::  p_par_int
  TYPE(t_gridref_state), TARGET, INTENT(IN)  ::  p_par_grf

  ! Other input fields (on reduced grid)
  REAL(wp), TARGET, INTENT(IN) ::                                               &
    rg_lwflxall(:,:,:), rg_trsolall(:,:,:)

  ! Corresponding output fields (on full grid)
  REAL(wp), INTENT(OUT) ::                                                &
    lwflxall(:,:,:), trsolall(:,:,:)

  ! Intermediate storage fields needed in the case of MPI parallelization
  REAL(wp), ALLOCATABLE, TARGET ::                                          &
    z_lwflxall(:,:,:), z_trsolall(:,:,:)

  ! Pointers to output fields (no MPI) or intermediate fields (MPI)
  REAL(wp), POINTER ::                                                      &
    p_lwflxall(:,:,:), p_trsolall(:,:,:)

  ! Pointers to types needed to minimize code duplication for MPI/no-MPI cases
  TYPE(t_grid_cells), POINTER     :: p_gcp => NULL()
  TYPE(t_int_state),  POINTER     :: p_int => NULL()
  TYPE(t_gridref_state), POINTER  :: p_grf => NULL()
  TYPE(t_patch),      POINTER     :: p_pp => NULL()

  ! Indices
  INTEGER :: jg, jgp, i_chidx, i_nchdom, nblks_c_lp,jk

  INTEGER :: nlev, nlevp1      !< number of full and half levels

  LOGICAL :: l_parallel, l_limit(2)
  REAL(wp) :: rlimval(2)
!-----------------------------------------------------------------------

  IF (msg_level >= 10) THEN
    WRITE(message_text,'(a,i2,a,i2)') 'Downscaling of radiation output fields',&
      p_par_patch%id,' =>',p_patch%id
    CALL message('downscale_rad_output',message_text)
  ENDIF

  IF (my_process_is_mpi_seq()) THEN
    l_parallel = .FALSE.
  ELSE
    l_parallel = .TRUE.
  ENDIF

  jgp = p_par_patch%id

  ! For the time being, the radiation grid is assumed to have the same levels
  ! as the full grid
  nlev   = p_patch%nlev
  nlevp1 = p_patch%nlevp1

  IF (l_parallel) THEN
    jg = p_patch%id
    p_grf => p_grf_state_local_parent(jg)
    p_int => p_int_state_local_parent(jg)
    p_gcp => p_patch_local_parent(jg)%cells
    p_pp  => p_patch_local_parent(jg)
  ELSE
    p_grf => p_par_grf
    p_int => p_par_int
    p_gcp => p_par_patch%cells
    p_pp  => p_par_patch
  ENDIF

  i_nchdom = MAX(1,p_patch%n_childdom)
  i_chidx  = p_patch%parent_child_index

  nblks_c_lp = p_pp%nblks_c

  ! Allocation of local storage fields at local parent level in MPI-case
  IF (l_parallel .AND. jgp == 0) THEN

    ALLOCATE(z_lwflxall(nproma,nlevp1,nblks_c_lp), z_trsolall(nproma,nlevp1,nblks_c_lp) )

  ENDIF

  ! Perform communication from parent to local parent grid in the MPI case,
  ! and set pointers such that further processing is the same for MPI / non-MPI cases
  IF (l_parallel .AND. jgp == 0) THEN

    CALL exchange_data_mult(p_pp%comm_pat_glb_to_loc_c, 2, 2*nlev+2, &
                            RECV1=z_lwflxall, SEND1=rg_lwflxall,     &
                            RECV2=z_trsolall, SEND2=rg_trsolall      )

     p_lwflxall   => z_lwflxall
     p_trsolall   => z_trsolall
     
   ELSE

     p_lwflxall   => rg_lwflxall
     p_trsolall   => rg_trsolall
     
   ENDIF

  ! Interpolate reduced-grid fields to full grid

  ! interpol_vec/scal_nudging needs all fields with full boundaries, so the arrays
  ! have to be sync'd before calling these routines.
  ! Please note that we cannot use sync_patch_array here (comparing parallel/non parallel results)
  ! since the arrays don't start with lower bound 1 in the non paralellel case!

  IF  (l_parallel) THEN

    CALL exchange_data_mult(p_pp%comm_pat_c, 2, 2*nlev+2, recv1=p_lwflxall, recv2=p_trsolall )

  ENDIF

  IF (p_test_run) THEN
    trsolall = 0._wp
    lwflxall = 0._wp
  ENDIF

  l_limit(1)   = .TRUE.      ! limit transmissivity to positive values
  l_limit(2)   = .FALSE.
  rlimval(:)   = 2.94e-37_wp ! seems to be the lower threshold for SW transmissivity
                             ! in the RRTM scheme

  CALL interpol_scal_nudging (p_pp, p_int, p_grf%p_dom(i_chidx), i_chidx, 0, 2, 1, &
    &                         f3din1=p_trsolall, f3dout1=trsolall,                 &
    &                         f3din2=p_lwflxall, f3dout2=lwflxall,                 &
    &                         llimit_nneg=l_limit, rlimval=rlimval                 )

  IF (msg_level >= 14) THEN
    DO jk = 1, nlevp1
      WRITE(message_text,'(a,i3,2E10.3)') 'max/min SW rg transmissivity  = ',&
           & jk, MAXVAL (p_trsolall(:,jk,:)), MINVAL(p_trsolall(:,jk,:))
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i3,2E10.3)') 'max/min SW fg transmissivity  = ',&
           & jk, MAXVAL (trsolall(:,jk,:)), MINVAL(trsolall(:,jk,:))
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i3,2E10.3)') 'max/min LW rg net flux  = ',&
           & jk, MAXVAL (p_lwflxall(:,jk,:)), MINVAL(p_lwflxall(:,jk,:))
      CALL message('', TRIM(message_text))
      WRITE(message_text,'(a,i3,2E10.3)') 'max/min LW fg net flux  = ',&
           & jk, MAXVAL (lwflxall(:,jk,:)), MINVAL(lwflxall(:,jk,:))
      CALL message('', TRIM(message_text))
    ENDDO
  ENDIF

  IF (l_parallel .AND. jgp == 0) THEN
    DEALLOCATE(z_lwflxall, z_trsolall)
  ENDIF

END SUBROUTINE downscale_rad_output_rg



!>
!! This routine optimizes the boundary interpolation of diagnostic physics fields for output
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2010-12-03
!!
SUBROUTINE interpol_phys_grf (ptr_pp,ptr_pc,ptr_int, ptr_grf, jg, jgc, jn )

  USE mo_nwp_phy_state,      ONLY: prm_diag

  ! Input:
  TYPE(t_patch),                INTENT(in) :: ptr_pp
  TYPE(t_patch),                INTENT(in) :: ptr_pc
  TYPE(t_gridref_single_state), INTENT(in) :: ptr_grf
  TYPE(t_int_state),            INTENT(in) :: ptr_int
  INTEGER,                      INTENT(in) :: jg,jgc,jn

  ! Local fields
  INTEGER, PARAMETER  :: nfields=9    ! Number of 2D fields for which boundary interpolation is needed
  INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, jb, jc

  ! Temporary storage to do boundary interpolation for all 2D fields in one step
  REAL(wp) :: z_aux3d_p(nproma,nfields,ptr_pp%nblks_c), &
              z_aux3d_c(nproma,nfields,ptr_pc%nblks_c)

  IF (p_test_run) THEN
     z_aux3d_p(:,:,:) = 0._wp
  ENDIF

  i_startblk = ptr_pp%cells%start_blk(1,1)
  i_endblk   = ptr_pp%nblks_c

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_pp, jb, i_startblk, i_endblk, i_startidx, i_endidx, 1)



    DO jc = i_startidx, i_endidx
      z_aux3d_p(jc,1,jb) = prm_diag(jg)%tot_prec(jc,jb)
      z_aux3d_p(jc,2,jb) = prm_diag(jg)%rain_gsp(jc,jb)
      z_aux3d_p(jc,3,jb) = prm_diag(jg)%snow_gsp(jc,jb)
      z_aux3d_p(jc,4,jb) = prm_diag(jg)%rain_con(jc,jb)
      z_aux3d_p(jc,5,jb) = prm_diag(jg)%snow_con(jc,jb)
      z_aux3d_p(jc,6,jb) = prm_diag(jg)%tracer_rate(jc,jb,1)
      z_aux3d_p(jc,7,jb) = prm_diag(jg)%tracer_rate(jc,jb,2)
      z_aux3d_p(jc,8,jb) = prm_diag(jg)%tracer_rate(jc,jb,3)
      z_aux3d_p(jc,9,jb) = prm_diag(jg)%tracer_rate(jc,jb,4)

    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

    ! Halo update is needed before interpolation
    CALL sync_patch_array(SYNC_C,ptr_pp,z_aux3d_p)

    CALL interpol_scal_grf (ptr_pp, ptr_pc, ptr_int, ptr_grf, jn, 1, &
      &                     z_aux3d_p, z_aux3d_c, llimit_nneg=.TRUE.,&
      &                     lnoshift=.TRUE.)


  i_startblk = ptr_pc%cells%start_blk(1,1)
  i_endblk   = ptr_pc%cells%end_blk(grf_bdywidth_c,1)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc)
  DO jb = i_startblk, i_endblk

    CALL get_indices_c(ptr_pc, jb, i_startblk, i_endblk,        &
                       i_startidx, i_endidx, 1, grf_bdywidth_c)

    DO jc = i_startidx, i_endidx

      prm_diag(jgc)%tot_prec(jc,jb)       = z_aux3d_c(jc,1,jb)
      prm_diag(jgc)%rain_gsp(jc,jb)       = z_aux3d_c(jc,2,jb)
      prm_diag(jgc)%snow_gsp(jc,jb)       = z_aux3d_c(jc,3,jb)
      prm_diag(jgc)%rain_con(jc,jb)       = z_aux3d_c(jc,4,jb)
      prm_diag(jgc)%snow_con(jc,jb)       = z_aux3d_c(jc,5,jb)
      prm_diag(jgc)%tracer_rate(jc,jb,1)  = z_aux3d_c(jc,6,jb)
      prm_diag(jgc)%tracer_rate(jc,jb,2)  = z_aux3d_c(jc,7,jb)
      prm_diag(jgc)%tracer_rate(jc,jb,3)  = z_aux3d_c(jc,8,jb)
      prm_diag(jgc)%tracer_rate(jc,jb,4)  = z_aux3d_c(jc,9,jb)

    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE interpol_phys_grf

!>
!! This routine performs the feedback of diagnostic physics fields for output
!!
!! @par Revision History
!! Developed  by Guenther Zaengl, DWD, 2010-12-03
!!
SUBROUTINE feedback_phys_diag(p_patch, p_grf_state, jg, jgp)

  TYPE(t_patch),       TARGET, INTENT(IN)    ::  p_patch(n_dom_start:n_dom)
  TYPE(t_gridref_state), TARGET, INTENT(IN)  ::  p_grf_state(n_dom_start:n_dom)

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
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc)
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
!$OMP END DO
!$OMP END PARALLEL


  IF (l_parallel) THEN

    CALL exchange_data(p_pp%comm_pat_loc_to_glb_c_fbk, RECV=z_aux3d_par, SEND=z_aux3d_lp)
    p_aux3d => z_aux3d_par

  ENDIF

  i_startblk = p_patch(jgp)%cells%start_blk(grf_fbk_start_c,i_chidx)
  i_endblk   = p_patch(jgp)%cells%end_blk(min_rlcell_int,i_chidx)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc)
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
!$OMP END DO
!$OMP END PARALLEL


  DEALLOCATE(z_aux3d_lp, z_aux3d_par)

END SUBROUTINE feedback_phys_diag


END MODULE mo_phys_nest_utilities
