!>
!! This module contains the I/O routines for prep_icon and the driver
!! routines for topography preprocessing
!!
!! @author Guenther Zaengl, DWD
!!
!!
!! @par Revision History
!! First version by Guenther Zaengl, DWD (2011-07-13)
!!
!! @par Copyright
!! 2002-2010 by DWD and MPI-M
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

MODULE mo_prepicon_utils

  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: filename_max
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_run_config,          ONLY: msg_level, nvclev, iqv, iqc, iqi, iqr, iqs
  USE mo_extpar_config,       ONLY: extpar_filename,     &
    &                               generate_extpar_filename => generate_filename
  USE mo_dynamics_config,     ONLY: nnow, nnow_rcf, nnew, nnew_rcf
  USE mo_nonhydrostatic_config,ONLY: ivctype
  USE mo_sleve_config,        ONLY: lread_smt
  USE mo_nonhydro_types,      ONLY: t_nh_state
  USE mo_nwp_lnd_types,       ONLY: t_lnd_state
  USE mo_prepicon_types,      ONLY: t_prepicon_state, t_pi_atm_in, t_pi_sfc_in, &
    &                               t_pi_atm, t_pi_sfc
  USE mo_prepicon_config,     ONLY: i_oper_mode, nlev_in, l_w_in, nlevsoil_in, &
    &                               l_sfc_in, l_hice_in, l_sst_in,             &
    &                               ifs2icon_filename, generate_filename
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH, max_dom
  USE mo_physical_constants,  ONLY: tf_salt
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_grid_config,         ONLY: n_dom, nroot, global_cell_type
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_grf_intp_data_strc,  ONLY: t_gridref_state
  USE mo_mpi,                 ONLY: p_pe, p_io, p_bcast, p_comm_work_test, p_comm_work
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_smooth_topo,         ONLY: smooth_topography
  USE mo_util_netcdf,         ONLY: read_netcdf_data, read_netcdf_data_single, nf
  USE mo_model_domain,        ONLY: p_patch
  USE mo_io_config,           ONLY: lkeep_in_sync
  USE mo_io_util,             ONLY: gather_array1, gather_array2, outvar_desc,    &
    &                               GATHER_C, GATHER_E, GATHER_V, num_output_vars
  USE mo_datetime,            ONLY: t_datetime
  USE mo_nh_init_utils,       ONLY: nflat, nflatlev, compute_smooth_topo, init_vert_coord,  &
                                    hydro_adjust
  USE mo_nh_init_nest_utils,  ONLY: topography_blending, topography_feedback
  USE mo_grid_config,         ONLY: lfeedback
  USE mo_ifs_coord,           ONLY: alloc_vct, init_vct, vct, vct_a, vct_b
  USE mo_lnd_nwp_config,      ONLY: nlev_soil, ntiles_total
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_master_nml,          ONLY: model_base_dir
  USE mo_phyparam_soil,       ONLY: csalb_snow_min, csalb_snow_max,crhosmin_ml,crhosmax_ml
  USE mo_seaice_nwp,          ONLY: frsi_min

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  PRIVATE

  CHARACTER(len=*), PARAMETER :: version = '$Id$'


  TYPE(t_prepicon_state), ALLOCATABLE, TARGET :: prepicon(:) 


  CHARACTER(LEN=10) :: psvar 
  CHARACTER(LEN=10) :: geop_ml_var  ! model level surface geopotential
  CHARACTER(LEN=10) :: geop_sfc_var ! surface-level surface geopotential
  CHARACTER(LEN=10) :: alb_snow_var ! snow albedo

  PUBLIC :: prepicon


  PUBLIC :: init_prepicon
  PUBLIC :: compute_coord_fields
  PUBLIC :: deallocate_prepicon
  PUBLIC :: copy_prepicon2prog


  CONTAINS

  !-------------
  !>
  !! SUBROUTINE init_prepicon
  !! Initialization routine of prep_icon:
  !! Reads in data and processes topography blending and feedback in the
  !! presence of nested domains
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE init_prepicon (p_int_state, p_grf_state, prepicon, extdata)

    TYPE(t_int_state),     TARGET, INTENT(IN) :: p_int_state(:)
    TYPE(t_gridref_state), TARGET, INTENT(IN) :: p_grf_state(:)

    TYPE(t_prepicon_state), INTENT(INOUT) :: prepicon(:)
    TYPE(t_external_data),  INTENT(INOUT), OPTIONAL :: extdata(:)

    INTEGER :: jg, jlev, jk
    LOGICAL :: l_exist

    INTEGER :: no_cells, no_levels
    INTEGER :: ncid, dimid, varid, mpi_comm

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_prepicon_utils:init_prepicon'
    CHARACTER(LEN=filename_max) :: topo_file(max_dom), ifs2icon_file(max_dom)

!-------------------------------------------------------------------------


    ! Allocate memory for prep_icon state
    CALL allocate_prepicon (prepicon)


    ! Copy the already smoothed extdata topography fields to prepicon
    !
    DO jg = 1, n_dom
      prepicon(jg)%topography_c(:,:) = extdata(jg)%atm%topography_c(:,:)
      prepicon(jg)%topography_v(:,:) = extdata(jg)%atm%topography_v(:,:)
    ENDDO




    SELECT CASE(i_oper_mode)
!      CASE(2)
!        CALL init_from_dwd()
      CASE(3)
        !
        ! read horizontally interpolated IFS analysis for atmosphere
        ! 
        CALL read_ifs_atm( prepicon )
        !
        ! read horizontally interpolated IFS analysis for surface
        !
        IF ( l_sfc_in ) THEN
          CALL read_ifs_sfc( prepicon )
        ENDIF
    END SELECT





    IF (n_dom > 1) CALL topo_blending_and_fbk(p_int_state, p_grf_state, prepicon, 1)

    IF (PRESENT(extdata)) THEN

      ! Copy blended topography fields back to the external parameter state
      DO jg = 1, n_dom

        extdata(jg)%atm%topography_c(:,:) = prepicon(jg)%topography_c(:,:)
        extdata(jg)%atm%topography_v(:,:) = prepicon(jg)%topography_v(:,:)
        IF (i_oper_mode > 1 .AND. l_sfc_in) THEN
          ! In addition, copy climatological deep-soil temperature to soil level nlev_soil
          ! These are limited to -60 deg C because less is definitely nonsense
          prepicon(jg)%sfc%tsoil(:,:,nlev_soil) = MAX(213.15_wp,extdata(jg)%atm%t_cl(:,:))
        ENDIF
      ENDDO

    ENDIF

  END SUBROUTINE init_prepicon





  !>
  !! Read in horizontally interpolated IFS analysis (atmosphere only)
  !!
  !! Reads in horizontally interpolated IFS analysis atmosphere data 
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !! Modification by Daniel Reinert, DWD (2012-12-18)
  !! - encapsulate reading of IFS analysis
  !!
  SUBROUTINE read_ifs_atm (prepicon)

    TYPE(t_prepicon_state), INTENT(INOUT) :: prepicon(:)

    INTEGER :: jg, jlev, jk
    LOGICAL :: l_exist

    INTEGER :: no_cells, no_levels
    INTEGER :: ncid, dimid, varid, mpi_comm

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_prepicon_utils:read_ifs_atm'

    CHARACTER(LEN=filename_max) :: ifs2icon_file(max_dom)


    !-------------------------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!
    !! oper_mode = 3  !!
    !!!!!!!!!!!!!!!!!!!!

    DO jg = 1, n_dom

      jlev = p_patch(jg)%level


      ! Skip reading the atmospheric input data if a model domain 
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE


      ! Read in data from IFS2ICON
      !
      IF(p_pe == p_io ) THEN 
        !
        ! generate file name
        !
        ifs2icon_file(jg) = generate_filename(ifs2icon_filename, model_base_dir, &
          &                                   nroot, jlev, jg)
        INQUIRE (FILE=ifs2icon_file(jg), EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(TRIM(routine),'IFS2ICON file is not found: '//TRIM(ifs2icon_file(jg)))
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(ifs2icon_file(jg)), NF_NOWRITE, ncid), routine)

        !
        ! get number of cells
        !
        CALL nf(nf_inq_dimid(ncid, 'ncells', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)

        !
        ! get number of vertical levels
        !
        CALL nf(nf_inq_dimid(ncid, 'lev', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_levels), routine)

        !
        ! check the number of cells and vertical levels
        !
        IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in IFS2ICON file do not match.')
        ENDIF

        IF(nlev_in /= no_levels) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'nlev_in does not match the number of levels in IFS2ICON file.')
        ENDIF


        !
        ! Check if surface pressure (PS) or its logarithm (LNPS) is provided as input
        !
        IF (nf_inq_varid(ncid, 'PS', varid) == nf_noerr) THEN
          psvar = 'PS'
        ELSE IF (nf_inq_varid(ncid, 'LNPS', varid) == nf_noerr) THEN
          psvar = 'LNPS'
        ENDIF
        !
        ! Check if model-level surface Geopotential is provided as GEOSP or GEOP_ML
        !
        IF (nf_inq_varid(ncid, 'GEOSP', varid) == nf_noerr) THEN
          geop_ml_var = 'GEOSP'
        ELSE IF (nf_inq_varid(ncid, 'GEOP_ML', varid) == nf_noerr) THEN
          geop_ml_var = 'GEOP_ML'
        ELSE
          CALL finish(TRIM(routine),'Could not find model-level sfc geopotential')
        ENDIF

      ENDIF

      IF (msg_level >= 10) THEN
        WRITE(message_text,'(a)') 'surface pressure variable: '//TRIM(psvar)
        CALL message('', TRIM(message_text))
        WRITE(message_text,'(a)') 'Model-level surface geopotential: '//TRIM(geop_ml_var)
        CALL message('', TRIM(message_text))
      ENDIF



      ! start reading atmospheric fields
      !
      CALL read_netcdf_data_single (ncid, 'T', p_patch(jg)%n_patch_cells_g,                  &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     nlev_in,prepicon(jg)%atm_in%temp)

      CALL read_netcdf_data_single (ncid, 'U', p_patch(jg)%n_patch_cells_g,                  &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     nlev_in,prepicon(jg)%atm_in%u)

      CALL read_netcdf_data_single (ncid, 'V', p_patch(jg)%n_patch_cells_g,                  &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     nlev_in,prepicon(jg)%atm_in%v)

      IF (l_w_in) THEN ! note: input vertical velocity is in fact omega (Pa/s)
        CALL read_netcdf_data_single (ncid, 'W', p_patch(jg)%n_patch_cells_g,                &
          &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
          &                     nlev_in,prepicon(jg)%atm_in%omega)
      ENDIF

      CALL read_netcdf_data_single (ncid, 'QV', p_patch(jg)%n_patch_cells_g,                 &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     nlev_in,prepicon(jg)%atm_in%qv)

      CALL read_netcdf_data_single (ncid, 'QC', p_patch(jg)%n_patch_cells_g,                 &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     nlev_in,prepicon(jg)%atm_in%qc)

      CALL read_netcdf_data_single (ncid, 'QI', p_patch(jg)%n_patch_cells_g,                 &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     nlev_in,prepicon(jg)%atm_in%qi)

      CALL read_netcdf_data_single (ncid, 'QR', p_patch(jg)%n_patch_cells_g,                 &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     nlev_in,prepicon(jg)%atm_in%qr)

      CALL read_netcdf_data_single (ncid, 'QS', p_patch(jg)%n_patch_cells_g,                 &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     nlev_in,prepicon(jg)%atm_in%qs)

      CALL read_netcdf_data (ncid, TRIM(psvar), p_patch(jg)%n_patch_cells_g,          &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%atm_in%psfc)

      CALL read_netcdf_data (ncid, TRIM(geop_ml_var), p_patch(jg)%n_patch_cells_g,    &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%atm_in%phi_sfc)




      ! Allocate and read in vertical coordinate tables
      !
      IF (jg == 1) THEN

        ALLOCATE(vct_a(nlev_in+1), vct_b(nlev_in+1), vct(2*(nlev_in+1)))

        IF(p_pe == p_io) THEN
          CALL nf(nf_inq_varid(ncid, 'hyai', varid), routine)
          CALL nf(nf_get_var_double(ncid, varid, vct_a), routine)

          CALL nf(nf_inq_varid(ncid, 'hybi', varid), routine)
          CALL nf(nf_get_var_double(ncid, varid, vct_b), routine)
        ENDIF

        IF(p_test_run) THEN
          mpi_comm = p_comm_work_test 
        ELSE
          mpi_comm = p_comm_work
        ENDIF

        CALL p_bcast(vct_a, p_io, mpi_comm)
        CALL p_bcast(vct_b, p_io, mpi_comm)


        vct(1:nlev_in+1)             = vct_a(:)
        vct(nlev_in+2:2*(nlev_in+1)) = vct_b(:)

        nvclev = 2*(nlev_in+1)

        CALL alloc_vct(nlev_in)
        CALL init_vct(nlev_in)

        IF (msg_level >= 15) THEN
          WRITE(message_text,'(a)') 'vct table values of input data:'
          CALL message('', TRIM(message_text))

          DO jk = 1, nlev_in+1
            WRITE(message_text,'(a,i4,F12.4,F12.8)') 'jk, vct_a, vct_b: ',jk, vct_a(jk), vct_b(jk)
            CALL message('', TRIM(message_text))
          ENDDO
        ENDIF

      ENDIF  ! jg=1


      ! close file
      !
      IF(p_pe == p_io) CALL nf(nf_close(ncid), routine)

    ENDDO ! loop over model domains



  END SUBROUTINE read_ifs_atm





  !>
  !! Read in horizontally interpolated IFS analysis (surface only)
  !!
  !! Reads in horizontally interpolated IFS analysis surface data
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !! Modification by Daniel Reinert, DWD (2012-12-18)
  !! - encapsulate reading of IFS analysis
  !!
  SUBROUTINE read_ifs_sfc (prepicon)

    TYPE(t_prepicon_state), INTENT(INOUT) :: prepicon(:)

    INTEGER :: jg, jlev
    LOGICAL :: l_exist

    INTEGER :: no_cells, no_levels
    INTEGER :: ncid, dimid, varid, mpi_comm

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_prepicon_utils:read_ifs_sfc'

    CHARACTER(LEN=filename_max) :: ifs2icon_file(max_dom)
    LOGICAL :: l_sst_present = .FALSE.     !TRUE if SST is present in the IFS input file


    !-------------------------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!
    !! oper_mode = 3  !!
    !!!!!!!!!!!!!!!!!!!!

    DO jg = 1, n_dom

      jlev = p_patch(jg)%level


      ! Skip reading the atmospheric input data if a model domain 
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE


      ! Read in data from IFS2ICON
      !
      IF(p_pe == p_io ) THEN 
        !
        ! generate file name
        !
        ifs2icon_file(jg) = generate_filename(ifs2icon_filename, model_base_dir, &
          &                                   nroot, jlev, jg)
        INQUIRE (FILE=ifs2icon_file(jg), EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(TRIM(routine),'IFS2ICON file is not found: '//TRIM(ifs2icon_file(jg)))
        ENDIF

        !
        ! open file
        !
        CALL nf(nf_open(TRIM(ifs2icon_file(jg)), NF_NOWRITE, ncid), routine)

        !
        ! get number of cells
        !
        CALL nf(nf_inq_dimid(ncid, 'ncells', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)

        !
        ! get number of vertical levels
        !
        CALL nf(nf_inq_dimid(ncid, 'lev', dimid), routine)
        CALL nf(nf_inq_dimlen(ncid, dimid, no_levels), routine)

        !
        ! check the number of cells and vertical levels
        !
        IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'Number of patch cells and cells in IFS2ICON file do not match.')
        ENDIF

        IF(nlev_in /= no_levels) THEN
          CALL finish(TRIM(ROUTINE),&
          & 'nlev_in does not match the number of levels in IFS2ICON file.')
        ENDIF




        ! Check, if the surface-level surface geopotential (GEOP_SFC) is available. 
        ! If GEOP_SFC is missing, a warning will be issued and the model-level surface 
        ! geopotential (GEOSP or GEOP_ML) will be used instead.
        IF (nf_inq_varid(ncid, 'GEOP_SFC', varid) == nf_noerr) THEN
          geop_sfc_var = 'GEOP_SFC'
        ELSE

          WRITE (message_text,'(a,a)')                            &
            &  'surface-level surface geopotential is missing. ', &
            &  'use model-level surface geopotential, instead.'
          CALL message(TRIM(routine),TRIM(message_text))

          ! use model level geopotential instead
          geop_sfc_var = geop_ml_var
        ENDIF

        ! Check, if the snow albedo ALB_SNOW is available. 
        ! If ALB_SNOW is missing, a warning will be issued and RHO_SNOW 
        ! will be used instead to determine FRESHSNOW.
        IF (nf_inq_varid(ncid, 'ALB_SNOW', varid) == nf_noerr) THEN
          WRITE (message_text,'(a,a)')                            &
            &  'snow albedo available, ', &
            &  'used to determine freshsnow.'
          alb_snow_var = 'ALB_SNOW'
        ELSE

          WRITE (message_text,'(a,a)')                            &
            &  'snow albedo is missing. ', &
            &  'use snow density value, instead.'
          CALL message(TRIM(routine),TRIM(message_text))

          alb_snow_var = 'RHO_SNOW'
        ENDIF


        ! Check, if sea-ice thickness field is provided as input
        ! IF H_ICE is missing, set l_hice_in=.FALSE.
        IF (nf_inq_varid(ncid, 'H_ICE', varid) == nf_noerr) THEN
          WRITE (message_text,'(a,a)')                            &
            &  'sea-ice thickness available'
          l_hice_in = .TRUE.
        ELSE

          WRITE (message_text,'(a,a)')                            &
            &  'sea-ice thickness not available. ', &
            &  'initialize with constant value (0.5 m), instead.'
          CALL message(TRIM(routine),TRIM(message_text))

          l_hice_in = .FALSE.
        ENDIF

        ! Check, if sea surface temperature field is provided as input
        ! IF SST is missing, set l_sst_in=.FALSE.
        IF (nf_inq_varid(ncid, 'SST', varid) == nf_noerr) THEN
          WRITE (message_text,'(a,a)')                            &
            &  'sea surface temperature available'
          l_sst_present = .TRUE.

        ELSE

          WRITE (message_text,'(a,a)')                            &
            &  'sea surface temperature not available. ', &
            &  'initialize with skin temperature, instead.'
          CALL message(TRIM(routine),TRIM(message_text))
          l_sst_present = .FALSE.
          l_sst_in = .FALSE.     !it has to be set to FALSE
        ENDIF

      ENDIF  ! p_io



      IF(p_test_run) THEN
        mpi_comm = p_comm_work_test 
      ELSE
        mpi_comm = p_comm_work
      ENDIF

      CALL p_bcast(l_hice_in, p_io, mpi_comm)

      CALL p_bcast(l_sst_in, p_io, mpi_comm)

      CALL p_bcast(l_sst_present, p_io, mpi_comm)

      CALL p_bcast(alb_snow_var, p_io, mpi_comm)


      ! start reading surface fields
      !
      CALL read_netcdf_data (ncid, TRIM(geop_sfc_var), p_patch(jg)%n_patch_cells_g,   &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%phi)

      CALL read_netcdf_data (ncid, 'SKT', p_patch(jg)%n_patch_cells_g,                &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%tskin)
      IF ( l_sst_present) THEN
       CALL read_netcdf_data (ncid, 'SST', p_patch(jg)%n_patch_cells_g,                &
         &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
         &                     prepicon(jg)%sfc_in%sst)
      ELSE 
       prepicon(jg)%sfc_in%sst(:,:)=0.0_wp
      END IF

      CALL read_netcdf_data (ncid, 'T_SNOW', p_patch(jg)%n_patch_cells_g,             &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%tsnow)

      CALL read_netcdf_data (ncid,TRIM(alb_snow_var), p_patch(jg)%n_patch_cells_g,    &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%snowalb)
 
      CALL read_netcdf_data (ncid, 'W_SNOW', p_patch(jg)%n_patch_cells_g,             &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%snowweq)

      CALL read_netcdf_data (ncid,'RHO_SNOW', p_patch(jg)%n_patch_cells_g,            &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%snowdens)

      CALL read_netcdf_data (ncid, 'W_I', p_patch(jg)%n_patch_cells_g,                &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%skinres)

      CALL read_netcdf_data (ncid, 'LSM', p_patch(jg)%n_patch_cells_g,                &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%ls_mask)

      CALL read_netcdf_data (ncid, 'CI', p_patch(jg)%n_patch_cells_g,                 &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%seaice)

      CALL read_netcdf_data (ncid, 'STL1', p_patch(jg)%n_patch_cells_g,               &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%tsoil(:,:,1))

      CALL read_netcdf_data (ncid, 'STL2', p_patch(jg)%n_patch_cells_g,               &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%tsoil(:,:,2))

      CALL read_netcdf_data (ncid, 'STL3', p_patch(jg)%n_patch_cells_g,               &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%tsoil(:,:,3))

      CALL read_netcdf_data (ncid, 'STL4', p_patch(jg)%n_patch_cells_g,               &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%tsoil(:,:,4))

      CALL read_netcdf_data (ncid, 'SMIL1', p_patch(jg)%n_patch_cells_g,              &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%wsoil(:,:,1))

      CALL read_netcdf_data (ncid, 'SMIL2', p_patch(jg)%n_patch_cells_g,              &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%wsoil(:,:,2))

      CALL read_netcdf_data (ncid, 'SMIL3', p_patch(jg)%n_patch_cells_g,              &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%wsoil(:,:,3))

      CALL read_netcdf_data (ncid, 'SMIL4', p_patch(jg)%n_patch_cells_g,              &
        &                     p_patch(jg)%n_patch_cells, p_patch(jg)%cells%glb_index, &
        &                     prepicon(jg)%sfc_in%wsoil(:,:,4))



      ! close file
      !
      IF(p_pe == p_io) CALL nf(nf_close(ncid), routine)

    ENDDO ! loop over model domains



  END SUBROUTINE read_ifs_sfc







  RECURSIVE SUBROUTINE topo_blending_and_fbk(p_int, p_grf, prepicon, jg)

    TYPE(t_int_state),     TARGET, INTENT(IN) :: p_int(:)
    TYPE(t_gridref_state), TARGET, INTENT(IN) :: p_grf(:)

    TYPE(t_prepicon_state), INTENT(INOUT) :: prepicon(:)

    INTEGER, INTENT(IN) :: jg

    INTEGER :: jgc, jn


    ! Loop over nested domains
    DO jn = 1, p_patch(jg)%n_childdom

      jgc = p_patch(jg)%child_id(jn)


      CALL topography_blending(p_patch(jg), p_patch(jgc), p_int(jg),  &
               p_int(jgc), p_grf(jg)%p_dom(jn), jn,                   &
               prepicon(jg)%topography_c, prepicon(jgc)%topography_c, &
               prepicon(jgc)%topography_v                             )

      IF (p_patch(jgc)%n_childdom > 0) &
        CALL topo_blending_and_fbk(p_int, p_grf, prepicon, jgc)

    ENDDO

    DO jn = 1, p_patch(jg)%n_childdom

      jgc = p_patch(jg)%child_id(jn)

      IF (lfeedback(jgc)) THEN
        CALL topography_feedback(p_patch(jg), p_int(jg), p_grf(jg), jn, &
          prepicon(jg)%topography_c, prepicon(jgc)%topography_c,        &
          prepicon(jg)%topography_v                                     )
      ENDIF

    ENDDO


  END SUBROUTINE topo_blending_and_fbk


  !-------------
  !>
  !! SUBROUTINE copy_prepicon2prog
  !! Copies atmospheric and surface fields interpolated by prep_icon to the
  !! prognostic model state variables if the respective parts of prep_icon 
  !! are called directly 
  !!
  !! Required input: prepicon state
  !! Output is written on fields of NH state and land state
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-28)
  !!
  !!
  SUBROUTINE copy_prepicon2prog(prepicon, p_nh_state, p_lnd_state, ext_data)

    TYPE(t_prepicon_state), INTENT(IN) :: prepicon(:)

    TYPE(t_nh_state),      INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_lnd_state),     INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_external_data), INTENT(   IN) :: ext_data(:)

    INTEGER :: jg, jb, jk, jc, je, jt, js, jp, ic
    INTEGER :: nblks_c, npromz_c, nblks_e, npromz_e, nlen, nlev, nlevp1, ntl, ntlr

!$OMP PARALLEL PRIVATE(jg,nblks_c,npromz_c,nblks_e,npromz_e,nlev,nlevp1,ntl,ntlr)
    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      nblks_c   = p_patch(jg)%nblks_c
      npromz_c  = p_patch(jg)%npromz_c
      nblks_e   = p_patch(jg)%nblks_e
      npromz_e  = p_patch(jg)%npromz_e
      nlev      = p_patch(jg)%nlev
      nlevp1    = p_patch(jg)%nlevp1
      ntl       = nnow(jg)
      ntlr      = nnow_rcf(jg)

!$OMP DO PRIVATE(jb,jk,je,nlen) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_e

        IF (jb /= nblks_e) THEN
          nlen = nproma
        ELSE
          nlen = npromz_e
        ENDIF

        ! Wind speed
        DO jk = 1, nlev
          DO je = 1, nlen
            p_nh_state(jg)%prog(ntl)%vn(je,jk,jb) = prepicon(jg)%atm%vn(je,jk,jb)
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO

!$OMP DO PRIVATE(jb,jk,jc,nlen,jt,js,jp,ic) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks_c

        IF (jb /= nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = npromz_c
        ENDIF

        ! 3D fields
        DO jk = 1, nlev
          DO jc = 1, nlen
            ! Dynamic prognostic variables on cell points
            p_nh_state(jg)%prog(ntl)%w(jc,jk,jb)       = prepicon(jg)%atm%w(jc,jk,jb)
            p_nh_state(jg)%prog(ntl)%theta_v(jc,jk,jb) = prepicon(jg)%atm%theta_v(jc,jk,jb)
            p_nh_state(jg)%prog(ntl)%exner(jc,jk,jb)   = prepicon(jg)%atm%exner(jc,jk,jb)
            p_nh_state(jg)%prog(ntl)%rho(jc,jk,jb)     = prepicon(jg)%atm%rho(jc,jk,jb)

            p_nh_state(jg)%prog(ntl)%rhotheta_v(jc,jk,jb) = &
              p_nh_state(jg)%prog(ntl)%rho(jc,jk,jb) * p_nh_state(jg)%prog(ntl)%theta_v(jc,jk,jb)

            ! Moisture variables
            p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqv) = prepicon(jg)%atm%qv(jc,jk,jb)
            p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqc) = prepicon(jg)%atm%qc(jc,jk,jb)
            p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqi) = prepicon(jg)%atm%qi(jc,jk,jb)
            p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqr) = prepicon(jg)%atm%qr(jc,jk,jb)
            p_nh_state(jg)%prog(ntlr)%tracer(jc,jk,jb,iqs) = prepicon(jg)%atm%qs(jc,jk,jb)
          ENDDO
        ENDDO

        ! w at surface level
        DO jc = 1, nlen
          p_nh_state(jg)%prog(ntl)%w(jc,nlevp1,jb)      = prepicon(jg)%atm%w(jc,nlevp1,jb)
          p_nh_state(jg)%prog(nnew(jg))%w(jc,nlevp1,jb) = prepicon(jg)%atm%w(jc,nlevp1,jb)
       ENDDO

        ! ground temperature
        IF (l_sfc_in) THEN
          DO jc = 1, nlen
            p_lnd_state(jg)%prog_lnd(ntlr)%t_g(jc,jb)         = prepicon(jg)%sfc%tskin(jc,jb)
            p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_g(jc,jb) = prepicon(jg)%sfc%tskin(jc,jb)
            p_lnd_state(jg)%diag_lnd%t_skin(jc,jb)            = prepicon(jg)%sfc%tskin(jc,jb)
          ENDDO
          ! Fill also SST and sea ice fraction fields over ocean points; SST is limited to 30 deg C
          ! Note: missing values of the sea ice fraction, which may occur due to differing land-sea masks, 
          ! are indicated with -999.9; non-ocean points are filled with zero for both fields
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, ext_data(jg)%atm%sp_count(jb)
            jc = ext_data(jg)%atm%idx_lst_sp(ic,jb)
            IF ( l_sst_in .AND. prepicon(jg)%sfc%sst(jc,jb) > 10._wp  ) THEN
              p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = prepicon(jg)%sfc%sst(jc,jb)              
            ELSE
             p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = MIN(303.15_wp,prepicon(jg)%sfc%tskin(jc,jb))
            ENDIF
            !
            ! In case of missing sea ice fraction values, we make use of the sea 
            ! surface temperature (tskin over ocean points). For tskin<=tf_salt, 
            ! we set the sea ice fraction to one. For tskin>tf_salt, we set it to 0.
            ! Note: tf_salt=271.45K is the salt-water freezing point
            !
            IF ( prepicon(jg)%sfc%seaice(jc,jb) > -999.0_wp ) THEN
              p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = prepicon(jg)%sfc%seaice(jc,jb) 
            ELSE    ! missing value
              IF ( prepicon(jg)%sfc%tskin(jc,jb) <= tf_salt ) THEN
                p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = 1._wp     ! sea ice point
              ELSE
                p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = 0._wp     ! water point
              ENDIF
            ENDIF

            ! For fr_seaice in ]0,frsi_min[, set fr_seaice to 0
            ! For fr_seaice in ]1-frsi_min,1[, set fr_seaice to 1. This will ensure in 
            ! init_sea_lists, that sea-ice and water fractions sum up exactly to the total 
            ! sea fraction.
            IF (p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) < frsi_min ) THEN
               p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = 0._wp
            ENDIF
            IF (p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) > (1._wp-frsi_min) ) THEN
               p_lnd_state(jg)%diag_lnd%fr_seaice(jc,jb) = 1._wp
            ENDIF

          ENDDO
          ! In addition, write skin temperature to lake points, limited to 33 deg C. These will
          ! be used to initialize lake points until something more reasonable becomes available
!CDIR NODEP,VOVERTAKE,VOB
          DO ic = 1, ext_data(jg)%atm%fp_count(jb)
            jc = ext_data(jg)%atm%idx_lst_fp(ic,jb)
            p_lnd_state(jg)%diag_lnd%t_seasfc(jc,jb) = MIN(306.15_wp,prepicon(jg)%sfc%tskin(jc,jb))
          ENDDO
        ENDIF

        IF ( l_sfc_in .AND. atm_phy_nwp_config(jg)%inwp_surface > 0 ) THEN
          DO jt = 1, ntiles_total
            DO jc = 1, nlen
               p_lnd_state(jg)%prog_lnd(ntlr)%t_snow_t(jc,jb,jt)           = &
                &                                                prepicon(jg)%sfc%tsnow   (jc,jb)

               ! Initialize freshsnow
               ! for seapoints, freshsnow is set to 0
               IF(alb_snow_var == 'ALB_SNOW') THEN
              p_lnd_state(jg)%diag_lnd%freshsnow_t(jc,jb,jt)    =  MAX(0._wp,MIN(1._wp,   &
            &                           (prepicon(jg)%sfc%snowalb (jc,jb)-csalb_snow_min) &
            &                          /(csalb_snow_max-csalb_snow_min)))                 &
            &                          * REAL(NINT(ext_data(jg)%atm%fr_land(jc,jb)),wp) 
              ELSE
              p_lnd_state(jg)%diag_lnd%freshsnow_t(jc,jb,jt)    =  MAX(0._wp,MIN(1._wp,   &
            &                     1._wp - ((prepicon(jg)%sfc%snowalb (jc,jb)-crhosmin_ml) &
            &                    /(crhosmax_ml-crhosmin_ml))))                            &
            &                    * REAL(NINT(ext_data(jg)%atm%fr_land(jc,jb)),wp)
               END IF

              p_lnd_state(jg)%prog_lnd(ntlr)%w_snow_t(jc,jb,jt)           = &
                &                                                prepicon(jg)%sfc%snowweq (jc,jb)
              p_lnd_state(jg)%prog_lnd(ntlr)%rho_snow_t(jc,jb,jt)         = &
                &                                                prepicon(jg)%sfc%snowdens(jc,jb) 
              p_lnd_state(jg)%prog_lnd(ntlr)%w_i_t(jc,jb,jt)              = &
                &                                                prepicon(jg)%sfc%skinres (jc,jb)
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_snow_t(jc,jb,jt)   = &
                &                                                prepicon(jg)%sfc%tsnow   (jc,jb)
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%w_snow_t(jc,jb,jt)   = &
                &                                                prepicon(jg)%sfc%snowweq (jc,jb)
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%rho_snow_t(jc,jb,jt) = &
                &                                                prepicon(jg)%sfc%snowdens(jc,jb)
              p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%w_i_t(jc,jb,jt)      = &
                &                                                prepicon(jg)%sfc%skinres (jc,jb)
            ENDDO
          ENDDO

          ! Multi-layer surface fields
          DO jt = 1, ntiles_total

            DO js = 0, nlev_soil
              jp = js+1 ! indexing for the ICON state field starts at 1
              DO jc = 1, nlen
                p_lnd_state(jg)%prog_lnd(ntlr)%t_so_t(jc,jp,jb,jt)        = &
                  &                                              prepicon(jg)%sfc%tsoil(jc,jb,js)
                p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%t_so_t(jc,jp,jb,jt)= &
                  &                                              prepicon(jg)%sfc%tsoil(jc,jb,js)
              ENDDO
            ENDDO

            ! For soil water, no comparable layer shift exists
            DO js = 1, nlev_soil
              DO jc = 1, nlen
                p_lnd_state(jg)%prog_lnd(ntlr)%w_so_t(jc,js,jb,jt) = &
                  &                                              prepicon(jg)%sfc%wsoil(jc,jb,js)
                p_lnd_state(jg)%prog_lnd(nnew_rcf(jg))%w_so_t(jc,js,jb,jt)= &
                  &                                              prepicon(jg)%sfc%wsoil(jc,jb,js)
              ENDDO
            ENDDO

          ENDDO
        ENDIF

      ENDDO
!$OMP END DO NOWAIT

    ENDDO
!$OMP END PARALLEL

    ! Finally, compute exact hydrostatic adjustment for thermodynamic fields
    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE
      ntl = nnow(jg)

      CALL hydro_adjust(p_patch(jg), p_nh_state(jg)%metrics,                                  &
                        p_nh_state(jg)%prog(ntl)%rho,     p_nh_state(jg)%prog(ntl)%exner,     &
                        p_nh_state(jg)%prog(ntl)%theta_v, p_nh_state(jg)%prog(ntl)%rhotheta_v )

    ENDDO

  END SUBROUTINE copy_prepicon2prog




  SUBROUTINE compute_coord_fields(p_int, prepicon)

    TYPE(t_int_state),  INTENT(IN)       :: p_int(:)
    TYPE(t_prepicon_state), INTENT(INOUT):: prepicon(:)

    INTEGER :: jg, jgp, nblks_c, npromz_c, nlev, nlevp1
    INTEGER :: i_nchdom
    INTEGER :: nshift_total(n_dom)       !< Total shift of model top w.r.t. global domain
    LOGICAL :: l_half_lev_centr

    !------------------------------------------------------------------------

    SELECT CASE (global_cell_type)
    CASE (6)
      l_half_lev_centr = .TRUE.
      ! The HALF LEVEL where the model layer are flat, moves one layer upward.
      ! there could also be a zero there
      nflat = nflat-1
    CASE DEFAULT
      l_half_lev_centr = .FALSE.
    END SELECT

    nshift_total(1) = 0

    DO jg = 1,n_dom

      nblks_c   = p_patch(jg)%nblks_c
      npromz_c  = p_patch(jg)%npromz_c

      i_nchdom   = MAX(1,p_patch(jg)%n_childdom)

      ! number of vertical levels
      nlev   = p_patch(jg)%nlev
      nlevp1 = p_patch(jg)%nlevp1

      ! total shift of model top with respect to global domain
      IF (jg > 1) THEN
        jgp = p_patch(jg)%parent_id
        nshift_total(jg) = nshift_total(jgp) + p_patch(jg)%nshift
        nflatlev(jg)     = nflatlev(1) - nshift_total(jg)
      ENDIF
      IF (global_cell_type == 6) nflatlev(jg) = nflat
      IF (jg > 1 .AND. nshift_total(jg) > 0 .AND. nflatlev(jg) < 1) THEN
        CALL finish ('mo_prepicon_utils:compute_coord_fields', &
                     'nflat must be higher than the top of the innermost nested domain')
      ENDIF

      ! Compute smooth topography when SLEVE coordinate is used
      IF (ivctype == 2  .AND. .NOT. lread_smt ) THEN
        CALL compute_smooth_topo(p_patch(jg), p_int(jg),            &
          prepicon(jg)%topography_c, prepicon(jg)%topography_c_smt, &
          prepicon(jg)%topography_v, prepicon(jg)%topography_v_smt  )
      ENDIF

      ! Compute 3D coordinate fields for cell points (for vertex points needed
      ! only temporarily in set_nh_metrics)
      CALL init_vert_coord(prepicon(jg)%topography_c, prepicon(jg)%topography_c_smt, &
                           prepicon(jg)%z_ifc, prepicon(jg)%z_mc,                    &
                           nlev, nblks_c, npromz_c, nshift_total(jg), nflatlev(jg),  &
                           l_half_lev_centr)


    ENDDO


  END SUBROUTINE compute_coord_fields




  !-------------
  !>
  !! SUBROUTINE allocate_prepicon
  !! Allocates the components of the prepicon data type
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE allocate_prepicon (prepicon)

    TYPE(t_prepicon_state), INTENT(INOUT) :: prepicon(:)

    ! Local variables: loop control and dimensions
    INTEGER :: jg, nlev, nlevp1, nblks_c, nblks_v, nblks_e

!-------------------------------------------------------------------------

    ! Loop over model domains
    DO jg = 1, n_dom

      nlev = p_patch(jg)%nlev
      nlevp1 = nlev + 1
      nblks_c = p_patch(jg)%nblks_c
      nblks_v = p_patch(jg)%nblks_v
      nblks_e = p_patch(jg)%nblks_e

      IF (i_oper_mode == 2) nlev_in = p_patch(jg)%nlev

      ! basic prep_icon data
      ALLOCATE(prepicon(jg)%topography_c    (nproma,nblks_c),        &
               prepicon(jg)%topography_c_smt(nproma,nblks_c) ,       &
               prepicon(jg)%topography_v    (nproma,nblks_v),        &
               prepicon(jg)%topography_v_smt(nproma,nblks_v) ,       &
               prepicon(jg)%z_ifc           (nproma,nlevp1,nblks_c), &
               prepicon(jg)%z_mc            (nproma,nlev  ,nblks_c) )

      IF (i_oper_mode >= 2) THEN
        ! Allocate atmospheric input data
        ALLOCATE(prepicon(jg)%atm_in%psfc(nproma,         nblks_c), &
                 prepicon(jg)%atm_in%phi_sfc(nproma,      nblks_c), &
                 prepicon(jg)%atm_in%pres (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%z3d  (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%temp (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%u    (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%v    (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%w    (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%omega(nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%qv   (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%qc   (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%qi   (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%qr   (nproma,nlev_in,nblks_c), &
                 prepicon(jg)%atm_in%qs   (nproma,nlev_in,nblks_c)  )

        ! Allocate surface input data
        ! The extra soil temperature levels are not read in; they are only used to simplify vertical interpolation
        ALLOCATE(prepicon(jg)%sfc_in%phi      (nproma,nblks_c                ), &
                 prepicon(jg)%sfc_in%tskin    (nproma,nblks_c                ), &
                 prepicon(jg)%sfc_in%sst      (nproma,nblks_c                ), &
                 prepicon(jg)%sfc_in%tsnow    (nproma,nblks_c                ), &
                 prepicon(jg)%sfc_in%snowalb  (nproma,nblks_c                ), &
                 prepicon(jg)%sfc_in%snowweq  (nproma,nblks_c                ), &
                 prepicon(jg)%sfc_in%snowdens (nproma,nblks_c                ), &
                 prepicon(jg)%sfc_in%skinres  (nproma,nblks_c                ), &
                 prepicon(jg)%sfc_in%ls_mask  (nproma,nblks_c                ), &
                 prepicon(jg)%sfc_in%seaice   (nproma,nblks_c                ), &
                 prepicon(jg)%sfc_in%tsoil    (nproma,nblks_c,0:nlevsoil_in+1), &
                 prepicon(jg)%sfc_in%wsoil    (nproma,nblks_c,0:nlevsoil_in+1)  )

        ! Allocate atmospheric output data
        ALLOCATE(prepicon(jg)%atm%vn        (nproma,nlev  ,nblks_e), &
                 prepicon(jg)%atm%u         (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%v         (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%w         (nproma,nlevp1,nblks_c), &
                 prepicon(jg)%atm%temp      (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%exner     (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%pres      (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%rho       (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%theta_v   (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%qv        (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%qc        (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%qi        (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%qr        (nproma,nlev  ,nblks_c), &
                 prepicon(jg)%atm%qs        (nproma,nlev  ,nblks_c)  )

        ! Allocate surface output data
        ALLOCATE(prepicon(jg)%sfc%tskin    (nproma,nblks_c             ), &
                 prepicon(jg)%sfc%sst      (nproma,nblks_c             ), &
                 prepicon(jg)%sfc%tsnow    (nproma,nblks_c             ), &
                 prepicon(jg)%sfc%snowalb  (nproma,nblks_c             ), &
                 prepicon(jg)%sfc%snowweq  (nproma,nblks_c             ), &
                 prepicon(jg)%sfc%snowdens (nproma,nblks_c             ), &
                 prepicon(jg)%sfc%skinres  (nproma,nblks_c             ), &
                 prepicon(jg)%sfc%ls_mask  (nproma,nblks_c             ), &
                 prepicon(jg)%sfc%seaice   (nproma,nblks_c             ), &
                 prepicon(jg)%sfc%tsoil    (nproma,nblks_c,0:nlev_soil ), &
                 prepicon(jg)%sfc%wsoil    (nproma,nblks_c,nlev_soil)     )


      ENDIF

    ENDDO ! loop over model domains

  END SUBROUTINE allocate_prepicon

  !-------------
  !>
  !! SUBROUTINE deallocate_prepicon
  !! Deallocates the components of the prepicon data type
  !!
  !! @par Revision History
  !! Initial version by Guenther Zaengl, DWD(2011-07-14)
  !!
  !!
  SUBROUTINE deallocate_prepicon (prepicon)

    TYPE(t_prepicon_state), INTENT(INOUT) :: prepicon(:)

    ! Local variables: loop control
    INTEGER :: jg

!-------------------------------------------------------------------------

    ! Loop over model domains
    DO jg = 1, n_dom

      ! basic prep_icon data
      DEALLOCATE(prepicon(jg)%topography_c,     &
                 prepicon(jg)%topography_c_smt, &
                 prepicon(jg)%topography_v,     &
                 prepicon(jg)%topography_v_smt, &
                 prepicon(jg)%z_ifc,            &
                 prepicon(jg)%z_mc              )

      IF (i_oper_mode >= 2) THEN
        ! atmospheric input data
        DEALLOCATE(prepicon(jg)%atm_in%psfc,    &
                   prepicon(jg)%atm_in%phi_sfc, &
                   prepicon(jg)%atm_in%pres,    &
                   prepicon(jg)%atm_in%z3d,     &
                   prepicon(jg)%atm_in%temp,    &
                   prepicon(jg)%atm_in%u,       &
                   prepicon(jg)%atm_in%v,       &
                   prepicon(jg)%atm_in%w,       &
                   prepicon(jg)%atm_in%omega,   &
                   prepicon(jg)%atm_in%qv,      &
                   prepicon(jg)%atm_in%qc,      &
                   prepicon(jg)%atm_in%qi,      &
                   prepicon(jg)%atm_in%qr,      &
                   prepicon(jg)%atm_in%qs )

        ! surface input data
        DEALLOCATE(prepicon(jg)%sfc_in%phi,      &
                   prepicon(jg)%sfc_in%tskin,    &
                   prepicon(jg)%sfc_in%sst,    &
                   prepicon(jg)%sfc_in%tsnow,    &
                   prepicon(jg)%sfc_in%snowalb,  &
                   prepicon(jg)%sfc_in%snowweq,  &
                   prepicon(jg)%sfc_in%snowdens, &
                   prepicon(jg)%sfc_in%skinres,  &
                   prepicon(jg)%sfc_in%ls_mask,  &
                   prepicon(jg)%sfc_in%seaice,   &
                   prepicon(jg)%sfc_in%tsoil,    &
                   prepicon(jg)%sfc_in%wsoil     )



        ! atmospheric output data
        DEALLOCATE(prepicon(jg)%atm%vn,      &
                   prepicon(jg)%atm%u,       &
                   prepicon(jg)%atm%v,       &
                   prepicon(jg)%atm%w,       &
                   prepicon(jg)%atm%temp,    &
                   prepicon(jg)%atm%exner,   &
                   prepicon(jg)%atm%pres,    &  
                   prepicon(jg)%atm%rho,     &
                   prepicon(jg)%atm%theta_v, &
                   prepicon(jg)%atm%qv,      &
                   prepicon(jg)%atm%qc,      &
                   prepicon(jg)%atm%qi,      &
                   prepicon(jg)%atm%qr,      &
                   prepicon(jg)%atm%qs       )

        ! surface output data
        DEALLOCATE(prepicon(jg)%sfc%tskin,    &
                   prepicon(jg)%sfc%sst,      &
                   prepicon(jg)%sfc%tsnow,    &
                   prepicon(jg)%sfc%snowalb,  &
                   prepicon(jg)%sfc%snowweq,  &
                   prepicon(jg)%sfc%snowdens, &
                   prepicon(jg)%sfc%skinres,  &
                   prepicon(jg)%sfc%ls_mask,  &
                   prepicon(jg)%sfc%seaice,   &
                   prepicon(jg)%sfc%tsoil,    &
                   prepicon(jg)%sfc%wsoil     )


      ENDIF


    ENDDO ! loop over model domains

  END SUBROUTINE deallocate_prepicon

END MODULE mo_prepicon_utils

