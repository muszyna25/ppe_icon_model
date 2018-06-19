!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/model-distribution.html and the
!! file COPYING in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the headers of the routines.
!!
MODULE mo_psrad_communication
#ifndef NOMPI
#ifdef HAVE_YAXT
#define USE_PSRAD_COMMUNICATION
#endif
#endif


  USE iso_c_binding,               ONLY: c_ptr
  USE mo_exception,                ONLY: finish, message


  USE mo_kind,                     ONLY: wp
  USE mo_master_control,           ONLY: atmo_process, ps_radiation_process
  USE mo_parallel_config,          ONLY: nproma
  USE mo_model_domain,             ONLY: t_patch, p_patch
  USE mo_grid_config,              ONLY: n_dom, n_dom_start
  USE mo_run_config,               ONLY: num_lev

  USE mo_psrad_general,            ONLY: nbndsw

#ifdef USE_PSRAD_COMMUNICATION
  USE mpi
  USE mo_mpi,                      ONLY: get_mpi_work_intercomm, p_pe_work, &
                                         get_my_global_mpi_communicator
  USE yaxt,                        ONLY: xt_initialize, xt_initialized, &
                                         xt_redist_msg, xt_redist, &
                                         xt_redist_single_array_base_new, &
                                         xt_redist_collection_new, &
                                         xt_redist_delete, xt_redist_s_exchange
#endif

  IMPLICIT NONE

  PRIVATE

  INTEGER :: psrad_atm_intercomm
#ifdef USE_PSRAD_COMMUNICATION
  TYPE(xt_redist), ALLOCATABLE :: exchange_redist_atmo_2_psrad(:)
  TYPE(xt_redist), ALLOCATABLE :: exchange_redist_psrad_2_atmo(:)
#endif

  PUBLIC :: setup_atmo_2_psrad_communication, &
            setup_psrad_2_atmo_communication, &
            exchange_data_atmo_2_psrad, &
            exchange_data_psrad_2_atmo, &
            free_atmo_psrad_communication

CONTAINS

#ifdef USE_PSRAD_COMMUNICATION
  SUBROUTINE generate_redists(is_psrad, p_patch, no_of_levels, &
                              exchange_redist_atmo_2_psrad, &
                              exchange_redist_psrad_2_atmo)

    LOGICAL, INTENT(IN)             :: is_psrad
    TYPE(t_patch), INTENT(IN)       :: p_patch
    INTEGER, INTENT(IN)             :: no_of_levels

    TYPE(xt_redist), INTENT(OUT) :: exchange_redist_atmo_2_psrad
    TYPE(xt_redist), INTENT(OUT) :: exchange_redist_psrad_2_atmo

    INTEGER :: alloc_cell_blocks
    INTEGER :: shape2d(2), shape3d(3), shape3d_layer_interfaces(3)
    INTEGER :: dt_int_2d, dt_logical_2d, dt_wp_nbndsw, dt_wp_2d, dt_wp_3d, &
               dt_wp_3d_cfc, dt_wp_3d_layer

    TYPE(xt_redist_msg) :: redist_msg_int(1)
    TYPE(xt_redist_msg) :: redist_msg_int_2d(1)
    TYPE(xt_redist_msg) :: redist_msg_logical_2d(1)
    TYPE(xt_redist_msg) :: redist_msg_wp(1)
    TYPE(xt_redist_msg) :: redist_msg_wp_nbndsw(1)
    TYPE(xt_redist_msg) :: redist_msg_wp_2d(1)
    TYPE(xt_redist_msg) :: redist_msg_wp_3d(1)
    TYPE(xt_redist_msg) :: redist_msg_wp_3d_cfc(1)
    TYPE(xt_redist_msg) :: redist_msg_wp_3d_layer(1)

    TYPE(xt_redist) :: atmo_2_psrad_redist_int
    TYPE(xt_redist) :: atmo_2_psrad_redist_int_2d
    TYPE(xt_redist) :: atmo_2_psrad_redist_logical_2d
    TYPE(xt_redist) :: atmo_2_psrad_redist_wp
    TYPE(xt_redist) :: atmo_2_psrad_redist_wp_nbndsw
    TYPE(xt_redist) :: atmo_2_psrad_redist_wp_2d
    TYPE(xt_redist) :: atmo_2_psrad_redist_wp_3d
    TYPE(xt_redist) :: atmo_2_psrad_redist_wp_3d_cfc
    TYPE(xt_redist) :: atmo_2_psrad_redist_wp_3d_layer
    TYPE(xt_redist) :: psrad_2_atmo_redist_int_2d
    TYPE(xt_redist) :: psrad_2_atmo_redist_logical_2d
    TYPE(xt_redist) :: psrad_2_atmo_redist_wp_2d
    TYPE(xt_redist) :: psrad_2_atmo_redist_wp_3d
    TYPE(xt_redist) :: psrad_2_atmo_redist_wp_3d_layer

    INTEGER :: atmo_2_psrad_send_count, atmo_2_psrad_recv_count
    INTEGER :: psrad_2_atmo_send_count, psrad_2_atmo_recv_count

    INTEGER :: ierr

    alloc_cell_blocks = p_patch%alloc_cell_blocks

    ! shapes of data arrays
    shape2d  = (/nproma,       alloc_cell_blocks/)
    shape3d  = (/nproma, no_of_levels, alloc_cell_blocks/)
    shape3d_layer_interfaces = (/nproma,no_of_levels+1,alloc_cell_blocks/)

    ! MPI data types of a single field
    CALL MPI_Type_contiguous(shape2d(1) * shape2d(2), MPI_INTEGER, &
                             dt_int_2d, ierr)
    CALL MPI_Type_contiguous(shape2d(1) * shape2d(2), MPI_LOGICAL, &
                             dt_logical_2d, ierr)
    CALL MPI_Type_contiguous(nbndsw, MPI_DOUBLE_PRECISION, &
                             dt_wp_nbndsw, ierr)
    CALL MPI_Type_contiguous(shape2d(1) * shape2d(2), MPI_DOUBLE_PRECISION, &
                             dt_wp_2d, ierr)
    CALL MPI_Type_contiguous(shape3d(1) * shape3d(2) * shape3d(3), &
                             MPI_DOUBLE_PRECISION, dt_wp_3d, ierr)
    CALL MPI_Type_contiguous(2 * shape3d(1) * shape3d(2) * shape3d(3), &
                             MPI_DOUBLE_PRECISION, dt_wp_3d_cfc, ierr)
    CALL MPI_Type_contiguous(shape3d_layer_interfaces(1) * &
                             shape3d_layer_interfaces(2) * &
                             shape3d_layer_interfaces(3), &
                             MPI_DOUBLE_PRECISION, dt_wp_3d_layer, ierr)

    redist_msg_int(1)%rank = p_pe_work
    redist_msg_int(1)%datatype = MPI_INTEGER
    redist_msg_int_2d(1)%rank = p_pe_work
    redist_msg_int_2d(1)%datatype = dt_int_2d
    redist_msg_logical_2d(1)%rank = p_pe_work
    redist_msg_logical_2d(1)%datatype = dt_logical_2d
    redist_msg_wp(1)%rank = p_pe_work
    redist_msg_wp(1)%datatype = MPI_DOUBLE_PRECISION
    redist_msg_wp_nbndsw(1)%rank = p_pe_work
    redist_msg_wp_nbndsw(1)%datatype = dt_wp_nbndsw
    redist_msg_wp_2d(1)%rank = p_pe_work
    redist_msg_wp_2d(1)%datatype = dt_wp_2d
    redist_msg_wp_3d(1)%rank = p_pe_work
    redist_msg_wp_3d(1)%datatype = dt_wp_3d
    redist_msg_wp_3d_cfc(1)%rank = p_pe_work
    redist_msg_wp_3d_cfc(1)%datatype = dt_wp_3d_cfc
    redist_msg_wp_3d_layer(1)%rank = p_pe_work
    redist_msg_wp_3d_layer(1)%datatype = dt_wp_3d_layer

    atmo_2_psrad_send_count = MERGE(0, 1, is_psrad)
    atmo_2_psrad_recv_count = MERGE(1, 0, is_psrad)
    psrad_2_atmo_send_count = MERGE(1, 0, is_psrad)
    psrad_2_atmo_recv_count = MERGE(0, 1, is_psrad)

    atmo_2_psrad_redist_int = &
      xt_redist_single_array_base_new( &
        atmo_2_psrad_send_count, atmo_2_psrad_recv_count, &
        redist_msg_int, redist_msg_int, psrad_atm_intercomm)
    atmo_2_psrad_redist_int_2d = &
      xt_redist_single_array_base_new( &
        atmo_2_psrad_send_count, atmo_2_psrad_recv_count, &
        redist_msg_int_2d, redist_msg_int_2d, psrad_atm_intercomm)
    psrad_2_atmo_redist_int_2d = &
      xt_redist_single_array_base_new( &
        psrad_2_atmo_send_count, psrad_2_atmo_recv_count, &
        redist_msg_int_2d, redist_msg_int_2d, psrad_atm_intercomm)
    atmo_2_psrad_redist_logical_2d = &
      xt_redist_single_array_base_new( &
        atmo_2_psrad_send_count, atmo_2_psrad_recv_count, &
        redist_msg_logical_2d, redist_msg_logical_2d, psrad_atm_intercomm)
    psrad_2_atmo_redist_logical_2d = &
      xt_redist_single_array_base_new( &
        psrad_2_atmo_send_count, psrad_2_atmo_recv_count, &
        redist_msg_logical_2d, redist_msg_logical_2d, psrad_atm_intercomm)
    atmo_2_psrad_redist_wp = &
      xt_redist_single_array_base_new( &
        atmo_2_psrad_send_count, atmo_2_psrad_recv_count, &
        redist_msg_wp, redist_msg_wp, psrad_atm_intercomm)
    atmo_2_psrad_redist_wp_nbndsw = &
      xt_redist_single_array_base_new( &
        atmo_2_psrad_send_count, atmo_2_psrad_recv_count, &
        redist_msg_wp_nbndsw, redist_msg_wp_nbndsw, psrad_atm_intercomm)
    atmo_2_psrad_redist_wp_2d = &
      xt_redist_single_array_base_new( &
        atmo_2_psrad_send_count, atmo_2_psrad_recv_count, &
        redist_msg_wp_2d, redist_msg_wp_2d, psrad_atm_intercomm)
    psrad_2_atmo_redist_wp_2d = &
      xt_redist_single_array_base_new( &
        psrad_2_atmo_send_count, psrad_2_atmo_recv_count, &
        redist_msg_wp_2d, redist_msg_wp_2d, psrad_atm_intercomm)
    atmo_2_psrad_redist_wp_3d = &
      xt_redist_single_array_base_new( &
        atmo_2_psrad_send_count, atmo_2_psrad_recv_count, &
        redist_msg_wp_3d, redist_msg_wp_3d, psrad_atm_intercomm)
    psrad_2_atmo_redist_wp_3d = &
      xt_redist_single_array_base_new( &
        psrad_2_atmo_send_count, psrad_2_atmo_recv_count, &
        redist_msg_wp_3d, redist_msg_wp_3d, psrad_atm_intercomm)
    atmo_2_psrad_redist_wp_3d_cfc = &
      xt_redist_single_array_base_new( &
        atmo_2_psrad_send_count, atmo_2_psrad_recv_count, &
        redist_msg_wp_3d_cfc, redist_msg_wp_3d_cfc, &
        psrad_atm_intercomm)
    atmo_2_psrad_redist_wp_3d_layer = &
      xt_redist_single_array_base_new( &
        atmo_2_psrad_send_count, atmo_2_psrad_recv_count, &
        redist_msg_wp_3d_layer, redist_msg_wp_3d_layer, &
        psrad_atm_intercomm)
    psrad_2_atmo_redist_wp_3d_layer = &
      xt_redist_single_array_base_new( &
        psrad_2_atmo_send_count, psrad_2_atmo_recv_count, &
        redist_msg_wp_3d_layer, redist_msg_wp_3d_layer, &
        psrad_atm_intercomm)

    exchange_redist_atmo_2_psrad = &
      xt_redist_collection_new((/atmo_2_psrad_redist_int, & !< aerosol control
                                 atmo_2_psrad_redist_int, & !< number of levels
                                 atmo_2_psrad_redist_int_2d, & ! convection_type
                                 atmo_2_psrad_redist_wp, & !< orbit and time dependent solar constant for radiation time step
                                 atmo_2_psrad_redist_wp_nbndsw, & !< fraction of TSI in the 14 RRTM SW bands
                                 atmo_2_psrad_redist_logical_2d, & ! loland
                                 atmo_2_psrad_redist_logical_2d, & ! loglac
                                 atmo_2_psrad_redist_wp_2d, & ! pcos_mu0
                                 atmo_2_psrad_redist_wp_2d, & ! daylght_frc
                                 atmo_2_psrad_redist_wp_2d, & ! albvisdir
                                 atmo_2_psrad_redist_wp_2d, & ! albnirdir
                                 atmo_2_psrad_redist_wp_2d, & ! albvisdif
                                 atmo_2_psrad_redist_wp_2d, & ! albnirdif
                                 atmo_2_psrad_redist_wp_3d, & ! zf
                                 atmo_2_psrad_redist_wp_3d_layer, & ! zh
                                 atmo_2_psrad_redist_wp_3d, & ! dz
                                 atmo_2_psrad_redist_wp_2d, & ! surface_pressure
                                 atmo_2_psrad_redist_wp_3d, & ! pressure
                                 atmo_2_psrad_redist_wp_2d, & ! surface_temperature
                                 atmo_2_psrad_redist_wp_3d, & ! temperature_fl
                                 atmo_2_psrad_redist_wp_3d_layer, & ! temperature_hl
                                 atmo_2_psrad_redist_wp_3d, & ! dry_air_mass
                                 atmo_2_psrad_redist_wp_3d, & ! water_vapor
                                 atmo_2_psrad_redist_wp_3d, & ! cloud_water
                                 atmo_2_psrad_redist_wp_3d, & ! cloud_ice
                                 atmo_2_psrad_redist_wp_3d, & ! cloud_nuclei
                                 atmo_2_psrad_redist_wp_3d, & ! fractional_cloud_cover
                                 atmo_2_psrad_redist_wp_3d, & ! co2
                                 atmo_2_psrad_redist_wp_3d, & ! ch4
                                 atmo_2_psrad_redist_wp_3d, & ! n2o
                                 atmo_2_psrad_redist_wp_3d_cfc, & ! cfc
                                 atmo_2_psrad_redist_wp_3d, & ! o3
                                 atmo_2_psrad_redist_wp_3d/), & ! o2
                               psrad_atm_intercomm)

    exchange_redist_psrad_2_atmo = &
      xt_redist_collection_new((/psrad_2_atmo_redist_wp_3d_layer, & ! Clear-sky downwelling_longwave_flux_in_air
                                 psrad_2_atmo_redist_wp_3d_layer, & ! Clear-sky upwelling_longwave_flux_in_air
                                 psrad_2_atmo_redist_wp_3d_layer, & ! Clear-sky downwelling_shortwave_flux_in_air_assuming_clear_sky
                                 psrad_2_atmo_redist_wp_3d_layer, & ! Clear-sky upwelling_shortwave_flux_in_air
                                 psrad_2_atmo_redist_wp_3d_layer, & ! All-sky downwelling_longwave_flux_in_air
                                 psrad_2_atmo_redist_wp_3d_layer, & ! All-sky upwelling_longwave_flux_in_air
                                 psrad_2_atmo_redist_wp_3d_layer, & ! All-sky downwelling_shortwave_flux_in_air_assuming_clear_sky
                                 psrad_2_atmo_redist_wp_3d_layer, & ! All-sky upwelling_shortwave_flux_in_air
                                 psrad_2_atmo_redist_wp_2d, & ! surface_downwelling_direct_visible_flux_in_air_at_rad_time
                                 psrad_2_atmo_redist_wp_2d, & ! surface_downwelling_direct_par_flux_in_air_at_rad_time
                                 psrad_2_atmo_redist_wp_2d, & ! surface_downwelling_direct_nearir_flux_in_air_at_rad_time
                                 psrad_2_atmo_redist_wp_2d, & ! surface_downwelling_diffuse_visible_flux_in_air_at_rad_time
                                 psrad_2_atmo_redist_wp_2d, & ! surface_downwelling_diffuse_par_flux_in_air_at_rad_time
                                 psrad_2_atmo_redist_wp_2d, & ! surface_downwelling_diffuse_nearir_flux_in_air_at_rad_time
                                 psrad_2_atmo_redist_wp_2d, & ! surface_upwelling_visible_flux_in_air_at_rad_time
                                 psrad_2_atmo_redist_wp_2d, & ! surface_upwelling_par_flux_in_air_at_rad_time
                                 psrad_2_atmo_redist_wp_2d/), & ! surface_upwelling_nearir_flux_in_air_at_rad_time
                               psrad_atm_intercomm)

    CALL xt_redist_delete(atmo_2_psrad_redist_int)
    CALL xt_redist_delete(atmo_2_psrad_redist_int_2d)
    CALL xt_redist_delete(atmo_2_psrad_redist_logical_2d)
    CALL xt_redist_delete(atmo_2_psrad_redist_wp)
    CALL xt_redist_delete(atmo_2_psrad_redist_wp_nbndsw)
    CALL xt_redist_delete(atmo_2_psrad_redist_wp_2d)
    CALL xt_redist_delete(atmo_2_psrad_redist_wp_3d)
    CALL xt_redist_delete(atmo_2_psrad_redist_wp_3d_cfc)
    CALL xt_redist_delete(atmo_2_psrad_redist_wp_3d_layer)
    CALL xt_redist_delete(psrad_2_atmo_redist_int_2d)
    CALL xt_redist_delete(psrad_2_atmo_redist_logical_2d)
    CALL xt_redist_delete(psrad_2_atmo_redist_wp_2d)
    CALL xt_redist_delete(psrad_2_atmo_redist_wp_3d)
    CALL xt_redist_delete(psrad_2_atmo_redist_wp_3d_layer)

    ! clean up
    CALL MPI_Type_free(dt_wp_3d_layer, ierr)
    CALL MPI_Type_free(dt_wp_3d_cfc, ierr)
    CALL MPI_Type_free(dt_wp_3d, ierr)
    CALL MPI_Type_free(dt_wp_2d, ierr)
    CALL MPI_Type_free(dt_wp_nbndsw, ierr)
    CALL MPI_Type_free(dt_logical_2d, ierr)
    CALL MPI_Type_free(dt_int_2d, ierr)

  END SUBROUTINE generate_redists
#endif

  SUBROUTINE setup_communication(is_psrad)

    LOGICAL, INTENT(IN) :: is_psrad

#ifdef USE_PSRAD_COMMUNICATION
    INTEGER :: i

    IF (.NOT. xt_initialized()) &
      CALL xt_initialize(get_my_global_mpi_communicator())

    ALLOCATE(exchange_redist_atmo_2_psrad(n_dom_start:n_dom), &
             exchange_redist_psrad_2_atmo(n_dom_start:n_dom))

    DO i = n_dom_start, n_dom

      CALL generate_redists(is_psrad, p_patch(i), num_lev(i), &
                            exchange_redist_atmo_2_psrad(i), &
                            exchange_redist_psrad_2_atmo(i))
    END DO
#endif

  END SUBROUTINE setup_communication

  SUBROUTINE setup_atmo_2_psrad_communication()

#ifdef USE_PSRAD_COMMUNICATION
    psrad_atm_intercomm = get_mpi_work_intercomm(ps_radiation_process)

    CALL setup_communication(.FALSE.)
#endif

  END SUBROUTINE setup_atmo_2_psrad_communication

  SUBROUTINE setup_psrad_2_atmo_communication()

#ifdef USE_PSRAD_COMMUNICATION
    psrad_atm_intercomm = get_mpi_work_intercomm(atmo_process)

    CALL setup_communication(.TRUE.)
#endif

  END SUBROUTINE setup_psrad_2_atmo_communication

  SUBROUTINE free_atmo_psrad_communication
 
   INTEGER :: ierr

#ifdef USE_PSRAD_COMMUNICATION
   INTEGER :: i

    DO i = n_dom_start, n_dom

      CALL xt_redist_delete(exchange_redist_atmo_2_psrad(i))
      CALL xt_redist_delete(exchange_redist_psrad_2_atmo(i))

    END DO

    DEALLOCATE(exchange_redist_atmo_2_psrad)
    DEALLOCATE(exchange_redist_psrad_2_atmo)

    CALL MPI_COMM_FREE(psrad_atm_intercomm, ierr)
#endif

  END SUBROUTINE free_atmo_psrad_communication

  SUBROUTINE exchange_data_atmo_2_psrad(idom, &
                                        irad_aero, &
                                        klev, &
                                        ktype, &
                                        psctm, &
                                        ssi_factor, &
                                        loland, &
                                        loglac, &
                                        pcos_mu0,     &
                                        daylght_frc,  &
                                        alb_vis_dir,  &
                                        alb_nir_dir,  &
                                        alb_vis_dif,  &
                                        alb_nir_dif,  &
                                        zf,         &
                                        zh,         &
                                        dz,         &
                                        pp_sfc,       &
                                        pp_fl,      &
                                        tk_sfc,       &
                                        tk_fl,      &
                                        tk_hl,      &
                                        xm_dry,     &
                                        xm_vap,     &
                                        xm_liq,     &
                                        xm_ice,     &
                                        cdnc,       &
                                        xc_frc, &
                                        xm_co2, &
                                        xm_ch4, &
                                        xm_n2o, &
                                        xm_cfc, &
                                        xm_o3, &
                                        xm_o2)

    INTEGER, INTENT(IN) :: idom

    !INTEGER,INTENT(INOUT)  :: &
    TYPE(c_ptr), INTENT(IN) :: &
         irad_aero,        & !< aerosol control
         klev,             & !< number of levels
         ktype               !< type of convection

    !REAL(wp),INTENT(INOUT) :: &
    TYPE(c_ptr), INTENT(IN) :: &
        psctm,                & !< orbit and time dependent solar constant for radiation time step
        ssi_factor              !< fraction of TSI in the 14 RRTM SW bands

    !LOGICAL,INTENT(INOUT) ::   &
    TYPE(c_ptr), INTENT(IN) ::  &
         loland,                & !< land sea mask, land=.true.
         loglac                   !< glacier mask, glacier=.true.

    !REAL(WP),INTENT(INOUT)  :: &
    TYPE(c_ptr), INTENT(IN)  :: &
         pcos_mu0,     & !< mu0 for solar zenith angle
         daylght_frc,  & !< daylight fraction; with diurnal cycle 0 or 1, with zonal mean in [0,1]
         alb_vis_dir,  & !< surface albedo for vis range and dir light
         alb_nir_dir,  & !< surface albedo for NIR range and dir light
         alb_vis_dif,  & !< surface albedo for vis range and dif light
         alb_nir_dif,  & !< surface albedo for NIR range and dif light
         zf,         & !< geometric height at full level in m
         zh,         & !< geometric height at half level in m
         dz,         & !< geometric height thickness in m
         pp_sfc,       & !< surface pressure in Pa
         pp_fl,      & !< full level pressure in Pa
         tk_sfc,       & !< surface temperature in K
         tk_fl,      & !< full level temperature in K
         tk_hl,      & !< half level temperature in K
         xm_dry,     & !< dry air     mass in kg/m2
         xm_vap,     & !< water vapor mass in kg/m2
         xm_liq,     & !< cloud water mass in kg/m2
         xm_ice,     & !< cloud ice   mass in kg/m2
         cdnc,       & !< cloud nuclei concentration
         xc_frc,     & !< fractional cloud cover
         xm_co2,     & !< co2 mass in kg/m2
         xm_ch4,     & !< ch4 mass in kg/m2
         xm_n2o,     & !< n2o mass in kg/m2
         xm_cfc,   & !< cfc mass in kg/m2
         xm_o3,      & !< o3  mass in kg/m2
         xm_o2       !< o2  mass in kg/m2

    TYPE(c_ptr) :: src_data_cptr(33), dst_data_cptr(33)

    src_data_cptr( 1) = irad_aero
    src_data_cptr( 2) = klev
    src_data_cptr( 3) = ktype
    src_data_cptr( 4) = psctm
    src_data_cptr( 5) = ssi_factor
    src_data_cptr( 6) = loland
    src_data_cptr( 7) = loglac
    src_data_cptr( 8) = pcos_mu0
    src_data_cptr( 9) = daylght_frc
    src_data_cptr(10) = alb_vis_dir
    src_data_cptr(11) = alb_nir_dir
    src_data_cptr(12) = alb_vis_dif
    src_data_cptr(13) = alb_nir_dif
    src_data_cptr(14) = zf
    src_data_cptr(15) = zh
    src_data_cptr(16) = dz
    src_data_cptr(17) = pp_sfc
    src_data_cptr(18) = pp_fl
    src_data_cptr(19) = tk_sfc
    src_data_cptr(20) = tk_fl
    src_data_cptr(21) = tk_hl
    src_data_cptr(22) = xm_dry
    src_data_cptr(23) = xm_vap
    src_data_cptr(24) = xm_liq
    src_data_cptr(25) = xm_ice
    src_data_cptr(26) = cdnc
    src_data_cptr(27) = xc_frc
    src_data_cptr(28) = xm_co2
    src_data_cptr(29) = xm_ch4
    src_data_cptr(30) = xm_n2o
    src_data_cptr(31) = xm_cfc
    src_data_cptr(32) = xm_o3
    src_data_cptr(33) = xm_o2

    dst_data_cptr = src_data_cptr

#ifdef USE_PSRAD_COMMUNICATION
    CALL xt_redist_s_exchange(exchange_redist_atmo_2_psrad(idom), &
                              src_data_cptr, dst_data_cptr)
#else
    CALL finish("exchange_data_atmo_2_psrad", " Requires the YAXT library")
#endif


  END SUBROUTINE exchange_data_atmo_2_psrad


  SUBROUTINE exchange_data_psrad_2_atmo(idom, &
                                        lw_upw, &
                                        lw_upw_clr, &
                                        lw_dnw, &
                                        lw_dnw_clr, &
                                        sw_upw, &
                                        sw_upw_clr, &
                                        sw_dnw, &
                                        sw_dnw_clr, &
                                        vis_dn_dir_sfc, &
                                        par_dn_dir_sfc, &
                                        nir_dn_dir_sfc, &
                                        vis_dn_dff_sfc, &
                                        par_dn_dff_sfc, &
                                        nir_dn_dff_sfc, &
                                        vis_up_sfc, &
                                        par_up_sfc, &
                                        nir_up_sfc)

    INTEGER, INTENT(IN) :: idom

    !REAL(wp), INTENT(INOUT)   :: &
    TYPE(c_ptr), INTENT(IN)  :: &
         lw_upw,    & !< All-sky   upward   longwave  at all levels
         lw_upw_clr,& !< Clear-sky upward   longwave  at all levels
         lw_dnw,    & !< All-sky   downward longwave  at all levels
         lw_dnw_clr,& !< Clear-sky downward longwave  at all levels
         sw_upw,    & !< All-sky   upward   shortwave at all levels
         sw_upw_clr,& !< Clear-sky upward   shortwave at all levels
         sw_dnw,    & !< All-sky   downward shortwave at all levels
         sw_dnw_clr   !< Clear-sky downward shortwave at all levels

    !REAL(wp), INTENT(INOUT) :: &
    TYPE(c_ptr), INTENT(IN)  :: &
         vis_dn_dir_sfc, & !< Diffuse downward flux surface visible radiation
         par_dn_dir_sfc, & !< Diffuse downward flux surface PAR
         nir_dn_dir_sfc, & !< Diffuse downward flux surface near-infrared radiation
         vis_dn_dff_sfc, & !< Direct  downward flux surface visible radiation
         par_dn_dff_sfc, & !< Direct  downward flux surface PAR
         nir_dn_dff_sfc, & !< Direct  downward flux surface near-infrared radiation
         vis_up_sfc,     & !< Upward  flux surface visible radiation
         par_up_sfc,     & !< Upward  flux surface PAR
         nir_up_sfc        !< Upward  flux surface near-infrared radiation

    TYPE(c_ptr) :: src_data_cptr(17), dst_data_cptr(17)

    src_data_cptr( 1) = lw_upw
    src_data_cptr( 2) = lw_upw_clr
    src_data_cptr( 3) = lw_dnw
    src_data_cptr( 4) = lw_dnw_clr
    src_data_cptr( 5) = sw_upw
    src_data_cptr( 6) = sw_upw_clr
    src_data_cptr( 7) = sw_dnw
    src_data_cptr( 8) = sw_dnw_clr
    src_data_cptr( 9) = vis_dn_dir_sfc
    src_data_cptr(10) = par_dn_dir_sfc
    src_data_cptr(11) = nir_dn_dir_sfc
    src_data_cptr(12) = vis_dn_dff_sfc
    src_data_cptr(13) = par_dn_dff_sfc
    src_data_cptr(14) = nir_dn_dff_sfc
    src_data_cptr(15) = vis_up_sfc
    src_data_cptr(16) = par_up_sfc
    src_data_cptr(17) = nir_up_sfc
    dst_data_cptr = src_data_cptr

#ifdef USE_PSRAD_COMMUNICATION
    CALL xt_redist_s_exchange(exchange_redist_psrad_2_atmo(idom), &
                              src_data_cptr, dst_data_cptr)
#else
    CALL finish("exchange_data_atmo_2_psrad", " Requires the YAXT library")
#endif

  END SUBROUTINE exchange_data_psrad_2_atmo

END MODULE mo_psrad_communication