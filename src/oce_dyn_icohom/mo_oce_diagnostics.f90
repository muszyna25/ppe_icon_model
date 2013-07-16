!>
!! Contains basic diagnostics for ICON ocean model.
!!
!!
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2011/02)
!!  Extended   by Stephan Lorenz,   MPI-M (2012)
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
MODULE mo_oce_diagnostics
  USE mo_kind,               ONLY: wp, dp, i8
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range, t_subset_indexed
  USE mo_grid_tools,         ONLY: get_oriented_edges_from_global_vertices
  USE mo_mpi,                ONLY: my_process_is_stdio, p_field_sum, get_my_mpi_work_id, &
    &                              p_comm_work_test, p_comm_work, p_io, p_bcast
  USE mo_sync,               ONLY: global_sum_array
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates
  USE mo_math_constants,     ONLY: rad2deg
  USE mo_impl_constants,     ONLY: sea_boundary,sea, &
    &                              min_rlcell, min_rledge, min_rlcell, &
    &                              max_char_length, MIN_DOLIC
  USE mo_ocean_nml,          ONLY: n_zlev, no_tracer, &
    &                              gibraltar, denmark_strait,drake_passage, indonesian_throughflow,&
    &                              scotland_iceland, &
    &                              ab_const, ab_beta, ab_gam, iswm_oce, idisc_scheme
  USE mo_dynamics_config,    ONLY: nold,nnew
  USE mo_parallel_config,    ONLY: nproma, p_test_run
  USE mo_run_config,         ONLY: dtime, nsteps
  USE mo_physical_constants, ONLY: grav, rho_ref
  USE mo_oce_state,          ONLY: t_hydro_ocean_state, t_hydro_ocean_diag,&
    &                              set_lateral_boundary_values
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D,t_patch_vert, t_grid_edges
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_exception,          ONLY: message, finish, message_text
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_oce_physics,        ONLY: t_ho_params
  USE mo_sea_ice_types,      ONLY: t_sfc_flx, t_sea_ice
  USE mo_datetime,           ONLY: t_datetime, datetime_to_string, date_len
  USE mo_linked_list,        ONLY: t_var_list
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_io_units,           ONLY: find_next_free_unit
  USE mo_util_file,          ONLY: util_symlink, util_rename, util_islink, util_unlink
  USE mo_statistics,         ONLY: subset_sum

IMPLICIT NONE

!PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'

INTEGER :: diag_unit = -1 ! file handle for the global timeseries output
INTEGER :: moc_unit = -1 ! file handle for the global timeseries output
CHARACTER(len=max_char_length)  :: diag_fname, moc_fname

!
! PUBLIC INTERFACE
!
PUBLIC :: calculate_oce_diagnostics
PUBLIC :: construct_oce_diagnostics
PUBLIC :: destruct_oce_diagnostics
PUBLIC :: t_oce_monitor
PUBLIC :: t_oce_timeseries
PUBLIC :: calc_moc
PUBLIC :: calc_psi

TYPE t_oce_monitor
    REAL(wp) :: volume
    REAL(wp) :: kin_energy
    REAL(wp) :: pot_energy
    REAL(wp) :: total_energy
    REAL(wp) :: vorticity
    REAL(wp) :: enstrophy
    REAL(wp) :: potential_enstrophy
    REAL(wp) :: absolute_vertical_velocity
    REAL(wp) :: forc_swflx    ! surface short wave heat flux                              [W/m2]
    REAL(wp) :: forc_lwflx    ! surface long wave heat flux                               [W/m2]
    REAL(wp) :: forc_ssflx    ! surface sensible heat flux                                [W/m2]
    REAL(wp) :: forc_slflx    ! surface latent heat flux                                  [W/m2]
    REAL(wp) :: forc_precip   ! total precipitation flux                                  [m/s]
    REAL(wp) :: forc_evap     ! evaporation flux                                          [m/s]
    REAL(wp) :: forc_runoff   ! river runoff flux                                         [m/s]
    REAL(wp) :: forc_fwbc     ! sum of forcing surface freshwater flux from BC            [m/s]
    REAL(wp) :: forc_fwrelax  ! diagnosed surface freshwater flux due to relaxation       [m/s]
    REAL(wp) :: forc_fwfx     ! diagnosed sum of forcing surface freshwater flux          [m/s]
    REAL(wp) :: forc_hfrelax  ! diagnosed surface heat flux due to relaxation             [m/s]
    REAL(wp) :: forc_hflx     ! diagnosed sum of forcing surface heat flux                [W/m2]
    REAL(wp) :: ice_volume_nh !                                                           [km3]
    REAL(wp) :: ice_volume_sh !                                                           [km3]
    REAL(wp) :: ice_extent_nh !                                                           [km2]
    REAL(wp) :: ice_extent_sh !                                                           [km2]
    REAL(wp) :: gibraltar     ! though flow                                               [Sv]
    REAL(wp) :: denmark_strait! though flow                                               [Sv]
    REAL(wp) :: drake_passage ! though flow                                               [Sv]
    REAL(wp) :: indonesian_throughflow !                                                  [Sv]
    REAL(wp) :: scotland_iceland !                                                        [Sv]
    REAL(wp), ALLOCATABLE :: tracer_content(:)

END TYPE t_oce_monitor

TYPE t_oce_timeseries

    TYPE(t_oce_monitor), ALLOCATABLE :: oce_diagnostics(:)    ! time array of diagnostic values
    CHARACTER(len=40), DIMENSION(31)  :: names = (/ &
      & "volume                                  ", &
      & "kin_energy                              ", &
      & "pot_energy                              ", &
      & "total_energy                            ", &
      & "vorticity                               ", &
      & "enstrophy                               ", &
      & "potential_enstrophy                     ", &
      & "absolute_vertical_velocity              ", &
      & "forc_swflx                              ", &
      & "forc_lwflx                              ", &
      & "forc_ssflx                              ", &
      & "forc_slflx                              ", &
      & "forc_precip                             ", &
      & "forc_evap                               ", &
      & "forc_runoff                             ", &
      & "forc_fwbc                               ", &
      & "forc_fwrelax                            ", &
      & "forc_fwfx                               ", &
      & "forc_hfrelax                            ", &
      & "forc_hflx                               ", &
      & "ice_volume_nh                           ", &
      & "ice_volume_sh                           ", &
      & "ice_extent_nh                           ", &
      & "ice_extent_sh                           ", &
      & "gibraltar                               ", &
      & "denmark_strait                          ", &
      & "drake_passage                           ", &
      & "indonesian_throughflow                  ", &
      & "scotland_iceland                        ", &
      & "total_temperature                       ", &
      & "total_salinity                          "/)

END TYPE t_oce_timeseries

TYPE t_oce_section
  TYPE(t_subset_indexed) :: subset
  REAL(wp), POINTER  :: orientation(:)
END TYPE t_oce_section

INTEGER, PARAMETER  :: oce_section_count = 5
PRIVATE             :: oce_section_count
TYPE(t_oce_section) :: oce_sections(oce_section_count)
PRIVATE             :: oce_sections

CONTAINS
!-------------------------------------------------------------------------
!
!
!  !The constructor of the types related to ocean diagnostics
!>
!!
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2011).
!!
SUBROUTINE construct_oce_diagnostics( p_patch_3D, p_os, oce_ts, datestring )
  TYPE(t_patch_3D ),TARGET, INTENT(INOUT) :: p_patch_3D
  TYPE(t_hydro_ocean_state), TARGET       :: p_os
  TYPE(t_oce_timeseries),POINTER          :: oce_ts
  CHARACTER(len=32)                       :: datestring

  !local variable
  INTEGER :: i,ist
  CHARACTER(len=max_char_length), PARAMETER :: &
    & routine = ('mo_oce_diagnostics:construct_oce_diagnostics')
  !-----------------------------------------------------------------------
  CHARACTER(len=max_char_length)      :: headerLine
  TYPE(t_patch), POINTER              :: p_patch
  CHARACTER(len=max_char_length)      :: listname
  INTEGER :: nblks_c,nblks_e,nblks_v
  !-----------------------------------------------------------------------
  p_patch   => p_patch_3D%p_patch_2D(1)

  CALL message (TRIM(routine), 'start')
  ALLOCATE(oce_ts)

  ALLOCATE(oce_ts%oce_diagnostics(0:nsteps))

  oce_ts%oce_diagnostics(0:nsteps)%volume                     = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%kin_energy                 = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%pot_energy                 = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%total_energy               = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%vorticity                  = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%enstrophy                  = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%potential_enstrophy        = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%absolute_vertical_velocity = 0.0_wp

  oce_ts%oce_diagnostics(0:nsteps)%forc_swflx                 = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%forc_lwflx                 = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%forc_ssflx                 = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%forc_slflx                 = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%forc_precip                = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%forc_evap                  = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%forc_runoff                = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%forc_fwbc                  = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%forc_fwrelax               = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%forc_fwfx                  = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%forc_hfrelax               = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%forc_hflx                  = 0.0_wp

  oce_ts%oce_diagnostics(0:nsteps)%ice_volume_nh              = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%ice_volume_sh              = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%ice_extent_nh              = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%ice_extent_sh              = 0.0_wp

  ! through flows
  oce_ts%oce_diagnostics(0:nsteps)%gibraltar                  = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%denmark_strait             = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%drake_passage              = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%indonesian_throughflow     = 0.0_wp
  oce_ts%oce_diagnostics(0:nsteps)%scotland_iceland           = 0.0_wp

  DO i=0,nsteps
    ALLOCATE(oce_ts%oce_diagnostics(i)%tracer_content(1:no_tracer))
    oce_ts%oce_diagnostics(i)%tracer_content(1:no_tracer) = 0.0_wp
  END DO

  ! open textfile for global timeseries
  diag_fname = 'oce_diagnostics-'//TRIM(datestring)//'.txt'
  diag_unit = find_next_free_unit(10,99)
  OPEN (unit=diag_unit,file=diag_fname,IOSTAT=ist)
  ! header of the text file
  headerLine = ''
  ! * add timestep columns
  write(headerLine,'(a)') 'step date time'
  ! * add columne for each monitored variable
  DO i=1,SIZE(oce_ts%names)
    WRITE(headerLine,'(a,a,a)')TRIM(headerLine),' ',TRIM(oce_ts%names(i))
  END DO
  write(diag_unit,'(a)')TRIM(headerLine)

  ! open file for MOC - extraordinary at this time
  moc_fname='MOC.'//TRIM(datestring)
  moc_unit = find_next_free_unit(10,99)
  OPEN (moc_unit,file=moc_fname,form='unformatted')
  WRITE(message_text,'(2a)') ' MOC-file opened successfully, filename=',TRIM(moc_fname)
  CALL message (TRIM(routine), message_text)

  ! compute subsets for given sections path allong edges
  CALL get_oriented_edges_from_global_vertices(    &
     &  edge_subset = oce_sections(1)%subset,      &
     &  orientation = oce_sections(1)%orientation, &
     &  patch_3D = p_patch_3D,                     &
     & global_vertex_array = gibraltar,            &
     & subset_name = 'gibraltar')
  CALL get_oriented_edges_from_global_vertices(    &
     &  edge_subset = oce_sections(2)%subset,      &
     &  orientation = oce_sections(2)%orientation, &
     &  patch_3D = p_patch_3D,                     &
     & global_vertex_array =denmark_strait,        &
     & subset_name = 'denmark_strait')
  CALL get_oriented_edges_from_global_vertices(    &
     &  edge_subset = oce_sections(3)%subset,      &
     &  orientation = oce_sections(3)%orientation, &
     &  patch_3D = p_patch_3D,                     &
     & global_vertex_array =drake_passage,         &
     & subset_name = 'drake_passage')
  CALL get_oriented_edges_from_global_vertices(    &
     &  edge_subset = oce_sections(4)%subset,      &
     &  orientation = oce_sections(4)%orientation, &
     &  patch_3D = p_patch_3D,                     &
     & global_vertex_array =indonesian_throughflow,&
     & subset_name = 'indonesian_throughflow')
  CALL get_oriented_edges_from_global_vertices(    &
     &  edge_subset = oce_sections(5)%subset,      &
     &  orientation = oce_sections(5)%orientation, &
     &  patch_3D = p_patch_3D,                     &
     & global_vertex_array =scotland_iceland,      &
     & subset_name = 'scotland_iceland')

  CALL message (TRIM(routine), 'end')
END SUBROUTINE construct_oce_diagnostics
!-------------------------------------------------------------------------
!
!
!  !The destructor of the types related to ocean diagnostics
!>
!!
!!
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2011).
!!
SUBROUTINE destruct_oce_diagnostics(oce_ts)
!
TYPE(t_oce_timeseries),POINTER         :: oce_ts
!
!local variables
INTEGER :: i,iret
CHARACTER(len=max_char_length)  :: linkname
CHARACTER(len=max_char_length)  :: message_text

CHARACTER(len=max_char_length), PARAMETER :: &
       & routine = ('mo_oce_diagnostics:destruct_oce_diagnostics')
!-----------------------------------------------------------------------
   DO i=0,nsteps
     DEALLOCATE(oce_ts%oce_diagnostics(i)%tracer_content)
   END DO
   DEALLOCATE(oce_ts%oce_diagnostics)
   DEALLOCATE(oce_ts)
   ! close the global diagnostics text file and the SRV MOC file
   CLOSE(unit=diag_unit)
   CLOSE(unit=moc_unit)
   ! create a link to the last diagnostics file
   linkname = 'oce_diagnostics.txt'
   IF (util_islink(TRIM(linkname))) THEN
     iret = util_unlink(TRIM(linkname))
   ENDIF
   iret = util_symlink(TRIM(diag_fname),TRIM(linkname))
   WRITE(message_text,'(t1,a,t50,a)') TRIM(diag_fname), TRIM(linkname)
   CALL message('',message_text)

CALL message (TRIM(routine), 'end')
END SUBROUTINE destruct_oce_diagnostics
!-------------------------------------------------------------------------
!>
! !  calculate_oce_diagnostics
!
! @par Revision History
! Developed  by  Peter Korn, MPI-M (2010).
!
SUBROUTINE calculate_oce_diagnostics(p_patch_3D, p_os, p_sfc_flx, p_ice, &
    & p_phys_param, timestep, datetime, oce_ts)
  TYPE(t_patch_3D ),TARGET, INTENT(IN)    :: p_patch_3D
  TYPE(t_hydro_ocean_state), TARGET       :: p_os
  TYPE(t_sfc_flx),    INTENT(IN)          :: p_sfc_flx
  TYPE (t_sea_ice),   INTENT(IN)          :: p_ice
  TYPE (t_ho_params)                      :: p_phys_param
  INTEGER                                 :: timestep
  TYPE(t_datetime), INTENT(IN)            :: datetime
  TYPE(t_oce_timeseries),POINTER          :: oce_ts

  !Local variables
  INTEGER :: i_startidx_c, i_endidx_c!,i_startblk_c, i_endblk_c,
  INTEGER :: jk,jc,jb!,je
  INTEGER :: i_no_t, i
  REAL(wp):: prism_vol, surface_height, prism_area, surface_area, z_w
  INTEGER :: reference_timestep
  TYPE(t_patch), POINTER     :: p_patch
  REAL(wp) :: sflux

  TYPE(t_subset_range), POINTER :: owned_cells
  TYPE(t_oce_monitor),  POINTER :: monitor
  CHARACTER(len=1024)           :: line, nvars
  CHARACTER(len=1024)           :: fmt_string, real_fmt
  CHARACTER(len=date_len)       :: datestring
  REAL(wp), PARAMETER           :: equator = 0.00001_wp

  !-----------------------------------------------------------------------
  p_patch        => p_patch_3D%p_patch_2D(1)
  owned_cells    => p_patch%cells%owned
  !-----------------------------------------------------------------------
  monitor        => oce_ts%oce_diagnostics(timestep)
  surface_area   = 0.0_wp
  surface_height = 0.0_wp
  prism_vol      = 0.0_wp
  prism_area     = 0.0_wp
  z_w            = 0.0_wp
  CALL datetime_to_string(datestring, datetime, plain=.TRUE.)

  !cell loop to calculate cell based monitored fields volume, kinetic energy and tracer content
  SELECT CASE (iswm_oce)
  CASE (1) ! shallow water mode
    !Potential energy in SW-casep_patch%patch_oce%del_zlev_m(1)
    DO jb = owned_cells%start_block,owned_cells%end_block
      CALL get_index_range(owned_cells, jb, i_startidx_c, i_endidx_c)

      DO jc =  i_startidx_c, i_endidx_c
        IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN
          prism_vol = p_patch%cells%area(jc,jb)*p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,1,jb)
            !prism_vol = p_op_coeff%fixed_vol_norm(jc,1,jb)*p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,1,jb)!p_os%p_prog(nold(1))%h(jc,jb)
          monitor%volume = monitor%volume + prism_vol


          monitor%pot_energy = monitor%pot_energy &
            &+ 0.5_wp*grav* prism_vol*p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,1,jb)

          monitor%kin_energy = monitor%kin_energy &
            &+p_os%p_diag%kin(jc,1,jb)*prism_vol!p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,1,jb)

          monitor%total_energy=monitor%kin_energy+monitor%pot_energy
          DO i_no_t=1, no_tracer
            monitor%tracer_content(i_no_t) = monitor%tracer_content(i_no_t)&
              & + prism_vol*p_os%p_prog(nold(1))%tracer(jc,1,jb,i_no_t)
          END DO
        ENDIF
      END DO
    END DO

  CASE DEFAULT !3D model
    DO jb = owned_cells%start_block, owned_cells%end_block
    CALL get_index_range(owned_cells, jb, i_startidx_c, i_endidx_c)
      !We are dealing with the surface layer first
      DO jc =  i_startidx_c, i_endidx_c

        ! area
        prism_area   = p_patch%cells%area(jc,jb)
        surface_area = surface_area + prism_area
        ! sum of top layer vertical velocities abolsute values
        monitor%absolute_vertical_velocity = &
          & monitor%absolute_vertical_velocity + abs(p_os%p_diag%w(jc,1,jb))*prism_area

        monitor%forc_swflx   = monitor%forc_swflx   + p_sfc_flx%forc_swflx(jc,jb)*prism_area
        monitor%forc_lwflx   = monitor%forc_lwflx   + p_sfc_flx%forc_lwflx(jc,jb)*prism_area
        monitor%forc_ssflx   = monitor%forc_ssflx   + p_sfc_flx%forc_ssflx(jc,jb)*prism_area
        monitor%forc_slflx   = monitor%forc_slflx   + p_sfc_flx%forc_slflx(jc,jb)*prism_area
        monitor%forc_precip  = monitor%forc_precip  + p_sfc_flx%forc_precip(jc,jb)*prism_area
        monitor%forc_evap    = monitor%forc_evap    + p_sfc_flx%forc_evap(jc,jb)*prism_area
        monitor%forc_runoff  = monitor%forc_runoff  + p_sfc_flx%forc_runoff(jc,jb)*prism_area
        monitor%forc_fwbc    = monitor%forc_fwbc    + p_sfc_flx%forc_fwbc(jc,jb)*prism_area
        monitor%forc_fwrelax = monitor%forc_fwrelax + p_sfc_flx%forc_fwrelax(jc,jb)*prism_area
        monitor%forc_fwfx    = monitor%forc_fwfx    + p_sfc_flx%forc_fwfx(jc,jb)*prism_area
        monitor%forc_hfrelax = monitor%forc_hfrelax + p_sfc_flx%forc_hfrelax(jc,jb)*prism_area
        monitor%forc_hflx    = monitor%forc_hflx    + p_sfc_flx%forc_hflx(jc,jb)*prism_area

        ! northern hemisphere
        IF (p_patch%cells%center(jc,jb)%lat > equator) THEN
          monitor%ice_volume_nh  = monitor%ice_volume_nh + prism_area*SUM(p_ice%vol(jc,:,jb)*p_ice%conc(jc,:,jb))
          monitor%ice_extent_nh  = monitor%ice_extent_nh + p_ice%concSum(jc,jb)*prism_area
        ELSE
          ! southern hemisphere
          monitor%ice_volume_sh  = monitor%ice_volume_sh + prism_area*SUM(p_ice%vol(jc,:,jb)*p_ice%conc(jc,:,jb))
          monitor%ice_extent_sh  = monitor%ice_extent_sh + p_ice%concSum(jc,jb)*prism_area
        END IF


        DO jk = 1,p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)

          !local volume
          surface_height = merge(p_os%p_prog(nnew(1))%h(jc,jb),0.0_wp, 1 == jk)
          prism_vol      = prism_area * (p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,jk,jb) + surface_height)

          !Fluid volume
          monitor%volume = monitor%volume + prism_vol

          !kinetic energy
          monitor%kin_energy = monitor%kin_energy + p_os%p_diag%kin(jc,jk,jb)*prism_vol

          !Potential energy
          IF(jk==1)THEN
            z_w = (p_os%p_diag%w(jc,jk,jb)*p_os%p_prog(nold(1))%h(jc,jb)&
              & +p_os%p_diag%w(jc,jk+1,jb)*0.5_wp*p_patch_3D%p_patch_1D(1)%del_zlev_i(jk))&
              &/(0.5_wp*p_patch_3D%p_patch_1D(1)%del_zlev_i(jk)+p_os%p_prog(nold(1))%h(jc,jb))
          ELSEIF(jk>1.AND.jk<n_zlev)THEN
            z_w = (p_os%p_diag%w(jc,jk,jb)*p_patch_3D%p_patch_1D(1)%del_zlev_i(jk)&
              & +p_os%p_diag%w(jc,jk+1,jb)*p_patch_3D%p_patch_1D(1)%del_zlev_i(jk+1))&
              &/(p_patch_3D%p_patch_1D(1)%del_zlev_i(jk)+p_patch_3D%p_patch_1D(1)%del_zlev_i(jk+1))
          ENDIF
          monitor%pot_energy = monitor%pot_energy + grav*z_w* p_os%p_diag%rho(jc,jk,jb)* prism_vol

          !Tracer content
          DO i_no_t=1, no_tracer
            monitor%tracer_content(i_no_t) = &
              & monitor%tracer_content(i_no_t) + prism_vol*p_os%p_prog(nold(1))%tracer(jc,jk,jb,i_no_t)
          END DO
        END DO
      END DO
    END DO

  END SELECT

  ! compute global sums {
  monitor%volume                     = global_sum_array(monitor%volume)
  surface_area                       = global_sum_array(surface_area)
  monitor%kin_energy                 = global_sum_array(monitor%kin_energy)/monitor%volume
  monitor%pot_energy                 = global_sum_array(monitor%pot_energy)/monitor%volume
  monitor%total_energy               = global_sum_array(monitor%total_energy)/monitor%volume
  monitor%vorticity                  = global_sum_array(monitor%vorticity)
  monitor%enstrophy                  = global_sum_array(monitor%enstrophy)
  monitor%potential_enstrophy        = global_sum_array(monitor%potential_enstrophy)
  monitor%absolute_vertical_velocity = global_sum_array(monitor%absolute_vertical_velocity)/surface_area
  monitor%forc_swflx                 = global_sum_array(monitor%forc_swflx)/surface_area
  monitor%forc_lwflx                 = global_sum_array(monitor%forc_lwflx)/surface_area
  monitor%forc_ssflx                 = global_sum_array(monitor%forc_ssflx)/surface_area
  monitor%forc_slflx                 = global_sum_array(monitor%forc_slflx)/surface_area
  monitor%forc_precip                = global_sum_array(monitor%forc_precip)/surface_area
  monitor%forc_evap                  = global_sum_array(monitor%forc_evap)/surface_area
  monitor%forc_runoff                = global_sum_array(monitor%forc_runoff)/surface_area
  monitor%forc_fwbc                  = global_sum_array(monitor%forc_fwbc)/surface_area
  monitor%forc_fwrelax               = global_sum_array(monitor%forc_fwrelax)/surface_area
  monitor%forc_fwfx                  = global_sum_array(monitor%forc_fwfx)/surface_area
  monitor%forc_hfrelax               = global_sum_array(monitor%forc_hfrelax)/surface_area
  monitor%forc_hflx                  = global_sum_array(monitor%forc_hflx)/surface_area
  monitor%ice_volume_nh              = global_sum_array(monitor%ice_volume_nh)/1.0e9_wp
  monitor%ice_volume_sh              = global_sum_array(monitor%ice_volume_sh)/1.0e9_wp
  monitor%ice_extent_nh              = global_sum_array(monitor%ice_extent_nh)/1.0e6_wp
  monitor%ice_extent_sh              = global_sum_array(monitor%ice_extent_sh)/1.0e6_wp
  DO i_no_t=1,no_tracer
    monitor%tracer_content(i_no_t) = global_sum_array(monitor%tracer_content(i_no_t))
  END DO
  ! fluxes through given paths
  IF (my_process_is_stdio()) &
    & write(0,*) "---------------  fluxes --------------------------------"
  DO i=1,oce_section_count
    sflux = section_flux(oce_sections(i), p_os%p_prog(nnew(1))%vn)
!
! #slo# disabled since subset%block is not allocated (#3759, HPC_sun_debug)
! #ifdef NOMPI
!     IF (my_process_is_stdio()) &
!       & write(0,*) oce_sections(i)%subset%name, ":", sflux, 'at edges:',oce_sections(i)%subset%block
! #else
    IF (my_process_is_stdio()) &
      & write(0,*) oce_sections(i)%subset%name, ":", sflux
! #endif

    SELECT CASE (i)
    CASE (1)
      monitor%gibraltar              = sflux*rho_ref
    CASE (2)
      monitor%denmark_strait         = sflux*rho_ref
    CASE (3)
      monitor%drake_passage          = sflux*rho_ref
    CASE (4)
      monitor%indonesian_throughflow = sflux*rho_ref
    CASE (5)
      monitor%scotland_iceland       = sflux*rho_ref
    END SELECT
  ENDDO
  IF (my_process_is_stdio()) &
    & write(0,*) "---------------  end fluxes ----------------------------"
  ! } dbg_print


  IF (my_process_is_stdio()) THEN
    ! write things to diagnostics output file
    real_fmt   = 'es26.18'
    ! * number of non-tracer diag. variables
    write(nvars,'(i3)') SIZE(oce_ts%names)-no_tracer
    write(fmt_string,'(a)') '(i5.5,1x,a,1x,'//TRIM(ADJUSTL(nvars))//TRIM(real_fmt)//')'
    ! create date and time string
    ! * non-tracer diags
    write(line,fmt_string) &
      & timestep, &
      & TRIM(datestring), &
      & monitor%volume, &
      & monitor%kin_energy, &
      & monitor%pot_energy, &
      & monitor%total_energy, &
      & monitor%vorticity, &
      & monitor%enstrophy, &
      & monitor%potential_enstrophy, &
      & monitor%absolute_vertical_velocity, &
      & monitor%forc_swflx, &
      & monitor%forc_lwflx, &
      & monitor%forc_ssflx, &
      & monitor%forc_slflx, &
      & monitor%forc_precip, &
      & monitor%forc_evap, &
      & monitor%forc_runoff, &
      & monitor%forc_fwbc, &
      & monitor%forc_fwrelax, &
      & monitor%forc_fwfx, &
      & monitor%forc_hfrelax, &
      & monitor%forc_hflx, &
      & monitor%ice_volume_nh, &
      & monitor%ice_volume_sh, &
      & monitor%ice_extent_nh, &
      & monitor%ice_extent_sh, &
      & monitor%gibraltar, &
      & monitor%denmark_strait, &
      & monitor%drake_passage,  &
      & monitor%indonesian_throughflow, &
      & monitor%scotland_iceland
    ! * tracers
    DO i_no_t=1,no_tracer
    write(line,'(a,'//TRIM(real_fmt)//')') TRIM(line),monitor%tracer_content(i_no_t)
    END DO

    write(diag_unit,'(a)') TRIM(line)
  END IF

END SUBROUTINE calculate_oce_diagnostics
!-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  REAL(wp) FUNCTION section_flux(in_oce_section, velocity_values)
    TYPE(t_oce_section) :: in_oce_section
    REAL(wp), POINTER :: velocity_values(:,:,:)

    INTEGER :: i, k, edge_idx, edge_block
    REAL(wp) :: oriented_length
    REAL(wp), ALLOCATABLE :: flux_weights(:,:)
    TYPE(t_grid_edges), POINTER ::  edges
    TYPE(t_patch_vert),POINTER :: patch_vertical

    CHARACTER(LEN=*), PARAMETER :: method_name='mo_oce_diagnostics:section_flux'

    edges          => in_oce_section%subset%patch%edges
    patch_vertical => in_oce_section%subset%patch_3D%p_patch_1D(1)

    ! calculate weights
    ! flux_weights can also be preallocated
    ALLOCATE(flux_weights(n_zlev, MAX(in_oce_section%subset%size, 1)))
    flux_weights(:,:) = 0.0_wp
    DO i=1, in_oce_section%subset%size

      edge_idx   = in_oce_section%subset%idx(i)
      edge_block = in_oce_section%subset%block(i)
      oriented_length = edges%primal_edge_length(edge_idx, edge_block) * &
        & in_oce_section%orientation(i) ! this can also be pre-calculated and stored in in_oce_section%orientation

      !write(0,*) "oriented_length:",  oriented_length

      DO k=1, n_zlev
        flux_weights(k, i) = patch_vertical%prism_thick_e(edge_idx, k, edge_block) * oriented_length ! maybe also use slm
       !write(0,*) i, k, in_oce_section%subset%name, " flux_weights:",  flux_weights(k, i), &
       !  & patch_vertical%prism_thick_e(edge_idx, k, edge_block)
       !write(0,*) i, k, in_oce_section%subset%name, " velocity_value:", velocity_values(edge_idx, k, edge_block)
      ENDDO

    ENDDO


    section_flux = subset_sum(                           &
      & values                 = velocity_values,        &
      & indexed_subset         = in_oce_section%subset,  &
      & subset_indexed_weights = flux_weights)

    DEALLOCATE(flux_weights)

    !write(0,*) get_my_mpi_work_id(), ": section_flux on subset ", in_oce_section%subset%name, ":", &
    !  & section_flux, in_oce_section%subset%size

  END FUNCTION section_flux
  !-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!
!
!!  Calculation of meridional overturning circulation (MOC)
!
!   Calculation of meridional overturning circulation for different basins
!   (Atlantic, Pacific, Indian, global)
!>
!!
!! @par Revision History
!! Developed  by  Stephan Lorenz, MPI-M (2012).
!!  based on code from MPIOM
!
! TODO: implement variable output dimension (1 deg resolution) and smoothing extent
! TODO: calculate the 1 deg resolution meridional distance
!!
SUBROUTINE calc_moc (p_patch, p_patch_3D, w, datetime)

  TYPE(t_patch), TARGET, INTENT(IN)  :: p_patch
  TYPE(t_patch_3D ),TARGET, INTENT(INOUT)  :: p_patch_3D
  REAL(wp), INTENT(in)               :: w(:,:,:)   ! vertical velocity at cell centers
                                                   ! dims: (nproma,nlev+1,nblks_c)
  TYPE(t_datetime), INTENT(IN)       :: datetime
  !
  ! local variables
  ! INTEGER :: i
  INTEGER, PARAMETER ::  jbrei=3   !  latitudinal smoothing area is 2*jbrei-1 rows of 1 deg
  INTEGER :: jb, jc, jk, i_startidx, i_endidx !, il_e, ib_e
  INTEGER :: lbrei, lbr, idate, itime
  INTEGER :: mpi_comm
  INTEGER(i8) :: i1,i2,i3,i4

  REAL(wp) :: z_lat, z_lat_deg, z_lat_dim
  REAL(wp) :: global_moc(180,n_zlev), atlant_moc(180,n_zlev), pacind_moc(180,n_zlev)
  REAL(dp) :: local_moc(180), res_moc(180)

  TYPE(t_subset_range), POINTER :: dom_cells

  CHARACTER(len=max_char_length), PARAMETER :: routine = ('mo_oce_diagnostics:calc_moc')

  !-----------------------------------------------------------------------

  IF(p_test_run) THEN
    mpi_comm = p_comm_work_test
  ELSE
    mpi_comm = p_comm_work
  ENDIF

  global_moc(:,:) = 0.0_wp
  pacind_moc(:,:) = 0.0_wp
  atlant_moc(:,:) = 0.0_wp

  ! set barrier:
  ! CALL MPI_BARRIER(0)

  ! with all cells no sync is necessary
  !owned_cells => p_patch%cells%owned
  dom_cells   => p_patch%cells%in_domain

  !write(81,*) 'MOC: datetime:',datetime

  DO jk = 1, n_zlev   !  not yet on intermediate levels
    DO jb = dom_cells%start_block, dom_cells%end_block
      CALL get_index_range(dom_cells, jb, i_startidx, i_endidx)
      DO jc = i_startidx, i_endidx

        !  could be replaced by vertical loop to bottom
        IF ( p_patch_3D%lsm_c(jc,jk,jb) <= sea_boundary ) THEN

          ! lbrei: corresponding latitude row of 1 deg extension
          !       1 south pole
          !     180 north pole
          z_lat = p_patch%cells%center(jc,jb)%lat
          z_lat_deg = z_lat*rad2deg
          lbrei = NINT(90.0_wp + z_lat_deg)
          lbrei = MAX(lbrei,1)
          lbrei = MIN(lbrei,180)

          ! get neighbor edge for scaling
      !   il_e = p_patch%cells%edge_idx(jc,jb,1)
      !   ib_e = p_patch%cells%edge_blk(jc,jb,1)

          ! z_lat_dim: scale to 1 deg resolution
          ! z_lat_dim: latitudinal extent of triangle divided by latitudinal smoothing extent
      !   z_lat_dim = p_patch%edges%primal_edge_length(il_e,ib_e) / &
      !     & (REAL(2*jbrei, wp) * 111111._wp*1.3_wp)
          z_lat_dim = 1.0_wp

          ! distribute MOC over (2*jbrei)+1 latitude rows
          !  - no weighting with latitudes done
          !  - lbrei: index of 180 X 1 deg meridional resolution
          DO lbr = -jbrei, jbrei
            lbrei = NINT(90.0_wp + z_lat_deg + REAL(lbr, wp) * z_lat_dim)
            lbrei = MAX(lbrei,1)
            lbrei = MIN(lbrei,180)

            global_moc(lbrei,jk) = global_moc(lbrei,jk) - &
              !  multiply with wet (or loop to bottom)
              &                    p_patch%cells%area(jc,jb) * rho_ref * w(jc,jk,jb) * &
              &                    p_patch_3D%wet_c(jc,jk,jb) / &
              &                    REAL(2*jbrei + 1, wp)

            IF (p_patch_3D%basin_c(jc,jb) == 1) THEN         !  1: Atlantic; 0: Land

              atlant_moc(lbrei,jk) = atlant_moc(lbrei,jk) - &
                &                    p_patch%cells%area(jc,jb) * rho_ref * w(jc,jk,jb) * &
                &                    p_patch_3D%wet_c(jc,jk,jb) / &
                &                    REAL(2*jbrei + 1, wp)
            ELSE IF (p_patch_3D%basin_c(jc,jb) >= 2) THEN   !  2: Indian; 4: Pacific
              pacind_moc(lbrei,jk) = pacind_moc(lbrei,jk) - &
                &                    p_patch%cells%area(jc,jb) * rho_ref * w(jc,jk,jb) * &
                &                    p_patch_3D%wet_c(jc,jk,jb) / &
                &                    REAL(2*jbrei + 1, wp)
            END IF

          END DO

        END IF
      END DO
    END DO

    ! test parallelization:
    ! function field_sum_all using mpi_allreduce and working precisions wp does not exist
  ! res_moc(:) = p_field_sum_all_wp(global_moc(:,jk))
  ! res_moc(:) = p_field_sum_all_wp(atlant_moc(:,jk))
  ! res_moc(:) = p_field_sum_all_wp(pacind_moc(:,jk))

    ! function field_sum using mpi_reduce, then broadcast
    local_moc(:)     = REAL(global_moc(:,jk),dp)
    res_moc(:)       = p_field_sum(local_moc, mpi_comm)
    CALL p_bcast(res_moc(:), p_io, mpi_comm)
    global_moc(:,jk) = REAL(res_moc(:),wp)

    local_moc(:)     = REAL(atlant_moc(:,jk),dp)
    res_moc(:)       = p_field_sum(local_moc, mpi_comm)
    CALL p_bcast(res_moc(:), p_io, mpi_comm)
    atlant_moc(:,jk) = REAL(res_moc(:),wp)

    local_moc(:)     = REAL(pacind_moc(:,jk),dp)
    res_moc(:)       = p_field_sum(local_moc, mpi_comm)
    CALL p_bcast(res_moc(:), p_io, mpi_comm)
    pacind_moc(:,jk) = REAL(res_moc(:),wp)

  END DO  ! n_zlev-loop

  IF (my_process_is_stdio()) THEN
    DO lbr=179,1,-1   ! fixed to 1 deg meridional resolution

        global_moc(lbr,:)=global_moc(lbr+1,:)+global_moc(lbr,:)
        atlant_moc(lbr,:)=atlant_moc(lbr+1,:)+atlant_moc(lbr,:)
        pacind_moc(lbr,:)=pacind_moc(lbr+1,:)+pacind_moc(lbr,:)

    END DO

    ! write out MOC in extra format, file opened in mo_hydro_ocean_run  - integer*8
    !  - correct date in extra format - i.e YYYYMMDD - no time info
    !idate=datetime%month*1000000+datetime%day*10000+datetime%hour*100+datetime%minute
    idate = datetime%year*10000+datetime%month*100+datetime%day
    itime = datetime%hour*100+datetime%minute
    WRITE(message_text,*) 'Write MOC at year =',datetime%year,', date =',idate,' time =', itime
    CALL message (TRIM(routine), message_text)

    DO jk = 1,n_zlev
      i1=INT(idate,i8)
      i2 = INT(777,i8)
      i3 = INT(p_patch_3D%p_patch_1D(1)%zlev_i(jk),i8)
      i4 = INT(180,i8)
      write(moc_unit) i1,i2,i3,i4
      write(moc_unit) (global_moc(lbr,jk),lbr=1,180)
      i2 = INT(778,i8)
      write(moc_unit) i1,i2,i3,i4
      write(moc_unit) (atlant_moc(lbr,jk),lbr=1,180)
      i2 = INT(779,i8)
      write(moc_unit) i1,i2,i3,i4
      write(moc_unit) (pacind_moc(lbr,jk),lbr=1,180)

    END DO
  END IF

END SUBROUTINE calc_moc
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!
!
!!  Calculation of horizontal stream function
!
!>
!!
!! @par Revision History
!! Developed  by  Stephan Lorenz, MPI-M (2012).
!!  based on code from MPIOM
!
! TODO: implement variable output dimension (1 deg resolution) and smoothing extent
!!
SUBROUTINE calc_psi (p_patch,p_patch_3D, u, h, u_vint, datetime)

  TYPE(t_patch), TARGET, INTENT(IN)  :: p_patch
  TYPE(t_patch_3D ),TARGET, INTENT(INOUT)  :: p_patch_3D
  REAL(wp), INTENT(IN)               :: u(:,:,:)     ! zonal velocity at cell centers
  REAL(wp), INTENT(IN)               :: h(:,:)       ! elevation on cell centers
                                                     ! dims: (nproma,nlev,nblks_c)
  REAL(wp), INTENT(OUT)              :: u_vint(:,:)  ! barotropic zonal velocity on icon grid
  TYPE(t_datetime), INTENT(IN)       :: datetime
  !
  ! local variables
  ! INTEGER :: i

  ! switch for writing stream function (not yet in namelist); 1: icon-grid; 2: regular grid output
  INTEGER, PARAMETER ::  idiag_psi = 1

  INTEGER, PARAMETER ::  nlat = 180                    ! meridional dimension of regular grid
  INTEGER, PARAMETER ::  nlon = 360                    ! zonal dimension of regular grid

  ! smoothing area is 2*jsmth-1 lat/lon areas of 1 deg
  INTEGER, PARAMETER ::  jsmth = 3
  INTEGER            :: jb, jc, jk, i_startidx, i_endidx
  INTEGER            :: jlat, jlon, jlt, jln, jltx, jlnx, jsmth2
  INTEGER(i8)        :: idate, iextra(4)


  REAL(wp) :: z_lat_deg, z_lon_deg, z_lat_dist, delta_z, rsmth
  REAL(wp) :: z_uint_reg(nlon,nlat)                     ! vertical integral on regular grid
  REAL(wp) :: psi_reg(nlon,nlat)                        ! horizontal stream function

  TYPE(t_subset_range), POINTER :: all_cells, dom_cells

  !CHARACTER(len=max_char_length), PARAMETER :: routine = ('mo_oce_diagnostics:calc_psi')

  !-----------------------------------------------------------------------

  psi_reg(:,:)    = 0.0_wp
  z_uint_reg(:,:) = 0.0_wp

  jsmth2          = 2*jsmth + 1
  rsmth           = REAL(jsmth2*jsmth2, wp)


  ! with all cells no sync is necessary
  all_cells => p_patch%cells%all
  dom_cells => p_patch%cells%in_domain

  ! (1) barotropic system:
  !     vertical integration of zonal velocity times vertical layer thickness [m/s*m]
  u_vint(:,:)     = 0.0_wp
  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
    DO jk = 1, n_zlev

      DO jc = i_startidx, i_endidx
        delta_z = p_patch_3D%p_patch_1D(1)%del_zlev_m(jk)
        IF (jk == 1) delta_z = p_patch_3D%p_patch_1D(1)%del_zlev_m(jk) + h(jc,jb)
        u_vint(jc,jb) = u_vint(jc,jb) - u(jc,jk,jb)*delta_z*p_patch_3D%wet_c(jc,jk,jb)
      END DO
    END DO
  END DO

  IF (idiag_psi == 1) RETURN

  ! (2) distribute integrated zonal velocity (u*dz) on 1x1 deg grid
  !     this code is not mature yet

  ! in domain: count all cells only once
  DO jb = dom_cells%start_block, dom_cells%end_block
    CALL get_index_range(dom_cells, jb, i_startidx, i_endidx)
    DO jc = i_startidx, i_endidx
      z_lat_deg = p_patch%cells%center(jc,jb)%lat * rad2deg
      z_lon_deg = p_patch%cells%center(jc,jb)%lon * rad2deg

   !  ! 0 <= lon <= 360 deg
   !  z_lon_deg = z_lon_deg + 180.0_wp

      ! jlat/jlon: corresponding latitude/longitude coordinates of 1 deg extension
      ! jlat: 1 = south of 89.0S; 89 = 1S-Eq.; 90 = Eq-1N;  180 = north of 89N
      ! jlon: 1 = 180W-179W; 180 = 1-0 deg west; 360 = 179E-180E

      jlat = NINT(91.0_wp + z_lat_deg)
      jlon = NINT(z_lon_deg + 180.5_wp)

      ! distribute stream function over rsmth=(2*jsmth+1)**2 lat/lon regular grid points
      !  - no weighting with latitudes done
      !  - no correction with regular lsm done
      DO jltx = jlat-jsmth, jlat+jsmth

        jlt = jltx
        IF (jlt <    1) jlt =      1-jlt  ! apply equatorwards
        IF (jlt > nlat) jlt = 2*nlat-jlt  ! apply equatorwards
        DO jlnx = jlon-jsmth, jlon+jsmth

          jln = jlnx
          IF (jln <    1) jln = jln+nlon  ! circular boundary
          IF (jln > nlon) jln = jln-nlon  ! circular boundary

          z_uint_reg(jln,jlt) = z_uint_reg(jln,jlt) + u_vint(jc,jb) / rsmth

 ! 99 format('J lat=',f8.2,' lon=',f8.2,' jlat=',i4,' jlon=',i4,' lsm=',i3, &
 !      &    ' jlt=',i4,  ' jln=',i4,' uint=',1p10e12.3)
 ! 98 format(' lat=',f8.2,' lon=',f8.2,' jlat=',i4,' jlon=',i4,' lsm=',i3, &
 !      &    ' uint=',1p10e12.3)
 !    if ((jlat==101 .and. jlon==270) &
 !      & write(82,99) z_lat_deg,z_lon_deg,jlat,jlon,v_base%lsm_c(jc,1,jb), &
 !      &              jlt,jln,z_uint_reg(jln,jlt)

        END DO
      END DO
 !    write(82,98) z_lat_deg,z_lon_deg,jlat,jlon,v_base%lsm_c(jc,1,jb),z_uint_reg(jlon,jlat)

    END DO
  END DO

  ! (3) calculate meridional integral on regular grid starting from south pole:

  DO jlt = nlat-1, 1, -1
    z_uint_reg(:,jlt) = z_uint_reg(:,jlt) + z_uint_reg(:,jlt+1)
  END DO

  ! (4) calculate stream function: scale with length of 1 deg*rho [m/s*m*m*kg/m3=kg/s]

  ! meridional distance of 1 deg
  ! ATTENTION - fixed 1 deg resolution should be related to icon-resolution
  z_lat_dist = 111111.0_wp  ! * 1.3_wp ??

  psi_reg(:,:) = z_uint_reg(:,:) * z_lat_dist * rho_ref

  ! stream function on icon grid without calculation of meridional integral
  !  - tbd after interpolation to regular grid externally
!  psi    (:,:) = u_vint    (:,:)              * rho_ref


  ! write out in extra format - integer*8
  idate = INT(datetime%month*1000000+datetime%day*10000+datetime%hour*100+datetime%minute,i8)
  write(0,*) 'write global PSI at iyear, idate:',datetime%year, idate

  iextra(1) = INT(idate,i8)
  iextra(2) = INT(780,i8)
  iextra(3) = INT(0,i8)
  iextra(4) = INT(nlon*nlat,i8)

  write(80) (iextra(jb),jb=1,4)
  write(80) ((psi_reg(jln,jlt),jln=1,nlon),jlt=1,nlat)

  do jlat=1,nlat
      write(82,*) 'jlat=',jlat
      write(82,'(1p10e12.3)') (psi_reg(jlon,jlat),jlon=1,nlon)
  enddo

END SUBROUTINE calc_psi
!-------------------------------------------------------------------------

END MODULE mo_oce_diagnostics
