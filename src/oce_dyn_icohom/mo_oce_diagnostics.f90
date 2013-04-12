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
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range
  USE mo_mpi,                ONLY: my_process_is_stdio, p_field_sum, &
    &                              p_comm_work_test, p_comm_work, p_io, p_bcast
  USE mo_sync,               ONLY: global_sum_array
  USE mo_math_utilities,     ONLY: t_cartesian_coordinates!, gc2cc
  USE mo_math_constants,     ONLY: rad2deg
  USE mo_impl_constants,     ONLY: sea_boundary,sea, &
    &                              min_rlcell, min_rledge, min_rlcell, &
    &                              max_char_length, MIN_DOLIC
  USE mo_ocean_nml,          ONLY: n_zlev, no_tracer, &
      &                            ab_const, ab_beta, ab_gam, iswm_oce, idisc_scheme
  USE mo_dynamics_config,    ONLY: nold,nnew
  USE mo_parallel_config,    ONLY: nproma, p_test_run
  USE mo_run_config,         ONLY: dtime, nsteps
  USE mo_physical_constants, ONLY: grav, rho_ref
  USE mo_oce_state,          ONLY: t_hydro_ocean_state, t_hydro_ocean_diag,&! v_base, &
    &                              set_lateral_boundary_values
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D,t_patch_vert
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_exception,          ONLY: message, finish, message_text
  USE mo_loopindices,        ONLY: get_indices_c!, get_indices_e
  USE mo_oce_physics,        ONLY: t_ho_params
  USE mo_sea_ice_types,      ONLY: t_sfc_flx
  USE mo_datetime,           ONLY: t_datetime
  USE mo_linked_list,        ONLY: t_var_list
  USE mo_var_list,           ONLY: add_var,                  &
    &                              new_var_list,             &
    &                              delete_var_list,          &
    &                              default_var_list_settings,&
    &                              add_ref
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi_constants
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff
  USE mo_io_units,           ONLY: find_next_free_unit
IMPLICIT NONE

!PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'

INTEGER :: diag_unit = -1 ! file handle for the global timeseries output

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
    REAL(wp), ALLOCATABLE :: tracer_content(:)

END TYPE t_oce_monitor

TYPE t_oce_timeseries

    TYPE(t_oce_monitor), ALLOCATABLE :: oce_diagnostics(:)    ! time array of diagnostic values
    CHARACTER(len=40),DIMENSION(10)  :: names = (/ &
    & "volume                                  ", &
    & "kin_energy                              ", &
    & "pot_energy                              ", &
    & "total_energy                            ", &
    & "vorticity                               ", &
    & "enstrophy                               ", &
    & "potential_enstrophy                     ", &
    & "absolute_vertical_velocity              ", &
    & "total_temperature                       ", &
    & "total_salinity                          "/)

END TYPE t_oce_timeseries

CONTAINS
!-------------------------------------------------------------------------
!>
! !  calculate_oce_diagnostics
! 
! @par Revision History
! Developed  by  Peter Korn, MPI-M (2010).
! 
SUBROUTINE calculate_oce_diagnostics(p_patch_3D, p_os, p_sfc_flx, p_phys_param, timestep, oce_ts)
  TYPE(t_patch_3D ),TARGET, INTENT(INOUT)    :: p_patch_3D
  TYPE(t_hydro_ocean_state), TARGET          :: p_os
  TYPE(t_sfc_flx), INTENT(INOUT)             :: p_sfc_flx
  TYPE (t_ho_params)                         :: p_phys_param
  INTEGER                                    :: timestep
  TYPE(t_oce_timeseries),POINTER             :: oce_ts

  !Local variables
  INTEGER :: i_startidx_c, i_endidx_c!,i_startblk_c, i_endblk_c,
  INTEGER :: jk,jc,jb!,je
  INTEGER :: i_no_t, i
  REAL(wp) :: prism_vol = 0.0_wp, surface_height = 0.0_wp, prism_area = 0.0_wp
  REAL(wp) :: z_w
  REAL(wp) :: surface_area = 0.0_wp
  INTEGER  :: reference_timestep
  TYPE(t_patch), POINTER     :: p_patch

  TYPE(t_subset_range), POINTER :: owned_cells
  TYPE(t_oce_monitor),  POINTER :: monitor
  CHARACTER(len=1024)           :: line, nvars
  CHARACTER(len=1024)           :: fmt_string, real_fmt

  !-----------------------------------------------------------------------
  p_patch     => p_patch_3D%p_patch_2D(1)
  owned_cells => p_patch%cells%owned
  !-----------------------------------------------------------------------
  monitor     => oce_ts%oce_diagnostics(timestep)

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
        DO jk = 1,p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)

          ! area
          prism_area     = p_patch%cells%area(jc,jb)

          !local volume
          surface_height = merge(p_os%p_prog(nnew(1))%h(jc,jb),0.0_wp, 1 == jk)
          prism_vol      = prism_area * (p_patch_3D%p_patch_1D(1)%prism_thick_c(jc,1,jb) + surface_height)

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

          ! top level diagnostics
          IF (jk /= 1) cycle
          surface_area = surface_area + prism_area
          ! sum of top layer vertical velocities abolsute values
          monitor%absolute_vertical_velocity = &
            & monitor%absolute_vertical_velocity + abs(p_os%p_diag%w(jc,jk,jb))*prism_area
        END DO
      END DO
    END DO

  END SELECT

  ! compute global sums {
  monitor%volume              = global_sum_array(monitor%volume)
  monitor%kin_energy          = global_sum_array(monitor%kin_energy)/monitor%volume
  monitor%pot_energy          = global_sum_array(monitor%pot_energy)/monitor%volume
  monitor%total_energy        = global_sum_array(monitor%total_energy)/monitor%volume
  monitor%vorticity           = global_sum_array(monitor%vorticity)
  monitor%enstrophy           = global_sum_array(monitor%enstrophy)
  monitor%potential_enstrophy = global_sum_array(monitor%potential_enstrophy)
  surface_area                = global_sum_array(surface_area)
  monitor%absolute_vertical_velocity = global_sum_array(monitor%absolute_vertical_velocity)/surface_area
  DO i_no_t=1,no_tracer
    monitor%tracer_content(i_no_t) = global_sum_array(monitor%tracer_content(i_no_t))
  END DO
  ! } dbg_print

  ! write things to diagnostics output file
  real_fmt   = 'es26.18'
  ! * number of non-tracer diag. variables
  write(nvars,'(i3)') SIZE(oce_ts%names)-no_tracer
  write(fmt_string,'(a)') '(i5.5,1x,'//TRIM(ADJUSTL(nvars))//TRIM(real_fmt)//')'
  ! * non-tracer diags
  write(line,fmt_string) &
    & timestep, &
    & monitor%volume, &
    & monitor%kin_energy, &
    & monitor%pot_energy, &
    & monitor%total_energy, &
    & monitor%vorticity, &
    & monitor%enstrophy, &
    & monitor%potential_enstrophy, &
    & monitor%absolute_vertical_velocity

  ! * tracers
  DO i_no_t=1,no_tracer
    write(line,'(a,'//TRIM(real_fmt)//')') TRIM(line),monitor%tracer_content(i_no_t)
  END DO
  write(diag_unit,'(a)') TRIM(line)

END SUBROUTINE calculate_oce_diagnostics
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
SUBROUTINE construct_oce_diagnostics( p_patch_3D, p_os, oce_ts )
  TYPE(t_patch_3D ),TARGET, INTENT(INOUT) :: p_patch_3D
  TYPE(t_hydro_ocean_state), TARGET       :: p_os
  TYPE(t_oce_timeseries),POINTER          :: oce_ts

  !local variables
  INTEGER :: i,ist
  REAL(wp) :: delta_z, prism_vol
  TYPE(t_oce_monitor), POINTER :: ptr_monitor
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

  DO i=0,nsteps
    ALLOCATE(oce_ts%oce_diagnostics(i)%tracer_content(1:no_tracer))
    oce_ts%oce_diagnostics(i)%tracer_content(1:no_tracer) = 0.0_wp
  END DO

  ! open textfile for global timeseries
  diag_unit = find_next_free_unit(10,99)
  OPEN (unit=diag_unit,file='oce_diagnostics.txt',IOSTAT=ist)
  !HEADER
  headerLine = ''
  ! * add timestep columns
  write(headerLine,'(a)') 'step'
  ! * add columne for each monitored variable
  DO i=1,SIZE(oce_ts%names)
    WRITE(headerLine,'(a,a,a)')TRIM(headerLine),' ',TRIM(oce_ts%names(i))
  END DO
  print *,headerLine
  write(diag_unit,'(a)')TRIM(headerLine)

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
INTEGER :: i

CHARACTER(len=max_char_length), PARAMETER :: &
       & routine = ('mo_oce_diagnostics:destruct_oce_diagnostics')
!-----------------------------------------------------------------------
   DO i=0,nsteps 
     DEALLOCATE(oce_ts%oce_diagnostics(i)%tracer_content)
   END DO
   DEALLOCATE(oce_ts%oce_diagnostics)
   DEALLOCATE(oce_ts)
   ! close the global diagnostics text file
   CLOSE(unit=diag_unit)

CALL message (TRIM(routine), 'end')
END SUBROUTINE destruct_oce_diagnostics
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
! TODO: calculate in parallel
!! 
SUBROUTINE calc_moc (p_patch,p_patch_3D, w, datetime)

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
          !  - not yet parallelized
          DO lbr = -jbrei, jbrei
            lbrei = NINT(90.0_wp + z_lat_deg + REAL(lbr, wp) * z_lat_dim)
            lbrei = MAX(lbrei,1)
            lbrei = MIN(lbrei,180)
      
            global_moc(lbrei,jk) = global_moc(lbrei,jk) - &
              &                    p_patch%cells%area(jc,jb) * rho_ref * w(jc,jk,jb) / &
              &                    REAL(2*jbrei + 1, wp)

            IF (p_patch_3D%basin_c(jc,jb) <2) THEN

              atlant_moc(lbrei,jk) = atlant_moc(lbrei,jk) - &
                &                    p_patch%cells%area(jc,jb) * rho_ref * w(jc,jk,jb) / &
                &                    REAL(2*jbrei + 1, wp)
            ELSE
              pacind_moc(lbrei,jk) = pacind_moc(lbrei,jk) - &
                &                    p_patch%cells%area(jc,jb) * rho_ref * w(jc,jk,jb) / &
                &                    REAL(2*jbrei + 1, wp)
            END IF

        !   if (jk==3) write(81,*) 'jb,jc,lbr,lbrei,zlatc,zlat,glb:', &
        !     & jb,jc,lbr,lbrei,lbr+int(z_lat_deg),z_lat_deg,global_moc(lbrei,jk)

          END DO
   
        END IF
      END DO
    END DO

    ! test parallelization:
    ! new function field_sum_all using mpi_allreduce and working precisions wp
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
      write(77) i1,i2,i3,i4
      write(77) (global_moc(lbr,jk),lbr=1,180)
      i2 = INT(778,i8)
      write(77) i1,i2,i3,i4
      write(77) (atlant_moc(lbr,jk),lbr=1,180)
      i2 = INT(779,i8)
      write(77) i1,i2,i3,i4
      write(77) (pacind_moc(lbr,jk),lbr=1,180)

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
