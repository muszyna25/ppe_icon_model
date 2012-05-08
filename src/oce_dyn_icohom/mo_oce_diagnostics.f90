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
!-------------------------------------------------------------------------  
!
!    ProTeX FORTRAN source: Style 2  
!    modified for ICON project, DWD/MPI-M 2006
!  
!-------------------------------------------------------------------------  
!  
!   
! 
USE mo_kind,                      ONLY: wp, i8
USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
USE mo_math_utilities,            ONLY: t_cartesian_coordinates!, gc2cc
USE mo_math_constants,            ONLY: rad2deg
USE mo_impl_constants,            ONLY: sea_boundary,sea, &
  &                                     min_rlcell, min_rledge, min_rlcell, &
  &                                     max_char_length, MIN_DOLIC
USE mo_ocean_nml,                 ONLY: n_zlev, no_tracer, &
  &                                     ab_const, ab_beta, ab_gam, iswm_oce, idisc_scheme
USE mo_dynamics_config,          ONLY: nold,nnew
USE mo_parallel_config,  ONLY: nproma
USE mo_run_config,                ONLY: dtime, nsteps
USE mo_physical_constants,        ONLY: grav, rho_ref
USE mo_oce_state,                 ONLY: t_hydro_ocean_state, t_hydro_ocean_diag, v_base, &
  &                                     set_lateral_boundary_values
USE mo_model_domain,              ONLY: t_patch
USE mo_ext_data_types,            ONLY: t_external_data
USE mo_exception,                 ONLY: message, finish!, message_text
USE mo_loopindices,               ONLY: get_indices_c, get_indices_e !, get_indices_v
USE mo_oce_math_operators,        ONLY: div_oce, grad_fd_norm_oce, grad_fd_norm_oce_2d,&
  &                                     height_related_quantities
USE mo_oce_physics,               ONLY: t_ho_params
USE mo_sea_ice,                   ONLY: t_sfc_flx
USE mo_intp_data_strc,            ONLY: t_int_state
USE mo_scalar_product,            ONLY: calc_scalar_product_veloc, calc_scalar_product_veloc_3D
USE mo_intp_data_strc,            ONLY: t_int_state
USE mo_intp_rbf,                  ONLY: rbf_vec_interpol_cell, rbf_vec_interpol_edge
USE mo_oce_ab_timestepping,       ONLY: calc_vert_velocity
USE mo_datetime,                  ONLY: t_datetime

IMPLICIT NONE

!PRIVATE

! !VERSION CONTROL:
CHARACTER(len=*), PARAMETER :: version = '$Id$'

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
    REAL(wp), ALLOCATABLE :: tracer_content(:)

END TYPE t_oce_monitor

TYPE t_oce_timeseries

    TYPE(t_oce_monitor), ALLOCATABLE :: oce_diagnostics(:)    ! time array of diagnostic values

END TYPE t_oce_timeseries


CONTAINS
!-------------------------------------------------------------------------  
!
!  
!>
!! !  Solves the free surface equation.
!! 
!! @par Revision History
!! Developed  by  Peter Korn, MPI-M (2010).
!! 
SUBROUTINE calculate_oce_diagnostics(p_patch, p_os, p_sfc_flx, p_phys_param, timestep, oce_ts)
!
TYPE(t_patch), TARGET, INTENT(in)             :: p_patch
TYPE(t_hydro_ocean_state), TARGET             :: p_os 
TYPE(t_sfc_flx), INTENT(INOUT)                :: p_sfc_flx
TYPE (t_ho_params)                            :: p_phys_param
INTEGER                                       :: timestep
TYPE(t_oce_timeseries),POINTER                :: oce_ts
!
!Local variables
!REAL(wp) :: z_h_c(nproma,p_patch%nblks_c)
!REAL(wp) :: z_h_e(nproma,p_patch%nblks_e)

!CHARACTER(len=max_char_length) :: string
INTEGER :: rl_start_c, rl_end_c, i_startblk_c, i_endblk_c,i_startidx_c, i_endidx_c
!INTEGER :: rl_start_e, rl_end_e, i_startblk_e, i_endblk_e!, i_startidx_e, i_endidx_e
INTEGER :: jk,jc,jb!,je
INTEGER :: i_no_t, i !, z_dolic

!REAL(wp) :: z_volume, z_volume_initial
!REAL(wp) :: z_kin_energy,z_kin_energy_initial
!REAL(wp) :: z_pot_energy,z_pot_energy_initial
!REAL(wp) :: z_total_energy,z_total_energy_initial
!REAL(wp) :: z_vorticity,z_vorticity_initial
!REAL(wp) :: z_enstrophy,z_enstrophy_initial
!REAL(wp) :: z_potential_enstrophy,z_potential_enstrophy_initial
!REAL(wp) :: z_tracer_content(1:no_tracer), z_tracer_content_initial(1:no_tracer)
REAL(wp) :: delta_z, prism_vol
REAL(wp) :: z_w
INTEGER  :: referenz_timestep
TYPE(t_oce_monitor), POINTER :: ptr_monitor
!LOGICAL :: l_first_timestep
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diagnostics:calculate_oce_diagnotics')
!-----------------------------------------------------------------------
rl_start_c = 1
rl_end_c   = min_rlcell
i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)
! rl_start_e = 1
! rl_end_e   = min_rledge
! i_startblk_e = p_patch%edges%start_blk(rl_start_e,1)
! i_endblk_e   = p_patch%edges%end_blk(rl_end_e,1)



!direct pointer to monitored quantitiy at actual timestep
ptr_monitor =>oce_ts%oce_diagnostics(timestep)

!cell loop to calculate cell based monitored fields volume, kinetic energy and tracer content
IF(iswm_oce/=1)THEN

  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
    &                             rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c

      !z_dolic = v_base%dolic_c(jc,jb)
      !IF ( z_dolic>=MIN_DOLIC)THEN 
        DO jk=1,n_zlev!z_dolic

          IF ( v_base%lsm_oce_c(jc,jk,jb) == sea ) THEN

!              delta_z=p_os%p_diag%cons_thick_c(jc,jk,jb)
             delta_z = v_base%del_zlev_m(jk)! 
             IF (jk == 1) THEN
               delta_z = v_base%del_zlev_m(jk)&
                      & + p_os%p_prog(nnew(1))%h(jc,jb)
             ENDIF

            prism_vol = p_patch%cells%area(jc,jb)*delta_z

            !Fluid volume 
            ptr_monitor%volume = ptr_monitor%volume + prism_vol

            !kinetic energy
            ptr_monitor%kin_energy = ptr_monitor%kin_energy+ p_os%p_diag%kin(jc,jk,jb)*prism_vol

            !Potential energy
            IF(jk==1)THEN
              z_w = (p_os%p_diag%w(jc,jk,jb)*p_os%p_prog(nold(1))%h(jc,jb)&
                 & +p_os%p_diag%w(jc,jk+1,jb)*0.5_wp*v_base%del_zlev_i(jk))&
                 &/(0.5_wp*v_base%del_zlev_i(jk)+p_os%p_prog(nold(1))%h(jc,jb))
            ELSEIF(jk>1.AND.jk<n_zlev)THEN
              z_w = (p_os%p_diag%w(jc,jk,jb)*v_base%del_zlev_i(jk)&
                 & +p_os%p_diag%w(jc,jk+1,jb)*v_base%del_zlev_i(jk+1))&
                 &/(v_base%del_zlev_i(jk)+v_base%del_zlev_i(jk+1))
            ENDIF 

            ptr_monitor%pot_energy = ptr_monitor%pot_energy&
            &+ grav*z_w* p_os%p_diag%rho(jc,jk,jb)* prism_vol

            !Tracer content
            DO i_no_t=1, no_tracer
              ptr_monitor%tracer_content(i_no_t) = ptr_monitor%tracer_content(i_no_t)&
              & + prism_vol*p_os%p_prog(nold(1))%tracer(jc,jk,jb,i_no_t)*v_base%wet_c(jc,jk,jb)
            END DO
          ENDIF
        END DO
      !ENDIF
    END DO
  END DO

!   !divide by volume
!   IF(ptr_monitor%volume/=0.0_wp)THEN
!     ptr_monitor%kin_energy = ptr_monitor%kin_energy/ptr_monitor%volume
!     ptr_monitor%pot_energy = ptr_monitor%pot_energy/ptr_monitor%volume
!   ELSE
!     ptr_monitor%kin_energy = 0.0_wp
!     ptr_monitor%pot_energy = 0.0_wp
!   ENDIF
!   ptr_monitor%total_energy = ptr_monitor%pot_energy + ptr_monitor%kin_energy
! 
!   IF(ptr_monitor%volume/=0.0_wp)THEN
!     DO i_no_t=1, no_tracer
!       ptr_monitor%tracer_content(i_no_t) = ptr_monitor%tracer_content(i_no_t)/ptr_monitor%volume
!     END DO
!   ELSE
!     DO i_no_t=1, no_tracer
!       ptr_monitor%tracer_content(i_no_t) = 0.0_wp
!     END DO
!   ENDIF

ELSEIF(iswm_oce==1)THEN
  !Potential energy in SW-casep_patch%patch_oce%del_zlev_m(1)
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
       &                             rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c
      IF ( v_base%lsm_oce_e(jc,1,jb) <= sea_boundary ) THEN
        prism_vol = p_patch%cells%area(jc,jb)*p_os%p_prog(nold(1))%h(jc,jb)
        ptr_monitor%volume = ptr_monitor%volume + prism_vol

        ptr_monitor%pot_energy = ptr_monitor%pot_energy &
        &+ 0.5_wp*grav*p_os%p_prog(nold(1))%h(jc,jb)*p_os%p_prog(nold(1))%h(jc,jb)

        ptr_monitor%kin_energy = ptr_monitor%kin_energy &
        &+p_os%p_diag%kin(jc,1,jb)*p_os%p_prog(nold(1))%h(jc,jb)

        DO i_no_t=1, no_tracer
          ptr_monitor%tracer_content(i_no_t) = ptr_monitor%tracer_content(i_no_t)&
          & + prism_vol*p_os%p_prog(nold(1))%tracer(jc,1,jb,i_no_t)
        END DO
      ENDIF
    END DO
  END DO
!   IF(ptr_monitor%volume/=0.0_wp)THEN
!     DO i_no_t=1, no_tracer
!       ptr_monitor%tracer_content(i_no_t) = ptr_monitor%tracer_content(i_no_t)/ptr_monitor%volume
!     END DO
!   ELSE
!     DO i_no_t=1, no_tracer
!       ptr_monitor%tracer_content(i_no_t) = 0.0_wp
!     END DO
!   ENDIF
!   IF(ptr_monitor%volume/=0.0_wp)THEN
!     ptr_monitor%pot_energy = ptr_monitor%pot_energy/ptr_monitor%volume
!     ptr_monitor%kin_energy = ptr_monitor%kin_energy/ptr_monitor%volume
!   ELSE
!     ptr_monitor%kin_energy = 0.0_wp
!     ptr_monitor%pot_energy = 0.0_wp
!   ENDIF 
!   ptr_monitor%total_energy = ptr_monitor%pot_energy + ptr_monitor%kin_energy
ENDIF

DO i=timestep,timestep
  write(0,*)'ACTUAL VALUES OF VOLUME NORMALIZED BY INITIAL VALUE:          ', i,&
  &oce_ts%oce_diagnostics(i)%volume/oce_ts%oce_diagnostics(0)%volume

  IF(oce_ts%oce_diagnostics(1)%kin_energy/=0.0_wp)THEN
    write(0,*)'ACTUAL VALUES OF KINETIC ENERGY NORMALIZED BY INITIAL VALUE:  ',  i,&
    &oce_ts%oce_diagnostics(i)%kin_energy/ oce_ts%oce_diagnostics(1)%kin_energy
  ENDIF

  IF(oce_ts%oce_diagnostics(0)%pot_energy/=0.0_wp)THEN
    write(0,*)'ACTUAL VALUES OF POTENTIAL ENERGY NORMALIZED BY INITIAL VALUE:',i,&
    &oce_ts%oce_diagnostics(i)%pot_energy/ oce_ts%oce_diagnostics(0)%pot_energy
  ENDIF

  IF(oce_ts%oce_diagnostics(0)%total_energy/=0.0_wp)THEN
    write(0,*)'ACTUAL VALUES OF TOTAL ENERGY NORMALIZED BY INITIAL VALUE:    ',i,&
    &oce_ts%oce_diagnostics(i)%total_energy/ oce_ts%oce_diagnostics(0)%total_energy
  ENDIF
  DO i_no_t=1, no_tracer
    IF(oce_ts%oce_diagnostics(0)%tracer_content(i_no_t)/=0.0_wp)THEN
      write(0,*)'ACTUAL VALUES OF TOTAL TRACER CONTENT:                        ',i_no_t,i,&
      &oce_ts%oce_diagnostics(i)%tracer_content(i_no_t)&
      &/ oce_ts%oce_diagnostics(0)%tracer_content(i_no_t)

  ENDIF
END DO
END DO


!IF(timestep==nsteps)THEN
  referenz_timestep=1
  !DO i=referenz_timestep,timestep
  DO i=timestep,timestep
    write(1234,*)'ACTUAL VALUES OF VOLUME NORMALIZED BY INITIAL VALUE:          ', i,&
    &oce_ts%oce_diagnostics(i)%volume&
    &/oce_ts%oce_diagnostics(referenz_timestep)%volume

!     IF(oce_ts%oce_diagnostics(1)%kin_energy/=0.0_wp)THEN
!       write(1234,*)'ACTUAL VALUES OF KINETIC ENERGY NORMALIZED BY INITIAL VALUE:  ',  i,&
!       &oce_ts%oce_diagnostics(i)%kin_energy&
!       &/ oce_ts%oce_diagnostics(referenz_timestep)%kin_energy
!     ENDIF
! 
!     IF(oce_ts%oce_diagnostics(0)%pot_energy/=0.0_wp)THEN
!       write(1234,*)'ACTUAL VALUES OF POTENTIAL ENERGY NORMALIZED BY INITIAL VALUE:',i,&
!       &oce_ts%oce_diagnostics(i)%pot_energy&
!       &/ oce_ts%oce_diagnostics(referenz_timestep)%pot_energy
!     ENDIF
! 
!     IF(oce_ts%oce_diagnostics(0)%total_energy/=0.0_wp)THEN
!       write(1234,*)'ACTUAL VALUES OF TOTAL ENERGY NORMALIZED BY INITIAL VALUE:    ',i,&
!       &oce_ts%oce_diagnostics(i)%total_energy&
!       &/ oce_ts%oce_diagnostics(referenz_timestep)%total_energy
!     ENDIF
     DO i_no_t=1, no_tracer
       IF(oce_ts%oce_diagnostics(0)%tracer_content(i_no_t)/=0.0_wp)THEN
         write(1234,*)'ACTUAL VALUES OF TOTAL TRACER CONTENT:                        ',i_no_t,i,&
         &oce_ts%oce_diagnostics(i)%tracer_content(i_no_t)&
         &/ oce_ts%oce_diagnostics(referenz_timestep)%tracer_content(i_no_t)
       ENDIF
    END DO
  END DO
!ENDIF
!CALL message (TRIM(routine), 'end')
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
SUBROUTINE construct_oce_diagnostics(p_patch, p_os, p_ext_data, oce_ts)
!
TYPE(t_patch), TARGET, INTENT(in)             :: p_patch
TYPE(t_hydro_ocean_state), TARGET             :: p_os
TYPE(t_external_data), TARGET, INTENT(IN)     :: p_ext_data 
TYPE(t_oce_timeseries),POINTER                :: oce_ts
!
!local variables
INTEGER :: i
!CHARACTER(len=max_char_length) :: string
INTEGER :: rl_start_c, rl_end_c, i_startblk_c, i_endblk_c,i_startidx_c, i_endidx_c
!INTEGER :: rl_start_e, rl_end_e, i_startblk_e, i_endblk_e, i_startidx_e, i_endidx_e
INTEGER :: jk,jc,jb!,je
INTEGER :: i_no_t
REAL(wp) :: delta_z, prism_vol
REAL(wp) :: z_w
!TYPE(t_cartesian_coordinates) :: z_u_cc(nproma,n_zlev,p_patch%nblks_c)
TYPE(t_oce_monitor), POINTER :: ptr_monitor
CHARACTER(len=max_char_length), PARAMETER :: &
       & routine = ('mo_oce_diagnostics:construct_oce_diagnostics')
!-----------------------------------------------------------------------
  CALL message (TRIM(routine), 'start')
  ALLOCATE(oce_ts)

  ALLOCATE(oce_ts%oce_diagnostics(0:nsteps))

   oce_ts%oce_diagnostics(0:nsteps)%volume              = 0.0_wp
   oce_ts%oce_diagnostics(0:nsteps)%kin_energy          = 0.0_wp
   oce_ts%oce_diagnostics(0:nsteps)%pot_energy          = 0.0_wp
   oce_ts%oce_diagnostics(0:nsteps)%total_energy        = 0.0_wp
   oce_ts%oce_diagnostics(0:nsteps)%vorticity           = 0.0_wp
   oce_ts%oce_diagnostics(0:nsteps)%potential_enstrophy = 0.0_wp

   DO i=0,nsteps 
     ALLOCATE(oce_ts%oce_diagnostics(i)%tracer_content(1:no_tracer))
     oce_ts%oce_diagnostics(i)%tracer_content(1:no_tracer) = 0.0_wp
   END DO


rl_start_c = 1
rl_end_c   = min_rlcell
i_startblk_c = p_patch%cells%start_blk(rl_start_c,1)
i_endblk_c   = p_patch%cells%end_blk(rl_end_c,1)


!calculate initial values
  ptr_monitor =>oce_ts%oce_diagnostics(0)
  IF(iswm_oce/=1)THEN

    DO jb = i_startblk_c, i_endblk_c
      CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
       &                             rl_start_c, rl_end_c)

      DO jk=1,n_zlev

        DO jc = i_startidx_c, i_endidx_c

          IF ( v_base%lsm_oce_c(jc,jk,jb) == sea ) THEN

            delta_z = v_base%del_zlev_m(jk)
            IF (jk == 1) THEN
             delta_z = v_base%del_zlev_m(jk)&
                   & + p_os%p_prog(nold(1))%h(jc,jb)
            ENDIF

            prism_vol = p_patch%cells%area(jc,jb)*delta_z

            !Fluid volume 
            ptr_monitor%volume = ptr_monitor%volume + prism_vol

            !kinetic energy
            ptr_monitor%kin_energy = ptr_monitor%kin_energy&
                              &+ p_os%p_diag%kin(jc,jk,jb)*prism_vol

            !Potential energy
            IF(jk==1)THEN
              z_w = (p_os%p_diag%w(jc,jk,jb)*p_os%p_prog(nold(1))%h(jc,jb)&
              & +p_os%p_diag%w(jc,jk+1,jb)*0.5_wp*v_base%del_zlev_i(jk))&
              &/(0.5_wp*v_base%del_zlev_i(jk)+p_os%p_prog(nold(1))%h(jc,jb))
            ELSEIF(jk<n_zlev)THEN
              z_w = (p_os%p_diag%w(jc,jk,jb)*v_base%del_zlev_i(jk)&
              & +p_os%p_diag%w(jc,jk+1,jb)*v_base%del_zlev_i(jk+1))&
              &/(v_base%del_zlev_i(jk)+v_base%del_zlev_i(jk+1))
          ENDIF 

          ptr_monitor%pot_energy = ptr_monitor%pot_energy&
          &+ grav*z_w* p_os%p_diag%rho(jc,jk,jb)* prism_vol

          !Tracer content
          DO i_no_t=1, no_tracer
            ptr_monitor%tracer_content(i_no_t) = ptr_monitor%tracer_content(i_no_t)&
            & + prism_vol*p_os%p_prog(nold(1))%tracer(jc,jk,jb,i_no_t)*v_base%wet_c(jc,jk,jb)
          END DO
        ENDIF
      END DO
    END DO
  END DO

  !divide by volume
  !ptr_monitor%kin_energy = ptr_monitor%kin_energy/ptr_monitor%volume
  !ptr_monitor%pot_energy = ptr_monitor%pot_energy/ptr_monitor%volume
  !ptr_monitor%total_energy = ptr_monitor%pot_energy + ptr_monitor%kin_energy

  !DO i_no_t=1, no_tracer
  !  ptr_monitor%tracer_content(i_no_t) = ptr_monitor%tracer_content(i_no_t)!/ptr_monitor%volume
  !END DO


ELSEIF(iswm_oce==1)THEN
  !Potential energy in SW-case
  DO jb = i_startblk_c, i_endblk_c
    CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
       &                             rl_start_c, rl_end_c)
    DO jc = i_startidx_c, i_endidx_c

      prism_vol = p_patch%cells%area(jc,jb)*p_os%p_prog(nold(1))%h(jc,jb)
      ptr_monitor%volume = ptr_monitor%volume + prism_vol

      ptr_monitor%pot_energy = ptr_monitor%pot_energy &
      &+ 0.5_wp*grav*p_os%p_prog(nold(1))%h(jc,jb)*p_os%p_prog(nold(1))%h(jc,jb)

      ptr_monitor%kin_energy = ptr_monitor%kin_energy &
      &+p_os%p_diag%kin(jc,1,jb)*p_os%p_prog(nold(1))%h(jc,jb)

          !Tracer content
          DO i_no_t=1, no_tracer
            ptr_monitor%tracer_content(i_no_t) = ptr_monitor%tracer_content(i_no_t)&
            & + prism_vol*p_os%p_prog(nold(1))%tracer(jc,1,jb,i_no_t)*v_base%wet_c(jc,1,jb)
          END DO
    END DO
  END DO


!ptr_monitor%pot_energy = ptr_monitor%pot_energy/ptr_monitor%volume
!ptr_monitor%kin_energy = ptr_monitor%kin_energy/ptr_monitor%volume
!ptr_monitor%total_energy = ptr_monitor%pot_energy + ptr_monitor%kin_energy
ENDIF

write(*,*)'INITIAL VALUES OF VOLUME          :',oce_ts%oce_diagnostics(0)%volume
write(*,*)'INIIAL VALUES OF KINETIC ENERGY   :',oce_ts%oce_diagnostics(0)%kin_energy
write(*,*)'INITIAL VALUES OF POTENTIAL ENERGY:',oce_ts%oce_diagnostics(0)%pot_energy
write(*,*)'INITIAL VALUES OF TOTAL ENERGY    :',oce_ts%oce_diagnostics(0)%total_energy


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
!! 
SUBROUTINE calc_moc (p_patch, w, datetime)

  TYPE(t_patch), TARGET, INTENT(IN)  :: p_patch
  REAL(wp), INTENT(in)               :: w(:,:,:)   ! vertical velocity at cell centers
                                                   ! dims: (nproma,nlev+1,nblks_c)
  TYPE(t_datetime), INTENT(IN)       :: datetime
  !
  ! local variables
  ! INTEGER :: i
  INTEGER, PARAMETER ::  jbrei=3   !  latitudinal smoothing area is 2*jbrei-1 rows of 1 deg
  INTEGER :: jb, jc, jk, i_startidx, i_endidx !, il_e, ib_e
  INTEGER :: lbrei, lbr, idate
  INTEGER(i8) :: i1,i2,i3,i4

  REAL(wp) :: z_lat, z_lat_deg, z_lat_dim
  REAL(wp) :: global_moc(180,n_zlev), atlant_moc(180,n_zlev), pacind_moc(180,n_zlev)

  TYPE(t_subset_range), POINTER :: dom_cells
  
  !CHARACTER(len=max_char_length), PARAMETER :: routine = ('mo_oce_diagnostics:calc_moc')

  !-----------------------------------------------------------------------

  global_moc(:,:) = 0.0_wp
  pacind_moc(:,:) = 0.0_wp
  atlant_moc(:,:) = 0.0_wp
    
  ! with all cells no sync is necessary
  !owned_cells => p_patch%cells%owned
  dom_cells   => p_patch%cells%in_domain

  !write(81,*) 'MOC: datetime:',datetime

  DO jk = 1, n_zlev   !  not yet on intermediate levels
    DO jb = dom_cells%start_block, dom_cells%end_block
      CALL get_index_range(dom_cells, jb, i_startidx, i_endidx)
      DO jc = i_startidx, i_endidx
        IF ( v_base%lsm_oce_c(jc,jk,jb) <= sea_boundary ) THEN
   
          ! lbrei: corresponding latitude row of 1 deg extension
          !       1 south pole
          !     180 north pole
          z_lat = p_patch%cells%center(jc,jb)%lat
          z_lat_deg = z_lat*rad2deg
          lbrei = NINT(90.0_wp + z_lat_deg)
          lbrei=MAX(lbrei,1)
          lbrei=MIN(lbrei,180)

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
            lbrei=MAX(lbrei,1)
            lbrei=MIN(lbrei,180)
      
            global_moc(lbrei,jk) = global_moc(lbrei,jk) - &
              &                    p_patch%cells%area(jc,jb) * rho_ref * w(jc,jk,jb) / &
              &                    REAL(2*jbrei + 1, wp)
      
         ! horizontal transports
         !  IF ( k == 1 ) THEN
         !    global_hfl(lbrei) = global_hfl(lbrei) &
         !         - area(i, j) * flum(i, j) / REAL(2*jbrei + 1, wp)
         !    global_wfl(lbrei) = global_wfl(lbrei) &
         !         - area(i, j) * pem(i, j) / REAL(2*jbrei + 1, wp)
         !  END IF

            IF (v_base%basin_c(jc,jb) <2) THEN

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
  END DO

  !IF (p_pe==p_io) THEN
    DO lbr=179,1,-1   ! fixed to 1 deg meridional resolution

        global_moc(lbr,:)=global_moc(lbr+1,:)+global_moc(lbr,:)
        pacind_moc(lbr,:)=pacind_moc(lbr+1,:)+pacind_moc(lbr,:)
        atlant_moc(lbr,:)=atlant_moc(lbr+1,:)+atlant_moc(lbr,:)

    END DO

    ! write out in extra format - integer*8
    idate=datetime%month*1000000+datetime%day*10000+datetime%hour*100+datetime%minute
    write(82,*) 'global MOC at iyear, idate:',datetime%year, idate
    DO jk=1,n_zlev
      i1=int(idate,8)
      i2=int(777,8)
      i3=int(v_base%zlev_i(jk),8)
      i4=int(180,8)
      write(77) i1,i2,i3,i4
      write(77) (global_moc(lbr,jk),lbr=1,180)
      i2=int(778,8)
      write(78) i1,i2,i3,i4
      write(78) (atlant_moc(lbr,jk),lbr=1,180)
      i2=int(779,8)
      write(79) i1,i2,i3,i4
      write(79) (pacind_moc(lbr,jk),lbr=1,180)

  !   write(82,*) 'jk=',jk
  !   write(82,'(1p10e12.3)') (global_moc(lbr,jk),lbr=1,180)
    END DO
  !END IF

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
SUBROUTINE calc_psi (p_patch, u, h, datetime)

  TYPE(t_patch), TARGET, INTENT(IN)  :: p_patch
  REAL(wp), INTENT(IN)               :: u(:,:,:)       ! zonal velocity at cell centers
  REAL(wp), INTENT(IN)               :: h(:,:)         ! elevation on cell centers
                                                       ! dims: (nproma,nlev,nblks_c)
! REAL(wp), INTENT(OUT)              :: psi(:,:)       ! horizontal stream function
  TYPE(t_datetime), INTENT(IN)       :: datetime
  !
  ! local variables
  ! INTEGER :: i

  INTEGER, PARAMETER ::  nlat = 180                    ! meridional dimension of regular grid
  INTEGER, PARAMETER ::  nlon = 360                    ! zonal dimension of regular grid

  ! smoothing area is 2*jsmth-1 lat/lon areas of 1 deg
  INTEGER, PARAMETER ::  jsmth = 3                  
  INTEGER            :: jb, jc, jk, i_startidx, i_endidx
  INTEGER            :: jlat, jlon, jlt, jln, jltx, jlnx, jsmth2
  INTEGER(i8)        :: idate, iextra(4)


  REAL(wp) :: z_lat_deg, z_lon_deg, z_lat_dist, delta_z, rsmth
  REAL(wp) :: z_uint(nproma,p_patch%nblks_c)            ! vertical integral of horizontal velocity
  REAL(wp) :: z_uint_reg(nlon,nlat)                     ! vertical integral on regular grid
  REAL(wp) :: psi(nlon,nlat)                            ! horizontal stream function

  TYPE(t_subset_range), POINTER :: all_cells, dom_cells
  
  !CHARACTER(len=max_char_length), PARAMETER :: routine = ('mo_oce_diagnostics:calc_psi')

  !-----------------------------------------------------------------------

  psi(:,:)        = 0.0_wp
  z_uint(:,:)     = 0.0_wp
  z_uint_reg(:,:) = 0.0_wp

  jsmth2          = 2*jsmth + 1
  rsmth           = REAL(jsmth2*jsmth2, wp)
   

  ! with all cells no sync is necessary
  all_cells => p_patch%cells%all
  dom_cells => p_patch%cells%in_domain

  ! (1) vertical integration of zonal velocity times vertical layer thickness [m/s*m]
  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
    DO jk = 1, n_zlev

      DO jc = i_startidx, i_endidx
        delta_z = v_base%del_zlev_m(jk)
        IF (jk == 1) delta_z = v_base%del_zlev_m(jk) + h(jc,jb)
        z_uint(jc,jb) = z_uint(jc,jb) - u(jc,jk,jb)*delta_z*v_base%wet_c(jc,jk,jb)
      END DO
    END DO
  END DO

  ! (2) distribute integrated zonal velocity (u*dz) on 1x1 deg grid

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

          z_uint_reg(jln,jlt) = z_uint_reg(jln,jlt) + z_uint(jc,jb) / rsmth

 ! 99 format('J lat=',f8.2,' lon=',f8.2,' jlat=',i4,' jlon=',i4,' lsm=',i3, &
 !      &    ' jlt=',i4,  ' jln=',i4,' uint=',1p10e12.3)
 ! 98 format(' lat=',f8.2,' lon=',f8.2,' jlat=',i4,' jlon=',i4,' lsm=',i3, &
 !      &    ' uint=',1p10e12.3)
 !    if ((jlat==101 .and. jlon==270) &
 !      & write(82,99) z_lat_deg,z_lon_deg,jlat,jlon,v_base%lsm_oce_c(jc,1,jb), &
 !      &              jlt,jln,z_uint_reg(jln,jlt)

        END DO
      END DO
 !    write(82,98) z_lat_deg,z_lon_deg,jlat,jlon,v_base%lsm_oce_c(jc,1,jb),z_uint_reg(jlon,jlat)

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

  psi(:,:) = z_uint_reg(:,:) * z_lat_dist * rho_ref


  ! write out in extra format - integer*8
  idate = int(datetime%month*1000000+datetime%day*10000+datetime%hour*100+datetime%minute,8)
  write(0,*) 'write global PSI at iyear, idate:',datetime%year, idate

  iextra(1) = int(idate,8)
  iextra(2) = int(780,8)
  iextra(3) = int(0,8)
  iextra(4) = int(nlon*nlat,8)

  write(80) (iextra(jb),jb=1,4)
  write(80) ((psi(jln,jlt),jln=1,nlon),jlt=1,nlat)

  do jlat=1,nlat
      write(82,*) 'jlat=',jlat
      write(82,'(1p10e12.3)') (psi(jlon,jlat),jlon=1,nlon)
  enddo

END SUBROUTINE calc_psi
!-------------------------------------------------------------------------  

END MODULE mo_oce_diagnostics
