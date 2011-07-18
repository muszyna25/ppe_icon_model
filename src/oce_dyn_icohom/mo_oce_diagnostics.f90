!>
!! Contains basic diagnostics for ICON ocean model.
!! 
!! 
!! @par Revision History
!!  Developed  by Peter Korn,       MPI-M (2011/02)
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
USE mo_kind,                      ONLY: wp
!USE mo_mpi,                       ONLY: p_pe, p_io
USE mo_math_utilities,            ONLY: t_cartesian_coordinates!, gc2cc
USE mo_impl_constants,            ONLY: sea_boundary,sea, &
  &                                     min_rlcell, min_rledge, min_rlcell, &
  &                                     max_char_length
USE mo_ocean_nml,                ONLY: n_zlev, no_tracer,&! toplev, &
                                    &   ab_const, ab_beta, ab_gam, iswm_oce, idisc_scheme
USE mo_dynamics_config,          ONLY: nold,nnew
USE mo_parallel_config,  ONLY: nproma
USE mo_run_config,                ONLY: dtime, nsteps
USE mo_physical_constants,        ONLY: grav!, re
USE mo_oce_state,                 ONLY: t_hydro_ocean_state, t_hydro_ocean_diag, &
  &                                     set_lateral_boundary_values
USE mo_model_domain,              ONLY: t_patch
USE mo_exception,                 ONLY: message, finish!, message_text
USE mo_loopindices,               ONLY: get_indices_c, get_indices_e !, get_indices_v
USE mo_oce_math_operators,        ONLY: div_oce, grad_fd_norm_oce, grad_fd_norm_oce_2d,&
                                  & height_related_quantities
USE mo_oce_physics,               ONLY: t_ho_params
USE mo_oce_forcing,               ONLY: t_ho_sfc_flx
USE mo_interpolation,             ONLY: t_int_state
USE mo_scalar_product,            ONLY: calc_scalar_product_for_veloc
USE mo_interpolation,             ONLY: t_int_state, rbf_vec_interpol_edge,       &
                                        rbf_vec_interpol_cell!verts2edges_scalar
USE mo_oce_ab_timestepping,       ONLY: calc_vert_velocity

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
TYPE(t_ho_sfc_flx), INTENT(INOUT)             :: p_sfc_flx
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
INTEGER :: i_no_t, i

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

  DO jk=1,n_zlev

    delta_z = p_patch%patch_oce%del_zlev_m(jk)

    DO jb = i_startblk_c, i_endblk_c
      CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
         &                             rl_start_c, rl_end_c)
      DO jc = i_startidx_c, i_endidx_c
        IF (jk == 1) THEN
         delta_z = p_patch%patch_oce%del_zlev_m(jk)&
                 & + p_os%p_prog(nold(1))%h(jc,jb)
        ENDIF

        prism_vol = p_patch%cells%area(jc,jb)*delta_z

        !Fluid volume 
        ptr_monitor%volume = ptr_monitor%volume + prism_vol

        !kinetic energy
        ptr_monitor%kin_energy = ptr_monitor%kin_energy+ p_os%p_diag%kin(jc,jk,jb)*prism_vol

        !Potential energy
        IF(jk==1)THEN
          z_w = (p_os%p_diag%w(jc,jk,jb)*p_os%p_prog(nold(1))%h(jc,jb)&
             & +p_os%p_diag%w(jc,jk+1,jb)*0.5_wp*p_patch%patch_oce%del_zlev_i(jk))&
             &/(0.5_wp*p_patch%patch_oce%del_zlev_i(jk)+p_os%p_prog(nold(1))%h(jc,jb))
        ELSEIF(jk<n_zlev)THEN
          z_w = (p_os%p_diag%w(jc,jk,jb)*p_patch%patch_oce%del_zlev_i(jk)&
             & +p_os%p_diag%w(jc,jk+1,jb)*p_patch%patch_oce%del_zlev_i(jk+1))&
             &/(p_patch%patch_oce%del_zlev_i(jk)+p_patch%patch_oce%del_zlev_i(jk+1))
        ENDIF 

        ptr_monitor%pot_energy = ptr_monitor%pot_energy&
        &+ grav*z_w* p_os%p_diag%rho(jc,jk,jb)* prism_vol

        !Tracer content
        DO i_no_t=1, no_tracer
          ptr_monitor%tracer_content(i_no_t) = ptr_monitor%tracer_content(i_no_t)&
          & + prism_vol*p_os%p_prog(nold(1))%tracer(jc,jk,jb,i_no_t)
        END DO

      END DO
    END DO
  END DO

  !divide by volume
  ptr_monitor%kin_energy = ptr_monitor%kin_energy/ptr_monitor%volume
  ptr_monitor%pot_energy = ptr_monitor%pot_energy/ptr_monitor%volume
  ptr_monitor%total_energy = ptr_monitor%pot_energy + ptr_monitor%kin_energy

  DO i_no_t=1, no_tracer
    ptr_monitor%tracer_content(i_no_t) = ptr_monitor%tracer_content(i_no_t)/ptr_monitor%volume
  END DO


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
    END DO
  END DO

ptr_monitor%pot_energy = ptr_monitor%pot_energy/ptr_monitor%volume
ptr_monitor%kin_energy = ptr_monitor%kin_energy/ptr_monitor%volume
ptr_monitor%total_energy = ptr_monitor%pot_energy + ptr_monitor%kin_energy

ENDIF

DO i=timestep,timestep
  write(*,*)'ACTUAL VALUES OF VOLUME NORMALIZED BY INITIAL VALUE:         ', i,&
  &oce_ts%oce_diagnostics(i)%volume/oce_ts%oce_diagnostics(0)%volume

IF(oce_ts%oce_diagnostics(1)%kin_energy/=0.0_wp)THEN
  write(*,*)'ACTUAL VALUES OF KINETIC ENERGY NORMALIZED BY INITIAL VALUE:  ',  i,&!timestep,&
  &oce_ts%oce_diagnostics(i)%kin_energy/ oce_ts%oce_diagnostics(1)%kin_energy
ENDIF

IF(oce_ts%oce_diagnostics(0)%pot_energy/=0.0_wp)THEN
  write(*,*)'ACTUAL VALUES OF POTENTIAL ENERGY NORMALIZED BY INITIAL VALUE:',i,&
  &oce_ts%oce_diagnostics(i)%pot_energy/ oce_ts%oce_diagnostics(0)%pot_energy
ENDIF

IF(oce_ts%oce_diagnostics(0)%total_energy/=0.0_wp)THEN
  write(*,*)'ACTUAL VALUES OF TOTAL ENERGY NORMALIZED BY INITIAL VALUE:     ',i,&
  &oce_ts%oce_diagnostics(i)%total_energy/ oce_ts%oce_diagnostics(0)%total_energy
ENDIF
END DO
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
SUBROUTINE construct_oce_diagnostics(p_patch, p_os, oce_ts)
!
TYPE(t_patch), TARGET, INTENT(in)             :: p_patch
TYPE(t_hydro_ocean_state), TARGET             :: p_os 
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
! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diagnostics:construct_oce_diagnostics')
!-----------------------------------------------------------------------
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


CALL height_related_quantities( p_patch, p_os)

IF(idisc_scheme==1)THEN
CALL calc_scalar_product_for_veloc( p_patch,                &
                                    & p_os%p_prog(nold(1))%vn,&
                                    & p_os%p_prog(nold(1))%vn,&
                                    & p_os%p_diag%h_e,        &
                                    & p_os%p_diag)
ELSE
! CALL rbf_vec_interpol_edge( p_os%p_prog(nold(1))%vn,&
!                           & p_patch,                &
!                           & p_int,                  &
!                           & p_os%p_diag%vt)
! CALL rbf_vec_interpol_cell( p_os%p_prog(nold(1))%vn,&
!                           & p_patch,&
!                           & p_int,&
!                           & p_os%p_diag%u,  &
!                           & p_os%p_diag%v)

!add calculation of kinetic energy

ENDIF

CALL calc_vert_velocity( p_patch, p_os)

!calculate initial values
  ptr_monitor =>oce_ts%oce_diagnostics(0)
  IF(iswm_oce/=1)THEN

    DO jk=1,n_zlev

      delta_z = p_patch%patch_oce%del_zlev_m(jk)

      DO jb = i_startblk_c, i_endblk_c
        CALL get_indices_c(p_patch, jb, i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c, &
         &                             rl_start_c, rl_end_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (jk == 1) THEN
           delta_z = p_patch%patch_oce%del_zlev_m(jk)&
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
             & +p_os%p_diag%w(jc,jk+1,jb)*0.5_wp*p_patch%patch_oce%del_zlev_i(jk))&
             &/(0.5_wp*p_patch%patch_oce%del_zlev_i(jk)+p_os%p_prog(nold(1))%h(jc,jb))
        ELSEIF(jk<n_zlev)THEN
          z_w = (p_os%p_diag%w(jc,jk,jb)*p_patch%patch_oce%del_zlev_i(jk)&
             & +p_os%p_diag%w(jc,jk+1,jb)*p_patch%patch_oce%del_zlev_i(jk+1))&
             &/(p_patch%patch_oce%del_zlev_i(jk)+p_patch%patch_oce%del_zlev_i(jk+1))
        ENDIF 

        ptr_monitor%pot_energy = ptr_monitor%pot_energy&
        &+ grav*z_w* p_os%p_diag%rho(jc,jk,jb)* prism_vol

        !Tracer content
        DO i_no_t=1, no_tracer
          ptr_monitor%tracer_content(i_no_t) = ptr_monitor%tracer_content(i_no_t)&
          & + prism_vol*p_os%p_prog(nold(1))%tracer(jc,jk,jb,i_no_t)
        END DO
      END DO
    END DO
  END DO

  !divide by volume
  ptr_monitor%kin_energy = ptr_monitor%kin_energy/ptr_monitor%volume
  ptr_monitor%pot_energy = ptr_monitor%pot_energy/ptr_monitor%volume
  ptr_monitor%total_energy = ptr_monitor%pot_energy + ptr_monitor%kin_energy

  DO i_no_t=1, no_tracer
    ptr_monitor%tracer_content(i_no_t) = ptr_monitor%tracer_content(i_no_t)/ptr_monitor%volume
  END DO


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
    END DO
  END DO

ptr_monitor%pot_energy = ptr_monitor%pot_energy/ptr_monitor%volume
ptr_monitor%kin_energy = ptr_monitor%kin_energy/ptr_monitor%volume
ptr_monitor%total_energy = ptr_monitor%pot_energy + ptr_monitor%kin_energy
ENDIF

write(*,*)'INITIAL VALUES OF VOLUME          :',oce_ts%oce_diagnostics(0)%volume
write(*,*)'INIIAL VALUES OF KINETIC ENERGY   :',oce_ts%oce_diagnostics(0)%kin_energy
write(*,*)'INITIAL VALUES OF POTENTIAL ENERGY:',oce_ts%oce_diagnostics(0)%pot_energy
write(*,*)'INITIAL VALUES OF TOTAL ENERGY    :',oce_ts%oce_diagnostics(0)%total_energy


!CALL message (TRIM(routine), 'end')
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

! CHARACTER(len=max_char_length), PARAMETER :: &
!        & routine = ('mo_oce_diagnostics:destruct_oce_diagnostics')
!-----------------------------------------------------------------------
   DO i=0,nsteps 
     DEALLOCATE(oce_ts%oce_diagnostics(i)%tracer_content)
   END DO
   DEALLOCATE(oce_ts%oce_diagnostics)
   DEALLOCATE(oce_ts)

!CALL message (TRIM(routine), 'end')
END SUBROUTINE destruct_oce_diagnostics
!-------------------------------------------------------------------------  



END MODULE mo_oce_diagnostics
