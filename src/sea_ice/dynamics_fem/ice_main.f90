! =====================
!  Standing alone sea ice
!  Version 2, based on version 1, with several new features added
!  Most important are true VP solver and FCT advection
!  Questions to S. Danilov (dynamics) and Q. Wang and R. Timmermann
! (thermodynamics). Many updates and corrections
!  to version 1 are due to R. Timmermann
! ======================

!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

! This file contails all modules, subroutine feim_step (ice model
! step) and main program. The modules here are derived from
! FESOM, the coupled version relies on full FESOM modules.


!=======================================================================
!subroutine feim_step
! 
! Sea ice model step
!
!use param
!use mesh
!use elements 
!use ice
!use iceparam
!use parsup
!use atm_forcing
!implicit none
!
!integer  :: time, globaltime, n, i
!integer  :: year, month, day, globalday
!
!REAL(wp) :: t0,t1, t2, t3
!t0=MPI_Wtime()
!
! ! ===== Dynamics
! if(ice_VP_rheology) then
! call VPdynamics
! else
! call EVPdynamics
! end if
! 
! t2=MPI_Wtime()     
!  
! ! ===== Advection part
! if(ice_advection==ice_FCT) then 
!         ! ===== FCT ice advection: 
!         call ice_TG_rhs
!         call fct_ice_solve    
! end if
! if(ice_advection==ice_BE) then 
!         call icestiff_fill2       ! Assemble the stiffness matrix of 
!                                   ! advection
!         call ice_rhs              ! compute rhs for the advection step
! 
!         call solveIce(12)         ! Advect m_ice
!         call solveIce(13)         ! Advect a_ice
!         call solveIce(14)         ! Advect m_snow
! end if  
!   
! call cut_off
! ! ===== Thermodynamic part
! !call thermodynamics
! t1=MPI_Wtime()
!if (mype==0) then
! write(*,*) 'FEIM step took ', t1-t0
! write(*,*) 'FEIM dynamics  ', t2-t0
!endif
!
!end subroutine feim_step
!==============================================================================
! Ice main routine 

!program ice_main
subroutine ice_init_fem_old

  use mo_ice_elements
  use mo_ice
  use mo_ice_parsup
  use mo_ice_init

!  USE mo_run_config,          ONLY: dtime
  USE mo_kind,    ONLY: wp
  USE mo_ice_mesh,            ONLY: set_fem_mesh

implicit none
!  integer      :: steps, step

  ! set the parameters

!  meshpath='/iblade/home/sdanilov/ice/mesh2/'
!  meshpath='/iblade/user/sdanilov/fv/mesh/box/'
!  meshpath='./'
!  dt=3600.0/2.0            ! time step
!  dt=dtime

  !
  ! ######## THIS piece was moved to namelist, mo_physical_constants
  !
!  evp_rheol_steps=120  ! the number of sybcycling steps in EVP
!  delta_min=2.0e-11_wp ! (1/s) Limit for minimum divergence (Hibler, Hunke
                       ! normally use 2.0e-9, which does much stronger
		       ! limiting; valid for both VP and EVP
!  Clim_evp= 61500.0_wp   !(kg/m^3) (see Hunke)
                       ! limits viscosity to be CFL stable in EVP
!  zeta_min= 4.0e+8_wp  !(kg/s), Minimum viscosity. Implemented in EVP, but
                       ! commented out as it does not lead to
                       ! correct physics.
!  theta_io=0._wp       !    0.436   ! ice/ocean rotation angle. Implemented in EVP, can be
                       ! added to VP
!  ice_gamma_fct=0.5_wp ! smoothing parameter in ice fct advection
!  ice_diff=10.0_wp     ! diffusion to stabilize ice advection
!  ice_VP_rheology=.false.    ! VP=.true., EVP=.false.
!  ice_VP_soltol=1.0e-4_wp
!  ice_advection=ice_BE      ! ice_FCT switches on FCT advection, and
                            ! ice_BE switches on Backward Euler advection
		       ! The real value of diffusion is scaled with element
		       ! area: diff=ice_diff*max(1, area/0.5e+8)
		       ! (10 km triangle); It is, however limited from
		       ! above: diff=min(2000,diff); One can change this scaling
		       ! if inapropriate in ice_stiff_fill2 (Backward Euler)
		       ! or ice_TG_rhs (FCT). Generally, one can use very small
		       ! values for FCT scheme to preserve sharp fronts in
		       ! ice thickness and area coverage.
!  steps= 24*10*2       ! Total number of steps
  !
  ! ########
  !

  ! ================ DO not change
!  Tevp_inv=3.0_wp/dtime
!  Clim_evp=Clim_evp*(evp_rheol_steps/dt)**2/Tevp_inv  ! This is combination
                                                       ! it always enters
  ! ================

! Einar: No parallelization yet
  !call par_init
  maxPENum=1
  npes=1

  ! ########
  ! Standing-alone mesh requires setting up the mesh
  ! In the coupled version this set of routines is called
  ! in the ocean model

  call set_fem_mesh
  !########
  call icestiff_matrix         ! We need it here because partit
                               ! called in set_par_support needs it
			       ! In coupled version, this is not needed
			       ! and can be moved to array_setup_ice
			       ! before call ice_mass_matrix_fill
  !########
! Einar: No parallelization yet
  call set_par_support

  ! This is the ice specific part. Instead of real io routines
  ! dynamical fields and forcing are now initialized and set by
  ! init_fields and update_forcing routines.
  ! ########
  call array_setup_ice
  ! ########
! Einar: We only need to initialize u and v
  u_ice = 0._wp; v_ice = 0._wp
!  call init_fields
! Einar: Handled by ICON
!  if(r_restart) then
!     call read_restart_ice
!     write(*,*) 'restart'
!  end if
!  write(*,*) 'fields are initialized'

! Einar: Initialization finished
!  do step=1, steps
!     if(mype==0)  write(*,*) 'DAY ', step*dt/(24*3600)
!      call update_forcing(step)
!
!     call feim_step
!     if((step/(24*3*5))*24*3*5==step) then
!     call ice_out
!     call write_restart_ice
!     end if
!     lfirst=.false.
!  end do
!  call ice_out
!  call write_restart_ice
!  call par_ex
end subroutine ice_init_fem_old
