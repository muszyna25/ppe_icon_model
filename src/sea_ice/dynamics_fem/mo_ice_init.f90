!
! array allocation + matrix setup routines+field initialization+io.
! Initialization and IO are just very simple versions (when coupled to 
! FESOM, ice uses FESOM machinery) that allow one to run test cases.
!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.

! Matrices are needed (i) if VP solver is used, (ii) in advection
! routines and (iii) for partitioning

! There are, in principle, 3 matrices which share the same 
! structure stored in icestiff. We use icestiff%values for both implicit 
! advection schemes and VP algorithms. Since they are re-assembled 
! at each time step, there is no problem with sharing the place for
! values.
! Additionally, there is mass matrix which is used for FCT advection 
! scheme.  It has the same structure, so we use mass_matrix 
! of the size of icestiff%values only to keep  the values.
!==============================================================================
module mo_ice_init

  USE mo_kind,         ONLY: wp
  USE mo_run_config,   ONLY: dtime
  USE mo_sea_ice_nml,  ONLY: Tevp_inv
  USE mo_ice_elements
  USE mo_ice
  USE mo_ice_mesh
  USE mo_ice_parsup
  USE mo_ice_evp_old, ONLY: init_evp_solver_coeffs_old
  USE mo_ice_evp, ONLY: init_evp_solver_coeffs
!  use mo_ice_atm_forcing
!  USE mo_physical_constants,         ONLY: grav, earth_radius

  IMPLICIT NONE

  PUBLIC :: ice_init_fem

  PUBLIC :: array_setup_ice
  PUBLIC :: icestiff_matrix
  PRIVATE :: ice_mass_matrix_fill

CONTAINS

!=================================================================
subroutine ice_init_fem

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
  !
  ! ########
  !

  ! ================ DO not change
  Tevp_inv=3.0_wp/dtime
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
  ! Vladimir: init some coeffs for EVP solver
  call init_evp_solver_coeffs_old
  call init_evp_solver_coeffs

  ! write(*,*) 'fields are initialized'

end subroutine ice_init_fem

!=================================================================
subroutine array_setup_ice
!
! inializing sea ice model 
!

integer   :: k

!---------------------------------------------------------------------------------------------------
! The Coriolis term is now calculated in ice_fem_grid_init (modified by Vladimir Lapin 11/08/2015)
!---------------------------------------------------------------------------------------------------
! allocate(coriolis_nod2D(nod2D))
! coriolis_nod2D=2._wp*omega*sin(coord_nod2D(2,:))

 allocate(col_pos(nod2d))


! Allocate memory for atmospheric forcing 
! allocate(shortwave(nod2D), longwave(nod2D), Tair(nod2D), Tdew(nod2D))
! allocate(Pair(nod2D), precipitation(nod2D), evaporation(nod2D))
! allocate(u_wind(nod2D), v_wind(nod2D))
! allocate(shum(nod2d), cloudiness(nod2d))


! Allocate memory for variables of ice model      
 allocate(u_ice(nod2D), v_ice(nod2D))
 allocate(m_ice(nod2D), a_ice(nod2D), m_snow(nod2D))
 allocate(sigma11(elem2D), sigma12(elem2D), sigma22(elem2D))
 allocate(rhs_m(nod2D), rhs_a(nod2D), rhs_u(nod2D), rhs_v(nod2D))
 allocate(rhs_mis(nod2D))

 rhs_m=0._wp
 rhs_mis=0._wp
 rhs_a=0._wp
 rhs_u=0._wp
 rhs_v=0._wp
 sigma11=0._wp
 sigma22=0._wp
 sigma12=0._wp

! Allocate memory used for coupling (Partly used for input, and partly
! for output of information) 
 allocate(S_oc_array(nod2D), T_oc_array(nod2D))
 allocate(fresh_wa_flux(nod2D), net_heat_flux(nod2d))
 allocate(stress_atmice_x(nod2D), stress_atmice_y(nod2D))    
 allocate(stress_atmoce_x(nod2D), stress_atmoce_y(nod2D))    
 allocate(stress_iceoce_x(nod2D), stress_iceoce_y(nod2D))    
 allocate(elevation(nod2D))         ! =ssh  of ocean        
 allocate(u_w(nod2D), v_w(nod2D))   ! =uf and vf of ocean at surface nodes

! Fill in  the mass matrix    
 call ice_mass_matrix_fill
! Create lumped mass matrix
 allocate(lmass_matrix(nod2D))
DO k=1,nod2D
 lmass_matrix(k)=sum(mass_matrix(icestiff%rowptr(k):icestiff%rowptr(k+1)-1))
END DO
 ! write(*,*) 'Stiffness matrix is allocated'   
! Einar: ICON does advection
!if(ice_advection==ice_FCT) then
!call fem_fct_ice_init      !##### FCT advection scheme allocations
!end if
end subroutine array_setup_ice
!==========================================================================

SUBROUTINE icestiff_matrix()
!
! Sets the structure of the stiffness matrix for BE ice advection/VP problems
!

INTEGER                       :: k


! a) We must compute the number of potentially non-zero entries
! b) Allocate space
! c) Fill in the column index array with proper values 


! a)
icestiff%dim=nod2D
allocate(icestiff%rowptr(icestiff%dim+1))
icestiff%rowptr(1)=1
DO k=1,nod2D
  icestiff%rowptr(k+1)=icestiff%rowptr(k)+nghbr_nod2D(k)%nmb
END DO
icestiff%nza=icestiff%rowptr(nod2D+1)-1
! b)
allocate(icestiff%colind(icestiff%nza))
allocate(icestiff%values(icestiff%nza))
allocate(mass_matrix(icestiff%nza))
icestiff%values=0.0_wp

! c)
        ! ===== FIND colind ======
DO k=1,nod2D
  icestiff%colind(icestiff%rowptr(k):icestiff%rowptr(k+1)-1)= &
    nghbr_nod2D(k)%addresses
END DO
END SUBROUTINE icestiff_matrix

!=========================================================================
SUBROUTINE ice_mass_matrix_fill
! Computes entries of mass matrix used to speed up stress2rhs and in the
! implementation of the FCT scheme

INTEGER                       :: col, elem, elnodes(3), offset
INTEGER                       :: ipos, row, q, n, i

  mass_matrix =0.0_wp
  DO i=1, myDim_elem2D   
     elem=myList_elem2D(i)
     elnodes=elem2D_nodes(:,elem)   
     DO n=1,3
        row=elnodes(n)
        IF(part(row).NE.mype) CYCLE             
        ! Global-to-local neighbourhood correspondence  
        DO q=1,nghbr_nod2D(row)%nmb
           col_pos(nghbr_nod2D(row)%addresses(q))=q
        END DO 
        offset=icestiff%rowptr(row)-1
        DO q=1,3 
           col=elnodes(q)
           ipos=offset+col_pos(col)
           mass_matrix(ipos)=mass_matrix(ipos)+voltriangle(elem)/12.0_wp
           IF(q==n) THEN                     
             mass_matrix(ipos)=mass_matrix(ipos)+voltriangle(elem)/12.0_wp
           END IF
        END DO
     END DO
  END DO

END SUBROUTINE ice_mass_matrix_fill

! ============================================================================ 
!subroutine init_fields
!!
!! Simple initialization for a box model to test the dynamical part.
!! No thermodinamics is initialized here
!!
!
!REAL(wp)  :: xmin, xmax, ymin, ymax, Lx, Ly, yc, xc, meanf
!integer       :: n
!
!  coriolis_nod2D=1.4e-4_wp  ! redefines Coriolis
!
!   ! Set initial thickness and area coverage:
!   m_ice=2.0_wp
!   m_snow=0.0_wp
!   u_ice=0.0_wp
!   v_ice=0.0_wp
!   stress_atmice_x=0.0_wp
!   stress_atmice_y=0.0_wp
!   ! a_ice is defined later
!
!
!   ! Set ocean velocity (stationary in time):
!   xmin=minval(coord_nod2d(1,:))
!   xmax=maxval(coord_nod2d(1,:))
!   ymin=minval(coord_nod2d(2,:))
!   ymax=maxval(coord_nod2d(2,:))
!   Lx=xmax-xmin
!   Ly=ymax-ymin
!   xc=0.5_wp*(xmax+xmin)
!   yc=0.5_wp*(ymax+ymin)
!   DO n=1, nod2D
!   a_ice(n)=(coord_nod2d(1,n)-xmin)/Lx
!   u_w(n)=0.1_wp*(2._wp*(coord_nod2d(2,n)-ymin)-Ly)/Ly
!   v_w(n)=-0.1_wp*(2._wp*(coord_nod2d(1,n)-xmin)-Lx)/Lx
!   END DO
!   m_ice=m_ice*a_ice
!
!   ! Elevation computed approximately, from the geostrophy:
!   meanf= 1.4e-4_wp*earth_radius   !2*omega*sin(yc)*r_earth
!   DO n=1, nod2d
!      elevation(n)=-0.1_wp*meanf/grav *((coord_nod2d(2,n)-ymin)**2/Ly- &
!                    (coord_nod2d(2,n)-ymin)+ &
!                    (coord_nod2d(1,n)-xmin)**2/Lx -&
!                    (coord_nod2d(1,n)-xmin))
!   END DO
!
!end subroutine init_fields
! =============================================================================
!Subroutine update_forcing(step)
!!
!! Here only simple wind variability is introduced
!!
!
!  USE mo_math_constants,         ONLY: pi
!  USE mo_run_config,             ONLY: dtime
!  USE mo_physical_constants,     ONLY: Cd_ia
!  USE mo_kind,                   ONLY: wp
!
!  use mo_ice_elements
!  use mo_ice_mesh
!  use mo_ice
!  use mo_ice_parsup
!
!  use mo_ice_atm_forcing
!
!IMPLICIT NONE
!REAL(wp), parameter :: rhoair=  1.3_wp    ! Air density
!
!REAL(wp)  :: xmin, xmax, ymin, ymax, Lx, Ly, td
!integer       :: step, n
!   ! Set wind velocity (stationary in time):
!   xmin=minval(coord_nod2d(1,:))
!   xmax=maxval(coord_nod2d(1,:))
!   ymin=minval(coord_nod2d(2,:))
!   ymax=maxval(coord_nod2d(2,:))
!   Lx=xmax-xmin
!   Ly=ymax-ymin
!   td=4._wp*3600._wp*24.0_wp
!
!   DO n=1, nod2D
!   u_wind(n)=5.0_wp+(sin(2._wp*pi*step*dtime/td)-3.0_wp)*sin(2._wp*pi*(coord_nod2d(1,n)-xmin)/Lx) &
!                *sin(pi*(coord_nod2d(2,n)-ymin)/Ly)
!
!   v_wind(n)=5.0_wp+(sin(2._wp*pi*step*dtime/td)-3.0_wp)*sin(2._wp*pi*(coord_nod2d(2,n)-ymin)/Ly) &
!                *sin(pi*(coord_nod2d(1,n)-xmin)/Lx)
!   END DO
!   ! wind to stress:
!
!   stress_atmice_x = rhoair*Cd_ia*sqrt(u_wind**2+v_wind**2)*u_wind
!   stress_atmice_y = rhoair*Cd_ia*sqrt(u_wind**2+v_wind**2)*v_wind
!
!end subroutine update_forcing
!===================================================================

! Einar: ICON takes care of output and restarts
!subroutine ice_out
!
!USE MESH
!USE ELEMENTS
!USE ICEPARAM
!USE ice
!use parsup
!
!IMPLICIT NONE
!character*20   :: out_type='rewind'  !'append'  or 'rewind'
!INTEGER :: i
!
! call broadcast_nod2d(u_ice)
! call broadcast_nod2d(v_ice)
! call broadcast_nod2d(m_ice)
! call broadcast_nod2d(a_ice)
! call broadcast_nod2d(m_snow)
!
! if (mype==0) then     ! Use PE 0 for output
!
! write(*,*) 'writing results'
!
! open(17,file='u_ice.out',position=trim(out_type))
! DO i=1,nod2D
! write(17, '(1e12.6)') u_ice(i)
! END DO 
! DO i=1,nod2D
! write(17, '(1e12.6)') v_ice(i)
! END DO 
! close(17)
! 
! open(15,file='m_ice.out',position=trim(out_type))  
! DO i=1,nod2D
!    write(15,'(1e9.3)')  m_ice(i)
! END DO  
! close(15)
!
!
! open(16,file='a_ice.out',position=trim(out_type))  
! DO i=1,nod2D
! write(16, '(1e9.3)') a_ice(i)
! END DO 
! close(16)
!return 
! open(16,file='m_snow.out',position=trim(out_type))
! DO i=1,nod2D
! !write(16, '(1f6.3)') m_snow(i)
! write(16, *)  m_snow(i)
! END DO 
! close(16)
! end if
!
! end subroutine ice_out
!======================================================================
!subroutine write_restart_ice
!  use ice
!  use PARSUP
!  implicit none
!  integer :: i
!
!  if(mype/=0) return
!  write(*,*) 'Writing restart'
!  !
!  ! Write restart
!  open(25,file='restart_ice.out', form='unformatted')
!     write(25) m_ice, m_snow, a_ice, u_ice, v_ice
!  close(25)
!end subroutine write_restart_ice
!
!----------------------------------------------------------------------------
!
!subroutine read_restart_ice
!  use ice
!  implicit none
!  integer :: i
!
!  write(*,*) 'Reading restart'
!  !
!  ! Read restart
!  open(25,file='restart_ice.out', form='unformatted')
!  read(25)  m_ice, m_snow, a_ice, u_ice, v_ice
!  close(25)
!end subroutine read_restart_ice
!
!----------------------------------------------------------------------------
!
end module mo_ice_init
