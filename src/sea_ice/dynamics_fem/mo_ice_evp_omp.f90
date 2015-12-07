!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!
! Contains: Three routines of EVP dynamics. The driving routine is EVPdynamics.
! 2D indices are used to distinguish between boundary and internal nodes.
!
!
! Vladimir: added omp parallelization; rhs_a, rhs_m, rhs_mis are precalculated now
! The "original" AWI version of this solver is available in mo_ice_evp_old
!----------------------------
#include "omp_definitions.inc"
!----------------------------
!============================================================================
module mo_ice_evp
  !
  ! Arrays defined here are used to keep mesh information
  !
  use mo_ice_elements
  use mo_ice_mesh
  use mo_ice
  use mo_ice_parsup

  USE mo_kind,    ONLY: wp

  IMPLICIT NONE

  PUBLIC :: init_evp_solver_coeffs
  PUBLIC :: EVPdynamics

  PRIVATE :: stress_tensor
  PRIVATE :: stress2rhs
  PRIVATE :: precalc4rhs
  PRIVATE :: index_si_elements

  PRIVATE
  ! some aggregated parameters used in the EVP solver
  REAL(wp):: val3=1.0_wp/3.0_wp
  REAL(wp):: vale, dte, det1, det2
  REAL(wp):: ax, ay

CONTAINS

!===================================================================

subroutine index_si_elements
! Replaces "if" checking of whether sea ice is actually present in a given cell/node

    IMPLICIT NONE
    INTEGER :: elem,row, elnodes(3), i
    REAL(wp):: aa
    ! Temporary variables/buffers
    INTEGER :: buffy_array(myDim_elem2D)

!    ! count nodes with ice
    buffy_array = 0._wp
    si_nod2D = 0

    DO i=1, myDim_nod2D
        row=myList_nod2D(i)

         aa = m_ice(row)*a_ice(row)

         IF (aa > 0._wp) THEN
            si_nod2D = si_nod2D + 1
            buffy_array(si_nod2D)=row
         ENDIF
    ENDDO

    if (allocated(si_idx_nodes)) deallocate(si_idx_nodes)
    allocate(si_idx_nodes(si_nod2D))
    si_idx_nodes=buffy_array(1:si_nod2D)

    ! count elements with ice
    buffy_array = 0._wp
    si_elem2D = 0

    DO i=1,myDim_elem2D
         elem=myList_elem2D(i)
         elnodes=elem2D_nodes(:,elem)

         aa=product(m_ice(elnodes))*product(a_ice(elnodes))

         IF (aa > 0._wp) THEN
            si_elem2D = si_elem2D + 1
            buffy_array(si_elem2D)=elem
         ENDIF
    ENDDO

    if (allocated(si_idx_elem)) deallocate(si_idx_elem)
    allocate(si_idx_elem(si_elem2D))
    si_idx_elem=buffy_array(1:si_elem2D)

end subroutine index_si_elements
!===================================================================

subroutine precalc4rhs
! Some of the quantities used in the solver do not change
! with the subcycles. Hence, can be precalculated and stored.
! Those are rhs_a, rhs_m, mass

  use mo_physical_constants,  ONLY: rhoi, rhos

IMPLICIT NONE
INTEGER      :: row, elem, elnodes(3), nodels(6), k, i
REAL(wp) :: mass, aa
REAL(wp) :: cluster_area,elevation_elem(3),dx(3),dy(3)
REAL(wp) :: da(myDim_elem2D), dm(myDim_elem2D) ! temp storage for element-wise contributions from SSH

!ICON_OMP_PARALLEL

!ICON_OMP_DO PRIVATE(i,row) SCHEDULE(static,4)
 DO i=1, myDim_nod2D
     row=myList_nod2D(i)
     rhs_a(row)=0.0_wp    ! these are used as temporal storage here
     rhs_m(row)=0.0_wp    ! for the contribution due to ssh
     rhs_mis(row)=0.0_wp

     rhs_u(row)=0.0_wp    ! these will be reinitialized at every subcycling iteration for non-ice-free nodes
     rhs_v(row)=0.0_wp    ! but fill all nodal vals with zeros here as well
 END DO
!ICON_OMP_END_DO

!ICON_OMP_DO PRIVATE(i,elem) SCHEDULE(static,4)
 DO i=1, myDim_elem2D
     elem=myList_elem2D(i)
     da(elem)=0.0_wp    ! initialize
     dm(elem)=0.0_wp    !
 END DO
!ICON_OMP_END_DO

!ICON_OMP_DO PRIVATE(i,elem,elnodes,aa,dx,dy,elevation_elem) SCHEDULE(static,4)
DO i=1,si_elem2D
     elem=si_idx_elem(i)
     elnodes=elem2D_nodes(:,elem)

     dx=bafux(:,elem)
     dy=bafuy(:,elem)
     elevation_elem=elevation(elnodes)

     ! use rhs_m and rhs_a for storing the contribution from elevation:
     aa=9.81_wp*voltriangle(elem)/3.0_wp
     da(elem)=-aa*sum(dx*elevation_elem)
     dm(elem)=-aa*sum(dy*elevation_elem)
 end do
!ICON_OMP_END_DO

!ICON_OMP_DO  PRIVATE(i,row,nodels,k,elem)  SCHEDULE(static,4)
 DO i=1,si_nod2D
    row=si_idx_nodes(i)
    nodels=nod2D_elems(:,row)

     DO k=1,6
        elem=nodels(k)
        IF (elem > 0) THEN
         ! use rhs_m and rhs_a for storing the contribution from elevation:
        rhs_a(row)=rhs_a(row)+da(elem)
        rhs_m(row)=rhs_m(row)+dm(elem)
        END IF
     END DO
END DO
!ICON_OMP_END_DO

!ICON_OMP_DO  PRIVATE(i,row,cluster_area,mass) SCHEDULE(static,4)
  DO i=1,si_nod2D
     row=si_idx_nodes(i)
     cluster_area=lmass_matrix(row)
     mass=cluster_area*(m_ice(row)*rhoi+m_snow(row)*rhos)

     rhs_a(row) = rhs_a(row)/cluster_area
     rhs_m(row) = rhs_m(row)/cluster_area
     rhs_mis(row) = mass
  END DO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

end subroutine precalc4rhs
!===================================================================

subroutine init_evp_solver_coeffs
! Calculates coefficients which are used for calculations in stress_tensor
! Called once during the initialization step at ice_init_fem

  USE mo_run_config,          ONLY: dtime
  USE mo_sea_ice_nml,         ONLY: evp_rheol_steps, Tevp_inv, theta_io
  USE mo_physical_constants,  ONLY: ellipse

  val3=1.0_wp/3.0_wp
  vale=1.0_wp/(ellipse**2)

  dte=dtime/(1.0_wp*REAL(evp_rheol_steps,wp))
  det1=1.0_wp+0.5_wp*Tevp_inv*dte
  det2=1.0_wp+0.5_wp*Tevp_inv*dte*ellipse**2    !RTSD corrected 8.3.2006
                                              ! There is error in CICE
                          ! manual.
  det1=1.0_wp/det1
  det2=1.0_wp/det2

  ! theta_io should be zero - set in ice_main.f90
  ax=cos(theta_io)
  ay=sin(theta_io)

end subroutine init_evp_solver_coeffs
!===================================================================

subroutine stress_tensor(elem)
! EVP rheology implementation. Computes stress tensor components based on ice 
! velocity field. They are stored as elemental arrays (sigma11, sigma22 and
! sigma12). 

  USE mo_sea_ice_nml,         ONLY: delta_min, Tevp_inv
  USE mo_physical_constants,  ONLY: Pstar, c_pressure

implicit none

    INTEGER,                  INTENT(IN)     :: elem ! element index for which the stress tensor is calculated

    REAL(wp)   :: eps11, eps12, eps22, pressure, delta, delta_inv!, aa
    integer    :: elnodes(3)
    REAL(wp)   :: asum, msum, dx(3), dy(3)
    REAL(wp)   :: r1, r2, r3, si1, si2
    REAL(wp)   :: zeta, usum, vsum

! ATTENTION: the rows commented with !metrics contain terms due to 
! differentiation of metrics. 

     elnodes=elem2D_nodes(:,elem)

     dx=bafux(:,elem)
     dy=bafuy(:,elem)

!    meancos=sin_elem2D(elem)/cos_elem2D(elem)/earth_radius !is precalculated and stored as metrics_elem2D(elem)
     vsum=sum(v_ice(elnodes))                           !metrics   
     usum=sum(u_ice(elnodes))                           !metrics

      ! ===== Deformation rate tensor on element elem:
     eps11=sum(dx*u_ice(elnodes))
     eps11=eps11-val3*vsum*metrics_elem2D(elem)                !metrics
     eps22=sum(dy*v_ice(elnodes))
     eps12=0.5_wp*sum(dy*u_ice(elnodes) + dx*v_ice(elnodes))
     eps12=eps12+0.5_wp*val3*usum*metrics_elem2D(elem)          !metrics
      ! ===== moduli:
     delta=(eps11**2+eps22**2)*(1.0_wp+vale)+4.0_wp*vale*eps12**2 + &
            2.0_wp*eps11*eps22*(1.0_wp-vale)
     delta=sqrt(delta)
     msum=sum(m_ice(elnodes))*val3
     asum=sum(a_ice(elnodes))*val3
     
      ! ===== Hunke and Dukowicz c*h*p*
     pressure=pstar*(msum)*exp(-c_pressure*(1.0_wp-asum))
      ! =======================================
      ! ===== Here the EVP rheology piece starts
      ! =======================================
     pressure=0.5_wp*pressure
      ! ===== viscosity zeta should exceed zeta_min
      ! (done via limiting delta from above)
     !if(delta>pressure/zeta_min) delta=pressure/zeta_min
      ! ===== if viscosity is too big, it is limited too
      ! (done via correcting delta_inv)
     delta_inv=1.0_wp/max(delta,delta_min)
       
     !pressure=pressure*delta*delta_inv    ! Limiting pressure --- may not 
     !                                     ! be needed. Should be tested   
				     
      ! ===== Limiting pressure/Delta  (zeta): still it may happen that zeta is too
      ! large in regions with fine mesh so that CFL criterion is violated.
        				     
      zeta=pressure*delta_inv
      ! This place was introduced by Hunke, but seemingly is not used 
      ! in the current CICE. We artificially increase Clim_evp so
      ! that almost no limiting happens here
      !if (zeta>Clim_evp*voltriangle(elem)) then
      !zeta=Clim_evp*voltriangle(elem)
      !end if 
      
      pressure=pressure*Tevp_inv
      zeta=zeta*Tevp_inv
      				     
     r1=zeta*(eps11+eps22) - pressure
     r2=zeta*(eps11-eps22)
     r3=zeta*eps12
     si1=sigma11(elem)+sigma22(elem)
     si2=sigma11(elem)-sigma22(elem)
     
     si1=det1*(si1+dte*r1)
     si2=det2*(si2+dte*r2)
     sigma12(elem)=det2*(sigma12(elem)+dte*r3)
     sigma11(elem)=0.5_wp*(si1+si2)
     sigma22(elem)=0.5_wp*(si1-si2)
   
end subroutine stress_tensor

!===================================================================
subroutine stress2rhs(row)
! EVP implementation:
! Computes the divergence of stress tensor and puts the result into the
! rhs vectors 

IMPLICIT NONE

INTEGER  :: row, elem, nodels(6), k
REAL(wp) :: dx(6), dy(6)

     nodels=nod2D_elems(:,row)

     dx=bafux_nod(:,row)
     dy=bafuy_nod(:,row)

     ! initialize
     rhs_u(row)=0.0_wp
     rhs_v(row)=0.0_wp

     DO k=1,6
        elem=nodels(k)
        IF (elem > 0) THEN

        rhs_u(row)=rhs_u(row) - voltriangle(elem) * &
             (sigma11(elem)*dx(k)+sigma12(elem)*(dy(k)) &
             +sigma12(elem)*val3*metrics_elem2D(elem))                          !metrics
        rhs_v(row)=rhs_v(row) - voltriangle(elem) * &
             (sigma12(elem)*dx(k)+sigma22(elem)*dy(k) &
             -sigma11(elem)*val3*metrics_elem2D(elem))
        END IF
     END DO

     !  Vladimir: rhs_a and rhs_m are calculated in precalc4rhs
     rhs_u(row)=rhs_u(row)/rhs_mis(row) + rhs_a(row)
     rhs_v(row)=rhs_v(row)/rhs_mis(row) + rhs_m(row)

end subroutine stress2rhs
!===================================================================

subroutine EVPdynamics
! EVP implementation. Does cybcycling and boundary conditions.  
  USE mo_sea_ice_nml,           ONLY: evp_rheol_steps
  USE mo_physical_constants,    ONLY: rhoi, rhos, Cd_io, rho_ref

  USE mo_ice_fem_init,         ONLY: exchange_nod2D

IMPLICIT NONE
integer     :: shortstep
REAL(wp)    ::  drag, inv_mass, det, umod, rhsu, rhsv
integer     ::  i,j

! index elements/nodes where sea ice is present for faster loops
    call index_si_elements
! precalculate several arrays that do not change during subcycling
    call precalc4rhs

 DO shortstep=1, evp_rheol_steps
     ! ===== Boundary conditions
     do j=1, myDim_nod2D+eDim_nod2D
        i=myList_nod2D(j)
        if(index_nod2D(i)==1) then
        u_ice(i)=0.0_wp
        v_ice(i)=0.0_wp
        end if
     end do

!ICON_OMP_PARALLEL
!ICON_OMP_DO        PRIVATE(i)  SCHEDULE(static,4)
     DO i=1,si_elem2D
        call stress_tensor(si_idx_elem(i))
     ENDDO
!ICON_OMP_END_DO

!ICON_OMP_DO        PRIVATE(i) ICON_OMP_DEFAULT_SCHEDULE
     DO i=1,si_nod2D
         call stress2rhs(si_idx_nodes(i))
     END DO
!ICON_OMP_END_DO

!ICON_OMP_DO        PRIVATE(j,i,inv_mass,umod,drag,rhsu,rhsv,det) ICON_OMP_DEFAULT_SCHEDULE
     DO j=1,si_nod2D
        i=si_idx_nodes(j)
      if (index_nod2D(i)>0) CYCLE          ! Skip boundary nodes
      if (a_ice(i) > 0.01_wp) then             ! If ice is present, update velocities
       inv_mass=(rhoi*m_ice(i)+rhos*m_snow(i))/a_ice(i)
       inv_mass=max(inv_mass, 9.0_wp)        ! Limit the weighted mass
                                           ! if it is too small
       inv_mass=1.0_wp/inv_mass

       umod=sqrt((u_ice(i)-u_w(i))**2+(v_ice(i)-v_w(i))**2)
       drag=Cd_io*umod*rho_ref*inv_mass

       rhsu=u_ice(i)+dte*(drag*(ax*u_w(i)-ay*v_w(i))+inv_mass*stress_atmice_x(i)+rhs_u(i))
       rhsv=v_ice(i)+dte*(drag*(ax*v_w(i)+ay*u_w(i))+inv_mass*stress_atmice_y(i)+rhs_v(i))

       det=(1._wp+ax*drag*dte)**2+(dte*coriolis_nod2D(i)+dte*ay*drag)**2
       det=1.0_wp/det
       u_ice(i)=det*((1.0_wp+ax*drag*dte)*rhsu+dte*(coriolis_nod2D(i)+ay*drag)*rhsv)
       v_ice(i)=det*((1.0_wp+ax*drag*dte)*rhsv-dte*(coriolis_nod2D(i)+ay*drag)*rhsu)
      ! else                           ! Set ice velocity equal to water velocity
      ! u_ice(i)=u_w(i)
      ! v_ice(i)=v_w(i)
       end if
     end do
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

     call exchange_nod2D(u_ice)
     call exchange_nod2D(v_ice)
 END DO

end subroutine EVPdynamics
!===================================================================

end module mo_ice_evp
