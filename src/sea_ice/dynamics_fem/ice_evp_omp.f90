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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

!============================================================================
subroutine stress_tensor_omp
! EVP rheology implementation. Computes stress tensor components based on ice 
! velocity field. They are stored as elemental arrays (sigma11, sigma22 and
! sigma12). 

  use mo_ice_elements
  use mo_ice_mesh
  use mo_ice
  use mo_ice_parsup

  USE mo_run_config,          ONLY: dtime
  USE mo_sea_ice_nml,         ONLY: delta_min, evp_rheol_steps, Tevp_inv, theta_io
  use mo_physical_constants,  ONLY: Pstar, ellipse, c_pressure, earth_radius

  USE mo_kind,    ONLY: wp

implicit none

REAL(wp)   :: eps11, eps12, eps22, pressure, delta, aa
integer        :: elem, elnodes(3),i
REAL(wp)   :: val3, asum, msum, vale, dx(3), dy(3)
REAL(wp)   :: det1, det2, r1, r2, r3, si1, si2, dte 
REAL(wp)   :: zeta, delta_inv, meancos, usum, vsum

! ATTENTION: the rows commented with !metrics contain terms due to 
! differentiation of metrics. 

  val3=1.0_wp/3.0_wp
  vale=1.0_wp/(ellipse**2)
   
  dte=dtime/(1.0_wp*REAL(evp_rheol_steps,wp))
  det1=1.0_wp+0.5_wp*Tevp_inv*dte
  det2=1.0_wp+0.5_wp*Tevp_inv*dte*ellipse**2    !RTSD corrected 8.3.2006
                                              ! There is error in CICE 
					      ! manual. 
  det1=1.0_wp/det1
  det2=1.0_wp/det2
  
!ICON_OMP_PARALLEL_DO PRIVATE(i,elem,elnodes,aa,dx,dy,meancos, usum, vsum,  &
!ICON_OMP       eps11,eps12,eps22,delta,msum,asum,pressure,delta_inv,zeta,  &
!ICON_OMP       r1, r2, r3, si1, si2) ICON_OMP_DEFAULT_SCHEDULE
  do i=1,myDim_elem2D
     elem=myList_elem2D(i)     
     elnodes=elem2D_nodes(:,elem)
      ! ===== Check if there is ice on elem
        aa=product(m_ice(elnodes))*product(a_ice(elnodes))
        if (aa==0._wp) CYCLE     ! There is no ice in elem
      ! =====	
     dx=bafux(:,elem)
     dy=bafuy(:,elem)

     meancos=sin_elem2D(elem)/cos_elem2D(elem)/earth_radius  !metrics
     vsum=sum(v_ice(elnodes))                           !metrics   
     usum=sum(u_ice(elnodes))                           !metrics

      ! ===== Deformation rate tensor on element elem:
     eps11=sum(dx*u_ice(elnodes))
     eps11=eps11-val3*vsum*meancos                !metrics
     eps22=sum(dy*v_ice(elnodes))
     eps12=0.5_wp*sum(dy*u_ice(elnodes) + dx*v_ice(elnodes))
     eps12=eps12+0.5_wp*val3*usum*meancos          !metrics        
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
    end do
!ICON_OMP_END_PARALLEL_DO
   
end subroutine stress_tensor_omp

!===================================================================
subroutine stress2rhs_omp
! EVP implementation:
! Computes the divergence of stress tensor and puts the result into the
! rhs vectors 

  use mo_ice_elements
  use mo_ice_mesh
  use mo_ice
  use mo_ice_parsup


  use mo_physical_constants,  ONLY: rhoi, rhos, earth_radius
  USE mo_kind,                ONLY: wp

IMPLICIT NONE
INTEGER      :: row, elem, nodels(6), k, i!, j, elnodes(3)
REAL(wp) :: mass, aa
REAL(wp) :: cluster_area,elevation_elem(3)
REAL(wp) :: dx, dy, meancos, val3 ! dx(3), dy(3)

val3=1._wp/3.0_wp
!ICON_OMP_PARALLEL_DO PRIVATE(i,row) ICON_OMP_DEFAULT_SCHEDULE
 DO i=1, myDim_nod2D
     row=myList_nod2D(i) 
     rhs_u(row)=0.0_wp
     rhs_v(row)=0.0_wp
!     rhs_a(row)=0.0_wp    ! these are used as temporal storage here
!     rhs_m(row)=0.0_wp    ! for the contribution due to ssh
 END DO
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(i,row,nodels,k,elem,dx,dy,meancos,&  !j,elnodes,&
!ICON_OMP       aa,elevation_elem) ICON_OMP_DEFAULT_SCHEDULE
DO i=1, myDim_nod2D
    row=myList_nod2D(i)
    nodels=nod2D_elems(:,row)
    ! ===== Skip if ice is absent
     if (m_ice(row)*a_ice(row)==0._wp) CYCLE
    ! =====

     DO k=1,6
        elem=nodels(k)
        IF (elem > 0) THEN
         meancos=sin_elem2D(elem)/cos_elem2D(elem)/earth_radius         !metrics
!         dx=bafux(:,elem)
!         dy=bafuy(:,elem)
         dx=bafux_nod(elem,row)
         dy=bafuy_nod(elem,row)

!         elnodes=elem2D_nodes(:,elem)
!         DO j=1,3
!            if (row==elnodes(j)) exit ! find the node for the give element elem, that matches row
!         END DO

!        rhs_u(row)=rhs_u(row) - voltriangle(elem) * &
!             (sigma11(elem)*dx(j)+sigma12(elem)*(dy(j)) &
!             +sigma12(elem)*val3*meancos)                          !metrics
!        rhs_v(row)=rhs_v(row) - voltriangle(elem) * &
!             (sigma12(elem)*dx(j)+sigma22(elem)*dy(j) &
!             -sigma11(elem)*val3*meancos)
        rhs_u(row)=rhs_u(row) - voltriangle(elem) * &
             (sigma11(elem)*dx+sigma12(elem)*dy &
             +sigma12(elem)*val3*meancos)                          !metrics
        rhs_v(row)=rhs_v(row) - voltriangle(elem) * &
             (sigma12(elem)*dx+sigma22(elem)*dy &
             -sigma11(elem)*val3*meancos)

         ! use rhs_m and rhs_a for storing the contribution from elevation:
         ! do not recalculate at every iteration
!         aa=9.81_wp*voltriangle(elem)/3.0_wp
!         elevation_elem=elevation(elnodes)

!         rhs_a(row)=rhs_a(row)-aa*sum(dx*elevation_elem)
!         rhs_m(row)=rhs_m(row)-aa*sum(dy*elevation_elem)

        END IF
     END DO
END DO
!ICON_OMP_END_PARALLEL_DO

!!!  !ICON_OMP_PARALLEL_DO PRIVATE(i,row,cluster_area,mass) ICON_OMP_DEFAULT_SCHEDULE
!  DO i=1, myDim_nod2D
!     row=myList_nod2D(i)             
!     cluster_area=lmass_matrix(row)
!     mass=cluster_area*(m_ice(row)*rhoi+m_snow(row)*rhos)
     
!     if (mass.ne.0._wp) then
!     rhs_u(row)=rhs_u(row)/mass + rhs_a(row)/cluster_area 
!     rhs_v(row)=rhs_v(row)/mass + rhs_m(row)/cluster_area 
!     else
!     rhs_u(row)=0._wp
!     rhs_v(row)=0._wp
!     end if
!  END DO
!!!  !ICON_OMP_END_PARALLEL_DO

!ICON_OMP_WORKSHARE
    WHERE (rhs_mis > 0._wp)
     rhs_u=rhs_u/rhs_mis + rhs_a
     rhs_v=rhs_v/rhs_mis + rhs_m
    ENDWHERE
!ICON_OMP_END_WORKSHARE
  
end subroutine stress2rhs_omp
!===================================================================
subroutine precalc4rhs_omp
! Some of the quantities used in the solver do not change
! with the subcycles. Hence, can be precalculated and stored.
! Those are rhs_a, rhs_m, mass

  use mo_ice_elements
  use mo_ice_mesh
  use mo_ice
  use mo_ice_parsup

  use mo_physical_constants,  ONLY: rhoi, rhos
  USE mo_kind,                ONLY: wp

IMPLICIT NONE
INTEGER      :: row, elem, elnodes(3), nodels(6), k, i
REAL(wp) :: mass, aa
REAL(wp) :: cluster_area,elevation_elem(3)
REAL(wp) :: dx(3), dy(3)

!ICON_OMP_PARALLEL_DO PRIVATE(i,row) ICON_OMP_DEFAULT_SCHEDULE
 DO i=1, myDim_nod2D
     row=myList_nod2D(i) 
     rhs_a(row)=0.0_wp    ! these are used as temporal storage here
     rhs_m(row)=0.0_wp    ! for the contribution due to ssh
     rhs_mis(row)=0.0_wp
 END DO
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(i,elem,elnodes,aa,dx,dy,elevation_elem,&
!ICON_OMP       k,row) ICON_OMP_DEFAULT_SCHEDULE
 do i=1,myDim_elem2D
     elem=myList_elem2D(i)    
     elnodes=elem2D_nodes(:,elem)
      ! ===== Skip if ice is absent
     aa=product(m_ice(elnodes))*product(a_ice(elnodes))
     if (aa==0._wp) CYCLE
      ! =====
     dx=bafux(:,elem)
     dy=bafuy(:,elem)
     elevation_elem=elevation(elnodes)

     ! use rhs_m and rhs_a for storing the contribution from elevation:
     aa=9.81_wp*voltriangle(elem)/3.0_wp
     DO k=1,3
        row=elnodes(k)
        rhs_a(row)=rhs_a(row)-aa*sum(dx*elevation_elem)
        rhs_m(row)=rhs_m(row)-aa*sum(dy*elevation_elem)
     END DO     
 end do 
!ICON_OMP_END_PARALLEL_DO

!ICON_OMP_PARALLEL_DO PRIVATE(i,row,cluster_area,mass) ICON_OMP_DEFAULT_SCHEDULE
  DO i=1, myDim_nod2D
     row=myList_nod2D(i)             
     cluster_area=lmass_matrix(row)
     mass=cluster_area*(m_ice(row)*rhoi+m_snow(row)*rhos)
     
     if (mass.ne.0._wp) then
     rhs_a(row) = rhs_a(row)/cluster_area 
     rhs_m(row) = rhs_m(row)/cluster_area
     rhs_mis(row) = mass
!!    already taken care of at the initialization step
!     else
!     rhs_u(row)=0._wp
!     rhs_v(row)=0._wp
!     rhs_mis(row)=0._wp
     end if
  END DO
!ICON_OMP_END_PARALLEL_DO
  
end subroutine precalc4rhs_omp
!===================================================================
subroutine EVPdynamics_omp
! EVP implementation. Does cybcycling and boundary conditions.  

  use mo_ice_elements
  use mo_ice_mesh
  use mo_ice
  use mo_ice_parsup

  USE mo_run_config,            ONLY: dtime
  USE mo_sea_ice_nml,           ONLY: evp_rheol_steps, theta_io
  USE mo_physical_constants,    ONLY: rhoi, rhos, Cd_io, rho_ref

  USE mo_kind,    ONLY: wp
  USE mo_ice_fem_utils,   ONLY: exchange_nod2D
  USE mo_kind,            ONLY: wp

IMPLICIT NONE
integer      :: steps, shortstep
REAL(wp) :: rdt
REAL(wp)    ::  drag, inv_mass, det, umod, rhsu, rhsv
integer         ::  i,j
REAL(wp)    :: ax, ay

    rdt=dtime/(1.0_wp*REAL(evp_rheol_steps,wp))
    steps=evp_rheol_steps
    ! theta_io should be zero - set in ice_main.f90
    ax=cos(theta_io)
    ay=sin(theta_io)

! precalculate several arrays that do not change during subcycling
 call precalc4rhs_omp

do shortstep=1, steps 
 ! ===== Boundary conditions
!ICON_OMP_PARALLEL_DO PRIVATE(j,i) ICON_OMP_DEFAULT_SCHEDULE
 do j=1, myDim_nod2D+eDim_nod2D   
    i=myList_nod2D(j) 
    if(index_nod2D(i)==1) then
    u_ice(i)=0.0_wp
    v_ice(i)=0.0_wp
    end if
 end do  
!ICON_OMP_END_PARALLEL_DO
 
 call stress_tensor_omp
 !write(*,*) 'stress', maxval(sigma11), minval(sigma11)
 !if(shortstep<3) write(*,*) mype,'stress ', &
 !minval(sigma11(myList_elem2D(1:myDim_elem2D)))
 
 call stress2rhs_omp
  !write(*,*) 'rhs', maxval(rhs_u), minval(rhs_v)
  !if(shortstep<3)  write(*,*) mype, 'rhs  ', &
  !   maxval(rhs_u(myList_nod2D(1:myDim_nod2D)))

!ICON_OMP_PARALLEL_DO PRIVATE(j,i,inv_mass,umod,drag,rhsu,rhsv,det) ICON_OMP_DEFAULT_SCHEDULE
 do j=1,myDim_nod2D 
    i=myList_nod2D(j)
  if (index_nod2D(i)>0) CYCLE          ! Skip boundary nodes
  if (a_ice(i) > 0.01_wp) then             ! If ice is present, update velocities
   inv_mass=(rhoi*m_ice(i)+rhos*m_snow(i))/a_ice(i)
   inv_mass=max(inv_mass, 9.0_wp)        ! Limit the weighted mass 
                                       ! if it is too small
   inv_mass=1.0_wp/inv_mass

   umod=sqrt((u_ice(i)-u_w(i))**2+(v_ice(i)-v_w(i))**2)
   drag=Cd_io*umod*rho_ref*inv_mass
   
   rhsu=u_ice(i)+rdt*(drag*(ax*u_w(i)-ay*v_w(i))+inv_mass*stress_atmice_x(i)+rhs_u(i))
   rhsv=v_ice(i)+rdt*(drag*(ax*v_w(i)+ay*u_w(i))+inv_mass*stress_atmice_y(i)+rhs_v(i))

   det=(1._wp+ax*drag*rdt)**2+(rdt*coriolis_nod2D(i)+rdt*ay*drag)**2
   det=1.0_wp/det
   u_ice(i)=det*((1.0_wp+ax*drag*rdt)*rhsu+rdt*(coriolis_nod2D(i)+ay*drag)*rhsv)
   v_ice(i)=det*((1.0_wp+ax*drag*rdt)*rhsv-rdt*(coriolis_nod2D(i)+ay*drag)*rhsu)
  ! else                           ! Set ice velocity equal to water velocity
  ! u_ice(i)=u_w(i)
  ! v_ice(i)=v_w(i)
   end if
 end do
 !ICON_OMP_END_PARALLEL_DO

 call exchange_nod2D(u_ice)
 call exchange_nod2D(v_ice)    
 END DO
 


end subroutine EVPdynamics_omp
!===================================================================
