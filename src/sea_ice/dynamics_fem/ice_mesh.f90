! A set of auxilliary routines required to run standing-alone ice.
! They are not needed if ice is used as a part of coupled system
! Their role: reading the mesh and preparing it to further usage 


subroutine set_mesh
implicit none
! Einar: We don't read the mesh in from file. coord_nod2D and elem2D_nodes are set directly in
! ice_init in mo_sea_ice.f90
  !call read_2Dmesh
  call standard_element_definition
  call mesh_scaling 
  call basisfunctions
  call build_nghbr_arrays
end subroutine set_mesh  

!----------------------------------------------------------------------------
!
subroutine standard_element_definition
  !
  !  DEFINITION OF :
  !    - 2D STANDARD ELEMENTS
  !    - BASISFUNCTIONS ON 2D STANDARD ELEMENTS (stdbafu)
  !     Zeta    : stdbafu(1)=1-x-y   stdbafu(2)=x       stdbafu(3)=y  
  !    - SCALARPRODUCTS
  !    - DERIVATIVES OF STD.BASISFUNCTIONS
  !
  use mo_ice_param
  use mo_ice_elements
  use mo_ice_mesh
  USE mo_kind,  ONLY: wp
  implicit none

  integer :: i,j
  Vol2D=1._wp/2._wp
  allocate(derivative_stdbf(2,3))
  !
  derivative_stdbf= 0._wp
  derivative_stdbf(:,1)= -1._wp
  derivative_stdbf(1,2)=  1._wp
  derivative_stdbf(2,3)=  1._wp
  !
end subroutine standard_element_definition
!
!----------------------------------------------------------------------------
!
subroutine basisfunctions
  use mo_ice_param
  use mo_ice_elements
  use mo_ice_mesh
  USE mo_kind,  ONLY: wp
  implicit none

  REAL(wp)                           :: DET2D
  REAL(wp), dimension(2, 3)          :: derivative_loczeta
  REAL(wp), dimension(2,2)           :: jacobian2D, jacobian2D_inv
  integer                                :: elem,i

  allocate(bafux(3,elem2d), bafuy(3,elem2d))
  allocate(voltriangle(elem2d))

  bafux    =0._wp
  bafuy    =0._wp
  voltriangle   =0._wp

  do elem=1,elem2d
     call local_element_def(elem,   & 
          jacobian2D, jacobian2D_inv, DET2D,  derivative_loczeta)
     do i=1,3
        bafux(i,elem) = derivative_loczeta(1,i)
        bafuy(i,elem) = derivative_loczeta(2,i)
     enddo
     voltriangle(elem) = abs(DET2D) * Vol2D
  enddo

end subroutine basisfunctions
!
!----------------------------------------------------------------------------
!
subroutine local_element_def(element, &
     jacobian2D, jacobian2D_inv, DET, derivative_loczeta)
  !--------------------------------------------------------------------------
  !                     2D LOCAL ELEMENT DEFINITIONS
  !     USES: 
  !       coord_nod2D(2,nod2D)
  !       elem2D_nodes(3, elem2D)
  !       derivative_stabafu_x_2D(2, 3)
  !
  !     OUTPUT:
  !       jacobian2D(2,2), jacobian2D_inv(2,2), DET
  !       derivative_locbafu_x_2D(2, 3)
  !--------------------------------------------------------------------------
  use mo_ice_param
  use mo_ice_elements
  use mo_ice_mesh
  USE mo_kind,  ONLY: wp
  implicit none
  !
  integer, intent(IN)                         :: element
  REAL(wp), dimension(2,2),  intent(OUT)  :: jacobian2D
  REAL(wp), dimension(2,2),  intent(OUT)  :: jacobian2D_inv
  REAL(wp), intent(OUT)                   :: DET
  REAL(wp), dimension(2, 3), intent(OUT)  :: derivative_loczeta  
  !
  REAL(wp), dimension(2, 3)               :: local_cart
  REAL(wp), dimension(3, 2)               :: der_transp
  integer                                     :: i
  integer                                     :: node
  REAL(wp)                                :: meancos
  !
  meancos=cos_elem2D(element)
  do i=1, 3
     node=elem2D_nodes(i,element)
     ! cartesian coordinates
     local_cart(1,i)=coord_nod2D(1,node)
     local_cart(2,i)=r_earth * coord_nod2D(2,node)
  end do
  !
  !  T R A N S F O R M A T I O N - M A T R I X   "jacobian"
  do i=1,2
     jacobian2D(:,i)= local_cart(:,i+1)-local_cart(:,1)
     if (cyclic_geometry) then
        if (jacobian2D(1,i)> domain_length/2.0_wp) jacobian2D(1,i)=jacobian2D(1,i)-domain_length
        if (jacobian2D(1,i)<-domain_length/2.0_wp) jacobian2D(1,i)=jacobian2D(1,i)+domain_length
     end if
  end do
  jacobian2D(1,:)=jacobian2D(1,:)*meancos*r_earth
  !
  !  I N V E R S E   O F   jacobian
  call matrix_inverse_2x2(jacobian2D, jacobian2D_inv, DET)
  !
  der_transp=matmul(transpose(derivative_stdbf), jacobian2D_inv)
  derivative_loczeta=transpose(der_transp)
  !
end subroutine local_element_def
!
!----------------------------------------------------------------------------
!
subroutine  matrix_inverse_2x2 (A, AINV, DET)
  USE mo_kind,  ONLY: wp
  !
  implicit none
  !
  REAL(wp), dimension(2,2), intent(IN)  :: A
  REAL(wp), dimension(2,2), intent(OUT) :: AINV
  REAL(wp), intent(OUT)                 :: DET
  !
  integer                                   :: i,j
  !
  DET  = A(1,1)*A(2,2) - A(1,2)*A(2,1)
  if ( DET .eq. 0.0_wp )  then
     do j=1,2
        write(*,*) (A(i,j),i=1,2)
     end do
     stop 'SINGULAR 2X2 MATRIX'
  else
     AINV(1,1) =  A(2,2)/DET
     AINV(1,2) = -A(1,2)/DET
     AINV(2,1) = -A(2,1)/DET
     AINV(2,2) =  A(1,1)/DET
  endif
end subroutine matrix_inverse_2x2
!
!----------------------------------------------------------------------------
!
!subroutine read_2Dmesh
!  !
!  !--------------------------------------------------------------------------
!  ! READS 2D MESH AND VARIABLES CONCERNING THE 2D MESH    
!  !
!  ! OUTPUT:
!  !   Number of 2D-Nodes: nod2D
!  !   Number of 2D-Elements: elem2D
!  !   Number of 2D-element-neighbours is always 3
!  !
!  !   Mesh-Arrays: 
!  !     coord_nod2D(2,nod2D)
!  !     elem2D_nodes(3,elem2D)
!  !     index_nod2D(nod2D)
!  !     elem2D_nghbrs(3,elem2D)
!  !--------------------------------------------------------------------------
!  !
!  use PARAM
!  use MESH
!  use ELEMENTS
!  implicit none
!  !
!  integer           :: i, n, ind
!  REAL(wp)      :: x1, x2
!
!  open (20,file=trim(meshpath)//'nod2d.out',  status='old')
!  open (21,file=trim(meshpath)//'elem2d.out', status='old')
!  write(*,*) '2D mesh is opened'
!  ! 
!  ! nodes
!  !
!  read(20,*) nod2D       
!  allocate(coord_nod2D(2,nod2D), index_nod2D(nod2D))
!  do i=1, nod2D
!     read(20,*) n, x1, x2, ind
!     coord_nod2D(1,n)=x1 ! x-coords in degrees
!     coord_nod2D(2,n)=x2 ! y-coords in degrees
!     index_nod2D(n) =ind ! 1 or 0 (land mask)
!  end do
!  close(20)
!  !
!  ! elements
!  !
!  read(21,*)  elem2D      
!  allocate(elem2D_nodes(3,elem2D))  
!  do n=1, elem2D
!     read(21,*) elem2D_nodes(:,n) ! Nodes that make up each element
!  end do
!  close(20)
!end subroutine read_2Dmesh
!
!---------------------------------------------------------------------------
!
subroutine mesh_scaling  
  !
  ! Transforms degrees in rad in coord_nod2D(2,nod2D)     
  ! Fills in the arrays cos_elem2D(elem2D) and sin_elem2D(elem2D)
  !
  use mo_ice_param
  use mo_ice_elements
  use mo_ice_mesh
  USE mo_kind,  ONLY: wp
  !
  implicit none
  !
  integer         :: i, j, n, node
  !
  !  degrees to radians
  !
  coord_nod2D(1,:)=coord_nod2D(1,:)*rad
  coord_nod2D(2,:)=coord_nod2D(2,:)*rad

  allocate(cos_elem2D(elem2D), sin_elem2D(elem2D))

  do i=1, elem2D  
     cos_elem2D(i)=sum(cos(coord_nod2D(2,elem2D_nodes(:,i))))/3.0_wp
     sin_elem2D(i)=sum(sin(coord_nod2D(2,elem2D_nodes(:,i))))/3.0_wp
  end do
  !
  if(cartesian) cos_elem2D=1.0_wp
  if(cartesian) sin_elem2D=0.0_wp
end subroutine mesh_scaling
!
!---------------------------------------------------------------------------
!
subroutine build_nghbr_arrays
  !
  use mo_ice_param
  use mo_ice_elements
  use mo_ice_mesh
 
  !
  implicit none
  !
  integer                            :: j, k, m, a, tr(3), pr(6), flag, el,ml
  integer                            :: i, eltype, dgnod(2), pos
  integer, dimension(100)            :: AUX=0
  integer, allocatable, dimension(:) :: ind, check
  !
  !--------------- 2D mesh-------------------------------
  ! Builds nod_in_elem2D  
  allocate(ind(nod2D), check(elem2D))
  ind=0
  do j=1,elem2D
     tr=elem2D_nodes(:,j)
     ind(tr)=ind(tr)+1
  end do
  allocate(nod_in_elem2D(nod2D))
  
  nod_in_elem2D%nmb=ind(1:nod2D)    
  do j=1,nod2D   
     allocate(nod_in_elem2D(j)%addresses(ind(j)))
     
  end do
  ind=0
  do i=1,elem2D   
     tr=elem2D_nodes(:,i)
     ind(tr)=ind(tr)+1
     do j=1, 3
        nod_in_elem2D(tr(j))%addresses(ind(tr(j)))=i
     end do
  end do
  ! The list of elements is ordered, and no sorting is needed
  !
  ! Builds nghbr_nod2D  
  allocate(nghbr_nod2D(nod2D))
  
  check=0
  do j=1, nod2D
     flag=0
     do m=1,nod_in_elem2D(j)%nmb
        el=nod_in_elem2D(j)%addresses(m)
        do k=1, 3
           a=elem2D_nodes(k,el)       
           if (check(a)==0) then  
              check(a)=1 
              flag=flag+1         
              aux(flag)=a
           end if
        end do
     end do
     nghbr_nod2D(j)%nmb=flag
     allocate(nghbr_nod2D(j)%addresses(flag))
     
     ! We need to sort array aux(1:flag)
     do m=flag,1,-1
        ml=maxloc(aux(1:flag), 1)
        nghbr_nod2D(j)%addresses(m)=aux(ml)
        check(aux(ml))=0
        aux(ml)=-999
     end do
  end do
  !
  deallocate(ind, check)
end subroutine build_nghbr_arrays
