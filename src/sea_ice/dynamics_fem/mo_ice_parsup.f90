!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!----------------------------------------------------------------------------!
module mo_ice_parsup
  !
  ! Arrays used to organize the parallel support
  !
  USE mo_kind,    ONLY: wp
  USE mo_ice_elements
  USE mo_ice_mesh

  IMPLICIT NONE

  PUBLIC :: set_par_support

  PRIVATE :: communication_nod
  PRIVATE :: mymesh

  PUBLIC
  SAVE

! Einar: No parallelization yet
!#ifdef PETSC
!#include "finclude/petsc.h"
!#else
!  include 'mpif.h'
!#endif
  integer      :: maxPEnum
  type com_struct
     integer    :: rPEnum
     integer, dimension(:), pointer :: rPE
     integer, dimension(:), pointer :: rptr
     integer, dimension(:), pointer :: rlist
     integer    :: sPEnum
     integer, dimension(:), pointer :: sPE
     integer, dimension(:), pointer :: sptr
     integer, dimension(:), pointer :: slist
  end type com_struct

  type(com_struct)   :: com_nod2D, com_nod3D
  ! Buffer arrays to store information to be communicated
  type com_array
     REAL(wp), dimension(:), pointer :: array
  end  type com_array
  type(com_array), allocatable             :: s_buff_ssh(:), r_buff_ssh(:)
  type(com_array), allocatable             :: s_buff_ts(:), r_buff_ts(:)

  ! general MPI part
  integer            :: MPIERR
  integer            :: npes=0
  integer            :: mype=0
  integer, allocatable, dimension(:)  :: part, part3D

  ! Mesh partition
  integer                             :: myDim_nod2D, eDim_nod2D
  integer, allocatable, dimension(:)  :: myList_nod2D
  integer                             :: myDim_nod3D, eDim_nod3D
  integer, allocatable, dimension(:)  :: myList_nod3D
  integer                             :: myDim_elem2D
  integer, allocatable, dimension(:)  :: myList_elem2D
  integer                             :: myDim_elem3Dr  !(regular)
  integer                             :: myDim_elem3D   !(total)
  integer, allocatable, dimension(:)  :: myList_elem3D


CONTAINS

!subroutine par_init    ! initializes MPI
!
!  use PARSUP
!  implicit none
!
!
!  integer :: i
!
!  call MPI_INIT(i);
!  call MPI_Comm_Size(MPI_COMM_WORLD,npes,i)
!  call MPI_Comm_Rank(MPI_COMM_WORLD,mype,i)
!  maxPEnum=npes
!  call PETSCInitialize(PETSC_NULL_CHARACTER,i);
!
!  write(*,*) 'MPI has been initialized'
!end subroutine par_init
!=================================================================
!subroutine par_ex       ! finalizes MPI
!  use PARSUP
!  implicit none
!
!  call  MPI_Barrier(MPI_COMM_WORLD,MPIerr)
!  call  MPI_Finalize(MPIerr)
!
!end subroutine par_ex
!=================================================================
subroutine set_par_support

  integer   n, j, k, nini, nend, count

  ! Define partitioning vector

  allocate(part(nod2D))
  ! Geometrical partitioning (do not use except for ordered mesh)
  do n=0, npes-1
     nini=(nod2D/npes)*n+1
     nend=(nod2D/npes)*(n+1)
     if (n==npes-1) nend=nod2D
     part(nini:nend)=n
  end do
  ! This is METIS-based partitioning
   !call partit(icestiff%dim,icestiff%rowptr,icestiff%colind,npes, part)


  call communication_nod
  call mymesh

  ! Allocate communication buffers:


  if (npes>1) then
     allocate(s_buff_ssh(com_nod2D%sPEnum),r_buff_ssh(com_nod2D%sPEnum))
     do n=1, com_nod2D%sPEnum
        count=com_nod2D%sptr(n+1) - com_nod2D%sptr(n)
        allocate(s_buff_ssh(n)%array(count))
     end do
     do n=1, com_nod2D%rPEnum
        count=com_nod2D%rptr(n+1) - com_nod2D%rptr(n)
        allocate(r_buff_ssh(n)%array(count))
     end do

  end if
!   write(*,*) 'Communication arrays are set'
end subroutine set_par_support
!=======================================================================
subroutine communication_nod

  integer n,np, nz, prank, elem, elnodes(3), epe(3), counter, nini, nend
  integer, allocatable :: aux(:,:), pnum(:,:)

  ! Assume we have 2D partitioning vector in part. Find communication
  ! rules
  allocate(aux(npes,nod2D))
  allocate(pnum(npes,npes))
  aux=0
  pnum=0
  do elem=1,elem2D
     elnodes=elem2D_nodes(:,elem)
     epe=part(elnodes)+1
     if(epe(1).ne.epe(2)) then
        aux(epe(1), elnodes(2))=1
        aux(epe(2), elnodes(1))=1
     end if

     if(epe(2).ne.epe(3)) then
        aux(epe(3), elnodes(2))=1
        aux(epe(2), elnodes(3))=1
     end if
     if(epe(1).ne.epe(3)) then
        aux(epe(1), elnodes(3))=1
        aux(epe(3), elnodes(1))=1
     end if
  end do

  do n=1, nod2D
     do np=1, npes
        if(aux(np,n).ne.0) then
           pnum(np,part(n)+1)=pnum(np,part(n)+1)+1
        end if
     end do
  end do


  ! We know how many external nodes each PE needs
  ! This is the 'receive' list
  ! com_nod2D for 2D nodes
  !

  ! The number of external PE I receive information from
  com_nod2D%rPEnum=0
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        com_nod2D%rPEnum=com_nod2D%rPEnum+1
     end if
  end do

  ! Their ranks (PE numbers)

  counter=0
  allocate(com_nod2D%rPE(com_nod2D%rPEnum))
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        com_nod2D%rPE(counter)=n-1
     end if
  end do


  ! Ptr to list of external nodes ordered by external PE ranks

  counter=0
  allocate(com_nod2D%rptr(com_nod2D%rPEnum+1))
  com_nod2D%rptr(1)=1
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        com_nod2D%rptr(counter+1)=com_nod2D%rptr(counter)+ pnum(mype+1,n)
     end if
  end do

  ! List itself

  counter=0
  allocate(com_nod2D%rlist(com_nod2D%rptr(com_nod2D%rPEnum+1)-1))
  do np=1,com_nod2D%rPEnum
     prank=com_nod2D%rPE(np)
     do n=1, nod2D
        if((aux(mype+1,n)==1).and.(part(n)==prank)) then
           counter=counter+1
           com_nod2D%rlist(counter)=n
        end if
     end do
  end do
  ! Summary of this piece: mype receives
  ! information on external 2D nodes from
  ! comm_nod2D%rPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_nod2D%rPE(:)
  ! Pointers to external node numbers are in
  ! comm_nod2D%rptr(:)
  ! The node numbers are in
  ! comm_nod2D%list(:)
  ! Putting everything into structure takes many operations, but
  ! without the structure we will need to many names and arrays
  ! Do not forget that we need also send part, and we need analogous
  ! information for 3D nodes.

  ! SENDING PART
  com_nod2D%sPEnum=0
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        com_nod2D%sPEnum=com_nod2D%sPEnum+1
     end if
  end do

  ! Their ranks (PE numbers)

  counter=0
  allocate(com_nod2D%sPE(com_nod2D%sPEnum))
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        counter=counter+1
        com_nod2D%sPE(counter)=n-1
     end if
  end do

  ! Ptr to list of external nodes ordered by external PE ranks
  counter=0
  allocate(com_nod2D%sptr(com_nod2D%sPEnum+1))
  com_nod2D%sptr(1)=1
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        counter=counter+1
        com_nod2D%sptr(counter+1)=com_nod2D%sptr(counter)+ pnum(n,mype+1)
     end if
  end do

  ! List itself

  counter=0
  allocate(com_nod2D%slist(com_nod2D%sptr(com_nod2D%sPEnum+1)-1))
  do np=1,com_nod2D%sPEnum
     prank=com_nod2D%sPE(np)
     do n=1, nod2D
        if((aux(prank+1,n)==1).and.(part(n)==mype)) then
           counter=counter+1
           com_nod2D%slist(counter)=n
        end if
     end do
  end do

  ! mype sends its data to
  ! comm_nod2D%sPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_nod2D%sPE(:)
  ! Pointers to external node numbers are in
  ! comm_nod2D%sptr(:)
  ! The node numbers are in
  ! comm_nod2D%list(:)

  deallocate(pnum,aux)

end subroutine communication_nod
!==========================================================================


subroutine mymesh

  integer     n, counter, q

  !======= NODES

  ! Owned nodes + external nodes which I need:
  counter=0
  do n=1, nod2D
     if (part(n)==mype) counter=counter+1
  end do
  myDim_nod2D=counter
  eDim_nod2D=com_nod2D%rptr(com_nod2D%rPEnum+1)-1
  allocate(myList_nod2D(myDim_nod2D+eDim_nod2D))
  counter=0
  do n=1, nod2D
     if (part(n)==mype) then
        counter=counter+1
        myList_nod2D(counter)=n
     end if
  end do
  myList_nod2D(myDim_nod2D+1:myDim_nod2D+eDim_nod2D)=&
       com_nod2D%rlist

  ! Summary:
  ! myList_nod2D(myDim_nod2D+1:myDim_nod2D+eDim_nod2D)
  ! contains external nodes which mype needs;
  ! myList_nod2D(1:myDim_nod2D) contains owned nodes

  !======= ELEMENTS
  ! 2D elements
  counter=0
  do n=1, elem2D
     do q=1,3
        if(part(elem2D_nodes(q,n))==mype) then
           counter=counter+1
           exit
        end if
     end do
  end do
  myDim_elem2D=counter
  allocate(myList_elem2D(myDim_elem2D))
  counter=0
  do n=1, elem2D
     do q=1,3
        if(part(elem2D_nodes(q,n))==mype) then
           counter=counter+1
           myList_elem2D(counter)=n
           exit
        end if
     end do
  end do


!  write(0,*) "max,min of part", maxval(part), minval(part)
!  write(0,*) mype, ":: elem2D=", elem2D, " myDim_elem2D=", myDim_elem2D
!  write(0,*) mype, ":: nod2D=", nod2D, " myDim_nod2D=", myDim_nod2D
!  write(0,*) "max(myList_nod2D)", maxval(myList_nod2D)
!  write(0,*) "max(myList_elem2D)", maxval(myList_elem2D)

end subroutine mymesh

!=======================================================================

!subroutine exchange_nod2D(nod_array2D)
!USE MESH
!USE ICE
!USE PARSUP
!IMPLICIT NONE
!
!
! INTEGER  :: sreq(maxPEnum)
! INTEGER  :: rreq(maxPEnum)
! INTEGER  :: sstat(MPI_STATUS_SIZE,maxPEnum)
! INTEGER  :: rstat(MPI_STATUS_SIZE,maxPEnum)
! integer    :: n, sn, rn, dest, nini, nend, count, source,tag
! real(kind=8) :: nod_array2D(nod2D)
!  sn=com_nod2D%sPEnum
!  rn=com_nod2D%rPEnum
!  ! Put data to be communicated into send buffer
!
!  Do n=1,com_nod2D%sPEnum
!   nini=com_nod2D%sptr(n)
!   nend=com_nod2D%sptr(n+1)-1
!   s_buff_ssh(n)%array=nod_array2D(com_nod2D%slist(nini:nend))
!  end do
!
!  DO n=1, sn
!     dest=com_nod2D%sPE(n)
!     nini=com_nod2D%sptr(n)
!     count=com_nod2D%sptr(n+1) - nini
!
!
!     call MPI_ISEND(s_buff_ssh(n)%array, count, MPI_DOUBLE_PRECISION, dest, mype, &
!               MPI_COMM_WORLD, sreq(n), MPIerr)
!  END DO
!  DO n=1,rn
!     source=com_nod2D%rPE(n)
!     nini=com_nod2D%rptr(n)
!     count=com_nod2D%rptr(n+1) - nini
!
!
!     call MPI_IRECV(r_buff_ssh(n)%array, count, MPI_DOUBLE_PRECISION, source, &
!               source, MPI_COMM_WORLD, rreq(n), MPIerr)
!
!  END DO
!
!     call MPI_WAITALL(sn,sreq,sstat, MPIerr)
!     call MPI_WAITALL(rn,rreq,rstat, MPIerr)
!
!  ! Put received data to their destination
!  Do n=1,com_nod2D%rPEnum
!   nini=com_nod2D%rptr(n)
!   nend=com_nod2D%rptr(n+1)-1
!   nod_array2D(com_nod2D%rlist(nini:nend))=r_buff_ssh(n)%array
!  end do
!
!END SUBROUTINE exchange_nod2D
!============================================================================
!subroutine broadcast_nod2D(arr2D)
!! A 2D version of the previous routine
!use PARSUP
!USE MESH
!USE ELEMENTS
!
!IMPLICIT NONE
!
!INTEGER :: ireals
!integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
!INTEGER, ALLOCATABLE, DIMENSION(:) ::  isendbuf, irecvbuf
!
!real(kind=8) ::  arr2D(nod2D)
!real(kind=8), ALLOCATABLE, DIMENSION(:) ::  sendbuf, recvbuf
!
!IF ( mype == 0 ) THEN
!
!    DO  n = 1, npes-1
!
!       CALL MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
!                     0, MPI_COMM_WORLD, status, MPIerr )
!       sender = status(MPI_SOURCE)
!       ALLOCATE( recvbuf(1:nTS), irecvbuf(1:nTS) )
!       CALL MPI_RECV( irecvbuf(1), nTS, MPI_INTEGER, sender, &
!                      1, MPI_COMM_WORLD, status, MPIerr )
!       CALL MPI_RECV( recvbuf(1), nTS, MPI_DOUBLE_PRECISION, sender, &
!                      2, MPI_COMM_WORLD, status, MPIerr )
!
!       DO i = 1, nTS
!          arr2D(irecvbuf(i)) = recvbuf(i)
!       ENDDO
!       DEALLOCATE( recvbuf, irecvbuf )
!
!    ENDDO
!
!ELSE
!
!    ALLOCATE( sendbuf(1:myDim_nod2d), isendbuf(1:myDim_nod2D) )
!    DO n = 1, myDim_nod2D
!       isendbuf(n) = myList_nod2D(n)
!       sendbuf(n)  = arr2D(myList_nod2D(n))
!    ENDDO
!    CALL MPI_SEND( myDim_nod2D, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPIerr )
!    CALL MPI_SEND( isendbuf(1), myDim_nod2D, MPI_INTEGER, 0, 1, &
!                   MPI_COMM_WORLD, MPIerr )
!    CALL MPI_SEND( sendbuf(1), myDim_nod2D, MPI_DOUBLE_PRECISION, &
!                   0, 2, MPI_COMM_WORLD, MPIerr )
!    DEALLOCATE( sendbuf, isendbuf )
!
!ENDIF
!
!CALL MPI_BCAST( arr2D, nod2d, MPI_DOUBLE_PRECISION, 0, &
!                 MPI_COMM_WORLD, MPIerr)
!
!end subroutine broadcast_nod2D
!===================================================================

end module mo_ice_parsup
