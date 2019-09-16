program test_index_list

  use mo_index_list
  use mo_kind

  implicit none

  integer,     parameter :: n  = 200000
  integer,     parameter :: nb = 7
  integer(i4), target    :: conditions(n,nb)
  integer,     target    :: indices(n,nb), dev_indices(n,nb)
  integer                :: nvalid(nb), dev_nvalid(nb)
  real                   :: harvest(n,nb)

  integer :: i, b

  call random_number(harvest)
  conditions = int(harvest * 2)

  !$acc data copyin(conditions) create(dev_indices, dev_nvalid)

  ! Test the non-batched version

  nvalid(1) = 0
  do i = 1, n
    if (conditions(i,1)) then
      nvalid(1) = nvalid(1) + 1
      indices(nvalid(1),1) = i
    end if
  end do


  call generate_index_list(conditions(:,1), dev_indices(:,1), 1, n, dev_nvalid(1), 1)
  !$acc wait(1)
  !$acc update host(dev_indices(:,1))

  print *, "CHECK NON-BATCHED: ", nvalid(1) == dev_nvalid(1), all(indices(:,1) == dev_indices(:,1))

  ! Test the batched version

  nvalid = 0
  do b = 1, nb
    do i = 1, n
      if (conditions(i,b)) then
        nvalid(b) = nvalid(b) + 1
        indices(nvalid(b),b) = i
      end if
    end do
  end do

  call generate_index_list_batched(conditions, dev_indices, 1, n, dev_nvalid, 1)
  !$acc wait(1)
  !$acc update host(dev_indices, dev_nvalid)

  print *, "CHECK BATCHED: ", all(nvalid == dev_nvalid), all(indices == dev_indices)

  !$acc end data

!  print *, "n=", nvalid, " data: ", indices
!
!  print *, "DEVICE:"
!  print *, "n=", dev_nvalid, " data: ", dev_indices


end program test_index_list
