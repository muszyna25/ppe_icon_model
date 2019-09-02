program test_index_list

  use mo_index_list
  use mo_kind

  implicit none

  integer,     parameter :: n = 200000
  integer(i1), target    :: conditions(n)
  integer,     target    :: indices(n), dev_indices(n)
  integer                :: nvalid, dev_nvalid
  real                   :: harvest(n)

  integer :: i

  call random_number(harvest)
  conditions = int(harvest * 2)

  nvalid = 0
  do i = 1, n
    if (conditions(i)) then
      nvalid = nvalid + 1
      indices(nvalid) = i
    end if
  end do

!  print *, "Predicates: ", conditions


  !$acc data copyin(conditions) copyout(dev_indices)
  call generate_index_list(conditions, dev_indices, dev_nvalid, 1)
  !$acc wait(1)
  !$acc end data

!  print *, "n=", nvalid, " data: ", indices
!
!  print *, "DEVICE:"
!  print *, "n=", dev_nvalid, " data: ", dev_indices

  print *, "CHECK:"
  print *, all(indices == dev_indices)


end program test_index_list
