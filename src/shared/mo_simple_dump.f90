module mo_simple_dump

    use mo_kind

    implicit none

    private

    public :: dump2text

    interface dump2text
        module procedure dump2text1_wp
        module procedure dump2text2_wp
        module procedure dump2text3_wp
    end interface dump2text

    contains

        subroutine dump2text1_wp(array, name, preserve_original)
            real(wp), intent(inout) :: array(:)
            character(len=*), intent(in) :: name
            logical, intent(in), optional :: preserve_original

            logical :: do_copy
            real(wp), allocatable :: backup(:)

            do_copy = .true.
            if (present(preserve_original)) do_copy = preserve_original

            if (do_copy) backup = array

            !$acc update host(array)
            open(17, file=name)
            write (17, *) array
            close(17)

            if (do_copy) array = backup
        end subroutine dump2text1_wp

        subroutine dump2text2_wp(array, name, preserve_original)
            real(wp), intent(inout) :: array(:, :)
            character(len=*), intent(in) :: name
            logical, intent(in), optional :: preserve_original

            logical :: do_copy
            real(wp), allocatable :: backup(:,:)

            do_copy = .true.
            if (present(preserve_original)) do_copy = preserve_original

            if (do_copy) backup = array

            !$acc update host(array)
            open(17, file=name)
            write (17, *) array
            close(17)

            if (do_copy) array = backup
        end subroutine dump2text2_wp

        subroutine dump2text3_wp(array, name, preserve_original)
            real(wp), intent(inout) :: array(:, :, :)
            character(len=*), intent(in) :: name
            logical, intent(in), optional :: preserve_original

            logical :: do_copy
            real(wp), allocatable :: backup(:,:,:)

            do_copy = .true.
            if (present(preserve_original)) do_copy = preserve_original

            if (do_copy) backup = array

            !$acc update host(array)
            open(17, file=name)
            write (17, *) array
            close(17)

            if (do_copy) array = backup
        end subroutine dump2text3_wp



end module mo_simple_dump
