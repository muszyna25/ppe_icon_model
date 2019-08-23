module mo_simple_dump

    use mo_kind

    implicit none

    private

    public :: dump2text

    interface dump2text
        module procedure dump2text_vp
        module procedure dump2text_wp
        module procedure dump2text_sp
    end interface dump2text

    contains

        subroutine dump2text_vp(array, name, preserve_original)
            real(vp), intent(inout) :: array(:, :, :)
            character(len=*), intent(in) :: name
            logical, intent(in), optional :: preserve_original

            logical :: do_copy
            real(vp), allocatable :: backup(:,:,:)

            do_copy = .true.
            if (present(preserve_original)) do_copy = preserve_original

            if (do_copy) backup = array

            !$acc update host(array)
            open(7, file=name)
            write (7, *) array
            close(7)

            if (do_copy) array = backup
        end subroutine dump2text_vp


        subroutine dump2text_wp(array, name, preserve_original)
            real(wp), intent(inout) :: array(:, :, :)
            character(len=*), intent(in) :: name
            logical, intent(in), optional :: preserve_original

            logical :: do_copy
            real(wp), allocatable :: backup(:,:,:)

            do_copy = .true.
            if (present(preserve_original)) do_copy = preserve_original

            if (do_copy) backup = array

            !$acc update host(array)
            open(7, file=name)
            write (7, *) array
            close(7)

            if (do_copy) array = backup
        end subroutine dump2text_wp


        subroutine dump2text_sp(array, name, preserve_original)
            real(sp), intent(inout) :: array(:, :, :)
            character(len=*), intent(in) :: name
            logical, intent(in), optional :: preserve_original

            logical :: do_copy
            real(sp), allocatable :: backup(:,:,:)

            do_copy = .true.
            if (present(preserve_original)) do_copy = preserve_original

            if (do_copy) backup = array

            !$acc update host(array)
            open(7, file=name)
            write (7, *) array
            close(7)

            if (do_copy) array = backup
        end subroutine dump2text_sp

end module mo_simple_dump
