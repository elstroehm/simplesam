module module_output

!!!Dieses Modul enthält einige IO-Funktionen

    use module_transform

    implicit none

contains

    subroutine write_file(unitnumber, f_values, io_error)
        integer, intent(in)              :: unitnumber
        real(p), dimension(:,:), intent(in) :: f_values
        integer, intent(in)              :: io_error

        if (io_error == 0) then
            write(unitnumber,rec=1) f_values
        else
            print*, "error occured during write, iostat value: ", io_error
        end if

    end subroutine write_file

    subroutine write_array(unitnumber, array, io_error)
        integer, intent(in)               :: unitnumber
        real(p), dimension(:), intent(in)    :: array
        integer, intent(in)               :: io_error

        if(io_error == 0) then
            write(unitnumber, rec=1) array
        else
            print*, "error occured during write, iostat value: ", io_error
        end if

    end subroutine write_array

    subroutine print_params(params)
        type(mapping), intent(in)      :: params
        integer                        :: i

        print*, "Transformation matrix: "
        do i=1, params%dim
            print*, params%a(i,:)
        end do
        print*, "Inverse of transformation matrix: "
        do i=1, params%dim
            print*, params%a_inv(i,:)
        end do
        print*, "Translation: "
        print*, params%b
        print*, "________________________"

    end subroutine print_params

    subroutine print_coords(coords_array)
        real(p), dimension(:,:), intent(in) :: coords_array
        integer                             :: i

        do i = 1, size(coords_array(:,1))
            print*, coords_array(i,:)
        end do
        print*, "______________________________________"

    end subroutine print_coords

end module module_output
