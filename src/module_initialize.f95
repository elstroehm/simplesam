module module_initialize

!!!Dieses Modul enthält einige Methoden zur Initialisierung der Transformation und Koordianten.

    use module_matrix
    use module_transform
    use module_function

    implicit none

contains

    !!!Zufällige Startparameter für die Transformation
    function init_params(dim)
        type(mapping)       :: init_params
        integer, intent(in) :: dim
        integer             :: i, j, seed
        real(p)                :: time, number

        call init_mapping_empty(init_params, dim)
        do i = 1, dim
            do j = 1, dim
                call cpu_time(time)
                call random_seed(seed)
                call random_number(number)
                init_params%a(i, j) = number * 10.0 - 5.0
            end do
        end do
        do i = 1, dim
            call cpu_time(time)
            call random_seed(seed)
            call random_number(number)
            init_params%b(i) = number * 3.0 - 1.5
        end do

        init_params%dim = dim
        init_params%a_inv = init_params%a
        call invert_matrix(init_params%a_inv, dim)

    end function init_params

    subroutine random_params(params)
        type(mapping), intent(inout) :: params
        integer                      :: i, j, seed
        real(p)                      :: time, number

        do i = 1, params%dim
            do j = 1, params%dim
                call cpu_time(time)
                call random_seed(seed)
                call random_number(number)
                params%a(i, j) = number * 10.0 - 5.0
            end do
        end do
        do i = 1, params%dim
            call cpu_time(time)
            call random_seed(seed)
            call random_number(number)
            params%b(i) = number * 3.0 - 1.5
        end do

        params%a_inv = params%a
        call invert_matrix(params%a_inv, params%dim)

    end subroutine random_params

    !!!Zufällige, aber homogen verteilte, Punktkoordinaten im W-Raum
    function init_coords(amount, dim)
        real(p), allocatable, dimension(:,:) :: init_coords
        integer, intent(in)               :: amount, dim
        integer                           :: i, j, seed
        real(p)                              :: time, number

        allocate(init_coords(1:amount,1:dim))
        do i = 1, amount
            do j = 1, dim
                call cpu_time(time)
                call random_seed(seed)
                call random_number(number)
                init_coords(i, j) = number * 2.0 - 1.0
            end do
        end do

    end function init_coords

    subroutine init_coords_random(w_coords_array)
        real(p), dimension(:,:), intent(inout) :: w_coords_array
        integer                                :: i, j, seed
        real(p)                                :: time, number

        do i = 1, size(w_coords_array(:,1))
            do j = 1, size(w_coords_array(1,:))
                call cpu_time(time)
                call random_seed(seed)
                call random_number(number)
                w_coords_array(i, j) = number * 2.0 - 1.0
            end do
        end do

    end subroutine init_coords_random

    subroutine init_array_random(array)
        real(p), dimension(:), intent(inout) :: array
        integer                              :: i, seed
        real(p)                              :: time, number

        do i = 1, size(array)
            call cpu_time(time)
            call random_seed(seed)
            call random_number(number)
            array(i) = number * 2.0 - 1.0
        end do

    end subroutine init_array_random

    subroutine init_number_random(number)
        real(p), intent(inout) :: number
        integer                :: seed
        real(p)                :: time

        call cpu_time(time)
        call random_seed(seed)
        call random_number(number)

    end subroutine init_number_random

    !!!Liest alle für das Programm nötigen Parameter ein
    !!!Zufällige oder vorgefertigte Starttransformation und homogene, zufällige Punkteverteilung im W-Raum
    !!!Funktionszeiger f_prt zuweisen
    !!!Werte für center, height und width definieren
    !!!und was man sonst noch so braucht
    ! subroutine init_config(params, w_coords_array, x_coords_array, amount, dim)
    !     type(mapping), intent(inout)        :: params
    !     real(p), dimension(:,:), intent(inout) :: w_coords_array, x_coords_array
    !     integer, intent(inout)              :: amount, dim
    !     integer, parameter                  :: unit_number = 25

    !     f_ptr => f_value_peak
    !     open(unit=unit_number, file='config.txt', status='old', action='read')

    ! end subroutine init_config

end module module_initialize
