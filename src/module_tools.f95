module module_tools

    use module_matrix
    use module_transform
    use module_function
    use module_optimize
    use module_initialize
    use module_output

    implicit none

contains

    !!!Wurde nicht in der BA verwendet. Ziel war hier, Transformationen zu generieren, welche
    !!!die Punktemenge W auf die Punkte in center_coords_array "zusammenquetschen", sodass die resultierende
    !!!Menge X einen "Radius" von size hat. Hiermit sollte das Integrationsgebiet "abgetastet" werden,
    !!!um peaks an bestimmten Stellen zu erkennen. Dies sollte eine optimierte Vorgehensweise zum Finden
    !!!eines Startpunktes für die Transformation darstellen.
    subroutine scan_raster_random(x_coords_array, w_coords_array, params, attempts_counter, boundary, size)
        real(p), dimension(:,:), intent(inout) :: x_coords_array
        real(p), dimension(:,:), intent(in)    :: w_coords_array
        type(mapping), intent(inout)           :: params
        integer, intent(in)                    :: attempts_counter
        real(p), intent(in)                    :: boundary, size
        type(mapping)                          :: dummy_params
        real(p)                                :: opt_value_min, opt_value_current, bound
        real(p), allocatable, dimension(:,:)   :: center_coords_array
        integer                                :: i, j, k

        allocate(center_coords_array(1:attempts_counter,1:params%dim))
        call init_mapping_empty(dummy_params, params%dim)

        bound = atanh(boundary)
        center_coords_array = init_coords(attempts_counter, params%dim)

        call compute_raw_data(x_coords_array, w_coords_array, params)
        opt_value_min = opt_val()

        do k = 1, attempts_counter
            do i = 1, params%dim
                do j = 1, params%dim
                    if (i == j) then
                        if (abs(boundary - center_coords_array(k, i)) <= size) then
                            dummy_params%a(i, j) = bound / (  atanh(center_coords_array(k, i) + &
                                                & abs(boundary - center_coords_array(k, i))) &
                                                & - atanh(center_coords_array(k, i))  )
                            dummy_params%a_inv(i, j) = 1 / params%a(i, j)
                        else
                            dummy_params%a(i, j) = bound / (  atanh(center_coords_array(k, i) + size) &
                                                & - atanh(center_coords_array(k, i))  )
                            dummy_params%a_inv(i, j) = 1 / params%a(i, j)
                        end if
                    else
                        dummy_params%a(i, j)          = 0
                        dummy_params%a_inv(i, j) = 0
                    end if
                end do
                dummy_params%b(i) = - dummy_params%a(i, i) * atanh(center_coords_array(k, i))
            end do

            x_coords_array = dummy_params%backtransform_array(w_coords_array)
            call compute_raw_data(x_coords_array, w_coords_array, dummy_params)
            opt_value_current = opt_val()

            if(opt_value_current < opt_value_min) then
                opt_value_min = opt_value_current
                params = dummy_params
            end if

            print*, "optimizer current and min: ", opt_value_current, opt_value_min
        end do

        x_coords_array = params%backtransform_array(w_coords_array)
        call compute_raw_data(x_coords_array, w_coords_array, params)

    end subroutine scan_raster_random

    !!!Methode, um "manuell" eine Transformation zu erstellen, welche die Punktemenge
    !!!W auf den Punkt center_coords zusammenquetscht, sodass die resultierende Punktemenge
    !!!den Radius size_values hat (Radius in jede Koordinatenrichtung").
    subroutine init_trafo_peak(params, center_coords, size_values, boundary)
        type(mapping), intent(inout)      :: params
        real(p), dimension(:), intent(in) :: center_coords, size_values
        real(p), intent(in)               :: boundary
        real(p)                           :: bound
        integer                           :: i, j

        bound = atanh(boundary)
        do i = 1, size(center_coords)
            do j = 1, size(center_coords)
                if (i == j) then
                    if (abs(boundary - center_coords(i)) <= size_values(i)) then
                        params%a(i, j) = bound / (atanh(center_coords(i) + abs(boundary - center_coords(i))) &
                                        & - atanh(center_coords(i)))
                        params%a_inv(i, j) = 1 / params%a(i, j)
                    else
                        params%a(i, j) = bound / (atanh(center_coords(i) + size_values(i)) - atanh(center_coords(i)))
                        params%a_inv(i, j) = 1 / params%a(i, j)
                    end if
                else
                    params%a(i, j)          = 0
                    params%a_inv(i, j) = 0
                end if
            end do
            params%b(i) = - params%a(i, i) * atanh(center_coords(i))
        end do

    end subroutine init_trafo_peak

!!!Es folgen Prozeduren zur Berechnung der Integrale der verschiedenen samplings und des q-factors

    function q_factor()
        real(p) :: q_factor
        integer :: i

        q_factor = 0.0
        do i = 1, size(f_vals)
            q_factor = q_factor + abs(abs(f_vals(i) * det_vals(i)) - c)
        end do
        q_factor = q_factor / size(f_vals)

    end function q_factor

    subroutine int_N_data_unweighted(int_N, variance, dim)
        real(p), intent(inout) :: int_N, variance
        integer, intent(in)    :: dim
        integer                :: i
        real(p)                :: mult, ratio

        int_N    = 0.0
        variance = 0.0
        ratio    = size(f_vals) / (size(f_vals) - 1)

        do i = 1, size(f_vals)
            mult = 1 / real(i, p)
            int_N    = int_N + mult * (f_vals(i) * abs(det_vals(i)) - int_N)
        end do
        do i = 1, size(f_vals)
            mult = 1 / real(i, p)
            variance = variance + mult * ((f_vals(i) * abs(det_vals(i)) - int_N)**2 - variance)
        end do
        int_N    = (2**dim) * int_N
        variance = (4**dim) * variance * ratio

    end subroutine int_N_data_unweighted

    subroutine int_N_data_weighted_online(int_N, variance, dim)
        real(p), intent(inout) :: int_N, variance
        integer, intent(in)    :: dim
        real(p)                :: total_old, total, ratio, weight
        integer                :: i

        int_N     = 0.0
        variance  = 0.0
        ratio     = size(f_vals) / (size(f_vals) - 1)
        total     = 0.0
        total_old = 0.0

        do i = 1, size(f_vals)
            weight = 1 / (abs(g_vals(i)) * abs(g_vals(i)))
            if (weight > 10000.0) then
                weight = 10000.0
            end if
            total = total + weight
            int_N = int_N + (weight / total) * (f_vals(i) * abs(det_vals(i)) - int_N)
        end do
        total     = 0.0
        do i = 1, size(f_vals)
            weight = 1 / (abs(g_vals(i)) * abs(g_vals(i)))
            if (weight > 10000.0) then
                weight = 10000.0
            end if
            total = total + weight
            variance = variance + (weight / total) * ((f_vals(i) * abs(det_vals(i)) - int_N)**2 - variance)
        end do
        int_N = (2**dim) * int_N
        variance = (4**dim) * ratio * variance

    end subroutine int_N_data_weighted_online

    subroutine int_I_data_online(w_coords_array, int_I, variance)
        real(p), dimension(:,:), intent(in) :: w_coords_array
        real(p), intent(inout)              :: int_I, variance
        integer                             :: i
        real(p)                             :: mult, ratio

        int_I    = 0.0
        variance = 0.0
        ratio    = size(w_coords_array(:,1)) / (size(w_coords_array(:,1)) - 1)

        do i = 1, size(w_coords_array(:,1))
            mult  = 1 / real(i, p)
            int_I = int_I + mult * (f_ptr(w_coords_array(i,:)) - int_I)
        end do
        int_I = int_I
        do i = 1, size(w_coords_array(:,1))
            mult = 1 / real(i, p)
            variance = variance + mult * ((f_ptr(w_coords_array(i,:)) - int_I)**2 - variance)
        end do
        int_I = (2**size(w_coords_array(1,:))) * int_I
        variance = (4**size(w_coords_array(1,:))) * variance * ratio

    end subroutine int_I_data_online

    function integral_I(amount, dim)
        real(p)                              :: integral_I
        integer, intent(in)                  :: amount, dim
        real(p), allocatable, dimension(:,:) :: w_coords_array
        integer                              :: i
        real(p)                              :: mult

        allocate(w_coords_array(amount,dim))
        w_coords_array = init_coords(amount, dim)
        integral_I = 0.0
        do i = 1, size(w_coords_array(:,1))
            mult = 1 / real(i, p)
            integral_I = integral_I + mult * (f_ptr(w_coords_array(i,:)) - integral_I)
        end do
        integral_I = integral_I * (2**size(w_coords_array(1,:)))

    end function integral_I

end module module_tools
