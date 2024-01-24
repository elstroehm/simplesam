program test_algorithm

!!!Das erste Testprogramm. Wurde nicht für die BA verwendet.

    use module_matrix
    use module_transform
    use module_function
    use module_optimize
    use module_algorithm
    use module_output
    use module_initialize
    use module_tools

    implicit none

    integer, parameter :: q = selected_real_kind(6, 6)

    real(p), allocatable, dimension(:,:) :: w_coords_array, x_coords_array, dummy_coords, a_change
    type(mapping)                        :: params, params_difference, start_params
    integer                              :: dim, amount, cycles1, cycles2
    real(p)                              :: threshold, stepwidth, opt_min, f_max, g_mean, g_mean_old
    logical                              :: reinitialize
    integer                              :: i, count_reinitialize, attempts_counter, k, l, pos, neg, pos_old, neg_old
    real(p), allocatable, dimension(:)   :: det_deriv_values, f_deriv_values, h_derivs
    real(p), allocatable, dimension(:)   :: g_derivs, size_array, b_change, h_values
    real(p)                              :: gamma_val, boundary, size_val, new_stepwidth, variance, int_N
    real(p)                              :: b, integral
    real(q), allocatable, dimension(:,:) :: w, x
    logical, allocatable, dimension(:)   :: blabla

    dim              = 3
    amount           = 1000
    threshold        = 100.0
    stepwidth        = 1
    cycles1          = 1000
    cycles2          = 1000
    size_array       = [0.17, 0.17]
    new_stepwidth    = 1


    allocate(center(dim))
    allocate(h_derivs(1:amount), g_derivs(1:amount), b_change(1:dim), a_change(1:dim,1:dim), h_values(1:amount))

    center = [0.5, 0.2, 0.7]
    height = 1000000
    f_ptr => f_value_peak

    allocate(w_coords_array(1:amount,1:dim), x_coords_array(1:amount,1:dim), dummy_coords(1:amount,1:dim))
    allocate(w(1:amount,1:dim), x(1:amount,1:dim))

    call init_module_optimize(amount)

    h_vals_ptr   => compute_h_vals_asym
    h_derivs_ptr => compute_h_derivs_asym

    integral = integral_I(10000000, dim)
    c = integral / (2**dim)

    w_coords_array = init_coords(amount, dim)

    call init_mapping_id(params, dim)
    ! call random_params(params)
    ! call init_trafo_peak(params, center, size_array, boundary)

    x_coords_array = params%backtransform_array(w_coords_array)
    x = real(x_coords_array, q)

    ! call init_number_random(height)
    ! height = (height * 100000.0) + 1000.0
    ! b = 1 + (1 / height)
    ! integral = (real(2.0, p)**(dim-1)) * log(real(2.0, p) * height + real(1.0, p))
    ! c = integral / (2**dim)

    call compute_raw_data(x_coords_array, w_coords_array, params)

    gamma_val = opt_val()

    call optimize_coordinate_descent_inv(x_coords_array, w_coords_array, params, stepwidth, cycles1)

    ! h_vals_ptr   => compute_h_vals_abs
    ! h_derivs_ptr => compute_h_derivs_abs
    ! call optimize_coordinate_descent(x_coords_array, w_coords_array, params, new_stepwidth, cycles2)

    call compute_raw_data(x_coords_array, w_coords_array, params)

    print*, "optimizer at start:"
    print*, gamma_val
    print*, "______________________________________"

    print*, "optimizer"
    print*, opt_val()
    print*, "______________________________________"

    print*, "difference"
    print*, gamma_val - opt_val()
    print*, "______________________________________"

    ! print*, "starting parameters:"
    ! call print_params(start_params)
    ! print*, "______________________________________"
    !
    ! print*, "end parameters:"
    ! call print_params(params)
    ! print*, "______________________________________"

    print*, "integral, c"
    print*, integral, "  ", c
    print*, "______________________________________"

    int_N = 0.0
    variance = 0.0
    call int_N_data_weighted_online(int_N, variance, dim)

    print*, "importance sampling integral, estimated error and real error, weighted:"
    print*, int_N, "  ", sqrt(variance / real(amount)), "  ", abs(1 - int_N / ((2**dim) * c))
    print*, "______________________________________"

    int_N = 0.0
    variance = 0.0
    call int_N_data_unweighted(int_N, variance, dim)

    print*, "importance sampling integral, estimated error and real error, unweighted:"
    print*, int_N, "  ", sqrt(variance / real(amount)), "  ", abs(1 - int_N / ((2**dim) * c))
    print*, "______________________________________"

    int_N = 0.0
    variance = 0.0
    call int_I_data_online(w_coords_array, int_N, variance)

    print*, "uniform sampling integral, estimated error and real error:"
    print*, int_N, "  ", sqrt(variance / real(amount)), "  ", abs(1 - int_N / ((2**dim) * c))
    print*, "______________________________________"

    print*, "q-factor:"
    print*, q_factor()
    print*, "______________________________________"

    w = real(w_coords_array, q)
    x = real(x_coords_array, q)

    open(unit=20, file='../data/w_coords.txt', status='old', action='write')
    open(unit=30, file='../data/x_coords.txt', status='old', action='write')

    do i = 1, amount
        write(20, *) w(i,:)
        write(30, *) x(i,:)
    end do

    close(unit=20)
    close(unit=30)

end program test_algorithm
