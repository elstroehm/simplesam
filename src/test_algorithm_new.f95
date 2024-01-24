program test_algorithm

!!!Mit diesem Modul wurden alle Tests für die BA durchgeführt.

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

    real(p), allocatable, dimension(:,:) :: w_coords_array, x_coords_array, dummy_coords
    type(mapping)                        :: params
    integer                              :: dim, amount, cycles1, cycles2
    real(p)                              :: stepwidth1, stepwidth2
    integer                              :: k, i
    real(p)                              :: int, variance, integral, gamma_val, quality, start, finish
    real(p)                              :: int1, int2, var1, var2, int3, var3, b, i1, i2, i3, i4, i5
    real(q), allocatable, dimension(:,:) :: w, x
    real, dimension(4)                   :: data_array
    real, dimension(10)                  :: data_arr

    dim        = 2
    amount     = 1000
    stepwidth1 = 1
    stepwidth2 = 0.1
    cycles1    = 3000
    cycles2    = 10000
    ! center     = [0.5, 0.5]
    height     = 10000
    f_ptr      => f_value_wall

    allocate(w_coords_array(1:amount,1:dim), x_coords_array(1:amount,1:dim), dummy_coords(1:amount,1:dim))
    allocate(w(1:amount,1:dim), x(1:amount,1:dim))
    call init_module_optimize(amount)

    h_vals_ptr   => compute_h_vals_asym
    h_derivs_ptr => compute_h_derivs_asym

    dummy_coords = init_coords(1000000, dim)
    call set_c_integral(dummy_coords)

    w_coords_array = init_coords(amount, dim)

    call init_mapping_id(params, dim)
    ! params%a_inv(1,:) = [100.0, 0.0]
    ! params%a_inv(2,:) = [0.0, 100.0]
    ! params%a = params%a_inv
    ! call invert_matrix(params%a, dim)

    x_coords_array = params%backtransform_array(w_coords_array)
    x = real(x_coords_array, q)

    call compute_raw_data(x_coords_array, w_coords_array, params)

    do i = 1, amount
        print*, det_vals(i)
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! open(unit=20, file='../data/q_factor.txt', status='old', action='write')
    ! call optimize_coordinate_descent_inv_record_q(x_coords_array, w_coords_array, params, stepwidth1, cycles1, 20)
    ! close(unit=20)

    ! call optimize_coordinate_descent_inv(x_coords_array, w_coords_array, params, stepwidth1, cycles1)
    ! h_vals_ptr   => compute_h_vals_abs
    ! h_derivs_ptr => compute_h_derivs_abs
    ! call optimize_coordinate_descent_inv(x_coords_array, w_coords_array, params, stepwidth2, cycles2)

    ! call compute_raw_data(x_coords_array, w_coords_array, params)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!Sourcecode für die variance reduction bei Erhöhung der Anzahl der grid points

    ! integral = 28.1708
    ! integral = 19.807

    ! open(unit=20, file='../data/I_data.txt', status='old', action='write')
    ! open(unit=30, file='../data/N_data.txt', status='old', action='write')
    ! open(unit=40, file='../data/U_data.txt', status='old', action='write')
    ! do k = 10, 1000, 1
    !     deallocate(f_vals, det_vals, g_vals)
    !     deallocate(w_coords_array, x_coords_array)
    !     allocate(w_coords_array(1:k**2,dim), x_coords_array(1:k**2,dim))
    !     call init_module_optimize(k**2)
    !     call init_coords_random(w_coords_array)
    !     x_coords_array = params%backtransform_array(w_coords_array)
    !     call compute_raw_data(x_coords_array, w_coords_array, params)
    !     int      = 0.0
    !     variance = 0.0
    !     call int_I_data_online(w_coords_array, int, variance)
    !     data_array = [real(k), real(int), real(sqrt(variance/(k**2))), real(abs(1 - (int / integral)))]
    !     write(20, *) data_array
    !     int      = 0.0
    !     variance = 0.0
    !     call int_N_data_weighted_online(int, variance, dim)
    !     data_array = [real(k), real(int), real(sqrt(variance/(k**2))), real(abs(1 - (int / integral)))]
    !     write(30, *) data_array
    !     int      = 0.0
    !     variance = 0.0
    !     call int_N_data_unweighted(int, variance, dim)
    !     data_array = [real(k), real(int), real(sqrt(variance/(k**2))), real(abs(1 - (int / integral)))]
    !     write(40, *) data_array
    !     int      = 0.0
    !     variance = 0.0
    !     print*, k    
    ! end do
    ! close(unit=20)
    ! close(unit=30)
    ! close(unit=40)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!Sourcecode für die "consecutive runs"

    ! open(unit=20, file='../data/consecutive_runs.txt', status='old', action='write')
    ! do k = 1, 100
    !     int1 = 0.0
    !     int2 = 0.0
    !     int3 = 0.0
    !     var1 = 0.0
    !     var2 = 0.0
    !     var3 = 0.0
    !     call init_coords_random(w_coords_array)
    !     call init_id(params)
    !     x_coords_array = w_coords_array 
    !     ! call init_array_random(center)
    !     ! integral = integral_I(10000000, dim)
    !     call init_number_random(height)
    !     height = (height * 100000.0) + 1000.0
    !     b = 1 + (1 / height)
    !     integral = real(2.0, p) * log(real(2.0, p) * height + real(1.0, p))
    !     c = integral / (2**dim)
    !     call compute_raw_data(x_coords_array, w_coords_array, params)
    !     call optimize_coordinate_descent_inv(x_coords_array, w_coords_array, params, real(1.0, p), 1000)
    !     call int_N_data_weighted_online(int1, var1, dim)
    !     call int_N_data_unweighted(int2, var2, dim)
    !     call int_I_data_online(w_coords_array, int3, var3)
    !     data_arr = [real(integral), real(int1), real(sqrt(var1/k)),  real(abs(1 - (int1 / integral))), &
    !     & real(int2), real(sqrt(var2/k)), real(abs(1 - (int2 / integral))), &
    !     & real(int3), real(sqrt(var3/k)), real(abs(1 - (int3 / integral)))]
    !     gamma_val = opt_val()
    !     quality   = q_factor()
    !     ! write(20, *) real(k), center, gamma_val, quality
    !     write(20, *) real(k), b, gamma_val, quality
    !     write(20, *) data_arr
    !     print*, k
    ! end do
    ! close(unit=20)

!!!Sourcecode für die Tests in höheren Dimensionen

    ! open(unit=20, file='../data/consecutive_runs.txt', status='old', action='write')
    ! center = [0.5, 0.2, 0.7]
    ! height = 1000000.0
    ! f_ptr  => f_value_peak
    ! b = 0.0
    ! i1 = integral_I(10000000, dim)
    ! i2 = integral_I(10000000, dim)
    ! i3 = integral_I(10000000, dim)
    ! i4 = integral_I(10000000, dim)
    ! i5 = integral_I(10000000, dim)
    ! integral = (i1 + i2 + i3 + i4 + i5) / 5
    ! c = integral / (2**dim)
    ! do k = 1, 5
    !     int1 = 0.0
    !     int2 = 0.0
    !     int3 = 0.0
    !     var1 = 0.0
    !     var2 = 0.0
    !     var3 = 0.0
    !     call init_coords_random(w_coords_array)
    !     call init_id(params)
    !     x_coords_array = w_coords_array
    !     call compute_raw_data(x_coords_array, w_coords_array, params)
    !     call optimize_coordinate_descent_inv(x_coords_array, w_coords_array, params, real(1.0, p), 1000)
    !     call int_N_data_weighted_online(int1, var1, dim)
    !     call int_N_data_unweighted(int2, var2, dim)
    !     call int_I_data_online(w_coords_array, int3, var3)
    !     data_arr = [real(integral), real(int1), real(sqrt(var1/k)),  real(abs(1 - (int1 / integral))), &
    !     & real(int2), real(sqrt(var2/k)), real(abs(1 - (int2 / integral))), &
    !     & real(int3), real(sqrt(var3/k)), real(abs(1 - (int3 / integral)))]
    !     gamma_val = opt_val()
    !     quality   = q_factor()
    !     write(20, *) real(k), gamma_val, quality
    !     write(20, *) data_arr
    !     print*, k
    ! end do
    ! f_ptr  => f_value_wall
    ! b = 1 + (1 / height)
    ! integral = (2**(dim-1)) * log(real(2.0, p) * height + real(1.0, p))
    ! c = integral / (2**dim)
    ! do k = 1, 5
    !     int1 = 0.0
    !     int2 = 0.0
    !     int3 = 0.0
    !     var1 = 0.0
    !     var2 = 0.0
    !     var3 = 0.0
    !     call init_coords_random(w_coords_array)
    !     call init_id(params)
    !     x_coords_array = w_coords_array 
    !     call compute_raw_data(x_coords_array, w_coords_array, params)
    !     call optimize_coordinate_descent_inv(x_coords_array, w_coords_array, params, real(1.0, p), 1000)
    !     call int_N_data_weighted_online(int1, var1, dim)
    !     call int_N_data_unweighted(int2, var2, dim)
    !     call int_I_data_online(w_coords_array, int3, var3)
    !     data_arr = [real(integral), real(int1), real(sqrt(var1/k)),  real(abs(1 - (int1 / integral))), &
    !     & real(int2), real(sqrt(var2/k)), real(abs(1 - (int2 / integral))), &
    !     & real(int3), real(sqrt(var3/k)), real(abs(1 - (int3 / integral)))]
    !     gamma_val = opt_val()
    !     quality   = q_factor()
    !     write(20, *) real(k), gamma_val, quality
    !     write(20, *) data_arr
    !     print*, k
    ! end do
    ! center = [0.5, 0.5, 0.5]
    ! f_ptr  => f_value_cylinder
    ! i1 = integral_I(10000000, dim)
    ! i2 = integral_I(10000000, dim)
    ! i3 = integral_I(10000000, dim)
    ! i4 = integral_I(10000000, dim)
    ! i5 = integral_I(10000000, dim)
    ! integral = (i1 + i2 + i3 + i4 + i5) / 5
    ! c = integral / (2**dim)
    ! do k = 1, 5
    !     int1 = 0.0
    !     int2 = 0.0
    !     int3 = 0.0
    !     var1 = 0.0
    !     var2 = 0.0
    !     var3 = 0.0
    !     call init_coords_random(w_coords_array)
    !     call init_id(params)
    !     x_coords_array = w_coords_array 
    !     call compute_raw_data(x_coords_array, w_coords_array, params)
    !     call optimize_coordinate_descent_inv(x_coords_array, w_coords_array, params, real(1.0, p), 1000)
    !     call int_N_data_weighted_online(int1, var1, dim)
    !     call int_N_data_unweighted(int2, var2, dim)
    !     call int_I_data_online(w_coords_array, int3, var3)
    !     data_arr = [real(integral), real(int1), real(sqrt(var1/k)),  real(abs(1 - (int1 / integral))), &
    !     & real(int2), real(sqrt(var2/k)), real(abs(1 - (int2 / integral))), &
    !     & real(int3), real(sqrt(var3/k)), real(abs(1 - (int3 / integral)))]
    !     gamma_val = opt_val()
    !     quality   = q_factor()
    !     write(20, *) real(k), gamma_val, quality
    !     write(20, *) data_arr
    !     print*, k
    ! end do
    ! close(unit=20)

    ! integral = 29.4958 !(0.2,0.1)
    ! integral = 28.1708 !(0.5,0.5)
    ! integral = 24.5790 !(0.8,0.8)
    ! integral = 24.0550 !(0.6,0.9)

    ! print*, "optimizer at start, optimizer at end, difference:"
    ! print*, gamma_val, opt_val(), gamma_val - opt_val()
    ! print*, "______________________________________"

    ! print*, "integral value:"
    ! print*, integral
    ! print*, "______________________________________"

    ! call int_N_data_weighted(int_N, variance_N, dim)

    ! print*, "optimized sampling integral, estimated error and real error, weighted:"
    ! print*, int_N, "  ", sqrt(variance_N / real(amount)), "  ", abs(1 - int_N / ((2**dim) * c))
    ! print*, "______________________________________"

    ! call int_N_data(int_N, variance_N, dim)

    ! print*, "optimized sampling integral, estimated error and real error, unweighted:"
    ! print*, int_N, "  ", sqrt(variance_N / real(amount)), "  ", abs(1 - int_N / ((2**dim) * c))
    ! print*, "______________________________________"

    ! call int_I_data(w_coords_array, int_N, variance_I)

    ! print*, "uniform sampling integral, estimated error and real error, unweighted:"
    ! print*, int_N, "  ", sqrt(variance_I / real(amount)), "  ", abs(1 - int_N / ((2**dim) * c))
    ! print*, "______________________________________"

    ! print*, "q-factor:"
    ! print*, q_factor()
    ! print*, "______________________________________"

!!!Sourcecode für die Laufzeitmessung in höheren Dimensionen

    ! center = [0.5, 0.2, 0.7]
    ! height = 1000000.0
    ! f_ptr  => f_value_peak
    ! integral = 12.37
    ! c = integral / (2**dim)

    ! center = [0.5, 0.2, 0.7]
    ! height = 1000000.0
    ! f_ptr  => f_value_wall
    ! integral = 58.03
    ! c = integral / (2**dim)

    ! center = [0.5, 0.5, 0.5]
    ! height = 1000000.0
    ! f_ptr  => f_value_cylinder
    ! integral = 84.86
    ! c = integral / (2**dim)

    ! call init_coords_random(w_coords_array)
    ! call init_id(params)
    ! call cpu_time(start)
    ! call compute_raw_data(x_coords_array, w_coords_array, params)
    ! call optimize_coordinate_descent_inv_num(x_coords_array, w_coords_array, params, real(0.001, p), 100)
    ! call cpu_time(finish)
    ! print*, finish - start

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