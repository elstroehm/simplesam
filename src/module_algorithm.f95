module module_algorithm

    use module_matrix
    use module_transform
    use module_function
    use module_initialize
    use module_output
    use module_optimize

contains

    !!!Diese Prozedur wurde zur Aufzeichnung der q-factor Konvergenz verwendet. Es ist die gleiche Prozedur
    !!!wie optimize_coordinate_descent_inv, nur dass hier zusätzlich nach jedem cycle der q-factor in
    !!!eine Datei geschrieben wird.
    subroutine optimize_coordinate_descent_inv_record_q(x_coords_array, w_coords_array, params, stepwidth, cycles, filenumber)
        real(p), dimension(:,:), intent(inout) :: x_coords_array
        real(p), dimension(:,:), intent(in)    :: w_coords_array
        type(mapping), intent(inout)           :: params
        real(p), intent(in)                    :: stepwidth
        integer, intent(in)                    :: cycles, filenumber
        integer                                :: i, j, k
        real(p)                                :: opt_value_old, opt_value, opt_deriv, q_factor
        real(p), allocatable, dimension(:,:)   :: x_coords_old, a_steps
        type(mapping)                          :: params_old
        real(p), allocatable, dimension(:)     :: f_vals_old, det_vals_old, g_vals_old, b_steps

        allocate(x_coords_old(1:size(x_coords_array(:,1)), params%dim))
        allocate(f_vals_old(1:size(w_coords_array(:,1))), det_vals_old(1:size(w_coords_array(:,1))), &
                & g_vals_old(1:size(w_coords_array(:,1))))
        allocate(a_steps(1:params%dim,1:params%dim), b_steps(1:params%dim))
        call init_mapping_empty(params_old, params%dim)

        do i = 1, params%dim
            do j = 1, params%dim
                a_steps(i, j) = stepwidth
            end do
            b_steps(i) = stepwidth
        end do

        call compute_raw_data(x_coords_array, w_coords_array, params)
        opt_value_old       = opt_val()
        params_old          = params
        x_coords_old        = x_coords_array
        f_vals_old          = f_vals
        det_vals_old        = det_vals
        g_vals_old          = g_vals

        print*, k, " optimizer: ", opt_val(), b_steps(i)

        do k = 1, cycles
            do i = 1, params%dim
                do j = 1, params%dim
                    opt_deriv = opt_deriv_a_inv(x_coords_array, w_coords_array, params, i, j, a_steps(i, j))

                    if (opt_deriv < 0) then
                        params%a_inv(i, j) = params%a_inv(i, j) + a_steps(i, j)
                    else
                        params%a_inv(i, j) = params%a_inv(i, j) - a_steps(i, j)
                    end if

                    x_coords_array = params%backtransform_array(w_coords_array)
                    call compute_raw_data(x_coords_array, w_coords_array, params)
                    opt_value = opt_val()

                    if (opt_value < opt_value_old) then     !!!Wenn die Transformation besser ist, übernehme diese
                        opt_value_old = opt_value
                        params%a      = params%a_inv
                        call invert_matrix(params%a, params%dim)
                        params_old    = params
                        x_coords_old  = x_coords_array
                        f_vals_old    = f_vals
                        det_vals_old  = det_vals
                        g_vals_old    = g_vals
                        a_steps(i, j) = 1.024 * a_steps(i, j)
                    else                                        !!!ansonsten lasse sie unverändert
                        params         = params_old
                        x_coords_array = x_coords_old
                        f_vals         = f_vals_old
                        det_vals       = det_vals_old
                        g_vals         = g_vals_old
                        a_steps(i, j)  = 0.98 * a_steps(i, j)
                    end if
                end do
            end do

            do i = 1, params%dim
                opt_deriv = opt_deriv_b(x_coords_array, w_coords_array, params, i, b_steps(i))

                if (opt_deriv < 0) then
                    params%b(i) = params%b(i) + b_steps(i)
                else
                    params%b(i) = params%b(i) - b_steps(i)
                end if

                x_coords_array  = params%backtransform_array(w_coords_array)
                call compute_raw_data(x_coords_array, w_coords_array, params)
                opt_value = opt_val()

                if (opt_value < opt_value_old) then
                    opt_value_old = opt_value
                    params_old    = params
                    x_coords_old  = x_coords_array
                    f_vals_old    = f_vals
                    det_vals_old  = det_vals
                    g_vals_old    = g_vals
                    b_steps(i)    = 1.024 * b_steps(i)
                else
                    params         = params_old
                    x_coords_array = x_coords_old
                    f_vals         = f_vals_old
                    det_vals       = det_vals_old
                    g_vals         = g_vals_old
                    b_steps(i)     = 0.98 * b_steps(i)
                end if
            end do

            print*, k, " optimizer: ", opt_val()

            q_factor = 0.0
            do l = 1, size(f_vals)
                q_factor = q_factor + abs(f_vals(l) * abs(det_vals(l)) - c)
            end do
            q_factor = q_factor / size(f_vals)
            write(filenumber, *) q_factor

        end do

    end subroutine optimize_coordinate_descent_inv_record_q

    !!!Diese Prozedur wurde für alle Tests in der Bachelorarbeit verwendet. Es handelt sich um einen
    !!!coordinate descent Algorithmus, der die Transformationsparameter a^{-1}_ij und b_i optimiert.
    subroutine optimize_coordinate_descent_inv(x_coords_array, w_coords_array, params, stepwidth, cycles)
        real(p), dimension(:,:), intent(inout) :: x_coords_array
        real(p), dimension(:,:), intent(in)    :: w_coords_array
        type(mapping), intent(inout)           :: params
        real(p), intent(in)                    :: stepwidth
        integer, intent(in)                    :: cycles
        integer                                :: i, j, k
        real(p)                                :: opt_value_old, opt_value, opt_deriv
        real(p), allocatable, dimension(:,:)   :: x_coords_old, a_steps
        type(mapping)                          :: params_old
        real(p), allocatable, dimension(:)     :: f_vals_old, det_vals_old, g_vals_old, b_steps

        allocate(x_coords_old(1:size(x_coords_array(:,1)), params%dim))
        allocate(f_vals_old(1:size(w_coords_array(:,1))), det_vals_old(1:size(w_coords_array(:,1))), &
                & g_vals_old(1:size(w_coords_array(:,1))))
        allocate(a_steps(1:params%dim,1:params%dim), b_steps(1:params%dim))
        call init_mapping_empty(params_old, params%dim)

        !!!a_steps und b_steps mit den Anfangswerten belegen
        do i = 1, params%dim
            do j = 1, params%dim
                a_steps(i, j) = stepwidth
            end do
            b_steps(i) = stepwidth
        end do

        !!!Berechne Funktionswerte f(x), Jacobideterminanten det(J(w)), Werte für g(w, a) und
        !!!den Wert der Zielfunktion \gamma(a)
        call compute_raw_data(x_coords_array, w_coords_array, params)
        opt_value_old       = opt_val()
        params_old          = params
        x_coords_old        = x_coords_array
        f_vals_old          = f_vals
        det_vals_old        = det_vals
        g_vals_old          = g_vals

        ! print*, k, " optimizer: ", opt_val()

        do k = 1, cycles    !!!Anzahl der cycles
            do i = 1, params%dim    !!!Optimiere Matrixparameter a^{-1}_ij
                do j = 1, params%dim

                    !!!Berechne Ableitung der Zielfunktion \gamma(a)
                    opt_deriv = opt_deriv_a_inv(x_coords_array, w_coords_array, params, i, j, a_steps(i, j))

                    !!!Schritt in Richtung des Abstiegs
                    if (opt_deriv < 0) then
                        params%a_inv(i, j) = params%a_inv(i, j) + a_steps(i, j)
                    else
                        params%a_inv(i, j) = params%a_inv(i, j) - a_steps(i, j)
                    end if

                    !!!Nun Funktionswerte f(x), Jacobideterminanten det(J(w)), Werte für g(w, a) und \gamma(a)
                    !!!An der neuen Stelle berechnen.
                    x_coords_array = params%backtransform_array(w_coords_array)
                    call compute_raw_data(x_coords_array, w_coords_array, params)
                    opt_value = opt_val()

                    if (opt_value < opt_value_old) then     !!!Wenn die Transformation besser ist, übernehme diese
                        opt_value_old = opt_value
                        params%a      = params%a_inv
                        call invert_matrix(params%a, params%dim)
                        params_old    = params
                        x_coords_old  = x_coords_array
                        f_vals_old    = f_vals
                        det_vals_old  = det_vals
                        g_vals_old    = g_vals
                        a_steps(i, j) = 1.024 * a_steps(i, j)
                    else                                    !!!ansonsten lasse sie unverändert
                        params         = params_old
                        x_coords_array = x_coords_old
                        f_vals         = f_vals_old
                        det_vals       = det_vals_old
                        g_vals         = g_vals_old
                        a_steps(i, j)  = 0.98 * a_steps(i, j)
                    end if
                end do
            end do

            do i = 1, params%dim    !!!Optimiere Translationsparameter b_i
                opt_deriv = opt_deriv_b(x_coords_array, w_coords_array, params, i, b_steps(i))

                if (opt_deriv < 0) then
                    params%b(i) = params%b(i) + b_steps(i)
                else
                    params%b(i) = params%b(i) - b_steps(i)
                end if

                x_coords_array  = params%backtransform_array(w_coords_array)
                call compute_raw_data(x_coords_array, w_coords_array, params)
                opt_value = opt_val()

                if (opt_value < opt_value_old) then
                    opt_value_old = opt_value
                    params_old    = params
                    x_coords_old  = x_coords_array
                    f_vals_old    = f_vals
                    det_vals_old  = det_vals
                    g_vals_old    = g_vals
                    b_steps(i)    = 1.024 * b_steps(i)
                else
                    params         = params_old
                    x_coords_array = x_coords_old
                    f_vals         = f_vals_old
                    det_vals       = det_vals_old
                    g_vals         = g_vals_old
                    b_steps(i)     = 0.98 * b_steps(i)
                end if
            end do

            ! print*, k, " optimizer: ", opt_val()

        end do

    end subroutine optimize_coordinate_descent_inv

    !!!Wie optimize_coordinate_descent_inv, jedoch werden hier nicht die a^{-1}_ij optimiert,
    !!!sondern die a_ij. Wurde in der BA nicht verwendet.
    subroutine optimize_coordinate_descent(x_coords_array, w_coords_array, params, stepwidth, cycles)
        real(p), dimension(:,:), intent(inout) :: x_coords_array
        real(p), dimension(:,:), intent(in)    :: w_coords_array
        type(mapping), intent(inout)           :: params
        real(p), intent(in)                    :: stepwidth
        integer, intent(in)                    :: cycles
        integer                                :: i, j, k
        real(p)                                :: opt_value_old, opt_value, opt_deriv
        real(p), allocatable, dimension(:,:)   :: x_coords_old, a_steps
        type(mapping)                          :: params_old
        real(p), allocatable, dimension(:)     :: f_vals_old, det_vals_old, g_vals_old, b_steps

        allocate(x_coords_old(1:size(x_coords_array(:,1)), params%dim))
        allocate(f_vals_old(1:size(w_coords_array(:,1))), det_vals_old(1:size(w_coords_array(:,1))), &
                & g_vals_old(1:size(w_coords_array(:,1))))
        allocate(a_steps(1:params%dim,1:params%dim), b_steps(1:params%dim))
        call init_mapping_empty(params_old, params%dim)

        do i = 1, params%dim
            do j = 1, params%dim
                a_steps(i, j) = stepwidth
            end do
            b_steps(i) = stepwidth
        end do

        call compute_raw_data(x_coords_array, w_coords_array, params)
        opt_value_old       = opt_val()
        params_old          = params
        x_coords_old        = x_coords_array
        f_vals_old          = f_vals
        det_vals_old        = det_vals
        g_vals_old          = g_vals

        ! print*, k, " optimizer: ", opt_val(), b_steps(i)

        do k = 1, cycles
            do i = 1, params%dim
                do j = 1, params%dim
                    opt_deriv = opt_deriv_a(x_coords_array, w_coords_array, params, i, j, a_steps(i, j))

                    if (opt_deriv < 0) then
                        params%a(i, j) = params%a(i, j) + a_steps(i, j)
                    else
                        params%a(i, j) = params%a(i, j) - a_steps(i, j)
                    end if

                    params%a_inv = params%a
                    call invert_matrix(params%a_inv, params%dim)
                    x_coords_array = params%backtransform_array(w_coords_array)
                    call compute_raw_data(x_coords_array, w_coords_array, params)
                    opt_value = opt_val()

                    if (opt_value < opt_value_old) then     !!!Wenn die Transformation besser ist, übernehme diese
                        opt_value_old = opt_value
                        params_old    = params
                        x_coords_old  = x_coords_array
                        f_vals_old    = f_vals
                        det_vals_old  = det_vals
                        g_vals_old    = g_vals
                        a_steps(i, j) = 1.024 * a_steps(i, j)
                    else                                        !!!ansonsten lasse sie unverändert
                        params         = params_old
                        x_coords_array = x_coords_old
                        f_vals         = f_vals_old
                        det_vals       = det_vals_old
                        g_vals         = g_vals_old
                        a_steps(i, j)  = 0.98 * a_steps(i, j)
                    end if
                end do
            end do

            do i = 1, params%dim
                opt_deriv = opt_deriv_b(x_coords_array, w_coords_array, params, i, b_steps(i))

                if (opt_deriv < 0) then
                    params%b(i) = params%b(i) + b_steps(i)
                else
                    params%b(i) = params%b(i) - b_steps(i)
                end if

                x_coords_array  = params%backtransform_array(w_coords_array)
                call compute_raw_data(x_coords_array, w_coords_array, params)
                opt_value = opt_val()

                if (opt_value < opt_value_old) then
                    opt_value_old = opt_value
                    params_old    = params
                    x_coords_old  = x_coords_array
                    f_vals_old    = f_vals
                    det_vals_old  = det_vals
                    g_vals_old    = g_vals
                    b_steps(i)    = 1.024 * b_steps(i)
                else
                    params         = params_old
                    x_coords_array = x_coords_old
                    f_vals         = f_vals_old
                    det_vals       = det_vals_old
                    g_vals         = g_vals_old
                    b_steps(i)     = 0.98 * b_steps(i)
                end if
            end do

            ! print*, k, " optimizer: ", opt_val()

        end do

    end subroutine optimize_coordinate_descent

    !!!Ein direktes Gradientenabstiegsverfahren. Wurde in der BA nicht verwendet, funktioniert aber.
    !!!Es konvergiert wider Erwartens nicht so gut wie das coordinate descent Verfahren.
    subroutine optimize_gradient_inv(x_coords_array, w_coords_array, params, stepwidth, cycles)
        real(p), dimension(:,:), intent(inout) :: x_coords_array
        real(p), dimension(:,:), intent(in)    :: w_coords_array
        type(mapping), intent(inout)           :: params
        real(p), intent(inout)                 :: stepwidth
        integer, intent(in)                    :: cycles
        integer                                :: i, j, k
        real(p)                                :: opt_value_old, opt_value, abs_gradient
        real(p), allocatable, dimension(:,:)   :: x_coords_old, opt_gradient_a
        type(mapping)                          :: params_old
        real(p), allocatable, dimension(:)     :: f_vals_old, det_vals_old, g_vals_old, opt_gradient_b

        allocate(opt_gradient_a(1:params%dim,1:params%dim), x_coords_old(1:size(x_coords_array(:,1)), params%dim))
        allocate(f_vals_old(1:size(w_coords_array(:,1))), det_vals_old(1:size(w_coords_array(:,1))), &
                & g_vals_old(1:size(w_coords_array(:,1))), opt_gradient_b(1:params%dim))
        call init_mapping_empty(params_old, params%dim)

        call compute_raw_data(x_coords_array, w_coords_array, params)
        opt_value_old  = opt_val()
        params_old     = params
        x_coords_old   = x_coords_array
        f_vals_old     = f_vals
        det_vals_old   = det_vals
        g_vals_old     = g_vals
        abs_gradient   = 0.0

        do k = 1, cycles
            do i = 1, params%dim
                do j = 1, params%dim
                    opt_gradient_a(i, j) = opt_deriv_a_inv(x_coords_array, w_coords_array, params, i, j, stepwidth)
                    abs_gradient = abs_gradient + opt_gradient_a(i, j)**2
                end do
            end do
            do i = 1, params%dim
                opt_gradient_b(i) = opt_deriv_b(x_coords_array, w_coords_array, params, i, stepwidth)
                abs_gradient = abs_gradient + opt_gradient_b(i)**2
            end do
            abs_gradient = sqrt(abs_gradient)

            params%a_inv = params%a_inv - (stepwidth * opt_gradient_a / abs_gradient)
            params%b     = params%b - (stepwidth * opt_gradient_b / abs_gradient)

            x_coords_array = params%backtransform_array(w_coords_array)
            call compute_raw_data(x_coords_array, w_coords_array, params)
            opt_value = opt_val()

            if (opt_value < opt_value_old) then
                opt_value_old = opt_value
                params%a      = params%a_inv
                call invert_matrix(params%a, params%dim)
                params_old    = params
                x_coords_old  = x_coords_array
                f_vals_old    = f_vals
                det_vals_old  = det_vals
                g_vals_old    = g_vals
            else
                params         = params_old
                x_coords_array = x_coords_old
                f_vals         = f_vals_old
                det_vals       = det_vals_old
                g_vals         = g_vals_old
            end if

            if (opt_value >= opt_value_old) then
                stepwidth = 0.98 * stepwidth
            else
                stepwidth = 1.024 * stepwidth
            end if

            print*, k, " optimizer: ", opt_val()

        end do

    end subroutine  optimize_gradient_inv

!!!Alle folgenden Prozeduren sind gleicher Bauart wie die vorherigen, nur dass hier auf die
!!!numerisch Brechnung der Ableitung von det(J(w)) zurückgegriffen wurde. Diese Methode ist,
!!!wie in der BA beschrieben, ineffizienter.

    !!!Wurde in der BA verwendet, um den Laufzeitvergleich in höheren Dimensionen anzustellen.
    subroutine optimize_coordinate_descent_inv_num(x_coords_array, w_coords_array, params, stepwidth, cycles)
        real(p), dimension(:,:), intent(inout) :: x_coords_array
        real(p), dimension(:,:), intent(in)    :: w_coords_array
        type(mapping), intent(inout)           :: params
        real(p), intent(in)                    :: stepwidth
        integer, intent(in)                    :: cycles
        integer                                :: i, j, k
        real(p)                                :: opt_value_old, opt_value, opt_deriv
        real(p), allocatable, dimension(:,:)   :: x_coords_old, a_steps
        type(mapping)                          :: params_old
        real(p), allocatable, dimension(:)     :: f_vals_old, det_vals_old, g_vals_old, b_steps

        allocate(x_coords_old(1:size(x_coords_array(:,1)), params%dim))
        allocate(f_vals_old(1:size(w_coords_array(:,1))), det_vals_old(1:size(w_coords_array(:,1))), &
                & g_vals_old(1:size(w_coords_array(:,1))))
        allocate(a_steps(1:params%dim,1:params%dim), b_steps(1:params%dim))
        call init_mapping_empty(params_old, params%dim)

        do i = 1, params%dim
            do j = 1, params%dim
                a_steps(i, j) = stepwidth
            end do
            b_steps(i) = stepwidth
        end do

        call compute_raw_data(x_coords_array, w_coords_array, params)
        opt_value_old       = opt_val()
        params_old          = params
        x_coords_old        = x_coords_array
        f_vals_old          = f_vals
        det_vals_old        = det_vals
        g_vals_old          = g_vals

        do k = 1, cycles
            do i = 1, params%dim
                do j = 1, params%dim
                    opt_deriv = opt_deriv_a_inv_num(x_coords_array, w_coords_array, params, i, j, a_steps(i, j))

                    if (opt_deriv < 0) then
                        params%a_inv(i, j) = params%a_inv(i, j) + a_steps(i, j)
                    else
                        params%a_inv(i, j) = params%a_inv(i, j) - a_steps(i, j)
                    end if

                    x_coords_array = params%backtransform_array(w_coords_array)
                    call compute_raw_data(x_coords_array, w_coords_array, params)
                    opt_value = opt_val()

                    if (opt_value < opt_value_old) then     !!!Wenn die Transformation besser ist, übernehme diese
                        opt_value_old = opt_value
                        params%a      = params%a_inv
                        call invert_matrix(params%a, params%dim)
                        params_old    = params
                        x_coords_old  = x_coords_array
                        f_vals_old    = f_vals
                        det_vals_old  = det_vals
                        g_vals_old    = g_vals
                        a_steps(i, j) = 1.024 * a_steps(i, j)
                    else                                        !!!ansonsten lasse sie unverändert
                        params         = params_old
                        x_coords_array = x_coords_old
                        f_vals         = f_vals_old
                        det_vals       = det_vals_old
                        g_vals         = g_vals_old
                        a_steps(i, j)  = 0.98 * a_steps(i, j)
                    end if
                end do
            end do

            do i = 1, params%dim
                opt_deriv = opt_deriv_b_num(x_coords_array, w_coords_array, params, i, b_steps(i))

                if (opt_deriv < 0) then
                    params%b(i) = params%b(i) + b_steps(i)
                else
                    params%b(i) = params%b(i) - b_steps(i)
                end if

                x_coords_array  = params%backtransform_array(w_coords_array)
                call compute_raw_data(x_coords_array, w_coords_array, params)
                opt_value = opt_val()

                if (opt_value < opt_value_old) then
                    opt_value_old = opt_value
                    params_old    = params
                    x_coords_old  = x_coords_array
                    f_vals_old    = f_vals
                    det_vals_old  = det_vals
                    g_vals_old    = g_vals
                    b_steps(i)    = 1.024 * b_steps(i)
                else
                    params         = params_old
                    x_coords_array = x_coords_old
                    f_vals         = f_vals_old
                    det_vals       = det_vals_old
                    g_vals         = g_vals_old
                    b_steps(i)     = 0.98 * b_steps(i)
                end if
            end do

            ! print*, k, " optimizer: ", opt_val()

        end do

    end subroutine optimize_coordinate_descent_inv_num

    !!!Nicht in BA verwendet.
    subroutine optimize_coordinate_descent_num(x_coords_array, w_coords_array, params, stepwidth, cycles)
        real(p), dimension(:,:), intent(inout) :: x_coords_array
        real(p), dimension(:,:), intent(in)    :: w_coords_array
        type(mapping), intent(inout)           :: params
        real(p), intent(in)                    :: stepwidth
        integer, intent(in)                    :: cycles
        integer                                :: i, j, k
        real(p)                                :: opt_value_old, opt_value, opt_deriv
        real(p), allocatable, dimension(:,:)   :: x_coords_old, a_steps
        type(mapping)                          :: params_old
        real(p), allocatable, dimension(:)     :: f_vals_old, det_vals_old, g_vals_old, b_steps

        allocate(x_coords_old(1:size(x_coords_array(:,1)), params%dim))
        allocate(f_vals_old(1:size(w_coords_array(:,1))), det_vals_old(1:size(w_coords_array(:,1))), &
                & g_vals_old(1:size(w_coords_array(:,1))))
        allocate(a_steps(1:params%dim,1:params%dim), b_steps(1:params%dim))
        call init_mapping_empty(params_old, params%dim)

        do i = 1, params%dim
            do j = 1, params%dim
                a_steps(i, j) = stepwidth
            end do
            b_steps(i) = stepwidth
        end do

        call compute_raw_data(x_coords_array, w_coords_array, params)
        opt_value_old       = opt_val()
        params_old          = params
        x_coords_old        = x_coords_array
        f_vals_old          = f_vals
        det_vals_old        = det_vals
        g_vals_old          = g_vals

        do k = 1, cycles
            do i = 1, params%dim
                do j = 1, params%dim
                    opt_deriv = opt_deriv_a_num(x_coords_array, w_coords_array, params, i, j, a_steps(i, j))

                    if (opt_deriv < 0) then
                        params%a(i, j) = params%a(i, j) + a_steps(i, j)
                    else
                        params%a(i, j) = params%a(i, j) - a_steps(i, j)
                    end if

                    params%a_inv = params%a
                    call invert_matrix(params%a_inv, params%dim)
                    x_coords_array = params%backtransform_array(w_coords_array)
                    call compute_raw_data(x_coords_array, w_coords_array, params)
                    opt_value = opt_val()

                    if (opt_value < opt_value_old) then     !!!Wenn die Transformation besser ist, übernehme diese
                        opt_value_old = opt_value
                        params_old    = params
                        x_coords_old  = x_coords_array
                        f_vals_old    = f_vals
                        det_vals_old  = det_vals
                        g_vals_old    = g_vals
                        a_steps(i, j) = 1.024 * a_steps(i, j)
                    else                                        !!!ansonsten lasse sie unverändert
                        params         = params_old
                        x_coords_array = x_coords_old
                        f_vals         = f_vals_old
                        det_vals       = det_vals_old
                        g_vals         = g_vals_old
                        a_steps(i, j)  = 0.98 * a_steps(i, j)
                    end if
                end do
            end do

            do i = 1, params%dim
                opt_deriv = opt_deriv_b(x_coords_array, w_coords_array, params, i, b_steps(i))

                if (opt_deriv < 0) then
                    params%b(i) = params%b(i) + b_steps(i)
                else
                    params%b(i) = params%b(i) - b_steps(i)
                end if

                x_coords_array  = params%backtransform_array(w_coords_array)
                call compute_raw_data(x_coords_array, w_coords_array, params)
                opt_value = opt_val()

                if (opt_value < opt_value_old) then
                    opt_value_old = opt_value
                    params_old    = params
                    x_coords_old  = x_coords_array
                    f_vals_old    = f_vals
                    det_vals_old  = det_vals
                    g_vals_old    = g_vals
                    b_steps(i)    = 1.024 * b_steps(i)
                else
                    params         = params_old
                    x_coords_array = x_coords_old
                    f_vals         = f_vals_old
                    det_vals       = det_vals_old
                    g_vals         = g_vals_old
                    b_steps(i)     = 0.98 * b_steps(i)
                end if
            end do

            ! print*, k, " optimizer: ", opt_val()

        end do

    end subroutine optimize_coordinate_descent_num

    !!!Nicht in BA verwendet.
    subroutine optimize_gradient_inv_num(x_coords_array, w_coords_array, params, stepwidth, cycles)
        real(p), dimension(:,:), intent(inout) :: x_coords_array
        real(p), dimension(:,:), intent(in)    :: w_coords_array
        type(mapping), intent(inout)           :: params
        real(p), intent(inout)                 :: stepwidth
        integer, intent(in)                    :: cycles
        integer                                :: i, j, k
        real(p)                                :: opt_value_old, opt_value, abs_gradient
        real(p), allocatable, dimension(:,:)   :: x_coords_old, opt_gradient_a
        type(mapping)                          :: params_old
        real(p), allocatable, dimension(:)     :: f_vals_old, det_vals_old, g_vals_old, opt_gradient_b

        allocate(opt_gradient_a(1:params%dim,1:params%dim), x_coords_old(1:size(x_coords_array(:,1)), params%dim))
        allocate(f_vals_old(1:size(w_coords_array(:,1))), det_vals_old(1:size(w_coords_array(:,1))), &
                & g_vals_old(1:size(w_coords_array(:,1))), opt_gradient_b(1:params%dim))
        call init_mapping_empty(params_old, params%dim)

        call compute_raw_data(x_coords_array, w_coords_array, params)
        opt_value_old  = opt_val()
        params_old     = params
        x_coords_old   = x_coords_array
        f_vals_old     = f_vals
        det_vals_old   = det_vals
        g_vals_old     = g_vals
        abs_gradient   = 0.0

        do k = 1, cycles
            do i = 1, params%dim
                do j = 1, params%dim
                    opt_gradient_a(i, j) = opt_deriv_a_inv_num(x_coords_array, w_coords_array, params, i, j, stepwidth)
                    abs_gradient = abs_gradient + opt_gradient_a(i, j)**2
                end do
            end do
            do i = 1, params%dim
                opt_gradient_b(i) = opt_deriv_b_num(x_coords_array, w_coords_array, params, i, stepwidth)
                abs_gradient = abs_gradient + opt_gradient_b(i)**2
            end do
            abs_gradient = sqrt(abs_gradient)

            params%a_inv = params%a_inv - (stepwidth * opt_gradient_a / abs_gradient)
            params%b     = params%b - (stepwidth * opt_gradient_b / abs_gradient)

            x_coords_array = params%backtransform_array(w_coords_array)
            call compute_raw_data(x_coords_array, w_coords_array, params)
            opt_value = opt_val()

            if (opt_value < opt_value_old) then
                opt_value_old = opt_value
                params%a      = params%a_inv
                call invert_matrix(params%a, params%dim)
                params_old    = params
                x_coords_old  = x_coords_array
                f_vals_old    = f_vals
                det_vals_old  = det_vals
                g_vals_old    = g_vals
            else
                params         = params_old
                x_coords_array = x_coords_old
                f_vals         = f_vals_old
                det_vals       = det_vals_old
                g_vals         = g_vals_old
            end if

            if (opt_value >= opt_value_old) then
                stepwidth = 0.98 * stepwidth
            else
                stepwidth = 1.024 * stepwidth
            end if

            print*, k, " optimizer: ", opt_val()

        end do

    end subroutine  optimize_gradient_inv_num

end module module_algorithm
