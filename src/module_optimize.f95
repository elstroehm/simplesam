module module_optimize

    use module_matrix
    use module_transform
    use module_function
    use module_initialize
    use module_output

    implicit none

    abstract interface
        subroutine compute_h_vals(h_vals)
            integer, parameter :: p = selected_real_kind(12, 30)
            real(p), dimension(:), intent(inout) :: h_vals

        end subroutine compute_h_vals
    end interface

    abstract interface
        subroutine compute_h_derivs(h_derivs, g_derivs)
            integer, parameter :: p = selected_real_kind(12, 30)
            real(p), dimension(:), intent(inout) :: h_derivs
            real(p), dimension(:), intent(in)    :: g_derivs

        end subroutine compute_h_derivs
    end interface

    !!!Mit dem Zeiger h_ptr kann auf verschiedene Prozeduren für die Berechnung von h(g(w, a)) verwiesen werden,
    !!!sodass diese direkt im Optimizer verwendet werden können. Hier können dann verschiedene Ansätze direkt
    !!!implementiert werden.
    procedure(compute_h_vals), pointer   :: h_vals_ptr
    procedure(compute_h_derivs), pointer :: h_derivs_ptr

    real(p), allocatable, dimension(:) :: f_vals, det_vals, g_vals
    real(p)                            :: c

contains

!!!_____________________modulinterne Variablen definieren_________________________

    !!!Alloziere Speicher für die Arrays f_vals, det_vals und g_vals
    subroutine init_module_optimize(amount)
        integer, intent(in) :: amount

        allocate(f_vals(1:amount), det_vals(1:amount), g_vals(1:amount))
    end subroutine init_module_optimize

    subroutine set_c_integral(w_coords_array)
        real(p), dimension(:,:), intent(in) :: w_coords_array
        integer                             :: i

        c = 0.0
        do i = 1, size(w_coords_array(:,1))
            c = c + f_ptr(w_coords_array(i,:))
        end do

        c = c / size(w_coords_array(:,1))

    end subroutine set_c_integral

!!!_____________________________________Verschiedene Prozeduren für h(g(w, a))_____________________________________

    subroutine compute_h_vals_abs(h_vals)
        real(p), dimension(:), intent(inout) :: h_vals
        integer                              :: i

        do i = 1, size(g_vals)
            h_vals(i)   = abs(g_vals(i))
        end do

    end subroutine compute_h_vals_abs

    subroutine compute_h_derivs_abs(h_derivs, g_derivs)
        real(p), dimension(:), intent(inout) :: h_derivs
        real(p), dimension(:), intent(in)    :: g_derivs
        integer                              :: i

        do i = 1, size(g_vals)
            if (g_vals(i) < 0) then
                h_derivs(i) = - g_derivs(i)
            else
                h_derivs(i) = g_derivs(i)
            end if
        end do

    end subroutine compute_h_derivs_abs

    subroutine compute_h_vals_square(h_vals)
        real(p), dimension(:), intent(inout) :: h_vals
        integer                              :: i

        do i = 1, size(g_vals)
            h_vals(i) = g_vals(i) * g_vals(i)
        end do

    end subroutine compute_h_vals_square

    subroutine compute_h_derivs_square(h_derivs, g_derivs)
        real(p), dimension(:), intent(inout) :: h_derivs
        real(p), dimension(:), intent(in)    :: g_derivs
        integer                              :: i

        do i = 1, size(g_vals)
            h_derivs(i) = 2 * g_vals(i) * g_derivs(i)
        end do

    end subroutine compute_h_derivs_square

    subroutine compute_h_vals_sqrt(h_vals)
        real(p), dimension(:), intent(inout) :: h_vals
        integer                              :: i

        do i = 1, size(g_vals)
            h_vals(i) = sqrt(abs(g_vals(i)))
        end do

    end subroutine compute_h_vals_sqrt

    subroutine compute_h_derivs_sqrt(h_derivs, g_derivs)
        real(p), dimension(:), intent(inout) :: h_derivs
        real(p), dimension(:), intent(in)    :: g_derivs
        integer                              :: i

        do i = 1, size(g_vals)
            if (g_vals(i) > 0) then
                h_derivs(i) = 0.5 / sqrt(g_vals(i)) * g_derivs(i)
            else
                h_derivs(i) = -0.5 / sqrt(-g_vals(i)) * g_derivs(i)
            end if
        end do

    end subroutine compute_h_derivs_sqrt

    subroutine compute_h_vals_log(h_vals)
        real(p), dimension(:), intent(inout) :: h_vals
        integer                              :: i

        do i = 1, size(g_vals)
            h_vals(i) = log(abs(g_vals(i)) + 1.0) - 1.0
        end do

    end subroutine compute_h_vals_log

    subroutine compute_h_derivs_log(h_derivs, g_derivs)
        real(p), dimension(:), intent(inout) :: h_derivs
        real(p), dimension(:), intent(in)    :: g_derivs
        integer                              :: i

        do i = 1, size(g_vals)
            h_derivs(i) = (1 / (abs(g_vals(i)) + 1)) * g_derivs(i)
        end do

    end subroutine compute_h_derivs_log

    subroutine compute_h_vals_asym(h_vals)
        real(p), dimension(:), intent(inout) :: h_vals
        integer                              :: i

        do i = 1, size(g_vals)
            h_vals(i) = (1 / (g_vals(i) + c)) * g_vals(i) * g_vals(i)
        end do

    end subroutine compute_h_vals_asym

    subroutine compute_h_derivs_asym(h_derivs, g_derivs)
        real(p), dimension(:), intent(inout) :: h_derivs
        real(p), dimension(:), intent(in)    :: g_derivs
        integer                              :: i

        do i = 1, size(g_vals)
            h_derivs(i) = ( g_vals(i) / (g_vals(i) + c) ) *( 1 + (c / (g_vals(i) + c)) ) * g_derivs(i)
        end do

    end subroutine compute_h_derivs_asym

    subroutine compute_h_vals_asym2(h_vals)
        real(p), dimension(:), intent(inout) :: h_vals
        integer                              :: i

        do i = 1, size(g_vals)
            h_vals(i) = (1 / (g_vals(i) + c)**2) * g_vals(i)**4
        end do

    end subroutine compute_h_vals_asym2

    subroutine compute_h_derivs_asym2(h_derivs, g_derivs)
        real(p), dimension(:), intent(inout) :: h_derivs
        real(p), dimension(:), intent(in)    :: g_derivs
        integer                              :: i

        do i = 1, size(g_vals)
            h_derivs(i) = ( g_vals(i)**3 / (g_vals(i) + c)**2 ) *( 2 + (2 * c / (g_vals(i) + c)) ) * g_derivs(i)
        end do

    end subroutine compute_h_derivs_asym2

!!!Prozeduren für die Berechnung von g(w, a) und deren Ableitung, an allen Punkten im W-Raum
!!!Diese schreiben allesamt in die Arrays f_vals, det_vals und g_vals
!!!Und noch ein paar weitere Funktionen

    !!!Schreibt die entsprechenden Werte für die Funktionswerte, Jacobideterminanten und g(w, a) in
    !!!die dafür vorgesehenen Arrays.
    subroutine compute_raw_data(x_coords_array, w_coords_array, params)
        real(p), dimension(:,:), intent(in)  :: w_coords_array, x_coords_array
        type(mapping), intent(in)            :: params
        integer                              :: k

        do k = 1, size(w_coords_array(:,1))
            f_vals(k)   = f_ptr(x_coords_array(k,:))
            det_vals(k) = determinant(params%jacobian_backtrafo(x_coords_array(k,:), w_coords_array(k,:)), params%dim)
            g_vals(k)   = abs(f_vals(k) * det_vals(k)) - c
        end do

    end subroutine compute_raw_data

!!!Namensschema für die kommenden Prozeduren für die Ableitung von g(w, a): ohne "num" bedeutet,
!!!dass zur Berechnung der Ableitung der Jacobideterminante die analytischen Ausdrücke aus Appendix A
!!!mittels der dafür vorgesehenen Prozeduren in module_transform verwendet werden, mit "num" bedeutet,
!!!dass die Berechnung numerisch vorgenommen wird.

    subroutine compute_g_derivs_a(g_derivs, x_coords_array, w_coords_array, params, i, j, stepwidth)
        real(p), dimension(:), intent(inout) :: g_derivs
        real(p), dimension(:,:), intent(in)  :: w_coords_array, x_coords_array
        type(mapping), intent(in)            :: params
        integer, intent(in)                  :: i, j
        real(p), intent(in)                  :: stepwidth
        real(p), allocatable, dimension(:)   :: f_derivs, det_derivs
        integer                              :: k

        allocate(f_derivs(1:size(f_vals)), det_derivs(1:size(det_vals)))
        f_derivs   = f_derivative_a_array(params, w_coords_array, i, j, stepwidth)
        det_derivs = params%det_derivative_a_array(x_coords_array, w_coords_array, det_vals, i, j)

        do k = 1, size(w_coords_array(:,1))

            if (det_vals(k) < 0) then
                if (f_vals(k) < 0) then
                    g_derivs(k) = f_derivs(k) * det_vals(k) + det_derivs(k) * f_vals(k)
                else
                    g_derivs(k) = - (f_derivs(k) * det_vals(k) + det_derivs(k) * f_vals(k))
                end if
            else
                if (f_vals(k) < 0) then
                    g_derivs(k) = - (f_derivs(k) * det_vals(k) + det_derivs(k) * f_vals(k))
                else
                    g_derivs(k) = f_derivs(k) * det_vals(k) + det_derivs(k) * f_vals(k)
                end if
            end if
        end do

    end subroutine compute_g_derivs_a

    subroutine compute_g_derivs_a_inv(g_derivs, x_coords_array, w_coords_array, params, i, j, stepwidth)
        real(p), dimension(:), intent(inout) :: g_derivs
        real(p), dimension(:,:), intent(in)  :: w_coords_array, x_coords_array
        type(mapping), intent(in)            :: params
        integer, intent(in)                  :: i, j
        real(p), intent(in)                  :: stepwidth
        real(p), allocatable, dimension(:)   :: f_derivs, det_derivs
        integer                              :: k

        allocate(f_derivs(1:size(f_vals)), det_derivs(1:size(det_vals)))
        f_derivs   = f_derivative_a_inv_array(params, w_coords_array, i, j, stepwidth)
        det_derivs = params%det_derivative_a_inv_array(x_coords_array, w_coords_array, det_vals, i, j)

        do k = 1, size(w_coords_array(:,1))

            if (det_vals(k) < 0) then
                if (f_vals(k) < 0) then
                    g_derivs(k) = f_derivs(k) * det_vals(k) + det_derivs(k) * f_vals(k)
                else
                    g_derivs(k) = - (f_derivs(k) * det_vals(k) + det_derivs(k) * f_vals(k))
                end if
            else
                if (f_vals(k) < 0) then
                    g_derivs(k) = - (f_derivs(k) * det_vals(k) + det_derivs(k) * f_vals(k))
                else
                    g_derivs(k) = f_derivs(k) * det_vals(k) + det_derivs(k) * f_vals(k)
                end if
            end if
        end do

    end subroutine compute_g_derivs_a_inv

    subroutine compute_g_derivs_b(g_derivs, x_coords_array, w_coords_array, params, i, stepwidth)
        real(p), dimension(:), intent(inout) :: g_derivs
        real(p), dimension(:,:), intent(in)  :: w_coords_array, x_coords_array
        type(mapping), intent(in)            :: params
        integer, intent(in)                  :: i
        real(p), intent(in)                  :: stepwidth
        real(p), allocatable, dimension(:)   :: f_derivs, det_derivs
        integer                              :: k

        allocate(f_derivs(1:size(f_vals)), det_derivs(1:size(det_vals)))
        f_derivs   = f_derivative_b_array(params, w_coords_array, i, stepwidth)
        det_derivs = params%det_derivative_b_array(x_coords_array, det_vals, i)

        do k = 1, size(w_coords_array(:,1))

            if (det_vals(k) < 0) then
                if (f_vals(k) < 0) then
                    g_derivs(k) = f_derivs(k) * det_vals(k) + det_derivs(k) * f_vals(k)
                else
                    g_derivs(k) = - (f_derivs(k) * det_vals(k) + det_derivs(k) * f_vals(k))
                end if
            else
                if (f_vals(k) < 0) then
                    g_derivs(k) = - (f_derivs(k) * det_vals(k) + det_derivs(k) * f_vals(k))
                else
                    g_derivs(k) = f_derivs(k) * det_vals(k) + det_derivs(k) * f_vals(k)
                end if
            end if
        end do

    end subroutine compute_g_derivs_b

    subroutine compute_g_derivs_a_num(g_derivs, x_coords_array, w_coords_array, params, i, j, stepwidth)
        real(p), dimension(:), intent(inout) :: g_derivs
        real(p), dimension(:,:), intent(in)  :: w_coords_array, x_coords_array
        type(mapping), intent(in)            :: params
        integer, intent(in)                  :: i, j
        real(p), intent(in)                  :: stepwidth
        type(mapping)                        :: params_1, params_2
        real(p), allocatable, dimension(:,:) :: x_coords_1, x_coords_2
        real(p)                              :: f_1, f_2, det_1, det_2
        integer                              :: k
        real(p)                              :: f_deriv, det_deriv

        allocate(x_coords_1(1:size(w_coords_array(:,1)),1:params%dim), x_coords_2(1:size(w_coords_array(:,1)),1:params%dim))

        call init_mapping_empty(params_1, params%dim)
        call init_mapping_empty(params_2, params%dim)
        params_1 = params
        params_2 = params
        params_1%a(i, j) = params_1%a(i, j) - stepwidth
        params_2%a(i, j) = params_2%a(i, j) + stepwidth
        params_1%a_inv = params_1%a
        params_2%a_inv = params_2%a
        call invert_matrix(params_1%a_inv, params%dim)
        call invert_matrix(params_2%a_inv, params%dim)

        x_coords_1 = params_1%backtransform_array(w_coords_array)
        x_coords_2 = params_2%backtransform_array(w_coords_array)

        do k = 1, size(w_coords_array(:,1))
            f_1   = f_ptr(x_coords_1(k,:))
            f_2   = f_ptr(x_coords_2(k,:))
            det_1 = determinant( params_1%jacobian_backtrafo(x_coords_1(k,:), w_coords_array(k,:)), params%dim)
            det_2 = determinant( params_2%jacobian_backtrafo(x_coords_2(k,:), w_coords_array(k,:)), params%dim)

            f_deriv   = (f_2 - f_1) / (2 * stepwidth)
            det_deriv = (det_2 - det_1) / (2 * stepwidth)

            if (det_vals(k) < 0) then
                if (f_vals(k) < 0) then
                    g_derivs(k) = f_deriv * det_vals(k) + det_deriv * f_vals(k)
                else
                    g_derivs(k) = - (f_deriv * det_vals(k) + det_deriv * f_vals(k))
                end if
            else
                if (f_vals(k) < 0) then
                    g_derivs(k) = - (f_deriv * det_vals(k) + det_deriv * f_vals(k))
                else
                    g_derivs(k) = f_deriv * det_vals(k) + det_deriv * f_vals(k)
                end if
            end if
        end do

    end subroutine compute_g_derivs_a_num

    subroutine compute_g_derivs_a_inv_num(g_derivs, x_coords_array, w_coords_array, params, i, j, stepwidth)
        real(p), dimension(:), intent(inout) :: g_derivs
        real(p), dimension(:,:), intent(in)  :: w_coords_array, x_coords_array
        type(mapping), intent(in)            :: params
        integer, intent(in)                  :: i, j
        real(p), intent(in)                  :: stepwidth
        type(mapping)                        :: params_1, params_2
        real(p), allocatable, dimension(:,:) :: x_coords_1, x_coords_2
        real(p)                              :: f_1, f_2, det_1, det_2
        integer                              :: k
        real(p)                              :: f_deriv, det_deriv

        allocate(x_coords_1(1:size(w_coords_array(:,1)),1:params%dim), x_coords_2(1:size(w_coords_array(:,1)),1:params%dim))

        call init_mapping_empty(params_1, params%dim)
        call init_mapping_empty(params_2, params%dim)
        params_1 = params
        params_2 = params
        params_1%a_inv(i, j) = params_1%a_inv(i, j) - stepwidth
        params_2%a_inv(i, j) = params_2%a_inv(i, j) + stepwidth

        x_coords_1 = params_1%backtransform_array(w_coords_array)
        x_coords_2 = params_2%backtransform_array(w_coords_array)

        do k = 1, size(w_coords_array(:,1))
            f_1   = f_ptr(x_coords_1(k,:))
            f_2   = f_ptr(x_coords_2(k,:))
            det_1 = determinant( params_1%jacobian_backtrafo(x_coords_1(k,:), w_coords_array(k,:)), params%dim)
            det_2 = determinant( params_2%jacobian_backtrafo(x_coords_2(k,:), w_coords_array(k,:)), params%dim)

            f_deriv   = (f_2 - f_1) / (2 * stepwidth)
            det_deriv = (det_2 - det_1) / (2 * stepwidth)

            if (det_vals(k) < 0) then
                if (f_vals(k) < 0) then
                    g_derivs(k) = f_deriv * det_vals(k) + det_deriv * f_vals(k)
                else
                    g_derivs(k) = - (f_deriv * det_vals(k) + det_deriv * f_vals(k))
                end if
            else
                if (f_vals(k) < 0) then
                    g_derivs(k) = - (f_deriv * det_vals(k) + det_deriv * f_vals(k))
                else
                    g_derivs(k) = f_deriv * det_vals(k) + det_deriv * f_vals(k)
                end if
            end if
        end do

    end subroutine compute_g_derivs_a_inv_num

    subroutine compute_g_derivs_b_num(g_derivs, x_coords_array, w_coords_array, params, i, stepwidth)
        real(p), dimension(:), intent(inout) :: g_derivs
        real(p), dimension(:,:), intent(in)  :: w_coords_array, x_coords_array
        type(mapping), intent(in)            :: params
        integer, intent(in)                  :: i
        real(p), intent(in)                  :: stepwidth
        type(mapping)                        :: params_1, params_2
        real(p), allocatable, dimension(:,:) :: x_coords_1, x_coords_2
        real(p)                              :: f_1, f_2, det_1, det_2
        integer                              :: k
        real(p)                              :: f_deriv, det_deriv

        allocate(x_coords_1(1:size(w_coords_array(:,1)),1:params%dim), x_coords_2(1:size(w_coords_array(:,1)),1:params%dim))

        call init_mapping_empty(params_1, params%dim)
        call init_mapping_empty(params_2, params%dim)
        params_1 = params
        params_2 = params
        params_1%b(i) = params_1%b(i) - stepwidth
        params_2%b(i) = params_2%b(i) + stepwidth

        x_coords_1 = params_1%backtransform_array(w_coords_array)
        x_coords_2 = params_2%backtransform_array(w_coords_array)

        do k = 1, size(w_coords_array(:,1))
            f_1   = f_ptr(x_coords_1(k,:))
            f_2   = f_ptr(x_coords_2(k,:))
            det_1 = determinant( params_1%jacobian_backtrafo(x_coords_1(k,:), w_coords_array(k,:)), params%dim)
            det_2 = determinant( params_2%jacobian_backtrafo(x_coords_2(k,:), w_coords_array(k,:)), params%dim)

            f_deriv   = (f_2 - f_1) / (2 * stepwidth)
            det_deriv = (det_2 - det_1) / (2 * stepwidth)

            if (det_vals(k) < 0) then
                if (f_vals(k) < 0) then
                    g_derivs(k) = f_deriv * det_vals(k) + det_deriv * f_vals(k)
                else
                    g_derivs(k) = - (f_deriv * det_vals(k) + det_deriv * f_vals(k))
                end if
            else
                if (f_vals(k) < 0) then
                    g_derivs(k) = - (f_deriv * det_vals(k) + det_deriv * f_vals(k))
                else
                    g_derivs(k) = f_deriv * det_vals(k) + det_deriv * f_vals(k)
                end if
            end if
        end do

    end subroutine compute_g_derivs_b_num

!!!________________________Prozeduren für die Berechnung des Optimizers und dessen Gradienten_________________________

    function opt_val()
        real(p)                            :: opt_val
        real(p), allocatable, dimension(:) :: h_vals
        integer                            :: i

        allocate(h_vals(1:size(g_vals)))
        call h_vals_ptr(h_vals)
        opt_val = 0.0
        do i = 1, size(g_vals)
            opt_val = opt_val + h_vals(i)
        end do
        opt_val = opt_val / size(g_vals)

    end function opt_val


    function opt_deriv_a(x_coords_array, w_coords_array, params, i, j, stepwidth)
        real(p)                             :: opt_deriv_a
        real(p), dimension(:,:), intent(in) :: x_coords_array, w_coords_array
        type(mapping), intent(in)           :: params
        integer, intent(in)                 :: i, j
        real(p), intent(in)                 :: stepwidth
        real(p), allocatable, dimension(:)  :: g_derivs, h_derivs
        integer                             :: k

        allocate(g_derivs(size(w_coords_array(:,1))), h_derivs(size(w_coords_array(:,1))))
        call compute_g_derivs_a(g_derivs, x_coords_array, w_coords_array, params, i, j, stepwidth)
        call h_derivs_ptr(h_derivs, g_derivs)

        opt_deriv_a = 0.0
        do k = 1, size(g_derivs)
            opt_deriv_a = opt_deriv_a + h_derivs(k)
        end do
        opt_deriv_a = opt_deriv_a / size(h_derivs)

    end function opt_deriv_a

    function opt_deriv_a_inv(x_coords_array, w_coords_array, params, i, j, stepwidth)
        real(p)                             :: opt_deriv_a_inv
        real(p), dimension(:,:), intent(in) :: x_coords_array, w_coords_array
        type(mapping), intent(in)           :: params
        integer, intent(in)                 :: i, j
        real(p), intent(in)                 :: stepwidth
        real(p), allocatable, dimension(:)  :: g_derivs, h_derivs
        integer                             :: k

        allocate(g_derivs(size(w_coords_array(:,1))), h_derivs(size(w_coords_array(:,1))))
        call compute_g_derivs_a_inv(g_derivs, x_coords_array, w_coords_array, params, i, j, stepwidth)
        call h_derivs_ptr(h_derivs, g_derivs)

        opt_deriv_a_inv = 0.0
        do k = 1, size(g_derivs)
            opt_deriv_a_inv = opt_deriv_a_inv + h_derivs(k)
        end do
        opt_deriv_a_inv = opt_deriv_a_inv / size(h_derivs)

    end function opt_deriv_a_inv

    function opt_deriv_b(x_coords_array, w_coords_array, params, i, stepwidth)
        real(p)                             :: opt_deriv_b
        real(p), dimension(:,:), intent(in) :: x_coords_array, w_coords_array
        type(mapping), intent(in)           :: params
        integer, intent(in)                 :: i
        real(p), intent(in)                 :: stepwidth
        real(p), allocatable, dimension(:)  :: g_derivs, h_derivs
        integer                             :: k

        allocate(g_derivs(size(w_coords_array(:,1))), h_derivs(size(w_coords_array(:,1))))
        call compute_g_derivs_b(g_derivs, x_coords_array, w_coords_array, params, i, stepwidth)
        call h_derivs_ptr(h_derivs, g_derivs)

        opt_deriv_b = 0.0
        do k = 1, size(g_derivs)
            opt_deriv_b = opt_deriv_b + h_derivs(k)
        end do
        opt_deriv_b = opt_deriv_b / size(h_derivs)

    end function opt_deriv_b

    function opt_deriv_a_num(x_coords_array, w_coords_array, params, i, j, stepwidth)
        real(p)                             :: opt_deriv_a_num
        real(p), dimension(:,:), intent(in) :: x_coords_array, w_coords_array
        type(mapping), intent(in)           :: params
        integer, intent(in)                 :: i, j
        real(p), intent(in)                 :: stepwidth
        real(p), allocatable, dimension(:)  :: g_derivs, h_derivs
        integer                             :: k

        allocate(g_derivs(size(w_coords_array(:,1))), h_derivs(size(w_coords_array(:,1))))
        call compute_g_derivs_a_num(g_derivs, x_coords_array, w_coords_array, params, i, j, stepwidth)
        call h_derivs_ptr(h_derivs, g_derivs)

        opt_deriv_a_num = 0.0
        do k = 1, size(g_derivs)
            opt_deriv_a_num = opt_deriv_a_num + h_derivs(k)
        end do
        opt_deriv_a_num = opt_deriv_a_num / size(h_derivs)

    end function opt_deriv_a_num

    function opt_deriv_a_inv_num(x_coords_array, w_coords_array, params, i, j, stepwidth)
        real(p)                             :: opt_deriv_a_inv_num
        real(p), dimension(:,:), intent(in) :: x_coords_array, w_coords_array
        type(mapping), intent(in)           :: params
        integer, intent(in)                 :: i, j
        real(p), intent(in)                 :: stepwidth
        real(p), allocatable, dimension(:)  :: g_derivs, h_derivs
        integer                             :: k

        allocate(g_derivs(size(w_coords_array(:,1))), h_derivs(size(w_coords_array(:,1))))
        call compute_g_derivs_a_inv_num(g_derivs, x_coords_array, w_coords_array, params, i, j, stepwidth)
        call h_derivs_ptr(h_derivs, g_derivs)

        opt_deriv_a_inv_num = 0.0
        do k = 1, size(g_derivs)
            opt_deriv_a_inv_num = opt_deriv_a_inv_num + h_derivs(k)
        end do
        opt_deriv_a_inv_num = opt_deriv_a_inv_num / size(h_derivs)

    end function opt_deriv_a_inv_num

    function opt_deriv_b_num(x_coords_array, w_coords_array, params, i, stepwidth)
        real(p)                             :: opt_deriv_b_num
        real(p), dimension(:,:), intent(in) :: x_coords_array, w_coords_array
        type(mapping), intent(in)           :: params
        integer, intent(in)                 :: i
        real(p), intent(in)                 :: stepwidth
        real(p), allocatable, dimension(:)  :: g_derivs, h_derivs
        integer                             :: k

        allocate(g_derivs(size(w_coords_array(:,1))), h_derivs(size(w_coords_array(:,1))))
        call compute_g_derivs_b_num(g_derivs, x_coords_array, w_coords_array, params, i, stepwidth)
        call h_derivs_ptr(h_derivs, g_derivs)

        opt_deriv_b_num = 0.0
        do k = 1, size(g_derivs)
            opt_deriv_b_num = opt_deriv_b_num + h_derivs(k)
        end do
        opt_deriv_b_num = opt_deriv_b_num / size(h_derivs)

    end function opt_deriv_b_num

!!!________________The same again, using online algorithms__________________

    ! function opt_val()
    !     real(p)                            :: opt_val
    !     real(p), allocatable, dimension(:) :: h_vals
    !     integer                            :: i
    !     real(p)                             :: mult

    !     allocate(h_vals(1:size(g_vals)))
    !     call h_vals_ptr(h_vals)
    !     opt_val = 0.0
    !     do i = 1, size(g_vals)
    !         mult = 1 / real(i, p)
    !         opt_val = opt_val + mult * (h_vals(i) - opt_val)
    !     end do

    ! end function opt_val


    ! function opt_deriv_a(x_coords_array, w_coords_array, params, i, j, stepwidth)
    !     real(p)                             :: opt_deriv_a
    !     real(p), dimension(:,:), intent(in) :: x_coords_array, w_coords_array
    !     type(mapping), intent(in)           :: params
    !     integer, intent(in)                 :: i, j
    !     real(p), intent(in)                 :: stepwidth
    !     real(p), allocatable, dimension(:)  :: g_derivs, h_derivs
    !     integer                             :: k
    !     real(p)                             :: mult

    !     allocate(g_derivs(size(w_coords_array(:,1))), h_derivs(size(w_coords_array(:,1))))
    !     call compute_g_derivs_a(g_derivs, x_coords_array, w_coords_array, params, i, j, stepwidth)
    !     call h_derivs_ptr(h_derivs, g_derivs)

    !     opt_deriv_a = 0.0
    !     do k = 1, size(g_derivs)
    !         mult = 1 / real(i, p)
    !         opt_deriv_a = opt_deriv_a + mult * (h_derivs(k) - opt_deriv_a)
    !     end do

    ! end function opt_deriv_a

    ! function opt_deriv_a_inv(x_coords_array, w_coords_array, params, i, j, stepwidth)
    !     real(p)                             :: opt_deriv_a_inv
    !     real(p), dimension(:,:), intent(in) :: x_coords_array, w_coords_array
    !     type(mapping), intent(in)           :: params
    !     integer, intent(in)                 :: i, j
    !     real(p), intent(in)                 :: stepwidth
    !     real(p), allocatable, dimension(:)  :: g_derivs, h_derivs
    !     integer                             :: k
    !     real(p)                             :: mult

    !     allocate(g_derivs(size(w_coords_array(:,1))), h_derivs(size(w_coords_array(:,1))))
    !     call compute_g_derivs_a_inv(g_derivs, x_coords_array, w_coords_array, params, i, j, stepwidth)
    !     call h_derivs_ptr(h_derivs, g_derivs)

    !     opt_deriv_a_inv = 0.0
    !     do k = 1, size(g_derivs)
    !         mult = 1 / real(i, p)
    !         opt_deriv_a_inv = opt_deriv_a_inv + mult * (h_derivs(k) - opt_deriv_a_inv)
    !     end do

    ! end function opt_deriv_a_inv

    ! function opt_deriv_b(x_coords_array, w_coords_array, params, i, stepwidth)
    !     real(p)                             :: opt_deriv_b
    !     real(p), dimension(:,:), intent(in) :: x_coords_array, w_coords_array
    !     type(mapping), intent(in)           :: params
    !     integer, intent(in)                 :: i
    !     real(p), intent(in)                 :: stepwidth
    !     real(p), allocatable, dimension(:)  :: g_derivs, h_derivs
    !     integer                             :: k
    !     real(p)                             :: mult

    !     allocate(g_derivs(size(w_coords_array(:,1))), h_derivs(size(w_coords_array(:,1))))
    !     call compute_g_derivs_b(g_derivs, x_coords_array, w_coords_array, params, i, stepwidth)
    !     call h_derivs_ptr(h_derivs, g_derivs)

    !     opt_deriv_b = 0.0
    !     do k = 1, size(g_derivs)
    !         mult = 1 / real(i, p)
    !         opt_deriv_b = opt_deriv_b + mult * (h_derivs(k) - opt_deriv_b)
    !     end do

    ! end function opt_deriv_b

    ! function opt_deriv_a_num(x_coords_array, w_coords_array, params, i, j, stepwidth)
    !     real(p)                             :: opt_deriv_a_num
    !     real(p), dimension(:,:), intent(in) :: x_coords_array, w_coords_array
    !     type(mapping), intent(in)           :: params
    !     integer, intent(in)                 :: i, j
    !     real(p), intent(in)                 :: stepwidth
    !     real(p), allocatable, dimension(:)  :: g_derivs, h_derivs
    !     integer                             :: k
    !     real(p)                             :: mult

    !     allocate(g_derivs(size(w_coords_array(:,1))), h_derivs(size(w_coords_array(:,1))))
    !     call compute_g_derivs_a_num(g_derivs, x_coords_array, w_coords_array, params, i, j, stepwidth)
    !     call h_derivs_ptr(h_derivs, g_derivs)

    !     opt_deriv_a_num = 0.0
    !     do k = 1, size(g_derivs)
    !         mult = 1 / real(i, p)
    !         opt_deriv_a_num = opt_deriv_a_num + mult * (h_derivs(k) - opt_deriv_a_num)
    !     end do

    ! end function opt_deriv_a_num

    ! function opt_deriv_a_inv_num(x_coords_array, w_coords_array, params, i, j, stepwidth)
    !     real(p)                             :: opt_deriv_a_inv_num
    !     real(p), dimension(:,:), intent(in) :: x_coords_array, w_coords_array
    !     type(mapping), intent(in)           :: params
    !     integer, intent(in)                 :: i, j
    !     real(p), intent(in)                 :: stepwidth
    !     real(p), allocatable, dimension(:)  :: g_derivs, h_derivs
    !     integer                             :: k
    !     real(p)                             :: mult

    !     allocate(g_derivs(size(w_coords_array(:,1))), h_derivs(size(w_coords_array(:,1))))
    !     call compute_g_derivs_a_inv_num(g_derivs, x_coords_array, w_coords_array, params, i, j, stepwidth)
    !     call h_derivs_ptr(h_derivs, g_derivs)

    !     opt_deriv_a_inv_num = 0.0
    !     do k = 1, size(g_derivs)
    !         mult = 1 / real(i, p)
    !         opt_deriv_a_inv_num = opt_deriv_a_inv_num + mult * (h_derivs(k) - opt_deriv_a_inv_num)
    !     end do

    ! end function opt_deriv_a_inv_num

    ! function opt_deriv_b_num(x_coords_array, w_coords_array, params, i, stepwidth)
    !     real(p)                             :: opt_deriv_b_num
    !     real(p), dimension(:,:), intent(in) :: x_coords_array, w_coords_array
    !     type(mapping), intent(in)           :: params
    !     integer, intent(in)                 :: i
    !     real(p), intent(in)                 :: stepwidth
    !     real(p), allocatable, dimension(:)  :: g_derivs, h_derivs
    !     integer                             :: k
    !     real(p)                             :: mult

    !     allocate(g_derivs(size(w_coords_array(:,1))), h_derivs(size(w_coords_array(:,1))))
    !     call compute_g_derivs_b_num(g_derivs, x_coords_array, w_coords_array, params, i, stepwidth)
    !     call h_derivs_ptr(h_derivs, g_derivs)

    !     opt_deriv_b_num = 0.0
    !     do k = 1, size(g_derivs)
    !         mult = 1 / real(i, p)
    !         opt_deriv_b_num = opt_deriv_b_num + mult * (h_derivs(k) - opt_deriv_b_num)
    !     end do

    ! end function opt_deriv_b_num

end module module_optimize
