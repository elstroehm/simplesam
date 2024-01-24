module module_function

    use module_matrix
    use module_transform

    implicit none

    abstract interface
        function f_value(x_coords)
            integer, parameter :: p = selected_real_kind(12, 30)
            real(p)                           :: f_value
            real(p), dimension(:), intent(in) :: x_coords
        end function f_value
    end interface

    !!!Der Funktionszeiger enhält eine Referenz auf diejenige Funktion, welche für den Algorithmus
    !!!verwendet werden soll (Phase-Space-Sampling eines Lorentz-Peaks, einer abfallenden Wand etc.).
    !!!Damit wird vermieden, dass für unterschiedliche Funktionen jeweils eigene Prozeduren zum Ableiten
    !!!implementiert werden müssen.
    !!!
    !!!center, height und width enthalten die Funktionsargumente, sodass diese nicht jedes mal
    !!!explizit an den Zeiger übergeben werden müssen.

    procedure(f_value), pointer        :: f_ptr
    real(p), allocatable, dimension(:) :: center
    real(p)                            :: height

    !!!Der Parameter height ist nicht genau gleich den Parametern für a und b in der Bachelorarbeit.
    !!!Er gibt die maximale Höhe des Peaks f_1/der Wand f_2 an und somit height = 1 / a für f_1 bzw.
    !!!b = 1 / height + 1 für f_2. center gibt das Zentrum des Peaks an.

contains

!!!_______________Initialisierung der Parameter des modules________________

    subroutine init_module_function(f_ptr_, center_, height_)
        procedure(f_value), pointer, intent(in) :: f_ptr_
        real(p), dimension(:), intent(in)       :: center_
        real(p), intent(in)                     :: height_

        allocate(center(1:size(center_)))
        f_ptr  => f_ptr_
        center = center_
        height = height_

    end subroutine init_module_function

!!!______________Prozeduren für die Berechnung von f(x)_______________

    !!!Lorentz peak f_1
    function f_value_peak(x_coords)
        real(p)                           :: f_value_peak
        real(p), dimension(:), intent(in) :: x_coords
        integer                           :: i
        real(p)                           :: dummy

        dummy = 1 / height
        do i = 1, size(x_coords)
            dummy = dummy + (x_coords(i) - center(i))**2
        end do

        f_value_peak = 1 / dummy

    end function f_value_peak

    !!!Funktion f_2
    function f_value_wall(x_coords)
        real(p)                           :: f_value_wall
        real(p), dimension(:), intent(in) :: x_coords
        real(p)                           :: dummy

        dummy = 1 / height + (x_coords(1) + 1.0)
        f_value_wall = 1 / dummy

    end function f_value_wall
    
    !!!Funktion g_3
    function f_value_cylinder(x_coords)
        real(p)                           :: f_value_cylinder
        real(p), dimension(:), intent(in) :: x_coords
        integer                           :: i
        real(p)                           :: dummy

        dummy = 1 / height
        do i = 1, 2
            dummy = dummy + (x_coords(i) - center(i))**2
        end do

        f_value_cylinder = 1 / dummy

    end function f_value_cylinder

!!!________________Prozeduren für die Berechnung der Ableitung von f(x) nach einem a_n_________________

!!!Die folgenden Prozeduren berechnen die Ableitung von f nach Transformationsparametern auf numerische Weise
!!!Die Benennung der Prozeduren erfolgt nach dem gleichen Schema wie in module_transform

    function f_derivative_a(params, w_coords, i, j, stepwidth)
        real(p)                            :: f_derivative_a
        real(p), dimension(:), intent(in)  :: w_coords
        type(mapping), intent(in)          :: params
        integer, intent(in)                :: i, j
        real(p), intent(in)                :: stepwidth
        type(mapping)                      :: params_1, params_2
        real(p), allocatable, dimension(:) :: x_coords_1, x_coords_2
        real(p)                            :: f_1, f_2

        allocate(x_coords_1(1:params%dim), x_coords_2(1:params%dim))
        call init_mapping_empty(params_1, params%dim)
        call init_mapping_empty(params_2, params%dim)
        params_1 = params
        params_2 = params
        params_1%a(i, j) = params_1%a(i, j) - stepwidth
        params_2%a(i, j) = params_2%a(i, j) + stepwidth
        params_1%a_inv = params_1%a
        params_2%a_inv = params_2%a
        call invert_matrix(params_1%a_inv, params%dim)
        call invert_matrix(params_1%a_inv, params%dim)

        x_coords_1 = params_1%backtransform(w_coords)
        x_coords_2 = params_2%backtransform(w_coords)

        f_1 = f_ptr(x_coords_1)
        f_2 = f_ptr(x_coords_2)
        f_derivative_a = (f_2 - f_1) / (2 * stepwidth)

    end function f_derivative_a

    function f_derivative_a_inv(params, w_coords, i, j, stepwidth)
        real(p)                            :: f_derivative_a_inv
        real(p), dimension(:), intent(in)  :: w_coords
        type(mapping), intent(in)          :: params
        integer, intent(in)                :: i, j
        real(p), intent(in)                :: stepwidth
        type(mapping)                      :: params_1, params_2
        real(p), allocatable, dimension(:) :: x_coords_1, x_coords_2
        real(p)                            :: f_1, f_2

        allocate(x_coords_1(1:params%dim), x_coords_2(1:params%dim))
        call init_mapping_empty(params_1, params%dim)
        call init_mapping_empty(params_2, params%dim)
        params_1 = params
        params_2 = params
        params_1%a_inv(i, j) = params_1%a_inv(i, j) - stepwidth
        params_2%a_inv(i, j) = params_2%a_inv(i, j) + stepwidth

        x_coords_1 = params_1%backtransform(w_coords)
        x_coords_2 = params_2%backtransform(w_coords)

        f_1 = f_ptr(x_coords_1)
        f_2 = f_ptr(x_coords_2)
        f_derivative_a_inv = (f_2 - f_1) / (2 * stepwidth)

    end function f_derivative_a_inv

    function f_derivative_b(params, w_coords, i, stepwidth)
        real(p)                            :: f_derivative_b
        real(p), dimension(:), intent(in)  :: w_coords
        type(mapping), intent(in)          :: params
        integer, intent(in)                :: i
        real(p), intent(in)                :: stepwidth
        type(mapping)                      :: params_1, params_2
        real(p), allocatable, dimension(:) :: x_coords_1, x_coords_2
        real(p)                            :: f_1, f_2

        allocate(x_coords_1(1:params%dim), x_coords_2(1:params%dim))
        call init_mapping_empty(params_1, params%dim)
        call init_mapping_empty(params_2, params%dim)
        params_1 = params
        params_2 = params
        params_1%b(i) = params_1%b(i) - stepwidth
        params_2%b(i) = params_2%b(i) + stepwidth

        x_coords_1 = params_1%backtransform(w_coords)
        x_coords_2 = params_2%backtransform(w_coords)

        f_1 = f_ptr(x_coords_1)
        f_2 = f_ptr(x_coords_2)

        f_derivative_b = (f_2 - f_1) / (2 * stepwidth)

    end function f_derivative_b

    function f_derivative_a_array(params, w_coords_array, i, j, stepwidth)
        real(p), allocatable, dimension(:)   :: f_derivative_a_array
        real(p), dimension(:,:), intent(in)  :: w_coords_array
        type(mapping), intent(in)            :: params
        integer, intent(in)                  :: i, j
        real(p), intent(in)                  :: stepwidth
        type(mapping)                        :: params_1, params_2
        real(p), allocatable, dimension(:,:) :: x_coords_1, x_coords_2
        real(p)                              :: f_1, f_2
        integer                              :: k

        allocate(x_coords_1(1:size(w_coords_array(:,1)),1:params%dim), x_coords_2(1:size(w_coords_array(:,1)),1:params%dim))
        allocate(f_derivative_a_array(1:size(w_coords_array(:,1))))
        call init_mapping_empty(params_1, params%dim)
        call init_mapping_empty(params_2, params%dim)
        params_1 = params
        params_2 = params
        params_1%a(i, j) = params_1%a(i, j) - stepwidth
        params_2%a(i, j) = params_2%a(i, j) + stepwidth
        params_1%a_inv = params_1%a
        params_2%a_inv = params_2%a
        call invert_matrix(params_1%a_inv, params%dim)
        call invert_matrix(params_1%a_inv, params%dim)

        x_coords_1 = params_1%backtransform_array(w_coords_array)
        x_coords_2 = params_2%backtransform_array(w_coords_array)

        do k = 1, size(w_coords_array(:,1))
            f_1 = f_ptr(x_coords_1(k,:))
            f_2 = f_ptr(x_coords_2(k,:))
            f_derivative_a_array(k) = (f_2 - f_1) / (2 * stepwidth)
        end do

    end function f_derivative_a_array

    function f_derivative_a_inv_array(params, w_coords_array, i, j, stepwidth)
        real(p), allocatable, dimension(:)   :: f_derivative_a_inv_array
        real(p), dimension(:,:), intent(in)  :: w_coords_array
        type(mapping), intent(in)            :: params
        integer, intent(in)                  :: i, j
        real(p), intent(in)                  :: stepwidth
        type(mapping)                        :: params_1, params_2
        real(p), allocatable, dimension(:,:) :: x_coords_1, x_coords_2
        real(p)                              :: f_1, f_2
        integer                              :: k

        allocate(x_coords_1(1:size(w_coords_array(:,1)),1:params%dim), x_coords_2(1:size(w_coords_array(:,1)),1:params%dim))
        allocate(f_derivative_a_inv_array(1:size(w_coords_array(:,1))))
        call init_mapping_empty(params_1, params%dim)
        call init_mapping_empty(params_2, params%dim)
        params_1 = params
        params_2 = params
        params_1%a_inv(i, j) = params_1%a_inv(i, j) - stepwidth
        params_2%a_inv(i, j) = params_2%a_inv(i, j) + stepwidth

        x_coords_1 = params_1%backtransform_array(w_coords_array)
        x_coords_2 = params_2%backtransform_array(w_coords_array)

        do k = 1, size(w_coords_array(:,1))
            f_1 = f_ptr(x_coords_1(k,:))
            f_2 = f_ptr(x_coords_2(k,:))
            f_derivative_a_inv_array(k) = (f_2 - f_1) / (2 * stepwidth)
        end do

    end function f_derivative_a_inv_array

    function f_derivative_b_array(params, w_coords_array, i, stepwidth)
        real(p), allocatable, dimension(:)   :: f_derivative_b_array
        real(p), dimension(:,:), intent(in)  :: w_coords_array
        type(mapping), intent(in)            :: params
        integer, intent(in)                  :: i
        real(p), intent(in)                  :: stepwidth
        type(mapping)                        :: params_1, params_2
        real(p), allocatable, dimension(:,:) :: x_coords_1, x_coords_2
        real(p)                              :: f_1, f_2
        integer                              :: k

        allocate(x_coords_1(1:size(w_coords_array(:,1)),1:params%dim), x_coords_2(1:size(w_coords_array(:,1)),1:params%dim))
        allocate(f_derivative_b_array(1:size(w_coords_array(:,1))))
        call init_mapping_empty(params_1, params%dim)
        call init_mapping_empty(params_2, params%dim)
        params_1 = params
        params_2 = params
        params_1%b(i) = params_1%b(i) - stepwidth
        params_2%b(i) = params_2%b(i) + stepwidth

        x_coords_1 = params_1%backtransform_array(w_coords_array)
        x_coords_2 = params_2%backtransform_array(w_coords_array)

        do k = 1, size(w_coords_array(:,1))
            f_1 = f_ptr(x_coords_1(k,:))
            f_2 = f_ptr(x_coords_2(k,:))
            f_derivative_b_array(k) = (f_2 - f_1) / (2 * stepwidth)
        end do

    end function f_derivative_b_array

end module module_function
