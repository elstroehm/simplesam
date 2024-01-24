module module_transform

    !!!Enthält die Klasse mapping, welche die Transformation Phi darstellt
    !!!Methoden zum Transformieren von X->W und W->X, Berechnung der Jacobimatrix und
    !!!Berechnung Ableitung der Jacobideterminante nach Transformationsparametern a_n
    !!!bei gegebenen Koordinaten w_coords

    use module_matrix

    implicit none

    type mapping
        real(p), allocatable, dimension(:,:) :: a       !Matrix a mit Einträgen a_ij
        real(p), allocatable, dimension(:,:) :: a_inv   !Inverse der Matrix a
        real(p), allocatable, dimension(:)   :: b       !Translationsvektor b
        integer                              :: dim     !Dimension n

    contains
        procedure :: transform
        procedure :: transform_array
        procedure :: backtransform
        procedure :: backtransform_array

        procedure :: jacobian_trafo
        procedure :: jacobian_backtrafo

        procedure :: det_derivative_a
        procedure :: det_derivative_a_inv
        procedure :: det_derivative_b
        procedure :: det_derivative_a_array
        procedure :: det_derivative_a_inv_array
        procedure :: det_derivative_b_array

    end type mapping

contains

!!!__________________________Konstruktoren__________________________

    !!!Konstruktor für die Klasse mapping. Erzeugt eine Abbildung mit den Parametern a (Matrix), b (Translation) und dim
    subroutine init_mapping_explicit(params, a, b, dim)
        integer, intent(in)                 :: dim
        real(p), dimension(:,:), intent(in) :: a
        real(p), dimension(:), intent(in)   :: b
        type(mapping), intent(inout)        :: params

        allocate(params%a(1:dim, 1:dim), params%a_inv(1:dim, 1:dim), &
                 & params%b(1:dim))
                 params%a = a
                 params%a_inv = a
        call invert_matrix(params%a_inv, dim)
        params%b = b
        params%dim = dim

    end subroutine init_mapping_explicit

    !!!Konstruktor für die Klasse mapping. Alloziert lediglich den nötigen Speicher für a und b,
    !!!ohne diese mit Werten zu belegen. Definiert dim
    subroutine init_mapping_empty(params, dim)
        type(mapping), intent(out) :: params
        integer, intent(in)        :: dim

        allocate(params%a(1:dim, 1:dim), params%a_inv(1:dim, 1:dim), &
                 & params%b(1:dim))
        params%dim = dim

    end subroutine init_mapping_empty

    !!!Konstruktor für die Klasse mapping. Erzeugt die identische Abbildung
    subroutine init_mapping_id(params, dim)
        type(mapping), intent(out) :: params
        integer, intent(in)        :: dim
        integer                    :: i, j

        allocate(params%a(1:dim, 1:dim), params%a_inv(1:dim, 1:dim), &
                 & params%b(1:dim))
        params%dim = dim

        do i = 1, dim
            do j = 1, dim
                if (i == j) then
                    params%a(i, j) = 1.0
                    params%a_inv(i, j) = 1.0
                    cycle
                end if
                params%a(i, j) = 0.0
                params%a_inv(i, j) = 0.0
            end do
        end do
        params%b = 0

    end subroutine init_mapping_id

    !!!Kein Konstruktor. Setzt params auf die identische Abbildung
    subroutine init_id(params)
        type(mapping), intent(inout) :: params
        integer                      :: i, j

        do i = 1, params%dim
            do j = 1, params%dim
                if (i == j) then
                    params%a(i, j) = 1.0
                    params%a_inv(i, j) = 1.0
                    cycle
                end if
                params%a(i, j) = 0.0
                params%a_inv(i, j) = 0.0
            end do
        end do
        params%b = 0

    end subroutine init_id

!!!___________________________________Methoden zur Koordinatentransformation______________________________________

!!!Zur Benennung der Prozeduren: ohne "array" bedeutet, dass nut ein einzelner Koordinatenvektor als Argument
!!!übergeben wird, mit "array" bedeutet, dass ein ganzes Tupel an Vektoren übergeben wird.

    !!!Transformation X->W
    function transform(this, x_coords)
        real(p), allocatable, dimension(:) :: transform
        class(mapping), intent(in)         :: this
        real(p), dimension(:), intent(in)  :: x_coords

        allocate(transform(this%dim))
        transform = tanh(matmul(this%a, atanh(x_coords)) + this%b)

    end function transform

    !!!Rücktransformation W->x
    function transform_array(this, x_coords_array)
        real(p), allocatable, dimension(:,:) :: transform_array
        class(mapping), intent(in)           :: this
        real(p), dimension(:,:), intent(in)  :: x_coords_array
        integer                              :: i

        allocate(transform_array(1:size(x_coords_array(:,1)), 1:this%dim))
        transform_array = matmul(atanh(x_coords_array), transpose(this%a))
        do i = 1, size(transform_array(:,1))
            transform_array(i,:) = tanh(transform_array(i,:) + this%b)
        end do

    end function transform_array

    !!!Transformation X->W eines arrays von Koordinaten
    !!!Zeilen = Vektoren, Spalten = Koordinaten der Vektoren
    function backtransform(this, w_coords)
        real(p), allocatable, dimension(:) :: backtransform
        class(mapping), intent(in)         :: this
        real(p), dimension(:), intent(in)  :: w_coords

        allocate(backtransform(1:this%dim))
        backtransform = tanh(matmul(this%a_inv, atanh(w_coords) - this%b))

    end function backtransform

    !!!Rücktransformation W->X eines arrays von Koordinaten
    !!!Zeilen = Vektoren, Spalten = Koordinaten der Vektoren
    function backtransform_array(this, w_coords_array)
        real(p), allocatable, dimension(:,:) :: backtransform_array
        class(mapping), intent(in)           :: this
        real(p), dimension(:,:), intent(in)  :: w_coords_array
        integer                              :: i

        allocate(backtransform_array(1:size(w_coords_array(:,1)), 1:this%dim))
        do i = 1, size(backtransform_array(:,1))
            backtransform_array(i,:) = atanh(w_coords_array(i,:)) - this%b
        end do
        backtransform_array = tanh(matmul(backtransform_array, transpose(this%a_inv)))

    end function backtransform_array

!!!_________________Methoden zur Berechnung der Jacobimatrix________________

    !!!Jacobideterminante der Transformation am Punkt x_coords in X
    !!!w_coords = Koordinaten von x_coords im W-Raum
    function jacobian_trafo(this, x_coords, w_coords)
        real(p), dimension(:), intent(in)    :: x_coords, w_coords
        class(mapping), intent(in)           :: this
        real(p), allocatable, dimension(:,:) :: jacobian_trafo
        integer                              :: i, j
        real(p), allocatable, dimension(:)   :: x_squared, w_squared

        allocate(jacobian_trafo(1:this%dim, 1:this%dim), &
                & x_squared(1:this%dim), w_squared(1:this%dim))
        x_squared = x_coords * x_coords
        w_squared = w_coords * w_coords

        do i = 1, this%dim
            do j = 1, this%dim
                jacobian_trafo(i, j) = this%a(i, j) * (1 - w_squared(i)) / (1 - x_squared(j))
            end do
        end do

    end function jacobian_trafo

    !!!Jacobideterminante der Rücktransformation am Punkt w_coords in W
    !!!x_coords = Koordinaten von w_coords im X-Raum
    function jacobian_backtrafo(this, x_coords, w_coords)
        real(p), dimension(:), intent(in)    :: x_coords, w_coords
        class(mapping), intent(in)           :: this
        real(p), allocatable, dimension(:,:) :: jacobian_backtrafo
        integer                              :: i, j
        real(p), allocatable, dimension(:)   :: x_squared, w_squared

        allocate(jacobian_backtrafo(1:this%dim, 1:this%dim), &
                & x_squared(1:this%dim), w_squared(1:this%dim))
        x_squared = x_coords * x_coords
        w_squared = w_coords * w_coords

        do i = 1, this%dim
            do j = 1, this%dim
                jacobian_backtrafo(i, j) = this%a_inv(i, j) * (1 - x_squared(i)) / (1 - w_squared(j))
            end do
        end do

    end function jacobian_backtrafo

!!!_________________Methoden zur Berechnung der Ableitung der Jacobideterminante_________________

!!!Die folgenden Prozeduren berechnen die Ableitung von det(J(w)) nach den Transformationsparametern mithilfe
!!!der in Appendix A gegebenen analytischen Ausdrücke.
!!!Zur Benennung der Prozeduren: a, a_inv, b bedeutet jeweils Ableitung nach Matrixparameter a_ij, a^{-1}_ij
!!!oder Translationsparameter b_i

    !!!Berechnet die Ableitung von det(J(w)) an der Stelle w_coords bezüglich eines Matrixparameters a_ij
    function det_derivative_a(this, x_coords, w_coords, det_val, i, j)
        real(p)                           :: det_derivative_a
        class(mapping), intent(in)        :: this
        real(p), dimension(:), intent(in) :: x_coords, w_coords
        real(p), intent(in)               :: det_val
        integer, intent(in)               :: i, j

        det_derivative_a = - det_val * (this%a_inv(j, i) - 2 * w_coords(i) * atanh(x_coords(j)))

    end function det_derivative_a

    !!!Berechnet die Ableitung von det(J(w)) an der Stelle w_coords bezüglich eines Matrixparameters a^{-1}_ij
    function det_derivative_a_inv(this, x_coords, w_coords, det_val, i, j)
        real(p)                           :: det_derivative_a_inv
        class(mapping), intent(in)        :: this
        real(p), dimension(:), intent(in) :: x_coords, w_coords
        real(p), intent(in)               :: det_val
        integer, intent(in)               :: i, j

        det_derivative_a_inv = det_val * (this%a(j, i) - &
        & 2 * x_coords(i) * (atanh(w_coords(j)) - this%b(j)) )

    end function det_derivative_a_inv

    !!!Berechnet die Ableitung von det(J(w)) an der Stelle w_coords bezüglich eines Translationsparameters b_ij
    function det_derivative_b(this, x_coords, det_val, i)
        real(p)                           :: det_derivative_b
        class(mapping), intent(in)        :: this
        real(p), dimension(:), intent(in) :: x_coords
        real(p), intent(in)               :: det_val
        integer, intent(in)               :: i
        integer                           :: l

        det_derivative_b = 0.0
        do l = 1, this%dim
            det_derivative_b = det_derivative_b + (x_coords(l) * this%a_inv(l, i))
        end do
        det_derivative_b = 2 * det_val * det_derivative_b

    end function det_derivative_b

    !!!Berechnet die Ableitung von det(J(w)) an allen Stellen in w_coords_array bezüglich eines Matrixparameters a^{-1}_ij
    function det_derivative_a_array(this, x_coords_array, w_coords_array, det_vals, i, j)
        real(p), allocatable, dimension(:)  :: det_derivative_a_array
        class(mapping), intent(in)          :: this
        real(p), dimension(:,:), intent(in) :: x_coords_array, w_coords_array
        real(p), dimension(:), intent(in)   :: det_vals
        integer, intent(in)                 :: i, j
        integer                             :: k

        allocate(det_derivative_a_array(1:size(w_coords_array(:,1))))

        do k = 1, size(det_vals)
            det_derivative_a_array(k) = - det_vals(k) * &
            & (this%a_inv(j, i) - 2 * w_coords_array(k, i) * atanh(x_coords_array(k, j)))
        end do

    end function det_derivative_a_array

    !!!Berechnet die Ableitung von det(J(w)) an allen Stellen in w_coords_array bezüglich eines Matrixparameters a_ij
    function det_derivative_a_inv_array(this, x_coords_array, w_coords_array, det_vals, i, j)
        real(p), allocatable, dimension(:)  :: det_derivative_a_inv_array
        class(mapping), intent(in)          :: this
        real(p), dimension(:,:), intent(in) :: x_coords_array, w_coords_array
        real(p), dimension(:), intent(in)   :: det_vals
        integer, intent(in)                 :: i, j
        integer                             :: k

        allocate(det_derivative_a_inv_array(1:size(w_coords_array(:,1))))

        do k = 1, size(det_vals)
            det_derivative_a_inv_array(k) = det_vals(k) * (this%a(j, i) - &
            & 2 * x_coords_array(k, i) * (atanh(w_coords_array(k, j)) - this%b(j)) )
        end do

    end function det_derivative_a_inv_array

    !!!Berechnet die Ableitung von det(J(w)) an allen Stellen in w_coords_array bezüglich eines Translationsparameters b_i
    function det_derivative_b_array(this, x_coords_array, det_vals, i)
        real(p), allocatable, dimension(:)  :: det_derivative_b_array
        class(mapping), intent(in)          :: this
        real(p), dimension(:,:), intent(in) :: x_coords_array
        real(p), dimension(:), intent(in)   :: det_vals
        integer, intent(in)                 :: i
        integer                             :: k, l
        real(p)                             :: dummy

        allocate(det_derivative_b_array(1:size(x_coords_array(:,1))))

        do k = 1, size(det_vals)
            dummy = 0.0
            do l = 1, this%dim
                dummy = dummy + ((x_coords_array(k, l) * this%a_inv(l, i)))
            end do
            det_derivative_b_array(k) = 2 * det_vals(k) * dummy
        end do

    end function det_derivative_b_array

end module module_transform
