module module_matrix

!!!contains the basic framework for matrix computations
!!!maybe promote to a class "matrix" for object oriented coding
!!!change the arguments from pointers to normal variables or arrays

!!!compute determinant of a matrix
!!!compute lu-decomposition of a matrix
!!!compute inverse of a matrix
!!!swap rows and columns of a matrix

    implicit none

    integer, parameter :: p = selected_real_kind(12, 30)

contains

    !!!performs an lu-decomposition of "matrix" and overwrites it with the result
    !!!"perm_matrix" denotes the permutation matrix used to rearrange the rows of "matrix", here commented out (will not be used)
    !!!"perm_sign" yields "-1" for an odd number of row permutations and "1" for an even number of row permutations
    !!!"dim" is the dimension of the matrix
    !!!returns zero matrix if the algorithm fails
    subroutine lu_decompose(matrix, dim, perm_sign)
        real(p), dimension(:,:), intent(inout) :: matrix
        integer, intent(in)                    :: dim
        integer, intent(out)                   :: perm_sign
        real(p)                                :: pivot, dummy
        integer                                :: pivot_row, i, j, k

        !!!set permutation_matrix and perm_ij to identity and "perm_sign" to "1"
        ! do i = 1, dim
        !     do j = 1, dim
        !         if (i == j) then
        !             perm_matrix(i, j) = 1
        !             perm_ij(i, j)     = 1
        !             cycle
        !         end if
        !         perm_matrix(i, j) = 0
        !         perm_ij(i, j)     = 0
        !     end do
        ! end do
        perm_sign = 1

        !!!Crout's algorithm for lu-decomposition, loop through j and i and distinguish between i<j, i=j and i>j
        !!!and calculate the corresponding matrix entries, overwriting the entries in "matrix"
        do j = 1, dim
            do i = 1, j-1
                do k = 1, i-1
                    matrix(i, j) = matrix(i, j) - (matrix(i, k) * matrix(k, j))
                end do
            end do
            pivot = 0.0
            do i = j, dim
                do k = 1, j-1   !!!no division by matrix(j, j) yet, could be zero -> row permutation might be necessary
                    matrix(i, j) = matrix(i, j) - (matrix(i, k) * matrix(k, j))
                end do
                if (abs(matrix(i, j)) > pivot) then !!!find pivot element for division
                    pivot = matrix(i, j)
                    pivot_row = i
                end if
            end do
            if (pivot == 0) then    !!!pivot element vanishes -> Crout's algorithm fails -> terminate and return zero matrix
                matrix = 0
                return
            else if (j /= pivot_row) then   !!!keep track of the row permutations, change perm_matrix and perm_sign accordingly
                do k = 1, dim
                    dummy                = matrix(j, k)
                    matrix(j, k)         = matrix(pivot_row, k)
                    matrix(pivot_row, k) = dummy
                end do
                ! perm_ij(j, j)                   = 0
                ! perm_ij(pivot_row, pivot_row)   = 0
                ! perm_ij(j, pivot_row)           = 1
                ! perm_ij(pivot_row, j)           = 1
                ! perm_matrix = matmul(perm_ij, perm_matrix)
                perm_sign = -perm_sign
            end if
            pivot = 1/pivot
            do i = j+1, dim   !!!divide the lower diagonal elements in the column by the pivot element
                matrix(i, j) = pivot * matrix(i, j)
            end do
        end do

    end subroutine lu_decompose

    !!!computes the determinant of matrix using the subroutine lu_decompose
    function determinant(matrix, dim)
        real(p), dimension(:,:), intent(in)  :: matrix
        integer, intent(in)                  :: dim
        real(p)                              :: determinant
        real(p), allocatable, dimension(:,:) :: copymatrix
        integer                              :: i, perm_sign

        allocate(copymatrix(1:dim, 1:dim))
        copymatrix = matrix
        call lu_decompose(copymatrix, dim, perm_sign)
        determinant = perm_sign

        do i = 1, dim
            determinant = determinant * copymatrix(i, i)
        end do

    end function determinant

    !!!computes the inverse "matrix"
    !!!"dim" is the dimension of "matrix"
    !!!returns zero matrix if "matrix" is not invertible
    !!!matrix gets overwritten in the process
    !!!with real(p) entries -> precision up to at least four or five digits
    subroutine invert_matrix(matrix, dim)
        real(p), dimension(:,:), intent(inout) :: matrix
        integer, intent(in)                    :: dim
        real(p), allocatable, dimension(:,:)   :: inverted
        real(p)                                :: dummy_1, dummy_2
        integer                                :: i, j, k

        allocate(inverted(1:dim, 1:dim))

        !!!set inverted to identity
        do i = 1, dim
            do j = 1, dim
                if (i == j) then
                    inverted(i, j) = 1
                    cycle
                end if
                inverted(i, j) = 0
            end do
        end do

        do i = 1, dim - 1
            if (matrix(i,i) == 0) then  !!!pivot element vanishes -> search for pivot elements in the other rows
                do j = i, dim
                    if (matrix(j, i) == 0) then     !!!new pivot element not found yet
                        cycle
                    end if
                    do k = 1, dim       !!!pivot element found -> swap rows
                        dummy_1          = matrix(i, k)
                        matrix(i, k)     = matrix(j, k)
                        matrix(j, k)     = dummy_1
                        dummy_1          = inverted(i, k)
                        inverted(i, k)   = inverted(j, k)
                        inverted(j, k)   = dummy_1
                    end do
                    exit    !!!since the pivot element has been found -> exit the loop
                end do
                if(matrix(i, i) == 0) then  !!!search for pivot element not successful -> no pivot element present -> singular matrix
                    matrix = 0
                    return
                end if
            end if
            dummy_1 = 1 / matrix(i, i)
            do j = i + 1, dim       !!!gauss-elimination on matrix and inverted -> upper triangle form of matrix
                dummy_2 = matrix(j, i) * dummy_1
                do k = 1, i-1
                    inverted(j, k)     = inverted(j, k) - dummy_2 * inverted(i, k)
                end do
                do k = i, dim
                    matrix(j, k)       = matrix(j, k) - dummy_2 * matrix(i, k)
                    inverted(j, k)     = inverted(j, k) - dummy_2 * inverted(i, k)
                end do
            end do
        end do
        do i = dim, 1, -1       !!!bring matrix to diagonal form
            if (matrix(i, i) == 0) then     !!!can't divide by zero -> singular matrix
                matrix = 0
                return
            end if
            dummy_1 = 1 / matrix(i, i)
            do j = i-1, 1, -1
                dummy_2 = dummy_1 * matrix(j, i)
                do k = 1, i-1
                    inverted(j, k) = inverted(j, k) - (dummy_2 * inverted(i, k))
                end do
                do k = i, dim
                    matrix(j, k)   = matrix(j, k) - (dummy_2 * matrix(i, k))
                    inverted(j, k) = inverted(j, k) - (dummy_2 * inverted(i, k))
                end do
            end do
        end do
        do i = 1, dim       !!!divide inverted by diagonal elements of matrix
            if(matrix(i, i) == 0) then
                matrix = 0
                return
            end if
            dummy_1 = 1 / matrix(i, i)
            do j = 1, dim
                inverted(i, j) = inverted(i, j) * dummy_1
            end do
        end do
        matrix = inverted

    end subroutine invert_matrix

end module module_matrix
