program solve_linear_system
    implicit none
    integer, parameter :: n = 3
    integer :: info, i, j
    double precision :: A(n, n), b(n)
    integer :: ipiv(n)
    double precision, allocatable :: A_copy(:,:)
  
    ! Define the matrix A in column-major order
    A = reshape([3.0d0,  2.0d0, -1.0d0, &
                -1.0d0,  1.0d0,  2.0d0, &
                 2.0d0,  3.0d0,  1.0d0], &
                shape(A))
  
    ! Define the right-hand side vector b
    b = [1.0d0, 2.0d0, 3.0d0]
  
    ! Allocate a copy of A to keep the original matrix for verification
    allocate(A_copy(n, n))
    A_copy = A
  
    ! Print the original matrix A
    print *, 'Original Matrix A:'
    do i = 1, n
      print *, (A(i, j), j = 1, n)
    end do
  
    ! Print the original vector b
    print *, 'Original Vector b:'
    print *, b
  
    ! Call LAPACK routine DGESV to solve the linear system
    call dgesv(n, 1, A, n, ipiv, b, n, info)
  
    ! Check for successful execution
    if (info /= 0) then
      print *, 'Error: DGESV did not converge. Info = ', info
    else
      ! Print the solution vector x
      print *, 'Solution vector x:'
      print *, b
  
      ! Compute A_copy * x to verify the solution
      print *, 'Reconstructed b (A_copy * x):'
      do i = 1, n
        b(i) = 0.0d0
        do j = 1, n
          b(i) = b(i) + A_copy(i, j) * b(j)
        end do
      end do
  
      ! Print the reconstructed b
      print *, b
  
      ! Print the original vector b again for comparison
      print *, 'Original Vector b for comparison:'
      print *, [1.0d0, 2.0d0, 3.0d0]
    end if
  
end program solve_linear_system
! program solve_linear_system

!     use ModData
    
!     implicit none
!     integer, parameter :: n = 3
!     integer :: info, i, j
!     real(WP) :: A(n, n), b(n)
!     integer :: ipiv(n)
  
!     ! Define the matrix A in column-major order
!     A = reshape([3.0, -1.0, 2.0, &
!                2.0,  1.0, 3.0, &
!               -1.0,  2.0, 1.0], &
!               shape(A))

!     A = reshape((/1., 0., 2., 0.,  1., 1., &
!                 2.,  3., 0./), shape(A))
  
!     ! Define the right-hand side vector b
!     b = [7., 5., 8.]
  
!     ! Print the matrix A
!     print *, 'Matrix A:'
!     do i = 1, n
!       print *, (A(i, j), j = 1, n)
!     end do
  
!     ! Print the vector b
!     print *, 'Vector b:'
!     print *, b
  
!     ! Call LAPACK routine DGESV to solve the linear system
!     call dgesv(n, 1, A, n, ipiv, b, n, info)
  
!     ! Check for successful execution
!     if (info /= 0) then
!       print *, 'Error: DGESV did not converge. Info = ', info
!     else
!       ! Print the solution vector x
!       print *, 'Solution vector x:'
!       print *, b

!       do i = 1, n
!         print *, (A(i, j), j = 1, n)
!       end do
!     end if
  
! end program solve_linear_system