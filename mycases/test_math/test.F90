program test

  use ModConf
  use ModBasicMath
  use MPI

  implicit none

  real(WP) :: A(4, 5), B(5, 4)
  integer :: i, j

  call InitMPI
  A = reshape( [1, 0, 0, 0, 0, 0, 0, 4, 0, 3, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0], shape(A))

  print *, '4 x 5 array A'
  do i = 1, 4
    print *, A(i, :)
  end do

  call Matrix_PseudoInvert(A, B)

  print *, '5 x 4 array B'
  do i = 1, 5
    print *, B(i, :)
  end do

  call FinalizeMPI
  stop


end program test