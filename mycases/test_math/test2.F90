program test

  use ModDataTypes
  use ModConf
  use f95_lapack, only: LA_POSV
  use MPI

  implicit none

  integer :: i, j
  real(WP) :: A(5, 5)
  real(WP) :: B(5, 3)
  A = reshape( [38, 3, 4, 6, 5, 3, 48, 6, 7, 7, 4, 6, 52, 1, 6, 6, 7, 1, 56, 10, 5, 7, 6, 10, 74], shape(A))
  B = reshape( [56, 71, 69, 80, 102, 112, 142, 138, 160, 204, 168, 213, 207, 240, 306], shape(B))

  call InitMPI
  print *, 'A before:'
  do i = 1, 5
    do j= 1, 5
      write (*, "(f8.3)", advance="no") A(i,j)
    end do
    print *, ""
  end do

  print *, 'B before:'
  do i = 1, 5
    do j= 1, 3
      write (*, "(f8.3)", advance="no") B(i,j)
    end do
    print *, ""
  end do

  print *, ""
  print *, "calling LA_POSV(A, B)"
  print *, ""

  ! call DPOTRS(A, B)
  call LA_POSV(A, B)

  print *, 'A after:'
  do i = 1, 5
    do j= 1, 5
      write (*, "(f8.3)", advance="no") A(i,j)
    end do
    print *, ""
  end do

  print *, 'B after:'
  do i = 1, 5
    do j= 1, 3
      write (*, "(f8.3)", advance="no") B(i,j)
    end do
    print *, ""
  end do



  call FinalizeMPI
  stop


end program test