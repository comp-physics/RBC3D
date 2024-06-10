program test

  use ModData
  use ModDataTypes
  use ModDataStruct
  use ModBasicMath

  implicit none

  integer :: i
  real(WP) :: c(3), c1(3), a(3), res(3), matrix(3, 3)

  call InitMPI

  c = (/0, 0, 1/)
  c1 = (/1, 0, 0/)
  a = (/0, 1, 0/)

  matrix = RotateMatrix(c1)

  print *, "rotation matrix: "
  do i = 1, 3
    print *, matrix(i, :)
  end do

  call FinalizeMPI
end program test

! program test

!     use ModData
!     use ModDataTypes
!     use ModDataStruct
!     use ModBasicMath

!     implicit none

!     call InitMPI

!     integer :: n_points = 4, i
!     real(WP) :: A(4, 2)
!     real(WP) :: b(4)
!     real(WP) :: x0 = 0., x1 = 0., x2 = 0., x11 = 0., x12 = 0., x22 = 0., timeit
!     real(WP) :: x, y

!     real(WP) :: c(3), c1(3), res(3)

!     c = (/0, 0, 1/)
!     c1 = (/1, 0, 0/)

!     res = CrossProd(c, c1)

!     call FinalizeMPI
! end program test