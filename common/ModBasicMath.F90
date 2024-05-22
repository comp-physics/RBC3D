! Collection of basic math operation
module ModBasicMath

  use ModDataTypes
  use f95_lapack, only: LA_POSV
  use MPI
  use ModConf
!  use mkl95_lapack, only : LA_POSV=>POSV

  implicit none

  private

  public :: VecNorm, &
            CrossProd, &
            InvMat2, &
            InvMat3, &
            RotateMatrix, &
            TriArea, &
            Matrix_PseudoInvert, &
            QuadFit_1D, &
            QuadFit_2D, &
            Min_Quad_2D, &
            MaskFunc_Exact, &
            MaskFunc, &
            BSplineFunc, &
            RandomNumber

contains

!**********************************************************************
! L2 norm of a vector
! |a|
  function VecNorm(a) result(c)
    real(WP) :: a(:), c

    ! print *, "a", a

    c = sqrt(sum(a*a))

  end function VecNorm

!**********************************************************************
! Cross product of two vectors
! a x b
  function CrossProd(a, b) result(c)
    real(WP) :: a(3), b(3), c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  end function CrossProd

!**********************************************************************
! Inverse of a 2x2 matrix
  function InvMat2(a) result(b)
    real(WP) :: a(2, 2), b(2, 2)

    real(WP) :: idetA

    idetA = 1./(a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1))
    b(1, 1) = idetA*a(2, 2)
    b(1, 2) = -idetA*a(1, 2)
    b(2, 1) = -idetA*a(2, 1)
    b(2, 2) = idetA*a(1, 1)

  end function InvMat2

!**********************************************************************
! Inverse of a 3 by 3 matrix
  function InvMat3(a) result(b)
    real(WP) :: a(3, 3), b(3, 3)

    real(WP) :: detA

    detA = a(1, 1)*a(2, 2)*a(3, 3) + a(1, 2)*a(2, 3)*a(3, 1) + a(1, 3)*a(2, 1)*a(3, 2) &
           - a(1, 1)*a(2, 3)*a(3, 2) - a(1, 2)*a(2, 1)*a(3, 3) - a(1, 3)*a(2, 2)*a(3, 1)

    b(1, 1) = a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)
    b(1, 2) = -(a(1, 2)*a(3, 3) - a(1, 3)*a(3, 2))
    b(1, 3) = a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2)

    b(2, 1) = -(a(2, 1)*a(3, 3) - a(2, 3)*a(3, 1))
    b(2, 2) = a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1)
    b(2, 3) = -(a(1, 1)*a(2, 3) - a(1, 3)*a(2, 1))

    b(3, 1) = a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1)
    b(3, 2) = -(a(1, 1)*a(3, 2) - a(1, 2)*a(3, 1))
    b(3, 3) = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)

    b = (1./detA)*b

  end function InvMat3

!**********************************************************************
! Compute the rotation matrix
! Arguments:
!  c1 -- the new z-axis
! Note:
!  1. c1 must have unit norm
!  2. The transform matrix rotates the z-axis (0, 0, 1) to c1
  function RotateMatrix(c1) result(mat)
    real(WP) :: c1(3), mat(3, 3)

    ! a -- the axis of rotation
    ! a = c x c1, where c is the z-axis
    ! b = c x a
    ! b1 = c1 x a
    ! mat = a * a + b1 * b + c1 * c
    real(WP) :: a(3), b(3), c(3), b1(3)
    integer :: i, j
    real(WP), parameter :: eps = 1.E-10

    c = (/0._WP, 0._WP, 1._WP/)
    a = (/-c1(2), c1(1), 0._WP/)

    if (sum(a*a) < eps) then
      mat = 0.
      mat(1, 1) = 1.
      mat(2, 2) = 1.
      mat(3, 3) = 1.
      return
    end if

    a = a/sqrt(sum(a**2))
    b = (/-a(2), a(1), 0._WP/)
    b1 = CrossProd(c1, a)

    forall (i=1:3, j=1:3) mat(i, j) = a(i)*a(j) + b1(i)*b(j) + c1(i)*c(j)

  end function RotateMatrix

!**********************************************************************
! The area of a triangle
! Argument:
!  x(i,:) -- i-th corner coordinate
  function TriArea(x) result(area)
    real(WP) :: x(3, 3), area

    real(WP) :: a(3), b(3), c(3)

    a = x(2, :) - x(1, :)
    b = x(3, :) - x(1, :)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

    area = 0.5*sqrt(sum(c*c))

  end function TriArea

!**********************************************************************
! Compute the psedo-inverse of a matrix
! Arguments:
!  A -- original matrix
!  B -- A^{-1}
  subroutine Matrix_PseudoInvert(A, B)
    real(WP) :: A(:, :), B(:, :)

    integer :: M, N
    real(WP), allocatable :: AtA(:, :)

    ! Allocate working arrays
    M = size(A, 1)
    N = size(A, 2)
    allocate (AtA(N, N))

    B = transpose(A)
    AtA = matmul(A, A)
    call LA_POSV(AtA, B)

    ! Deallocate working arrays
    deallocate (AtA)

  end subroutine Matrix_PseudoInvert

!**********************************************************************
! 1D Quadratic fit
! Arguments:
!  x, f -- coordinates and function values
!  a0, a1, a2 -- parameters fitted
! Note:
!   f(x) = a0 + a1*x + a2*(x**2)
  subroutine QuadFit_1D(x, f, a0, a1, a2)
    real(WP) :: x(:), f(:), a0, a1, a2

    real(WP) :: lhs(3, 3), rhs(3)
    real(WP) :: xi, xi2, xi3, xi4
    integer :: i

    lhs = 0.
    rhs = 0.
    do i = 1, size(x)
      xi = x(i)
      xi2 = xi*xi
      xi3 = xi*xi2
      xi4 = xi*xi3

      lhs(1, 1) = lhs(1, 1) + 1.
      lhs(1, 2) = lhs(1, 2) + xi
      lhs(1, 3) = lhs(1, 3) + xi2

      lhs(2, 2) = lhs(2, 2) + xi2
      lhs(2, 3) = lhs(2, 3) + xi3

      lhs(3, 3) = lhs(3, 3) + xi4

      rhs = rhs + (/1._WP, xi, xi2/)*f(i)
    end do ! i

    lhs(2, 1) = lhs(1, 2)
    lhs(3, 1) = lhs(1, 3)
    lhs(3, 2) = lhs(2, 3)

    call LA_POSV(lhs, rhs)

    a0 = rhs(1)
    a1 = rhs(2)
    a2 = rhs(3)

  end subroutine QuadFit_1D

!**********************************************************************
! 2D Quadratic fit
! Arguments:
!  x(i,:), f(i) -- coordinates and functional value of the i-th point
! Note:
!   f(x) = a0 + a1*x + a2*y + a11*(x**2) + a12*(x*y) + a22*(y**2)
  subroutine QuadFit_2D(x, f, a0, a1, a2, a11, a12, a22)
    real(WP) :: x(:, :), f(:), a0, a1, a2, a11, a12, a22

    real(WP) :: lhs(6, 6), rhs(6), u(6)
    real(WP) :: xi, yi
    integer :: i, ii, jj
    integer :: ierr

    real(WP) :: clockBgn, clockEnd

    clockBgn = MPI_WTime()

    ! Compute the upper half of lhs and rhs
    lhs = 0.
    rhs = 0.

    do i = 1, size(x, 1)
      xi = x(i, 1)
      yi = x(i, 2)

      u(1) = 1.
      u(2) = xi
      u(3) = yi
      u(4) = xi*xi
      u(5) = xi*yi
      u(6) = yi*yi

      do ii = 1, 6
        do jj = ii, 6
          lhs(ii, jj) = lhs(ii, jj) + u(ii)*u(jj)
        end do ! jj
      end do ! ii

      rhs = rhs + u(:)*f(i)
    end do ! i

    ! Compute the lower half of the symmetric lhs
    do ii = 2, 6
      do jj = 1, ii - 1
        lhs(ii, jj) = lhs(jj, ii)
      end do ! jj
    end do ! ii
    ! lapack call
    call LA_POSV(lhs, rhs, INFO=ierr)

    a0 = rhs(1)
    a1 = rhs(2)
    a2 = rhs(3)
    a11 = rhs(4)
    a12 = rhs(5)
    a22 = rhs(6)

    clockEnd = MPI_WTime()
    if (rootWorld) then
      write(*, '(A, I3, A, F12.2)') &
        "time cost = ", clockEnd - clockBgn
    end if

  end subroutine QuadFit_2D

!**********************************************************************
! Find the minimum of a quadratic function
! Arguments:
!  a0 to a22 -- coefficients of the qudratic function
!  xmin -- coordinate of the minimum point
!  fmin -- minimum function value
!
! Note:
!  -- f = a0 + a1*x + a2*y + a11*x^2 + a12*x*y + a22*y^2
!  -- It is up to the user to find whether fmin is really the minimum,
!     since xmin could also be a maximum or saddle point
  subroutine Min_Quad_2D(a0, a1, a2, a11, a12, a22, xmin, fmin)
    real(WP) :: a0, a1, a2, a11, a12, a22
    real(WP) :: xmin(2)
    real(WP), optional :: fmin

    real(WP) :: lhs(2, 2), rhs(2), det, idet

    lhs(1, 1) = 2.*a11
    lhs(1, 2) = a12
    lhs(2, 1) = a12
    lhs(2, 2) = 2.*a22
    det = lhs(1, 1)*lhs(2, 2) - lhs(1, 2)*lhs(2, 1)

    rhs = -(/a1, a2/)

    if (det > 0) then
      ! Exclude saddle point
      idet = 1./det
      xmin(1) = idet*(lhs(2, 2)*rhs(1) - lhs(1, 2)*rhs(2))
      xmin(2) = idet*(-lhs(2, 1)*rhs(1) + lhs(1, 1)*rhs(2))
    else
      xmin = 0.
    end if

    if (present(fmin)) then
      fmin = a0 + a1*xmin(1) + a2*xmin(2) &
             + a11*xmin(1)**2 + a12*xmin(1)*xmin(2) + a22*xmin(2)**2
    end if

  end subroutine Min_Quad_2D

!**********************************************************************
! An infinitely smooth mask function
!  f(0) = 1
!  f(1) = 0
  function MaskFunc_Exact(x) result(f)
    real(WP) :: x, f

    real(WP) :: t

    t = abs(x)

    if (t < 0.01) then
      f = 1.
    else if (t > 0.99) then
      f = 0.
    else
      f = exp(2*(exp(-1./t))/(t - 1))
    end if

  end function MaskFunc_Exact

!**********************************************************************
! Mask function by linear interpolation
  function MaskFunc(x) result(f)
    real(WP) :: x, f

    integer, parameter :: N = 8192
    real(WP) :: ftab(0:N)
    logical, save :: table_inited = .false.
    real(WP) :: s
    integer :: i

    ! Build look up table
    if (.not. table_inited) then
      do i = 0, N
        s = real(i)/N
        ftab(i) = MaskFunc_Exact(s)
      end do ! i
      table_inited = .true.
    end if

    ! Linear interpolation
    s = abs(x)*N
    i = floor(s)

    if (i >= N) then
      f = 0.
    else
      f = ftab(i)*(i + 1 - s) + ftab(i + 1)*(s - i)
    end if

  end function MaskFunc

!**********************************************************************
! Evaluate the B-spline function on a uniform mesh
! Argument:
!  xc -- position of the delta-function singularity
!  P -- number of B-spline points
!  imin -- start point where the B-spline function value is non-zero
!  w(:) -- B-spline function values, it should be mapped to the physical
!          axis as w(1) -> w(imin), w(2) -> w(imin+1), and so on
! Note:
!  -- The order of B-spline function is (P-1)
!  -- The dimension of w must be equal or bigger than P
  subroutine BsplineFunc(xc, P, imin, w)
    real(WP) :: xc
    integer :: P
    integer :: imin
    real(WP) :: w(:)

    real(WP) :: u(1:P)
    integer :: pp, j

    imin = floor(xc) - (P - 1)
    u(1) = imin - (xc - P)
    do j = 2, P
      u(j) = u(j - 1) + 1
    end do ! i

    w(1) = 1.
    if (P >= 2) w(2:P) = 0.

    do pp = 2, P
      do j = pp, 2, -1
        w(j) = u(j)/(pp - 1.)*w(j) + (pp - u(j))/(pp - 1.)*w(j - 1)
      end do ! j
      w(1) = u(1)/(pp - 1.)*w(1)
    end do ! pp

  end subroutine BsplineFunc

!**********************************************************************
! Random number generator, copied from Numerical Recipies
  function RandomNumber(idum)
    integer, intent(IN), optional :: idum
    real(WP)        :: RandomNumber

    integer, parameter :: K4B = selected_int_kind(9)
    integer(K4B) :: idum1
    integer(K4B), parameter :: IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836
    real(WP), save :: am
    integer(K4B), save :: ix = -1, iy = -1, k
    !$omp threadprivate(am, ix, iy, k)

    ! Initialize
    if (iy < 0) then
      if (present(idum)) then
        idum1 = idum
      else
        write (*, *) 'No random number seed provided, use default value 1'
        idum1 = 1
      end if

      am = nearest(1.0, -1.0)/IM
      iy = ior(ieor(888889999_K4B, abs(idum1)), 1_K4B)
      ix = ieor(777755555_K4B, abs(idum1))
    end if

    ix = ieor(ix, ishft(ix, 13))
    ix = ieor(ix, ishft(ix, -17))
    ix = ieor(ix, ishft(ix, 5))
    k = iy/IQ
    iy = IA*(iy - k*IQ) - IR*k
    if (iy < 0) iy = iy + IM

    RandomNumber = am*ior(iand(IM, ieor(ix, iy)), 1_K4B)

  end function RandomNumber

!**********************************************************************

end module ModBasicMath
