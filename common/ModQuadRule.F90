! Quadrature for numerical integration
module ModQuadRule

  use ModDataTypes

  implicit none

  ! 1D Gauss-Legendre quadrature
  ! n -- number of quadrature points
  ! rst(i) -- the i-th abscissa
  ! w(i) -- the i-th weight
  type t_GaussQuad1D
    integer :: n
    real(WP), dimension(:), pointer :: rst
    real(WP), dimension(:), pointer :: w(:)
  end type t_GaussQuad1D

  ! 2D Gauss-Legendre quadrature
  ! n -- number of quadrature points
  ! rst(i,:) -- the i-th abscissa
  ! w(i) -- the i-th weight
  type t_GaussQuad2D
    integer :: n
    real(WP), dimension(:, :), pointer :: rst
    real(WP), dimension(:), pointer :: w(:)
  end type t_GaussQuad2D

  type(t_GaussQuad2D), target :: gqTri7, gqTri3
  type(t_GaussQuad2D), target :: gqQuad9

  public

  public :: GaussQuad_Init, &
            GaussQuad_Finalize, &
            GauLeg, &
            GauLeg_Sinh

  private :: MyAsinh

contains

!**********************************************************************
! Purpose:
!  Set up abscissa and weights of Gauss quadrature rules
  subroutine GaussQuad_Init

    integer :: n
    real(WP) :: r, s, t
    real(WP), allocatable :: xtmp(:), wtmp(:)
    integer :: i, j, m

    !======================================================================
    ! Create 3-point quadrature for triangle
    n = 3
    gqTri3%n = n
    allocate (gqTri3%rst(n, 2))
    allocate (gqTri3%w(n))

    gqTri3%w = 0.5/3.

    r = 2./3.
    s = 1./6.
    t = 1./6.

    gqTri3%rst(:, 1) = (/r, s, t/)
    gqTri3%rst(:, 2) = (/s, t, r/)

    !======================================================================
    ! Create 7-point quadrature for triangle
    n = 7
    gqTri7%n = n
    allocate (gqTri7%rst(n, 2))
    allocate (gqTri7%w(n))

    r = (6.-sqrt(15.))/21.
    s = r
    t = 1 - r - s
    gqTri7%rst(1:3, 1) = (/r, s, t/)
    gqTri7%rst(1:3, 2) = (/s, t, r/)
    gqTri7%w(1:3) = (155.-sqrt(15.))/2400.

    r = (6.+sqrt(15.))/21.
    s = r
    t = 1 - r - s
    gqTri7%rst(4:6, 1) = (/r, s, t/)
    gqTri7%rst(4:6, 2) = (/s, t, r/)
    gqTri7%w(4:6) = (155.+sqrt(15.))/2400.

    r = 1./3.
    s = 1./3.
    t = 1 - r - s
    gqTri7%rst(7, 1) = r
    gqTri7%rst(7, 2) = s
    gqTri7%w(7) = 9./80.

    !======================================================================
    ! Create 9-point quadratuer for quadrilaterals
    n = 9
    gqQuad9%n = 9
    allocate (gqQuad9%rst(n, 2), gqQuad9%w(n))

    allocate (xtmp(3), wtmp(3))
    call GauLeg(-1._WP, 1._WP, 3, xtmp, wtmp)

    do i = 1, 3
    do j = 1, 3
      m = i + (j - 1)*3
      gqQuad9%rst(m, 1:2) = (/xtmp(i), xtmp(j)/)
      gqQuad9%w(m) = wtmp(i)*wtmp(j)
    end do ! j
    end do ! i

    deallocate (xtmp, wtmp)

  end subroutine GaussQuad_Init

!**********************************************************************
  subroutine GaussQuad_Finalize

    deallocate (gqTri3%rst, gqTri3%w)
    deallocate (gqTri7%rst, gqTri7%w)
    deallocate (gqQuad9%rst, gqQuad9%w)

  end subroutine GaussQuad_Finalize

!**********************************************************************
! Given the lower and upper limits of integration x1 and x2, this routine
! returns arrays x and w of length N containing the abscissas and weights
! of the Gauss-Legendre N-point quadrature formula. The parameter EPS is
! the relative precision. Note that internal computations are done in
! double precision.
!
! Note:
!  Copied from "Numerical Recipes in Fortran 90"
  subroutine GauLeg(x1, x2, n, x, w)
    real(WP) :: x1, x2
    integer :: n
    real(WP) :: x(:), w(:)

    character(*), parameter :: func_name = 'GauLeg'
    double precision, parameter :: EPS = 3.D-14
    integer :: its, j, m
    integer, parameter :: MAXIT = 10
    double precision :: xl, xm
    double precision, dimension((n + 1)/2) :: p1, p2, p3, pp, z, z1
    logical, dimension((n + 1)/2) :: unfinished

    ! The roots are symmetric in the interval, so we
    ! only have to find half of them
    m = (n + 1)/2
    xm = 0.5*(x2 + x1)
    xl = 0.5*(x2 - x1)

    ! Initial approximations to the roots.
    z(1) = 1
    do j = 2, size(z)
      z(j) = z(j - 1) + 1
    end do ! j

    z = cos(PI*(z - 0.25)/(n + 0.5))
    unfinished = .true.

    ! Newton’s method carried out simultaneously on the roots
    do its = 1, MAXIT
      where (unfinished)
        p1 = 1.0
        p2 = 0.0
      end where

      ! Loop up the recurrence relation to get the Legendre
      ! polynomials evaluated at z
      do j = 1, n
        where (unfinished)
          p3 = p2
          p2 = p1
          p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j
        end where
      end do

      ! p1 now contains the desired Legendre polynomials.
      ! We next compute pp, the derivatives, by a standard
      ! relation involving also p2, the polynomials of one
      ! lower order.
      where (unfinished)
        pp = n*(z*p1 - p2)/(z*z - 1.0)
        z1 = z
        z = z1 - p1/pp          ! Newton’s method.
        unfinished = (abs(z - z1) > EPS)
      end where

      if (.not. any(unfinished)) exit
    end do

    if (its == MAXIT + 1) then
      write (*, *) 'Subroutine ', func_name
      write (*, *) 'Error: too many iterations'
      stop
    end if

    ! Scale the root to the desired interval, and
    ! put in its symmetric counterpart.
    x(1:m) = xm - xl*z
    x(n:n - m + 1:-1) = xm + xl*z
    ! Compute the weight and its symmetric counterpart.
    w(1:m) = 2.0*xl/((1.0 - z**2)*pp**2)
    w(n:n - m + 1:-1) = w(1:m)

  end subroutine GauLeg

!**********************************************************************
! Create a Gauss-Legendre N-point quadrature formula by a sinh transformation
! Arguments:
!  xmin, xmax -- integration interval
!  (a, b) -- the target point coordinate.  b is very small (but non-zero),
!            which will cause a near singularity for the integration around
!            (a,0)
!  n, x, w -- same as those in subroutine GauLeg
!
! References:
!   K. Road and H. Tasmania, "A sinh transformation for evaluating nearly
!   singular boundary element integrals", Int. J. Numer. Engng 2005;
!   62:564-578.
  subroutine GauLeg_Sinh(xmin, xmax, a, b, n, x, w)
    real(WP) :: xmin, xmax, a, b
    integer :: n
    real(WP) :: x(:), w(:)

    real(WP) :: xm, xl, a0, b0
    real(WP) :: u1, u2, mu, eta
    real(WP) :: s(n)
    integer :: i

    ! Scale the integration interval
    !  x = xm + xl*y, with y in [-1, 1]
    xm = 0.5*(xmin + xmax)
    xl = 0.5*(xmax - xmin)

    a0 = (a - xm)/xl
    b0 = b*xl

    ! Do sinh transformation such that
    !  y = a + b*sinh(mu*s - eta), with s in [-1,1]
    u1 = MyAsinh((1 + a0)/b0)
    u2 = MyAsinh((1 - a0)/b0)
    mu = 0.5*(u1 + u2)
    eta = 0.5*(u1 - u2)

    ! Compute the un-stretched Gauss quadrature rule for s in [-1, 1]
    call GauLeg(-1._WP, 1._WP, n, s, w)

    ! Multiply quadrature weight by the Jacobi and do the coordinate
    ! transform
    ! After this, we obtain the stretched quadrature rule for y in [-1,1]
    do i = 1, n
      x(i) = a0 + b0*sinh(mu*s(i) - eta)
      w(i) = w(i)*b0*mu*cosh(mu*s(i) - eta)
    end do ! i

    ! Stretch the interval back to [xmin,xmax]
    w(1:n) = xl*w(1:n)
    x(1:n) = xm + xl*x(1:n)

  end subroutine GauLeg_Sinh

!**********************************************************************
! Arcsin function
! Arguments:
!  f = sinh(z) = (exp(z) + exp(-z))/2
  function MyAsinh(f) result(z)
    real(WP) :: f, z

    z = log(f + sqrt(1 + f*f))

  end function MyAsinh

!**********************************************************************

end module ModQuadRule
