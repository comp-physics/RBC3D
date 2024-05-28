module ModSpline

  use ModDataTypes
  use ModDataStruct
  use ModBasicMath
  use ModFFT
  use ModPolarPatch

  implicit none

  private

  public :: Spline_Create, &
            Spline_Destroy, &
            Spline_Build, &
            Spline_Build_on_Sphere, &
            Spline_Interp, &
            Spline_FindProjection

contains

!**********************************************************************
! Create a spline object
  subroutine Spline_Create(spln, m, n, nvar)
    type(t_spline) :: spln
    integer :: m, n, nvar

    spln%m = m
    spln%n = n
    spln%nvar = nvar

    allocate (spln%u(0:m - 1, 0:n - 1, nvar), &
              spln%u1(0:m - 1, 0:n - 1, nvar), &
              spln%u2(0:m - 1, 0:n - 1, nvar), &
              spln%u12(0:m - 1, 0:n - 1, nvar), &
              spln%kx(nvar), spln%ky(nvar))

    spln%hx = TWO_PI/m
    spln%hy = TWO_PI/n

    spln%ihx = 1./spln%hx
    spln%ihy = 1./spln%hy

    spln%kx = 0.
    spln%ky = 0.

  end subroutine Spline_Create

!**********************************************************************
! Destroy a spline object
  subroutine Spline_Destroy(spln)
    type(t_spline) :: spln

    spln%m = 0
    spln%n = 0
    spln%nvar = 0

    deallocate (spln%u, spln%u1, spln%u2, spln%u12)
    deallocate (spln%kx, spln%ky)

  end subroutine Spline_Destroy

!**********************************************************************
! Build a spline representation of a periodic function
! Arguments:
!  spln -- a spline object
!  f -- input array
!  dfx(:), dfy(:) -- the change of f over one period in each direction
! Note:
!  The shape of f must the same as that of spln%u
  subroutine Spline_Build(spln, f, dfx, dfy)
    type(t_spline) :: spln
    real(WP), dimension(spln%m, spln%n, spln%nvar) :: f
    real(WP), dimension(spln%nvar), optional :: dfx, dfy

    real(WP), allocatable :: f1(:, :, :)
    integer :: m, n, nvar, i, j

    m = spln%m
    n = spln%n
    nvar = spln%nvar

    spln%u = f

    if (present(dfx)) then
      spln%kx = i_2PI*dfx

      do i = 1, m
      do j = 1, n
        spln%u(i, j, :) = f(i, j, :) - dfx*real(i - 1)/m
      end do ! j
      end do ! i
    else
      spln%kx = 0.
    end if

    if (present(dfy)) then
      spln%ky = i_2PI*dfy

      do i = 1, m
      do j = 1, n
        spln%u(i, j, :) = f(i, j, :) - dfy*real(j - 1)/n
      end do ! j
      end do ! i
    else
      spln%ky = 0.
    end if

    call FFT_Diff(m, n, nvar, spln%u, spln%u1, spln%u2)
    call FFT_Diff(m, n, nvar, spln%u1, u2=spln%u12)

  end subroutine Spline_Build

!**********************************************************************
! Build a cubic spline formula sphere
! Arguments:
!  spline -- sphereical spline
!  f --  variable defined on a Gaussian mesh on a sphere
! Note:
!  The first dimension of f must be half of that of spln%u
  subroutine Spline_Build_on_Sphere(spln, f)
    type(t_spline) :: spln
    real(WP), dimension(0:spln%m/2, 0:spln%n - 1, spln%nvar) :: f

    integer :: nlat, nlon, i, j

    nlat = spln%m/2
    nlon = spln%n

    do j = 0, nlon - 1
    do i = 0, nlat - 1
      spln%u(i, j, :) = f(i, j, :)
      spln%u(nlat + i, j, :) = f(nlat - i, mod(j + nlon/2, nlon), :)
    end do ! i
    end do ! j

    spln%kx = 0.
    spln%ky = 0.
    call FFT_Diff(spln%m, spln%n, spln%nvar, spln%u, spln%u1, spln%u2)
    call FFT_Diff(spln%m, spln%n, spln%nvar, spln%u1, u2=spln%u12)

  end subroutine Spline_Build_on_Sphere

!**********************************************************************
! Interpolate a function from the spline representation
! Arguments:
!  spln -- splines
!  x, y -- coordinates
!  f -- interpolated values
  subroutine Spline_Interp(spln, x, y, f)
    type(t_spline) :: spln
    real(WP) :: x, y, f(spln%nvar)

    integer :: i1, j1, i2, j2, l
    real(WP) :: s, t, cx(4), cy(4), u(4, 4)

    ! Interpolate in phi direction
    i1 = floor(x*spln%ihx)
    j1 = floor(y*spln%ihy)

    s = x*spln%ihx - i1
    t = y*spln%ihy - j1

    i1 = modulo(i1, spln%m)
    j1 = modulo(j1, spln%n)

    i2 = modulo(i1 + 1, spln%m)
    j2 = modulo(j1 + 1, spln%n)

    cx = (/1 + s*s*(-3.+2.*s), s*s*(3.-2.*s), &
           spln%hx*s*(1.+s*(-2.+s)), spln%hx*s*s*(-1.+s)/)
    cy = (/1 + t*t*(-3.+2.*t), t*t*(3.-2.*t), &
           spln%hy*t*(1.+t*(-2.+t)), spln%hy*t*t*(-1.+t)/)

    do l = 1, size(f)
      u(1, :) = (/spln%u(i1, j1, l), spln%u(i1, j2, l), &
                  spln%u2(i1, j1, l), spln%u2(i1, j2, l)/)
      u(2, :) = (/spln%u(i2, j1, l), spln%u(i2, j2, l), &
                  spln%u2(i2, j1, l), spln%u2(i2, j2, l)/)
      u(3, :) = (/spln%u1(i1, j1, l), spln%u1(i1, j2, l), &
                  spln%u12(i1, j1, l), spln%u12(i1, j2, l)/)
      u(4, :) = (/spln%u1(i2, j1, l), spln%u1(i2, j2, l), &
                  spln%u12(i2, j1, l), spln%u12(i2, j2, l)/)

      f(l) = dot_product(cx, matmul(u, cy))
    end do ! l

    ! Add the linearly increasing part
    f = f + spln%kx*x + spln%ky*y

  end subroutine Spline_Interp

!**********************************************************************
! Find the projection of a point on a spline surface
!
! Arguments:
!  spln -- spline
!  xTar -- some arbitrary point
!  th0, phi0 -- the local coordinate of x0
!  x0 -- the physical coordinate of the surface point closest to xTar
! Note:
!  -- A good initial guess is needed for (th0, phi0)
  subroutine Spline_FindProjection(spln, xtar, th0, phi0, x0)
    type(t_Spline) :: spln
    real(WP) :: xtar(3), th0, phi0, x0(3)

    integer :: iter
    integer, parameter :: iterMax = 3
    integer, parameter :: nth = 2   ! number of points on the patch arond (th0, phi0)
    integer, parameter :: nphi = 8
    real(WP) :: h
    real(WP) :: thPat_L(nth), phiPat_L(nphi), thPat(nth, nphi), phiPat(nth, nphi)
    real(WP) :: xyGq_L(0:nth*nphi, 2), dist2Gq(0:nth*nphi)
    real(WP) :: c0, c1, c2, c11, c12, c22, xyMin_L(2), dist2MinEst
    real(WP) :: xmin(3), thMin_L, phiMin_L, thMin, phiMin, dist2min
    real(WP) :: xx(3)
    integer :: ith, iphi, i

    ! Initialize x0
    call Spline_Interp(spln, th0, phi0, x0)

    ! Build a local patch around x0
    h = max(spln%hx, spln%hy)

    do iter = 1, iterMax
      thPat_L = (/(ith*h/nth, ith=1, nth)/)
      phiPat_L = (/((iphi - 1)*two_pi/nphi, iphi=1, nphi)/)
      call PolarPatch_Build(th0, phi0, thPat_L, phiPat_L, thPat, phiPat)

      ! Compute the local 2D Cartesian coordinate of patch points
      xyGq_L(0, :) = 0.
      call Spline_Interp(spln, th0, phi0, xx)
      xx = xx - xtar
      dist2Gq(0) = sum(xx**2)

      do iphi = 1, nphi
        do ith = 1, nth
          i = ith + nth*(iphi - 1)

          xyGq_L(i, 1) = thPat_L(ith)*cos(phiPat_L(iphi))
          xyGq_L(i, 2) = thPat_L(ith)*sin(phiPat_L(iphi))

          call Spline_Interp(spln, thPat(ith, iphi), phiPat(ith, iphi), xx)
          xx = xx - xtar
          dist2Gq(i) = sum(xx**2)
        end do ! ith
      end do ! iphi

      ! Find a quadartic interpolation of distance to xtar
      call QuadFit_2D(xyGq_L, dist2Gq, c0, c1, c2, c11, c12, c22)

      call Min_Quad_2D(c0, c1, c2, c11, c12, c22, xyMin_L, dist2MinEst)

      thMin_L = sqrt(sum(xyMin_L**2))
      phiMin_L = atan2(xyMin_L(2), xyMin_L(1))
      call PolarPatch_Map(th0, phi0, thMin_L, phiMin_L, thMin, phiMin)
      call Spline_Interp(spln, thMin, phiMin, xMin)
      xx = xMin - xtar
      dist2Min = sum(xx**2)

      if (dist2Min > dist2Gq(0)) exit    ! Fail to converge

      ! Reduce the search radius
      th0 = thMin
      phi0 = phiMin
      x0 = xMin
      h = 0.5*h
    end do ! iter

  end subroutine Spline_FindProjection

!**********************************************************************

end module ModSpline
