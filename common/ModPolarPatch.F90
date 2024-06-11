! Construction of local polar coordinate patch
module ModPolarPatch

  use ModDataTypes
  use ModDataStruct
  use ModBasicMath
  use ModQuadRule

  implicit none

  private

  public :: RbcPolarPatch_Create, &
            RbcPolarPatch_Destroy, &
            PolarPatch_Build, &
            PolarPatch_FindPoints, &
            PolarPatch_Map, &
            DistOnSphere

contains

!**********************************************************************
! Build a RBC surface polar patch
! Argument:
!   patch --
!   rbc -- a cell
  subroutine RbcPolarPatch_Create(patch, rbc)
    type(t_RbcPolarPatch) :: patch
    type(t_Rbc) :: rbc

    integer :: nlat, nlon, ilat, ilon
    integer :: nrad, nazm, irad, iazm
    real(WP) :: radius, h

    ! Set up parameters
    nlat = rbc%nlat
    nlon = rbc%nlon

    patch%nlat = nlat
    patch%nlon = nlon

    radius = PI/sqrt(real(nlat))
    h = PI/rbc%nlat
    nrad = 2*nint(radius/h)
    nazm = 2*nrad

    ! Copy values
    patch%nlat = nlat
    patch%nlon = nlon
    patch%radius = radius
    patch%nrad = nrad
    patch%nazm = nazm

    ! Allocate arrays
    allocate (patch%thL(nrad), patch%phiL(nazm), patch%w(nrad))
    allocate (patch%thG(nrad, nazm, nlat, nlon), &
              patch%phiG(nrad, nazm, nlat, nlon)); 
    ! Compute local patch
    call GauLeg(0._WP, radius, nrad, patch%thL, patch%w)
    do irad = 1, nrad
      patch%w(irad) = patch%w(irad)*sin(patch%thL(irad))*(TWO_PI/nazm)
      patch%w(irad) = patch%w(irad)*Maskfunc(patch%thL(irad)/radius)
    end do ! irad

    patch%phiL = (/((iazm - 1)*TWO_PI/nazm, iazm=1, nazm)/)

    ! compute global patch
    do ilon = 1, nlon
    do ilat = 1, nlat
      call PolarPatch_Build(rbc%th(ilat), rbc%phi(ilon), patch%thL, patch%phiL, &
                            patch%thG(:, :, ilat, ilon), patch%phiG(:, :, ilat, ilon))
    end do ! ilat
    end do ! ilon

  end subroutine RbcPolarPatch_Create

!**********************************************************************
! Destroy a RBC polar patch
  subroutine RbcPolarPatch_Destroy(patch)
    type(t_RbcPolarPatch) :: patch

    if (associated(patch%thL)) then
      deallocate (patch%thL, patch%phiL, patch%w)
    end if

    if (associated(patch%thG)) then
      deallocate (patch%thG, patch%phiG)
    end if

  end subroutine RbcPolarPatch_Destroy

!**********************************************************************
! Build up a local polar patch
! Arguments:
!  (th0, pi0) -- origin of the patch
!  thL(:), phiL(:) -- local polar coordinates
!  thG(:,:), phiG(:,:) -- mesh coordinates after the patch is
!       mapped to the surface
  subroutine PolarPatch_Build(th0, phi0, thL, phiL, thG, phiG)
    real(WP) :: th0, phi0
    real(WP) :: thL(:), phiL(:)
    real(WP) :: thG(:, :), phiG(:, :)

    integer :: nth, nphi, i, j
    real(WP) :: sin_th0, cos_th0
    real(WP), allocatable :: sin_th(:), cos_th(:), sin_phi(:), cos_phi(:)
    real(WP) :: A(3, 3), x(3)

    ! Allocate working arrays
    nth = size(thL)
    nphi = size(phiL)
    allocate (sin_th(nth), cos_th(nth), sin_phi(nphi), cos_phi(nphi))

    ! Pre-compute the expensive trigonomic functions
    sin_th = sin(thL)
    cos_th = cos(thL)

    sin_phi = sin(phiL)
    cos_phi = cos(phiL)

    ! Compute the rotation matrix, assuming phi0 = 0
    sin_th0 = sin(th0)
    cos_th0 = cos(th0)

    A = 0.
    A(1, 1) = cos_th0; A(1, 3) = sin_th0
    A(2, 2) = 1.
    A(3, 1) = -sin_th0; A(3, 3) = cos_th0

    do i = 1, nth
    do j = 1, nphi
      x = (/sin_th(i)*cos_phi(j), sin_th(i)*sin_phi(j), cos_th(i)/)
      x = matmul(A, x)

      x(3) = max(-1.0_WP, min(1.0_WP, x(3)))
      thG(i, j) = acos(x(3))

      ! Add back phi0
      phiG(i, j) = atan2(x(2), x(1)) + phi0
    end do ! i
    end do ! j

    phiG = phiG - floor(phiG*I_2PI)*TWO_PI

    ! Deallocate working arrays
    deallocate (sin_th, cos_th, sin_phi, cos_phi)

  end subroutine PolarPatch_Build

!**********************************************************************
! Find the Cartesian mesh points within a polar patch
! Arguments:
!  (th0, phi0) -- origin of the patch
!  r0 -- radius of the patch
!  ths(:), phis(:) -- the global mesh coordinates
!  n -- number of mesh points on the patch
!  ijs -- indices of mesh points on the patch
! Note:
!  -- phis(:) is uniform in [0, 2*pi]
!  -- For a planar surface, ths(:) is also uniform in [0, 2*pi]
  subroutine PolarPatch_FindPoints(th0, phi0, r0, ths, phis, n, ijs)
    real(WP) :: th0, phi0, r0
    integer :: type_surf
    real(WP) :: ths(:), phis(:)
    integer :: n, ijs(:, :)

    real(WP), parameter :: eps = 1.E-10
    integer :: nth, nphi
    real(WP) :: h_th, ih_th, h_phi, ih_phi, dth, dphi
    real(WP) :: cosr0
    integer :: i, j, imin, imax, jmin, jmax

    ! Set up parameters
    nth = size(ths)
    nphi = size(phis)

    h_phi = TWO_PI/nphi
    ih_phi = 1.0/h_phi

    n = 0

    cosr0 = cos(r0)

    do i = 1, nth
      if (ths(i) <= th0 - r0) cycle
      if (ths(i) >= th0 + r0) exit

      dth = ths(i) - th0

      dphi = 1.-(cos(dth) - cosr0)/((sin(ths(i)) + eps)*(sin(th0) + eps))
      dphi = max(-1._WP, min(1._WP, dphi))
      dphi = acos(dphi)

      if (dphi > pi - eps) then
        jmin = 1
        jmax = nphi
      else
        jmin = ceiling((phi0 - dphi - phis(1))*ih_phi) + 1
        jmax = floor((phi0 + dphi - phis(1))*ih_phi) + 1
      end if

      do j = jmin, jmax
        n = n + 1
        ijs(n, 1) = i
        ijs(n, 2) = modulo(j - 1, nphi) + 1
      end do ! j
    end do ! i

  end subroutine PolarPatch_FindPoints

!**********************************************************************
! Transform a local polar coordinate to the global coordinate on the
! reference surface
! Arguments:
!  th0, phi0 -- reference coordinate of the origin
!  dth, dphi -- relative polar coordinate
!  th, phi -- transformed coordinate on the reference surface
  subroutine PolarPatch_Map(th0, phi0, dth, dphi, th, phi)
    real(WP) :: th0, phi0, dth, dphi, th, phi

    real(WP) :: s1(3), s2(3), x0(3), x(3)

    s1 = (/cos(th0)*cos(phi0), cos(th0)*sin(phi0), -sin(th0)/)
    s2 = (/-sin(phi0), cos(phi0), 0._WP/)
    x0 = (/sin(th0)*cos(phi0), sin(th0)*sin(phi0), cos(th0)/)

    ! Compute the coordinate before rotation
    x = (/sin(dth)*cos(dphi), sin(dth)*sin(dphi), cos(dth)/)

    ! The mapping, which is a rotation transforms
    !  (1,0,0) to s1
    !  (0,1,0) to s2
    !  (0,0,1) to x0
    x = x(1)*s1 + x(2)*s2 + x(3)*x0
    x = x/sqrt(dot_product(x, x))

    x(3) = max(-1._WP, min(1._WP, x(3)))
    th = acos(x(3))
    phi = atan2(x(2), x(1))
    if (phi < 0) phi = phi + TWO_PI

  end subroutine PolarPatch_Map

!**********************************************************************
! The shortest distance between two points on the surface
! Arguments:
!  th0, phi0 -- first point
!  th1, phi1 -- second point
  function DistOnSphere(th0, phi0, th1, phi1) result(dth)
    real(WP) :: th0, phi0, th1, phi1
    real(WP) :: dth

    dth = cos(th0 - th1) - sin(th0)*sin(th1)*(1.-cos(phi0 - phi1))
    dth = min(1._WP, max(-1._WP, dth))
    dth = acos(dth)

  end function DistOnSphere

!**********************************************************************

end module ModPolarPatch
