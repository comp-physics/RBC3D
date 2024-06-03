! Subroutines for rigid wall boundary conditions
module ModWall

  use ModDataTypes
  use ModBasicMath
  use ModDataStruct

  implicit none

  private

  public :: Wall_Create, &
            Wall_Destroy, &
            Wall_Build_V2V, &
            Wall_ComputeGeometry

contains

!**********************************************************************
! Create a wall
! Arguments:
!  wall --
!  nvert, nele -- number of vertices and elements
  subroutine Wall_Create(wall, nvert, nele)
    type(t_wall) :: wall
    integer :: nvert, nele

    wall%nvert = nvert
    wall%nele = nele

    allocate (wall%x(nvert, 3))
    allocate (wall%e2v(nele, 3), wall%v2v(nvert))

    allocate (wall%a3(nele, 3), wall%area(nele), wall%epsDist(nele))
    allocate (wall%f(nvert, 3), wall%g(nvert, 3))

    allocate (wall%indxVertGlb(nvert))

    ! Default ID
    wall%id = -1

  end subroutine Wall_Create

!**********************************************************************
! Destroy a wall
! Arguments:
!   wall --
  subroutine Wall_Destroy(wall)
    type(t_wall) :: wall

    deallocate (wall%x)
    deallocate (wall%e2v, wall%v2v)

    deallocate (wall%a3, wall%area, wall%epsDist)
    deallocate (wall%f, wall%g)

    deallocate (wall%indxVertGlb)

  end subroutine Wall_Destroy

!**********************************************************************
! Build the v2v(:) array
! Arguments:
!  wall --
!  Lb -- periodic box size
  subroutine Wall_Build_V2v(wall, Lb)
    type(t_wall) :: wall
    real(WP) :: Lb(3)

    integer :: nvertBd
    integer, allocatable :: ivertBd(:)
    real(WP) :: iLb(3), xmin(3), xmax(3), eps, dx1(3), dx2(3), xx(3)
    integer :: ivert, p1, p2, ivert1, ivert2

    ! Identify all candidate points
    iLb = 1./Lb
    eps = minval(Lb)*1.E-5

    nvertBd = 0
    allocate (ivertBd(wall%nvert))

    xmin = minval(wall%x, dim=1)
    xmax = maxval(wall%x, dim=1)
    do ivert = 1, wall%nvert
      dx1 = wall%x(ivert, :) - xmin
      dx2 = wall%x(ivert, :) - xmax
      if (minval(abs(dx1)) < eps .or. minval(abs(dx2)) < eps) then
        nvertBd = nvertBd + 1
        ivertBd(nvertBd) = ivert
      end if
    end do ! ivert

    ! Build v2v list
    wall%v2v = 0
    do p1 = 1, nvertBd - 1
      ivert1 = ivertBd(p1)
      if (wall%v2v(ivert1) > 0) cycle

      do p2 = p1 + 1, nvertBd
        ivert2 = ivertBd(p2)

        xx = wall%x(ivert2, :) - wall%x(ivert1, :)
        xx = xx - nint(xx*iLb)*Lb

        if (maxval(abs(xx)) < eps) then
          wall%v2v(ivert2) = ivert1
        end if
      end do ! p2
    end do ! p1

    ! Deallocate working arrays
    deallocate (ivertBd)

  end subroutine Wall_Build_V2V

!**********************************************************************
! Compute wall geometry
  subroutine Wall_ComputeGeometry(wall)
    type(t_wall) :: wall

    integer :: iele, ivert, l
    real(WP) :: xele(3, 3), x12(3), x13(3)
    real(WP) :: a3(3), a3N

    do iele = 1, wall%nele
      do l = 1, 3
        ivert = wall%e2v(iele, l)
        xele(l, :) = wall%x(ivert, :)
      end do ! l

      x12 = xele(2, :) - xele(1, :)
      x13 = xele(3, :) - xele(1, :)

      a3 = CrossProd(x12, x13)
      a3N = sqrt(sum(a3*a3))
      a3 = a3/a3N

      wall%a3(iele, :) = a3
      wall%area(iele) = 0.5*a3N
      wall%epsDist(iele) = sqrt(wall%area(iele))
    end do ! iele

    wall%areaTot = sum(wall%area)

  end subroutine Wall_ComputeGeometry

!**********************************************************************

end module ModWall
