module ModHashTable

  use ModDataTypes
  use ModDataStruct
  use ModConf

  implicit none

  private

  public :: HashTable_Build, &
            HashTable_Index, &
            HashTable_ComputeNumBlocks

contains

!**********************************************************************
! Update cell list
! Arguments:
!  Nc -- number of cells
!  iLbNnc --
!  x -- point coordinates
  subroutine HashTable_Build(Nc, iLbNc, nPoint, x, hoc, next)
    integer :: Nc(3)
    real(WP) :: iLbNc(3)
    integer :: nPoint
    real(WP) :: x(:, :)
    integer :: hoc(0:Nc(1) + 1, 0:Nc(2) + 1, 0:Nc(3) + 1)
    integer :: next(:)

    integer :: i, i1, i2, i3

    ! Initialize
    hoc = -1
    next = -1

    do i = 1, nPoint
      call HashTable_Index(Nc, iLbNc, x(i, :), i1, i2, i3)

      if (i3 >= 0 .and. i3 <= Nc(3) + 1) then
        next(i) = hoc(i1, i2, i3)
        hoc(i1, i2, i3) = i
      end if
    end do ! i

    ! Ghost cells
    hoc(0, :, :) = hoc(Nc(1), :, :)
    hoc(Nc(1) + 1, :, :) = hoc(1, :, :)

    hoc(:, 0, :) = hoc(:, Nc(2), :)
    hoc(:, Nc(2) + 1, :) = hoc(:, 1, :)

    if (SINGLE_NODE) then
      hoc(:, :, 0) = hoc(:, :, Nc(3))
      hoc(:, :, Nc(3) + 1) = hoc(:, :, 1)
    end if

  end subroutine HashTable_Build

!**********************************************************************
! Compute the Hash table index of a point
! Argument:
!  Nc -- number of cell blocks
!  iLbNc -- Nc/(local domain size)
!  x -- point coordinate
!  i1, i2, i3 -- Hash index of x(:)
  subroutine HashTable_Index(Nc, iLbNc, x, i1, i2, i3)
    integer :: Nc(3)
    real(WP) :: iLbNc(3)
    real(WP) :: x(3)
    integer :: i1, i2, i3

    real(WP) :: z, zc

    i1 = modulo(floor(x(1)*iLbNc(1)), Nc(1)) + 1
    i2 = modulo(floor(x(2)*iLbNc(2)), Nc(2)) + 1

    ! The z-direction needs special care
    if (SINGLE_NODE) then
      i3 = modulo(floor(x(3)*iLbNc(3)), Nc(3)) + 1
    else
      ! Translate along z-direction to make the point as close to the local domain
      ! as possible
      zc = 0.5*(nodeZmin + nodeZmax)
      z = x(3) - zc
      z = zc + (z - nint(z*iLb(3))*Lb(3))    ! z is the translated x(3)

      i3 = floor((z - nodeZmin)*iLbNc(3)) + 1
    end if

  end subroutine HashTable_Index

!**********************************************************************
! Compute the number of Hash list cells
! Arguments:
!  rc -- cut-off distance
!  Nc -- number of cells
!  iLbNc -- Nc/Lb
  subroutine HashTable_ComputeNumBlocks(rc, Nc, iLbNc)
    real(WP) :: rc
    integer :: Nc(3)
    real(WP) :: iLbNc(3)

    integer :: numNodes, ierr

    call MPI_Comm_Size(MPI_Comm_Ewald, numNodes, ierr)

    Nc(1) = floor(Lb(1)/rc)
    Nc(2) = floor(Lb(2)/rc)
    Nc(3) = floor((nodeZmax - nodeZmin)/rc)

    ! Constraint
    Nc(1) = max(Nc(1), 3)
    Nc(2) = max(Nc(2), 3)
    if (numNodes == 1) then
      Nc(3) = max(Nc(3), 3)
    else if (numNodes == 2) then
      Nc(3) = max(Nc(3), 2)
    else
      Nc(3) = max(Nc(3), 1)
    end if

    iLbNc(1:2) = iLb(1:2)*Nc(1:2)
    iLbNc(3) = Nc(3)/(nodeZmax - nodeZmin)

  end subroutine HashTable_ComputeNumBlocks

!**********************************************************************

end module ModHashTable
