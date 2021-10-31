module ModTargetList
  
  use ModDataTypes
  use ModDataStruct
  use ModConf

  implicit none

  private

  private :: SetActiveFlag

  public :: TargetList_Create, &
    TargetList_Destroy, &
    TargetList_Update, &
    TargetList_CreateFromRaw, &
    TargetList_CollectArray

contains

!**********************************************************************
! Create target list from cells and walls
! Note:
!   We include every points in the target list to make MPI_AllGather
!   easier
  subroutine TargetList_Create(list, rbcs, walls)
    type(t_TargetList) :: list
    type(t_Rbc),target,optional :: rbcs(:)
    type(t_Wall),target,optional :: walls(:)

    type(t_Rbc),pointer :: rbc
    type(t_Wall),pointer :: wall
    integer :: nPoint, p, irbc, ilat, ilon, iwall, ivert

    nPoint = 0
    if (present(rbcs)) then
      do irbc = 1, size(rbcs)
        rbc => rbcs(irbc)
        nPoint = nPoint + rbc%nlat*rbc%nlon
      end do ! i
    end if
    if (present(walls)) then
      do iwall = 1, size(walls)
        wall => walls(iwall)
        nPoint = nPoint + wall%nvert
      end do ! i
    end if

    ! Allocate arrays
    list%nPoint = nPoint

    allocate(list%x(nPoint,3), list%indx(nPoint,0:2), list%active(nPoint) )
    allocate(list%lam(nPoint),list%Acoef(nPoint))

    ! Set up indx
    p = 0;

    if (present(rbcs)) then
      do irbc = 1, size(rbcs)
        rbc => rbcs(irbc)

    do ilon = 1, rbc%nlon
        do ilat = 1, rbc%nlat
      p = p + 1
      list%indx(p,0:2) = (/ rbc%ID, ilat, ilon /)
    end do ! ilat
    end do ! ilon
      end do ! irbc
    end if

    if (present(walls)) then
      do iwall = 1, size(walls)
        wall => walls(iwall)

    do ivert = 1, wall%nvert
      p = p + 1
      list%indx(p,0:2) = (/ wall%Id, ivert, -1 /)
    end do ! ivert
      end do ! i
    end if

  end subroutine TargetList_Create

!**********************************************************************
! Destroy target list
  subroutine TargetList_Destroy(list)
    type(t_TargetList) :: list

    deallocate(list%x,list%lam, list%Acoef)
    deallocate(list%indx, list%active)

  end subroutine TargetList_Destroy

!**********************************************************************
! Update target list
  subroutine TargetList_Update(list, rbcs, walls)
    type(t_TargetList) :: list
    type(t_Rbc),target,optional :: rbcs(:)
    type(t_Wall),target,optional :: walls(:)

    type(t_Rbc),pointer :: rbc
    type(t_Wall),pointer :: wall
    integer :: p, irbc, nvert, iwall, i
    integer,allocatable :: whichNode(:)
    integer :: nodeNum, numNodes
    integer :: ierr

    p = 0

    if (present(rbcs)) then
      do irbc = 1, size(rbcs)
        rbc => rbcs(irbc)

    nvert = rbc%nlat * rbc%nlon
    list%lam(p+1:p+nvert) = viscRat(rbc%celltype)   !COEF
    list%Acoef(p+1:p+nvert) = Acoef(rbc%celltype)   !COEF
    list%x(p+1:p+nvert,:) = reshape(rbc%x, (/nvert,3/))
    p = p + nvert
      end do ! irbc
    end if

    if (present(walls)) then
      do iwall = 1, size(walls)
        wall => walls(iwall)

    nvert = wall%nvert
    list%lam(p+1:p+nvert) = 1. !COEF
    list%Acoef(p+1:p+nvert) = 2.  !COEF
    list%x(p+1:p+nvert,:) = wall%x
    p = p + nvert
      end do ! iwall
    end if

    call SetActiveFlag(list)

  end subroutine TargetList_Update

!**********************************************************************
! Create target list from raw points
  subroutine TargetList_CreateFromRaw(list, x)
    type(t_TargetList) :: list
    real(WP) :: x(:,:)

    integer :: nPoint, i

!!$    nPoint = clist%nPoint+size(x,1)

    nPoint = size(x,1)
    list%nPoint = nPoint
    allocate(list%x(nPoint,3), list%indx(nPoint,0:2), list%active(nPoint))
    allocate(list%lam(nPoint),list%Acoef(nPoint))

!!$    list%x(1:clist%nPoint,:) = clist%x
!!$    list%indx(1:clist%nPoint,:) = clist%indx
!!$    list%lam(1:clist%nPoint) = clist%lam
!!$    list%Acoef(1:clist%nPoint) = clist%Acoef
!!$
!!$    list%x(clist%nPoint+1:nPoint,:) = x
!!$    list%indx(clist%nPoint+1:nPoint,:) = -1
!!$    list%lam(clist%nPoint+1:nPoint) = 1.
!!$    list%Acoef(clist%nPoint+1:nPoint) = 2.

    list%x = x
    list%indx = -1
    list%lam = 1.
    list%Acoef = 2.

    call SetActiveFlag(list)

  end subroutine TargetList_CreateFromRaw

!**********************************************************************
  subroutine TargetList_CollectArray(list, nvar, a, COMM)
    type(t_TargetList) :: list
    integer :: nvar
    real(WP) :: a(list%npoint,nvar)
    integer :: COMM

    integer :: nsend
    integer,allocatable :: ia(:)
    real(WP),allocatable :: asend(:,:)
    integer :: p, i

    ! Allocate working arrays
    nsend = count(list%active)
    allocate(ia(nsend), asend(nsend,nvar))

    p = 0
    do i = 1, list%npoint
      if (list%active(i)) then
        p = p + 1
    ia(p) = i
    asend(p,:) = a(i,:)
      end if
    end do ! i

    a = 0.
    call CollectArray(nsend, nvar, ia, asend, a, COMM, MPI_SUM)

    ! Deallocate working arrays
    deallocate(ia, asend)

  end subroutine TargetList_CollectArray

!**********************************************************************
  subroutine SetActiveFlag(list)
    type(t_TargetList) :: list

    integer,allocatable :: whichNode(:)
    integer :: nodeNum, numNodes
    integer :: i
    integer :: ierr

    ! The head node determines which processor every target point belongs to
    ! and broadcast the information
    allocate(whichNode(list%nPoint))

    call MPI_Comm_Size(MPI_Comm_Ewald, numNodes, ierr)
    call MPI_Comm_Rank(MPI_Comm_Ewald, nodeNum, ierr)
    if (rootEwald) then
      do i = 1, list%nPoint
        whichNode(i) = floor(list%x(i,3)*iLb(3)*numNodes)
        whichNode(i) = modulo(whichNode(i), numNodes)
      end do ! i
    end if
    call MPI_Bcast(whichNode, list%nPoint, MPI_Integer, 0, MPI_Comm_Ewald, ierr)

    do i = 1, list%nPoint
      list%active(i) = (whichNode(i) == nodeNum)
    end do ! i

    deallocate(whichNode)

  end subroutine SetActiveFlag

!**********************************************************************

end module ModTargetList
