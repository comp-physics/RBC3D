! Data
module ModData

  use MPI
  use ModDataTypes
  use ModDataStruct
  use ModConf
  use ModSourceList
  use ModTargetList

  implicit none

  integer :: nrbc, nwall
  type(t_rbc),allocatable,target :: rbcRefs(:)
  type(t_RbcPolarPatch),target :: rbcPatch
  type(t_rbc),allocatable,target :: rbcs(:)
  type(t_wall),allocatable,target :: walls(:)

  type(t_SourceList),target :: slist_rbc, slist_wall
  type(t_TargetList),target :: tlist_rbc, tlist_wall

  public :: GlobData_Init, &
    GlobData_Finalize, &
    SyncSurfCoord

contains

!**********************************************************************
! Set up some global data record
! Note:
!  -- The allocation of rbcs(:) and walls(:) is done by user
  subroutine GlobData_Init

    integer :: irbc, iwall, ivert, p
    type(t_wall),pointer :: wall

    ! Assign surface Id
    p = 0

    do irbc = 1, nrbc
      p = p + 1
      rbcs(irbc)%id = p
    end do ! irbc

    do iwall = 1, nwall
      p = p + 1
      walls(iwall)%id = p
    end do ! iwall

    ! Assign global index to wall mesh points
    p = 0
    do iwall = 1, nwall
      wall => walls(iwall)

      do ivert = 1, wall%nvert
        if (wall%v2v(ivert) == 0) then
      p = p + 1
      wall%indxVertGlb(ivert) = p
    else
      wall%indxVertGlb(ivert) = wall%indxVertGlb(wall%v2v(ivert))
    end if
      end do ! ivert
    end do ! iwall

    ! Create global source and target point list
    call SourceList_Create(slist_rbc, Lb, rc_Ewd, rbcs=rbcs)
    call TargetList_Create(tlist_rbc, rbcs=rbcs)

    call SourceList_Create(slist_wall, Lb, rc_Ewd, walls=walls)
    call TargetList_Create(tlist_wall, walls=walls)


  end subroutine GlobData_Init

!**********************************************************************
! Finalize some global data record
! Note:
!  -- The deallocation of rbcs(:) and walls(:) is done by user
  subroutine GlobData_Finalize

    call SourceList_Destroy(slist_rbc)
    call TargetList_Destroy(tlist_rbc)

    call SourceList_Destroy(slist_wall)
    call TargetList_Destroy(tlist_wall)

  end subroutine GlobData_Finalize

!**********************************************************************
! Synchronize cells and walls
  subroutine SyncSurfCoord

    integer :: irbc, iwall
    type(t_Rbc),pointer :: rbc
    type(t_Wall),pointer :: wall
    integer :: ierr

    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      call MPI_BCast(rbc%x, size(rbc%x), MPI_WP, 0, MPI_Comm_World, ierr)
    end do ! irbc

    do iwall = 1, nwall
      wall => walls(iwall)
      call MPI_BCast(wall%f, size(wall%f), MPI_INTEGER, 0, MPI_Comm_World, ierr)
    end do ! iwall

  end subroutine SyncSurfCoord

!**********************************************************************

end module ModData
