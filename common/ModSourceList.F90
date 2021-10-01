module ModSourceList

  use ModDataTypes
  use ModDataStruct
  use ModHashTable
  use ModConf

  implicit none

  private

  public :: SourceList_Create, &
    SourceList_Destroy, &
    SourceList_UpdateCoord, &
    SourceList_UpdateDensity

contains

!**********************************************************************
! Init the source list 
! Arguments:
!  list --
!  Lb -- size of the box
!  rc -- cut-off distance
!  rbcs -- 
!  walls --
  subroutine SourceList_Create(list, Lb, rc, rbcs, walls)
    type(t_SourceList) :: list
    real(WP) :: Lb(3), rc
    type(t_Rbc),target,optional :: rbcs(:)
    type(t_Wall),target,optional :: walls(:)

    integer :: Np, Nc(3)
    real(WP) :: iLbNc(3)
    integer :: irbc, iwall
    type(t_Rbc),pointer :: rbc
    type(t_Wall),pointer :: wall

    ! Count the number of points
    Np = 0

    if (present(rbcs)) then
      do irbc = 1, size(rbcs)
        rbc => rbcs(irbc)
        Np = Np + rbc%nlat * rbc%nlon
      end do ! irbc
    end if

    if (present(walls)) then
      do iwall = 1, size(walls)
        wall => walls(iwall)
        Np = Np + wall%nele
      end do ! iwall
    end if

    ! Allocate arrays
    list%nPoint = Np
    allocate(list%x(Np,3), list%indx(Np,0:2) )
    allocate(list%a3(Np,3), list%f(Np,3), list%g(Np,3) )
    allocate(list%lam(Np), list%Bcoef(Np) )

    ! Create the Hash table
    call HashTable_ComputeNumBlocks(rc, Nc, iLbNc)
    list%Nc = Nc
    list%iLbNc = iLbNc
    allocate(list%hoc(0:Nc(1)+1, 0:Nc(2)+1, 0:Nc(3)+1), list%next(Np)) 

  end subroutine SourceList_Create

!**********************************************************************
! Destroy a source list
  subroutine SourceList_Destroy(list)
    type(t_SourceList) :: list

    integer :: ierr

    deallocate(list%x, list%indx)
    deallocate(list%a3, list%f, list%g)
    deallocate(list%lam, list%Bcoef)
    deallocate(list%hoc, list%next)

  end subroutine SourceList_Destroy

!**********************************************************************
! Update the point coordinates of a source list
! Arguments:
!  list --
!  rbcs -- 
!  walls --
! Note:
!  Only active source points are included
  subroutine SourceList_UpdateCoord(list, rbcs, walls)
    type(t_SourceList) :: list
    type(t_Rbc),target,optional :: rbcs(:)
    type(t_Wall),target,optional :: walls(:)

    type(t_Rbc),pointer :: rbc
    type(t_Wall),pointer :: wall
    integer :: p, irbc, ilat, ilon, iwall, iele, l
    real(WP) :: xele(3,3)

    ! Update points in the list
    p = 0

    if (present(rbcs)) then
      do irbc = 1, size(rbcs)
        rbc => rbcs(irbc)

    if (Cell_Has_Source(rbc) ) then
      do ilon = 1, rbc%nlon
      do ilat = 1, rbc%nlat
        p = p + 1
        list%x(p,:) = rbc%x(ilat,ilon,:)
        list%a3(p,:) = rbc%a3(ilat,ilon,:)

        list%indx(p,0:2) = (/ rbc%id, ilat, ilon /)

            list%lam(p) =  viscRat(rbc%celltype)  !COEF
            list%Bcoef(p) = Bcoef(rbc%celltype)  !COEF
 
      end do ! ilat
      end do ! ilon
        end if
      end do ! irbc
    end if

    if (present(walls)) then
      do iwall = 1, size(walls)
        wall => walls(iwall)

    do iele = 1, wall%nele
      do l = 1, 3
        xele(l,:) = wall%x(wall%e2v(iele,l),:)
      end do ! l

      if (Tri_Has_Source(xele)) then
        p = p + 1
        list%x(p,:) = THRD*sum(xele,dim=1)
        list%indx(p,0:2) = (/ wall%ID, iele, -1 /)
            list%lam(p) = 1.
            list%Bcoef(p) = 0.  !COEF
      end if
    end do ! iele
      end do ! iwall
    end if

    ! Rebuild the Hash table
    list%nPoint = p
    call HashTable_Build(list%Nc, list%iLbNc, list%nPoint, list%x, list%hoc, list%next)

  end subroutine SourceList_UpdateCoord

!**********************************************************************
! Update sources
! Arguments:
!  UpdateF, UpdateG -- whether to update F and G
! Note:
!  -- We don't need to update density for wall points, because the surface
!     integration is done by Gauss quadrature on each element.
  subroutine SourceList_UpdateDensity(list, rbcs, UpdateF, UpdateG)
    type(t_SourceList) :: list
    type(t_Rbc),target :: rbcs(:)
    logical,optional :: UpdateF, UpdateG

    logical :: copyF, copyG
    type(t_Rbc),pointer :: rbc
    integer :: p, irbc, ilat, ilon
    real(WP) :: dS

    CopyF = .false.;    if (present(UpdateF)) CopyF = UpdateF
    CopyG = .false.;    if (present(UpdateG)) CopyG = UpdateG

    ! Update f and g arrays
    do p = 1, list%nPoint
      irbc = list%indx(p,0) - rbcs(1)%ID + 1
      ilat = list%indx(p,1)
      ilon = list%indx(p,2)

      rbc => rbcs(irbc)

      dS = rbc%detj(ilat,ilon)*rbc%w(ilat)

      if (CopyF) list%f(p,:) = rbc%f(ilat,ilon,:)*dS
      if (CopyG) list%g(p,:) = rbc%g(ilat,ilon,:)*dS
    end do ! p

  end subroutine SourceList_UpdateDensity

!**********************************************************************

end module ModSourceList
