! Add repulsion to prevent cell membranes from pentretrating each other
! and self-penetration
module ModRepulsion

  use ModDataTypes
  use ModDataStruct
  use ModConf
  use ModData
  use ModHashTable
  use ModPolarPatch
  use ModSpline
  use ModIntOnWalls, only: MinDistToTri

  implicit none

  private

  public :: AddIntraCellRepulsionForce
  public :: InterCellRepulsion
  public :: LeukWallRepulsion
  public :: AddR0GravityForce
  public :: AddR0Motion

  private :: Closest_Neighbor_Cell, &
             Closest_Neighbor_Wall

contains

!**********************************************************************
! Add repulsion force to prevent cells from self-overlapping
  subroutine AddR0Motion

    type(t_rbc), pointer :: rbc
    integer :: i, j, irbc, ilat, ilon
    integer :: ierr

    real(WP), dimension(3) :: nc

    if (rootworld) then
!       print *,"ADDING R0 GRAVITY FORCE ",fo
      print *, "MOVING CELL 1!!!!! +0.0005"
    end if

    nc(3) = 0.
    do irbc = 1, 1
      rbc => rbcs(irbc)

!       print *,rbc%Xc(1:2),Lb(1:2); stop

      nc(1:2) = rbc%Xc(1:2) - Lb(1:2)/2.
      nc = nc/SQRT(SUM(nc*nc))

      do ilon = 1, rbc%nlon
        do ilat = 1, rbc%nlat
          rbc%x(ilat, ilon, :) = rbc%x(ilat, ilon, :) + nc*0.0005
        end do
      end do

    end do ! irbc

  end subroutine AddR0Motion

!**********************************************************************
! Add repulsion force to prevent cells from self-overlapping
  subroutine AddR0GravityForce

    type(t_rbc), pointer :: rbc
    integer :: i, j, irbc, ilat, ilon
    integer :: ierr

    real(WP), dimension(3) :: nc
    real(WP), parameter    :: fo = 40000.

    if (rootworld) then
      print *, "ADDING R0 GRAVITY FORCE ", fo
    end if

    ! rbc => rbcs(1)
    ! do ilon = 1,rbc%nlon
    !       do ilat = 1,rbc%nlat

    !         rbc%f(ilat,ilon,3) = rbc%f(ilat,ilon,3) + fo*rbc%detj(ilat,ilon)*rbc%w(ilat)

    !   end do
    !end do

    nc(3) = 0.
    do irbc = 1, 1

      nc(1:2) = rbc%Xc(1:2) - Lb(1:2)/2.
      nc = nc/SQRT(SUM(nc*nc))

      rbc => rbcs(irbc)
      do ilon = 1, rbc%nlon
        do ilat = 1, rbc%nlat

          nc(1:2) = rbc%Xc(1:2) - Lb(1:2)/2.
          nc = nc/SQRT(SUM(nc*nc))
          !     rbc%f(ilat,ilon,:) = rbc%f(ilat,ilon,:) - nc*fo*rbc%detj(ilat,ilon)*rbc%w(ilat)
          rbc%x(ilat, ilon, :) = rbc%x(ilat, ilon, :) - nc*0.2
        end do
      end do

    end do ! irbc

    ! nc(3) = 0.
    ! do irbc = 1 , nrbc

    !    rbc => rbcs(irbc)
    !    do ilon = 1,rbc%nlon
    !       do ilat = 1,rbc%nlat

    !          nc(1:2) = rbc%x(ilat,ilon,1:2)-Lb(1:2)/2.
    !          nc = nc/SQRT(SUM(nc*nc))
    !          rbc%f(ilat,ilon,:) = rbc%f(ilat,ilon,:) - nc*fo*rbc%detj(ilat,ilon)*rbc%w(ilat)

    !       end do
    !    end do

    ! end do ! irbc

  end subroutine AddR0GravityForce

!**********************************************************************
! Add repulsion force to prevent cells from self-overlapping
  subroutine AddIntraCellRepulsionForce

    type(t_TargetList), pointer :: tlist
    type(t_SourceList), pointer :: slist
    type(t_rbc), pointer :: rbc
    integer :: cntDf, cntDf_Glb
    integer, allocatable :: indxDf(:)
    real(WP), allocatable :: df(:, :), dfGlb(:, :)
    real(WP) :: epsDist, epsDistRef, dist, distMin, distMin_Glb
    integer :: i, j, irbc, ilat, ilon, p, NN
    integer :: i1, i2, i3, j1, j2, j3
    real(WP) :: xi(3), a3i(3), thi, phii, xj(3), a3j(3), thj, phij
    real(WP) :: xx(3), rr, rrRef, FF, xxn, fi(3)
    integer :: ierr

    tlist => tlist_rbc
    slist => slist_rbc

    ! Allocate working arrays
    cntDf = count(tlist%active)
    allocate (indxDf(cntDf), df(cntDf, 3))
    allocate (dfGlb(tlist%npoint, 3))

    ! Initialize
    cntDf = 0
    distMin = huge(distMin)

    if (PhysEwald) then
      do i = 1, tlist%nPoint
        if (.not. tlist%active(i)) cycle

        irbc = tlist%indx(i, 0)
        ilat = tlist%indx(i, 1)
        ilon = tlist%indx(i, 2)

        rbc => rbcs(irbc)

        xi = rbc%x(ilat, ilon, :)
        a3i = rbc%a3(ilat, ilon, :)
        thi = rbc%th(ilat)
        phii = rbc%phi(ilon)

        ! Set distance threshold
        epsDist = 2*sqrt(rbc%area)/rbc%nlat
        epsDistRef = rbc%patch%radius

        !  Find the shortest distance of xi to the same surface
        !
        !  For a point xj on the same surface, let xx = xi - xj and xxn = -xx*a3j
        !  xj is not considered for distance calculation if any of the following happens
        !  -- The distance between xi and xj on the reference sphere is
        !     smaller than epsDistRef
        !  -- |xx| > 2*epsDist
        !  -- The angle between xx and -a3j is greater than 45 degree
        !
        !  Other wise, the distance between xi and xj is xxn

        ! Initialize
        dist = huge(dist)

        call HashTable_Index(slist%Nc, slist%iLbNc, xi, i1, i2, i3)
        do j1 = max(i1 - 1, 0), min(i1 + 1, slist%Nc(1) + 1)
        do j2 = max(i2 - 1, 0), min(i2 + 1, slist%Nc(2) + 1)
        do j3 = max(i3 - 1, 0), min(i3 + 1, slist%Nc(3) + 1)
          j = slist%hoc(j1, j2, j3)

          do while (j > 0)
            if (slist%indx(j, 0) .ne. irbc) goto 999

            ilat = slist%indx(j, 1)
            ilon = slist%indx(j, 2)

            xj = rbc%x(ilat, ilon, :)
            a3j = rbc%a3(ilat, ilon, :)
            thj = rbc%th(ilat)
            phij = rbc%phi(ilon)

            rrRef = DistOnSphere(thi, phii, thj, phij)
            if (rrRef < epsDistRef) goto 999

            xx = xi - xj
            rr = sqrt(sum(xx*xx))
            if (rr > epsDist) goto 999

            xxn = -dot_product(xx, a3j)
            if (xxn/rr < 0.707) goto 999

            dist = min(dist, rr)

999         j = slist%next(j)
          end do ! j
        end do ! j3
        end do ! j2
        end do ! j1

        if (dist < epsDist) then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Hard coded, future change needed
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          FF = ForceCoef*exp(1 - dist/epsDist)
          fi = FF*a3i

          cntDf = cntDf + 1
          indxDf(cntDf) = i
          df(cntDf, :) = fi
        end if

        distMin = min(rr, distMin)
      end do ! i
    end if

    ! Output information
    call MPI_AllReduce(distMin, distMin_Glb, 1, MPI_WP, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_AllReduce(cntDf, cntDf_Glb, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    if (rootWorld) then
      if (cntDf_Glb > 0) then
        write (*, '(A)') 'Intra-cellular repulsion:'
        write (*, '(A,ES12.2)') '  min dist = ', distMin_Glb
        write (*, '(A,I5,A)') '  repulsion forces added on ', cntDf_Glb, ' points'
      end if
    end if

    ! Synchronize
    if (cntDf_Glb > 0) then
      dfGlb = 0.0
      call CollectArray(cntDf, 3, indxDf, df, dfGlb, MPI_COMM_WORLD)

      p = 0
      do irbc = 1, nrbc
        rbc => rbcs(irbc)
        NN = rbc%nlat*rbc%nlon
        rbc%f = rbc%f + reshape(dfGlb(p + 1:p + NN, :), shape(rbc%f))
        p = p + NN
      end do ! irbc
    end if

    ! Deallocate working arrays
    deallocate (indxDf, df, dfGlb)

  end subroutine AddIntraCellRepulsionForce

!**********************************************************************
! Move cell mesh points when they are too close to another cell or wall
  subroutine InterCellRepulsion

    type(t_TargetList), pointer :: tlist
    type(t_Rbc), pointer :: rbc
    integer :: i, surfId_i
    real(WP) :: xi(3), dist, dist1, dist2, xj(3), x1(3), x2(3), xx(3), rr
    integer :: cntDX, cntDx_Glb
    integer, allocatable :: indxDx(:)
    real(WP), allocatable :: dx(:, :), dxGlb(:, :)
    real(WP) :: distMin, distMin_Glb

    real(WP) :: drig, dxrig(3)

    !    real(WP), parameter  :: viscRatThresh = 10.0
    integer :: irbc, NN, Ndx, p, ip
    integer :: ierr

    ! Update cell target and source list
    call TargetList_Update(tlist_rbc, rbcs)
    call SourceList_UpdateCoord(slist_rbc, rbcs)

    ! Allocate working arrays
    tlist => tlist_rbc
    cntDX = count(tlist%active)
    allocate (indxDx(cntDx), dx(cntDx, 3), dxGlb(tlist%npoint, 3))

    ! Initialize
    cntDx = 0
    distMin = huge(distMin)

    if (PhysEwald) then
      do i = 1, tlist%nPoint
        if (.not. tlist%active(i)) cycle

        xi = tlist%x(i, :)
        surfId_i = tlist%indx(i, 0)

        ! Set up distance threshold
        rbc => rbcs(surfId_i)

        call Closest_Neighbor_Cell(xi, surfId_i, epsDist, x1, dist1)
        call Closest_Neighbor_Wall(xi, surfId_i, epsDist, x2, dist2)

        if (dist1 < dist2) then
          xj = x1
          dist = dist1
        else
          xj = x2
          dist = dist2
        end if

        distMin = min(distMin, dist)

        if (dist < epsDist) then
          xx = xi - xj
          xx = xx - nint(xx*iLb)*Lb
          rr = sqrt(sum(xx*xx))

          cntDx = cntDx + 1
          indxDx(cntDx) = i
          dx(cntDx, :) = 0.5*xx*(epsDist - rr)/rr
        end if
      end do ! i
    end if

    ! Output information
    call MPI_AllReduce(distMin, distMin_Glb, 1, MPI_WP, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_AllReduce(cntDx, cntDx_Glb, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (rootWorld) then
      if (cntDx_Glb > 0) then
        write (*, '(A,F9.4)') 'Inter-cellular repulsion, epsDist = ', epsDist
        write (*, '(A,ES12.2)') '  min surface separation = ', distMin_Glb
        write (*, '(A,I5,A)') '  move ', cntDx_Glb, ' points'
        write (*, '(A,F6.2)') '  viscRatThresh = ', viscRatThresh
      end if
    end if

    ! Synchronize the changed coordinates
    if (cntDx_Glb > 0) then
      dxGlb = 0.0
      call CollectArray(cntDx, 3, indxDx, dx, dxGlb, MPI_COMM_WORLD)
!!$       if (rootWorld) print *,"MAXVAL dx",MAXVAL(dxGlb(:,3))
!!$       if (rootWorld) print *,"MINVAL dx",MINVAL(dxGlb(:,3))

      p = 0
      do irbc = 1, nrbc
        rbc => rbcs(irbc)
!!$          print *,"cell = ",irbc," MIN/MAX z= ",MINVAL(rbc%x(:,:,3)),MAXVAL(rbc%x(:,:,3))
        NN = rbc%nlat*rbc%nlon
        if (ABS(viscRat(rbc%celltype)) .lt. viscRatThresh) then

          rbc%x = rbc%x + reshape(dxGlb(p + 1:p + NN, :), shape(rbc%x))

        else if (rigidsep) then

          ! find max suggested displacement for cell i
!!$             drig = 0.; dxrig = 0.
!!$             do ip = 1,NN
!!$                if (SUM(dxGlb(p+ip,:)**2).gt.drig) then
!!$                   dxrig(:) = dxGlb(p+ip,:)
!!$                   drig = SUM(dxrig**2)
!!$                end if
!!$             end do

          dxrig = 0.
          Ndx = 0
          do ip = 1, NN
            if (SUM(dxGlb(p + ip, :)**2) .gt. 0.) then
              dxrig(:) = dxrig(:) + dxGlb(p + ip, :)
              Ndx = Ndx + 1
            end if
          end do
          if (Ndx .gt. 0) dxrig = dxrig/REAL(Ndx)

          ! apply that mx to all points on cell
          do i = 1, 3
            rbc%x(:, :, i) = rbc%x(:, :, i) + dxrig(i)
          end do

          if (Ndx .gt. 0 .and. rootWorld) print *, "RIGIDSEP AVERAGE CELL =", irbc, dxrig

        end if
!!$ print *,"AFTER MIN/MAX z= ",MINVAL(rbc%x(:,:,3)),MAXVAL(rbc%x(:,:,3))

        p = p + NN
      end do ! irbc

    end if

    ! Deallocate working arrays
    deallocate (indxDx, dx, dxGlb)

  end subroutine InterCellRepulsion

!**********************************************************************
! Move stiff-cell mesh points as rigid when they are too close to wall
  subroutine LeukWallRepulsion

    type(t_Rbc), pointer :: rbc
    integer :: i
    real(WP) :: xi(3), dist, xj(3), xx(3), rr
    real(WP) :: dx(3)
    real(WP) :: distMin, distMin_Glb

!    real(WP), parameter  :: viscRatThresh = 10.0
    integer :: irbc, ilat, ilon
    integer :: ierr

    if (PhysEwald) then

      do irbc = 1, nrbc
        rbc => rbcs(irbc)
        if (ABS(viscRat(rbc%celltype)) .lt. viscRatThresh) cycle

        dx = 0.
        distMin = huge(distMin)

        do ilon = 1, rbc%nlon
          do ilat = 1, rbc%nlat

            xi = rbc%x(ilat, ilon, :)
            call Closest_Neighbor_Wall(xi, irbc, epsDist, xj, dist)

            if (dist .lt. distMin) then
              distMin = dist
!                   print *,ilon,ilat,distMin
              if (dist .lt. epsDist) then

                xx = xi - xj
                xx = xx - nint(xx*iLb)*Lb
                rr = sqrt(sum(xx*xx))

                dx = 0.5*xx*(epsDist - rr)/rr

              end if
            end if

          end do
        end do

!          print "('DISTMIN1: ',4F10.4)",distMin,dx
        call MPI_AllReduce(distMin, distMin_Glb, 1, MPI_WP, MPI_MIN, MPI_COMM_WORLD, ierr)
        if (distMin .ne. distMin_Glb) dx = HUGE(dx)
        call MPI_AllReduce(dx, dx, 3, MPI_WP, MPI_MIN, MPI_COMM_WORLD, ierr)

        if (rootWorld .and. distMin .lt. epsDist) then
          write (*, '(A,F9.4)') 'LEUK-WALL repulsion, epsDist = ', epsDist
          write (*, '(A,ES12.2)') '  min surface separation = ', distMin_Glb
          write (*, '(A,F6.2)') '  viscRatThresh = ', viscRatThresh
          write (*, '(3F8.4)') dx
        end if

        do i = 1, 3
          rbc%x(:, :, i) = rbc%x(:, :, i) + dx(i)
        end do

      end do

    end if

  end subroutine LeukWallRepulsion

!**********************************************************************
! Find the cell that is closest to a given point
! Argument:
!  xi -- the point
!  surfId_i -- ID of the surface that xi(:) lies on
!  epsDist -- distance threshold
!  x0 -- projection of x0 on the closest neighbor surface
!  dist0 -- ||xi - x0||
  subroutine Closest_Neighbor_Cell(xi, surfId_i, epsDist, x0, dist0)
    real(WP) :: xi(3)
    integer :: surfId_i
    real(WP) :: epsDist
    real(WP) :: x0(3), dist0

    type(t_SourceList), pointer :: slist
    real(WP) :: xj(3), xx(3), rr, xtar(3)
    integer :: i, i1, i2, i3, j, j1, j2, j3
    integer :: surfId_j, irbc0, ilat0, ilon0
    type(t_Rbc), pointer :: rbc0
    real(WP) :: th0, phi0

    ! Initialize
    dist0 = huge(dist0)

    if (nrbc == 0) return   ! Need to actually have cells

    ! Loop over neighboring cells
    slist => slist_rbc
    call HashTable_Index(slist%Nc, slist%iLbNc, xi, i1, i2, i3)

    do j1 = max(i1 - 1, 0), min(i1 + 1, slist%Nc(1) + 1)
    do j2 = max(i2 - 1, 0), min(i2 + 1, slist%Nc(2) + 1)
    do j3 = max(i3 - 1, 0), min(i3 + 1, slist%Nc(3) + 1)
      j = slist%hoc(j1, j2, j3)

      do while (j > 0)
        xj = slist%x(j, :)
        surfId_j = slist%indx(j, 0)

        if (surfId_j == surfId_i) goto 999

        xx = xi - xj
        xx = xx - nint(xx*iLb)*Lb
        rr = sqrt(sum(xx*xx))

        if (rr < dist0) then
          irbc0 = surfId_j - rbcs(1)%id + 1
          ilat0 = slist%indx(j, 1)
          ilon0 = slist%indx(j, 2)
          dist0 = rr
        end if

999     j = slist%next(j)
      end do ! while
    end do ! j3
    end do ! j2
    end do ! j1

    if (dist0 <= 2*epsDist) then
      rbc0 => rbcs(irbc0)

      x0 = rbc0%x(ilat0, ilon0, :)
      th0 = rbc0%th(ilat0)
      phi0 = rbc0%phi(ilon0)

      xx = xi - x0
      xx = xx - nint(xx*iLb)*Lb
      xtar = x0 + xx
      call Spline_FindProjection(rbc0%spln_x, xtar, th0, phi0, x0)

      xx = xtar - x0
      dist0 = sqrt(sum(xx*xx))
    end if

  end subroutine Closest_Neighbor_Cell

!**********************************************************************
! Find the wall that is closest to a point
! Arguments:
!  xi -- the point
!  surfId_i -- ID of the surface that xi lies on
!  epsDist -- distance threshold
!  x0 -- the projection of xi on the closest neighboring wall
!  dist0 -- ||xi -x0||
  subroutine Closest_Neighbor_Wall(xi, surfId_i, epsDist, x0, dist0)
    real(WP) :: xi(3)
    integer :: surfId_i
    real(WP) :: epsDist
    real(WP) :: x0(3), dist0

    type(t_SourceList), pointer :: slist
    type(t_Wall), pointer :: wall
    integer :: i1, i2, i3, j1, j2, j3, j
    real(WP) :: xj(3), xx(3), rr
    integer :: surfId_j, iwall, iele, l
    real(WP) :: xele(3, 3), s, t, xtar(3)

    ! Initialize
    dist0 = huge(dist0)

    if (nwall == 0) return  ! Need to have a wall

    ! Loop over all neighboring elements
    slist => slist_wall
    call HashTable_Index(slist%Nc, slist%iLbNc, xi, i1, i2, i3)

    do j1 = max(i1 - 1, 0), min(i1 + 1, slist%Nc(1) + 1)
    do j2 = max(i2 - 1, 0), min(i2 + 1, slist%Nc(2) + 1)
    do j3 = max(i3 - 1, 0), min(i3 + 1, slist%Nc(3) + 1)
      j = slist%hoc(j1, j2, j3)

      do while (j > 0)
        surfId_j = slist%indx(j, 0)
        if (surfId_j == surfId_i) goto 999

        iwall = surfId_j - walls(1)%id + 1
        wall => walls(iwall)
        iele = slist%indx(j, 1)

        do l = 1, 3
          xele(l, :) = wall%x(wall%e2v(iele, l), :)
        end do ! l

        ! Translate xi as close to the triangle as possible
        xx = xi - xele(1, :)
        xx = xx - nint(xx*iLb)*Lb
        xtar = xele(1, :) + xx

        rr = MinDistToTri(xtar, xele, x0=xj)

        if (rr < dist0) then
          dist0 = rr
          x0 = xj
        end if

999     j = slist%next(j)
      end do ! while
    end do ! j3
    end do ! j2
    end do ! j1

  end subroutine Closest_Neighbor_Wall

!**********************************************************************

end module ModRepulsion
