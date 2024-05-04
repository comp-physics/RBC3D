program InitCond
    use ModDataTypes
    use ModDataStruct
    use ModWall
    use ModConf
    use ModData
    use ModIO
    use ModBasicMath
    use ModRbc

    integer, parameter :: ranseed = 112

    ! initial condition setup parameters
    real(WP), parameter :: hematocrit = 0.17
    real(WP) :: tuber = 4
    real(WP), parameter :: tubelen = 20

    integer :: nrbcMax ! how many cells

    type(t_Rbc), pointer :: rbc
    type(t_Wall), pointer :: wall
    character(CHRLEN) :: fn
    integer :: zmax
    integer :: i
    real(WP) :: th, actlen
    real(WP) :: clockBgn, clockEnd, totTime
    real(WP) :: offset(3)

    call InitMPI

    ! calculate number of cells for the defined hematocrit, assuming all blood cells are regularly shaped RBCs for volume
    ! hematocrit = 4 * nrbc / (3 * tube_radius^2 * tube_length)
    nrbcMax = int((3*(tubelen*tuber**2*hematocrit)) / 4)

    write (*, *) "Num RBCs in simulation is ", nrbcMax

    ! set other initialization params
    vbkg(1:2) = 0.
    vbkg(3) = 10.
    Nt = 0
    time = 0

    ! create wall
    nwall = 1
    allocate (walls(nwall))
    wall => walls(1)
    call ReadWallMesh("Input/new_cyl_D6_L13_33_hires.e", wall)
    actlen = 13.33
    wall%f = 0.
    do i = 1, wall%nvert
        th = ATAN2(wall%x(i, 1), wall%x(i, 2))
        wall%x(i, 1) = (tuber)*COS(th)
        wall%x(i, 2) = (tuber)*SIN(th)   
        wall%x(i, 3) = tubelen/actlen*wall%x(i, 3)
    end do

    ! recenter walls to be strictly positive
    call recenterWalls

    Lb(1) = maxval(wall%x(:, 1)) + 0.5
    Lb(2) = maxval(wall%x(:, 2)) + 0.5
    Lb(3) = tubelen

    allocate(rbcs(nrbcMax))

    nrbc = 0
    do i = 1, nrbcMax
        if (rootWorld) write(*, *) "Adding Cell #", nrbc + 1
        clockBgn = MPI_Wtime()
        if (i .le. 10) then
          call place_cell(2)
        else
          call place_cell(1)
        end if
        clockEnd = MPI_Wtime()

        nrbc = nrbc + 1
        totTime = clockEnd - clockBgn
        if (rootWorld) write (*, *) "Cell #", nrbc, "Placed; Time Cost = ", totTime
    end do

    ! if (rootWorld)
    !Write Cells Out to Cell TecPlot File
    write (fn, FMT=fn_FMT) 'D/', 'x', 0, '.dat'
    call WriteManyRBCs(fn, nrbc, rbcs)
    
    !Write Walls Out to Wall TecPlot File
    write (fn, FMT=fn_FMT) 'D/', 'wall', 0, '.dat'
    call WriteManyWalls(fn, nwall, walls)

    !Write restart blob
    write (fn, FMT=fn_FMT) 'D/', 'restart', 0, '.dat'
    call WriteRestart(fn, Nt0, time)
    fn = 'D/restart.LATEST.dat'
    call WriteRestart(fn, Nt0, time)

    if (nrbcMax .ne. 0) deallocate (rbcs)
    deallocate (walls)
    call FinalizeMPI
    stop


contains
    ! offset everything to be positive
    subroutine recenterWalls
        integer :: i, ii
        type(t_wall), pointer :: wall
        real(WP) :: offset(3)

        wall => walls(1)

        offset(1) = -1 * minval(wall%x(:, 1))
        offset(2) = -1 * minval(wall%x(:, 2))
        offset(3) = -1 * minval(wall%x(:, 3))

        ! shift wall
        do ii = 1, 3
            wall%x(:, ii) = wall%x(:, ii) + offset(ii)
        end do
    end subroutine recenterWalls

    subroutine choose_point(wall, pt, type)
        real(WP) :: pt(3)
        type(t_wall) :: wall

        real(WP) :: rad, len, maxN, minN
        real(WP) :: thresh
        integer :: i, type

        len = RandomNumber(ranseed)*tubelen
        thresh = 0.2
        maxN = -1E9
        minN = 1E9

        ! figure out the radius based on the location in the tube
        do i = 1, wall%nvert
            if (abs(wall%x(i, 3) - len) .le. thresh) then
                maxN = maxval((/ maxN, wall%x(i, 1), wall%x(i, 2) /))
                minN = minval((/ minN, wall%x(i, 1), wall%x(i, 2)/))
            end if
        end do

        rad = (maxN - minN) / 2

        pt(2) = RandomNumber(ranseed)*2*PI
        if (type .eq. 4) then
          print *, "placed celltype 2 closer to center"
          ! pt(3) = sqrt(RandomNumber(ranseed))*(rad / 3.0)
          pt(3) = .1
        else
          pt(3) = sqrt(RandomNumber(ranseed))*(rad)
        end if
        pt(1) = pt(3)*cos(pt(2)) + rad + minN
        pt(2) = pt(3)*sin(pt(2)) + rad + minN

        pt(3) = len
    end subroutine choose_point

    ! subroutine place_cell_hard(celltype)
    !   type(t_Rbc) :: newcell
    !   type(t_Rbc), pointer :: cell
    !   real(WP) :: xc(3)

    !   integer :: celli, ii, ierr
    !   integer :: nlat0, dealias
    !   logical :: place_success_loc, place_success_glb

    !   integer :: celltype

    !   nlat0 = 12
    !   dealias = 3
    !   xc(1:2) = tuber
    !   xc(3) = 5

    !   ! do i = 1, num
    !   ! end do

    !   !create template newcell
    !   select case (celltype)
    !   case (1)
    !       call RBC_Create(newcell, nlat0, dealias)
    !       call RBC_MakeBiconcave(newcell, 1., xc)
    !   case (2)
    !       call RBC_Create(newcell, nlat0, dealias)
    !       call RBC_MakeLeukocyte(newcell, 1., xc)
    !   case (3)
    !       call ImportReadRBC('Input/SickleCell.dat', newcell)
    !   case default
    !       stop "bad cellcase"
    !   end select
    !   newcell%celltype = celltype

    !   rbcs(nrbc + 1) = newcell
    ! end subroutine place_cell_hard

    ! find an open spot in the simulation to place a cell
    subroutine place_cell(celltype)
        type(t_Rbc) :: newcell
        type(t_Rbc), pointer :: cell
        real(WP) :: tmp_xc(3)

        integer :: celli, ii, ierr
        integer :: nlat0, dealias
        logical :: place_success_loc, place_success_glb

        integer :: celltype

        nlat0 = 12
        dealias = 3

        !create template newcell
        select case (celltype)
        case (1)
            call RBC_Create(newcell, nlat0, dealias)
            call RBC_MakeBiconcave(newcell, 1.)
        case (2)
            call RBC_Create(newcell, nlat0, dealias)
            call RBC_MakePlatelet(newcell, 1.)
        case (3)
            call ImportReadRBC('Input/SickleCell.dat', newcell)
        case default
            stop "bad cellcase"
        end select
        newcell%celltype = celltype

        place_success_glb = .false.
        tmp_xc = 0.

        !choose a location & orientation for newcell to fit into the simulation
        do while (.not. place_success_glb)

          !reset cell position
          do ii = 1,3
              newcell%x(:,:,ii) = newcell%x(:,:,ii) - tmp_xc(ii)
          end do

          !randomly rotate cell
          call rotate_cell(newcell)

          !randomly select a tmp_xc
          call choose_point(walls(1), tmp_xc, celltype)

          !shift newcell by the xc
          do ii = 1, 3
              newcell%x(:, :, ii) = newcell%x(:, :, ii) + tmp_xc(ii)
          end do

          place_success_glb = .true.
          newcell%x(:,:,1) = newcell%x(:,:,1) - tuber
          newcell%x(:,:,2) = newcell%x(:,:,2) - tuber
          if (any(newcell%x(:, :, 3) .ge. tubelen .or. &
                  newcell%x(:,:,3) .lt. 0 .or. &
                  (newcell%x(:, :, 1)**2 + newcell%x(:, :, 2)**2) .ge. (tuber+1)**2)) then
              place_success_loc = .false.
              place_success_glb = .false.
          end if
          newcell%x(:,:,1) = newcell%x(:,:,1) + tuber
          newcell%x(:,:,2) = newcell%x(:,:,2) + tuber
          if (.not. place_success_glb) cycle

          !check if cell collides with wall
          place_success_loc = .not. check_wall_collision(newcell)

          !if no wall collision, check if cell collides with other cells
          celli = 1
          do while ( place_success_loc .and. celli .le. nrbc )
              cell => rbcs(celli)
              place_success_loc = .not.check_cell_collision(cell, newcell)
              celli = celli + 1
          end do 

          place_success_glb = .true.
          call MPI_AllReduce(place_success_loc, place_success_glb, 1, MPI_Logical, MPI_LAND, MPI_COMM_WORLD, ierr)
        
        end do !while not place_success

        !the cell placement location was successful
        rbcs(nrbc + 1) = newcell
    end subroutine place_cell

  !helper for place cell, chooses a random orientation and rotates the cell
  subroutine rotate_cell(cell)
    type(t_Rbc) :: cell

    real(WP) :: v1, v2, vsq
    real(WP) :: zv(3), rc(3), rotmat(3, 3)

    integer :: i, j

    !randomly generate unit vector
    vsq = 10
    do while (vsq .ge. 1)
      v1 = RandomNumber(ranseed)*2 - 1
      v2 = RandomNumber(ranseed)*2 - 1
      vsq = v1**2 + v2**2
    end do
    zv(1) = v1*2*sqrt(1 - vsq)
    zv(2) = v2*2*sqrt(1 - vsq)
    zv(3) = 1 - (2*vsq)

    !generate rotation matrix
    rotmat = RotateMatrix(zv) !(/ 0., 1., 0./))

    !rotate cell with rotmat
    !each point: x = Rx
    forall (i=1:cell%nlat, j=1:cell%nlon)
      cell%x(i, j, :) = reshape(matmul(rotmat, cell%x(i, j, :)), (/3/))
    end forall

  end subroutine rotate_cell

  !helper for place_cell, checks if cell1 intersects/collides with wall
  logical function check_wall_collision(cell)
    type(t_Rbc) :: cell

    integer :: i, j, i2, j2
    real(WP) :: c1p(3), c2p(3), sp(3)

    integer :: nodeNum, numNodes
    real(WP) :: threshold
    threshold = 0.1
    call MPI_Comm_Rank(MPI_COMM_WORLD, nodeNum, ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, numNodes, ierr)
    nodeNum = nodeNum + 1

    check_wall_collision = .false.

    !each point in cell1
    do i = 1, cell%nlat
    do j = 1, cell%nlon

      !define 2 diagonal points c1p & c2p for collision detection
      c1p = cell%x(i, j, :)
      c2p = cell%x(modulo(i, cell%nlat) + 1, modulo(j, cell%nlon) + 1, :)

      !each wall
      do i2 = 1, nwall
        !each point in the wall
        do j2 = nodeNum, walls(i2)%nvert, numNodes

          sp = walls(i2)%x(j2, :)

          if (VecNorm(c1p - sp) .lt. threshold) then
            check_wall_collision = .true.
            return
          end if
          !check if the wall point lies inside the cube delineated by c1p & c2p
          if ( &
            ((c1p(1) .le. sp(1) .and. sp(1) .le. c2p(1)) .or. (c2p(1) .le. sp(1) .and. sp(1) .le. c1p(1))) .and. &  !x direction overlap
            ((c1p(2) .le. sp(2) .and. sp(2) .le. c2p(2)) .or. (c2p(2) .le. sp(2) .and. sp(2) .le. c1p(2))) .and. &  !y direction overlap
            ((c1p(3) .le. sp(3) .and. sp(3) .le. c2p(3)) .or. (c2p(3) .le. sp(3) .and. sp(3) .le. c1p(3))) &        !z direction overlap
            ) then
            ! since we found a collision, return true
            check_wall_collision = .true.
            return
          end if

        end do !j2
      end do !i2

    end do !j
    end do !i
  end function check_wall_collision

  ! helper for place_cell, checks if cell1 and cell2 intersect/collide
  ! creates miniature bounding boxes along cell meshes to do AABB-AABB collision detection
  ! should not have any false-negatives (return false even if the cells are intersecting)
  ! may have false-positives (might return true if cells aren't intersecting, but are very close)
  logical function check_cell_collision(cell1, cell2)
    type(t_Rbc) :: cell1, cell2

    integer :: i, j, i2, j2, ii
    integer :: p_it

    real(WP) :: p1(3), p2(3)
    real(WP) :: b1(3, 2), b2(3, 2)

    integer :: numNodes, nodeNum

    check_cell_collision = .false.

    !first check the centers; if distance > 4 then we are sure that they don't collide
    if (VecNorm(cell1%xc - cell2%xc) .ge. 4) return

    call MPI_Comm_Rank(MPI_COMM_WORLD, nodeNum, ierr)
    call MPI_Comm_Size(MPI_COMM_WORLD, numNodes, ierr)
    nodeNum = nodeNum + 1

    !each point in cell1
    do i = 1, cell1%nlat
    do j = 1, cell1%nlon

      ! create first AABB from patch of cell1
      p1 = cell1%x(i, j, :)
      p2 = cell1%x(modulo(i, cell1%nlat) + 1, modulo(j, cell1%nlon) + 1, :)

      do ii = 1, 3
        b1(ii, 1) = min(p1(ii), p2(ii))
        b1(ii, 2) = max(p1(ii), p2(ii))
      end do

      !each point in cell2
      do p_it = nodeNum, cell2%nlat * cell2%nlon, numNodes
      ! do i2 = 1, cell2%nlat
      ! do j2 = 1, cell2%nlon
        i2 = (p_it / cell2%nlon) + 1
        j2 = modulo(p_it, cell2%nlon) + 1

        !create second AABB from patch of cell2
        p1 = cell2%x(i2, j2, :)
        p2 = cell2%x(modulo(i2, cell2%nlat) + 1, modulo(j2, cell2%nlon) + 1, :)
        do ii = 1, 3
          b2(ii, 1) = min(p1(ii), p2(ii))
          b2(ii, 2) = max(p1(ii), p2(ii))
        end do

        !check for AABB collision, we find AABB intersection iff X, Y, and Z intervals all overlap
        if ( &
          b1(1, 1) .le. b2(1, 2) &
          .and. b1(1, 2) .ge. b2(1, 1) &
          .and. b1(2, 1) .le. b2(2, 2) &
          .and. b1(2, 2) .ge. b2(2, 1) &
          .and. b1(3, 1) .le. b2(3, 2) &
          .and. b1(3, 2) .ge. b2(3, 1) &
          ) then
            check_cell_collision = .true.
            return
        end if

      ! end do !j2
      ! end do !i2
      end do !p_it
    end do !j
    end do !i
  end function check_cell_collision

end program InitCond