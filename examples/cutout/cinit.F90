program randomized_cell_gen

  use ModDataTypes
  use ModDataStruct
  use ModRbc
  use ModWall
  use ModConf
  use ModData
  use ModIO
  use ModBasicMath

  integer, parameter :: ranseed = 5607

  !initial condition setup parameters
  real(WP), parameter :: hcrit = 0.05
  real(WP) :: tuber = 3
  real(WP), parameter :: tubelen = 20

  integer :: nrbcMax ! how many cells

  type(t_Wall), pointer :: wall
  character(CHRLEN) :: fn
  integer :: i
  real(WP) :: actlen
  real(WP) :: clockBgn, clockEnd
  integer :: nodeNum, numNodes
  real(WP) :: radEqv = 1.

  call InitMPI
  call MPI_Comm_Rank(MPI_COMM_WORLD, nodeNum, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, numNodes, ierr)
  nodeNum = nodeNum + 1

  ! calculate number of cells for the defined hematocrit, assuming all blood cells are healthy RBCs for volume
  ! hematocrit = 4 * nrbc / (3 * tube_radius^2 * tube_length)
  ! nrbcMax = ((3*(tubelen*tuber**2*hcrit))/4)
  nrbcMax = 38

  if (rootWorld) write (*, *) "Num RBCs in simulation is ", nrbcMax

  !set other initialization params
  vBkg(1:2) = 0.; vBkg(3) = 8.
  Nt = 0; time = 0.

  !Create wall
  nwall = 1
  allocate (walls(nwall))
  wall => walls(1)
  call ReadWallMesh('Input/cutout.e', wall)
  wall%f = 0.

  if (rootWorld) write (*, *) walls(1)%nvert
  actlen = tubelen

  !recenter walls to be strictly positive
  call recenterWalls

  !set periodic boundary box based on tube shape
  Lb(1) = maxval(wall%x(:, 1)) + 0.5
  Lb(2) = maxval(wall%x(:, 2)) + 0.5
  Lb(3) = tubelen

  !Write Walls Out to Wall TecPlot File
  write (fn, FMT=fn_FMT) 'D/', 'wall', 0, '.dat'
  call WriteManyWalls(fn, nwall, walls)

  !for each cell to add
  allocate (rbcs(nrbcMax))
  nrbc = 0

  do i = 1, nrbcMax
    if (rootWorld) write (*, *) "Adding Cell #", nrbc + 1

    clockBgn = MPI_Wtime()
    call place_cell(1, i)
    clockEnd = MPI_Wtime()

    nrbc = nrbc + 1
    if (rootWorld) write (*, *) "Cell #", nrbc, " Placed; Time Cost = ", clockEnd - clockBgn

    ! Write Cells Out to Cell TecPlot File
    if (modulo(i, 10) .eq. 0) then
      write (fn, FMT=fn_FMT) 'D/', 'x', 0, '.dat'
      call WriteManyRBCs(fn, nrbc, rbcs)

      ! Write restart blob
      write (fn, FMT=fn_FMT) 'D/', 'restart', 0, '.dat'
      call WriteRestart(fn, Nt0, time)
      fn = 'D/restart.LATEST.dat'
      call WriteRestart(fn, Nt0, time)
    end if
  end do

  write (fn, FMT=fn_FMT) 'D/', 'x', 0, '.dat'
  call WriteManyRBCs(fn, nrbc, rbcs)
  ! Write restart blob
  write (fn, FMT=fn_FMT) 'D/', 'restart', 0, '.dat'
  call WriteRestart(fn, Nt0, time)
  fn = 'D/restart.LATEST.dat'
  call WriteRestart(fn, Nt0, time)

  if (nrbcMax .ne. 0) deallocate (rbcs)
  if (nwall .ne. 0) deallocate (walls)
  call FinalizeMPI
  stop

contains

! offset everything to be positive
  subroutine recenterWalls
    integer :: i
    type(t_wall), pointer :: wall
    real(WP) :: offset(3)

    wall => walls(1)

    offset(1) = -1*minval(wall%x(:, 1))
    offset(2) = -1*minval(wall%x(:, 2))
    offset(3) = -1*minval(wall%x(:, 3))

    ! shift all walls
    do iwall = 1, nwall
      wall => walls(iwall)
      do i = 1, 3
        wall%x(:, i) = wall%x(:, i) + offset(i)
      end do
    end do

  end subroutine recenterWalls

! find an open spot in the simulation to place a cell
  subroutine place_cell(celltype, num)

    type(t_Rbc) :: newcell
    type(t_Rbc), pointer :: cell
    real(WP) :: tmp_xc(3)

    integer :: celli, ii, ierr
    integer :: nlat0
    logical :: place_success_loc, place_success_glb

    integer :: celltype, num
    real(WP) :: sphere_center(3), sphere_rad, d = 100.

    nlat0 = 12

    ! create template newcell
    select case (celltype)
    case (1)
      call RBC_Create(newcell, nlat0)
      call RBC_MakeBiconcave(newcell, radEqv)
    case (2)
      call RBC_Create(newcell, nlat0)
      call RBC_MakeLeukocyte(newcell, radEqv)
    case (3)
      call RBC_Create(newcell, nlat0)
      call RBC_MakeSphere(newcell, radEqv)
    case default
      stop "bad cellcase"
    end select
    newcell%celltype = celltype

    place_success_glb = .false.
    tmp_xc = 0.

    !choose a location & orientation for newcell to fit into the simulation
    do while (.not. place_success_glb)

      !reset cell position
      do ii = 1, 3
        newcell%x(:, :, ii) = newcell%x(:, :, ii) - tmp_xc(ii)
      end do

      call rotate_cell(newcell)

      ! randomly select a tmp_xc
      ! place first 22 cells in sphere
      if (num .le. 22) then
        sphere_center = (/4.5, 5., 9./)
        sphere_rad = 4.5
        do while (d .gt. sphere_rad*sphere_rad)
          tmp_xc(1) = (RandomNumber(ranseed)*sphere_rad*2) - sphere_rad
          tmp_xc(2) = (RandomNumber(ranseed)*sphere_rad*2) - sphere_rad
          tmp_xc(3) = (RandomNumber(ranseed)*sphere_rad*2) - sphere_rad

          d = sum(tmp_xc*tmp_xc)
        end do
        d = 100.
        tmp_xc = tmp_xc + sphere_center
      ! place rest of the cells in cylinder
      else
        tmp_xc(2) = RandomNumber(ranseed)*2*PI
        tmp_xc(3) = sqrt(RandomNumber(ranseed))*(tuber)
        ! shift by 3 + 1.5 to account for cylinder with rad 3 being centered 4.5 units away in x dir
        tmp_xc(1) = tmp_xc(3)*cos(tmp_xc(2)) + 4.5
        tmp_xc(2) = tmp_xc(3)*sin(tmp_xc(2)) + 3
        tmp_xc(3) = RandomNumber(ranseed)*tubelen
      end if

      ! shift newcell by the xc
      do ii = 1, 3
        newcell%x(:, :, ii) = newcell%x(:, :, ii) + tmp_xc(ii)
      end do

      place_success_glb = .true.
      newcell%x(:, :, 1) = newcell%x(:, :, 1) - tuber
      newcell%x(:, :, 2) = newcell%x(:, :, 2) - tuber
      ! if (rootWorld) print *, "TUBELEN CHECK"
      if (any(newcell%x(:, :, 3) .ge. tubelen .or. newcell%x(:, :, 3) .lt. 0)) then
        place_success_loc = .false.
        place_success_glb = .false.
      end if
      newcell%x(:, :, 1) = newcell%x(:, :, 1) + tuber
      newcell%x(:, :, 2) = newcell%x(:, :, 2) + tuber
      if (.not. place_success_glb) cycle

      !check if cell collides with wall
      place_success_loc = .not. check_wall_collision(newcell)
      place_success_glb = .true.
      call MPI_AllReduce(place_success_loc, place_success_glb, 1, MPI_Logical, MPI_LAND, MPI_COMM_WORLD, ierr)
      if (.not. place_success_glb) cycle

      ! if no wall collision, check if cell collides with other cells
      celli = 1
      do while (place_success_loc .and. celli .le. nrbc)
        cell => rbcs(celli)
        place_success_loc = .not. check_cell_collision(cell, newcell)
        celli = celli + 1
      end do

      place_success_glb = .true.
      call MPI_AllReduce(place_success_loc, place_success_glb, 1, MPI_Logical, MPI_LAND, MPI_COMM_WORLD, ierr)

    end do ! while not place_success

    ! the cell placement location was successful
    rbcs(nrbc + 1) = newcell
  end subroutine place_cell

! helper for place cell, chooses a random orientation and rotates the cell
  subroutine rotate_cell(cell)
    type(t_Rbc) :: cell

    real(WP) :: v1, v2, vsq
    real(WP) :: zv(3), rotmat(3, 3)

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
    rotmat = RotateMatrix(zv)

    !rotate cell with rotmat
    !each point: x = Rx
    forall (i=1:cell%nlat, j=1:cell%nlon)
      ! reshape is unnecessary
      cell%x(i, j, :) = matmul(rotmat, cell%x(i, j, :))
    end forall

  end subroutine rotate_cell

! helper for place_cell, checks if cell1 intersects/collides with wall
  logical function check_wall_collision(cell)
    type(t_Rbc) :: cell

    integer :: i, j, i2, j2
    real(WP) :: c1p(3), c2p(3), sp(3)

    real(WP) :: threshold
    threshold = 0.3

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

          if (norm2(c1p - sp) .lt. threshold) then
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

    check_cell_collision = .false.

    !first check the centers; if distance > 4 then we are sure that they don't collide
    if (NORM2(cell1%xc - cell2%xc) .ge. 4) return

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
      do p_it = nodeNum, cell2%nlat*cell2%nlon, numNodes
        ! do i2 = 1, cell2%nlat
        ! do j2 = 1, cell2%nlon
        i2 = (p_it/cell2%nlon) + 1
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

end program randomized_cell_gen
