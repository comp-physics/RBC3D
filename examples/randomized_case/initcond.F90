program randomized_cell_gen

  use ModDataTypes
  use ModDataStruct
  use ModRbc
  use ModWall
  use ModConf
  use ModData
  use ModIO
  use ModBasicMath

  integer, parameter :: ranseed = 161269

  !initial condition setup parameters
  real(WP), parameter :: hematocrit = 0.10
  real(WP), parameter :: tuber = 7
  real(WP), parameter :: tubelen = 20

  integer :: nrbcMax ! how many cells

  type(t_Rbc), pointer :: rbc
  type(t_Wall), pointer :: wall
  character(CHRLEN) :: fn
  integer :: zmax
  integer :: i
  real(WP) :: th, actlen
  real(WP) :: clockBgn, clockEnd

  integer :: numNodes, nodeNum

  call InitMPI
  call MPI_Comm_Rank(MPI_COMM_WORLD, nodeNum, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, numNodes, ierr)
  nodeNum = nodeNum + 1

  !calculate number of cells for the defined hematocrit, assuming all blood cells are healthy RBCs for volume
  !hematocrit = 4 * nrbc / (3 * tube_radius^2 * tube_length)
  nrbcMax = ((3*(tubelen*tuber**2*hematocrit))/4)

  if (rootWorld) write (*, *) "Num RBCs in simulation is ", nrbcMax

  !set periodic boundary box based on tube shape
  Lb(1) = tuber*2 + 0.5
  Lb(2) = tuber*2 + 0.5
  Lb(3) = tubelen

  !set other initialization params
  vBkg(1:2) = 0.; vBkg(3) = 8.
  Nt = 0; time = 0.

  !Create wall
  nwall = 1
  allocate (walls(nwall))
  wall => walls(1)
  call ReadWallMesh('Input/new_cyl_D6_L13_33_hires.e', wall)
  actlen = 13.33

  wall%f = 0.
  do i = 1, wall%nvert
    th = ATAN2(wall%x(i, 1), wall%x(i, 2))
    wall%x(i, 1) = tuber*COS(th)    !!!!!!!!!
    wall%x(i, 2) = tuber*SIN(th)    !!!!!!!!!
    wall%x(i, 3) = tubelen/actlen*wall%x(i, 3)
  end do

  !for each cell to add
  allocate (rbcs(nrbcMax))
  nrbc = 0
  do i = 1, nrbcMax

    clockBgn = MPI_Wtime()
    call place_cell(1)
    clockEnd = MPI_Wtime()

    nrbc = nrbc + 1
    if (rootWorld) write (*, *) "Cell #", nrbc, " Placed; Time Cost = ", clockEnd - clockBgn
  end do

  !recenter everything to be all positive
  call recenter

  !Write Cells Out to Cell TecPlot File
  write (fn, FMT=fn_FMT) 'D/', 'x', 0, '.dat'
  call WriteManyRBCs(fn, nrbc, rbcs)

  !Uncomment these lines to write cells to separate files by their celltype
  !write(fn, FMT=fn_FMT) 'D/', '1x', 0, '.dat'
  !call WriteManyRBCsByType(fn, nrbc, rbcs, 1)
  !write(fn, FMT=fn_FMT) 'D/', '2x', 0, '.dat'
  !call WriteManyRBCsByType(fn, nrbc, rbcs, 2)
  !write(fn, FMT=fn_FMT) 'D/', '3x', 0, '.dat'
  !call WriteManyRBCsByType(fn, nrbc, rbcs, 3)

  !Write Walls Out to Wall TecPlot File
  write (fn, FMT=fn_FMT) 'D/', 'wall', 0, '.dat'
  call WriteManyWalls(fn, nwall, walls)
  !Write restart blob
  write (fn, FMT=fn_FMT) 'D/', 'restart', 0, '.dat'
  call WriteRestart(fn, Nt0, time)
  fn = 'D/restart.LATEST.dat'
  call WriteRestart(fn, Nt0, time)

  deallocate (rbcs)
  deallocate (walls)
  call FinalizeMPI
  stop

contains

!offset everything to be positive
  subroutine recenter
    integer :: i, ii
    type(t_rbc), pointer :: cell
    type(t_wall), pointer :: wall
    real(WP) :: offset(3)

    offset(1) = tuber
    offset(2) = tuber
    offset(3) = tubelen/2

    !shift wall
    do i = 1, nwall
      wall => walls(i)

      do ii = 1, 3
        wall%x(:, ii) = wall%x(:, ii) + offset(ii)
      end do
    end do

    !shift cells
    do i = 1, nrbc
      cell => rbcs(i)

      do ii = 1, 3
        cell%x(:, :, ii) = cell%x(:, :, ii) + offset(ii)
      end do
    end do

  end subroutine recenter

!find an open spot in the simulation to place a cell
  subroutine place_cell(celltype)

    type(t_Rbc) :: newcell
    type(t_Rbc), pointer :: cell
    real(WP) :: tmp_xc(3)

    integer :: celli, ii, ierr
    integer :: nlat0
    logical :: place_success_loc, place_success_glb

    integer :: celltype

    nlat0 = 12

    !create template newcell
    select case (celltype)
    case (1)
      call RBC_Create(newcell, nlat0)
      call RBC_MakeBiconcave(newcell, 1.)
    case (2)
      call RBC_Create(newcell, nlat0)
      call RBC_MakeLeukocyte(newcell, 1.)
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
      do ii = 1, 3
        newcell%x(:, :, ii) = newcell%x(:, :, ii) - tmp_xc(ii)
      end do

      !randomly rotate cell
      call rotate_cell(newcell)

      !randomly select a tmp_xc
      tmp_xc(2) = RandomNumber(ranseed)*2*PI
      tmp_xc(3) = sqrt(RandomNumber(ranseed))*(tuber)
      tmp_xc(1) = tmp_xc(3)*cos(tmp_xc(2))
      tmp_xc(2) = tmp_xc(3)*sin(tmp_xc(2))
      tmp_xc(3) = RandomNumber(ranseed)*tubelen - tubelen/2

      !shift newcell by the xc
      do ii = 1, 3
        newcell%x(:, :, ii) = newcell%x(:, :, ii) + tmp_xc(ii)
      end do

      !check if newcell collides with wall, assume cylinder wall here
      if (any(abs(newcell%x(:, :, 3)) .ge. tubelen/2 .or. &
              newcell%x(:, :, 1)**2 + newcell%x(:, :, 2)**2 .ge. tuber**2)) then

        !placement collides with wall, try again
        cycle
      end if

      !if your wall isn't a cylinder, you can do this instead
      ! place_success_loc = .not. check_wall_collision(newcell)
      ! place_success_glb = .true.
      ! call MPI_AllReduce(place_success_loc, place_success_glb, 1, MPI_Logical, MPI_LAND, MPI_COMM_WORLD, ierr)
      ! if (.not. place_success_glb) cycle

      place_success_loc = .true.
      celli = 1
      do while (place_success_loc .and. celli .le. nrbc)
        cell => rbcs(celli)
        place_success_loc = .not. check_cell_collision(cell, newcell)
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
    real(WP), parameter :: threshold = 0.2 !'magic' adjustable number cell-to-wall distance tolerance

    integer :: i, j, i2, j2
    real(WP) :: cp(3), sp(3)

    check_wall_collision = .false.

    !each point in cell1
    do i = 1, cell%nlat
    do j = 1, cell%nlon

      !each point in cell
      cp = cell%x(i, j, :)

      !each wall
      do i2 = 1, nwall
        !each point in the wall
        do j2 = nodeNum, walls(i2)%nvert, numNodes

          sp = walls(i2)%x(j2, :)

          !just do a distance check
          if (VecNorm(sp - cp) .le. threshold) then
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
    if (VecNorm(cell1%xc - cell2%xc) .ge. 4) return

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

        !AABB collision check, we find AABB intersection iff X, Y, and Z intervals all overlap
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
