! takes a RESTART file
! Reads RESTART file
! Adds spaces between the RBCs
! and adds specified # of Sickle Cells
! Then writes the new cells-setup out to separate restart-file

program randomized_cell_gen

  use ModDataTypes
  use ModDataStruct
  use ModRbc
  use ModWall
  use ModConf
  use ModData
  use ModIO
  use ModBasicMath

  integer, parameter :: ranseed = 484

  type(t_Rbc), pointer :: rbc
  type(t_Wall), pointer :: wall
  character(CHRLEN) :: fn

  integer :: zmax
  integer :: i
  real(WP) :: th, actlen

  real(WP), parameter :: hematocrit = 0.25
  integer :: nrbcMax ! how many cells
  real(WP) :: tuber, tubelen

  real(WP) :: clockBgn, clockEnd

  call InitMPI

  !initialize parameter variables, ensuring for the defined hematocrit
  ! hematocrit = 4 * nrbc / (tube_radius^2 * tube_length)
  ! assuming all blood cells are healthy for the calculation
  tuber = 7
  nrbcMax = 50
  tubelen = 4*nrbcMax/(tuber**2*hematocrit)

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

  ! wall=>walls(2)

  ! call ReadWallMesh('Input/new_cyl_D6_L13_33_hires.e',wall)
  ! actlen = 13.33

  ! wall%f = 0.
  ! do i = 1,wall%nvert
  !     th = ATAN2(wall%x(i, 1), wall%x(i, 2))
  !     wall%x(i,1) = tuber*COS(th)    !!!!!!!!!
  !     wall%x(i,2) = tuber*SIN(th)    !!!!!!!!!
  !     wall%x(i,3) = tubelen/actlen*wall%x(i,3)
  ! end do

  !for each cell to add
  allocate (rbcs(nrbcMax))
  nrbc = 0
  do i = 1, nrbcMax
    clockBgn = MPI_Wtime()
    call place_cell
    clockEnd = MPI_Wtime()

    nrbc = nrbc + 1
    if (rootWorld) write (*, *) "Cell #", nrbc, " Placed; Time Cost = ", clockEnd - clockBgn
  end do

  !recenter everything to be all positive
  call recenter

  !Write Cells Out to Cell TecPlot File
  write (fn, FMT=fn_FMT) 'D/', 'x', 0, '.dat'
  call WriteManyRBCs(fn, nrbc, rbcs)
  write (fn, FMT=fn_FMT) 'D/', '1x', 0, '.dat'
  call WriteManyRBCsByType(fn, nrbc, rbcs, 1)
  write (fn, FMT=fn_FMT) 'D/', '2x', 0, '.dat'
  call WriteManyRBCsByType(fn, nrbc, rbcs, 2)
  write (fn, FMT=fn_FMT) 'D/', '3x', 0, '.dat'
  call WriteManyRBCsByType(fn, nrbc, rbcs, 3)
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

  subroutine place_cell

    type(t_Rbc) :: newcell
    type(t_Rbc), pointer :: cell
    real(WP) :: tmp_xc(3)

    integer :: celli, ii, i, j, i2, j2
    real(WP) :: c1p(3), c2p(3), sp(3)

    !create new template rbc
    call Rbc_Create(newcell, 12, 3)
    call Rbc_MakeBiconcave(newcell, 1.)
    newcell%celltype = 1

    do

      !randomly rotate cell
888   call rotate_cell(newcell)

      !randomly select a tmp_xc
      tmp_xc(2) = RandomNumber(ranseed)*2*PI
      tmp_xc(3) = sqrt(RandomNumber(ranseed))*(tuber)
      tmp_xc(1) = tmp_xc(3)*cos(tmp_xc(2))
      tmp_xc(2) = tmp_xc(3)*sin(tmp_xc(2))
      tmp_xc(3) = RandomNumber(ranseed)*tubelen - tubelen/2

      !shift sickle by the xc
      do ii = 1, 3
        newcell%x(:, :, ii) = newcell%x(:, :, ii) + tmp_xc(ii)
      end do

      !check if cell sticks outside the wall
      !assume only cylinder walls for now
      !for other wall-shapes, run collision detection algorithm between newcell & wall-mesh
      if (any(abs(newcell%x(:, :, 3)) .ge. tubelen/2 .or. &
              newcell%x(:, :, 1)**2 + newcell%x(:, :, 2)**2 .ge. tuber**2)) then
        do ii = 1, 3
          newcell%x(:, :, ii) = newcell%x(:, :, ii) - tmp_xc(ii)
        end do

        goto 888
      end if

      do celli = 1, nrbc
        cell => rbcs(celli)

        !all combinations of 2 vertices on the cell
        ! forall (i = 1:cell%nlat - 1, j=1:cell%nlon - 1)
        do i = 1, cell%nlat - 1
        do j = 1, cell%nlon - 1

          c1p = cell%x(i, j, :)
          c2p = cell%x(i + 1, j + 1, :)

          !each point on the new cell
          do i2 = 1, newcell%nlat
          do j2 = 1, newcell%nlon

            sp = newcell%x(i2, j2, :)

            if ( &
              ((c1p(1) .le. sp(1) .and. sp(1) .le. c2p(1)) .or. (c2p(1) .le. sp(1) .and. sp(1) .le. c1p(1))) .and. &  !x direction overlap
              ((c1p(2) .le. sp(2) .and. sp(2) .le. c2p(2)) .or. (c2p(2) .le. sp(2) .and. sp(2) .le. c1p(2))) .and. &  !y direction overlap
              ((c1p(3) .le. sp(3) .and. sp(3) .le. c2p(3)) .or. (c2p(3) .le. sp(3) .and. sp(3) .le. c1p(3))) &        !z direction overlap
              ) then
              !shift new cell back to center by the xc
              do ii = 1, 3
                newcell%x(:, :, ii) = newcell%x(:, :, ii) - tmp_xc(ii)
              end do !ii

              ! if (rootWorld) write(*,*) "Cell # ", nrbc + 1, " attempt failed"
              goto 888
            end if

          end do
          end do
        end do
        end do

      end do !celli

      !the cell is fine
      !replace an existing rbc in rbcs[] with the new sickle
      !exit out
      rbcs(nrbc + 1) = newcell
      exit

    end do

  end subroutine place_cell

  subroutine place_wbc

    type(t_Rbc) :: newcell
    type(t_Rbc), pointer :: cell
    real(WP) :: xc(3), radEqv

    integer :: nlat0, dealias

    nlat0 = 12
    dealias = 3
    radEqv = 1.0

    !create new template rbc
    call Rbc_Create(newcell, 12, 3)
    call Rbc_MakeLeukocyte(newcell, radEqv, xc)
    newcell%celltype = 2

    rbcs(nrbc + 1) = newcell

  end subroutine place_wbc

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

end program randomized_cell_gen
