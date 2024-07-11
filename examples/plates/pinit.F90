! Create restart file with initial condition
program InitCond

  use ModDataTypes
  use ModDataStruct
  use ModRbc
  use ModWall
  use ModConf
  use ModData
  use ModIO
  use ModBasicMath

  implicit none

  integer, parameter :: nrbcMax = 128
  integer, parameter :: nwallMax = 4
  type(t_rbc), pointer :: rbc
  type(t_rbc)         :: rbcRef
  type(t_wall), pointer :: wall
  real(WP) :: radEqv, szCell(3), radPlat
  integer :: nlat0, ii, nlatp
  real(WP) :: th, xc(3)
  real(WP) :: xmin, xmax, ymin, ymax, zmin, zmax
  integer :: iz, i, dealias
  integer, parameter :: ranseed = 161269
  character(CHRLEN) :: fn
  real :: lengtube, lengspacing, phi, actlen

  ! Initialize
  call InitMPI
  ! Allocate working arrays
  allocate (rbcs(nrbcMax))
  allocate (walls(nwallMax))

  ! Wall
  nwall = 1
  wall => walls(1)

  call ReadWallMesh('Input/new_cyl_D6_L13_33.e', wall)
  actlen = 13.33

  nrbc = 1
  nlat0 = 12
  nlatp = 4
  dealias = 3
  phi = 70/real(100)
  ! lengtube = nrbc/real(phi) !XXLL
  lengtube = 11

  ! lengspacing = lengtube/Real(nrbc)

  wall%f = 0.

  do i = 1, wall%nvert
    th = ATAN2(wall%x(i, 1), wall%x(i, 2))
    ! 10 is the new tube diameter
    wall%x(i, 1) = 10/2.0*COS(th)    !!!!!!!!!!!!!!!!!!!!!!
    wall%x(i, 2) = 10/2.0*SIN(th)    !!!!!!!!!!!!!!!!!!!!!!
    wall%x(i, 3) = lengtube/actlen*wall%x(i, 3)   !!!!!!!!!!!!!!!!!!!
  end do
  xmin = minval(wall%x(:, 1))
  xmax = maxval(wall%x(:, 1))

  ymin = minval(wall%x(:, 2))
  ymax = maxval(wall%x(:, 2))

  zmin = minval(wall%x(:, 3))
  zmax = maxval(wall%x(:, 3))

  ! size of the periodic box
  Lb(1) = xmax - xmin + 0.5
  Lb(2) = Lb(1)
  Lb(3) = zmax - zmin
  lengtube = Lb(3)
  lengspacing = lengtube/real(nrbc)
  print *, "lengtube: ", lengtube

  ! reference cell (unnecessary)
  xc = 0.
  radEqv = 1.0
  radPlat = .4

  call Rbc_Create(rbcRef, nlatp, dealias)
  call Rbc_MakePlatelet(rbcRef, radPlat, xc)

  ! dimensions of starting rbc
  do ii = 1, 3
    szCell(ii) = maxval(rbcRef%x(:, :, ii)) - minval(rbcRef%x(:, :, ii))
  end do

  print *, "szCell: ", szCell

  ! place 1 platelet at beginning of tube
  iz = 1
  xc(1:2) = 0.
  xc(3) = 1.
  print *, 'Xc', iz, xc
  rbc => rbcs(iz)
  rbc%celltype = 3
  call Rbc_Create(rbc, nlatp, dealias)
  call Rbc_MakePlatelet(rbc, radPlat, xc)

  ! do iz = 1, nrbc
  !   xc(1:2) = 0.
  !   xc(3) = lengspacing*(iz - 0.5)
  !   print *, 'Xc', iz, xc
  !   rbc => rbcs(iz)
  !   rbc%celltype = 1
  !   call Rbc_Create(rbc, nlat0, dealias)
  !   call Rbc_MakeBiConcave(rbc, radEqv, xc)
  ! end do

  ! Put things in the middle of the periodic box
  call Recenter_Cells_and_Walls

  ! Output
  write (*, '(A,3F10.3)') 'Periodic domain size = ', Lb

  Nt0 = 0; time = 0.
  vBkg(1:2) = 0.; vBkg(3) = 4.

  ! Write initial conditions
  if (nrbc > 0) then
    write (fn, FMT=fn_FMT) 'D/', 'x', 0, '.dat'
    call WriteManyRBCs(fn, nrbc, rbcs)
    write (*, '(A,A)') 'Cell file: ', trim(fn)

    fn = 'D/restart.LATEST.dat'
    call WriteRestart(fn, Nt0, time)
    write (*, '(A,A)') 'Binary restart file: ', trim(fn)
  end if

  if (nwall > 0) then
    write (fn, FMT=fn_FMT) 'D/', 'wall', 0, '.dat'
    call WriteManyWalls(fn, nwall, walls)
    write (*, '(A,A)') 'Wall file: ', trim(fn)
  end if

  write (fn, FMT=fn_FMT) 'D/', 'restart', Nt0, '.dat'
  call WriteRestart(fn, Nt0, time)
  write (*, '(A,A)') 'Binary restart file: ', trim(fn)

  ! Deallocate working arrays
  deallocate (rbcs)

  print *, "after deallocate"

  ! Finalize
  call FinalizeMPI

  print *, "FinalizeMPI"

  stop

contains

  subroutine Recenter_Cells_and_Walls
    integer :: irbc, iwall, ii, gen, repeat, j
    real :: x, y, a
    real, parameter :: PI = 3.14159265359
    type(t_rbc), pointer :: rbc
    type(t_wall), pointer :: wall
    real(WP) :: xmin(3), xmax(3), xc(3)

    ! cells
    ! translate cells
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      rbc%x(:, :, 1) = rbc%x(:, :, 1) + 0.5*Lb(1)
      rbc%x(:, :, 2) = rbc%x(:, :, 2) + 0.5*Lb(2)
      print *, 'Zc', sum(rbc%x(:, :, 3))/real(rbc%nlat*rbc%nlon)
    end do

    ! walls
    ! find the bounding box
    xmin = 1.D10
    xmax = -1.D10
    do iwall = 1, nwall
      wall => walls(iwall)
      do ii = 1, 3
        xmin(ii) = min(xmin(ii), minval(wall%x(:, ii)))
        xmax(ii) = max(xmax(ii), maxval(wall%x(:, ii)))
      end do
    end do
    xc = 0.5*(xmin + xmax)

    ! translate walls
    do iwall = 1, nwall
      wall => walls(iwall)
      do ii = 1, 3
        wall%x(:, ii) = wall%x(:, ii) + 0.5*Lb(ii) - xc(ii)
      end do
    end do

  end subroutine Recenter_Cells_and_Walls

end program InitCond
