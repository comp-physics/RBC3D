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
  use ModPostProcess
  ! use MPI

  implicit none

  integer, parameter :: nrbcMax = 128
  integer, parameter :: nwallMax = 4
  type(t_rbc), pointer :: rbc
  type(t_rbc)         :: rbcRef
  type(t_wall), pointer :: wall
  real(WP) :: radEqv, szCell(3)
  integer :: nlat0, nlat1, ii
  real(WP) :: theta, th, rc, xc(3), ztemp
  real(WP) :: xmin, xmax, ymin, ymax, zmin, zmax
  integer :: iz, i, p, l, dealias
  integer, parameter :: ranseed = 161269
  character(CHRLEN) :: fn
  ! real = double, fp
  real :: lengtube, lengspacing, phi, actlen
  real(WP) :: rand(27, 3), minDist, platRad, platDiam
  real :: tubeRad

  ! Initialize
  call InitMPI
  ! Allocate working arrays
  allocate (rbcs(nrbcMax))
  ! allocate(walls(nwallMax))

  tubeRad = 4.0

  nrbc = 2
  nlat0 = 12
  nlat1 = 4
  dealias = 3
  phi = 70/real(100)
  ! 11.34
  lengtube = 5 ! nrbc/real(phi) !XXLL

  lengspacing = (lengtube - ((2.62/2.82)*9))/9 ! lengtube/Real(nrbc)
  ! print *, "1"
  nwall = 1
  allocate (walls(nwall))
  wall => walls(1)

  call ReadWallMesh('Input/new_cyl_D6_L13_33.e', wall)
  actlen = 13.33
  wall%f = 0.
  do i = 1, wall%nvert
    th = ATAN2(wall%x(i, 1), wall%x(i, 2))
    wall%x(i, 1) = (tubeRad)*COS(th)    !!!!!!!!!
    wall%x(i, 2) = (tubeRad)*SIN(th)    !!!!!!!!!
    wall%x(i, 3) = lengtube/actlen*wall%x(i, 3)
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

  print *, "Lb", Lb(:)
  !  Lb   8.3013429849047267        8.3013429849047267        11.347517730496456
  ! lengspacing = lengtube/real(nrbc)
  print *, 'lengtube 2', lengtube
  print *, 'lengspacing 2', lengspacing

  ! Reference cell i
  xc = 0.

  radEqv = 1.0
  platRad = .3

  print *, "platRad", platRad

  call Rbc_Create(rbcRef, nlat0, dealias)
  call Rbc_MakePlatelet(rbcRef, platRad, xc)

  do ii = 1, 3
    ! dimensions of cell !
    szCell(ii) = maxval(rbcRef%x(:, :, ii)) - minval(rbcRef%x(:, :, ii))
  end do
  if (rootWorld) then
    print *, "szCell", szCell
  end if

  ! place cells
  xc(1:2) = 0.
  xc(3) = 1.
  ! print *, 'rbc iz:', iz, 'xc:', xc
  rbc => rbcs(1)
  rbc%celltype = 1
  call Rbc_Create(rbc, nlat0, dealias)
  call RBC_MakeBiConcave(rbc, radEqv, xc)

  xc(1:2) = 0.
  xc(3) = 4.
  ! print *, 'rbc iz:', iz, 'xc:', xc
  rbc => rbcs(2)
  rbc%celltype = 3
  call Rbc_Create(rbc, nlat1, dealias)
  call Rbc_MakePlatelet(rbc, platRad, xc)
  

  ! Put things in the middle of the periodic box
  call Recenter_Cells_and_Walls

  ! platDiam = GetDiameter(3)

  ! Output
  write (*, '(A,3F10.3)') 'Periodic domain size = ', Lb

  Nt0 = 0; time = 0.
  vBkg(1:2) = 0.; vBkg(3) = 8

  ! Write initial conditions
  if (nrbc > 0) then
    write (fn, FMT=fn_FMT) 'D/', 'x', 0, '.dat'
    call WriteManyRBCs(fn, nrbc, rbcs)

    write(fn, FMT=fn_FMT) 'D/', '3x', 0, '.dat'
    call WriteManyRBCsByType(fn, nrbc, rbcs, 3)

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

  ! Finalize
  call FinalizeMPI

  stop

contains
  subroutine Recenter_Cells_and_Walls
    integer :: irbc, iwall, ii, gen, repeat, j
    real :: x, y, a
    real, parameter :: PI = 3.14159265359
    type(t_rbc), pointer :: rbc
    type(t_wall), pointer :: wall
    real(WP) :: xmin(3), xmax(3), xc(3), xmax2, ymax2

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
      print *, "iwall", iwall
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

    xmax2 = maxval(wall%x(:, 1))
    print *, "xmax2: ", xmax2

    ymax2 = maxval(wall%x(:, 2))
    print *, "ymax2: ", ymax2

  end subroutine Recenter_Cells_and_Walls

end program InitCond