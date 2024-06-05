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
  use MPI

  implicit none

  integer, parameter :: nrbcMax = 128
  integer, parameter :: nwallMax = 4
  type(t_rbc), pointer :: rbc
  type(t_rbc)         :: rbcRef
  type(t_wall), pointer :: wall
  real(WP) :: radEqv, szCell(3)
  integer :: nlat0, ii
  real(WP) :: th, xc(3)! , rand(16, 3), thetas(16)
  real(WP) :: xmin, xmax, ymin, ymax, zmin, zmax
  integer :: iz, i, dealias
  integer, parameter :: ranseed = 49765
  character(CHRLEN) :: fn
  real :: lengtube, lengspacing, phi, actlen
  integer :: ierr
  real(WP), allocatable :: rand(:,:), thetas(:)
  real(WP) :: xs(12), zs(12)

  ! Initialize
  call InitMPI
  ! Allocate working arrays
  allocate (rbcs(nrbcMax))
  allocate (walls(nwallMax))

  ! Wall
  nwall = 1
  wall => walls(1)

  call ReadWallMesh('Input/fullbifurc.e', wall)
  actlen = 13.33

  nrbc = 26
  nlat0 = 12
  dealias = 3
  phi = 70/real(100)
  lengtube = nrbc/real(phi)

  allocate (rand(nrbc, 3))
  allocate (thetas(nrbc))

  ! lengspacing = lengtube/Real(nrbc)
  lengspacing = 2

  wall%f = 0.

  ! do i = 1, wall%nvert
  !   th = ATAN2(wall%x(i, 1), wall%x(i, 2))
  !   ! 10 is the new tube diameter
  !   wall%x(i, 1) = 10/2.0*COS(th)    !!!!!!!!!!!!!!!!!!!!!!
  !   wall%x(i, 2) = 10/2.0*SIN(th)    !!!!!!!!!!!!!!!!!!!!!!
  !   wall%x(i, 3) = lengtube/actlen*wall%x(i, 3)   !!!!!!!!!!!!!!!!!!!
  ! end do

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
  ! lengspacing = lengtube/real(nrbc)

  print *, "lengtube: ", lengtube

  ! reference cell
  xc = 0.
  radEqv = 1.0

  call Rbc_Create(rbcRef, nlat0, dealias)
  call Rbc_MakeBiconcave(rbcRef, radEqv, xc)

  ! dimensions of starting rbc
  do ii = 1, 3
    szCell(ii) = maxval(rbcRef%x(:, :, ii)) - minval(rbcRef%x(:, :, ii))
  end do

  if (rootWorld) then
    print *, "rootWorld populating random # array"
    call random_number(rand)
    do i = 1, size(rand, 1)
      rand(i, 1) = 2.5*(rand(i, 1) - 0.5)
      rand(i, 2) = 2.5*(rand(i, 2) - 0.5)
      rand(i, 3) = rand(i, 3)
    end do
    ! rand = 2*(rand - 0.5)

    call random_number(thetas)
    thetas = (PI*thetas) - PI/2
  end if

  call MPI_Bcast(rand, size(rand), MPI_WP, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(thetas, size(thetas), MPI_WP, 0, MPI_COMM_WORLD, ierr)

  print *, "16x3 rand num array"
  do i = 1, size(rand, 1)
    print *, rand(i, :)
  end do

  print *, "16x1 thetas array", thetas(:)

  xs = (/6., 5.196, 3.527, -3.527, -5.196, -6., 6., 5.196, 3.527, -3.527, -5.196, -6./)
  zs = (/22.625, 25.625, 27.821, 27.821, 25.625, 22.625, 20.375, 17.374, 15.179, 15.179, 17.375, 20.375/)

  ! place rbcs in a line along z-axis
  do iz = 1, nrbc
    if (iz .le. 14) then
      if ((iz .eq. 6) .or. (iz .eq. 9)) then
        print *, "iz", iz, "moved 1.5x dist in x direction"
        xc(1) = 1.4*rand(iz, 1)
      else
        xc(1) = rand(iz, 1)
      end if
      ! xc(1) = rand(iz, 1)
      xc(2) = rand(iz, 2)
      if (iz .le. 7) then
        ! xc(3) = lengspacing*1.2*(iz - 0.5) + (rand(iz, 3) / 3)
        if (iz .eq. 7) then
          xc(3) = lengspacing*1.2*(iz - 0.5) - (rand(iz, 3) / 3)
        else
          xc(3) = lengspacing*1.2*(iz - 0.5) + (rand(iz, 3) / 3)
        end if
      else
        xc(3) = 9 + lengspacing*1.2*(iz - 0.5) + (rand(iz, 3) / 3)
      end if
    else
      xc(1) = xs(iz - 14)
      xc(2) = 0
      xc(3) = zs(iz - 14)
    end if
    
    print *, 'Xc', iz, xc

    rbc => rbcs(iz)
    rbc%celltype = 1
    call Rbc_Create(rbc, nlat0, dealias)
    call Rbc_MakeBiConcave(rbc, radEqv, xc)
    rbc%xc = xc

    ! call Rotate2(rbc, thetas(iz))
    call Rotate(rbc)
  end do

  ! Put things in the middle of the periodic box
  call Recenter_Cells_and_Walls

  ! Output
  write (*, '(A,3F10.3)') 'Periodic domain size = ', Lb

  Nt0 = 0; time = 0.
  vBkg(1:2) = 0.; vBkg(3) = 10.

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

  subroutine Rotate(cell)
    type(t_Rbc) :: cell
    integer :: ilat, ilon, i, j, ii
    real(WP) :: vec(3), rotmat(3, 3), vecog(3)
    real(WP) :: vsq, v1, v2, zv(3)

    vsq = 10
    do while (vsq .ge. 1)
      v1 = RandomNumber(ranseed)*2 - 1
      v2 = RandomNumber(ranseed)*2 - 1
      vsq = v1**2 + v2**2
    end do
    zv(1) = v1*2*sqrt(1 - vsq)
    zv(2) = v2*2*sqrt(1 - vsq)
    zv(3) = 1 - (2*vsq)

    ! print *, "cell%xc(3) before: ", cell%xc(:)

    ! vec = (/6, 9, 3/)
    ! print *, "VecNorm(vec): ", VecNorm(vec)
    ! vec = vec / VecNorm(vec)
    ! print *, "vec: ", vec
    rotmat = RotateMatrix(zv)

    do ii = 1, 3
      cell%x(:, :, ii) = cell%x(:, :, ii) - cell%xc(ii)
    end do ! ii

    forall (i=1:cell%nlat, j=1:cell%nlon)
      cell%x(i, j, :) = matmul(rotmat, cell%x(i, j, :))
    end forall

    do ii = 1, 3
      cell%x(:, :, ii) = cell%x(:, :, ii) + cell%xc(ii)
    end do ! ii

    ! do ilat = 1, cell%nlat
    !   do ilon = 1, cell%nlon
    !     cell%x(ilat, ilon, :) = reshape(matmul(cell%x(ilat, ilon, :), mat), (/3/))
    !   end do
    ! end do

    ! print *, "cell%xc(3) after: ", cell%xc(:)

  end subroutine Rotate

  subroutine Rotate2(cell, theta)
    type(t_Rbc) :: cell
    real(WP) :: theta, rotmat(3, 3), first(3), second(3), third(3), v(3), k(3)
    real(WP) :: khat, temp(3)
    integer :: ilat, ilon

    k = (/1, 1, 0/)
    k = k / VecNorm(k)

    print *, "theta", theta

    do ii = 1, 3
      cell%x(:, :, ii) = cell%x(:, :, ii) - cell%xc(ii)
    end do ! ii

    do ilat = 1, cell%nlat
      do ilon = 1, cell%nlon
        v = cell%x(ilat, ilon, :)
        first = v * COS(theta) 
        second = CrossProd(k, v) * SIN(theta)
        third = k*dot_product(k, v) * (1 - COS(theta))
        cell%x(ilat, ilon, :) = first + second + third
      end do         
    end do

    do ii = 1, 3
      cell%x(:, :, ii) = cell%x(:, :, ii) + cell%xc(ii)
    end do ! ii
    
  end subroutine Rotate2

end program InitCond
