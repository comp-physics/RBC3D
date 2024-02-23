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

  integer,parameter :: nrbcMax = 128
  integer,parameter :: nwallMax = 4
  type(t_rbc),pointer :: rbc
  type(t_rbc)         :: rbcRef
  type(t_wall),pointer :: wall
  real(WP) :: radEqv, szCell(3)
  integer :: nlat0, ii
  real(WP) :: theta, th, rc, xc(3), ztemp
  real(WP) :: xmin, xmax, ymin, ymax, zmin, zmax
  integer :: iz, i, p, l, dealias
  integer,parameter :: ranseed = 161269
  character(CHRLEN) :: fn
  real :: lengtube,lengspacing, phi, actlen
  real(WP) :: rand(40, 3)
  integer :: j, ierr, half2, index, layer, newiz
  real :: wbcs(2), tubeDiam, layerx(4), layery(4)

    ! Initialize
    call InitMPI
    ! Allocate working arrays
    allocate(rbcs(nrbcMax))
    allocate(walls(nwallMax))

    tubeDiam = 8

    ! Wall
    nwall = 1
    wall=>walls(1)

    call ReadWallMesh('Input/cylwithhole.e',wall)
    actlen = 13.33

    nrbc = 40
    nlat0 = 12
    dealias = 3
    phi = 70/real(100)
    lengtube = 17.73 ! nrbc/real(phi) !XXLL

    lengspacing = .51 ! lengtube/Real(nrbc)

    print*, 'lengtube 1', lengtube
    print*, 'lengspacing 1', lengspacing

    wall%f = 0.

    do i = 1,wall%nvert
        th = ATAN2(wall%x(i,1),wall%x(i,2))
        ! 10 is the diameter
        wall%x(i,1) = tubeDiam/2.0*COS(th)    !!!!!!!!!!!!!!!!!!!!!!
        wall%x(i,2) = tubeDiam/2.0*SIN(th)    !!!!!!!!!!!!!!!!!!!!!!
        wall%x(i,3) = lengtube/actlen*wall%x(i,3)   !!!!!!!!!!!!!!!!!!!
    end do
    xmin = minval(wall%x(:,1))
    xmax = maxval(wall%x(:,1))

    ymin = minval(wall%x(:,2))
    ymax = maxval(wall%x(:,2))

    zmin = minval(wall%x(:,3))
    zmax = maxval(wall%x(:,3))

    ! size of the periodic box
    Lb(1) = xmax - xmin + 0.5
    Lb(2) = Lb(1)
    Lb(3) = zmax - zmin
    lengtube = Lb(3)
    ! lengspacing = lengtube/real(nrbc)
    print*, 'lengtube 2', lengtube
    print*, 'lengspacing 2', lengspacing

    ! Reference cell i
    xc = 0.
    ! radius of cell gets trans to ~1.4
    ! if you want rad of 8, pass in 5.7 for 8/1.4 (5.7*1.4=8)
    radEqv = 1.0

    call Rbc_Create(rbcRef, nlat0, dealias)
    call Rbc_MakeBiconcave(rbcRef, radEqv, xc)

    do ii = 1,3
        szCell(ii) = maxval(rbcRef%x(:,:,ii)) - minval(rbcRef%x(:,:,ii))
    end do

    if (rootWorld) then
        call random_number(rand)
    end if
    call MPI_Bcast(rand, 40*3, MPI_WP, 0, MPI_COMM_WORLD, ierr)
    ! call MPI_Bcast(wbcs, 2, MPI_WP, 0, MPI_COMM_WORLD, ierr)

    wbcs = 1 + FLOOR(40*wbcs)

    print*, '40 x 3 rand num array'
    do j = 1, nrbc
        print *, rand(j, :)
    end do

    layerx = (/-1, -1, 1, 1/)
    layery = (/-1, 1, -1, 1/)

    do iz = 1, nrbc
        index = (modulo(iz, 4))
        layer = (iz - 1) / 4
        xc(1) = layerx(index)
        xc(2) = layery(index)
        xc(3) = (layer + .5) + (lengspacing * layer)
        print*, 'Xc', iz, xc
        rbc => rbcs(iz)
        rbc%celltype = 1
        call Rbc_Create(rbc, nlat0, dealias)
        call Rbc_MakeBiConcave(rbc, radEqv, xc)
    end do


    ! Put things in the middle of the periodic box
    call Recenter_Cells_and_Walls

    ! Output
    write(*, '(A,3F10.3)') 'Periodic domain size = ', Lb

    Nt0 = 0; time = 0.
    vBkg(1:2) = 0.; vBkg(3) = 8.

    ! Write intial conditions
    if (nrbc > 0) then
        write(fn, FMT=fn_FMT) 'D/', 'x', 0, '.dat'
        call WriteManyRBCs(fn, nrbc, rbcs )
        write(*, '(A,A)') 'Cell file: ', trim(fn)

        fn = 'D/restart.LATEST.dat'
        call WriteRestart(fn, Nt0, time)
        write(*, '(A,A)') 'Binary restart file: ', trim(fn)
    end if

    if (nwall > 0) then
        write(fn, FMT=fn_FMT) 'D/', 'wall', 0, '.dat'
        call WriteManyWalls(fn, nwall, walls )
        write(*, '(A,A)') 'Wall file: ', trim(fn)
    end if

    write(fn, FMT=fn_FMT) 'D/', 'restart', Nt0, '.dat'
    call WriteRestart(fn, Nt0, time)
    write(*, '(A,A)') 'Binary restart file: ', trim(fn)

    ! Deallocate working arrays
    deallocate(rbcs)

    ! Finalize
    call FinalizeMPI

    stop

contains

subroutine Recenter_Cells_and_Walls
    integer :: irbc, iwall, ii, gen, repeat, j
    real :: x, y, a
    real, parameter :: PI = 3.14159265359
    type(t_rbc),pointer :: rbc
    type(t_wall),pointer :: wall
    real(WP) :: xmin(3), xmax(3), xc(3)

    ! cells
    ! translate cells
    do irbc = 1, nrbc
        rbc => rbcs(irbc)
        rbc%x(:,:,1) = rbc%x(:,:,1) + 0.5*Lb(1)
        rbc%x(:,:,2) = rbc%x(:,:,2) + 0.5*Lb(2)
        print*,'Zc', sum(rbc%x(:,:,3))/real(rbc%nlat*rbc%nlon)
    end do

    ! walls
    ! find the bounding box
    xmin = 1.D10
    xmax = -1.D10
    do iwall = 1,nwall
        wall => walls(iwall)
        do ii = 1,3
            xmin(ii) = min(xmin(ii), minval(wall%x(:,ii)) )
            xmax(ii) = max(xmax(ii), maxval(wall%x(:,ii)) )
        end do
    end do
    xc = 0.5*(xmin + xmax)

    ! translate walls
    do iwall = 1,nwall
        wall => walls(iwall)
        do ii = 1,3
            wall%x(:,ii) = wall%x(:,ii) + 0.5*Lb(ii) - xc(ii)
        end do
    end do

end subroutine Recenter_Cells_and_Walls

end program InitCond