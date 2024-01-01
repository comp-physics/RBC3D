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

  integer,parameter :: nrbcMax = 128
  integer,parameter :: nwallMax = 4
  integer :: rbcl
  type(t_rbc),pointer :: rbc
  type(t_rbc)         :: rbcRef
  type(t_wall),pointer :: wall
  real(WP) :: radEqv, szCell(3)
  integer :: nlat0, ii
  real(WP) :: theta, th, rc, xc(3)
  real(WP) :: xmin, xmax, ymin, ymax, zmin, zmax
  integer :: iz, i, p, l, dealias
  integer,parameter :: ranseed = 161269
  character(CHRLEN) :: fn
  real :: lengtube,lengspacing, phi, actlen
  

    ! Initialize
    call InitMPI
    ! Allocate working arrays
    allocate(rbcs(nrbcMax))
    allocate(walls(nwallMax))

    ! Wall
    nwall = 1
    wall=>walls(1)

    ! if (10.ge.12) then
    !     call ReadMyWallMesh('Input/cylinder_D_8.01_L_13.39.dat', wall)
    !     actlen = 13.3921
    ! else
        call ReadWallMesh('Input/valve.e', wall) !('Input/cylinder_D_6.0_L_13.33.g',wall)
        actlen = 10 !13.33 !/ sqrt(2.) ! minval(wall%x(:,3)) - maxval(wall%x(:,3))
    ! ! end if
    ! wall=>walls(2)
    ! call ReadWallMesh('Input/new_cyl_D2_L13_33.e', wall)

    nrbc = 1
    nlat0 = 12
    dealias = 3
    phi = 70/real(100)
    ! lengtube = nrbc/real(phi) !XXLL

    lengspacing = lengtube/Real(nrbc)

    wall%f = 0.

    !Use the larger wall for periodic boundary box set up
    wall=>walls(1)

    xmin = minval(wall%x(:,1))
    xmax = maxval(wall%x(:,1))

    ymin = minval(wall%x(:,2))
    ymax = maxval(wall%x(:,2))

    zmin = minval(wall%x(:,3))
    zmax = maxval(wall%x(:,3))

    ! size of the periodic box
    Lb(1) = (xmax - xmin + 0.5) !* 10
    Lb(2) = (ymax - ymin + 0.5) !* 10 !make really far away periodic box, for debugging ! Lb(1)
    Lb(3) = zmax - zmin
    lengtube = Lb(3)
    lengspacing = lengtube/real(nrbc)

    ! Reference cell 
    xc = 0.
    radEqv = 1.0

    call Rbc_Create(rbcRef, nlat0, dealias)
    call Rbc_MakeBiconcave(rbcRef, radEqv, xc)

    do ii = 1,3
        szCell(ii) = maxval(rbcRef%x(:,:,ii)) - minval(rbcRef%x(:,:,ii))
    end do

    ! place cells
    do iz = 1,nrbc
        ! do i = 1,9
        ! xc(1:2) = 0.
        xc(1) = 0 !3
        xc(2) = -5 !3
        ! xc(3) = 0 !10 
        xc(3) = lengspacing*(iz-0.5) ! - lengtube / 2
        print*, 'Xc', iz, xc
            ! xc(1) = layerx(i) * 20 ! + 2*(randoms(1, i, iz) - 0.5)
            ! xc(2) = layery(i) * 20 ! + 2*(randoms(2, i, iz) - 0.5) ! multiply by tube radius
            ! xc(3) = iz * 6 - 3   ! the Z is too close, will have to re-arrange   + randoms(3, i, iz) - 0.5

        rbc => rbcs(iz)
        rbc%celltype = 1
        call Rbc_Create(rbc, nlat0, dealias)
        call Rbc_MakeBiConcave(rbc, radEqv, xc)
    ! end do
    end do

    ! Put things in the middle of the periodic box
    call Recenter_Cells_and_Walls

    ! Output
    write(*, '(A,3F10.3)') 'Periodic domain size = ', Lb

    Nt0 = 0; time = 0.
    vBkg(1:2) = 0.; vBkg(3) = 8.
    ! vBkg(1) = 0.; vBkg(2:3) = 5.

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
