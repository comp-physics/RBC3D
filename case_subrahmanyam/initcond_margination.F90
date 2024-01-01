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
  integer,parameter :: rbcl = 10
  integer,parameter :: rbcpl = 12
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
  real(WP) :: layerx(rbcpl), layery(rbcpl), randoms(3, rbcpl, rbcl)
  integer :: ierr

    ! Initialize
    call InitMPI
    ! Allocate working arrays
    allocate(rbcs(nrbcMax))
    allocate(walls(nwallMax))

    !note that perturbation is at most ~ 0.1*radius, so scale the perturbation to [-0.14, 0.14]
    !experimenting with larger perturbation of .2*radius
    if (rootWorld) then
        call random_number(randoms);
        randoms(:,:,:) = (randoms(:,:,:) - 0.5) * 2 * 0.14 * 2
    end if
    call MPI_Bcast(randoms, 3 * rbcpl * rbcl, MPI_WP, 0, MPI_COMM_WORLD, ierr);

    !Define layer points (hardcoded for now)
    layerx = (/0.0, -0.4685099164220228, 0.4685099164220228, &
        -0.24816347057168683, 0.24816347057168683, -0.7434871929505062, &
        0.7434871929505062, 0.0, -0.6511095339780459, 0.6511095339780459, &
        -0.27497727652848347, 0.27497727652848347/)

    layery = (/-0.7518365294283131, -0.5880107356137642, -0.5880107356137642, &
        -0.14327724653759516, -0.14327724653759516, -0.11173612173951383, &
        -0.11173612173951383, 0.2865544930751903, 0.37591826471415657, &
        0.37591826471415657, 0.6997468573532779, 0.6997468573532779/)

    ! Wall
    nwall = 1
    wall=>walls(1)

    ! if (10.ge.12) then
    !     call ReadMyWallMesh('Input/cylinder_D_8.01_L_13.39.dat', wall)
    !     actlen = 13.3921
    ! else
        call ReadWallMesh('Input/new_cyl_D10_L13_33.e', wall) !('Input/cylinder_D_6.0_L_13.33.g',wall)
        actlen = 13.33 !13.33 !/ sqrt(2.) ! minval(wall%x(:,3)) - maxval(wall%x(:,3))
    ! ! end if
    ! wall=>walls(2)
    ! call ReadWallMesh('Input/new_cyl_D2_L13_33.e', wall)
    ! rbcl = 10

    nrbc = rbcl * rbcpl
    nlat0 = 12
    dealias = 3
    phi = 70/real(100)

    lengspacing = (1.4 / 4) * (90 * rbcpl) / (400 * PI * 0.2) !lengtube/real(nrbc)

    lengtube = rbcl * lengspacing !nrbc/real(phi) !XXLL

    wall%f = 0.

    do i = 1,wall%nvert
        th = ATAN2(wall%x(i,1),wall%x(i,2))
        wall%x(i,1) = 7*COS(th)    !!!!!!!!!!!!!!!!!!!!!!
        wall%x(i,2) = 7*SIN(th)    !!!!!!!!!!!!!!!!!!!!!!
        wall%x(i,3) = lengtube/actlen*wall%x(i,3)   !!!!!!!!!!!!!!!!!!!
    end do

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

    ! 1.4 / 4 is the micron -> simulation units conversion (4 microns = 1.4 sim units)
    ! 90 is ~ um^3 volume of a RBC
    ! 400 * PI is the surface area of the tube's cross-section, times hematocrit of 0.2

    ! Reference cell 
    xc = 0.
    radEqv = 1.0

    call Rbc_Create(rbcRef, nlat0, dealias)
    call Rbc_MakeBiconcave(rbcRef, radEqv, xc)

    do ii = 1,3
        szCell(ii) = maxval(rbcRef%x(:,:,ii)) - minval(rbcRef%x(:,:,ii))
    end do

    ! place cells
    do iz = 1,rbcl
        do i = 1,rbcpl
            xc(1) = layerx(i) * 6 + randoms(3, i, iz)
            xc(2) = layery(i) * 6 + randoms(3, i, iz) ! multiply by tube radius
            xc(3) = iz * lengspacing + randoms(3, i, iz)  ! the Z is too close, will have to re-arrange + randoms(3, i, iz) - 0.5

            rbc => rbcs((iz - 1) * rbcpl + i)

            ! set 1/10 cells to be sickles
            if (modulo(((iz - 1) * rbcpl + i), 10) .ne. 0) then
                rbc%celltype = 1
                call Rbc_Create(rbc, nlat0, dealias)
                call Rbc_MakeBiConcave(rbc, radEqv, xc)
            else
                rbc%celltype = 3
                call ImportReadRBC('Input/SickleCell_12nlat.dat', rbc, xc)
            end if 

        end do
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

        !Write out only type-1 cells (healthy RBC cells)
        write(fn, FMT=fn_FMT) 'D/', '1x', 0, '.dat'
        call WriteManyRBCsByType(fn, nrbc, rbcs, 1)
        ! Write out only type-3 cells (sickle cells)
        write(fn, FMT=fn_FMT) 'D/', '3x', 0, '.dat'
        call WriteManyRBCsByType(fn, nrbc, rbcs, 3)

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
