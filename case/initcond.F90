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
  use ModTimeInt
  use ModQuadRule
  use ModSphpk
  use ModPolarPatch
  implicit none

  integer,parameter :: nrbcMax = 128
  integer,parameter :: nwallMax = 4
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
  real :: lengtube, lengspacing, phi

  integer :: nlat, nlon
  real, allocatable, dimension(:,:,:) :: base_cell

    call InitAll
    call GlobData_Finalize
    call TimeInt_Finalize
    call FinalizeSystem

    call InitSystem
    call SetEwaldPrms
    call DomainDecomp
    call GlobData_Init
    call TimeInt_Init

    !get base_cell shape here
    rbc => rbcs(1)
    nlat = rbc%nlat; nlon = rbc%nlon
    allocate( base_cell(nlat,nlon,3) )

    base_cell(:,:,:) = rbc%x(:,:,:)
    do i = 1,3
        base_cell(:,:,i) = base_cell(:,:,i) - ( sum( base_cell(:,:,i) )/real(nlat*nlon) )
    end do

    print*, 'system finalize'
    call GlobData_Finalize
    call TimeInt_Finalize
    call FinalizeSystem

    print*, 'finalize all'
    call FinalizeAll

    ! Initialize
    print*, 'initmpi'
    call InitMPI

    ! Allocate working arrays
    allocate(rbcs(nrbcMax))
    allocate(walls(nwallMax))


    ! Wall
    nwall = 1
    wall=>walls(1)

    call ReadMyWallMesh('Input/cylinder_D_8.01_L_13.39.dat', wall)

    nrbc = 1
    nlat0 = 12
    dealias = 3
    phi = 70/real(100)
    lengtube = nrbc/real(phi) !XXLL

    lengspacing = lengtube/Real(nrbc)

    wall%f = 0.

    do i = 1,wall%nvert
        th = ATAN2(wall%x(i,1),wall%x(i,2))
        wall%x(i,1) = 10./2.0*COS(th)    !!!!!!!!!!!!!!!!!!!!!!
        wall%x(i,2) = 10./2.0*SIN(th)    !!!!!!!!!!!!!!!!!!!!!!
        wall%x(i,3) = lengtube/13.3921*wall%x(i,3)   !!!!!!!!!!!!!!!!!!!
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

    ! Reference cell
    xc = 0.
    radEqv = 1.0

    print*, 'get rbcs'
    call Rbc_Create(rbcRef, nlat0, dealias)
    call Rbc_MakeBiconcave(rbcRef, radEqv, xc)

    do ii = 1,3
        szCell(ii) = maxval(rbcRef%x(:,:,ii)) - minval(rbcRef%x(:,:,ii))
    end do

    ! place cells
    do iz = 1,nrbc
        xc(1:2) = 0.
        xc(3) = lengspacing*iz
        print*, 'Xc', iz, xc

        rbc => rbcs(iz)
        rbc%celltype = 1
        call Rbc_Create(rbc, nlat0, dealias)
        call Rbc_MakeBiConcave(rbc, radEqv, xc)
        rbc%x(:,:,:) = base_cell(:,:,:)
        rbc%x(:,:,3) = rbc%x(:,:,3) + xc(3)
    end do

    ! Put things in the middle of the periodic box
    call Recenter_Cells_and_Walls
    call ReboxRbcs

    ! Output
    write(*, '(A,3F10.3)') 'Periodic domain size = ', Lb

    Nt0 = 0; time = 0.
    vBkg(1:2) = 0.; vBkg(3) = 8.

    print*, 'write rbcs'
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
!    call FinalizeMPI

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
        !rbc%x(:,:,3) = rbc%x(:,:,3) - 0.5*Lb(3)
        print*,'Zc', sum(rbc%x(:,:,3))/real(rbc%nlat*rbc%nlon)
   end do

    !print*, sum(rbc%x(:,:,1))/real(rbc%nlat*rbc%nlon)
    !print*, sum(rbc%x(:,:,2))/real(rbc%nlat*rbc%nlon)

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

    !print*, minval(wall%x(:,1)), minval(wall%x(:,3))
    !print*, maxval(wall%x(:,1)), maxval(wall%x(:,3))

end subroutine Recenter_Cells_and_Walls



 subroutine InitAll

    ! System intialization
    call InitMPI()
    call GaussQuad_Init
    call ReadConfig('Input/tube_base.in')
    call InitSystem

    ! Note: SetEwaldPrms and DomainDecomp must be called after InitSystem
    call SetEwaldPrms
    call DomainDecomp
    call GlobData_Init

    ! Prepare for time integration
    call IO_Init
    call TimeInt_Init

  end subroutine InitAll

 subroutine FinalizeAll

    call IO_Finalize
    !call TimeInt_Finalize

    !call GlobData_Finalize
    !call FinalizeSystem

    call GaussQuad_Finalize
!    call FinalizeMPI

 end subroutine FinalizeAll

subroutine InitSystem

    integer :: irbc, iwall
    type(t_Rbc),pointer :: rbc,rbcRef
    type(t_Wall),pointer :: wall
    integer :: nlat0
    real(WP) :: radEqv
    integer :: ierr, nwalls

    nwalls = 1

    call ReadRestart('restart_basecell')


    ! Reference cells
    allocate(rbcRefs(2))

    if (nrbc > 0) then
      radEqv = 1.
      nlat0 = rbcs(1)%nlat0

      rbcRef => rbcRefs(1)
      call RBC_Create(rbcRef, nlat0)
      call RBC_MakeBiconcave(rbcRef, radEqv)
      call RBC_ComputeGeometry(rbcRef)

      rbcRef => rbcRefs(2)
      call RBC_Create(rbcRef, nlat0)
      call RBC_MakeSphere(rbcRef, radEqv)
      call RBC_ComputeGeometry(rbcRef)

    end if

    ! Wall periodic boundary condition
    do iwall = 1, nwall
      wall => walls(iwall)
      call Wall_Build_V2V(wall, Lb)
    end do ! iwall

    ! Mechanical properties of cells and walls
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      select case(rbc%celltype)
      case(1)
       rbc%ES = 12.4
       rbc%ED = 200.
       rbc%EB = 6.69D-2
      case(2)
       print *,"CASE 2 --- celltype"
       rbc%ES = 10.
       rbc%ED = 50.
       rbc%EB = 6.D-2
      case default
       stop "bad cellcase"
      end select
    end do ! irbc

    do iwall = 1, nwall
      wall => walls(iwall)
    end do ! iwall

    vbkg(1:2) = 0. ; vbkg(3) = 8.
  end subroutine InitSystem

  subroutine FinalizeSystem

    integer :: irbc, iwall, ierr
    type(t_Rbc),pointer :: rbc
    type(t_Wall),pointer :: wall

    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      call RbcPolarPatch_Destroy(rbc%patch)
      call RBC_Destroy(rbc)
    end do ! irbc

    do iwall = 1, nwall
      wall => walls(iwall)
      call MatDestroy(wall%lhs,ierr)
      print*, 'MatDestroy ierr =', ierr
      call Wall_Destroy(wall)
    end do ! iwall

    if (nrbc > 0) deallocate(rbcs)
    if (nwall > 0) deallocate(walls)

    call RbcPolarPatch_Destroy(rbcRefs(1)%patch)
    call RbcPolarPatch_Destroy(rbcRefs(2)%patch)
    call RBC_Destroy(rbcRefs(1))
    call RBC_Destroy(rbcRefs(2))

    deallocate(rbcRefs)

  end subroutine FinalizeSystem


  subroutine ReboxRbcs

    integer :: irbc, ii
    type(t_rbc),pointer :: rbc
    real(WP) :: xc
    real :: iLb(3)

    iLb(:) = 1/Lb(:)

    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      do ii = 1, 3
        xc = 0.5*(maxval(rbc%x(:,:,ii)) + minval(rbc%x(:,:,ii)))
        rbc%x(:,:,ii) = rbc%x(:,:,ii) - floor(xc*iLb(ii))*Lb(ii)
      end do ! ii
    end do ! irbc

  end subroutine ReboxRbcs

end program InitCond
