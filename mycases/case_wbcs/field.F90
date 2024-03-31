! cells in cylindrical tubes
program cells_in_tube

  use ModDataTypes
  use ModDataStruct
  use ModRBC
  use ModWall
  use ModConf
  use ModData
  use ModTimeInt
  use ModIO
  use ModBasicMath
  use ModQuadRule
  use ModPostProcess
  use ModNoSlip
  use ModIntOnWalls
  use ModVelSolver
  use ModRepulsion

  implicit none
  character(CHRLEN) :: fn

! #include "../../petsc_include.h"

  call InitAll

  call GetVelocity

  call FinalizeAll

  stop

contains

  SUBROUTINE GetVelocity

    integer :: irbc, ii, iwall
    type(t_RBC), pointer :: rbc
    type(t_TargetList) :: tlist
    real(WP), allocatable :: xs(:, :), vs(:, :)
    real(WP), allocatable :: x(:, :, :), v(:, :, :)
    real(WP) :: xmin(3), xmax(3), szcell(3)
    integer ::  N1, N2, N3, npoint, i1, i2, i3, idx
    character(CHRLEN) :: fn
    integer, parameter :: funit = 90
    integer :: ierr

    real(WP) :: Xo
    integer  :: Ny, Nz
    real(WP), parameter :: dx0 = 0.3
    real(WP), dimension(3) :: vl

    !======================================================================
    ! Set up background velocity and RBC surface density
    vBkg(1) = 0.
    vBkg(2) = 0.
    vBkg(3) = 8.
    call MPI_Bcast(vBkg, 3, MPI_WP, 0, MPI_COMM_WORLD, ierr); 
    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      call RBC_ComputeGeometry(rbc)

      rbc%f = 0.0
      do ii = 1, 3
        rbc%g(:, :, ii) = vBkg(ii)
      end do

      call Rbc_BuildSurfaceSource(rbc, xFlag=.true., fFlag=.true., gFlag=.true.)
    end do ! irbc
    call SourceList_UpdateCoord(slist_rbc, rbcs)
    call SourceList_UpdateDensity(slist_rbc, rbcs, UpdateF=.true., UpdateG=.true.)

    if (nwall > 0) then
      do iwall = 1, nwall
        call Wall_ComputeGeometry(walls(iwall))
      end do ! iwall

      ! The source list and target list only need to be created and updated once
      call SourceList_UpdateCoord(slist_wall, walls=walls)
      call TargetList_Update(tlist_wall, walls=walls)

      if (PhysEwald) then
        do iwall = 1, nwall
          call PrepareSingIntOnWall(walls(iwall))
        end do ! iwall
      end if
    end if

    Xo = Lb(1)/2.
    Ny = NINT(Lb(2)/dx0); print *, "Ny = ", Ny
    Nz = NINT(Lb(3)/dx0); print *, "Nz = ", Nz
    allocate (x(Ny, Nz, 3))

    call MakeMesh(Xo, Ny, Nz, x)
    npoint = Ny*Nz
    allocate (xs(npoint, 3), vs(npoint, 3))

    ! Evolve cells
    call Compute_Rbc_Vel

    ! Enforce no-slip condition on the wall
    call NoSlipWall

    ! call leukVel(rbcs(20), vl)

    xs = RESHAPE(x, SHAPE(xs))

    print *, "A"
    call TargetList_CreateFromRaw(tlist, xs)

    print *, "B"
    call CalcVelocityField(tlist, vs)

    allocate (v(Ny, Nz, 3))

    v = RESHAPE(vs, SHAPE(v))
    v(:, :, 3) = v(:, :, 3) - vBkg(3)
    !    v(:,:,:,3) = v(:,:,:,3) - vl(3)  ! substract leuk velocity

    call writeVel(Ny, Nz, v, x)

    call TargetList_Destroy(tlist)
    deallocate (xs, vs)

  END SUBROUTINE GetVelocity

  SUBROUTINE writeVel(Ny, Nz, v, x)
    integer                          :: Ny, Nz
    real(WP), dimension(Ny, Nz, 3)     :: v, x

    integer   :: i, j, k, m
    ! can reformat file by changing these
    open (1, file='v.f')
    write (1, *) Ny, Nz, 3
    write (1, *) v
    close (1)

    open (1, file='v.g')
    write (1, *) Ny, Nz
    write (1, *) x(:, :, 2:3)
    close (1)

  END SUBROUTINE writeVel

  SUBROUTINE MakeMesh(Xo, Ny, Nz, x)
    real(WP)                         :: Xo
    integer                          :: Ny, Nz
    real(WP), dimension(:, :, :)      :: x

    real(WP)  :: y, z, dy, dz
    integer   :: i, j, k

    dy = Lb(2)/REAL(Ny); print *, "dy = ", dy
    dz = Lb(3)/REAL(Nz); print *, "dz = ", dz

    do k = 1, Nz
      z = 0.5*dz + REAL(k - 1)*dz
      do j = 1, Ny
        y = 0.5*dy + REAL(j - 1)*dy
        x(j, k, 1) = Lb(1)/2
        x(j, k, 2) = y
        x(j, k, 3) = z
        !  do k = 1,Nz
        !     z = 0.5*dz+REAL(k-1)*dz
        !     do j = 1,Ny
        !        y = 0.5*dy+REAL(j-1)*dy
        !       !  x(j,k,1) = Lb(1)/2.
        !       !  x(j,k,2) = y
        !        x(j,k,1) = x
        !        x(j,k,3) = z

        !  print *, "COORDINATE", x(j, k, 1), x(j, k, 2), x(j, k, 3)
      end do
    end do

    open (1, file='v.g')
    write (1, *) Ny, Nz
    write (1, *) x(:, :, 2:3)
    close (1)

  END SUBROUTINE MakeMesh

  SUBROUTINE leukVel(rbc, vl)
    type(t_RBC)            :: rbc
    real(WP), dimension(3) :: vl

    real(WP)               :: ds
    integer                :: ilat, ilon

    call RBC_ComputeGeometry(rbc)

    vl = 0.
    do ilat = 1, rbc%nlat
      do ilon = 1, rbc%nlon
        ds = rbc%detj(ilat, ilon)*rbc%w(ilat)
        vl = vl + ds*rbc%g(ilat, ilon, :)
      end do ! ilon
    end do ! ilat
    vl = vl/rbc%area

    print *, "LEUK VEL", vl

  END SUBROUTINE leukVel

  !**********************************************************************
  subroutine InitAll

    ! System intialization
    call InitMPI()
    call GaussQuad_Init
    call ReadConfig('Input/tube.in')
    call InitSystem

    ! Note: SetEwaldPrms and DomainDecomp must be called after InitSystem
    call SetEwaldPrms
    call DomainDecomp
    call GlobData_Init

    ! Prepare for time integration
    call IO_Init
    call TimeInt_Init

  end subroutine InitAll

  !**********************************************************************
  subroutine FinalizeAll

    call IO_Finalize
    call TimeInt_Finalize

    call GlobData_Finalize
    call FinalizeSystem

    call GaussQuad_Finalize
    call FinalizeMPI

  end subroutine FinalizeAll

  !**********************************************************************
  subroutine InitSystem

    integer :: irbc, iwall
    type(t_Rbc), pointer :: rbc, rbcRef
    type(t_Wall), pointer :: wall
    integer :: nlat0
    real(WP) :: radEqv
    integer :: ierr

    call ReadRestart(restart_file)

    ! Reference cells
    allocate (rbcRefs(3))

    if (nrbc > 0) then
      radEqv = 1.
      nlat0 = rbcs(1)%nlat0

      rbcRef => rbcRefs(1)
      call RBC_Create(rbcRef, nlat0)
      call RBC_MakeBiconcave(rbcRef, radEqv)
      call RBC_ComputeGeometry(rbcRef)

      rbcRef => rbcRefs(2)
      call RBC_Create(rbcRef, nlat0)
      call RBC_MakeLeukocyte(rbcRef, radEqv)
      call RBC_ComputeGeometry(rbcRef)

      rbcRef => rbcRefs(3)
      call ImportReadRBC('Input/SickleCell.dat', rbcRef)
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
      select case (rbc%celltype)
      case (1)
        rbc%ES = 12.4
        rbc%ED = 200.
        rbc%EB = 6.69D-2
      case (2)
        ! print *,"CASE 2 --- celltype"
        rbc%ES = 887
        rbc%ED = 200.
        rbc%EB = 2.44D-2
      case (3)
        rbc%ES = 12.4*20/7.1
        rbc%ED = 200*49.4/15.4
        rbc%EB = 6.69D-2*19.5/5.7
      case default
        stop "bad cellcase"
      end select
    end do ! irbc

    do iwall = 1, nwall
      wall => walls(iwall)
    end do ! iwall

    ! Background velocity
    !    if (Nt0 == 0) then
    vbkg(1:2) = 0.
    vbkg(3) = 8.
    !    end if
    print *, vbkg

  end subroutine InitSystem

  !**********************************************************************
  subroutine FinalizeSystem

    integer :: irbc, iwall
    type(t_Rbc), pointer :: rbc
    type(t_Wall), pointer :: wall

    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      call RBC_Destroy(rbc)
    end do ! irbc

    do iwall = 1, nwall
      wall => walls(iwall)
      call Wall_Destroy(wall)
    end do ! iwall

    if (nrbc > 0) deallocate (rbcs)
    if (nwall > 0) deallocate (walls)

  end subroutine FinalizeSystem

end program
