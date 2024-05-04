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
    integer :: npoint
    character(CHRLEN) :: fn
    integer, parameter :: funit = 90
    integer :: ierr

    integer  :: Nx, Ny, Nz
    real(WP), parameter :: dx0 = 0.3
    real(WP), dimension(3) :: vl

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

    Nx = NINT(Lb(1)/dx0)
    Ny = NINT(Lb(2)/dx0)
    Nz = NINT(Lb(3)/dx0)

    npoint = Nx*Ny*Nz
    allocate (xs(npoint, 3), vs(npoint, 3))

    ! Evolve cells
    call Compute_Rbc_Vel

    ! Enforce no-slip condition on the wall
    call NoSlipWall

    call MakeMesh(Nx, Ny, Nz, xs)

    print *, "A"
    call TargetList_CreateFromRaw(tlist, xs)

    print *, "B"
    call CalcVelocityField(tlist, vs)

    call WriteVelCsv(npoint, vs, xs)

    call TargetList_Destroy(tlist)
    deallocate (xs, vs)

  END SUBROUTINE GetVelocity

  subroutine WriteVelCsv(npoint, vels, pts)
    integer :: npoint
    real(WP), dimension(npoint, 3) :: vels, pts

    integer :: i

    fn = "./D/field.csv"
    open (1, file=fn)
    write (1, *) "X,Y,Z,Vx,Vy,Vz"
    do i = 1, npoint
      write (1, *) pts(i, 1), ",", pts(i, 2), ",", pts(i, 3), ",", vels(i, 1), ",", vels(i, 2), ",", vels(i, 3)
    end do !i
    close (1)
  end subroutine writeVelCsv

  subroutine MakeMesh(Nx, Ny, Nz, pts)
    integer :: Nx, Ny, Nz
    real(WP), dimension(:, :) :: pts ! has dimensions (Nx * Ny * Nz , 3)

    real(WP) :: x, y, z, dx, dy, dz
    integer :: x_it, y_it, z_it, p_it

    dx = Lb(1)/REAL(Nx); 
    dy = Lb(2)/REAL(Ny); 
    dz = Lb(3)/REAL(Nz); 
    p_it = 1
    do x_it = 1, Nx
      x = 0.5*dx + REAL(x_it - 1)*dx
      do y_it = 1, Ny
        y = 0.5*dy + REAL(y_it - 1)*dy
        do z_it = 1, Nz
          z = 0.5*dz + REAL(z_it - 1)*dz

          pts(p_it, :) = (/x, y, z/)
          p_it = p_it + 1

        end do !z_it
      end do !y_it
    end do !x_it
  end subroutine MakeMesh

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
    allocate (rbcRefs(2))

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
      select case (rbc%celltype)
      case (1)
        rbc%ES = 12.4
        rbc%ED = 200.
        rbc%EB = 6.69D-2
      case (2)
        print *, "CASE 2 --- celltype"
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

    ! Background velocity
    !    if (Nt0 == 0) then
!      vbkg(1:2) = 0.
!      vbkg(3) = 10.
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
