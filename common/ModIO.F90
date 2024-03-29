module ModIO

  use MPI
  use NETCDF    ! For reading wall mesh
  use ModDataTypes
  use ModDataStruct
  use ModConf
  use ModData
  use ModRBC
  use ModWall
  use ModSphpk
  use ModPostProcess

  implicit none

  private

  public :: IO_Init, &
            IO_Finalize, &
            WriteAll, &
            WriteRBC, &
            WriteManyRBCs, &
            WriteManyRBCsNoFilt, &
            WriteWall, &
            WriteManyWalls, &
            WritePressGrad, &
            WriteFlowRate, &
            !   WriteVelField, &
            !   WriteViscStress, &
            ReadMyWallMesh, &
            ReadWallMesh, &
            WriteRestart, &
            ReadRestart, &
            ReadRestart_NoWalls, &
            ExportWriteRBC, &
            ImportReadRBC, &
            WriteManyRBCsByType

contains

!**********************************************************************
  subroutine IO_Init

    character(CHRLEN) :: fn

    if (rootWorld) then
      if (pgrad_out > 0) then
        write (fn, FMT=fn_FMT), 'D/', 'pgrad', Nt0, '.dat'
        open (pGrad_unit, file=trim(fn), action='write')
      end if

      if (flow_out > 0) then
        write (fn, FMT=fn_FMT), 'D/', 'flow', Nt0, '.dat'
        open (flow_unit, file=trim(fn), action='write')
      end if

      if (ftot_out > 0) then
        write (fn, FMT=fn_FMT), 'D/', 'ftot', Nt0, '.dat'
        open (ftot_unit, file=trim(fn), action='write')
      end if

    end if

  end subroutine IO_Init

!**********************************************************************
  subroutine IO_Finalize

    logical :: file_opened

    if (rootWorld) then
      inquire (pgrad_unit, opened=file_opened)
      if (file_opened) close (pgrad_unit)

      inquire (flow_unit, opened=file_opened)
      if (file_opened) close (flow_unit)
    end if

  end subroutine IO_Finalize

!**********************************************************************
! Write everything
! Arguments:
!  lt -- time step
!  time -- current time
  subroutine WriteAll(lt, time)
    integer :: lt
    real(WP) :: time

    character(CHRLEN) :: fn

    if (rootWorld) then
      if (cell_out > 0 .and. mod(lt, cell_out) == 0) then
        write (fn, FMT=fn_FMT) 'D/', 'x', lt, '.dat'
        call WriteManyRBCs(fn, nrbc, rbcs)
        write (fn, FMT=fn_FMT) 'D/', 'xe', lt, '.dat'
        call WriteExactPts(fn, nrbc, rbcs)

        ! Comment these 9 lines if you're only generating cells of 1 type
        ! Write out only type-1 cells (healthy RBCs)
        ! write (fn, FMT=fn_FMT) 'D/', '1x', lt, '.dat'
        ! call WriteManyRBCsByType(fn, nrbc, rbcs, 1)
        ! ! Write out only type-2 cells (WBCs)
        ! write (fn, FMT=fn_FMT) 'D/', '2x', lt, '.dat'
        ! call WriteManyRBCsByType(fn, nrbc, rbcs, 2)
        ! ! Write out only type-3 cells (sickle cells)
        ! write (fn, FMT=fn_FMT) 'D/', '3x', lt, '.dat'
        ! call WriteManyRBCsByType(fn, nrbc, rbcs, 3)
      end if

      if (wall_out > 0 .and. mod(lt, wall_out) == 0) then
        write (fn, FMT=fn_FMT) 'D/', 'wall', lt, '.dat'
        call WriteManyWalls(fn, nwall, walls)
      end if

      if (pgrad_out > 0 .and. mod(lt, pGrad_out) == 0) then
        call WritePressGrad(lt, time)
      end if

      if (flow_out > 0 .and. mod(lt, flow_out) == 0) then
        call WriteFlowRate(lt, time)
      end if

      if (restart_out > 0 .and. mod(lt, restart_out) == 0) then
        write (fn, FMT=fn_FMT) 'D/', 'restart', lt, '.dat'
        call WriteRestart(fn, lt, time)
        fn = 'D/restart.LATEST.dat'
        call WriteRestart(fn, lt, time)
      end if
    end if

  end subroutine WriteAll

!**********************************************************************
! Write a rbc mesh to file
  subroutine WriteRBC(fn, rbc)
    character(*) :: fn
    type(t_rbc) :: rbc

    real(WP), dimension(:, :, :), allocatable :: x, xa, xb
    integer :: nlat, nlon, ilat, ilon

    nlat = rbc%nlat
    nlon = rbc%nlon

    ! Allocate working array
    allocate (x(0:nlat, nlon, 3), xa(nlon/2 + 1, nlat + 1, 3), xb(nlon/2 + 1, nlat + 1, 3))

    ! Convert the mesh from Gaussian to uniform
    call ShAnalGau(nlat, nlon, 3, rbc%x, size(rbc%x, 1), size(rbc%x, 2), &
                   xa, xb, size(xa, 1), size(xa, 2), rbc%wshags)
    call ShSynthEqu(nlat + 1, nlon, 3, x, size(x, 1), size(x, 2), &
                    xa, xb, size(xa, 1), size(xa, 2), rbc%wshses)

    open (cell_unit, file=trim(fn), action='write')
    write (cell_unit, '(A)') 'VARIABLES = X, Y, Z'
    write (cell_unit, '(A,I9,A,I9,A)') 'ZONE I=', nlat + 1, ' J=', nlon + 1, ' F=POINT'
    do ilon = 1, nlon
    do ilat = 0, nlat
      write (cell_unit, '(3F20.10)') x(ilat, ilon, :)
    end do ! ilat
    end do ! ilon
    do ilat = 0, nlat
      write (cell_unit, '(3F20.10)') x(ilat, 1, :)
    end do ! ilat
    close (cell_unit)

    ! Deallocate working arrays
    deallocate (x, xa, xb)

  end subroutine WriteRBC

!**********************************************************************
! Write the shape of many red blood cells to file
! Arguments:
!  fn -- file name
!  nrbc -- number of cells
!  rbcs -- red blood cells
  subroutine WriteManyRBCs(fn, nrbc, rbcs)
    character(*) :: fn
    integer :: nrbc
    type(t_rbc), target :: rbcs(:)

    type(t_Rbc), pointer :: rbc
    real(WP), dimension(:, :, :), allocatable :: x, xa, xb
    integer :: irbc, nlat, nlon, ilat, ilon

    ! Argument check
    if (nrbc <= 0) return

    ! Write file head
    open (cell_unit, file=trim(fn), action='write')
    write (cell_unit, '(A)') 'VARIABLES = X, Y, Z'

    ! Write the record of every cell
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      nlat = rbc%nlat; nlon = rbc%nlon

      ! Allocate working arrays
      allocate (x(0:nlat, nlon, 3), xa(nlon/2 + 1, nlat + 1, 3), xb(nlon/2 + 1, nlat + 1, 3))

      call ShAnalGau(nlat, nlon, 3, rbc%x, size(rbc%x, 1), size(rbc%x, 2), &
                     xa, xb, size(xa, 1), size(xa, 2), rbc%wshags)
      call ShFilter(nlat, nlon, 3, xa, xb, size(xa, 1), size(xa, 2), rbc%nlat0, rbc%nlon0)
      call ShSynthEqu(nlat + 1, nlon, 3, x, size(x, 1), size(x, 2), &
                      xa, xb, size(xa, 1), size(xa, 2), rbc%wshses)

      write (cell_unit, '(A,I9,A,I9,A)') 'ZONE I=', nlat + 1, '  J=', nlon + 1, '  F=POINT'

      do ilon = 1, nlon
      do ilat = 0, nlat
        write (cell_unit, '(3F20.10)') x(ilat, ilon, :)
      end do ! ilat
      end do ! ilon
      do ilat = 0, nlat
        write (cell_unit, '(3F20.10)') x(ilat, 1, :)
      end do ! ilat

      ! Deallocate working arrays
      deallocate (x, xa, xb)
    end do ! irbc

    close (cell_unit)

  end subroutine WriteManyRBCs

!*********************************************************************
! Write the shape of blood cells of specified type to file
! Arguments:
!  fn -- file suffix name
!  nrbc -- nubmer of cells
!  rbcs -- blood cells
!  type -- type filter (1: rbc, 2: leukocyte, 3: sickle cell)
  subroutine WriteManyRBCsByType(fn, nrbc, rbcs, type)
    character(*) :: fn
    integer :: nrbc
    type(t_rbc), target :: rbcs(:)
    integer :: type

    type(t_Rbc), pointer :: rbc
    real(WP), dimension(:, :, :), allocatable :: x, xa, xb
    integer :: irbc, nlat, nlon, ilat, ilon

    ! Argument check
    if (nrbc <= 0) return

    ! Write file head
    open (cell_unit, file=trim(fn), action='write')
    write (cell_unit, '(A)') 'VARIABLES = X, Y, Z'

    ! Write the record of every cell
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      nlat = rbc%nlat; nlon = rbc%nlon

      if (rbc%celltype .ne. type) then
        cycle
      end if

      ! Allocate working arrays
      allocate (x(0:nlat, nlon, 3), xa(nlon/2 + 1, nlat + 1, 3), xb(nlon/2 + 1, nlat + 1, 3))

      call ShAnalGau(nlat, nlon, 3, rbc%x, size(rbc%x, 1), size(rbc%x, 2), &
                     xa, xb, size(xa, 1), size(xa, 2), rbc%wshags)
      call ShFilter(nlat, nlon, 3, xa, xb, size(xa, 1), size(xa, 2), rbc%nlat0, rbc%nlon0)
      call ShSynthEqu(nlat + 1, nlon, 3, x, size(x, 1), size(x, 2), &
                      xa, xb, size(xa, 1), size(xa, 2), rbc%wshses)

      write (cell_unit, '(A,I9,A,I9,A)') 'ZONE I=', nlat + 1, '  J=', nlon + 1, '  F=POINT'

      do ilon = 1, nlon
      do ilat = 0, nlat
        write (cell_unit, '(3F20.10)') x(ilat, ilon, :)
      end do ! ilat
      end do ! ilon
      do ilat = 0, nlat
        write (cell_unit, '(3F20.10)') x(ilat, 1, :)
      end do ! ilat

      ! Deallocate working arrays
      deallocate (x, xa, xb)
    end do ! irbc

    close (cell_unit)

  end subroutine WriteManyRBCsByType

  subroutine WriteManyRBCsNoFilt(fn, nrbc, rbcs)
    character(*) :: fn
    integer :: nrbc
    type(t_rbc), target :: rbcs(:)

    type(t_Rbc), pointer :: rbc
    !real(WP),dimension(:,:,:),allocatable :: x, xa, xb
    integer :: irbc, nlat, nlon, ilat, ilon

    ! Argument check
    if (nrbc <= 0) return

    ! Write file head
    open (cell_unit, file=trim(fn), action='write')
    write (cell_unit, '(A)') 'VARIABLES = X, Y, Z'

    ! Write the record of every cell
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      nlat = rbc%nlat; nlon = rbc%nlon

      write (cell_unit, '(A,I9,A,I9,A)') 'ZONE I=', nlat, '  J=', nlon + 1, '  F=POINT'

      do ilon = 1, nlon
      do ilat = 1, nlat
        write (cell_unit, '(3F20.10)') rbc%x(ilat, ilon, :)
      end do ! ilat
      end do ! ilon
      do ilat = 1, nlat
        write (cell_unit, '(3F20.10)') rbc%x(ilat, 1, :)
      end do ! ilat

    end do ! irbc

    close (cell_unit)

  end subroutine WriteManyRBCsNoFilt

  subroutine WriteExactPts(fn, nrbc, rbcs)
    character(*) :: fn
    integer :: nrbc
    type(t_rbc), target :: rbcs(:)

    type(t_Rbc), pointer :: rbc
    real(WP), dimension(:, :, :), allocatable :: x, xa, xb
    integer :: irbc, nlat, nlon, ilat, ilon, iz

    ! Argument check
    if (nrbc <= 0) return

    ! Write file head
    open (cell_unit, file=trim(fn), action='write')
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      nlat = rbc%nlat; nlon = rbc%nlon

      do ilon = 1, nlon
      do ilat = 1, nlat
      do iz = 1, 3
        write (cell_unit, '(F25.15)') rbc%x(ilat, ilon, iz)
      end do
      end do
      end do
    end do

    close (cell_unit)

  end subroutine WriteExactPts

!**********************************************************************
  subroutine WriteWall(fn, wall)
    character(*) :: fn
    type(t_wall) :: wall

    integer :: ivert, iele

    open (wall_unit, file=trim(fn), action='write')
    write (wall_unit, '(A)') 'VARIABLES = X, Y, Z'
    write (wall_unit, '(A,I9,A,I9,A)') &
      'ZONE N = ', wall%nvert, ' E = ', wall%nele, ' F=FEPOINT ET=TRIANGLE'

    do ivert = 1, wall%nvert
      write (wall_unit, '(3F20.10)') wall%x(ivert, :)
    end do ! ivert

    do iele = 1, wall%nele
      write (wall_unit, '(3I9)') wall%e2v(iele, :)
    end do ! iele

    close (wall_unit)

  end subroutine WriteWall

!**********************************************************************
! Write walls
! Arguments:
!  fn -- file name
!  nwall -- number of walls
!  walls --
  subroutine WriteManyWalls(fn, nwall, walls)
    character(*) :: fn
    integer :: nwall
    type(t_wall), target :: walls(:)

    integer :: iwall, nvert, nele, ivert, iele
    type(t_wall), pointer :: wall

    ! Argument check
    if (nwall <= 0) return

    open (wall_unit, file=trim(fn), action='write')
    write (wall_unit, '(A)') 'VARIABLES = X, Y, Z'

    do iwall = 1, nwall
      wall => walls(iwall)

      nvert = wall%nvert
      nele = wall%nele

      write (wall_unit, '(A,I9,A,I9,A)') &
        'ZONE N = ', nvert, ' E = ', nele, ' F=FEPOINT ET=TRIANGLE'

      do ivert = 1, nvert
        write (wall_unit, '(3F20.10)') wall%x(ivert, :)
      end do ! ivert

      do iele = 1, nele
        write (wall_unit, '(3I9)') wall%e2v(iele, :)
      end do ! iele
    end do ! iwall

    close (wall_unit)

  end subroutine WriteManyWalls

!**********************************************************************
! Write pressure gradient
! Arguments:
!  lt -- time step
!  time -- current time
  subroutine WritePressGrad(lt, time)
    integer :: lt
    real(WP) :: time

    real(WP) :: pGrad(3)

    pgrad = WallShearForce()
    pgrad = pGrad*product(iLb)

    write (pGrad_unit, '(F15.5, 3ES15.5)') time, pGrad

  end subroutine WritePressGrad

!**********************************************************************
! Write cellular flow rate
! Arguments:
!  lt -- time step
!  time -- current time
  subroutine WriteFlowRate(lt, time)
    integer :: lt
    real(WP) :: time

    real(WP) :: flow

    flow = CellFlowRate()

    write (flow_unit, '(F15.5, ES15.5)') time, flow

  end subroutine WriteFlowRate

!!**********************************************************************
!! Write velocity field
!! Argument:
!!  fn -- file name
!  subroutine WriteVelField(fn)
!    character(*) :: fn
!
!    real(WP) :: rmax, xc(3), r, th
!    real(WP) :: C1, C2
!    integer,parameter :: Nx = 30, Ny = 30
!    real(WP),allocatable :: x(:,:), v(:,:)
!    integer :: i, j, n
!
!    ! Allocate working arrays
!    allocate(x(Nx*Ny,3), v(Nx*Ny,3))
!
!    !======================================================================
!    ! Set up mesh mesh
!    ! A circular domain
!    rmax = 1.
!    xc = (/ 0.5*Lb(1), 0.5*Lb(2), 0.5*Lb(3) /)
!
!    do j = 1, Ny
!    do i = 1, Nx
!      r = (i - 1)*rmax/real(Nx - 1)
!      th = (j - 1)*two_Pi/real(Ny - 1)
!
!      n = i + (j-1)*Nx
!      x(n,1) = xc(1) + r*cos(th)
!      x(n,2) = xc(2) + r*sin(th)
!      x(n,3) = xc(3)
!    end do ! i
!    end do ! j
!
!!    ! A rectangular domain
!!    do j = 1, Ny
!!    do i = 1, Nx
!!      n = i + (j-1)*Nx
!!      x(n,1) = 0.25 + real(i-1)/(Nx-1)
!!      x(n,2) = 0.25 + real(j-1)/(Ny-1)
!!      x(n,3) = 0.5*Lb(3)
!!    end do ! i
!!    end do ! j
!    !======================================================================
!
!    C1 = 1/(8*PI)
!    C2 = (1 - viscRat)/(8*PI)
!    call VelocityField(C1, C2, x, v)
!
!    if (rootWorld) then
!      ! Write the velocity field
!      open(vel_unit, file=trim(fn), action='write')
!
!      write(vel_unit, '(A)') 'VARIABLES = X, Y, Z, U, V, W'
!      write(vel_unit, '(A,I0,A,I9,A)') 'ZONE I=', Nx, ' J=', Ny, ' F=POINT'
!      do j = 1, Ny
!      do i = 1, Nx
!   n = i + Nx*(j-1)
!   write(vel_unit, '(6ES15.5)') x(n,:), v(n,:)
!      end do ! j
!      end do ! i
!
!      close(vel_unit)
!    end if
!
!    ! Deallocate working arrays
!    deallocate(x, v)
!
!  end subroutine WriteVelField
!
!!**********************************************************************
!! Output the extra viscous stress due to RBCs
!  subroutine WriteViscStress(lt, time)
!    integer :: lt
!    real(WP) :: time
!
!    real(WP) :: tau(3,3)
!
!    call ComputeParticleStress(tau)
!    tau = tau*product(iLb)
!
!    write(visc_unit, '(F15.5,2ES15.5)') time, tau(1,3), tau(3,1)
!
!  end subroutine WriteViscStress
!
!**********************************************************************
! Read a wall mesh file in exodus-II finite-element format
  subroutine ReadWallMesh(fn, wall)
    character(*) :: fn
    type(t_wall) :: wall

    integer :: ncid, dimid, varid, ierr
    integer :: nvert, nele, num_nod_per_el
    integer, allocatable :: connect(:, :)
    character(*), parameter :: func_name = "ReadWallMesh"

    ! Check whether the file exists
    ierr = NF90_OPEN(trim(fn), NF90_NOWRITE, ncid)
    if (ierr .ne. 0) then
      write (*, *) 'Subroutine ', func_name
      write (*, *) 'Error: file ', trim(fn), ' does not exist'
      stop
    end if

    ! Check whether the mesh is Tri3
    ierr = nf90_inq_dimid(ncid, "num_nod_per_el1", dimid)
    ierr = nf90_inquire_dimension(ncid, dimid, len=num_nod_per_el)
    if (num_nod_per_el .ne. 3) then
      write (*, *) 'Subroutine ', func_name
      write (*, *) 'Error: the input mesh is not of Tri3 type'
      stop
    end if

    ! Find the size of the mesh
    ierr = NF90_INQ_DIMID(ncid, "num_nodes", dimid)
    ierr = NF90_INQUIRE_DIMENSION(ncid, dimid, len=nvert)
    ierr = NF90_INQ_DIMID(ncid, "num_elem", dimid)
    ierr = NF90_INQUIRE_DIMENSION(ncid, dimid, len=nele)
    call Wall_Create(wall, nvert, nele)

    ! Read node coordinates
    ! ierr = NF90_INQ_VARID(ncid, "coord", varid)
    ! ierr = NF90_GET_VAR(ncid, varid, wall%x)

    ierr = NF90_INQ_VARID(ncid, "coordx", varid)
    ierr = NF90_GET_VAR(ncid, varid, wall%x(:, 1))

    ierr = NF90_INQ_VARID(ncid, "coordy", varid)
    ierr = NF90_GET_VAR(ncid, varid, wall%x(:, 2))

    ierr = NF90_INQ_VARID(ncid, "coordz", varid)
    ierr = NF90_GET_VAR(ncid, varid, wall%x(:, 3))

    ! Read mesh connectivity
    allocate (connect(3, nele))

    ierr = NF90_INQ_VARID(ncid, "connect1", varid)
    ierr = NF90_GET_VAR(ncid, varid, connect)
    wall%e2v = transpose(connect)

    deallocate (connect)

    ! Close file
    ierr = nf90_close(ncid)

  end subroutine ReadWallMesh

!*********************************************************************
  ! Given a cell (rbc) and a file name (fn), this writes the rbc's
  ! plain coordinate values, celltype, and resolution (nlat & nlon) to the external file.
  ! This file can later be reread into a later simulation to replicate the cell-shape of the inputted rbc.
  ! This method was used after inducing an "UncurvedSickleSpheroid" (see ModRbc RBC_MakeUncurvedSickleSpheroid)
  ! inside a flow, so that we can export the final curved sickle cell shape to a file for later usage.
  subroutine ExportWriteRBC(fn, rbc)
    character(*) :: fn
    type(t_rbc) :: rbc

    open (cell_unit, file=trim(fn), action='write')

    write (cell_unit, *) rbc%nlat0, rbc%nlon0
    write (cell_unit, *) rbc%nlat, rbc%nlon
    write (cell_unit, *) rbc%celltype
    write (cell_unit, *) rbc%x

    close (cell_unit)
  end subroutine ExportWriteRBC

  ! Works in-tandem with ModIO ExportWriteRBC:
  ! Given a cell (rbc), file name (fn), and optionally center (xc), import the data from the file
  ! into the cell. The cell will be created with the resolution defined from the file data (and cannot be modified).
  ! Used to import Sickle Cell model from the input file.
  subroutine ImportReadRBC(fn, rbc, xc)
    character(*) :: fn
    type(t_rbc) :: rbc
    Real(WP), optional :: xc(3)

    integer :: nlat0, nlon0, nlat, nlon, dealias_fac, celltype, ilat, ilon, ierr, ii
    if (rootWorld) then
      open (unit=cell_unit, file=trim(fn), status='old', action='read')

      read (cell_unit, *) nlat0, nlon0
      read (cell_unit, *) nlat, nlon
      read (cell_unit, *) celltype
    end if !root world

    call MPI_Bcast(nlat0, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(nlon0, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(nlat, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(nlon, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(celltype, 1, MPI_Integer, 0, MPI_Comm_World, ierr)

    if (nlat/real(nlat0) .lt. 1.5) then
      dealias_fac = 100
    else
      dealias_fac = nlat/nlat0
    end if !dealias factor

    rbc%celltype = celltype
    call Rbc_Create(rbc, nlat0, dealias_fac)
    rbc%nlat = nlat
    rbc%nlon = nlon

    if (rootWorld) then
      read (cell_unit, *) rbc%x
    end if !root world
    close (cell_unit)
    call MPI_Bcast(rbc%x, size(rbc%x), MPI_WP, 0, MPI_Comm_World, ierr)

    ! recenter to zero - the exported cell was offset to fit in the box
    do ilon = 1, nlon
    do ilat = 1, nlat
      rbc%x(ilat, ilon, 1) = rbc%x(ilat, ilon, 1) - 5.25
      rbc%x(ilat, ilon, 2) = rbc%x(ilat, ilon, 2) - 5.25
      rbc%x(ilat, ilon, 3) = rbc%x(ilat, ilon, 3) - (5/7.0)
    end do
    end do

    if (present(xc)) then
      do ii = 1, 3
        rbc%x(:, :, ii) = rbc%x(:, :, ii) + xc(ii)
      end do ! ii
    end if

  end subroutine ImportReadRBC

!**********************************************************************
! Read a wall mesh file in dat format
  subroutine ReadMyWallMesh(fn, wall)
    character(*) :: fn
    type(t_wall) :: wall

    integer :: ncid, dimid, varid, ierr
    integer :: nvert, nele, num_nod_per_el
    integer, allocatable :: connect(:, :)
    character(*), parameter :: func_name = "ReadMyWallMesh"

    LOGICAL :: file_exists

    INQUIRE (FILE=trim(fn), EXIST=file_exists)

    open (unit=99, file=trim(fn), status='old', action='read')

    write (*, *) 'Mesh File: ', trim(fn)

    ! Check whether the file exists
    if (file_exists .neqv. .TRUE.) then
      write (*, *) 'Subroutine ', func_name
      write (*, *) 'Error: file ', trim(fn), ' does not exist'
      stop
    end if

    read (99, *) num_nod_per_el
    read (99, *) nvert
    read (99, *) nele

    ! Check whether the mesh is Tri3
    if (num_nod_per_el .ne. 3) then
      write (*, *) 'Subroutine ', func_name
      write (*, *) 'Error: the input mesh is not of Tri3 type'
      stop
    end if

    write (*, *) 'Number of Nodes per Element: ', num_nod_per_el
    write (*, *) 'Number of Nodes: ', nvert
    write (*, *) 'Number of Elements: ', nele

    call Wall_Create(wall, nvert, nele)

    read (99, *) wall%x
    read (99, *) wall%e2v

    ! Close file
    close (99)
    ierr = nf90_close(ncid)

  end subroutine ReadMyWallMesh

!**********************************************************************
! Write restart file
! Arguments:
!  lt --
!  time
  subroutine WriteRestart(fn, lt, time)
    character(*) :: fn
    integer :: lt
    real(WP) :: time

    integer :: irbc, iwall
    type(t_rbc), pointer :: rbc
    type(t_wall), pointer :: wall

    open (restart_unit, file=trim(fn), form='unformatted', action='write')

    write (restart_unit) Lb

    write (restart_unit) lt
    write (restart_unit) time
    write (restart_unit) vBkg

    ! Cells
    write (restart_unit) nrbc
    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      write (restart_unit) rbc%nlat0, rbc%nlon0
      write (restart_unit) rbc%nlat, rbc%nlon
      write (restart_unit) rbc%celltype
      write (restart_unit) rbc%starting_area
      write (restart_unit) rbcs(irbc)%x
    end do ! irbc

    ! Walls
    write (restart_unit) nwall
    do iwall = 1, nwall
      wall => walls(iwall)

      write (restart_unit) wall%nvert, wall%nele
      write (restart_unit) wall%x
      write (restart_unit) wall%f
      write (restart_unit) wall%e2v
    end do ! iwall

    close (restart_unit)

  end subroutine WriteRestart

!**********************************************************************
! Read restart file
! Arguments:
!  fn -- restart file name
!  Nt0 -- time step when the restart file is created
!  time -- time when the restart file is created
  subroutine ReadRestart(fn)
    character(*) :: fn

    integer :: irbc, iwall
    integer :: nlat0, nlon0, nlat, nlon, nvert, nele, celltype
    integer :: dealias_fac
    real(WP) :: starting_area
    type(t_Rbc), pointer :: rbc
    type(t_Wall), pointer :: wall
    integer :: ierr
    character(*), parameter :: func_name = 'ReadRestart'

    ! global parameters
    if (rootWorld) then
      write (*, *) 'Read restart file ', TRIM(fn)
      print *, 'error 0'
      open (restart_unit, file=trim(fn), form='unformatted', action='read')

      print *, 'error 1'
      read (restart_unit) Lb
      iLb = 1./Lb
      print *, 'error 2'
      print *, 'z'
      read (restart_unit) Nt0
      print *, 'a'
      read (restart_unit) time0
      print *, 'b'
      read (restart_unit) vBkg
      print *, '0'
      write (*, *) 'Lb = ', Lb
      write (*, *) 'Nt0 = ', Nt0
      write (*, *) 'time0 = ', time0
      write (*, *) 'vBkg = ', vBkg
    end if

    print *, '1'
    call MPI_Bcast(Lb, 3, MPI_WP, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(iLb, 3, MPI_WP, 0, MPI_Comm_World, ierr)

    call MPI_Bcast(Nt0, 1, MPI_INTEGER, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(time0, 1, MPI_WP, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(vBkg, 3, MPI_WP, 0, MPI_Comm_World, ierr)
    print *, '2'
    ! cells
    if (rootWorld) then
      read (restart_unit) nrbc

      write (*, *) 'nrbc = ', nrbc
    end if
    call MPI_Bcast(nrbc, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
    print *, '3'
    allocate (rbcs(nrbc))

    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      if (rootWorld) then
        read (restart_unit) nlat0, nlon0
        read (restart_unit) nlat, nlon
        ! celltype = 1; print *,"NO READ CELL TYPE"
        read (restart_unit) celltype
        read (restart_unit) starting_area
        write (*, *) 'irbc : ', irbc, ' nlat0 = ', nlat0, 'type = ', celltype, 'starting_area = ', starting_area
        write (*, *) 'nlat =', nlat
      end if
      call MPI_Bcast(nlat0, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
      call MPI_Bcast(nlon0, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
      call MPI_Bcast(nlat, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
      call MPI_Bcast(nlon, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
      call MPI_Bcast(celltype, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
      call MPI_Bcast(starting_area, 1, MPI_WP, 0, MPI_Comm_World, ierr)

      rbc%celltype = celltype
      rbc%starting_area = starting_area

      if (nlat/real(nlat0) .lt. 1.5) then
        dealias_fac = 100
      else
        dealias_fac = nlat/nlat0
      end if

      call Rbc_Create(rbc, nlat0, dealias_fac)

      ! Check array dimension
      if (rootWorld) then
        if (rbc%nlat .ne. nlat .or. rbc%nlon .ne. nlon) then
          print *, 'nlat,nlon =', nlat, nlon
          print *, 'rbc nlat,nlon =', rbc%nlat, rbc%nlon
          write (*, *) 'Subroutine ', func_name
          write (*, *) 'Error: invalid array dimension'
          stop
        end if
      end if

      if (rootWorld) then
        read (restart_unit) rbc%x
      end if
      call MPI_Bcast(rbc%x, size(rbc%x), MPI_WP, 0, MPI_Comm_World, ierr)
    end do ! irbc
    print *, '4'
    ! Walls
    if (rootWorld) then
      read (restart_unit) nwall

      write (*, *) 'nwall = ', nwall
    end if
    call MPI_Bcast(nwall, 1, MPI_Integer, 0, MPI_Comm_World, ierr)

    allocate (walls(nwall))

    do iwall = 1, nwall
      wall => walls(iwall)

      if (rootWorld) then
        read (restart_unit) nvert, nele

        write (*, *) 'iwall : ', iwall, ' nvert = ', nvert, ' nele = ', nele
      end if
      call MPI_Bcast(nvert, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
      call MPI_Bcast(nele, 1, MPI_Integer, 0, MPI_Comm_World, ierr)

      call Wall_Create(wall, nvert, nele)

      if (rootWorld) then
        read (restart_unit) wall%x
        read (restart_unit) wall%f
        read (restart_unit) wall%e2v
      end if
      call MPI_Bcast(wall%x, size(wall%x), MPI_WP, 0, MPI_Comm_World, ierr)
      call MPI_Bcast(wall%f, size(wall%f), MPI_WP, 0, MPI_Comm_World, ierr)
      call MPI_Bcast(wall%e2v, size(wall%e2v), MPI_Integer, 0, MPI_Comm_World, ierr)
    end do ! iwall

    if (rootWorld) then
      close (restart_unit)
      write (*, *)
    end if

  end subroutine ReadRestart

  subroutine ReadRestart_basecell(fn)
    character(*) :: fn

    integer :: irbc, iwall
    integer :: nlat0, nlon0, nlat, nlon, nvert, nele, celltype
    integer :: dealias_fac
    type(t_Rbc), pointer :: rbc
    type(t_Wall), pointer :: wall
    integer :: ierr
    character(*), parameter :: func_name = 'ReadRestart'

    if (rootWorld) then
      write (*, *) 'Read restart file ', TRIM(fn)
      open (restart_unit, file=trim(fn), form='unformatted', action='read')

      read (restart_unit) Lb
      iLb = 1./Lb
      read (restart_unit) Nt0
      read (restart_unit) time0
      read (restart_unit) vBkg

      write (*, *) 'Lb = ', Lb
      write (*, *) 'Nt0 = ', Nt0
      write (*, *) 'time0 = ', time0
      write (*, *) 'vBkg = ', vBkg
    end if
    call MPI_Bcast(Lb, 3, MPI_WP, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(iLb, 3, MPI_WP, 0, MPI_Comm_World, ierr)

    call MPI_Bcast(Nt0, 1, MPI_INTEGER, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(time0, 1, MPI_WP, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(vBkg, 3, MPI_WP, 0, MPI_Comm_World, ierr)

    ! cells
    if (rootWorld) then
      read (restart_unit) nrbc
      write (*, *) 'nrbc = ', nrbc
    end if
    call MPI_Bcast(nrbc, 1, MPI_Integer, 0, MPI_Comm_World, ierr)

    allocate (rbcs(nrbc))

    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      if (rootWorld) then
        read (restart_unit) nlat0, nlon0
        read (restart_unit) nlat, nlon
        read (restart_unit) celltype
        write (*, *) 'irbc : ', irbc, ' nlat0 = ', nlat0, 'type = ', celltype
        write (*, *) 'nlat =', nlat
      end if
      call MPI_Bcast(nlat0, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
      call MPI_Bcast(nlon0, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
      call MPI_Bcast(nlat, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
      call MPI_Bcast(nlon, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
      call MPI_Bcast(celltype, 1, MPI_Integer, 0, MPI_Comm_World, ierr)

      rbc%celltype = celltype
      dealias_fac = nlat/nlat0

      call Rbc_Create(rbc, nlat0, dealias_fac)

      if (rootWorld) then
        read (restart_unit) rbc%x
      end if
      call MPI_Bcast(rbc%x, size(rbc%x), MPI_WP, 0, MPI_Comm_World, ierr)
    end do

    if (rootWorld) then
      close (restart_unit)
    end if

  end subroutine

  subroutine ReadRestart_NoWalls(fn)
    character(*) :: fn

    integer :: irbc, iwall
    integer :: nlat0, nlon0, nlat, nlon, nvert, nele, celltype
    integer :: dealias_fac
    type(t_Rbc), pointer :: rbc
    type(t_Wall), pointer :: wall
    integer :: ierr
    character(*), parameter :: func_name = 'ReadRestart'

    ! global parameters
    if (rootWorld) then
      write (*, *) 'Read restart file ', TRIM(fn)
      print *, 'error 0'
      open (restart_unit, file=trim(fn), form='unformatted', action='read')

      print *, 'error 1'
      read (restart_unit) Lb
      iLb = 1./Lb
      print *, 'error 2'
      read (restart_unit) Nt0
      read (restart_unit) time0
      read (restart_unit) vBkg

      write (*, *) 'Lb = ', Lb
      write (*, *) 'Nt0 = ', Nt0
      write (*, *) 'time0 = ', time0
      write (*, *) 'vBkg = ', vBkg
    end if

    call MPI_Bcast(Lb, 3, MPI_WP, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(iLb, 3, MPI_WP, 0, MPI_Comm_World, ierr)

    call MPI_Bcast(Nt0, 1, MPI_INTEGER, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(time0, 1, MPI_WP, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(vBkg, 3, MPI_WP, 0, MPI_Comm_World, ierr)

    ! cells
    if (rootWorld) then
      read (restart_unit) nrbc

      write (*, *) 'nrbc = ', nrbc
    end if
    call MPI_Bcast(nrbc, 1, MPI_Integer, 0, MPI_Comm_World, ierr)

    allocate (rbcs(nrbc))

    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      if (rootWorld) then
        read (restart_unit) nlat0, nlon0
        read (restart_unit) nlat, nlon
        ! celltype = 1; print *,"NO READ CELL TYPE"
        read (restart_unit) celltype
        write (*, *) 'irbc : ', irbc, ' nlat0 = ', nlat0, 'type = ', celltype
        write (*, *) 'nlat =', nlat
      end if
      call MPI_Bcast(nlat0, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
      call MPI_Bcast(nlon0, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
      call MPI_Bcast(nlat, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
      call MPI_Bcast(nlon, 1, MPI_Integer, 0, MPI_Comm_World, ierr)
      call MPI_Bcast(celltype, 1, MPI_Integer, 0, MPI_Comm_World, ierr)

      rbc%celltype = celltype

      dealias_fac = nlat/nlat0
      call Rbc_Create(rbc, nlat0, dealias_fac)

      ! Check array dimension
      if (rootWorld) then
        if (rbc%nlat .ne. nlat .or. rbc%nlon .ne. nlon) then
          write (*, *) 'Subroutine ', func_name
          write (*, *) 'Error: invalid array dimension'
          stop
        end if
      end if

      if (rootWorld) then
        read (restart_unit) rbc%x
      end if
      call MPI_Bcast(rbc%x, size(rbc%x), MPI_WP, 0, MPI_Comm_World, ierr)
    end do ! irbc

    ! Walls
    if (rootWorld) then
      read (restart_unit) nwall

      write (*, *) 'nwall = ', nwall
    end if
    print *, 'sending walls'
    call MPI_Bcast(nwall, 1, MPI_Integer, 0, MPI_Comm_World, ierr)

    if (rootWorld) then
      close (restart_unit)
      write (*, *)
    end if

  end subroutine ReadRestart_NoWalls
!**********************************************************************

end module ModIO
