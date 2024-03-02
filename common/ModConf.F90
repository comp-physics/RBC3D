! Configuration of domain decomposition
module ModConf

  use MPI
  use ModDataTypes
  use ModDataStruct

  implicit none

  ! Domain size
  real(WP) :: Lb(3), iLb(3)

  ! PhysEwald, FourierEwald -- whether the node do physical and/or Fourier part
  !             of the Ewald sum
  logical :: PhysEwald, FourierEwald
  ! MPI_COMM_Ewald -- MPI communicator for computing Ewald sum
  integer :: MPI_COMM_Ewald
  ! Whether the node is the root of MPI_Comm_World and MPI_Comm_Ewald
  logical :: rootWorld, rootEwald
  ! only one node
  logical :: SINGLE_NODE

  ! Domain decomposition
  real(WP) :: nodeZmin, nodeZmax
  real(WP) :: nodeZminBuf, nodeZmaxBuf

  ! Ewald sum parameters
  ! alpha --
  ! eps -- error tolerance
  ! rc -- cut-off distance
  ! Nb -- number of mesh points for PME
  ! NbC -- number of Fourier coefficients for PME
  ! PBspln -- order of B-spline interpolation
  real(WP) :: alpha_Ewd, eps_Ewd
  real(WP) :: rc_Ewd
  integer :: Nb_Ewd(3), NbC_Ewd(3)
  integer :: PBspln_Ewd

  ! nCellTypes --- number of cell types
  integer :: nCellTypes

  ! viscRat -- viscosity ratio
  real(WP), allocatable, dimension(:)  :: viscRat
  real(WP), allocatable, dimension(:)  :: Acoef, Bcoef
  logical  :: Deflate

  real(WP), allocatable, dimension(:)  :: refRad  ! reference radius

  ! ForceCoef -- coefficient for intracell force
  ! epsDis -- intercell repulsion factor (not for intracell, yet!)
  ! viscRatThresh -- less than this to apply intercell repulsion
  ! rigidsep -- true for whole cell motion if viscRat > viscRatThres
  real(WP) :: ForceCoef
  real(WP) :: epsDist
  real(WP) :: viscRatThresh
  logical  :: rigidsep

  ! pGradTar -- target pressure gradient
  ! vBkg -- background velocity
  real(WP) :: pGradTar(3), vBkg(3)

  ! fmag -- magnetic force in xyz and toward tube center
  real(WP), dimension(4) :: fmags

  ! Time integration (from Nt0 to Nt)
  ! time0 -- initial time
  ! time -- current time
  ! Ts -- time step length
  integer :: Nt0, Nt
  real(WP) :: time0, time, Ts

  ! Input/output
  ! File name format: dir + prefix + number + suffix
  character(*), parameter :: fn_FMT = '(A,A,I9.9,A)'

  integer, parameter :: cell_unit = 21
  integer, parameter :: wall_unit = 22
  integer, parameter :: pgrad_unit = 23
  integer, parameter :: flow_unit = 24
  integer, parameter :: restart_unit = 25
  integer, parameter :: ftot_unit = 26

  integer :: cell_out
  integer :: wall_out
  integer :: pgrad_out
  integer :: flow_out
  integer :: restart_out
  integer :: ftot_out

  character(CHRLEN) :: restart_file

  public

  public :: InitMPI, &
            FinalizeMPI, &
            ReadConfig, &
            SetEwaldPrms, &
            DomainDecomp, &
            Is_Source, &
            Cell_Has_Source, &
            Tri_Has_Source, &
            CollectArray

contains

!**********************************************************************
! Initialize MPI and PETSC
  subroutine InitMPI(split_comm)
    logical, optional :: split_comm

#include "../petsc_include.h"
    integer :: numNodes, nodeNum
    character(MPI_Max_Processor_Name) :: machinename
    integer :: lenname
    integer :: stat(MPI_Status_Size)
    logical :: do_split_comm
    integer :: color, i
    integer :: ierr

    ! Init PETSc and MPI
    ! Note: MPI_Init is called by PETSc
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

    ! Get MPI information
    call MPI_Comm_Size(MPI_Comm_World, numNodes, ierr)
    call MPI_Comm_Rank(MPI_Comm_World, nodeNum, ierr)
    call MPI_Get_Processor_Name(machinename, lenname, ierr)

    ! Root of MPI_Comm_World indicator
    rootWorld = (nodeNum == 0)

    if (rootWorld) then
      print *, 'Node ', nodeNum, ' of ', numNodes, ' running on ', trim(machinename)
      do i = 1, numNodes - 1
        call MPI_Recv(lenname, 1, MPI_Integer, i, 1, MPI_Comm_World, stat, ierr)
        call MPI_Recv(machinename, lenname, MPI_Character, i, 1, MPI_Comm_World, stat, ierr)
        print *, 'Node ', i, ' of ', numNodes, ' running on ', trim(machinename)
      end do ! i
    else
      call MPI_Send(lenname, 1, MPI_Integer, 0, 1, MPI_Comm_World, stat, ierr)
      call MPI_Send(machinename, lenname, MPI_Character, 0, 1, MPI_Comm_World, stat, ierr)
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! Split communicator
    if (.not. present(split_comm)) then
      do_split_comm = .false.
    else if (.not. split_comm) then
      do_split_comm = .false.
    else if (numNodes > 16) then
      do_split_comm = .true.
    end if

    if (do_split_comm) then
      ! Generate two communicators for physical and Fourier Ewald sum
      if (nodeNum < numNodes - 16) then
        PhysEwald = .true.
        FourierEwald = .false.
        color = 1
      else
        PhysEwald = .false.
        FourierEwald = .true.
        color = 2
      end if

      call MPI_Comm_Split(MPI_COMM_WORLD, color, 0, MPI_COMM_Ewald, ierr)
    else
      ! Use the same communicator for physical and Fourier Ewald sum
      PhysEwald = .true.
      FourierEwald = .true.

      MPI_COMM_Ewald = MPI_COMM_WORLD
    end if

    ! root indicator for MPI_Comm_Ewald
    call MPI_Comm_Rank(MPI_Comm_Ewald, nodeNum, ierr)
    rootEwald = (nodeNum == 0)

    ! SINGLE_NODE indicator
    call MPI_Comm_Size(MPI_Comm_Ewald, numNodes, ierr)
    SINGLE_NODE = (numNodes == 1)

    if (PhysEwald) then
      call MPI_Comm_Size(MPI_COMM_Ewald, numNodes, ierr)
      call MPI_Comm_Rank(MPI_COMM_Ewald, nodeNum, ierr)
      if (rootEwald) then
        print *, 'Number of nodes for physical Ewald sum = ', numNodes
      end if
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    if (FourierEwald) then
      call MPI_Comm_Size(MPI_COMM_Ewald, numNodes, ierr)
      call MPI_Comm_Rank(MPI_COMM_Ewald, nodeNum, ierr)
      if (rootEwald) then
        print *, 'Number of nodes for Fourier Ewald sum = ', numNodes
      end if
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    if (rootWorld) then
      print *
    end if

  end subroutine InitMPI

!**********************************************************************
! Finalize MPI and PETSc
  subroutine FinalizeMPI

#include "../petsc_include.h"
    integer :: ierr

    if (MPI_COMM_Ewald /= MPI_COMM_WORLD) then
      call MPI_Comm_Free(MPI_COMM_Ewald, ierr)
    end if

    ! Note: MPI_Finalize is called by PETSc
    call PetscFinalize(ierr)

  end subroutine FinalizeMPI

!**********************************************************************
! Configure global parameters
! Arguments:
!  fn -- name of the configure file
  subroutine ReadConfig(fn)
    character(*) :: fn

    integer, parameter :: MAXCELLTYPES = 128
    real(WP), dimension(MAXCELLTYPES)  :: viscRatIN
    real(WP), dimension(MAXCELLTYPES)  :: refRadIN

    integer, parameter :: funit = 1
    integer :: ierr

    ! Root node reads configuration file
    if (rootWorld) then
      open (funit, file=trim(fn), action='read')

      read (funit, *) alpha_Ewd; print *, 'alpha_Ewd = ', alpha_Ewd
      read (funit, *) eps_Ewd; print *, 'eps_Ewd = ', eps_Ewd
      read (funit, *) PBspln_Ewd; print *, 'PBspln_Ewd = ', PBspln_Ewd

      read (funit, *) nCellTypes; print *, 'nCellTypes = ', nCelltypes
      read (funit, *) viscRatIN(1:nCellTypes)
      print *, 'viscRat = ', viscRatIN(1:nCellTypes)
      read (funit, *) refRadIN(1:nCellTypes)
      print *, 'refRad = ', refRadIN(1:nCellTypes)

      read (funit, *) Deflate; print *, 'Deflate = ', Deflate
      !print *, 'before error'
      read (funit, *) pGradTar(1); print *, 'pGradTar = ', pGradTar(1)
      read (funit, *) pGradTar(2); print *, 'pGradTar = ', pGradTar(2)
      read (funit, *) pGradTar(3); print *, 'pGradTar = ', pGradTar(3)
      !print *, 'after error'
      read (funit, *) Nt; print *, 'Nt = ', Nt
      read (funit, *) Ts; print *, 'Ts = ', Ts

      read (funit, *) cell_out; print *, 'cell_out = ', cell_out
      read (funit, *) wall_out; print *, 'wall_out = ', wall_out
      read (funit, *) pgrad_out; print *, 'pgrad_out = ', pgrad_out
      read (funit, *) flow_out; print *, 'flow_out = ', flow_out
      read (funit, *) ftot_out; print *, 'flow_out = ', ftot_out
      read (funit, *) restart_out; print *, 'restart_out = ', restart_out

      read (funit, *) restart_file; print *, 'restart file = ', trim(restart_file)

      read (funit, *) epsDist; print *, 'epsDist = ', epsDist
      read (funit, *) ForceCoef; print *, 'ForceCoef = ', ForceCoef
      read (funit, *) viscRatThresh; print *, 'viscRatThresh = ', viscRatThresh
      read (funit, *) rigidsep; print *, 'rigidsep = ', rigidsep
      read (funit, *) fmags; print *, 'fmags = ', fmags
      close (funit)

      print *
    end if

    ! Broadcast values
    call MPI_Bcast(alpha_Ewd, 1, MPI_WP, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(eps_Ewd, 1, MPI_WP, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(PBspln_Ewd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    call MPI_Bcast(nCellTypes, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    allocate (viscRat(nCellTypes))
    viscRat = viscRatIN(1:nCellTypes)
    allocate (refRad(nCellTypes))
    refRad = refRadIN(1:nCellTypes)

    call MPI_Bcast(viscRat, nCellTypes, MPI_WP, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(refRad, nCellTypes, MPI_WP, 0, MPI_COMM_WORLD, ierr)
    call initCOEFs

    call MPI_Bcast(Deflate, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    call MPI_Bcast(pgradTar, 3, MPI_WP, 0, MPI_Comm_World, ierr)

    call MPI_Bcast(Nt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(Ts, 1, MPI_WP, 0, MPI_COMM_WORLD, ierr)

    call MPI_Bcast(cell_out, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(wall_out, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(pgrad_out, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(flow_out, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(ftot_out, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(restart_out, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    call MPI_Bcast(restart_file, CHRLEN, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

    call MPI_Bcast(epsDist, 1, MPI_WP, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(ForceCoef, 1, MPI_WP, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(viscRatThresh, 1, MPI_WP, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(rigidsep, 1, MPI_LOGICAL, 0, MPI_Comm_World, ierr)
    call MPI_Bcast(fmags, 4, MPI_WP, 0, MPI_Comm_World, ierr)

  end subroutine ReadConfig

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This is initializing coefficients in the velcoity solver
!
!  The boundary integral equations can be put in the form
!
!   A q - B K q/4 Pi = R
!
!  where K is the double layer operator.  Acoef and Bcoef are A and B.
  subroutine initCOEFs

    integer  :: l

    allocate (Acoef(nCellTypes), Bcoef(nCellTypes))

    do l = 1, nCellTypes
      if (viscRat(l) .gt. 0) then  ! Finite viscosity cell
        Acoef(l) = 1.+viscRat(l)
        Bcoef(l) = 1.-viscRat(l)
      else  ! Rigid cell
        Acoef(l) = 1.   ! +1.
        Bcoef(l) = -1.    ! -1.
      end if
      print *, l, Acoef(l), Bcoef(l), viscRat(l)
    end do

  end subroutine initCOEFs

!**********************************************************************
! Compute Ewald sum parameters
  subroutine SetEwaldPrms

    integer :: numNodes, ierr
    integer :: iter
    integer, parameter :: iterMax = 10
    real(WP) :: s

    ! Cut-off distance for
    if (PhysEwald) then
      if (rootEwald) then
        ! Physical sum parameters
        ! Note: cut-off radius is such that
        !  (s^3  + 1.5*s^2 + 0.75/s)*exp(-s^2) < 0.75*sqrt(PI)*eps
        !  at the cut-off radius,
        !  where s = r*sqrt(PI/alpha) and eps is the error tolerance.
        !
        !  Here, the error tolerance is defined to be
        !  the double-layer Ewald physical sum divided by the Green's
        !  function of the infinite space at the cut-off radius.

        ! Use iteration to determine the cut-off radius
        s = 1.
        do iter = 1, iterMax
          s = 0.75*sqrt(PI)*eps_Ewd/(s**3 + 1.5*s + 0.75/s)
          s = sqrt(-log(s))
        end do ! iter
        rc_Ewd = sqrt(alpha_Ewd/PI)*s

        ! Need at least three cells in each direction
        rc_Ewd = min(minval(Lb/3.001), rc_Ewd)

        write (*, '(A)') 'Derived Ewald sum parameters:'
        write (*, '(A,F15.5)') 'rc = ', rc_Ewd
        write (*, '(A,3I5)') 'Nc = ', floor(Lb/rc_Ewd)
      end if

      call MPI_Bcast(rc_Ewd, 1, MPI_WP, 0, MPI_Comm_Ewald, ierr)
    end if
    call MPI_Barrier(MPI_Comm_World, ierr)

    if (FourierEwald) then
      if (rootEwald) then
        call MPI_Comm_Size(MPI_Comm_Ewald, numNodes, ierr)

        ! Fourier sum parameters
        Nb_Ewd = 2*ceiling(sqrt(-log(eps_Ewd)/(pi*alpha_Ewd))*Lb)
        Nb_Ewd(3) = max(Nb_Ewd(3), numNodes*PBspln_Ewd)
        Nb_Ewd(3) = ceiling(real(Nb_Ewd(3))/numNodes)*numNodes
        NbC_Ewd = (/Nb_Ewd(1)/2 + 1, Nb_Ewd(2), Nb_Ewd(3)/)

        write (*, '(A,3I5)') 'Nb = ', Nb_Ewd
        write (*, '(A,3I5)') 'NbC = ', NbC_Ewd
        write (*, *)
      end if

      call MPI_Bcast(Nb_Ewd, 3, MPI_Integer, 0, MPI_Comm_Ewald, ierr)
      call MPI_Bcast(NbC_Ewd, 3, MPI_Integer, 0, MPI_Comm_Ewald, ierr)
    end if
    call MPI_Barrier(MPI_Comm_World, ierr)

  end subroutine SetEwaldPrms

!**********************************************************************
! Initialize domain decompsotion
  subroutine DomainDecomp

    integer :: numNodes, nodeNum
    integer :: ierr
    real(WP) :: hz, hbuf

    ! Determine the domain boundary in z-direction
    call MPI_Comm_Size(MPI_COMM_Ewald, numNodes, ierr)
    call MPI_Comm_Rank(MPI_Comm_Ewald, nodeNum, ierr)
    hz = Lb(3)/numNodes

    hbuf = 0.
    if (PhysEwald) then
      hbuf = max(hbuf, rc_Ewd)
    end if
    if (FourierEwald) then
      hbuf = max(hbuf, PBspln_Ewd*Lb(3)/Nb_Ewd(3))
    end if

    nodeZmin = nodeNum*hz
    nodeZmax = (nodeNum + 1)*hz

    nodeZminBuf = nodeZmin - hbuf
    nodeZmaxBuf = nodeZmax + hbuf

  end subroutine DomainDecomp

!**********************************************************************
! Wheather a source point make contribution to the local domain
  function Is_Source(x)
    real(WP) :: x(3)
    logical :: Is_Source

    integer :: numNodes
    real(WP) :: z
    real(WP), parameter :: EPS = 1.D-5
    integer :: ierr

    if (SINGLE_NODE) then
      Is_Source = .true.
      return
    end if

    z = x(3) - nodeZminBuf
    z = z - floor(z*iLb(3))*Lb(3)

    if (z < nodeZmaxBuf - nodeZminBuf + EPS) then
      Is_Source = .true.
    else
      Is_Source = .false.
    end if

  end function Is_Source

!**********************************************************************
  function Cell_Has_Source(cell) result(hasSource)
    type(t_Rbc) :: cell
    logical :: hasSource

    real(WP) :: zmin, zmax
    real(WP), parameter :: EPS = 1.D-5

    if (SINGLE_NODE) then
      hasSource = .true.
      return
    end if

    zmin = minval(cell%x(:, :, 3)) - nodeZminBuf
    zmin = zmin - floor(zmin*iLb(3))*Lb(3)

    zmax = maxval(cell%x(:, :, 3)) - nodeZminBuf
    zmax = zmax - floor(zmax*iLb(3))*Lb(3)

    if (zmin < nodeZmaxBuf - nodeZminBuf + EPS) then
      hasSource = .true.
    else if (zmax < zmin) then
      hasSource = .true.
    else
      hasSource = .false.
    end if

  end function Cell_Has_Source

!**********************************************************************
! Wheather a triangle has active source points
! Arguments:
!  x(i,:) -- the coordinates of the i-th vertex
  function Tri_Has_Source(x) result(hasSource)
    real(WP) :: x(3, 3)
    logical :: hasSource

    logical :: T1, T2, T3

    if (SINGLE_NODE) then
      hasSource = .true.
      return
    end if

    T1 = Is_Source(x(1, :))
    T2 = Is_Source(x(2, :))
    T3 = Is_Source(x(3, :))
    hasSource = T1 .and. T2 .and. T3

  end function Tri_Has_Source

!**********************************************************************
! Collect array
! Arguments:
!  np -- number of points to send
!  ia -- index of elements of a(:) in the global array
!  a -- array to send
!  aGlb -- global array
!  COMM -- the communicator
!  OP -- MPI_OP
! Note:
!  -- For now, if OP is not present, aGlb(:) is set to the values of a(:)
!     i.e. aGlb(ia(i)) = a(i)
!  -- If OP does exists, the elements in a(:) are added to those of aGlb(:)
!     i.e. aGlb(ia(i)) += a(i)
  subroutine CollectArray(np, nvar, ia, a, aGlb, COMM, OP)
    integer :: np, nvar
    integer :: ia(:)
    real(WP) :: a(:, :)
    real(WP) :: aGlb(:, :)
    integer :: COMM
    integer, optional :: OP

    integer :: numNodes, nodeNum
    integer :: npTot, sendcount
    integer, allocatable :: nps(:), recvcounts(:), displs(:), iaRecv(:)
    real(WP), allocatable :: aSend(:), aRecv(:)
    integer :: i, p, ierr

    call MPI_Comm_Size(comm, numNodes, ierr)
    call MPI_Comm_Rank(comm, nodeNum, ierr)

    ! Allocate working arrays
    allocate (nps(0:numNodes - 1), recvcounts(0:numNodes - 1), displs(0:numNodes - 1))

    call MPI_AllGather(np, 1, MPI_Integer, nps, 1, MPI_Integer, comm, ierr)
    npTot = sum(nps)
    allocate (iaRecv(npTot))
    allocate (aSend(np*nvar), aRecv(npTot*nvar))

    ! Send and receive ia(:)
    sendCount = np

    recvcounts = nps
    displs(0) = 0
    do i = 1, numNodes - 1
      displs(i) = displs(i - 1) + recvCounts(i - 1)
    end do ! i

    call MPI_AllGatherV(ia, sendcount, MPI_INTEGER, &
                        iaRecv, recvcounts, displs, MPI_INTEGER, COMM, ierr)

    ! Send and receive a(:,:)
    sendCount = np*nvar
    do i = 1, np
      aSend((i - 1)*nvar + 1:i*nvar) = a(i, :)
    end do ! i

    recvcounts = nps*nvar
    displs(0) = 0
    do i = 1, numNodes - 1
      displs(i) = displs(i - 1) + recvcounts(i - 1)
    end do ! i
    call MPI_AllGatherV(aSend, sendcount, MPI_WP, &
                        aRecv, recvcounts, displs, MPI_WP, COMM, ierr)

    ! Re-assemble data
    do i = 1, npTot
      p = iaRecv(i)
      if (.not. present(OP)) then
        aGlb(p, :) = aRecv((i - 1)*nvar + 1:i*nvar)
      else
        aGlb(p, :) = aGlb(p, :) + aRecv((i - 1)*nvar + 1:i*nvar)
      end if
    end do ! i

    ! Deallocate working arrays
    deallocate (nps, recvcounts, displs)
    deallocate (iaRecv)
    deallocate (aSend, aRecv)

  end subroutine CollectArray

!**********************************************************************

end module ModConf
