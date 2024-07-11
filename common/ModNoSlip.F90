! Solve the wall force density
module ModNoSlip

  use ModDataTypes
  use ModDataStruct
  use ModData
  use ModSpline
  use ModRbc
  use ModIntOnRbcs
  use ModIntOnWalls
  use ModPME
  use ModSourceList
  use ModTargetList
  use ModSphPk

#include "petsc/finclude/petsc.h"
  use petsc

  implicit none

  ! npoint -- total number of wall mesh points, with redundancy
  ! dof -- degrees of freedom
  integer :: npoint, dof

  Mat :: mat_lhs
  Vec :: vec_rhs, vec_sol
  KSP :: ksp_lhs
  PC :: pc_lhs

  private

  public :: NoSlipWall, &
            Compute_Wall_Residual_Vel, &
            WallBuildMat   !! TESTING

  private :: MyMatMult, &
             AssembleArray, &
             My_VecGetValues, &
             My_VecSetValues

contains

!**********************************************************************
! Update the wall force density to enforce no-slip condition on wall
  subroutine NoSlipWall

    logical, save :: solver_inited = .false.
    integer :: iwall, ivert, p
    type(t_wall), pointer :: wall
    real(WP), allocatable :: v(:, :), v1D(:)
    real(WP), allocatable :: f0(:, :), finc(:, :), finc1D(:)
    integer :: niter
    real(WP) :: residual, vmax, vzmax
    integer :: ierr

    if (nwall == 0) return

    ! Initialize the solver
    if (.not. solver_inited) then
      ! Calculate system size
      npoint = 0
      dof = 0
      do iwall = 1, nwall
        wall => walls(iwall)
        npoint = npoint + wall%nvert
        dof = dof + 3*count(wall%v2v <= 0)
      end do ! iwall

      ! Set up petsc matrix-free GMRES solver
      call VecCreateSeq(PETSC_COMM_SELF, dof, vec_rhs, ierr)
      call VecCreateSeq(PETSC_COMM_SELF, dof, vec_sol, ierr)

      call MatCreateShell(PETSC_COMM_SELF, dof, dof, dof, dof, PETSC_NULL_INTEGER, mat_lhs, ierr)
      call MatShellSetOperation(mat_lhs, MATOP_MULT, MyMatMult, ierr)

      call KSPCreate(PETSC_COMM_SELF, ksp_lhs, ierr)
      ! call KSPSetOperators(ksp_lhs, mat_lhs, mat_lhs, SAME_NONZERO_PATTERN, ierr)
      call KSPSetOperators(ksp_lhs, mat_lhs, mat_lhs, ierr)
      call KSPSetType(ksp_lhs, KSPGMRES, ierr)
      call KSPGetPC(ksp_lhs, pc_lhs, ierr)
!      call KSPSetInitialGuessNonzero(ksp_lhs, PETSC_TRUE, ierr)  ! TRIAL
      call PCSetType(pc_lhs, PCNONE, ierr)

      !SHB modify from 1.D-3 to 1.D-6 and now to eps_Ewd
      call KSPSetTolerances(ksp_lhs, eps_Ewd, &
                            PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, &
                            80, ierr)
!!$      print *
!!$      print *, "SMALL WALL ITERATIONS"
!!$      print *, "SMALL WALL ITERATIONS"
!!$      print *, "SMALL WALL ITERATIONS"
!!$      print *, "SMALL WALL ITERATIONS"
!!$      print *

      solver_inited = .true.
    end if

    ! Allocate working arrays
    allocate (v(npoint, 3), v1D(dof))
    allocate (f0(npoint, 3), finc(npoint, 3), finc1D(dof))

    ! Save current wall force density
    p = 0
    do iwall = 1, nwall
      wall => walls(iwall)
      f0(p + 1:p + wall%nvert, :) = wall%f
      p = p + wall%nvert
    end do ! iwall

    ! Compute rhs
    call Compute_Wall_Residual_Vel(v)
    v = -v
    call AssembleArray(v, v1D, 1)
    call My_VecSetValues(vec_rhs, v1D)

    ! Solve the equation
    call VecZeroEntries(vec_sol, ierr)
    call KSPSolve(ksp_lhs, vec_rhs, vec_sol, ierr)

    call KSPGetIterationNumber(ksp_lhs, niter, ierr)
    call KSPGetResidualNorm(ksp_lhs, residual, ierr)

    if (rootWorld) then
      write (*, '(A)') 'Wall iteration:'
      write (*, '(A,I5,A,ES12.2)') 'niter = ', niter, ' residual = ', residual
    end if

    call My_VecGetValues(vec_sol, finc1D)
    call AssembleArray(finc, finc1D, -1)

    ! Update the wall force density
    p = 0
    do iwall = 1, nwall
      wall => walls(iwall)
      wall%f = f0(p + 1:p + wall%nvert, :) + finc(p + 1:p + wall%nvert, :)
      p = p + wall%nvert
    end do ! iwall

    ! Monitor residual velocity
    call Compute_Wall_Residual_Vel(v)
    if (rootWorld) then
      vmax = maxval(abs(v))
      vzmax = maxval(abs(v(:, 3)))
      write (*, '(A,ES12.2,A,ES12.2)') 'Residual velocity: vmax = ', vmax, ' vzmax = ', vzmax
    end if

    ! Deallocate working arrays
    deallocate (v, f0, finc, finc1D)

  end subroutine NoSlipWall

!**********************************************************************
  subroutine Compute_Wall_Residual_Vel(v)
    real(WP) :: v(:, :)

    type(t_targetlist), pointer :: tlist
    real(WP) :: c1, c2
    integer :: i, ii
    real(WP) :: vmax, vzmax
    integer :: ierr
    integer :: iwall, iele, ivert, l, offset
    type(t_Wall), pointer :: wall
    real(WP) :: vele(3, 3), vc(3)

    ! Allocate working arrays
    tlist => tlist_wall

    ! Compute the velocity contribution from cell surfaces to the wall
    ! and from the background velocity
    v = 0.

    c1 = 1./(4.*PI)  ! now 4 Pi because (tlist_wall%lam + 1) = 2 included later
    c2 = 1./(4.*Pi) !(1 - viscRat)/(8*PI)

    if (PhysEwald) then
      call AddIntOnRbcs(c1, c2, tlist, v)
      call AddIntOnWalls(c1, tlist, v)
    end if

    if (FourierEwald) then
      call PME_Distrib_Source(c1, c2, rbcs, walls)
      call PME_Transform
      call PME_Add_Interp_Vel(tlist, v)
    end if

    call TargetList_CollectArray(tlist, 3, v, MPI_COMM_WORLD)

    ! Add background velocity
    do ii = 1, 3
      v(:, ii) = v(:, ii) + vBkg(ii)   ! 2/(lam + 1) target = 1
    end do ! ii

  end subroutine Compute_Wall_Residual_Vel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Testing Routine for assembling the A matrix
!
! u = alpha A x + b
!
! It is the eigensystem of A that is to be deflated

  SUBROUTINE WallBuildMat(dof, AA)
    integer                   :: dof
    real(WP), dimension(dof, dof)  :: AA

    real(WP), allocatable :: g(:, :), v(:, :)

    real(WP), dimension(dof)      :: u, b

    Mat :: mat_t
    Vec :: vec_u, vec_b

    integer :: nPoint
    type(t_SourceList), pointer :: slist
    type(t_TargetList), pointer :: tlist

    type(t_wall), pointer :: wall
    integer :: ierr
    integer :: i, j, l, m, p, q, iwall

    ! do iwall = 1 , nwall
    !    wall => walls(iwall)
    !    call Wall_ComputeGeometry(wall)
    ! end do ! iwall

    call VecCreateSeq(PETSC_COMM_SELF, dof, vec_u, ierr)
    call VecCreateSeq(PETSC_COMM_SELF, dof, vec_b, ierr)

    call MatCreateShell(PETSC_COMM_SELF, dof, dof, dof, &
                        dof, PETSC_NULL_INTEGER, mat_t, ierr)
    call MatShellSetOperation(mat_t, MATOP_MULT, MyMatMult, ierr)

    do i = 1, dof
      u = 0.
      u(i) = 1.

      if (rootWorld) print *, i, " of ", dof

      call My_VecSetValues(vec_u, u)
      call MatMult(mat_t, vec_u, vec_b, ierr)
      call My_VecGetValues(vec_b, b)

      AA(:, i) = b

    end do

  end SUBROUTINE WallBuildMat

!**********************************************************************
! Matmul for matrix-free GMRES solver
! Arguments:
!  lhs_petsc -- lhs matrix
!  vec_u, vec_b -- b = lhs*u
!  ierr -- error code
  subroutine MyMatMult(lhs_petsc, vec_u, vec_b, ierr)
    Mat :: lhs_petsc
    Vec :: vec_u, vec_b
    integer :: ierr

    real(WP), allocatable :: f1D(:), f(:, :), v(:, :), v1D(:)
    type(t_TargetList), pointer :: tlist
    type(t_Wall), pointer :: wall
    integer :: iwall, p
    real(WP) :: c1, c2

    ! Set up parameters
    tlist => tlist_wall

    ! Allocate working arrays
    allocate (f1D(dof), f(npoint, 3), v(npoint, 3), v1D(dof))

    ! Transform vec_u to single-layer density on walls
    call My_VecGetValues(vec_u, f1D)
    call AssembleArray(f, f1D, -1)

    p = 0
    do iwall = 1, nwall
      wall => walls(iwall)
      wall%f = f(p + 1:p + wall%nvert, :)
      p = p + wall%nvert
    end do ! iwall

    ! Calculate the induced velocity field by the single-layer potential
    c1 = 1./(4.*PI)  ! switched to 4 Pi
    c2 = 0.

    v = 0.

    if (PhysEwald) then
      call AddIntOnWalls(c1, tlist, v)
    end if

    if (FourierEwald) then
      call PME_Distrib_Source(c1, c2, walls=walls)
      call PME_Transform
      call PME_Add_Interp_Vel(tlist, v)
    end if

    call TargetList_CollectArray(tlist, 3, v, MPI_COMM_WORLD)

    ! Convert to 1D array
    call AssembleArray(v, v1D, 1)
    call My_VecSetValues(vec_b, v1D)

    ! Deallocate working arrays
    deallocate (f1D, f, v, v1D)

  end subroutine MyMatMult

!**********************************************************************
  subroutine My_VecGetValues(x, a)
    Vec x
    real(WP) :: a(:)

    integer, allocatable :: ix(:)
    integer :: n, i, ierr

    ! Allocate working arrays
    call VecGetSize(x, n, ierr)
    allocate (ix(n))
    ix = (/(i, i=0, n - 1)/)

    call VecGetValues(x, n, ix, a, ierr)

    ! Deallocate working arrays
    deallocate (ix)

  end subroutine My_VecGetValues

!**********************************************************************
  subroutine My_VecSetValues(x, a)
    Vec x
    real(WP) :: a(:)

    integer, allocatable :: ix(:)
    integer :: n, i, ierr

    ! Allocate working arrays
    call VecGetSize(x, n, ierr)
    allocate (ix(n))
    ix = (/(i, i=0, n - 1)/)

    call VecSetValues(x, n, ix, a, INSERT_VALUES, ierr)
    call VecAssemblyBegin(x, ierr)
    call VecAssemblyEnd(x, ierr)

    ! Deallocate working arrays
    deallocate (ix)

  end subroutine My_VecSetValues

!**********************************************************************
! Assemble to and from a 1-dimensional array that has no redundancy
! Arguments:
!  u(:,:) --
!  u1D(:)
!  direction --
! Note:
!  direction = 1: u -> u1D
!           = -1, u <= u1D
  subroutine AssembleArray(u, u1D, direction)
    real(WP) :: u(:, :), u1D(:)
    integer :: direction

    type(t_wall), pointer :: wall
    integer :: iwall, ivert, pu, pu1D

    pu = 0

    do iwall = 1, nwall
      wall => walls(iwall)
      do ivert = 1, wall%nvert
        pu = pu + 1
        pu1D = wall%indxVertGlb(ivert)

        if (direction == 1) then
          u1D(3*pu1D - 2:3*pu1D) = u(pu, :)
        else
          u(pu, :) = u1D(3*pu1D - 2:3*pu1D)
        end if
      end do !ivert
    end do ! iwall

  end subroutine AssembleArray

!**********************************************************************

end module ModNoSlip
