! Solve the cell surface velocities for mismatched viscosity
module ModVelSolver

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

  Mat :: mat_lhs
  Vec :: vec_rhs, vec_sol
  KSP :: ksp_lhs
  PC :: pc_lhs

  private

  public :: Solve_RBC_Vel, &
            MyMatMult, &
            Glob_Sph_Trans, &
            DeflateGetSBModes, & !TEMP
            BuildMat !TEMP

  private :: Compute_Rhs, &
             My_VecGetValues, &
             My_VecSetValues

contains

!**********************************************************************
! Solve RBC surface velocities
! Note:
!   v -- rbc surface veocities
  subroutine Solve_RBC_Vel(v)
    real(WP) :: v(:, :)

    logical, save :: solver_inited = .false.
    integer :: dof
    integer :: irbc
    type(t_RBC), pointer :: rbc
    integer :: nlat, nlon
    integer :: ii, p
    real(WP), allocatable, save :: rhs(:), sol(:)
    integer :: niter
    real(WP) :: residual
    integer :: ierr

    if (nrbc == 0) return

    ! Initialize the solver
    if (.not. solver_inited) then
      ! Calculate matrix size
      dof = 0
      do irbc = 1, nrbc
        rbc => rbcs(irbc)
        dof = dof + 3*(rbc%nlat0**2)
      end do ! irbc
      write (*, *) "Just calculated dof, it is ", dof

      ! Set up petsc matrix-free GMRES solver
      allocate (rhs(dof), sol(dof))
      write (*, *) "Set up GMRES solver"

      call VecCreateSeq(PETSC_COMM_SELF, dof, vec_rhs, ierr)
      call VecCreateSeq(PETSC_COMM_SELF, dof, vec_sol, ierr)

      call MatCreateShell(PETSC_COMM_SELF, dof, dof, dof, dof, PETSC_NULL_INTEGER, mat_lhs, ierr)
      call MatShellSetOperation(mat_lhs, MATOP_MULT, MyMatMult, ierr)

      call KSPCreate(PETSC_COMM_SELF, ksp_lhs, ierr)
      call KSPSetOperators(ksp_lhs, mat_lhs, mat_lhs, ierr) ! call KSPSetOperators(ksp_lhs, mat_lhs, mat_lhs, SAME_NONZERO_PATTERN, ierr)
      call KSPSetType(ksp_lhs, KSPGMRES, ierr)
      call KSPGetPC(ksp_lhs, pc_lhs, ierr)
      call KSPSetInitialGuessNonzero(ksp_lhs, PETSC_TRUE, ierr)  ! TRIAL
      call PCSetType(pc_lhs, PCNONE, ierr)

      call KSPSetTolerances(ksp_lhs, 1.D-11, &
                            PETSC_DEFAULT_REAL, &
                            PETSC_DEFAULT_REAL, &
                            200, ierr); print *, "SETTING HIGH MAX ITER/ERR"

      solver_inited = .true.
    end if
!write(*,*) "solver initiated"
    ! Compute rhs
    call Compute_Rhs(rhs)
!write(*,*) "rhs computed"
    if (Deflate) then
      !if (rootWorld) print *,"DEFLATING --- GET SB MODES"
      call DeflateGetSBModes
    end if

    ! Solve the equation
    if (MAXVAL(viscRat) == 1 .and. MINVAL(viscRat) == 1) then
      sol = rhs
    else

      call My_VecSetValues(vec_rhs, rhs)

!      call My_VecSetValues(vec_sol, sol); print *, "SETTING INITIAL GUESS"

      call KSPSolve(ksp_lhs, vec_rhs, vec_sol, ierr)
      call My_VecGetValues(vec_sol, sol)

      call KSPGetIterationNumber(ksp_lhs, niter, ierr)
      call KSPGetResidualNorm(ksp_lhs, residual, ierr)

      print *, 'iterations'
      if (rootWorld) then
        write (*, '(A,I5,A,ES12.2)') 'niter = ', niter, ' residual = ', residual
      end if
    end if

    ! Get cell surface velocity
    call Glob_Sph_Trans(v, sol, FOUR_TO_PHYS)

    if (Deflate) then
      if (rootWorld) print *, "DEFLATING -- SOLUTION RECOVERY"
      call DeflateRecoverSol(size(v, 1), v)
    end if

  end subroutine Solve_RBC_Vel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Deflation as laid out by Pozrikis in his book
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Recover the original solution from the solution of the deflated system by
! adding back in the necessary eigenvalues
!
! x' = alpha A x' + b   ===>   recover x from x'
!
  subroutine DeflateRecoverSol(npoint, v)
    integer                       :: npoint
    real(WP), dimension(:, :)      :: v

    real(WP), allocatable, dimension(:, :, :) :: vv, vvw

    integer  :: irbc
    type(t_RBC), pointer :: rbc
    integer :: nlat0, nlon0, nlat, nlon, ilat, ilon, nvert

    real(WP) :: alpha, afac

    real(WP) :: ds, fac
    integer  :: p, m

    nlat = rbcs(1)%nlat
    nlon = rbcs(1)%nlon
    allocate (vv(nlat, nlon, 3), vvw(nlat, nlon, 3))

    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      alpha = Bcoef(rbc%celltype)/Acoef(rbc%celltype) !COEF
!       alpha = (1.-viscRat(rbc%celltype))/(1.+viscRat(rbc%celltype)) !COEF

      nlat = rbc%nlat
      nlon = rbc%nlon
      p = 1 + (irbc - 1)*(nlat*nlon)

      do ilon = 1, nlon
        do ilat = 1, nlat
          vv(ilat, ilon, :) = v(p, :)
          p = p + 1
        end do
      end do

      if (viscRat(rbc%celltype) .gt. 0) then  ! FINITE VISCOSITY
        vvw = vv
        afac = alpha/(alpha + 1.)
      else  ! RIGID
        vvw = 0.
        rbc%g = vv     ! store potential for rigid-body source
        ! this is needed for no slip wall and post process
        afac = -1./2.   ! -1./2.
        print *, "NEGATIVE RIGID RECOVERY"
      end if

      do m = 1, 6
        fac = 0.
        do ilon = 1, nlon
          do ilat = 1, nlat
            ds = rbc%detj(ilat, ilon)*rbc%w(ilat)
            fac = fac + SUM(vv(ilat, ilon, :)*rbc%qq(ilat, ilon, :, m))*ds
          end do ! ilon
        end do ! ilat
        vvw = vvw - afac*rbc%qq(:, :, :, m)*fac
      end do

      p = 1 + (irbc - 1)*(nlat*nlon)
      do ilon = 1, nlon
        do ilat = 1, nlat
          v(p, :) = vvw(ilat, ilon, :)
          p = p + 1
        end do
      end do

    end do

    deallocate (vv, vvw)

  end subroutine DeflateRecoverSol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modify the "forcing term" in x = A x + b:   b -> b'
!
  subroutine DeflateRHS(npoint, v)
    integer                       :: npoint
    real(WP), dimension(npoint, 3) :: v

    real(WP), allocatable, dimension(:, :, :) :: FF

    integer  :: irbc
    type(t_RBC), pointer :: rbc
    integer :: nlat0, nlon0, nlat, nlon, ilat, ilon, nvert

    real(WP) :: alpha

    real(WP) :: ds, fac
    integer  :: p

    nlat = rbcs(1)%nlat
    nlon = rbcs(1)%nlon
    allocate (FF(nlat, nlon, 3))

    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      alpha = Bcoef(rbc%celltype)/Acoef(rbc%celltype) !COEF
!       alpha = (1.-viscRat(rbc%celltype))/(1.+viscRat(rbc%celltype)) !COEF

      nlat = rbc%nlat
      nlon = rbc%nlon
      p = 1 + (irbc - 1)*(nlat*nlon)

      do ilon = 1, nlon
        do ilat = 1, nlat
          FF(ilat, ilon, :) = v(p, :)
          p = p + 1
        end do
      end do

      fac = 0.
      do ilon = 1, nlon
        do ilat = 1, nlat
          ds = rbc%detj(ilat, ilon)*rbc%w(ilat)
          fac = fac + SUM(FF(ilat, ilon, :)*rbc%a3(ilat, ilon, :))*ds
        end do ! ilon
      end do ! ilat
      FF = FF + alpha/(1.-alpha)*rbc%a3*fac/rbc%area

      p = 1 + (irbc - 1)*(nlat*nlon)
      do ilon = 1, nlon
        do ilat = 1, nlat
          v(p, :) = FF(ilat, ilon, :)
          p = p + 1
        end do
      end do

    end do

  end subroutine DeflateRHS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Modify the mat-vec multiply, accounting for SB (-1) and mass (+1) eigenvalues
!
!   x = A x + b  ==>  A -> A'
!
  subroutine DeflateMatMult(npoint, v)
    integer                       :: npoint
    real(WP), dimension(npoint, 3) :: v

    real(WP), allocatable, dimension(:, :, :) :: vv

    integer  :: irbc
    type(t_RBC), pointer :: rbc
    integer :: nlat0, nlon0, nlat, nlon, ilat, ilon, nvert

    real(WP) :: alpha

    real(WP) :: ds, fac
    integer  :: p, m

    nlat = rbcs(1)%nlat
    nlon = rbcs(1)%nlon
    allocate (vv(nlat, nlon, 3))

    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      alpha = Bcoef(rbc%celltype)/Acoef(rbc%celltype) !COEF
!       alpha = (1.-viscRat(rbc%celltype))/(1.+viscRat(rbc%celltype)) !COEF

      nlat = rbc%nlat
      nlon = rbc%nlon
      p = 1 + (irbc - 1)*(nlat*nlon)

      do ilon = 1, nlon
        do ilat = 1, nlat
          vv(ilat, ilon, :) = v(p, :)
          p = p + 1
        end do
      end do

      fac = 0.
      do ilon = 1, nlon
        do ilat = 1, nlat
          ds = rbc%detj(ilat, ilon)*rbc%w(ilat)
          fac = fac + SUM(rbc%g(ilat, ilon, :)*rbc%a3(ilat, ilon, :))*ds
        end do ! ilon
      end do ! ilat
      vv = vv + alpha*rbc%a3*fac/rbc%area

      do m = 1, 6
        fac = 0.
        do ilon = 1, nlon
          do ilat = 1, nlat
            ds = rbc%detj(ilat, ilon)*rbc%w(ilat)
            fac = fac + SUM(rbc%g(ilat, ilon, :)*rbc%qq(ilat, ilon, :, m))*ds
          end do ! ilon
        end do ! ilat
        vv = vv - alpha*rbc%qq(:, :, :, m)*fac
      end do

      p = 1 + (irbc - 1)*(nlat*nlon)
      do ilon = 1, nlon
        do ilat = 1, nlat
          v(p, :) = vv(ilat, ilon, :)
          p = p + 1
        end do
      end do

    end do

    deallocate (vv)

  end subroutine DeflateMatMult

!!!!!!!!!!!!!
! Build the solid-body motion eigenvectors -- stored as rbc%qq
!
  subroutine DeflateGetSBModes

    integer  :: irbc
    type(t_RBC), pointer :: rbc
    integer :: nlat0, nlon0, nlat, nlon, ilat, ilon, nvert
    real(WP), dimension(3) :: XX
    real(WP), allocatable, dimension(:, :, :, :)  :: pp, qq
    real(WP), allocatable, dimension(:, :, :)    :: projpq
    real(WP)                                   :: ds

    integer :: n, m

    real(WP)  :: fac

    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      nlat0 = rbc%nlat0
      nlon0 = rbc%nlon0

      nlat = rbc%nlat
      nlon = rbc%nlon

      allocate (pp(nlat, nlon, 3, 6), qq(nlat, nlon, 3, 6), projpq(nlat, nlon, 3))

      call RBC_ComputeGeometry(rbc)

      ! the basic vectors are taken to be the coordinate directions for both
      ! the translation directions and the vectors omega which define rotations
      ! w \cross (x-xc)
      !
      ! All are includeded in the GS orthonomalization for simplicity...
      !  ... it would be straightforward to normalize the translation
      do ilon = 1, nlon
        do ilat = 1, nlat

          pp(ilat, ilon, :, 1) = (/1., 0., 0./)
          pp(ilat, ilon, :, 2) = (/0., 1., 0./)
          pp(ilat, ilon, :, 3) = (/0., 0., 1./)

          XX = rbc%x(ilat, ilon, :) - rbc%xc(:)

          pp(ilat, ilon, :, 4) = (/0., -XX(3), XX(2)/)
          pp(ilat, ilon, :, 5) = (/XX(3), 0., -XX(1)/)
          pp(ilat, ilon, :, 6) = (/-XX(2), XX(1), 0./)

        end do
      end do

      ! Gram-Schmit orthonomalization
      do m = 1, 6

        qq(:, :, :, m) = pp(:, :, :, m)

        do n = 2, m
          call Proj(rbc, nlat, nlon, pp(1, 1, 1, m), qq(1, 1, 1, n - 1), projpq)
          qq(:, :, :, m) = qq(:, :, :, m) - projpq
        end do

      end do

      ! normalize
      do m = 1, 6
        fac = 0.
        do ilat = 1, nlat
          do ilon = 1, nlon
            ds = rbc%detj(ilat, ilon)*rbc%w(ilat)
            fac = fac + SUM(qq(ilat, ilon, :, m)*qq(ilat, ilon, :, m))*ds
          end do ! ilon
        end do ! ilat
        rbc%qq(:, :, :, m) = qq(:, :, :, m)/SQRT(fac)
      end do

      deallocate (pp, qq, projpq)

    end do

    ! Test orthonormalilty
!!$   do irbc = 1,nrbc
!!$      rbc => rbcs(irbc)
!!$      do n = 1,6
!!$         do m = 1,6
!!$            fac = 0.
!!$            do ilat = 1, nlat
!!$               do ilon = 1, nlon
!!$                  ds = rbc%detj(ilat,ilon)*rbc%w(ilat)
!!$                  fac = fac + SUM(rbc%qq(ilat,ilon,:,n)*rbc%qq(ilat,ilon,:,m))*ds
!!$               end do ! ilon
!!$            end do ! ilat
!!$            print *,n,m,fac
!!$         end do
!!$      end do
!!$   end do
!!$   stop "orthonomality test"

  end subroutine DeflateGetSBModes

!**********************************************************************
! Compute RHS of the cell BIE
! Argument:
!  rhs -- RHS
  subroutine Compute_Rhs(rhs)
    real(WP) :: rhs(:)

    type(t_TargetList), pointer :: tlist
    real(WP), allocatable :: v(:, :)
    real(WP) :: c1, c2
    integer :: nPoint, ii
    integer :: ierr

    tlist => tlist_rbc
    c1 = 1./(4.*pi)   !COEF
    c2 = 0.           !COEF

    ! Allocate working arrays
    npoint = tlist%nPoint
    allocate (v(npoint, 3))

    ! Compute the boundary integrals
    v = 0.
!print*, PhysEwald, 'PhysEwald'
!print*, FourierEwald, 'F Ewald'
    if (PhysEwald) then
      call AddIntOnRbcs(c1, c2, tlist, v)
      !print*, 'add int on rbc'
      call AddIntOnWalls(c1, tlist, v)
      !print*, 'add int on walls'
    end if

    !print *,"NO RHS SMOOTH"
    if (FourierEwald) then
      call PME_Distrib_Source(c1, c2, rbcs, walls)
      !print*, 'dist'
      call PME_Transform
      !print*, 'transform'; stop
      call PME_Add_Interp_Vel(tlist, v)
    end if
    ! print*, 'got ewald'

    call TargetList_CollectArray(tlist, 3, v, MPI_COMM_WORLD)

    !print*, 'got targ list'
!    v = 0. ;  print *,"ZERO RHS"
    ! Add background velocities
    c1 = 2.       !COEF
    do ii = 1, npoint
      v(ii, :) = v(ii, :) + c1*vBkg(:)/tlist%Acoef(ii)   !COEF
!       v(ii,:) = v(ii,:) + c1*vBkg(:)/(1.+tlist%lam(ii))   !COEF
    end do
    !print*, 'got vels'
    if (Deflate) then
      !if (rootWorld) print *,"DEFLATING --- RHS"
      call DeflateRHS(npoint, v)
    end if

    ! Convert to spherical harmonic coefficients
    call Glob_Sph_Trans(v, rhs, PHYS_TO_FOUR)

    ! Deallocate working arrays
    deallocate (v)

  end subroutine Compute_Rhs

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

    integer :: nPoint, dof, nvert
    real(WP), allocatable :: u(:), b(:)
    real(WP), allocatable :: g(:, :), v(:, :)
    integer :: irbc, p, ii
    type(t_RBC), pointer :: rbc
    real(WP), allocatable :: GdetJ(:, :, :)
    real(WP) :: c1, c2
    type(t_SourceList), pointer :: slist
    type(t_TargetList), pointer :: tlist

    ! Set up parameters
    slist => slist_rbc
    tlist => tlist_rbc

    ! Allocate working arrays
    call VecGetSize(vec_u, dof, ierr)
    allocate (u(dof), b(dof))

    nPoint = tlist%nPoint
    allocate (g(nPoint, 3), v(nPoint, 3))

    ! Transform u-array to double-layer density on cell surfaces
    call My_VecGetValues(vec_u, u)

    call Glob_Sph_Trans(g, u, FOUR_TO_PHYS)

    ! Build the double-layer density on cell surfaces
    p = 0
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      nvert = rbc%nlat*rbc%nlon

      rbc%g = reshape(g(p + 1:p + nvert, :), shape(rbc%g))
      p = p + nvert

      call Rbc_BuildSurfaceSource(rbc, gFlag=.true.)
    end do ! irbc
    call SourceList_UpdateDensity(slist, rbcs=rbcs, updateG=.true.)

    ! Compute the off-diagonal term
    c1 = 0.          !COEF
    c2 = -1./(4.*Pi) !COEF

    v = 0.

    if (PhysEwald) then
      call AddIntOnRbcs(c1, c2, tlist, v)
    end if

!!$    print *,"NO MATMUL SMOOTH"
    if (FourierEwald) then
      call PME_Distrib_Source(c1, c2, cells=rbcs)
      call PME_Transform
      call PME_Add_Interp_Vel(tlist, v)
    end if

    call TargetList_CollectArray(tlist, 3, v, MPI_COMM_WORLD)

    ! Add the diagonal term
    v = v + g

    if (Deflate) then
      !if (rootWorld) print *,"DEFLATING --- MATMULT"
      call DeflateMatMult(npoint, v)
    end if

    call Glob_Sph_Trans(v, b, PHYS_TO_FOUR)
    call My_VecSetValues(vec_b, b)

    ! Deallocate working arrays
    deallocate (u, b)
    deallocate (g, v)

  end subroutine MyMatMult

  subroutine Proj(cell, nlat, nlon, pp, qq, projpq)
    type(t_RBC) :: cell
    integer     :: nlat, nlon
    real(WP)    :: pp(nlat, nlon, 3), qq(nlat, nlon, 3), projpq(nlat, nlon, 3)

    real(WP)    :: fn, fd, ds
    integer     :: ilat, ilon

    fn = 0.; fd = 0.
    do ilat = 1, nlat
    do ilon = 1, nlon
      ds = cell%detj(ilat, ilon)*cell%w(ilat)
      fn = fn + SUM(pp(ilat, ilon, :)*qq(ilat, ilon, :))*ds
      fd = fd + SUM(qq(ilat, ilon, :)*qq(ilat, ilon, :))*ds
    end do ! ilon
    end do ! ilat
    projpq = fn/fd*qq

  end subroutine Proj

!**********************************************************************
! Global spherical harmonic transform
!
! Arguments:
!  v -- RBC surface velocities
!  c -- spherical harmonic coefficients
!  direction -- PHYS_TO_FOUR -- v to coeff
!               FOUR_TO_PHYS -- coeff to v
! Note:
!  Data alignment:
!   v(i,:) -- velocity of the i-th global point
!   c -- spherical harmonic coefficients of v
!        stored in one-dimensional form
!
! Important:
!   For the spherical harmonic coefficiets array, the first dimension
!   is for longitudinal direction, while the second is for the colatitudinal
!   direction.
  subroutine Glob_Sph_Trans(v, c, direction)
    real(WP) :: v(:, :), c(:)
    integer :: direction

    integer :: pv, pc
    integer :: irbc
    type(t_RBC), pointer :: rbc
    integer :: nlat0, nlon0, nlat, nlon, ilat, ilon, nvert
    real(WP), allocatable :: vtmp(:, :, :), va(:, :, :), vb(:, :, :)

    pv = 0
    pc = 0

    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      nlat0 = rbc%nlat0
      nlon0 = rbc%nlon0

      nlat = rbc%nlat
      nlon = rbc%nlon

      nvert = nlat*nlon

      ! Allocate working arrays
      allocate (vtmp(nlat, nlon, 3), va(nlat, nlat, 3), vb(nlat, nlat, 3))

      if (direction == PHYS_TO_FOUR) then
        ! v to c
        vtmp = reshape(v(pv + 1:pv + nvert, :), shape(vtmp))
        pv = pv + nvert

        call ShAnalGau(nlat, nlon, 3, vtmp, size(vtmp, 1), size(vtmp, 2), &
                       va, vb, size(va, 1), size(va, 2), rbc%wshags)

        do ilon = 1, nlat0
        do ilat = ilon, nlat0
          c(pc + 1:pc + 3) = (/va(ilon, ilat, :)/)
          pc = pc + 3
        end do ! ilat
        end do ! ilon

        do ilon = 2, nlat0
        do ilat = ilon, nlat0
          c(pc + 1:pc + 3) = (/vb(ilon, ilat, :)/)
          pc = pc + 3
        end do ! ilat
        end do ! ilon

      else if (direction == FOUR_TO_PHYS) then
        ! c to v
        va = 0.
        vb = 0.
        do ilon = 1, nlat0
        do ilat = ilon, nlat0
          va(ilon, ilat, :) = c(pc + 1:pc + 3)
          pc = pc + 3
        end do ! ilat
        end do ! ilon

        do ilon = 2, nlat0
        do ilat = ilon, nlat0
          vb(ilon, ilat, :) = c(pc + 1:pc + 3)
          pc = pc + 3
        end do ! ilat
        end do ! ilon

        call ShSynthGau(nlat, nlon, 3, vtmp, size(vtmp, 1), size(vtmp, 2), &
                        va, vb, size(va, 1), size(va, 2), rbc%wshsgs)

        v(pv + 1:pv + nvert, :) = reshape(vtmp, (/nvert, 3/))
        pv = pv + nvert
      end if

      ! Deallocate working arrays
      deallocate (vtmp, va, vb)
    end do !irbc

  end subroutine Glob_Sph_Trans

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
    Vec :: x
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Testing Routine for assembling the A matrix
!
! u = alpha A x + b
!
! It is the eigensystem of A that is to be deflated

  SUBROUTINE BuildMat(dof, AA)
    integer                   :: dof
    real(WP), dimension(dof, dof)  :: AA

    real(WP), allocatable :: g(:, :), v(:, :)

    real(WP), dimension(dof)      :: u, b

    Mat :: mat_t
    Vec :: vec_u, vec_b

    integer :: nPoint
    type(t_SourceList), pointer :: slist
    type(t_TargetList), pointer :: tlist

    type(t_rbc), pointer :: rbc
    integer :: ierr
    integer :: i, j, l, m, p, q, irbc

    ! Update rbc surface geometry
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      call RBC_ComputeGeometry(rbc)
      call Rbc_BuildSurfaceSource(rbc, xFlag=.true.)
    end do ! irbc

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

  end SUBROUTINE BuildMat

end module ModVelSolver
