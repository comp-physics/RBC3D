! Solve the electric field on the surface of the drop (MHF)
module ModElectSolver

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

  implicit none

#include "../petsc_include.h"

  Mat :: mat_lhs_elect
  Vec :: vec_rhs_elect, vec_sol_elect
  KSP :: ksp_lhs_elect
  PC :: pc_lhs_elect

  private

  public :: Solve_RBC_Elect, &
  	MyMatMult_Elect, &
	  Glob_Sph_Trans_Elect, &
    Compute_ElectPotential, &
    Compute_ElectFieldPlus, &
    BuildMat_Elect !TEMP

  private :: Compute_Rhs_Elect, &
	My_VecGetValues, &
	My_VecSetValues

contains

!**********************************************************************
! Solve RBC surface velocities
! Note:
!   EE -- delta E.n on drop surface
  subroutine Solve_RBC_Elect(EE)
    real(WP) :: EE(:)

    logical,save :: solver_inited = .false.
    integer :: dof
    integer :: irbc
    type(t_RBC),pointer :: rbc
    integer :: nlat, nlon
    integer :: ii, p
    real(WP),allocatable,save :: rhs(:), sol(:)
    integer :: niter
    real(WP) :: residual
    integer :: ierr

    if (nrbc == 0) return

    ! Initialize the solver
    if (.not. solver_inited) then
      ! Calcuate matrix size
      dof = 0
      do irbc = 1, nrbc
        rbc => rbcs(irbc)
        dof = dof + (rbc%nlat0**2)
      end do ! irbc
	write(*,*) "Just calculated ELECTRIC dof, it is ",dof

      ! Set up petsc matrix-free GMRES solver
      allocate(rhs(dof), sol(dof) )
	write(*,*) "Set up GMRES solver Electric"

 ! MHF comments
      call VecCreateSeq(PETSC_COMM_SELF, dof, vec_rhs_elect, ierr) ! creat vec_rhs(Nlat0**2)
      call VecCreateSeq(PETSC_COMM_SELF, dof, vec_sol_elect, ierr) ! creat vec_sol_elect(Nlat0**2)

      call MatCreateShell(PETSC_COMM_SELF, dof, dof, dof, dof, PETSC_NULL_INTEGER, mat_lhs_elect, ierr)! creat matrix "mat_lhs_elect" for matrix-free methods
      call MatShellSetOperation(mat_lhs_elect, MATOP_MULT, MyMatMult_Elect, ierr) !defines matrix operation "MyMatMult" for the shell matrix "mat_lhs_elect".

      ! set up the Krylov Space method (GMRES)
      call KSPCreate(PETSC_COMM_SELF, ksp_lhs_elect, ierr) !create an interative Krylov Space context
      call KSPSetOperators(ksp_lhs_elect, mat_lhs_elect, mat_lhs_elect, SAME_NONZERO_PATTERN, ierr) ! "mat_lhs_elect" defines the linear system
      call KSPSetType(ksp_lhs_elect, KSPGMRES, ierr) ! type of the Krylov method: GMRES
      call KSPGetPC(ksp_lhs_elect, pc_lhs_elect, ierr) ! get a preconditioner
      call KSPSetInitialGuessNonzero(ksp_lhs_elect, PETSC_TRUE, ierr)  ! TRIAL
      call PCSetType(pc_lhs_elect, PCNONE, ierr)

      call KSPSetTolerances(ksp_lhs_elect, 1.D-9, &
      		PETSC_DEFAULT_DOUBLE_PRECISION, &
                PETSC_DEFAULT_DOUBLE_PRECISION, &
		500, ierr); print *,"SETTING HIGH MAX ITER/ERR"

      solver_inited = .true.
    end if
!write(*,*) "solver initiated"
    ! Compute rhs
    call Compute_Rhs_Elect(rhs)

!write(*,*) "rhs computed"
!    if (Deflate) then
!      print *,"DEFLATION for electric field not avilable (MHF)"
       !if (rootWorld) print *,"DEFLATING --- GET SB MODES"
  !     call DeflateGetSBModes
!    end if

    ! Solve the equation
    if (MAXVAL(QQRatio) == 1 .and. MINVAL(QQratio) == 1)then
      sol = rhs
    else

    ! Solve the equation
      call My_VecSetValues(vec_rhs_elect, rhs)

!      call My_VecSetValues(vec_sol_elect, sol); print *, "SETTING INITIAL GUESS"

      call KSPSolve(ksp_lhs_elect, vec_rhs_elect, vec_sol_elect, ierr)
      call My_VecGetValues(vec_sol_elect, sol)

      call KSPGetIterationNumber(ksp_lhs_elect, niter, ierr)
      call KSPGetResidualNorm(ksp_lhs_elect, residual, ierr)

      print*, 'iterations'
      if (rootWorld) then
	write (*, '(A,I5,A,ES12.2)') '(Elect) niter = ', niter, ' residual = ', residual
      end if
    end if

    ! Get drop surface EE
    call Glob_Sph_Trans_Elect(EE, sol, FOUR_TO_PHYS)

!    if (Deflate) then
!       if (rootWorld) print *,"DEFLATING for electric problem unavailable"
  !     call DeflateRecoverSol(size(EE,1),EE)
!    end if

  end subroutine Solve_RBC_Elect

!**********************************************************************
! Compute RHS of the cell BIE
! Argument:
!  EE_rhs -- RHS
  subroutine Compute_Rhs_Elect(EE_rhs)
    real(WP) :: EE_rhs(:)

    type(t_TargetList),pointer :: tlist
    real(WP),allocatable :: EE(:)
    real(WP) :: c1, c2
    integer :: nPoint, ii
    integer :: ierr

    tlist => tlist_rbc
    c1 = 1.   !COEF
    c2 = 0.           !COEF

    ! Allocate working arrays
    npoint = tlist%nPoint
    allocate(EE(npoint))


    ! Compute the boundary integrals
    EE = 0.
!print*, PhysEwald, 'PhysEwald'
!print*, FourierEwald, 'F Ewald'
    if (PhysEwald) then
      call AddIntOnRbcs_ElectFld(c1, c2, tlist, EE)
      !print*, 'add int on rbc'
!      call AddIntOnWalls(c1, tlist, EE)
      !print*, 'add int on walls'
    end if

    !print *,"NO RHS SMOOTH"
  !  if (FourierEwald) then
  !!    print *,"NO RHS SMOOTH"
  !    call PME_Distrib_Source(c1, c2, rbcs, walls)
      !print*, 'dist'
  !    call PME_Transform
      !print*, 'transform'; stop
  !    call PME_Add_Interp_Vel(tlist, EE)
  !  end if
   ! print*, 'got ewald'


    call TargetList_CollectArray(tlist, 1, EE, MPI_COMM_WORLD)


  !  ! Add background velocities
  !  c1 = 2.       !COEF
  !  do ii = 1,npoint
  !     v(ii,:) = v(ii,:) + c1*vBkg(:)/tlist%Acoef(ii)   !COEF
!  !     v(ii,:) = v(ii,:) + c1*vBkg(:)/(1.+tlist%lam(ii))   !COEF
    !end do
    !print*, 'got vels'
  !  if (Deflate) then
  !     !if (rootWorld) print *,"DEFLATING --- RHS"
  !     call DeflateRHS(npoint,EE)
  !  end if

    ! Convert to spherical harmonic coefficients
    call Glob_Sph_Trans_Elect(EE, EE_rhs, PHYS_TO_FOUR)

    ! Deallocate working arrays
    deallocate(EE)

  end subroutine Compute_Rhs_Elect




!**********************************************************************
! Matmul for matrix-free GMRES solver
! Arguments:
!  lhs_petsc -- lhs matrix
!  vec_u, vec_b -- b = lhs*u
!  ierr -- error code
!  Note:
!  u and b are in spherical harmonic wave space. (MHF)

  subroutine MyMatMult_Elect(lhs_petsc, vec_u, vec_b, ierr)
    Mat :: lhs_petsc
    Vec :: vec_u, vec_b
    integer :: ierr

    integer :: nPoint, dof, nvert
    real(WP),allocatable :: u(:), b(:)
    real(WP),allocatable :: delEn(:), EE(:)
    integer :: irbc, p, ii
    type(t_RBC),pointer :: rbc
    real(WP),allocatable :: GdetJ(:,:,:)
    real(WP) :: c1, c2
    type(t_SourceList),pointer :: slist
    type(t_TargetList),pointer :: tlist

    ! Set up parameters
    slist => slist_rbc
    tlist => tlist_rbc

    ! Allocate working arrays
    call VecGetSize(vec_u, dof, ierr)
    allocate(u(dof), b(dof) )

    nPoint = tlist%nPoint
    allocate(delEn(nPoint), EE(nPoint))

    ! Transform u-array to double-layer density on drop surfaces (electric problem)
    call My_VecGetValues(vec_u, u) ! puts vec_u into u

    call Glob_Sph_Trans_Elect(delEn, u, FOUR_TO_PHYS)!given u, find delEn in physical space

    ! Build the double-layer density on drop surfaces (electric problem)
    p = 0
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      nvert = rbc%nlat*rbc%nlon

      rbc%delEn = reshape(delEn(p+1:p+nvert), shape(rbc%delEn))

      call Rbc_SphProject(rbc, 1, rbc%delEn) ! new 22-Jan-22

      p = p + nvert

      call Rbc_BuildSurfaceSource(rbc, delEnFlag=.true.)
    end do ! irbc
    call SourceList_UpdateDensity(slist, rbcs=rbcs, UpdateDelEn=.true.)

    ! Compute the off-diagonal term
    c1 = 0.          !COEF
    c2 = 1.          !COEF

    EE = 0.

    if (PhysEwald) then
      call AddIntOnRbcs_ElectFld(c1, c2, tlist, EE)
    end if

!!$    print *,"NO MATMUL SMOOTH"
  !  if (FourierEwald) then
  !    print *,"SMOOTH Ewald in MyMATMUL_Elect!!! (MHF)"
  !    call PME_Distrib_Source(c1, c2, cells=rbcs)
  !    call PME_Transform
  !    call PME_Add_Interp_Vel(tlist, v)
  !  end if

    call TargetList_CollectArray(tlist, 1, EE, MPI_COMM_WORLD)

    ! Add the diagonal term
    EE = EE + delEn

  !  if (Deflate) then
  !     !if (rootWorld) print *,"DEFLATING --- MATMULT"
  !     call DeflateMatMult(npoint,EE)
  !  end if

    call Glob_Sph_Trans_Elect(EE, b, PHYS_TO_FOUR)
    call My_VecSetValues(vec_b, b) !puts b into vec_b

    ! Deallocate working arrays
    deallocate(u, b)
    deallocate(delEn, EE)

  end subroutine MyMatMult_Elect



  subroutine Proj(cell,nlat,nlon,pp,qq,projpq)
    type(t_RBC) :: cell
    integer     :: nlat,nlon
    real(WP)    :: pp(nlat,nlon,3),qq(nlat,nlon,3),projpq(nlat,nlon,3)

    real(WP)    :: fn,fd,ds
    integer     :: ilat,ilon

    fn = 0.; fd = 0.
    do ilat = 1, nlat
    do ilon = 1, nlon
      ds = cell%detj(ilat,ilon)*cell%w(ilat)
      fn = fn + SUM(pp(ilat,ilon,:)*qq(ilat,ilon,:))*ds
      fd = fd + SUM(qq(ilat,ilon,:)*qq(ilat,ilon,:))*ds
    end do ! ilon
    end do ! ilat
    projpq = fn/fd*qq

  end subroutine Proj


!**********************************************************************
! Global spherical harmonic transform (MHF)
!
! Arguments:
!  EE -- n.delta E on drop surface
!  c -- spherical harmonic coefficients
!  direction -- PHYS_TO_FOUR -- EE to coeff
!               FOUR_TO_PHYS -- coeff to EE
! Note:
!  Data alignment:
!   EE(i,:) -- Electric field jumo of the i-th global point
!   c -- spherical harmonic coefficients of EE
!        stored in one-dimensional form
!
! Important:
!   For the spherical harmonic coefficiets array, the first dimension
!   is for longitudinal direciton, while the second is for the colatitudinal
!   direction.
  subroutine Glob_Sph_Trans_Elect(EE, c, direction)
    real(WP) :: EE(:), c(:)
    integer :: direction

    integer :: pv, pc
    integer :: irbc
    type(t_RBC),pointer :: rbc
    integer :: nlat0, nlon0, nlat, nlon, ilat, ilon, nvert
    real(WP),allocatable :: EEtmp(:,:), EEa(:,:), EEb(:,:)

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
      allocate(EEtmp(nlat,nlon), EEa(nlat,nlat), EEb(nlat,nlat) )

      if (direction == PHYS_TO_FOUR) then
        ! EE to c
	EEtmp = reshape(EE(pv+1:pv+nvert), shape(EEtmp) )
	pv = pv + nvert

	call ShAnalGau(nlat, nlon, 1, EEtmp, size(EEtmp,1), size(EEtmp,2), &
			EEa, EEb, size(EEa,1), size(EEa,2), rbc%wshags )

	do ilon = 1, nlat0
	do ilat = ilon, nlat0
	  c(pc+1) =  EEa(ilon,ilat)
	  pc = pc + 1
	end do ! ilat
	end do ! ilon

	do ilon = 2, nlat0
	do ilat = ilon, nlat0
	  c(pc+1) =  EEb(ilon,ilat)
	  pc = pc + 1
	end do ! ilat
	end do ! ilon

      else if (direction == FOUR_TO_PHYS) then
        ! c to EE
	EEa = 0.
	EEb = 0.
	do ilon = 1, nlat0
	do ilat = ilon, nlat0
          EEa(ilon,ilat) = c(pc+1)
	  pc = pc + 1
	end do ! ilat
	end do ! ilon

	do ilon = 2, nlat0
	do ilat = ilon, nlat0
          EEb(ilon,ilat) = c(pc+1)
	  pc = pc + 1
	end do ! ilat
	end do ! ilon

	call ShSynthGau(nlat, nlon, 1, EEtmp, size(EEtmp,1), size(EEtmp,2), &
			EEa, EEb, size(EEa,1), size(EEa,2), rbc%wshsgs )

	EE(pv+1:pv+nvert) = reshape(EEtmp, (/nvert/))
        pv = pv + nvert
      end if

      ! Deallocate working arrays
      deallocate(EEtmp, EEa, EEb)
    end do !irbc

  end subroutine Glob_Sph_Trans_Elect

!**********************************************************************
  subroutine My_VecGetValues(x, a)
    Vec x
    real(WP) :: a(:)

    integer,allocatable :: ix(:)
    integer :: n, i, ierr

    ! Allocate working arrays
    call VecGetSize(x, n, ierr) !PETSC: dtermines the global number of elements
    allocate(ix(n))
    ix = (/ (i, i=0,n-1) /)

!PETSC: get n elements from x starting from index ix, and put them into a
    call VecGetValues(x, n, ix, a, ierr)

    ! Deallocate working arrays
    deallocate(ix)

  end subroutine My_VecGetValues

!**********************************************************************
  subroutine My_VecSetValues(x, a)
    Vec :: x
    real(WP) :: a(:)

    integer,allocatable :: ix(:)
    integer :: n, i, ierr

    ! Allocate working arrays
    call VecGetSize(x, n, ierr)
    allocate(ix(n))
    ix = (/ (i, i=0,n-1) /)

    call VecSetValues(x, n, ix, a, INSERT_VALUES, ierr)
    call VecAssemblyBegin(x,ierr)
    call VecAssemblyEnd(x, ierr)

    ! Deallocate working arrays
    deallocate(ix)

  end subroutine My_VecSetValues

!**********************************************************************


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Testing Routine for assembling the A matrix
!
! u = alpha A x + b
!
! It is the eigensystem of A that is to be deflated

  SUBROUTINE BuildMat_Elect(dof, AA)
    integer                   :: dof
    real(WP), dimension(dof,dof)  :: AA

    real(WP),allocatable :: delEn(:), EE(:)!g(:,:), v(:,:)

    real(WP), dimension(dof)      :: u,b

    Mat :: mat_t
    Vec :: vec_u, vec_b

    integer :: nPoint
    type(t_SourceList),pointer :: slist
    type(t_TargetList),pointer :: tlist

    type(t_rbc),pointer :: rbc
    integer :: ierr
    integer :: i,j,l,m,p,q, irbc

   ! Update rbc surface geometry
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      call RBC_ComputeGeometry_MHF2(rbc)
      call Rbc_BuildSurfaceSource(rbc, xFlag=.true., delEnFlag=.true.)
    end do ! irbc

    call VecCreateSeq(PETSC_COMM_SELF, dof, vec_u, ierr)
    call VecCreateSeq(PETSC_COMM_SELF, dof, vec_b, ierr)

    call MatCreateShell(PETSC_COMM_SELF, dof, dof, dof, &
                        dof, PETSC_NULL_INTEGER, mat_t, ierr)
    call MatShellSetOperation(mat_t, MATOP_MULT, MyMatMult_Elect, ierr)

    do i = 1,dof
       u = 0.
       u(i) = 1.

       if (rootWorld) print *,i, " of ", dof

       call My_VecSetValues(vec_u, u)
       call MatMult(mat_t, vec_u, vec_b, ierr)
       call My_VecGetValues(vec_b, b)

       AA(:,i) = b

    end do

  end SUBROUTINE BuildMat_Elect

!**********************************************************
! This subroutine computes electric potential phi based on
! the BIE for delEn:
! phi = -x0.E_ext +\int_c delEn*G
!
! Note: PPhi   ---- in physical space
!
! by MHF

  subroutine Compute_ElectPotential(PPhi)

    real(WP) ::  PPhi(:)

    type(t_TargetList),pointer :: tlist
    ! type(t_RBC),pointer :: rbc
    integer :: nPoint, ii
    integer :: ierr

    tlist => tlist_rbc

    ! Compute the boundary integrals
    PPhi = 0.

    if (PhysEwald) then
      call AddIntOnRbcs_ElectPot(tlist,PPhi)

    end if

    call TargetList_CollectArray(tlist, 1, PPhi, MPI_COMM_WORLD)



    ! Convert to spherical harmonic coefficients
  !  call Glob_Sph_Trans_Elect(PPhi_p, PPhi, PHYS_TO_FOUR)

    ! Deallocate working arrays
!    deallocate(PPhi_p)

  end subroutine Compute_ElectPotential


  !**********************************************************
  ! This subroutine computes electric potential phi based on
  ! the BIE for delEn:
  ! EF_plus = E_ext +\int_c delEn* grad0G + 0.5*delEn *n
  !
  ! Note: PPhi   ---- in physical space
  !
  ! by MHF

    subroutine Compute_ElectFieldPlus(EF_plus)

      real(WP) ::  EF_plus(:,:)

      type(t_TargetList),pointer :: tlist
      ! type(t_RBC),pointer :: rbc
      integer :: nPoint, ii
      integer :: ierr

      tlist => tlist_rbc

      ! Compute the boundary integrals
      EF_plus = 0.

      if (PhysEwald) then
        call AddIntOnRbcs_ElectFieldPlus(tlist,EF_plus)

      end if


      call TargetList_CollectArray(tlist, 3, EF_plus, MPI_COMM_WORLD)



    end subroutine Compute_ElectFieldPlus

end module ModElectSolver
