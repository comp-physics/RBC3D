! Surface integration
module ModIntOnWalls

  use ModDataTypes
  use ModDataStruct
  use ModConf, rc => rc_Ewd
  use ModData
  use ModQuadRule
  use ModHashTable
  use ModBasicMath
  use ModEwaldFunc

  implicit none

  private

  public :: AddIntOnWalls, &
            SingIntOnWall, &
            PrepareSingIntOnWall, &
            Tri_Int_Regular, &
            Tri_Int_Duffy, &
            MinDistToTri

contains

!**********************************************************************
! Add integration from the walls
! Arguments:
!  c1 -- coefficient in front of the surface integrals
!  tlist -- the target list
! Note:
!  -- It is required that either all or none of the target points
!     are wall mesh points.
  subroutine AddIntOnWalls(c1, tlist, v)
    real(WP) :: c1
    type(t_TargetList) :: tlist
    real(WP) :: v(:, :)

    integer :: first_surfId
    type(t_SourceList), pointer :: slist
    type(t_Wall), pointer :: wall
    integer :: i, j, i1, i2, i3, j1, j2, j3
    integer :: p, iwall, iele, ivert, l
    real(WP) :: xele(3, 3), fele(3, 3)
    real(WP) :: xi(3), xx(3), rr, dv(3), s0, t0
    real(WP), allocatable :: vtmp(:, :)
    character(*), parameter :: func_name = 'AddIntOnWalls'
    integer :: ierr

    if (nwall == 0) return

    ! Self-interactions
    first_surfId = tlist%indx(1, 0)
    if (first_surfId == walls(1)%ID) then
      p = 0

      do iwall = 1, nwall
        wall => walls(iwall)
        allocate (vtmp(wall%nvert, 3))

        call SingIntOnWall(c1, wall, vtmp)
        do i = 1, wall%nvert
          v(p + i, :) = v(p + i, :) + vtmp(i, :)/tlist%Acoef(p + i)  !COEF
          ! v(p+i,:) = v(p+i,:) + vtmp(i,:)/(1.+tlist%lam(p+i))  !COEF
        end do

        deallocate (vtmp)

        p = p + wall%nvert
      end do ! iwall

      ! If there is only one wall, then all are self-interactions
      if (nwall == 1) return
    end if

    ! Non-self-interactions
    slist => slist_wall

    do i = 1, tlist%nPoint
      if (.not. tlist%active(i)) cycle
      xi = tlist%x(i, :)

      call HashTable_Index(slist%Nc, slist%iLbNc, xi, i1, i2, i3)
      do j1 = max(i1 - 1, 0), min(i1 + 1, slist%Nc(1) + 1)
      do j2 = max(i2 - 1, 0), min(i2 + 1, slist%Nc(2) + 1)
      do j3 = max(i3 - 1, 0), min(i3 + 1, slist%Nc(3) + 1)
        j = slist%hoc(j1, j2, j3)

        do while (j > 0)
          if (tlist%indx(i, 0) == slist%indx(j, 0)) goto 999

          iwall = slist%indx(j, 0) - walls(1)%ID + 1
          wall => walls(iwall)
          iele = slist%indx(j, 1)
          do l = 1, 3
            ivert = wall%e2v(iele, l)
            xele(l, :) = wall%x(ivert, :)
            fele(l, :) = wall%f(ivert, :)
          end do ! l

          ! Translate xele close to xi
          xx = nint((xi - xele(1, :))*iLb)*Lb
          do l = 1, 3
            xele(l, :) = xele(l, :) + xx
          end do ! l

          rr = MinDistToTri(xi, xele, s0, t0)
          if (rr > rc) goto 999

          if (rr < wall%epsDist(iele)) then
            call Tri_Int_Duffy(xele, fele, xi, s0, t0, dv)
          else
            call Tri_Int_Regular(xele, fele, xi, dv)
          end if
          v(i, :) = v(i, :) + c1*dv/tlist%Acoef(i)  !COEF
!     v(i,:) = v(i,:) + c1*dv/(1.+tlist%lam(i))  !COEF

999       j = slist%next(j)
        end do ! while
      end do ! j3
      end do ! j2
      end do ! j1
    end do ! i

  end subroutine AddIntOnWalls

!**********************************************************************
! Compute the self-interaction on the wall
! Arguments:
!  c1 -- coefficient before the single layer potential
!  wall -- the wall
!  v -- the velocities
  subroutine SingIntOnWall(c1, wall, v)
    real(WP) :: c1
    type(t_Wall) :: wall
    real(WP) :: v(:, :)

#include "../petsc_include.h"
    integer :: nrow, i
    integer, allocatable :: irows(:)
    Vec :: f_vec, lhsf_vec
    integer :: ierr

    nrow = 3*wall%nvert
    allocate (irows(nrow))
    call VecCreateSeq(PETSC_COMM_SELF, nrow, f_vec, ierr)
    call VecCreateSeq(PETSC_COMM_SELF, nrow, lhsf_vec, ierr)

    irows = (/(i, i=0, nrow - 1)/)
    call VecSetValues(f_vec, nrow, irows, wall%f, INSERT_VALUES, ierr)
    call VecAssemblyBegin(f_vec, ierr)
    call MatMult(wall%lhs, f_vec, lhsf_vec, ierr)
    call VecGetValues(lhsf_vec, nrow, irows, v, ierr)

    v = c1*v

    call VecDestroy(f_vec, ierr)
    call vecDestroy(lhsf_vec, ierr)
    deallocate (irows)

  end subroutine SingIntOnWall

! !**********************************************************************
! ! Compute the self-interaction on the wall
! ! Arguments:
! !  c1 -- coefficient before the single layer potential
! !  wall -- the wall
! !  v -- the velocities
!   subroutine SingIntOnWall(c1, wall, v)
!     real(WP) :: c1
!     type(t_Wall) :: wall
!     real(WP) :: v(:, :)
!     real(WP), allocatable :: v1D(:), f1D(:)

! #include "../petsc_include.h"
!     integer :: nrow, i, j
!     integer, allocatable :: irows(:)
!     Vec :: f_vec, lhsf_vec
!     integer :: ierr
!     integer, save :: count
!     integer :: Mat_m, Mat_n

!     count = count + 1

!     ! if (rootWorld) print *, "count: ", count

!     nrow = 3*wall%nvert
!     ! if (rootWorld) then
!     !   print *, "SINGINTONWALL nrow: ", nrow
!     !   print *, "size(wall%f, 1), size(wall%f, 2)", size(wall%f, 1), size(wall%f, 2)
!     !   do j = 1, 3
!     !     print *, wall%f(j, :)
!     !   end do
!     !   print *, "size(v, 1), size(v, 2)", size(v, 1), size(v, 2)
!     !   call MatGetSize(wall%lhs, Mat_m, Mat_n, ierr)
!     !   print *, "Mat_m, Mat_n", Mat_m, Mat_n
!     ! end if
    
!     allocate (irows(nrow))
!     allocate (f1D(nrow))
!     allocate (v1D(nrow))

!     f1D = reshape(wall%f, (/nrow/))

!     call VecCreateSeq(PETSC_COMM_SELF, nrow, f_vec, ierr)
!     call VecCreateSeq(PETSC_COMM_SELF, nrow, lhsf_vec, ierr)

!     irows = (/(i, i=0, nrow - 1)/)
!     call VecSetValues(f_vec, nrow, irows, f1D, INSERT_VALUES, ierr)
!     call VecAssemblyBegin(f_vec, ierr)

!     ! if (rootWorld) call VecView(f_vec, PETSC_VIEWER_STDOUT_SELF, ierr)

!     call MatMult(wall%lhs, f_vec, lhsf_vec, ierr)
!     call VecGetValues(lhsf_vec, nrow, irows, v1D, ierr)

!     v = reshape(v1D, (/wall%nvert, 3/))

!     v = c1*v

!     call VecDestroy(f_vec, ierr)
!     call vecDestroy(lhsf_vec, ierr)
!     deallocate (irows, f1D, v1D)

!   end subroutine SingIntOnWall

!**********************************************************************
! Prepare for the self-integration on the wall
! Arguments:
!  wall --
!
! Note:
!  The slist_wall and tlist_wall must be built before calling this subroutine.
  subroutine PrepareSingIntOnWall(wall)
    type(t_wall) :: wall

#include "../petsc_include.h"
    type(t_targetlist), pointer :: tlist
    type(t_sourcelist), pointer :: slist
    integer :: nvert, nrow
    logical, allocatable :: active(:)
    integer, allocatable :: nnz(:)
    integer :: p, i, j, i1, i2, i3, j1, j2, j3
    integer :: iele, ivert, l
    real(WP) :: xele(3, 3), fele(3, 3), s0, t0
    real(WP) :: xi(3), xj(3), xx(3), rr
    real(WP) :: rhs(3), lhs(3, 3, 3), values(9)
    integer :: irows(3), icols(3)
    integer :: ierr

    tlist => tlist_wall
    slist => slist_wall

    ! Allocate working arrays
    nvert = wall%nvert
    nrow = 3*nvert
    allocate (active(nvert))
    allocate (nnz(nrow))

    ! Set the active flag
    do p = 1, tlist%nPoint
      if (tlist%indx(p, 0) == wall%id) exit
    end do ! p
    active = tlist%active(p:p + nvert - 1)

    ! Estimate the number of non-zero elements of each row
    nnz = 0
    do i = 1, nvert
      if (.not. active(i)) cycle

      xi = wall%x(i, :)
      call HashTable_Index(slist%Nc, slist%iLbNc, xi, i1, i2, i3)

      do j1 = max(i1 - 1, 0), min(i1 + 1, slist%Nc(1) + 1)
      do j2 = max(i2 - 1, 0), min(i2 + 1, slist%Nc(2) + 1)
      do j3 = max(i3 - 1, 0), min(i3 + 1, slist%Nc(3) + 1)
        j = slist%hoc(j1, j2, j3)

        do while (j > 0)
          if (slist%indx(j, 0) .ne. wall%Id) goto 888

          xj = slist%x(j, :)
          xx = xi - xj; xx = xx - nint(xx*iLb)*Lb
          rr = sqrt(sum(xx*xx))

          if (rr < rc) then
            nnz(i) = nnz(i) + 3
            nnz(i + nvert) = nnz(i)
            nnz(i + 2*nvert) = nnz(i)
          end if

888       j = slist%next(j)
        end do ! j
      end do ! j3
      end do ! j2
      end do ! j1
    end do ! i

    nnz = nnz*2    ! a very conservative estimate
    call MatCreateSeqAIJ(PETSC_COMM_SELF, nrow, nrow, 0, nnz, wall%lhs, ierr)

    ! Assemble the matrix
    do i = 1, nvert
      if (.not. active(i)) cycle

      xi = wall%x(i, :)
      call HashTable_Index(slist%Nc, slist%iLbNc, xi, i1, i2, i3)

      do j1 = max(i1 - 1, 0), min(i1 + 1, slist%Nc(1) + 1)
      do j2 = max(i2 - 1, 0), min(i2 + 1, slist%Nc(2) + 1)
      do j3 = max(i3 - 1, 0), min(i3 + 1, slist%Nc(3) + 1)
        ! j is the source element
        j = slist%hoc(j1, j2, j3)

        do while (j > 0)
          if (slist%indx(j, 0) .ne. wall%id) goto 999

          iele = slist%indx(j, 1)

          do l = 1, 3
            ivert = wall%e2v(iele, l)
            xele(l, :) = wall%x(ivert, :)
          end do ! l

          ! Translate the triangle close to xi
          xx = nint((xi - xele(1, :))*iLb)*Lb
          do l = 1, 3
            xele(l, :) = xele(l, :) + xx
            fele(l, :) = 0.0     ! dummy
          end do ! l

          rr = MinDistToTri(xi, xele, s0, t0)
          if (rr > rc) goto 999

          if (rr > wall%epsDist(iele)) then
            call Tri_Int_Regular(xele, fele, xi, rhs, lhs)
          else
            call Tri_Int_Duffy(xele, fele, xi, s0, t0, rhs, lhs)
          end if

          irows = (/i, i + nvert, i + 2*nvert/) - 1
          do l = 1, 3
            ivert = wall%e2v(iele, l)
            icols = (/ivert, ivert + nvert, ivert + 2*nvert/) - 1
            values = reshape(transpose(lhs(l, :, :)), (/9/))
            ! PETSc InsertMode param wrong?
            call MatSetValues(wall%lhs, 3, irows, 3, icols, values, ADD_VALUES, ierr)
          end do ! l

999       j = slist%next(j)
        end do ! while
      end do ! j3
      end do ! j2
      end do ! j1
    end do ! i

    call MatAssemblyBegin(wall%lhs, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(wall%lhs, MAT_FINAL_ASSEMBLY, ierr)

    deallocate (active, nnz)

  end subroutine PrepareSingIntOnWall

!**********************************************************************
! Compute the regular surface integral over a triangle
! Arguments:
!  x(i,:), f(i,:) -- triangular vertex coordinates and force densities
!  xtar -- target point
!  rhs -- integration
!  lhs -- influence matrix
! Note:
!  -- rhs = \sum_i lhs(i,:,:)*f(i,:)
  subroutine Tri_Int_Regular(x, f, xtar, rhs, lhs)
    real(WP) :: x(3, 3), f(3, 3)
    real(WP) :: xtar(3)
    real(WP) :: rhs(3)
    real(WP), optional :: lhs(3, 3, 3)

    real(WP) :: detJ
    type(t_GaussQuad2D), pointer :: GQ
    integer :: iGq, ii, jj
    real(WP) :: s, t, xGq(3), fGq(3), dsGq, xx(3), rr, EA, EB
    real(WP) :: lhsTmp(3, 3)

    detJ = 2*TriArea(x)

    GQ => gqTri7
    rhs = 0.
    if (present(lhs)) lhs = 0.

    do iGq = 1, GQ%n
      s = GQ%rst(iGq, 1)
      t = GQ%rst(iGQ, 2)

      xGq = (1.-s - t)*x(1, :) + s*x(2, :) + t*x(3, :)
      fGq = (1.-s - t)*f(1, :) + s*f(2, :) + t*f(3, :)
      dsGq = GQ%w(iGQ)*detJ
      fGq = dsGq*fGq

      xx = xTar - xGq
      rr = sqrt(sum(xx**2))
      call EwaldCoeff_SL(rr, EA, EB)

      rhs = rhs + (EA*xx*dot_product(xx, fGq) + EB*fGq)

      if (present(lhs)) then
        forall (ii=1:3, jj=1:3) lhsTmp(ii, jj) = EA*xx(ii)*xx(jj)
        forall (ii=1:3) lhsTmp(ii, ii) = lhsTmp(ii, ii) + EB
        lhsTmp = dsGq*lhsTmp

        lhs(1, :, :) = lhs(1, :, :) + (1 - s - t)*lhsTmp
        lhs(2, :, :) = lhs(2, :, :) + s*lhsTmp
        lhs(3, :, :) = lhs(3, :, :) + t*lhsTmp
      end if
    end do ! iGQ

  end subroutine Tri_Int_Regular

!**********************************************************************
! Use Duffy's rule to compute the surface integral over a triangle
! Arguments:
!  x(i,:), f(i,:) -- triangular vertex coordinates and force densities
!  xtar -- target point
!  (s0, t0) -- reference coordinate of the nearest point to xtar
!  rhs -- rhs
!  lsh -- lhs matrix
  subroutine Tri_Int_Duffy(x, f, xtar, s0, t0, rhs, lhs)
    real(WP) :: x(3, 3), f(3, 3)
    real(WP) :: xtar(3), s0, t0
    real(WP) :: rhs(3)
    real(WP), optional :: lhs(3, 3, 3)

    real(WP) :: x0(3), f0(3), x1(3), f1(3), x2(3), f2(3), detJ0
    ! The quadrature rule is for integrating over [0,1]x[0,1]
    ! nGq, rGq, wGq -- 1-D quadrature rule
    logical, save :: INITED = .false.
    ! integer,save :: nGq
    integer, parameter :: nGq = 4
    real(WP), save :: rGq(4), wGq(4)
    !$omp threadprivate(INITED, nGq, rGq, wGq)
    real(WP) :: normal(3), detJ
    integer :: n, i, j, ii, jj
    real(WP) :: s, t, xGq(3), fGq(3), dsGq, xx(3), rr, EA, EB
    real(WP) :: s1, t1, s2, t2, sGlb, tGlb
    real(WP) :: lhsTmp(3, 3)

    ! Init the quadrature rule
    if (.not. INITED) then
      !  nGq = 4
      !  allocate(rGq(nGq), wGq(nGq))
      write (*, *) "***Initializing the quadrature rule in PSIOW***"
      call GauLeg(0._WP, 1._WP, nGq, rGq, wGq)
      INITED = .true.
    end if

    detJ0 = 2*TriArea(x)
    x0 = (1 - s0 - t0)*x(1, :) + s0*x(2, :) + t0*x(3, :)
    f0 = (1 - s0 - t0)*f(1, :) + s0*f(2, :) + t0*f(3, :)

    rhs = 0.
    if (present(lhs)) lhs = 0.

    ! Loop over the 3 sub-triangles
    do n = 1, 3
      x1 = x(n, :)
      x2 = x(mod(n, 3) + 1, :)

      f1 = f(n, :)
      f2 = f(mod(n, 3) + 1, :)

      normal = CrossProd(x1 - x0, x2 - x0)
      detJ = sqrt(sum(normal**2))

      do i = 1, nGq
      do j = 1, nGq
        s = rGq(i)
        t = s*rGq(j)

        xGq = (1.-s)*x0 + (s - t)*x1 + t*x2
        fGq = (1.-s)*f0 + (s - t)*f1 + t*f2
        dsGq = wGq(i)*wGq(j)*detJ*s
        fGq = dsGq*fGq

        xx = xTar - xGq
        rr = sqrt(sum(xx**2))
        call EwaldCoeff_SL(rr, EA, EB)

        rhs = rhs + (EA*xx*dot_product(xx, fGq) + EB*fGq)

        if (present(lhs)) then
          forall (ii=1:3, jj=1:3) lhsTmp(ii, jj) = EA*xx(ii)*xx(jj)
          forall (ii=1:3) lhsTmp(ii, ii) = lhsTmp(ii, ii) + EB
          lhsTmp = dsGq*lhsTmp

          ! Find the reference coordinate in the original triangle
          select case (n)
          case (1)
            s1 = 0.; t1 = 0.
            s2 = 1.; t2 = 0.
          case (2)
            s1 = 1.; t1 = 0.
            s2 = 0.; t2 = 1.
          case (3)
            s1 = 0.; t1 = 1.
            s2 = 0.; t2 = 0.
          end select
          sGlb = (1.-s)*s0 + (s - t)*s1 + t*s2
          tGlb = (1.-s)*t0 + (s - t)*t1 + t*t2

          lhs(1, :, :) = lhs(1, :, :) + (1.-sGlb - tGlb)*lhsTmp
          lhs(2, :, :) = lhs(2, :, :) + sGlb*lhsTmp
          lhs(3, :, :) = lhs(3, :, :) + tGlb*lhsTmp
        end if ! present lhs

      end do ! j
      end do ! i
    end do ! n

  end subroutine Tri_Int_Duffy

!**********************************************************************
! Find the shortest distance between a point and a triangle
! Arguments:
!  xTar -- the target point
!  x(i,:) -- i-th triangule vertex
!  s0, t0 -- the reference coordinate of the closest point
!  x0 -- the physical coordinate of the closest point
!
! Note:
!  x0 = (1 - s0 - t0)*x(1,:) + s0*x(2,:) + t0*x(3,:)
!
! Reference:
!  David Eberly, "Distance Between Point and Triangle in 3D"
  function MinDistToTri(xTar, x, s0, t0, x0)
    real(WP) :: MinDistToTri
    real(WP) :: xTar(3), x(3, 3)
    real(WP), optional :: s0, t0, x0(3)

    real(WP) :: s, t
    real(WP) :: x12(3), x13(3), a, b, c, det, invDet
    real(WP) :: x1Tar(3), d, e, f
    integer :: whichRegion

    ! Compute the distance of x0 to the triangle, and its projection
    ! on the triangle
    x12 = x(2, :) - x(1, :)
    x13 = x(3, :) - x(1, :)
    a = dot_product(x12, x12)
    b = dot_product(x12, x13)
    c = dot_product(x13, x13)

    det = a*c - b*b
    invDet = 1./det

    x1Tar = x(1, :) - xTar
    d = dot_product(x12, x1Tar)
    e = dot_product(x13, x1Tar)
    f = dot_product(x1Tar, x1Tar)

    s = b*e - c*d
    t = b*d - a*e

    if (s + t <= det) then
      if (s < 0) then
        if (t < 0) then
          whichRegion = 4
        else
          whichRegion = 3
        end if
      else if (t < 0) then
        whichRegion = 5
      else
        whichRegion = 0
      end if
    else
      if (s < 0) then
        whichRegion = 2
      else if (t < 0) then
        whichRegion = 6
      else
        whichRegion = 1
      end if
    end if

    if (whichRegion == 2) then
      if (-(c + e) < 0) then
        whichRegion = 3
      else
        whichRegion = 1
      end if
    else if (whichRegion == 4) then
      if (d < 0) then
        whichRegion = 5
      else
        whichRegion = 3
      end if
    else if (whichRegion == 6) then
      if (b + e - a - d < 0) then
        whichRegion = 1
      else
        whichRegion = 5
      end if
    end if

    select case (whichRegion)
    case (0)
      s = invDet*s
      t = invDet*t

    case (1)
      s = (c + e - b - d)/(a - 2*b + c)
      s = min(1._WP, max(0._WP, s))
      t = 1.-s

    case (3)
      s = 0.
      t = -e/c
      t = min(1._WP, max(0._WP, t))

    case (5)
      t = 0.
      s = -d/a
      s = min(1._WP, max(0._WP, s))
    end select

    MinDistToTri = sqrt(a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f); 
    if (present(s0)) s0 = s
    if (present(t0)) t0 = t
    if (present(x0)) x0 = (1 - s - t)*x(1, :) + s*x(2, :) + t*x(3, :)

  end function MinDistToTri

!**********************************************************************

end module ModIntOnWalls
