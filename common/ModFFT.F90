! Utilities using FFTW
! Note:
!  All arrays are assumped to be mapped from a 2*PI by 2*PI periodic domain
module ModFFT

  use ModDataTypes

  implicit none

#include "fftw3.f"

  private

  public :: FFT_Diff, &
  	FFT_Project

contains

!**********************************************************************
! Compute the derivative
! Arguments:
!  u(m,n,nvar) -- input array
!  m, n, nvar -- dimensions of u
!  u1, u2 -- derivative of u in the two dimensions
  subroutine FFT_Diff(m, n, nvar, u, u1, u2)
    integer :: m, n, nvar
    real(WP),dimension(m,n,nvar) :: u
    real(WP),dimension(m,n,nvar),optional :: u1, u2

    integer :: i, j
    double precision,allocatable :: w(:,:)
    double complex,allocatable :: wc(:,:)
    integer*8 :: fftw_fplan, fftw_bplan

    ! Compute u1
    if (present(u1)) then
      allocate(w(m,nvar), wc(m/2+1,nvar) )

      call dfftw_plan_many_dft_r2c(fftw_fplan, 1, m, nvar, &
		    w, m, 1, m, &
		    wc, m/2+1, 1, m/2+1, FFTW_ESTIMATE)
      call dfftw_plan_many_dft_c2r(fftw_bplan, 1, m, nvar, &
		    wc, m/2+1, 1, m/2+1, &
		    w, m, 1, m, FFTW_ESTIMATE)

      do j = 1, size(u,2)
	w = u(:,j,:)
	call dfftw_execute(fftw_fplan)

	do i = 1, m/2
	  wc(i,:) = (iota*(i-1)/real(m))*wc(i,:)
	end do ! j
	wc(m/2+1,:) = 0.

	call dfftw_execute(fftw_bplan)
	u1(:,j,:) = w
      end do ! j

      call dfftw_destroy_plan(fftw_fplan)
      call dfftw_destroy_plan(fftw_bplan)
      deallocate(w, wc)
    end if

    ! Compute u2
    if (present(u2)) then
      allocate(w(n,nvar), wc(n/2+1,nvar) )

      call dfftw_plan_many_dft_r2c(fftw_fplan, 1, n, nvar, &
		    w, n, 1, n, &
		    wc, n/2+1, 1, n/2+1, FFTW_ESTIMATE)
      call dfftw_plan_many_dft_c2r(fftw_bplan, 1, n, nvar, &
		    wc, n/2+1, 1, n/2+1, &
		    w, n, 1, n, FFTW_ESTIMATE)

      do i = 1, size(u,1)
	w = u(i,:,:)
	call dfftw_execute(fftw_fplan)

	do j = 1, n/2
	  wc(j,:) = (iota*(j-1)/real(n))*wc(j,:)
	end do ! j
	wc(n/2+1,:) = 0.

	call dfftw_execute(fftw_bplan)
	u2(i,:,:) = w
      end do ! i

      call dfftw_destroy_plan(fftw_fplan)
      call dfftw_destroy_plan(fftw_bplan)
      deallocate(w, wc)
    end if

  end subroutine FFT_Diff

!**********************************************************************
! Project one periodic array to another periodic one
! Arguments:
!  f -- input array
!  g -- projected array
  subroutine FFT_Project(mf, nf, nvar, f, mg, ng, g)
    integer :: mf, nf, nvar
    real(WP),dimension(mf,nf,nvar) :: f
    integer :: mg, ng
    real(WP),dimension(mg,ng,nvar) :: g

    integer :: l, m, n, i, j
    double precision,allocatable :: ftmp(:,:), gtmp(:,:)
    double complex,allocatable :: fh(:,:), gh(:,:)
    integer*8 :: fftw_fplan, fftw_bplan

    ! Trivial case
    if (mf == mg .and. nf == ng) then
      g = f
      return
    end if

    ! Allocate working arrays
    allocate(ftmp(mf,nf), gtmp(mg,ng) )
    allocate(fh(0:mf/2,0:nf-1), gh(0:mg/2,0:ng-1) )

    ! Prepare fftw
    call dfftw_plan_dft_r2c_2d(fftw_fplan, mf, nf, ftmp, fh, FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_2d(fftw_bplan, mg, ng, gh, gtmp, FFTW_ESTIMATE)

    do l = 1, nvar
      ftmp = f(:,:,l)
      call dfftw_execute(fftw_fplan)

      m = min(mf,ng)
      n = min(nf,ng)

      gh = 0.
      do i = 0, m/2-1
      do j = -n/2-1, n/2-1
        if (j >= 0) then
	  gh(i,j) = fh(i,j)
        else
	  gh(i,j+ng) = fh(i,j+nf)
	end if
      end do ! j
      end do ! i

      call dfftw_execute(fftw_bplan)

      g(:,:,l) = 1./real(mf*nf) * gtmp
    end do ! l

    ! Finalize fftw
    call dfftw_destroy_plan(fftw_fplan)
    call dfftw_destroy_plan(fftw_bplan)

    ! Deallocate working arrays
    deallocate(ftmp, gtmp)
    deallocate(fh, gh)

  end subroutine FFT_Project

!**********************************************************************

end module ModFFT
