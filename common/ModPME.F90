! Particle-mesh method for Ewald sum
module ModPME

  use MPI
  use ModDataTypes
  use ModDataStruct
  use ModConf, alpha => alpha_Ewd, Nb => Nb_Ewd, NbC => NbC_Ewd, PBspln => PBspln_Ewd
  use ModData
  use ModPFFTW
  use ModQuadRule
  use ModBasicMath

  implicit none

  ! Domain decomposition
  ! ixBgn, ixEnd -- range of local mesh index in physical space
  ! iqBgn, iqEnd -- range of local mesh index in Fourier space
  ! Note:
  !  -- The physical mesh is divided in z direction
  !  -- The Fourier mesh is divided in y direction
  integer, dimension(3) :: ixBgn, ixEnd
  integer, dimension(3) :: iqBgn, iqEnd

  ! ff -- distributed single-layer density
  ! tt -- distributed double-layer density
  ! vv -- distributed velocity
  ! ffC -- Fourier coefficients of ff
  ! ttC -- Fourier coefficients of tt
  ! vvC -- Fourier coefficnets of vv
  ! Note:
  !  -- The index of the Fourier coefficient arrays are defined as
  !       ffC(qz, qy, qx, :)
  double precision, allocatable :: ff(:, :, :, :), tt(:, :, :, :, :), vv(:, :, :, :)
  double complex, allocatable :: ffC(:, :, :, :), ttC(:, :, :, :, :), vvC(:, :, :, :)
  ! bb -- B-factor
  ! qq -- wave numbers
  real(WP), allocatable :: bb(:, :, :), qq(:, :, :, :)

  ! Indicators of non-zero single- and/or double-layer potential coefficiets
  logical :: flag_sing_lay, flag_doub_lay

  private

  public :: PME_Init, &
            PME_Finalize, &
            PME_Distrib_Source, &
            PME_Transform, &
            PME_Add_Interp_Vel

  private :: Update_Buff_Vel, &
             Distrib_Source, &
             Interp_Vel, &
             flag_Sing_Lay, flag_Doub_Lay

contains

!**********************************************************************
  subroutine PME_Distrib_Source(c1, c2, cells, walls)
    real(WP) :: c1, c2
    type(t_rbc), target, optional :: cells(:)
    type(t_wall), target, optional :: walls(:)

    type(t_SourceList), pointer :: slist
    integer :: p, ii, jj
    integer :: iwall, iele, ivert, l
    type(t_wall), pointer :: wall
    real(WP) :: xtmp(3), ftmp(3), gtmp(3), a3tmp(3), ttmp(3, 3)
    real(WP) :: xele(3, 3), fele(3, 3)

    ! Set up flags
    flag_sing_lay = (abs(c1) > 1.E-10)
    flag_doub_lay = (abs(c2) > 1.E-10)

    ! Initialize
    ff = 0.
    tt = 0.

    ! cells
    if (present(cells)) then
      slist => slist_rbc

      do p = 1, slist%nPoint
        xtmp = slist%x(p, :)

        if (flag_sing_lay) then
          ftmp = slist%f(p, :)
        else
          ftmp = 0.
        end if

        if (flag_doub_lay) then
          gtmp = slist%g(p, :)
          a3tmp = slist%a3(p, :)*slist%Bcoef(p)    !COEF
!     a3tmp = slist%a3(p,:)*(1.-slist%lam(p))    !COEF
          forall (ii=1:3, jj=1:3) ttmp(ii, jj) = gtmp(ii)*a3tmp(jj)
        else
          ttmp = 0.
        end if

        call Distrib_Source(xtmp, c1, c2, ftmp, ttmp)
      end do ! p
    end if

    ! walls
    if (present(walls)) then
      slist => slist_wall

      do p = 1, slist%npoint
        iwall = slist%indx(p, 0) - walls(1)%id + 1
        iele = slist%indx(p, 1)
        wall => walls(iwall)

        do l = 1, 3
          ivert = wall%e2v(iele, l)
          xele(l, :) = wall%x(ivert, :)
          fele(l, :) = wall%f(ivert, :)
        end do ! l

        xtmp = THRD*sum(xele, dim=1)

        if (flag_sing_lay) then
          ftmp = THRD*sum(fele, dim=1)*wall%area(iele)
        else
          ftmp = 0.
        end if

        ttmp = 0.   ! no double-layer density on wall

        call Distrib_Source(xtmp, c1, c2, ftmp, ttmp)
      end do ! p
    end if

  end subroutine PME_Distrib_Source

!**********************************************************************
! Forward FFT -> multiply weight -> backward FFT
  subroutine PME_Transform

    integer :: i, j, k, ii, jj
    real(WP) :: q(3), qt(3), q2, q2t
    real(WP) :: expq2t, phi0, phi1
    complex(WP) :: vvCtmp(3), ffCtmp(3), ttCtmp(3, 3), trtmp
    real(WP) :: vol

    ! Forward FFT
    if (flag_sing_lay) then
      do ii = 1, 3
        call Get_That(ff(:, :, :, ii), ffC(:, :, :, ii))
      end do ! ii
    end if

    if (flag_doub_lay) then
      do ii = 1, 3
      do jj = 1, 3
        call Get_That(tt(:, :, :, ii, jj), ttC(:, :, :, ii, jj))
      end do ! jj
      end do ! ii
    end if

    ! Compute vvC
    vvC = 0.

    vol = product(Lb)

    do i = iqBgn(1), iqEnd(1)
    do j = iqBgn(2), iqEnd(2)
    do k = iqBgn(3), iqEnd(3)
      if (i == 0 .and. j == 0 .and. k == 0) then
        vvC(k, j, i, :) = 0.
        cycle
      end if

      q = qq(k, j, i, :)
      qt = sqrt(pi*alpha)*q
      q2 = sum(q*q)
      q2t = pi*alpha*q2
      expq2t = exp(-q2t)

      phi0 = expq2t/q2t
      phi1 = (expq2t + phi0)/q2t

      ! Single layer potential
      if (flag_sing_lay) then
        ffCtmp = ffC(k, j, i, :)

        vvCtmp = 2*alpha/vol*phi1*(q2t*ffCtmp - qt*sum(qt*ffCTmp))
        vvC(k, j, i, :) = vvC(k, j, i, :) + vvCtmp
      end if

      ! Double layer potential
      if (flag_doub_lay) then
        ttCtmp = ttC(k, j, i, :, :)
        trtmp = ttCtmp(1, 1) + ttCtmp(2, 2) + ttCtmp(3, 3)

        vvCtmp = 0.
        vvCtmp = vvCtmp + iota*4*pi*alpha/vol*phi0*( &
                 q*trtmp + matmul(q, ttCtmp) + matmul(ttCtmp, q))
        vvCtmp = vvCtmp - iota*8*pi**2*alpha**2/vol*phi1* &
                 sum(q*matmul(ttCtmp, q))*q
        ! Note:
        !  The sign is tricky, see notes
        vvC(k, j, i, :) = vvC(k, j, i, :) - vvCtmp
      end if
    end do ! k
    end do ! j
    end do ! i

    do ii = 1, 3
      vvC(:, :, :, ii) = bb*vvC(:, :, :, ii)
    end do ! ii

    ! Backward FFT
    ! Note:
    !  -- After the transform, "ff" array stores transformed arrays
    !  -- We do not directly use "vv" array in the transform because
    !     its third dimension is bigger than the logical dimension
    do ii = 1, 3
      call Get_R(vvC(:, :, :, ii), ff(:, :, :, ii))
    end do ! ii
    vv(:, :, lbound(ff, 3):ubound(ff, 3), :) = ff

  end subroutine PME_Transform

!**********************************************************************
! Add the interpolated velocities
!  tlist -- target list
!  v -- result
  subroutine PME_Add_Interp_Vel(tlist, v)
    type(t_TargetList) :: tlist
    real(WP) :: v(:, :)

    integer :: i
    real(WP) :: dv(3)
    !print*, 'enter pme'
    ! Update buffer velocities
    call Update_Buff_Vel
!print*, '1'
!print*, 'tlist =', tlist%nPoint
    ! Add interpolated velocities
    do i = 1, tlist%nPoint
      if (tlist%active(i)) then
        call Interp_Vel(tlist%x(i, :), dv)
        v(i, :) = v(i, :) + dv/tlist%Acoef(i)  !COEF
!       v(i,:) = v(i,:) + dv/(1.+tlist%lam(i))  !COEF
      end if
    end do ! i

  end subroutine PME_Add_Interp_Vel

!**********************************************************************
! Initialize PME solver
  subroutine PME_Init

    integer :: i, j, k, ii, jj
    real(WP) :: q(3), q2, t
    ! Variables needed to compute b(:) factor
    integer :: imin, m
    complex(WP) :: b
    real(WP) :: MP(PBspln)
    integer :: ierr

    ! Initialize parallel fftw
    call Init_PFFTW(ixBgn, ixEnd, iqBgn, iqEnd)
    call Init_FFTWR
    call Init_FFTWT

    ! Allocate working arrays
    allocate (ff(ixBgn(1):ixEnd(1), ixBgn(2):ixEnd(2), ixBgn(3):ixEnd(3), 3))
    allocate (tt(ixBgn(1):ixEnd(1), ixBgn(2):ixEnd(2), ixBgn(3):ixEnd(3), 3, 3))
    ! Note:
    !  We need to extend the domain of velocity mesh in z-direction
    !  since velocity interplation needs data in the buffer zon
    allocate (vv(ixBgn(1):ixEnd(1), ixBgn(2):ixEnd(2), ixBgn(3) - PBspln:ixEnd(3), 3))

    allocate (ffC(iqBgn(3):iqEnd(3), iqBgn(2):iqEnd(2), iqBgn(1):iqEnd(1), 3))
    allocate (ttC(iqBgn(3):iqEnd(3), iqBgn(2):iqEnd(2), iqBgn(1):iqEnd(1), 3, 3))
    allocate (vvC(iqBgn(3):iqEnd(3), iqBgn(2):iqEnd(2), iqBgn(1):iqEnd(1), 3))

    allocate (bb(iqBgn(3):iqEnd(3), iqBgn(2):iqEnd(2), iqBgn(3):iqEnd(1)))
    allocate (qq(iqBgn(3):iqEnd(3), iqBgn(2):iqEnd(2), iqBgn(3):iqEnd(1), 3))

    ! Compute wave numbers
    do i = iqBgn(1), iqEnd(1)
    do j = iqBgn(2), iqEnd(2)
    do k = iqBgn(3), iqEnd(3)
      q(1) = i*iLb(1)

      if (j < Nb(2)/2) then
        q(2) = j*iLb(2)
      else
        q(2) = (j - Nb(2))*iLb(2)
      end if

      if (k < Nb(3)/2) then
        q(3) = k*iLb(3)
      else
        q(3) = (k - Nb(3))*iLb(3)
      end if

      ! Note: the Fourier coefficients computed by Adam's PFFTW code
      ! differ from those by FFTW
!      qq(k,j,i,:) = q
      qq(k, j, i, :) = (/q(1), q(2), -q(3)/)
    end do ! k
    end do ! j
    end do ! i

    ! Compute the B-factor
    ! Initialize
    bb = 1.
    if (0 >= iqBgn(2) .and. 0 <= iqEnd(2)) then
      bb(0, 0, 0) = 0.
    end if

    ! imin should be 1, so that MP(1) = M_p(1.0), MP(2) = M_p(2.0) and so on
    call BSplinefunc(PBSpln + epsilon(1._WP), PBSpln, imin, MP)

    do ii = 1, 3
      do k = iqBgn(ii), iqEnd(ii)
        b = 0.
        do m = 0, PBSpln - 2
          b = b + MP(m + 1)*exp(two_PI*iota*k*m/real(Nb(ii)))
        end do ! m
        b = exp(two_PI*iota*k*(PBSpln - 1.)/real(Nb(ii)))/b
        b = abs(b)**2

        select case (ii)
        case (1)
          bb(:, :, k) = bb(:, :, k)*b
        case (2)
          bb(:, k, :) = bb(:, k, :)*b
        case (3)
          bb(k, :, :) = bb(k, :, :)*b
        end select
      end do ! k
    end do ! ii

  end subroutine PME_Init

!**********************************************************************
! Finalize the PME solver
  subroutine PME_Finalize

    deallocate (ff, tt, vv)
    deallocate (ffC, ttC, vvC)
    deallocate (qq, bb)

    call Finalize_PFFTW

  end subroutine PME_Finalize

!**********************************************************************
! Update the buffer velocity
  subroutine Update_Buff_Vel

    integer :: numNodes, nodeNum
    integer :: sizeBuf
    integer :: nodeLeft, nodeRight
    real(WP), allocatable :: bufSend(:), bufRecv(:)
    integer :: stat(MPI_Status_Size)
    integer :: ierr

    ! Set up parameters
    call MPI_Comm_Size(MPI_COMM_Ewald, numNodes, ierr)
    call MPI_Comm_Rank(MPI_COMM_Ewald, nodeNum, ierr)

    sizeBuf = 3*Nb(1)*Nb(2)*PBspln
    nodeLeft = modulo(nodeNum - 1, numNodes)
    nodeRight = modulo(nodeNum + 1, numNodes)

!    print*, 'sizebuff =', sizebuf

    ! Allocate working arrays
    allocate (bufSend(sizeBuf), bufRecv(sizeBuf))
!print*, 'allocated'
!    ! Send data to the node to the left
!    bufSend = reshape( vv(:,:,ixBgn(3):ixBgn(3)+(PBspln-1),:), (/sizeBuf/) )
!    call MPI_SendRecv( bufSend, sizeBuf, MPI_WP, nodeLeft, nodeNum, &
!           bufRecv, sizeBuf, MPI_WP, nodeRight, nodeRight, &
!       MPI_COMM_Ewald, stat, ierr)
!    vv(:,:,ixEnd(3)+1:ixEnd(3)+PBspln,:) = &
!       reshape( bufRecv, (/Nb(1), Nb(2), PBspln, 3/) )

    ! Send data to the node to the right
    bufSend = reshape(vv(:, :, ixEnd(3) - PBspln + 1:ixEnd(3), :), (/sizeBuf/))
!    print*, 'bufsend'
    call MPI_SendRecv(bufSend, sizebuf, MPI_WP, nodeRight, nodeNum, &
                      bufRecv, sizeBuf, MPI_WP, nodeLeft, nodeLeft, &
                      MPI_COMM_Ewald, stat, ierr)
!   print*,' callled mpi'
    vv(:, :, ixBgn(3) - PBspln:ixBgn(3) - 1, :) = reshape(bufRecv, (/Nb(1), Nb(2), PBspln, 3/))

    ! Deallocate working arrays
    deallocate (bufSend, bufRecv)

  end subroutine Update_Buff_Vel

!**********************************************************************
! Distribute a source to Cartesian mesh
! Argument:
!  x(:) -- position
!  c1,c2 -- coefficients before the single- and double-layer potentials
!  f(:) -- single-layer source
!  t(:,:) -- double-layer source
  subroutine Distrib_Source(x, c1, c2, f, t)
    real(WP) :: x(3)
    real(WP) :: c1, c2, f(3), t(3, 3)

    real(WP) :: ih(3)
    real(WP) :: wx(PBspln), wy(PBspln), wz(PBspln), wxyz
    real(WP) :: ic, jc, kc
    integer :: imin, jmin, kmin, i0, j0, k0, i, j, k

    ! Calculate weights
    ih = Nb/Lb

    ic = x(1)*ih(1)
    call BSplinefunc(ic, PBspln, imin, wx)

    jc = x(2)*ih(2)
    call BSplinefunc(jc, PBspln, jmin, wy)

    kc = x(3)*ih(3)
    call BSplinefunc(kc, PBspln, kmin, wz)

    ! Do the distribution
    do k0 = 1, PBspln
      k = modulo(kmin + k0 - 1, Nb(3))
      if (k < ixBgn(3) .or. k > ixEnd(3)) cycle

      do j0 = 1, PBspln
      do i0 = 1, PBspln
        j = modulo(jmin + j0 - 1, Nb(2))
        i = modulo(imin + i0 - 1, Nb(1))

        wxyz = wx(i0)*wy(j0)*wz(k0)
        if (flag_sing_lay) ff(i, j, k, :) = ff(i, j, k, :) + (c1*wxyz)*f
        if (flag_doub_lay) tt(i, j, k, :, :) = tt(i, j, k, :, :) + (c2*wxyz)*t
      end do ! i0
      end do ! j0
    end do ! k0

  end subroutine Distrib_Source

!**********************************************************************
! Interpolate velocity from Cartesian mesh
! Arguments:
!  x -- position
!  v -- vector
  subroutine Interp_Vel(x, v)
    real(WP) :: x(3), v(3)

    real(WP) :: ih(3), wx(PBspln), wy(PBspln), wz(PBspln), wxyz
    real(WP) :: ic, jc, kc
    integer :: n, imin, jmin, kmin, i0, j0, k0, i, j, k

    ! Initialize
    v = 0.

    ! Calculate weights
    ih = Nb/Lb

    kc = x(3)*ih(3)
    call BSplinefunc(kc, PBspln, kmin, wz)

    jc = x(2)*ih(2)
    call BSplinefunc(jc, PBspln, jmin, wy)

    ic = x(1)*ih(1)
    call BSplinefunc(ic, PBspln, imin, wx)

    ! Do interpolation
    do k0 = 1, PBspln
      k = kmin + k0 - 1
      k = lbound(vv, 3) + modulo(k - lbound(vv, 3), Nb(3))
      if (k > ubound(vv, 3)) cycle

      do j0 = 1, PBspln
      do i0 = 1, PBspln
        j = modulo(jmin + j0 - 1, Nb(2))
        i = modulo(imin + i0 - 1, Nb(1))

        wxyz = wx(i0)*wy(j0)*wz(k0)
        v = v + wxyz*vv(i, j, k, :)
      end do ! i0
      end do ! j0
    end do ! k0

  end subroutine Interp_Vel

!**********************************************************************

end module ModPME
