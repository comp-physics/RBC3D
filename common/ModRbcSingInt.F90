module ModRbcSingInt

  use ModDataTypes
  use ModDataStruct
  use ModConf, rc=>rc_Ewd
  use ModEwaldFunc
  use ModBasicMath
  use ModPolarPatch
  use ModSpline
  use ModQuadRule

  implicit none

  public :: RBC_SingInt, &
  	RBC_NearSingInt, &
  	RBC_NearSingInt_Subtract, &
  	RBC_NearSingInt_ReAdd

  public :: NbrRbcList_Create, &
  	NbrRbcList_Destroy, &
	NbrRbcList_Insert

contains

!**********************************************************************
! Comptue the correction due to singular integration
! Arguments:
!  dv -- correction to the simple point-point sum
  subroutine RBC_SingInt(c1, c2, rbc, ilat0, ilon0, dv)
    real(WP) :: c1, c2
    type(t_Rbc) :: rbc
    integer :: ilat0, ilon0
    real(WP) :: dv(3)

    type(t_RbcPolarPatch),pointer :: patch
    integer :: p, ilat, ilon, irad, iazm
    real(WP),dimension(3) :: xi, xj, fj, gj, a3j, xx
    real(WP) :: th_j, phi_j, rr, EA, EB, mask

    patch => rbc%patch

    ! Initialize
    dv = 0.

    ! Add integrals on the local patch
    ! Note:
    !  The spline interpolation of surface coordinates does not exactly
    !  coincide with the orignal ones on mesh points. 
    !
    !  Since the coordinates of quadrature poins are by interpolation, we
    !  should also use interpolation for the target point position, even
    !  though the target point.
    !
    !  Fail to do that would cause a small but non-zero deviation of the
    !  target point from the center of the polar coordinate patch, and 
    !  results in divergence of the integration with increasing number 
    !  of qudrature points in radial direction.
    call Spline_Interp(rbc%spln_x, rbc%th(ilat0), rbc%phi(ilon0), xi)

    do irad = 1, patch%nrad
    do iazm = 1, patch%nazm
      th_j = patch%thG(irad,iazm,ilat0,ilon0)
      phi_j = patch%phiG(irad,iazm,ilat0,ilon0)
      call Spline_Interp(rbc%spln_x, th_j, phi_j, xj)

      ! xi is the target point
      xx = xj - xi
      rr = sqrt(sum(xx*xx))
      if (rr >= rc) cycle

      if (c1 .ne. 0) then
	call Spline_Interp(rbc%spln_FdetJ, th_j, phi_j, fj)
	fj = patch%w(irad)*fj

	call EwaldCoeff_SL(rr, EA, EB)
	dv = dv + c1*(EA*xx*dot_product(xx,fj) + EB*fj)
      end if

      if (c2 .ne. 0) then
	call Spline_Interp(rbc%spln_GdetJ, th_j, phi_j, gj)
	gj = patch%w(irad)*gj
	call Spline_Interp(rbc%spln_a3, th_j, phi_j, a3j)

	call EwaldCoeff_DL(rr, EA)
	dv = dv + c2*(EA*xx*dot_product(xx,gj)*dot_product(xx,a3j))
      end if
    end do ! iazm
    end do ! irad

  end subroutine RBC_SingInt


!**********************************************************************
! Compute the correction due to near singular integration
! Arguments:
!  c1, c2, rbc
!  xi -- target point
!  x0, th0, phi0 -- projection of xi on the surface
!  dv -- correction to the simple point-point sums
! 
! Note:
!    -- The target point xi must have been translated to the surface 
!       as close as possible
  subroutine RBC_NearSingInt(C1, C2, rbc, xi, x0, th0, phi0, dv)
    real(WP) :: c1, c2
    type(t_Rbc) :: rbc
    real(WP) :: xi(3), x0(3), th0, phi0
    real(WP) :: dv(3)

    type(t_RbcPolarPatch),pointer :: patch
    real(WP) :: sizePat, dist, dist1, dist2
    real(WP) :: a30(3), detJ0(1), g0(3)
    real(WP) :: xi0(3), xi1(3)
    real(WP) :: dvtmp(3), dv0(3), dv1(3)

    patch => rbc%patch

    ! Find distance to the surface
    call Spline_Interp(rbc%spln_a3, th0, phi0, a30)
    dist = dot_product(a30, xi-x0)

    ! Distance check
    if (dist > 2*rbc%meshSize) then
      dv = 0.
      return
    end if

    ! Set up patch size 
    sizePat = patch%radius*sqrt(rbc%area/(4*PI))
    dist1 = sign(0.01*sizePat, dist)

    ! Initialize 
    dv = 0.

    ! Subtract
    call RBC_NearSingInt_Subtract(c1, c2, rbc, xi, x0, th0, phi0, patch%radius, dvtmp)
    dv = dv + dvtmp

    ! Re-add
    if (abs(dist) >= abs(dist1)) then
      call RBC_NearSingInt_ReAdd(c1, c2, rbc, xi, x0, th0, phi0, patch%radius, dvtmp)
      dv = dv + dvtmp
    else
      xi1 = x0 + dist1*a30
      call RBC_NearSingInt_ReAdd(c1, c2, rbc, xi1, x0, th0, phi0, patch%radius, dv1)

      xi0 = x0
      call RBC_NearSingInt_ReAdd(c1, c2, rbc, xi0, x0, th0, phi0, patch%radius, dv0)

      ! Take care of the jump conditions when there is double-layer integral
      if (c2 .ne. 0) then
	call Spline_Interp(rbc%spln_detJ, th0, phi0, detJ0)
	call Spline_Interp(rbc%spln_GdetJ, th0, phi0, g0)
	g0 = g0/detJ0(1)

	if (dist > 0) then
	    dv0 = dv0 + c2*4*PI*g0
        else
	    dv0 = dv0 - c2*4*PI*g0
	end if
      end if

      dvtmp = dv0 + dist/dist1*(dv1 - dv0);
      dv = dv + dvtmp
    end if

  end subroutine RBC_NearSingInt


!**********************************************************************
! Compute the correction because of subtracting the masked
! surface integral near the projection point
! Arguments:
!  dv -- the correction
! Note:
!  -- dv on exit already contains the "-" sign, so it should be added to the original integral
  subroutine RBC_NearSingInt_Subtract(C1, C2, rbc, xi, x0, th0, phi0, radPat, dv)
    real(WP) :: c1, c2
    type(t_Rbc) :: rbc
    real(WP) :: xi(3), x0(3), th0, phi0, radPat
    real(WP) :: dv(3)

    integer :: ilat, ilon, p
    integer :: nptPat
    integer,allocatable,dimension(:,:) :: ijsPat
    real(WP),dimension(3) :: xj, fj, gj, a3j, xx
    real(WP) :: th_j, phi_j, dth, mask, rr, EA, EB

    ! Initialize
    dv = 0.

    ! Substract the imporperly added contribution from the surface
    allocate(ijsPat(rbc%nlat*rbc%nlon,2))
    call PolarPatch_FindPoints(th0, phi0, radPat, rbc%th, rbc%phi, nptPat, ijsPat)

    ! Subtract the integration improperly added
    do p = 1, nptPat
      ilat = ijsPat(p,1)
      ilon = ijsPat(p,2)

      xj = rbc%x(ilat,ilon,:)

      ! xi is the target point
      xx = xj - xi
      rr = sqrt(sum(xx**2))
      if (rr > rc) cycle

      dth = DistOnSphere(th0, phi0, rbc%th(ilat), rbc%phi(ilon))
      mask = MaskFunc(dth/radPat)

      if (C1 /= 0) then
	fj = rbc%f(ilat,ilon,:)*rbc%detJ(ilat,ilon)*rbc%w(ilat)
	call EwaldCoeff_SL(rr, EA, EB)
	dv = dv - c1*mask*(EA*xx*dot_product(xx,fj) + EB*fj)
      end if

      if (C2 /= 0) then
	gj = rbc%g(ilat,ilon,:)*rbc%detJ(ilat,ilon)*rbc%w(ilat)
	a3j = rbc%a3(ilat,ilon,:)
	call EwaldCoeff_DL(rr, EA)
	dv = dv - c2*mask*(EA*xx*dot_product(xx,gj)*dot_product(xx,a3j))
      end if
    end do ! p

    ! Deallocate working arrays
    deallocate(ijsPat)

  end subroutine RBC_NearSingInt_Subtract

!**********************************************************************
! Add the integral on the local patch around the projection point
! Arguments:
!  dv -- the correction term
  subroutine RBC_NearSingInt_ReAdd(C1, C2, rbc, xi, x0, th0, phi0, radPat, dv)
    real(WP) :: c1, c2
    type(t_Rbc) :: rbc
    real(WP) :: xi(3), x0(3), th0, phi0, radPat
    real(WP) :: dv(3)

    ! dist -- distance to the surface
    ! sizePat -- (approximate) physical radius of the patch
    real(WP) :: dist, sizePat
    integer :: nrad, nazm, irad, iazm
    real(WP),allocatable,dimension(:) :: thPat, phiPat, wtPat
    real(WP),allocatable,dimension(:,:) :: thPatG, phiPatG
    real(WP),dimension(3) :: xj, fj, gj, a3j, xx, dvi
    real(WP) :: th_j, phi_j, rr, EA, EB

    ! Calculate length scales
    dist = sqrt( dot_product(xi-x0, xi-x0) )
    sizePat = radPat*sqrt(rbc%area/(4*PI))

    ! Initialize
    dv = 0.

    ! Build the patch locally
    nrad = 16
    nazm = 2*nrad
    allocate(thPat(nrad), phiPat(nazm), wtPat(nrad) )
    allocate(thPatG(nrad,nazm), phiPatG(nrad,nazm) )

    if (dist > tiny(dist)) then
      ! Scale the distance to the reference unit sphere
      dist = dist*(radPat/sizePat)    
      call GauLeg_Sinh(0._WP, radPat, 0._WP, dist, nrad, thPat, wtPat)
    else
      call GauLeg(0._WP, radPat, nrad, thPat, wtPat)
    end if

    do irad = 1, nrad
      wtPat(irad) = wtPat(irad) * sin(thPat(irad)) * (two_PI/nazm)
      wtPat(irad) = wtPat(irad) * MaskFunc( thPat(irad)/radPat)
    end do ! irad
    phiPat = (/ ((iazm-1)*two_pi/nazm, iazm=1,nazm) /)

    ! Map the patch to the global coordinate
    call PolarPatch_Build(th0, phi0, thPat, phiPat, thPatG, phiPatG)

    ! Add the contribution from patch
    do irad = 1, nrad
    do iazm = 1, nazm
      th_j = thPatG(irad,iazm)
      phi_j = phiPatG(irad,iazm)
      call Spline_Interp(rbc%spln_x, th_j, phi_j, xj)

      ! xi is the target point in BIE
      xx = xj - xi
      rr = sqrt(sum(xx**2))
      if (rr > rc) cycle

      if (C1 .ne. 0) then
	call Spline_Interp(rbc%spln_FdetJ, th_j, phi_j, fj)
	fj = wtPat(irad)*fj

	call EwaldCoeff_SL(rr, EA, EB)
	dv = dv + c1*(EA*xx*dot_product(xx, fj) + EB*fj)
      end if

      if (C2 .ne. 0) then
	call Spline_Interp(rbc%spln_a3, th_j, phi_j, a3j)
	call Spline_Interp(rbc%spln_GdetJ, th_j, phi_j, gj)
	gj = wtPat(irad)*gj

	call EwaldCoeff_DL(rr, EA)
	dv = dv + c2*(EA*xx*dot_product(xx,gj)*dot_product(xx,a3j))
      end if
    end do ! iazm
    end do ! irad

    deallocate(thPat, phiPat, wtPat, thPatG, phiPatG)

  end subroutine RBC_NearSingInt_ReAdd



!**********************************************************************
! RBC neighbor list subroutines

!**********************************************************************
  subroutine NbrRbcList_Create(list)
    type(t_NbrRbcList) :: list

    integer,parameter :: NMax = 32 !32 !SHB edit

    allocate(list%indx(NMax,0:2), list%dist(NMax))
    list%N = 0

  end subroutine NbrRbcList_Create

!**********************************************************************
  subroutine NbrRbcList_Destroy(list)
    type(t_NbrRbcList) :: list

    deallocate(list%indx, list%dist)
    list%N = 0

  end subroutine NbrRbcList_Destroy

!**********************************************************************
! Add a element to the list
! Arguments:
!  list --
!  indx --
!  dist --
  subroutine NbrRbcList_Insert(list, indx, dist)
    type(t_NbrRbcList) :: list
    integer :: indx(0:2)
    real(WP) :: dist

    integer :: i

    do i = 1, list%N
      if (list%indx(i,0) == indx(0)) then
	if (dist < list%dist(i)) then
	  list%indx(i,1:2) = indx(1:2)
	  list%dist(i) = dist
	end if
	return
      end if
    end do ! i

    ! Expand the list if needed
    list%N = list%N + 1
    i = list%N
    list%indx(i,0:2) = indx
    list%dist(i) = dist

  end subroutine NbrRbcList_Insert

!**********************************************************************

end module ModRbcSingInt
