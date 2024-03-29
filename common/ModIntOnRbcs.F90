! Numerical integration on RBCs
module ModIntOnRbcs

  use ModDataTypes
  use ModDataStruct
  use ModConf, rc => rc_Ewd
  use ModData
  use ModHashTable, only: HashTable_Index
  use ModEwaldFunc
  use ModSpline
  use ModPolarPatch
  use ModRbcSingInt

  implicit none

  private

  public :: AddIntOnRbcs

  private :: AddLinearInt

contains

!**********************************************************************
  subroutine AddIntOnRbcs(c1, c2, tlist, v)
    real(WP) :: c1, c2
    type(t_TargetList) :: tlist
    real(WP) :: v(:, :)

    type(t_SourceList), pointer :: slist
    integer :: i, j, i1, i2, i3, j1, j2, j3, p
    integer :: surfId_i, surfId_j, iRbc_i, iRbc_j
    type(t_Rbc), pointer :: rbc_i, rbc_j
    real(WP) :: xi(3), xj(3), xx(3), rr
    real(WP) :: th_i, phi_i, th_j, phi_j, dth, mask
    type(t_NbrRbcList) :: rbcList
    integer :: ilat0, ilon0
    real(WP) :: th0, phi0, x0(3), a30(3)
    real(WP) :: EA, EB, dv(3)
    real(WP) :: c2MOD
    character(*), parameter :: func_name = 'AddIntOnRbcs'
    integer :: ierr

    if (nrbc == 0) return

    slist => slist_rbc
    do i = 1, tlist%nPoint
      ! Only calculate velocity at active target points
      if (.not. tlist%active(i)) cycle

      xi = tlist%x(i, :)
      surfId_i = tlist%indx(i, 0)
      irbc_i = surfId_i - rbcs(1)%Id + 1
      if (irbc_i >= 1 .and. irbc_i <= nrbc) then
        rbc_i => rbcs(irbc_i)
        th_i = rbc_i%th(tlist%indx(i, 1))
        phi_i = rbc_i%phi(tlist%indx(i, 2))
      end if

      call NbrRbcList_Create(rbcList)

      call HashTable_Index(slist%Nc, slist%iLbNc, xi, i1, i2, i3)
      do j1 = max(i1 - 1, 0), min(i1 + 1, slist%Nc(1) + 1)
      do j2 = max(i2 - 1, 0), min(i2 + 1, slist%Nc(2) + 1)
      do j3 = max(i3 - 1, 0), min(i3 + 1, slist%Nc(3) + 1)
        j = slist%hoc(j1, j2, j3)

        do while (j > 0)
          xj = slist%x(j, :)

          xx = xj - xi
          xx = xx - nint(xx*iLb)*Lb
          rr = sqrt(sum(xx*xx))
          if (rr > rc) goto 999

          surfId_j = slist%indx(j, 0)
          irbc_j = surfId_j - rbcs(1)%Id + 1
          rbc_j => rbcs(irbc_j)
          th_j = rbc_j%th(slist%indx(j, 1))
          phi_j = rbc_j%phi(slist%indx(j, 2))

          if (surfId_i == surfId_j) then
            dth = DistOnSphere(th_i, phi_i, th_j, phi_j)
            mask = MaskFunc(dth/rbc_i%patch%radius)
          else
            mask = 0.
            ! Update the neighboring cell list
            call NbrRbcList_Insert(rbcList, slist%indx(j, :), rr)
          end if

          if (c1 /= 0) then
            call EwaldCoeff_SL(rr, EA, EB)
            v(i, :) = v(i, :) + &
                      (1.-mask)*c1/tlist%Acoef(i) &  !COEF
                      !           (1. - mask)*c1/(1.+tlist%lam(i)) &  !COEF
                      *(EA*xx*dot_product(xx, slist%f(j, :)) + EB*slist%f(j, :))
          end if

          if (c2 /= 0) then
            call EwaldCoeff_DL(rr, EA)
            v(i, :) = v(i, :) + &
                      (1.-mask)*c2*slist%Bcoef(j)/tlist%Acoef(i) &  !COEF
                      !       (1. - mask)*c2*(1.-slist%lam(j))/(1.+tlist%lam(i)) &  !COEF
                      *(EA*xx*dot_product(xx, slist%g(j, :)) &
                        *dot_product(xx, slist%a3(j, :)))
          end if

999       j = slist%next(j)
        end do ! while
      end do ! j3
      end do ! j2
      end do ! j1

      ! Singular integration
      if (irbc_i >= 1 .and. irbc_i <= nrbc) then
        c2Mod = c2*Bcoef(rbc_i%celltype)  !COEF
!         c2Mod = c2*(1.-viscRat(rbc_i%celltype))  !COEF
        call RBC_SingInt(c1, c2Mod, rbc_i, tlist%indx(i, 1), tlist%indx(i, 2), dv)
        v(i, :) = v(i, :) + dv/tlist%Acoef(i)  !COEF
!   v(i,:) = v(i,:) + dv/(1.+tlist%lam(i))  !COEF
      end if

      ! Near-singular integration
      do p = 1, rbcList%N
        irbc_j = rbcList%indx(p, 0) - rbcs(1)%Id + 1
        ilat0 = rbcList%indx(p, 1)
        ilon0 = rbcList%indx(p, 2)

        rbc_j => rbcs(irbc_j)

        ! Translate xi close to the surface
        ! And then find the projection of xi on the surface
        th0 = rbc_j%th(ilat0)
        phi0 = rbc_j%phi(ilon0)
        x0 = rbc_j%x(ilat0, ilon0, :)

        xi = tlist%x(i, :)
        xx = xi - x0; xx = xx - nint(xx*iLb)*Lb
        xi = x0 + xx
        call Spline_FindProjection(rbc_j%spln_x, xi, th0, phi0, x0)

        ! Add correction
        c2Mod = c2*Bcoef(rbc_j%celltype) !COEF
!        c2Mod = c2*(1.-viscRat(rbc_j%celltype)) !COEF
        call RBC_NearSingInt(c1, c2Mod, rbc_j, xi, x0, th0, phi0, dv)
        v(i, :) = v(i, :) + dv/tlist%Acoef(i)   !COEF
!        v(i,:) = v(i,:) + dv/(1.+tlist%lam(i))   !COEF
      end do ! l

      call NbrRbcList_Destroy(rbcList)
    end do ! i

    ! Add the contribution from the linear part of the Ewald sum
!!$    print *,"NO NEAR-SING INT"
!!$    print *,"NO LIN INT"
    call AddLinearInt(c1, c2, tlist, v)

  end subroutine AddIntOnRbcs

!**********************************************************************
! Add the linear part of the double-layer potential
  subroutine AddLinearInt(c1, c2, tlist, v)
    real(WP) :: c1, c2
    type(t_TargetList) :: tlist
    real(WP) :: v(:, :)

    type(t_RBC), pointer :: rbc
    integer :: irbc, ilat, ilon, ii, i
    real(WP) :: ds, vn, xvint(3)

    if (c2 .eq. 0) return

    ! Compute the linear part of the double layer integration
    xvint = 0.

    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      do ilon = 1, rbc%nlon
        do ilat = 1, rbc%nlat
          ds = rbc%w(ilat)*rbc%detj(ilat, ilon)
          vn = dot_product(rbc%g(ilat, ilon, :), rbc%a3(ilat, ilon, :))
          xvint = xvint + Bcoef(rbc%celltype)*rbc%x(ilat, ilon, :)*vn*ds  !COEF
          ! xvint = xvint + (1.-viscRat(rbc%celltype))*rbc%x(ilat,ilon,:)*vn*ds  !COEF
        end do ! ilat
      end do ! ilon
    end do ! irbc

    xvint = -8*PI*product(iLb)*xvint

    ! Add the linear part
    ! xvint = c2 * xvint

    do i = 1, size(v, 1)
      if (tlist%active(i)) then
        v(i, :) = v(i, :) + c2*xvint/tlist%Acoef(i) !COEF
!          v(i,:) = v(i,:) + c2*xvint/(1.+tlist%lam(i)) !COEF
      end if
    end do ! i

  end subroutine AddLinearInt

!**********************************************************************

end module ModIntOnRbcs
