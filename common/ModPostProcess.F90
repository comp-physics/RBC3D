module ModPostProcess

  use ModDataTypes
  use ModDataStruct
  use ModData
  use ModTargetList
  use ModIntOnRbcs
  use ModIntOnWalls
  use ModPME
  use ModBasicMath

  implicit none

  private

  public :: CalcVelocityField, &
            WallShearForce, &
            CellFlowRate, &
            ComputeParticleStress, &
            DistFromWall

contains

!**********************************************************************
! Compute the velocity on an abitrary point sets
! Arguments:
!  tlist -- the target list
!  v(:,:) -- velocities
  subroutine CalcVelocityField(tlist, v)
    type(t_TargetList) :: tlist
    real(WP) :: v(:, :)

    real(WP) :: c1, c2
    integer :: ii

    ! Initialize
    v = 0.

    if (PhysEwald) then
      c1 = 1.0/(4.*PI)
      c2 = 1.0/(4.*PI)  ! (1.0 - viscRat)/(8*PI)
      call AddIntOnRbcs(c1, c2, tlist, v)
      call AddIntOnWalls(c1, tlist, v)
    end if

    if (FourierEwald) then
      c1 = 1.0/(4.*PI)
      c2 = 1.0/(4.*PI)  !(1.0 - viscRat)/(8*PI)
      call PME_Distrib_Source(c1, c2, rbcs, walls)
      call PME_Transform
      call PME_Add_Interp_Vel(tlist, v)
    end if

    call TargetList_CollectArray(tlist, 3, v, MPI_Comm_World)

    ! Add the background velocity
    do ii = 1, 3
      v(:, ii) = v(:, ii) + vBkg(ii)
    end do ! ii

  end subroutine CalcVelocityField

!**********************************************************************
! Compute the total wall shear force
  function WallShearForce() result(sf)
    real(WP) :: sf(3)

    type(t_wall), pointer :: wall
    integer :: iwall, iele, ivert, l
    real(WP) :: x_ele(3, 3), f_ele(3, 3)

    sf = 0.

    do iwall = 1, nwall
      wall => walls(iwall)

      do iele = 1, wall%nele
        do l = 1, 3
          ivert = wall%e2v(iele, l)
          f_ele(l, :) = wall%f(ivert, :)
        end do ! l

        sf = sf + thrd*sum(f_ele, 1)*wall%area(iele)
      end do ! iele
    end do ! iwall

  end function WallShearForce

!**********************************************************************
! Compute the total cellular flow rate
  function CellFlowRate() result(flowrate)
    real(WP) :: flowrate

    type(t_rbc), pointer :: rbc
    integer :: irbc, ilat, ilon
    real(WP) :: zc, ds, vn

    flowrate = 0.

    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      zc = 0.5*(minval(rbc%x(:, :, 3)) + maxval(rbc%x(:, :, 3)))

      do ilon = 1, rbc%nlon
      do ilat = 1, rbc%nlat
        ds = rbc%detJ(ilat, ilon)*rbc%w(ilat)
        vn = dot_product(rbc%g(ilat, ilon, :), rbc%a3(ilat, ilon, :))
        flowrate = flowrate + (rbc%x(ilat, ilon, 3) - zc)*vn*ds
      end do ! ilat
      end do ! ilon
    end do ! irbc

    flowrate = iLb(3)*flowrate

  end function CellFlowRate

!**********************************************************************
! Compute the excessive stress due to the exisitence of particles
  subroutine ComputeParticleStress(tau)
    real(WP) :: tau(3, 3)

    type(t_rbc), pointer :: rbc
    integer :: irbc, ilat, ilon, ii, jj
    real(WP) :: dS, x(3), xc(3), f(3)

    tau = 0.

    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      do ii = 1, 3
        xc(ii) = 0.5*(maxval(rbc%x(:, :, ii)) + minval(rbc%x(:, :, ii)))
      end do ! ii

      do ilon = 1, rbc%nlon
      do ilat = 1, rbc%nlat
        x = rbc%x(ilat, ilon, :) - xc
        f = rbc%f(ilat, ilon, :)
        dS = rbc%w(ilat)*rbc%detJ(ilat, ilon)
        forall (ii=1:3, jj=1:3) tau(ii, jj) = tau(ii, jj) - dS*f(ii)*x(jj)
      end do ! ilat
      end do ! ilon
    end do ! irbc

  end subroutine ComputeParticleStress

  ! i don't like this function declaration
  function DistFromWall(type) result(minDist)
    integer :: type
    real(WP) :: minDist, cellPoint(3), wallPoint(3), currDist, clockBgn, clockEnd
    type(t_rbc), pointer :: rbc
    type(t_wall), pointer :: wall
    integer :: irbc, ilat, ilon, iwall, ivert

    minDist = huge(minDist)

    if (rootWorld) then
      clockBgn = MPI_Wtime()
    end if

    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      if (rbc%celltype .eq. type) then
        do iwall = 1, nwall
          wall => walls(iwall)

          do ivert = 1, wall%nvert
            do ilat = 1, rbc%nlat
              do ilon = 1, rbc%nlon
                currDist = VecNorm((rbc%x(ilat, ilon, :)) - (wall%x(ivert, :)))
                if (currDist .le. minDist) then
                  minDist = currDist
                  cellPoint = rbc%x(ilat, ilon, :)
                  wallPoint = wall%x(ivert, :)
                end if
              end do !ilon
            end do !ilat
          end do !ivert

        end do !iwall
      end if
    end do

    if (rootWorld) then
      clockEnd = MPI_Wtime()
      print *, "minDist: ", minDist
      print *, "cellPoint: ", cellPoint(:)
      print *, "wallPoint: ", wallPoint(:)
      print *, "total time: ", clockEnd - clockBgn
    end if
  end function DistFromWall

!**********************************************************************

end module ModPostProcess
