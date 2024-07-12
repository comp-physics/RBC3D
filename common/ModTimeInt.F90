! Time integration by forward Euler scheme
module ModTimeInt

  use ModDataTypes
  use ModDataStruct
  use ModSphPk
  use ModRBC
  use ModWall
  use ModConf
  use ModData
  use ModIntOnRbcs
  use ModIntOnWalls
  use ModRbcSingInt
  use ModPME
  use ModVelSolver
  use ModNoSlip
  use ModRepulsion
  use ModIO
  use ModPostProcess

  implicit none

  private

  public :: TimeInt_Euler, &
            TimeInt_AxiSymm, &
            TimeIntModVC, &
            OneTimeInt, &
            OneTimeInt_SHB, &
            OneTimeIntModVC, &
            OneTimeInt_Verif, &
            TimeInt_Init, &
            TimeInt_Finalize, &
            Compute_Rbc_Vel, &
            Compute_Rbc_Vel_SHB, &
            TimeInt_AB2, &
            TimeInt_RK2

  private :: FilterRbcs, &
             ReboxRbcs, &
             PeriodicWall, &
             AdjustBkgVel, &
             t_rbcv

  ! To store cell velocities for higher-order time integrator method
  type t_rbcv
    real(WP), dimension(:, :, :), pointer :: v
  end type t_rbcv

contains

!**********************************************************************
! Initialize the time integrator
  subroutine TimeInt_Init

    integer :: nptRbc, nptWall
    integer :: irbc, iwall

    ! Pre-compute reference cell and wall geometries
    ! Prepare for the physical Ewald sum computation
    if (nrbc > 0) then
      call RbcPolarPatch_Create(rbcPatch, rbcs(1))

      ! call RBC_ComputeGeometry(rbcRef)  JBF:  not needed???
      rbcRefs(1)%patch => rbcPatch
      rbcRefs(2)%patch => rbcPatch
      rbcRefs(3)%patch => rbcPatch

      if (PhysEwald) then
      do irbc = 1, nrbc
        rbcs(irbc)%patch => rbcPatch
      end do ! irbc
      end if
    end if

    if (nwall > 0) then
      do iwall = 1, nwall
        call Wall_ComputeGeometry(walls(iwall))
      end do ! iwall

      ! The source list and target list only need to be created and updated once
      call SourceList_UpdateCoord(slist_wall, walls=walls)
      call TargetList_Update(tlist_wall, walls=walls)

      if (PhysEwald) then
        do iwall = 1, nwall
          call PrepareSingIntOnWall(walls(iwall))
        end do ! iwall
      end if
    end if

    ! Prepare for the Fourier Ewald sum
    if (FourierEwald) call PME_Init

  end subroutine TimeInt_Init

!**********************************************************************
! Finalize the time integrator
  subroutine TimeInt_Finalize

    if (FourierEwald) call PME_Finalize

  end subroutine TimeInt_Finalize

!**********************************************************************
! Time integrate using
! Note:
!  Time steps from Nt0+1 to Nt
  subroutine TimeInt_Euler

    integer :: lt
    integer :: irbc, iwall
    type(t_rbc), pointer :: rbc
    type(t_wall), pointer :: wall
    real(WP) :: clockBgn, clockEnd, areaExp
    integer :: ierr

    ! Time integration
    time = time0

    do lt = Nt0 + 1, Nt
      clockBgn = MPI_WTime() ! Start timing

      ! Evolve cells
      call Compute_Rbc_Vel

      ! Log area expansion of cells every 10 ts
      do irbc = 1, nrbc
        rbc => rbcs(irbc)
        if (rootWorld) then
          if (lt == 1) then
            rbc%starting_area = rbc%area
            print *, "STARTING AREA: ", rbc%area
            print *, "STARTING VOLUME: ", rbc%area, " for celltype ", rbc%celltype
          end if
          if (modulo(lt, 10) == 0) then
            areaExp = RBC_AreaExpansion(rbc)
            write (*, '(A, I3, A, F10.5, A)') &
              "area expansion of cell ", irbc, ": ", areaExp, "%"
          end if
        end if
      end do

      ! Enforce no-slip condition on the wall
      call NoSlipWall

      ! Evolve RBC
      do irbc = 1, nrbc
        rbc => rbcs(irbc)
        ! call RBC_ComputeGeometry(rbc);  print *,"UNNEEDED GEOMETRY"
        rbc%x = rbc%x + Ts*rbc%v
        ! rbc%x = rbc%x + Ts*rbc%g  ! old "NOTATION" --- pre-rigid-cell
      end do ! irbc

      ! call FilterRbcs
      call ReboxRbcs

      ! call AddR0Motion

      ! call LeukWallRepulsion

      call VolConstrainRbcs
      call InterCellRepulsion
      call FilterRbcs

      call VolConstrainRbcs
      call InterCellRepulsion
      call FilterRbcs

      call VolConstrainRbcs
      call InterCellRepulsion
      call FilterRbcs

      call VolConstrainRbcs
      call InterCellRepulsion
      call FilterRbcs

      ! Adjust background velocity
      ! call AdjustBkgVel

      ! Update time
      time = time + Ts

      clockEnd = MPI_WTime()

      ! Output results
      call WriteAll(lt, time)
      if (rootWorld) then
        write (*, '(A, I9, A, F15.5, A, F12.2)') &
          'lt = ', lt, '  T = ', time, ' time cost = ', clockEnd - clockBgn
        write (*, *)
      end if
    end do ! lt

  end subroutine TimeInt_Euler

!***********************************************************************
! Time Integrate using generic Runge-Kutta 2nd Order method
! Implementation from Chapter 5.1 of "Numerical Methods for Ordinary Differential Systems, the Initial Value Problem" by JD Lambert
  subroutine TimeInt_RK2(alpha)
    real(WP) :: alpha

    integer :: lt
    integer :: irbc, iwall
    type(t_rbc), pointer :: rbc
    type(t_wall), pointer :: wall
    real(WP) :: clockBgn, clockEnd
    integer :: ierr

    real(WP) :: b1, b2, c1, c2
    type(t_rbcv), dimension(:), allocatable :: k1 ! temporary velocities storage

    ! Set up RK2 constants from alpha
    c1 = 0
    c2 = alpha
    b2 = 1/(2*alpha)
    b1 = 1 - b2

    !set up time
    time = time0

    !allocate for k1
    allocate (k1(nrbc))
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      allocate (k1(irbc)%v, MOLD=rbc%v)
    end do

    !evaluate with RK2 generalized
    do lt = Nt0 + 1, Nt
      clockBgn = MPI_Wtime()

      !compute velocity (k1)
      call Compute_Rbc_Vel
      call NoSlipWall

      !save current velocity into k1 tmp-var
      do irbc = 1, nrbc
        rbc => rbcs(irbc)
        k1(irbc)%v = rbc%v
      end do

      !step rbcs forward so that cells are at (x_n + c2*h*U_k)
      do irbc = 1, nrbc
        rbc => rbcs(irbc)
        rbc%x = rbc%x + c2*Ts*(k1(irbc)%v)
      end do

      !calculate velocity (k2)
      call Compute_Rbc_Vel
      call NoSlipWall

      !revert cell positions back to x_n
      do irbc = 1, nrbc
        rbc => rbcs(irbc)
        rbc%x = rbc%x - c2*Ts*(k1(irbc)%v)
      end do

      !Use RK2 formula to find x_n+1
      !note that k1 is stored in the k1 tmp-var, and k2 is in each rbc%v
      do irbc = 1, nrbc
        rbc => rbcs(irbc)
        rbc%x = rbc%x + Ts*(b1*(k1(irbc)%v) + b2*rbc%v)
      end do

      call ReboxRbcs

      !update time
      time = time + Ts

      !output results
      clockEnd = MPI_Wtime()
      call WriteAll(lt, time)
      if (rootWorld) then
        write (*, '(A,I9,A,F15.5,A,F12.2)') &
          'lt = ', lt, '  T = ', time, ' time cost = ', clockEnd - clockBgn
        write (*, *)
      end if
    end do

    !deallocate k1 tmp var
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      deallocate (k1(irbc)%v)
    end do
    deallocate (k1)

  end subroutine TimeInt_RK2

!***********************************************************************
! Time Integrate using Adams-Bashforth 2nd Order
  subroutine TimeInt_AB2

    integer :: lt
    integer :: irbc, iwall
    type(t_rbc), pointer :: rbc
    type(t_wall), pointer ::wall
    real(WP) :: clockBgn, clockEnd
    integer :: ierr
    !Velocity at k-1
    type(t_rbcv), dimension(:), allocatable :: v_1

    !evaluate first timestep with Euler-Forward
    call Compute_Rbc_Vel
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      rbc%x = rbc%x + Ts*rbc%v
    end do

    !allocate space for each cell's previous velocity in v_1
    allocate (v_1(nrbc))
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      allocate (v_1(irbc)%v, MOLD=rbc%v)
    end do

    !evaluate X_2...n with AB2
    do lt = Nt0 + 2, Nt
      clockBgn = MPI_WTime() !Start timing

      !save previous velocity (U_n-1)
      do irbc = 1, nrbc
        rbc => rbcs(irbc)
        v_1(irbc)%v = rbc%v
      end do

      !compute current velocity (U_n)
      call Compute_Rbc_Vel
      call NoSlipWall

      !Evolve RBCs with AB2:
      do irbc = 1, nrbc
        rbc => rbcs(irbc)
        rbc%x = rbc%x + Ts*(3*rbc%v - v_1(irbc)%v)/2
      end do

      call ReboxRbcs

      !update time
      time = time + Ts
      clockEnd = MPI_WTime()

      !Output results
      call WriteAll(lt, time)
      if (rootWorld) then
        write (*, '(A,I9,A,F15.5,A,F12.2)') &
          'lt = ', lt, '  T = ', time, ' time cost = ', clockEnd - clockBgn
        write (*, *)
      end if
    end do

    !deallocate tmp-var
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      deallocate (v_1(irbc)%v)
    end do
    deallocate (v_1)
  end subroutine TimeInt_AB2

  subroutine TimeInt_AxiSymm

    integer :: lt
    integer :: irbc, iwall
    type(t_rbc), pointer :: rbc
    type(t_wall), pointer :: wall
    real(WP) :: clockBgn, clockEnd
    integer :: ierr

    time = time0

    do lt = Nt0 + 1, Nt
      clockBgn = MPI_WTime()

      call Compute_Rbc_Vel
      call NoSlipWall

      do irbc = 1, nrbc
        rbc => rbcs(irbc)
        rbc%x = rbc%x + Ts*rbc%v
      end do

      call ReboxRbcs
      call makeAxiSymm
      call recenter

      time = time + Ts
      clockEnd = MPI_WTime()

      call WriteAll(lt, time)
      if (rootWorld) then
        write (*, '(A,I9,A,F15.5,A,F12.2)') &
          'lt = ', lt, '  T = ', time, ' time cost = ', clockEnd - clockBgn
        write (*, *)
      end if
    end do

  end subroutine TimeInt_AxiSymm

!**********************************************************************

  subroutine TimeIntModVC

    integer :: lt
    integer :: irbc, iwall
    type(t_rbc), pointer :: rbc
    type(t_wall), pointer :: wall
    real(WP) :: clockBgn, clockEnd
    integer :: ierr, it

    ! Time integration
    time = time0

    do lt = Nt0 + 1, Nt
      clockBgn = MPI_WTime() ! Start timing

      ! Evolve cells
!      print *,"NO VEL"

      if (lt == Nt0 + 1) then
        do it = 1, 3
          call Compute_Rbc_Vel

          ! Enforce no-slip condition on the wall
!      print *,"NO NO SLIP"
          call NoSlipWall
        end do !it

      else

        call Compute_Rbc_Vel
        call NoSlipWall
      end if

!SHB    !call ModifiedVolConstrainRbcs

      ! Evolve RBC
      do irbc = 1, nrbc
        rbc => rbcs(irbc)
!        call RBC_ComputeGeometry(rbc);  print *,"UNNEEDED GEOMETRY"
        rbc%x = rbc%x + Ts*rbc%v
!   rbc%x = rbc%x + Ts*rbc%g  ! old "NOTATION" --- pre-rigid-cell
      end do ! irbc

      call FilterRbcs
      call ReboxRbcs

!      print *,"MULTIVOL"

!      call AddR0Motion

!      call LeukWallRepulsion

!      call InterCellRepulsion
!      call FilterRbcs

!      call InterCellRepulsion
!      call FilterRbcs

!      call InterCellRepulsion
!      call FilterRbcs

!      call InterCellRepulsion
!      call FilterRbcs

      ! Adjust background velocity
      call AdjustBkgVel

      ! Update time
      time = time + Ts

      clockEnd = MPI_WTime()

      ! Output results
      call WriteAll(lt, time)
      if (rootWorld) then
        write (*, '(A,I9,A,F15.5,A,F12.2)') &
          'lt = ', lt, '  T = ', time, ' time cost = ', clockEnd - clockBgn
        write (*, *)
      end if
    end do ! lt

  end subroutine TimeIntModVC

!**********************************************************************
  subroutine OneTimeInt_SHB(rbcvel)

    integer :: lt
    integer :: irbc, iwall
    type(t_rbc), pointer :: rbc
    type(t_wall), pointer :: wall
    real(WP) :: clockBgn, clockEnd
    integer :: ierr, ii
    real(WP), dimension(:, :, :, :) :: rbcvel

    time = time0
    lt = Nt0 + 1

    do ii = 1, 4
      call Compute_Rbc_Vel_SHB
      ! do we actually need this no slip call for one time step? I don't think so
      call NoSlipWall
    end do

    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      !rbc%x = rbc%x + Ts*rbc%v
      rbcvel(irbc, :, :, :) = rbc%v(:, :, :)
    end do

    !time = time + Ts

  end subroutine OneTimeInt_SHB

  subroutine OneTimeInt_Verif(modvec, xv, Ntime)

    integer :: lt, irbc, nlat0, nlat, nlon, r, c, z, ia, ib
    integer, intent(in) :: Ntime
    type(t_rbc), pointer :: rbc
    real, dimension(:, :) :: modvec, xv
    real, dimension(:, :, :), allocatable :: as, bs, myx

    modvec = 0.; xv = 0.

    time = time0
    rbc => rbcs(1)
    nlat = rbc%nlat
    nlon = rbc%nlon
    nlat0 = rbc%nlat0

    allocate (as(nlat, nlat, 3), bs(nlat, nlat, 3), myx(nlat, nlon, 3))

    do lt = 1, Ntime
      print *, 'lt =', lt

      !get spherical modes of cells
      ia = 0; ib = 0
      do irbc = 1, nrbc
        rbc => rbcs(irbc)
        do r = 1, nlat
          do c = 1, nlon
            do z = 1, 3
              ib = ib + 1
              xv(lt, ib) = rbc%x(r, c, z)
            end do
          end do
        end do
        as = 0.; bs = 0.; myx = rbc%x
        call ShAnalGau(nlat, nlon, 3, myx, size(myx, 1), size(myx, 2), as, bs, size(as, 1), size(as, 2), rbc%wshags)
        do r = 1, nlat0
          do c = r, nlat0
            do z = 1, 3
              ia = ia + 1
              modvec(lt, ia) = as(r, c, z)
            end do
            if (r .gt. 1) then
              do z = 1, 3
                ia = ia + 1
                modvec(lt, ia) = bs(r, c, z)
              end do
            end if
          end do
        end do
      end do

      call Compute_Rbc_Vel_SHB
      call NoSlipWall !new addition for test cases w/ walls

      do irbc = 1, nrbc
        rbc => rbcs(irbc)
        rbc%x = rbc%x + Ts*rbc%v
        !call WriteAll(lt,time)
      end do
      time = time + Ts
    end do

    deallocate (as, bs, myx)

  end subroutine OneTimeInt_Verif

  subroutine OneTimeInt

    !Same as  TimeInt, but for only one time step.

    integer :: lt, ii
    integer :: irbc, iwall
    type(t_rbc), pointer :: rbc
    type(t_wall), pointer :: wall
    real(WP) :: clockBgn, clockEnd
    integer :: ierr
!    real, allocatable, dimension(:,:,:) :: xsave

    ! Time integration
    time = time0

    lt = Nt0 + 1
    clockBgn = MPI_WTime() ! Start timing

!      rbc => rbcs(1)
!      allocate(xsave(rbc%nlat,rbc%nlon,3))
!      xsave = rbc%x

    do ii = 1, 6
      print *, 'ii =', ii
      call Compute_Rbc_Vel_SHB
      call NoSlipWall
    end do

!      rbc => rbcs(1)
!      print*, 'Max diff cell =', MaxVal(Abs(xsave - rbc%x))
!      deallocate(xsave)

    do irbc = 1, nrbc
      rbc => rbcs(irbc)
!        call RBC_ComputeGeometry(rbc);  print *,"UNNEEDED GEOMETRY"
      rbc%x = rbc%x + Ts*rbc%v
    end do ! irbc

!      call FilterRbcs
    call ReboxRbcs

!      print *,"MULTIVOL"

!      call AddR0Motion
!      call LeukWallRepulsion

    ! SHB commenting out the four series of calls to vol constrain and filtering
    !   call VolConstrainRbcs
    !   call InterCellRepulsion
    !   call FilterRbcs

    !   call VolConstrainRbcs
    !   call InterCellRepulsion
    !   call FilterRbcs

    !   call VolConstrainRbcs
    !   call InterCellRepulsion
    !   call FilterRbcs

    !   call VolConstrainRbcs
    !   call InterCellRepulsion
    !   call FilterRbcs
!!$
!!$  call VolConstrainRbcs
!!$  call InterCellRepulsion
!!$  call FilterRbcs

    ! Adjust background velocity
    !SHB!   call AdjustBkgVel

    ! Update time
    time = time + Ts

    clockEnd = MPI_WTime()

     !!! ! Output results
    ! call WriteAll(lt, time)
    ! if (rootWorld) then
    !   write(*, '(A,I9,A,F15.5,A,F12.2)')  &
    !           'lt = ', lt, '  T = ', time, ' time cost = ', clockEnd - clockBgn
    !   write(*, *)
    ! end if

  end subroutine OneTimeInt

!**********************************************************************

  subroutine OneTimeIntModVC

    !Same as  TimeInt, but for only one time step, and with modified volume
    !constraint effects.

    integer :: lt
    integer :: irbc, iwall
    type(t_rbc), pointer :: rbc
    type(t_wall), pointer :: wall
    real(WP) :: clockBgn, clockEnd
    integer :: ierr, i

    ! Time integration
    time = time0

    lt = Nt0 + 1
    clockBgn = MPI_WTime() ! Start timing

    ! Evolve cells
    ! print *,"NO VEL"

    !!testing for convergence here...
    do i = 1, 3
      call Compute_Rbc_Vel

      ! Enforce no-slip condition on the wall
      ! print *,"NO NO SLIP"
      call NoSlipWall

      call WriteAll(lt, time)
    end do !i

    ! Evolve RBC

    !Compute volume constraint effects
    !SHB  !call ModifiedVolConstrainRbcs

    do irbc = 1, nrbc
      rbc => rbcs(irbc)
!        call RBC_ComputeGeometry(rbc);  print *,"UNNEEDED GEOMETRY"
      rbc%x = rbc%x + Ts*rbc%v
!   rbc%x = rbc%x + Ts*rbc%g  ! old "NOTATION" --- pre-rigid-cell
    end do ! irbc

    call ReboxRbcs
    call FilterRbcs

    ! Adjust background velocity !!should this be here?
    !SHB! call AdjustBkgVel

    ! Update time
    time = time + Ts

    clockEnd = MPI_WTime()

     !!! ! Output results
    ! call WriteAll(lt, time)
    ! if (rootWorld) then
    !   write(*, '(A,I9,A,F15.5,A,F12.2)')  &
    !           'lt = ', lt, '  T = ', time, ' time cost = ', clockEnd - clockBgn
    !   write(*, *)
    ! end if

  end subroutine OneTimeIntModVc

!**********************************************************************
  subroutine Compute_Rbc_Vel

    real(WP), allocatable :: v(:, :)
    integer :: irbc, iwall, npoint, p, NN, ii
    type(t_rbc), pointer :: rbc, rbcRef
    type(t_wall), pointer :: wall
    real(WP), dimension(:, :, :), allocatable :: detJ, fdetJ, gdetJ
    integer :: ilat, ilon
    real(WP) :: vmax
    real(WP) :: vxmax, vxmin
    real(WP) :: vymax, vymin
    real(WP) :: vzmax, vzmin
    integer :: ierr

    if (nrbc == 0) return

    ! Allocate working arrays
    npoint = tlist_rbc%npoint
    allocate (v(npoint, 3))

    ! Update rbc surface geometry
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      call RBC_ComputeGeometry(rbc)
      call Rbc_BuildSurfaceSource(rbc, xFlag=.true.)
    end do ! irbc

    ! Update source and target coordinates
    call SourceList_UpdateCoord(slist_rbc, rbcs=rbcs)
    call TargetList_Update(tlist_rbc, rbcs=rbcs)

    ! Compute surface residual force
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      rbcRef => rbcRefs(rbc%celltype)
      call Shell_ResForce(rbc, rbcRef, rbc%f)
      rbc%f = -rbc%f

      ! The construction of splne_FdetJ is deferred till the
      ! intra-cellular repulsion forces are added
    end do !irbc

    ! Add intra-cell repulsion force
    call AddIntraCellRepulsionForce

!    call AddR0GravityForce

    ! Filter the surface force density and build spln_FdetJ
    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      call Rbc_SphProject(rbc, 3, rbc%f)
      call Rbc_BuildSurfaceSource(rbc, fFlag=.true.)
    end do ! irbc
    call SourceList_UpdateDensity(slist_rbc, rbcs=rbcs, updateF=.true.)

    ! Solve the RBC surface velocity and save to the local list
    call Solve_RBC_Vel(v)

    ! Save rbc velocity andn update the double-layer density sources
    p = 0
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      NN = rbc%nlat*rbc%nlon
      rbc%v = reshape(v(p + 1:p + NN, :), shape(rbc%v))  ! velocities
      if (viscRat(rbc%celltype) .gt. 0) then  ! FINITE VISCOSITY
        rbc%g = reshape(v(p + 1:p + NN, :), shape(rbc%g))
        ! for rigid, this g is the potential saved in the deflation solution
        ! recovery routines
      end if
      call Rbc_BuildSurfaceSource(rbc, gFlag=.true.)

      p = p + NN
    end do ! irbc
    call SourceList_UpdateDensity(slist_rbc, rbcs=rbcs, updateG=.true.)

    if (rootWorld) then
      vmax = maxval(abs(v))
      vxmax = maxval(v(:, 1))
      vxmin = minval(v(:, 1))
      vymax = maxval(v(:, 2))
      vymin = minval(v(:, 2))
      vzmax = maxval(v(:, 3))
      vzmin = minval(v(:, 3))
      write (*, '(A,3ES12.2)') 'cells: vmax = ', vmax
      write (*, '(A,3ES12.2)') 'cells: vxmax vxmin = ', vxmax, vxmin
      write (*, '(A,3ES12.2)') 'cells: vymax vymin = ', vymax, vymin
      write (*, '(A,3ES12.2)') 'cells: vzmax vzmin = ', vzmax, vzmin
    end if
    ! Deallocate working arrays
    deallocate (v)

  end subroutine Compute_Rbc_Vel

  subroutine Compute_Rbc_Vel_SHB

    real(WP), allocatable :: v(:, :)
    integer :: irbc, iwall, npoint, p, NN, ii
    type(t_rbc), pointer :: rbc, rbcRef
    type(t_wall), pointer :: wall
    real(WP), dimension(:, :, :), allocatable :: detJ, fdetJ, gdetJ
    integer :: ilat, ilon
    real(WP) :: vmax
    real(WP) :: vxmax, vxmin
    real(WP) :: vymax, vymin
    real(WP) :: vzmax, vzmin
    integer :: ierr

    if (nrbc == 0) return

    ! Allocate working arrays
    npoint = tlist_rbc%npoint
    allocate (v(npoint, 3))

    ! Update rbc surface geometry
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      call RBC_ComputeGeometry(rbc)
      call Rbc_BuildSurfaceSource(rbc, xFlag=.true.)
    end do ! irbc

    ! Update source and target coordinates
    call SourceList_UpdateCoord(slist_rbc, rbcs=rbcs)
    call TargetList_Update(tlist_rbc, rbcs=rbcs)

    ! Compute surface residual force
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      rbcRef => rbcRefs(rbc%celltype)
      call Shell_ResForce(rbc, rbcRef, rbc%f)
      rbc%f = -rbc%f

      ! The construction of splne_FdetJ is deferred till the
      ! intra-cellular repulsion forces are added
    end do !irbc

    ! Add intra-cell repulsion force
    ! SHB !call AddIntraCellRepulsionForce

!    call AddR0GravityForce

    ! Filter the surface force density and build spln_FdetJ
    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      call Rbc_SphProject(rbc, 3, rbc%f)
      call Rbc_BuildSurfaceSource(rbc, fFlag=.true.)
    end do ! irbc
    call SourceList_UpdateDensity(slist_rbc, rbcs=rbcs, updateF=.true.)

    ! Solve the RBC surface velocity and save to the local list
    call Solve_RBC_Vel(v)

    ! Save rbc velocity andn update the double-layer density sources
    p = 0
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      NN = rbc%nlat*rbc%nlon
      rbc%v = reshape(v(p + 1:p + NN, :), shape(rbc%v))  ! velocities
      if (viscRat(rbc%celltype) .gt. 0) then  ! FINITE VISCOSITY
        rbc%g = reshape(v(p + 1:p + NN, :), shape(rbc%g))
        ! for rigid, this g is the potential saved in the deflation solution
        ! recovery routines
      end if
      call Rbc_BuildSurfaceSource(rbc, gFlag=.true.)

      p = p + NN
    end do ! irbc
    call SourceList_UpdateDensity(slist_rbc, rbcs=rbcs, updateG=.true.)

    if (rootWorld) then
      vmax = maxval(abs(v))
      vxmax = maxval(v(:, 1))
      vxmin = minval(v(:, 1))
      vymax = maxval(v(:, 2))
      vymin = minval(v(:, 2))
      vzmax = maxval(v(:, 3))
      vzmin = minval(v(:, 3))
      write (*, '(A,3ES12.2)') 'cells: vmax = ', vmax
      write (*, '(A,3ES12.2)') 'cells: vxmax vxmin = ', vxmax, vxmin
      write (*, '(A,3ES12.2)') 'cells: vymax vymin = ', vymax, vymin
      write (*, '(A,3ES12.2)') 'cells: vzmax vzmin = ', vzmax, vzmin
    end if
    ! Deallocate working arrays
    deallocate (v)

  end subroutine Compute_Rbc_Vel_SHB

!**********************************************************************
! Filter red blood cells shapes
! Note:
!  -- It must be called every time the blood cells are moved
  subroutine FilterRbcs

    integer :: irbc
    type(t_rbc), pointer :: rbc
    integer :: nlat, nlon
    real(WP), dimension(:, :, :), allocatable :: xa, xb

    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      nlat = rbc%nlat
      nlon = rbc%nlon

      ! Allocate working arrays
      allocate (xa(nlat, nlat, 3), xb(nlat, nlat, 3))

      call ShAnalGau(nlat, nlon, 3, rbc%x, size(rbc%x, 1), size(rbc%x, 2), &
                     xa, xb, size(xa, 1), size(xa, 2), rbc%wshags)
      call ShFilter(nlat, nlon, 3, xa, xb, size(xa, 1), size(xa, 2), rbc%nlat0, rbc%nlon0)
      call ShSynthGau(nlat, nlon, 3, rbc%x, size(rbc%x, 1), size(rbc%x, 2), &
                      xa, xb, size(xa, 1), size(xa, 2), rbc%wshsgs)

      ! Deallocate working arrays
      deallocate (xa, xb)
    end do ! irbc

  end subroutine FilterRbcs

!**********************************************************************
! constrain volume of leukocytes
  subroutine VolConstrainRbcs

    integer :: irbc, ii, ilat, ilon
    type(t_rbc), pointer :: rbc
    real(WP) :: xc
    real     :: fac

    do irbc = 1, nrbc

      rbc => rbcs(irbc)

      call RBC_ComputeGeometry(rbc)
      fac = (rbcRefs(rbc%celltype)%vol - rbc%vol)/rbc%area
      fac = SIGN(MIN(ABS(fac), epsDist/20.), fac)
      if (rootWorld .and. ABS(fac) .gt. epsDist/40.) then
        print *, "VOL: ", rbc%celltype, rbcRefs(rbc%celltype)%vol, rbc%vol, fac
      end if
      do ilat = 1, rbc%nlat
        do ilon = 1, rbc%nlon
          rbc%x(ilat, ilon, :) = rbc%x(ilat, ilon, :) + fac*rbc%a3(ilat, ilon, :)
        end do
      end do

    end do ! irbc

  end subroutine VolConstrainRbcs

!**********************************************************************

  subroutine ModifiedVolConstrainRbcs

    integer :: irbc, ii, ilat, ilon
    type(t_rbc), pointer :: rbc
    real(WP) :: xc
    real     :: fac
    real(WP) :: leukVol, platVol, rbcVol

    real(WP) :: cellVol(3)

    leukVol = 4./3.*PI*(1.4)**3
    rbcVol = 4./3.*PI*(1.0)**3
    platVol = 4./3.*PI*(4.*3./(2.*2.82)/2.)**3  ! about 3 microns
    cellVol = (/rbcVol, leukVol, platVol/)

    do irbc = 1, nrbc

      rbc => rbcs(irbc)

      call RBC_ComputeGeometry(rbc)
      fac = (cellVol(rbc%celltype) - rbc%vol)/rbc%area
      fac = SIGN(MIN(ABS(fac), epsDist/20.), fac)
      if (rootWorld) then
        print *, "VOL: ", cellVol(rbc%celltype), rbc%vol, fac
      end if
      do ilat = 1, rbc%nlat
        do ilon = 1, rbc%nlon
          rbc%x(ilat, ilon, :) = rbc%x(ilat, ilon, :) + fac*rbc%a3(ilat, ilon, :)
             !!!Here's the modified part: effective velocity of VC
          rbc%v(ilat, ilon, :) = rbc%v(ilat, ilon, :) + (fac*rbc%a3(ilat, ilon, :))/Ts
        end do
      end do

    end do ! irbc

  end subroutine ModifiedVolConstrainRbcs

!**********************************************************************
! Put cells back into the box
  subroutine ReboxRbcs

    integer :: irbc, ii
    type(t_rbc), pointer :: rbc
    real(WP) :: xc

    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      do ii = 1, 3
        xc = 0.5*(maxval(rbc%x(:, :, ii)) + minval(rbc%x(:, :, ii)))
        rbc%x(:, :, ii) = rbc%x(:, :, ii) - floor(xc*iLb(ii))*Lb(ii)
      end do ! ii
    end do ! irbc

  end subroutine ReboxRbcs

!**********************************************************************
! Enforce wall periodicity
  subroutine PeriodicWall

    integer :: iwall
    type(t_wall), pointer :: wall
    integer :: ivert

    do iwall = 1, nwall
      wall => walls(iwall)
      do ivert = 1, wall%nvert
        if (wall%v2v(ivert) > 0) then
          wall%f(ivert, :) = wall%f(wall%v2v(ivert), :)
        end if
      end do ! ievert
    end do ! iwall

  end subroutine PeriodicWall

!**********************************************************************
! Adjust the background velocties to achieve the desired pressure gradient
  subroutine AdjustBkgVel

    real(WP) :: vol, pGrad(3), B
    integer :: ierr

! debug
    return
! end debug

    if (nwall == 0) return

    if (rootWorld) then
      vol = product(Lb)
!      B = 0.002*100.*(vol**2)/Lb(3)/sum(walls%areaTot)
!      B = min(B, 0.5/Ts)
      B = 100.

      pGrad = WallShearForce()/vol

      vBkg(1:2) = 0.
      vBkg(3) = vBkg(3) - B*(pGradTar(3) - pGrad(3))*Ts
    end if
    call MPI_Bcast(vBkg, 3, MPI_WP, 0, MPI_Comm_World, ierr)

  end subroutine AdjustBkgVel

!**********************************************************************
  subroutine makeAxiSymm
    integer :: irbc, i
    type(t_rbc), pointer :: rbc
    integer :: nlat, nlon
    real(WP), dimension(:, :, :), allocatable :: as, bs, myx

    do irbc = 1, nrbc
      rbc => rbcs(irbc)

      nlat = rbc%nlat
      nlon = rbc%nlon

      allocate (as(nlat, nlat, 3), bs(nlat, nlat, 3), myx(nlat, nlon, 3))
      as = 0.; bs = 0.; myx = 0.
      myx = rbc%x

      call ShAnalGau(nlat, nlon, 3, myx, size(myx, 1), size(myx, 2), as, bs, size(as, 1), size(as, 2), rbc%wshags)
      do i = 3, nlat, 2
        as(i, :, 1:2) = 0.
        bs(i, :, 1:2) = 0.
      end do
      call ShSynthGau(nlat, nlon, 3, myx, size(myx, 1), size(myx, 2), as, bs, size(as, 1), size(as, 2), rbc%wshsgs)
      if (rootWorld) print *, 'Max diff', maxval(abs(rbc%x - myx))
      rbc%x = myx

      deallocate (as, bs, myx)
    end do
  end subroutine makeAxiSymm

  subroutine recenter
    type(t_rbc), pointer :: rbc
    integer :: irbc
    real :: Xcc

    Xcc = Lb(2)/2.
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      !recenter (or get centers)
      call RBC_ComputeGeometry(rbc)
      if (rootworld) print *, 'Xc 1 =', rbc%xc(1) - Xcc, rbc%xc(2) - Xcc
      rbc%x(:, :, 1) = rbc%x(:, :, 1) - (rbc%xc(1) - Xcc)
      rbc%x(:, :, 2) = rbc%x(:, :, 2) - (rbc%xc(2) - Xcc)
      call RBC_ComputeGeometry(rbc)
      if (rootworld) print *, 'Xc 2 =', rbc%xc(1) - Xcc, rbc%xc(2) - Xcc
    end do
  end subroutine recenter

end module ModTimeInt
