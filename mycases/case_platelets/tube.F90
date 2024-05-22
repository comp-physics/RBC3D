! cells in cylindrical tubes
program cells_in_tube

  use ModDataTypes
  use ModDataStruct
  use ModRBC
  use ModWall
  use ModConf
  use ModData
  use ModTimeInt
  use ModIO
  use ModBasicMath
  use ModQuadRule

  implicit none
  integer :: cutoff
  character(CHRLEN) :: fn

#include "../../petsc_include.h"

  call InitAll

  call TimeInt_Euler

  call FinalizeAll

  stop

contains

!**********************************************************************
  subroutine InitAll

    ! System intialization
    call InitMPI()
    call GaussQuad_Init
    call ReadConfig('Input/tube.in')
    call InitSystem

    ! Note: SetEwaldPrms and DomainDecomp must be called after InitSystem
    call SetEwaldPrms
    call DomainDecomp
    call GlobData_Init

    ! Prepare for time integration
    call IO_Init
    call TimeInt_Init

  end subroutine InitAll

!**********************************************************************
  subroutine FinalizeAll

    call IO_Finalize
    call TimeInt_Finalize

    call GlobData_Finalize
    call FinalizeSystem

    call GaussQuad_Finalize
    call FinalizeMPI

  end subroutine FinalizeAll

!**********************************************************************
  subroutine InitSystem

    integer :: irbc, iwall
    type(t_Rbc), pointer :: rbc, rbcRef
    type(t_Wall), pointer :: wall
    integer :: nlat0, dealias
    real(WP) :: radEqv, radEqv2, radEqv3, platExp
    real(WP) :: shearRate, plasmaVisc, shearMod, bendingMod
    integer :: ierr

    shearMod = 4.2D-6
    bendingMod = 1.8D-19
    shearRate = 100
    plasmaVisc = 1.2D-3
    radEqv2 = 1.4
    dealias = 5

    call ReadRestart(restart_file)

    ! Reference cells
    allocate (rbcRefs(3))

    if (nrbc > 0) then
      radEqv = 1.
      radEqv2 = 1.4
      radEqv3 = .3
      platExp = 4.*3./(2.*2.82)/2.*radEqv3
      nlat0 = rbcs(1)%nlat0
      print *, "nlat0 tube", nlat0

      rbcRef => rbcRefs(1)
      call RBC_Create(rbcRef, nlat0)
      call RBC_MakeBiconcave(rbcRef, radEqv)
      call RBC_ComputeGeometry(rbcRef)

      rbcRef => rbcRefs(2)
      call RBC_Create(rbcRef, nlat0)
      call RBC_MakeLeukocyte(rbcRef, radEqv)
      call RBC_ComputeGeometry(rbcRef)

      rbcRef => rbcRefs(3)
      call RBC_Create(rbcRef, nlat0/3, dealias)
      call RBC_MakePlatelet(rbcRef, radEqv3)
      call RBC_ComputeGeometry(rbcRef)

    end if

    ! Wall periodic boundary condition
    do iwall = 1, nwall
      wall => walls(iwall)
      call Wall_Build_V2V(wall, Lb)
    end do ! iwall

    ! Mechanical properties of cells and walls
    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      select case (rbc%celltype)
      case (1) ! RBCs
        rbc%ES = shearMod / ((radEqv * 2.82 * 1D-6) * shearRate * plasmaVisc)
        ! rbc%ES = 12.4
        rbc%ED = 200.
        rbc%EB = bendingMod / ((radEqv * 2.82 * 1D-6)**3 * shearRate * plasmaVisc)
        ! rbc%EB = 6.69D-2
        if (rootWorld) then
          print *, "case 1: rbc%EB", rbc%EB, "rbc%ED", rbc%ED, "rbc%ES", rbc%ES
        end if
      case (2) ! WBCs
        rbc%ES = (shearMod * 100) / ((radEqv2 * 2.82 * 1D-6) * shearRate * plasmaVisc)
        rbc%ED = 200.
        rbc%EB = bendingMod / ((radEqv2 * 2.82 * 1D-6)**3 * shearRate * plasmaVisc)
        if (rootWorld) then
          print *, "case 2: rbc%EB", rbc%EB, "rbc%ED", rbc%ED, "rbc%ES", rbc%ES
        end if
      case (3) ! Platelets
        rbc%ES =  (shearMod / ((platExp * 2.82 * 1D-6) * shearRate * plasmaVisc)) / 5.
        ! rbc%ES   38.89 without scaling
        rbc%ED = 200. / 10.
        rbc%EB = (bendingMod / ((platExp * 2.82 * 1D-6)**3 * shearRate * plasmaVisc)) / 10.
        ! rbc%EB   2.0576 without scaling
        if (rootWorld) then
          print *, "case 3: rbc%EB", rbc%EB, "rbc%ED", rbc%ED, "rbc%ES", rbc%ES
        end if
      case default
        stop "bad cellcase"
      end select
    end do ! irbc

    do iwall = 1, nwall
      wall => walls(iwall)
    end do ! iwall

    ! Background velocity
!    if (Nt0 == 0) then
    vbkg(1:2) = 0.
    vbkg(3) = 8.
!    end if
    print *, vbkg

  end subroutine InitSystem

!**********************************************************************
  subroutine FinalizeSystem

    integer :: irbc, iwall
    type(t_Rbc), pointer :: rbc
    type(t_Wall), pointer :: wall

    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      call RBC_Destroy(rbc)
    end do ! irbc

    do iwall = 1, nwall
      wall => walls(iwall)
      call Wall_Destroy(wall)
    end do ! iwall

    if (nrbc > 0) deallocate (rbcs)
    if (nwall > 0) deallocate (walls)

  end subroutine FinalizeSystem

!**********************************************************************

end program
