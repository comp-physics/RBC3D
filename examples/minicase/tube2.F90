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

#include "petsc/finclude/petsc.h"
  use petsc

  implicit none
  integer :: cutoff
  character(CHRLEN) :: fn

  call InitAll

  call TimeInt_Euler

  call FinalizeAll

  stop

contains

!**********************************************************************
  subroutine InitAll

    ! System initialization
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
    integer :: nlat0, nlatp
    real(WP) :: radEqv, radPlat
    integer :: ierr

    call ReadRestart(restart_file)

    ! Reference cells
    allocate (rbcRefs(3))

    if (nrbc > 0) then
      radEqv = 1.
      radPlat = .4
      ! assumes first rbc is a rbc
      nlat0 = rbcs(1)%nlat0
      print *, "NLAT0", nlat0
      ! nlat0 = 12
      nlatp = 4

      rbcRef => rbcRefs(1)
      call RBC_Create(rbcRef, nlat0)
      call RBC_MakeBiconcave(rbcRef, radEqv)
      call RBC_ComputeGeometry(rbcRef)

      rbcRef => rbcRefs(2)
      call RBC_Create(rbcRef, nlat0)
      call RBC_MakeSphere(rbcRef, radEqv)
      call RBC_ComputeGeometry(rbcRef)

      rbcRef => rbcRefs(3)
      call RBC_Create(rbcRef, nlatp)
      call RBC_MakePlatelet(rbcRef, radPlat)
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
      case (1)
        rbc%ES = 12.4
        rbc%ED = 200.
        rbc%EB = 6.69D-2
      case (2)
        rbc%ES = 10.
        rbc%ED = 50.
        rbc%EB = 6.D-2
      case (3)
        rbc%ES = 10.
        rbc%ED = 50.
        rbc%EB = 6.D-2
      case default
        stop "bad cellcase"
      end select
    end do ! irbc

    do iwall = 1, nwall
      wall => walls(iwall)
    end do ! iwall

    ! Background velocity
    vbkg(1:2) = 0.
    vbkg(3) = 8.
    print *, "vbkg: ", vbkg

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
