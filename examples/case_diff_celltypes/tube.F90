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

  call TimeInt

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
    type(t_Rbc),pointer :: rbc,rbcRef
    type(t_Wall),pointer :: wall
    integer :: nlat0
    real(WP) :: radEqv
    integer :: ierr

    call ReadRestart(restart_file)

    ! Reference cells
    allocate(rbcRefs(3))

    if (nrbc > 0) then
      radEqv = 1.
      nlat0 = rbcs(1)%nlat0
      
      rbcRef => rbcRefs(1)
      call RBC_Create(rbcRef, nlat0)
      call RBC_MakeBiconcave(rbcRef, radEqv)
      call RBC_ComputeGeometry(rbcRef)

      rbcRef => rbcRefs(2)
      call RBC_Create(rbcRef, nlat0)
      call RBC_MakeLeukocyte(rbcRef, radEqv)
      call RBC_ComputeGeometry(rbcRef)

      rbcRef => rbcRefs(3)
      call ImportReadRBC('Input/SickleCell.dat', rbcRef)
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
      select case(rbc%celltype)
      case(1) 
       rbc%ES = 12.4
       rbc%ED = 200.
       rbc%EB = 6.69D-2
      ! Mechanical properties of leukocytes
      ! Calculated according to section 6 of reference paper with WBC radius of 1.4 sim units
      ! and Es* scaled by 10^2
      ! Reference:
      !   Zhao, H., Isfahani, A. H., Olson, L. N., & Freund, J. B. (2010). 
      !   A spectral boundary integral method for flowing blood cells. 
      !   Journal of Computational Physics, 229(10), 3726-3744.
      case(2)
       rbc%ES = 887
       rbc%ED = 200.
       rbc%EB = 2.44D-2 ! check = .024

      ! Mechanical properties for sickle cell roughly determined to be
      ! Es = (20 / 7.1) * Es for a healthy RBC (case(1) cell)
      ! Ed = (49.4 / 15.4) * Ed for healthy RBC (case(1) cell)
      ! Eb = (19.5 / 5.7) * Eb for a healthy RBC (case(1) cell)
      ! Reference:
      !   HeeSu Byun, Timothy R. Hillman, John M. Higgins, Monica Diez-Silva, Zhangli Peng, Ming Dao, Ramachandra R. Dasari, Subra Suresh, YongKeun Park
      !   Optical measurement of biomechanical properties of individual erythrocytes
      !   Acta Biomaterialia, Volume 8, Issue 11, 2012, Pages 4130-4138,
      !   https://doi.org/10.1016/j.actbio.2012.07.011
      case(3)
       rbc%ES = 12.4 * 20 / 7.1
       rbc%ED = 200 * 49.4 / 15.4
       rbc%EB = 6.69D-2 * 19.5 / 5.7
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
     print *,vbkg

  end subroutine InitSystem

!**********************************************************************
  subroutine FinalizeSystem
    
    integer :: irbc, iwall
    type(t_Rbc),pointer :: rbc
    type(t_Wall),pointer :: wall

    do irbc = 1, nrbc
      rbc => rbcs(irbc)
      call RBC_Destroy(rbc)
    end do ! irbc

    do iwall = 1, nwall
      wall => walls(iwall)
      call Wall_Destroy(wall)
    end do ! iwall

    if (nrbc > 0) deallocate(rbcs)
    if (nwall > 0) deallocate(walls)

  end subroutine FinalizeSystem

!**********************************************************************

end program 
