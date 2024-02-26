! Basic functions needed for evaluating Ewald sum
module ModEwaldFunc

  use ModDataTypes
  use ModConf, alpha => alpha_Ewd, rc => rc_Ewd

  implicit none

  integer, parameter :: Ntable = 8192

  private :: Ntable

  public :: EwaldCoeff_SL_Exact, &
            EwaldCoeff_DL_Exact, &
            EwaldCoeff_SL, &
            EwaldCoeff_DL

contains

!**********************************************************************
! The Ewlad coefficient of single layer potential
! Argument:
!  r -- separation distance
!  A, B -- the velocity induced is A*x*dot_product(x,f) + B*f,
  subroutine EwaldCoeff_SL_Exact(r, alpha, A, B)
    real(WP) :: r, alpha, A, B

    real(WP) :: r_t, c1, c2, ir, ir2

    if (alpha <= 0) then
      A = 1/(r*r*r)
      B = 1/r
      return
    end if

    r_t = sqrt(PI/alpha)*r

    if (r_t < 1.D-3) then
      A = 0.
      B = 0.
    else
      c1 = erfc(r_t)
      c2 = 2./sqrt(alpha)*exp(-r_t**2)

      ir = 1./r
      ir2 = ir*ir

      A = c1*ir*ir2 + c2*ir2
      B = c1*ir - c2
    end if

  end subroutine EwaldCoeff_SL_Exact

!**********************************************************************
! The exact Ewald sum coefficient for double layer potential
! Arguments:
!  r -- separation distance
!  A -- the resulting stress tensor is A*x*x*(x*f)/r^5
  subroutine EwaldCoeff_DL_Exact(r, alpha, A)
    real(WP) :: r, alpha, A

    real(WP) :: r_t

    if (alpha <= 0) then
      A = -6/(r**5)
      return
    end if

    r_t = sqrt(PI/alpha)*r

    if (r_t < 1.D-3) then
      A = 0.
    else
      A = exp(-r_t**2)*(1.5*r_t + r_t**3) + 0.75*sqrt(PI)*erfc(r_t)
      A = -8/sqrt(PI)*A    ! Now we get A*(r**5)
      A = A/(r**5)
    end if

  end subroutine EwaldCoeff_DL_Exact

!**********************************************************************
! Ewlad coefficient of single layer potential
! Argument:
!  r -- separation distance
!  A, B -- the velocity induced is A*x*dot_product(x,f) + B*f
  subroutine EwaldCoeff_SL(r, A, B)
    real(WP) :: r, A, B

    integer, parameter :: N = Ntable
    logical, save :: table_inited = .false.
    real(WP), save :: c1_tab(0:N), c2_tab(0:N)
    real(WP), save :: r_eps
    integer :: i
    real(WP) :: r_t, s, c1, c2, ir, ir2

    ! Build look up table
    if (.not. table_inited) then
      do i = 0, N
        r_t = sqrt(PI/alpha)*(i*rc/N)

        c1_tab(i) = erfc(r_t)
        c2_tab(i) = 2/sqrt(alpha)*exp(-r_t**2)
      end do ! i

      r_eps = 1.D-3*sqrt(alpha/PI)
      table_inited = .true.
    end if

    if (r < r_eps) then
      A = 0.
      B = 0.
    else
      s = N*r/rc
      i = floor(s)

      if (i >= N) then
        A = 0.
        B = 0.
      else
        c1 = c1_tab(i)*(i + 1 - s) + c1_tab(i + 1)*(s - i)
        c2 = c2_tab(i)*(i + 1 - s) + c2_tab(i + 1)*(s - i)

        ir = 1./r
        ir2 = ir*ir

        A = c1*ir*ir2 + c2*ir2
        B = c1*ir - c2
      end if
    end if

  end subroutine EwaldCoeff_SL

!**********************************************************************
! Ewald sum coefficient for double layer potential
! Argument:
!  r -- separation distance
!  A -- the resulting stress tensor is A*x*x*(x*f)
! Note:
!  Use lookup table, whose index is linear about r
!  and the value is A*(r**5)
  subroutine EwaldCoeff_DL(r, A)
    real(WP) :: r, A

    real(WP), save :: r_eps
    integer, parameter :: N = Ntable
    logical, save :: table_inited = .false.
    real(WP), save :: c1_tab(0:N)
    integer :: i
    real(WP) :: r_t, s, c1

    ! Build look up table
    if (.not. table_inited) then
      do i = 0, N
        r_t = sqrt(PI/alpha)*(i*rc/N)

        ! Now we get A*(r**5)
        c1_tab(i) = -8/sqrt(PI)*(exp(-r_t**2)*(1.5*r_t + r_t**3) + 0.75*sqrt(pi)*erfc(r_t))
      end do ! i

      r_eps = 1.D-3*sqrt(alpha/PI)
      table_inited = .true.
    end if

    if (r < r_eps) then
      A = 0.
    else
      s = N*r/rc
      i = floor(s)

      if (i >= N) then
        A = 0.
      else
        c1 = c1_tab(i)*(i + 1 - s) + c1_tab(i + 1)*(s - i)
        A = c1/(r**5)
      end if
    end if

  end subroutine EwaldCoeff_DL

!**********************************************************************

end module ModEwaldFunc
