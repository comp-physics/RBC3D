module ModDataTypes

  use MPI, only : MPI_DOUBLE_PRECISION

  implicit none

  ! WP -- real number precision
  ! MPI_WP -- real number precision for MPI
  integer,parameter :: WP = kind(1.D0)
  integer,parameter :: MPI_WP = MPI_DOUBLE_PRECISION

  ! CHRLEN -- default character string length
  integer,parameter :: CHRLEN = 1024

  ! Math constants
  real(WP),parameter :: HALF = 0.5
  real(WP),parameter :: THRD = 1.D0/3
  real(WP),parameter :: FOURTH = 1.D0/4
  real(WP),parameter :: PI = 3.14159265358979323846
  real(WP),parameter :: TWO_PI = 2*PI
  real(WP),parameter :: I_2PI = 1./TWO_PI
  complex(WP),parameter :: IOTA = ( 0._WP, 1._WP )

  integer,parameter :: PHYS_TO_FOUR = 1
  integer,parameter :: FOUR_TO_PHYS = -1

  public 

end module ModDataTypes
