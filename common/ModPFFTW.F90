! Parallel fftw by Adam Willis
! Edited by Hong Zhao
module ModPFFTW

  use MPI
  use ModDataTypes
  use ModConf, GN => Nb_Ewd, MPI_Comm_PFFTW => MPI_Comm_Ewald

  implicit none

#include "fftw3.f"

  integer :: NN     ! number of nodes
  integer :: my_rank

  integer, allocatable, dimension(:)  :: chunkR
  integer, allocatable, dimension(:, :)    :: pointsR
  integer, allocatable, dimension(:)  :: chunkT
  integer, allocatable, dimension(:, :)    :: pointsT

  double precision, allocatable, dimension(:, :, :) :: rhoR
  double complex, allocatable, dimension(:, :, :)   :: rhoRh, rhoT, rhoTh
  integer(8)    :: plan_R_Rh, plan_Rh_R
  integer(8)    :: plan_T_Th, plan_Th_T

  private

  public :: Init_PFFTW, &
            Init_FFTWT, &
            Init_FFTWR, &
            Get_That, &
            Get_R, &
            Finalize_PFFTW

  private :: getRdata_a2a, &
             getTdata_a2a, &
             get_rhoR, &
             get_rhoRh, &
             get_rhoT, &
             get_rhoTh

CONTAINS

!**********************************************************************
  SUBROUTINE Init_PFFTW(ixBgn, ixEnd, iqBgn, iqEnd)
    integer, dimension(3) :: ixBgn, ixEnd, iqBgn, iqEnd

    integer :: i, j, k
    integer :: ierr

    call MPI_Comm_Size(MPI_Comm_PFFTW, NN, ierr)
    call MPI_Comm_Rank(MPI_Comm_PFFTW, my_rank, ierr)

    allocate (chunkR(0:NN - 1), chunkT(0:NN - 1), pointsR(0:NN - 1, 2), pointsT(0:NN - 1, 2))

    chunkR(:) = INT(Gn(3)/NN)
    chunkT(:) = INT(Gn(2)/NN)

    if (MOD(Gn(3), NN) .ne. 0) then
      do i = 1, MOD(Gn(3), NN)
        chunkR(i) = chunkR(i) + 1
      end do
    end if

    if (MOD(Gn(2), NN) .ne. 0) then
      do i = 1, MOD(Gn(2), NN)
        chunkT(i) = chunkT(i) + 1
      end do
    end if

    pointsR(0, 1) = 1  !indicies of R on processor
    pointsR(0, 2) = chunkR(0)
    pointsT(0, 1) = 1
    pointsT(0, 2) = chunkT(0)
    do i = 1, NN - 1
      pointsR(i, 1) = pointsR(i - 1, 2) + 1
      pointsR(i, 2) = pointsR(i, 1) + chunkR(i) - 1
    end do

    do i = 1, NN - 1
      pointsT(i, 1) = pointsT(i - 1, 2) + 1
      pointsT(i, 2) = pointsT(i, 1) + chunkT(i) - 1
    end do

    ixBgn = (/0, 0, pointsR(my_rank, 1) - 1/)
    ixEnd = (/GN(1) - 1, GN(2) - 1, pointsR(my_rank, 2) - 1/)

    iqBgn = (/0, pointsT(my_rank, 1) - 1, 0/)
    iqEnd = (/GN(1)/2, pointsT(my_rank, 2) - 1, GN(3) - 1/)

  END SUBROUTINE Init_PFFTW

!**********************************************************************
  SUBROUTINE Init_FFTWT
    integer                     :: nembed, onembed
    integer                     :: idist, odist, howmany, rank, n, istride, ostride

    allocate (rhoT(Gn(3), chunkT(my_rank), Gn(1)/2 + 1), rhoTh(Gn(3), chunkT(my_rank), Gn(1)/2 + 1))

    rank = 1
    nembed = Gn(3)
    onembed = Gn(3)
    n = Gn(3)
    howmany = chunkT(my_rank)*(Gn(1)/2 + 1)
    idist = Gn(3)
    odist = Gn(3)
    istride = 1
    ostride = 1

    call dfftw_plan_many_dft(plan_T_Th, rank, n, howmany, &
                             rhoT, nembed, istride, idist, rhoTh, onembed, ostride, odist, +1, FFTW_ESTIMATE)

    call dfftw_plan_many_dft(plan_Th_T, rank, n, howmany, &
                             rhoTh, onembed, ostride, odist, rhoT, nembed, istride, idist, -1, FFTW_ESTIMATE)

  END SUBROUTINE Init_FFTWT

!**********************************************************************
  SUBROUTINE Init_FFTWR
    integer, dimension(2)    :: nembed, onembed
    integer :: idist, odist, howmany, rank, n, istride, ostride
    integer, dimension(2):: nList

    allocate (rhoR(Gn(1), Gn(2), chunkR(my_rank)), rhoRh(Gn(1)/2 + 1, Gn(2), chunkR(my_rank)))

    nembed(1) = Gn(1)           !x-dimension in
    nembed(2) = Gn(2)           !y-dimension in
    onembed(1) = Gn(1)/2 + 1     !x_dimension out
    onembed(2) = Gn(2)          !y-dimension out
    nList(1) = Gn(1)
    nList(2) = Gn(2)                !size of data for
    rank = 2                    !number of dimensions of transform
    howmany = chunkR(my_rank)   !number of transforms of rank and size n
    idist = Gn(1)*Gn(2)         !how far to move each for each in data
    odist = (Gn(1)/2 + 1)*Gn(2)  !how far to move for each out data
    istride = 1
    ostride = 1

    call dfftw_plan_many_dft_r2c(plan_R_Rh, rank, nList, howmany, &
                                 rhoR, nembed, istride, idist, rhoRh, onembed, ostride, odist, FFTW_ESTIMATE)

    call dfftw_plan_many_dft_c2r(plan_Rh_R, rank, nList, howmany, &
                                 rhoRh, onembed, ostride, odist, rhoR, nembed, istride, idist, FFTW_ESTIMATE)

  END SUBROUTINE Init_FFTWR

!**********************************************************************
  SUBROUTINE Get_That(rhoR_in, rhoTh_out)
    double precision, dimension(:, :, :), intent(in) :: rhoR_in
    double complex, dimension(:, :, :), intent(out) :: rhoTh_out
    character(20)                           :: fn
    integer                                 :: i, j, k

    !integer(8)                              :: ll_plan
    !double complex,dimension(Gn(1)/2+1,Gn(2))      :: buffer
    !double complex,dimension(Gn(3))                :: buffer_z

    rhoR = rhoR_in
    call get_rhoRh   !2-d fourier transform over different slabs in z

    call getTdata_a2a(rhoRh, rhoT)

    call get_rhoTh      !1-d fouier transforms in z over all x,y of process

    rhoTh_out = rhoTh

  END SUBROUTINE Get_That

!**********************************************************************
  SUBROUTINE Get_R(rhoTh_in, rhoR_out)
    double complex, dimension(:, :, :), intent(in) :: rhoTh_in
    double precision, dimension(:, :, :), intent(out) :: rhoR_out

    rhoTh = rhoTh_in

    call get_rhoT   !1-d inverse fft across z in all x,y

    call getRdata_a2a(rhoT, rhoRh)   !get full x,y data

    call get_rhoR   ! 2d inverse fouier transform over slabs in z

    ! Modified by Hong Zhao to remove the normalization
    rhoR_out = rhoR

  END SUBROUTINE Get_R

!**********************************************************************
  SUBROUTINE getRdata_a2a(dataT, dataR)
    !same at getRdata, but it uses mpi_alltoall...better be more optimized
    double complex, dimension(Gn(1)/2 + 1, Gn(2), chunkR(my_rank)), intent(out)   :: dataR
    double complex, dimension(Gn(3), chunkT(my_rank), Gn(1)/2 + 1), intent(in)    :: dataT

    character(20) :: fn
    integer                                                 :: i, j, k, l, m, xi, yi, zi, point

    integer, dimension(0:NN - 1)             :: sdispl, rdispl, send_counts, recv_counts

    double complex, allocatable :: send_buff(:, :, :), recv_buff2(:), recv_buff3(:)
    integer :: ierr

    ! Edit by Hong Zhao
    allocate (send_buff(Gn(1)/2 + 1, chunkT(my_rank), Gn(3)))
    ! End edit by Hong Zhao

    dataR = 0.

    do i = 0, NN - 1
      send_counts(i) = (Gn(1)/2 + 1)*chunkR(i)*chunkT(my_rank)
      recv_counts(i) = (Gn(1)/2 + 1)*chunkR(my_rank)*chunkT(i)
    end do

    sdispl(0) = 0
    rdispl(0) = 0
    do i = 1, NN - 1
      sdispl(i) = sdispl(i - 1) + send_counts(i - 1)
      rdispl(i) = rdispl(i - 1) + recv_counts(i - 1)
    end do

    !rearrange for sending...distributing  data
    !(z,y,x) --> (x,y,z)

    do zi = 1, Gn(3)
    do yi = 1, chunkT(my_rank)
      send_buff(:, yi, zi) = dataT(zi, yi, :)
    end do
    end do

    allocate (recv_buff3(2*sum(recv_counts(:))))

    recv_buff3 = 0.

    call MPI_alltoallv(send_buff, send_counts, sdispl, MPI_DOUBLE_COMPLEX, &
                       recv_buff3, recv_counts, rdispl, MPI_DOUBLE_COMPLEX, MPI_Comm_PFFTW, ierr)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!takes input data and loops over it to create the Transpose data
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    do i = 0, NN - 1  !loops over the data received from other processors
    do xi = 1, Gn(1)/2 + 1   !loop over x locations
      do zi = 1, ChunkR(my_rank)  !loop over z locations
      do yi = 1, ChunkT(i)  !loop over y locations
        point = rdispl(i) + (Gn(1)/2 + 1)*(yi - 1) + (Gn(1)/2 + 1)*chunkT(i)*(zi - 1) + xi
        dataR(xi, pointsT(i, 1) - 1 + yi, zi) = recv_buff3(point)
      end do
      end do
    end do
    end do

    ! Edit by Hong Zhao
    deallocate (send_buff)
    ! End edit by Hong Zhao

  END SUBROUTINE getRdata_a2a

!**********************************************************************
  SUBROUTINE getTdata_a2a(dataR, dataT)
    !same at getTdata, but it uses mpi_alltoall...better be more optimized
    double complex, dimension(Gn(1)/2 + 1, Gn(2), chunkR(my_rank)), intent(in)    :: dataR
    double complex, dimension(Gn(3), chunkT(my_rank), Gn(1)/2 + 1), intent(out)   :: dataT

    character(20)                                           :: fn
    integer                                                 :: i, j, k, l, m, xi, yi, zi, point

    integer, dimension(0:NN - 1)             :: sdispl, rdispl, send_counts, recv_counts

    double complex, allocatable :: send_buff(:, :, :), recv_buff2(:), recv_buff3(:)
    integer :: ierr

    ! Edit by Hong Zhao
    allocate (send_buff(Gn(1)/2 + 1, chunkR(my_Rank), Gn(2)))
    ! End edit by Hong Zhao

    dataT = 0.

    do i = 0, NN - 1
      send_counts(i) = (Gn(1)/2 + 1)*chunkR(my_rank)*chunkT(i)
      recv_counts(i) = (Gn(1)/2 + 1)*chunkR(i)*chunkT(my_rank)
    end do

    sdispl(0) = 0
    rdispl(0) = 0
    do i = 1, NN - 1
      sdispl(i) = sdispl(i - 1) + send_counts(i - 1)
      rdispl(i) = rdispl(i - 1) + recv_counts(i - 1)
    end do

    !rearrange for sending...distributing y data
    !(x,y,z) -> (x,z,y)
    do j = 1, Gn(2)
    do k = 1, chunkR(my_rank)
      send_buff(:, k, j) = dataR(:, j, k)
    end do
    end do

    allocate (recv_buff3(2*sum(recv_counts(:))))
    recv_buff3 = 0.

    call MPI_alltoallv(send_buff, send_counts, sdispl, MPI_DOUBLE_COMPLEX, &
                       recv_buff3, recv_counts, rdispl, MPI_DOUBLE_COMPLEX, MPI_Comm_PFFTW, ierr)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !takes input data and loops over it to create the Transpose data
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 0, NN - 1  !loops over the data received from other processors
    do xi = 1, Gn(1)/2 + 1   !loop over x locations
      do zi = 1, ChunkR(i)  !loop over z locations
      do yi = 1, ChunkT(my_rank)  !loop over y locations
        point = rdispl(i) + (Gn(1)/2 + 1)*(zi - 1) + (Gn(1)/2 + 1)*chunkR(i)*(yi - 1) + xi
        dataT(pointsR(i, 1) - 1 + zi, yi, xi) = recv_buff3(point)
      end do
      end do
    end do
    end do

    ! Edit by Hong Zhao
    deallocate (send_buff)
    ! End edit by Hong Zhao

  END SUBROUTINE getTdata_a2a

!**********************************************************************
  SUBROUTINE get_rhoR

    call dfftw_execute(plan_Rh_R)

  END SUBROUTINE get_rhoR

!**********************************************************************
  SUBROUTINE get_rhoT

    call dfftw_execute(plan_Th_T)

  END SUBROUTINE get_RhoT

!**********************************************************************
  SUBROUTINE get_rhoRh

    call dfftw_execute(plan_R_Rh)

  END SUBROUTINE get_rhoRh

!**********************************************************************
  SUBROUTINE get_rhoTh
    integer         :: i, j, k

    call dfftw_execute(plan_T_Th)

  END SUBROUTINE get_rhoTh

!**********************************************************************

  SUBROUTINE Finalize_PFFTW

    call dfftw_destroy_plan(plan_R_Rh)
    call dfftw_destroy_plan(plan_Rh_R)
    call dfftw_destroy_plan(plan_T_Th)
    call dfftw_destroy_plan(plan_Th_T)

    deallocate (chunkR, chunkT, pointsR, pointsT)
    deallocate (rhoR, rhoRh, rhoT, rhoTh)

  END SUBROUTINE Finalize_PFFTW

!***********************************************************************

end module ModPFFTW
