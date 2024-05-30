program test
    
    use ModData
    use ModDataTypes
    use ModDataStruct
    use ModBasicMath

    ! use ModDataTypes
    ! use ModDataStruct
    ! use ModRbc
    ! use ModWall
    ! use ModConf
    ! use ModData
    ! use ModIO
    ! use ModBasicMath
    ! use ModPostProcess

    implicit none

    integer :: n_points = 17, i
    real(WP) :: A(17, 2)
    real(WP) :: b(17)
    real(WP) :: x0 = 0., x1 = 0., x2 = 0., x11 = 0., x12 = 0., x22 = 0., timeit
    real(WP) :: x, y

    call InitMPI

    print *, "n_points", n_points

    A = reshape((/0.000E+00, 2.182E-02, 4.363E-02, 1.543E-02, &
        3.085E-02, 0.000E+00, 0.000E+00, -1.543E-02, -3.085E-02, & 
        -2.182E-02, -4.363E-02, -1.543E-02, -3.085E-02, -4.008E-18, &
        -8.015E-18, 1.543E-02, 3.085E-02, 0.000E+00, 0.000E+00, 0.000E+00, &
        1.543E-02, 3.085E-02, 2.182E-02, 4.363E-02, 1.543E-02, 3.085E-02, &
        2.672E-18, 5.344E-18, -1.543E-02, -3.085E-02, -2.182E-02, -4.363E-02, &
        -1.543E-02, -3.085E-02/), shape(A))
    b = (/1.1284, 1.1303, 1.1343, 1.1304, 1.1350, 1.1300, 1.1348, 1.1290, &
            1.1322, 1.1284, 1.1303, 1.1289, 1.1321, 1.1299, 1.1347, 1.1303, &
            1.1349/)

    print *, "Matrix A: "
    do i = 1, n_points
        print *, A(i, :)
    end do

    print *, "b: ", b

    call QuadFit_2D(A, b, x0, x1, x2, x11, x12, x22, timeit)

    print *, 'Computed coefficients:'
    print *, 'x0  = ', x0
    print *, 'x1  = ', x1
    print *, 'x2  = ', x2
    print *, 'x11 = ', x11
    print *, 'x12 = ', x12
    print *, 'x22 = ', x22

    print *, "time taken = ", timeit

    ! deallocate(x, f)

    call FinalizeMPI
end program test