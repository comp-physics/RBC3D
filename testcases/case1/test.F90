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

    integer :: n_points = 4, i
    real(WP) :: A(4, 2)
    real(WP) :: b(4)
    real(WP) :: x0 = 0., x1 = 0., x2 = 0., x11 = 0., x12 = 0., x22 = 0., timeit
    real(WP) :: x, y

    call InitMPI
    ! allocate arrays
    ! allocate(x(n_points, 2))
    ! allocate(f(n_points))

    print *, "n_points", n_points

    ! at each point x, y, b(x, y) = 1 + 2x + 3y + 4x^2 + 5xy + 6y^2
    do i = 1, n_points
        A(i, 1) = real(i, WP) / 4. ! x
        A(i, 2) = real(i, WP) / 4. ! y
        x = A(i, 1)
        y = A(i, 2)
        b(i) = .1 + .2*x + .3*y &
            + .4*(x**2) + .5*x*y + .6*(y**2)
    end do

    print *, "Matrix A: "
    do i = 1, n_points
        print *, "(", A(i, 1), ",", A(i, 2), ",", b(i), ")"
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