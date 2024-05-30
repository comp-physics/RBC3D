program test
    
    use ModData
    use ModDataTypes
    use ModDataStruct
    use ModBasicMath

    implicit none

    real(WP) :: x_vals(3)
    real(WP) :: f_vals(3)
    real(WP) :: c0, c1, c2

    x_vals = (/-2, 1, 6/)
    f_vals = (/5, 8, 4/)

    call QuadFit_1D(x_vals, f_vals, c0, c1, c2)

    print *, "c0: ", c0
    print *, "c1: ", c1
    print *, "c2: ", c2

end program test