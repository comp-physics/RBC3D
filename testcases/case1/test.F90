program test
    
    use ModDataTypes
    use ModBasicMath


    real(WP) :: mat2d(1: 17, 2), vec1d(1: 17)
    real(WP) :: b0, b1, b2, b11, b12, b22
    integer :: mat2drows, mat2dcols, j

    mat2d = reshape((/ 0., .0218166, .043633, 4., 5., 6., 7., 8., 9., 10., 11., 12., &
        13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., &
        27., 28., 29., 30., 31., 32., 33., 34. /), shape(mat2d))
    
    mat2drows = size(mat2d, 1)
    mat2dcols = size(mat2d, 2)
    print *, "mat2drows:", mat2drows, "mat2dcols:", mat2dcols

    do j = 1, mat2drows
        ! mat2d(j, :) = (/j, j*2/)
        print *, mat2d(j, :)
    end do

    



end program test