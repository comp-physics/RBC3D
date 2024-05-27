program test
    
    use ModData
    use ModDataTypes
    use ModDataStruct
    use ModBasicMath

    implicit none

    real(WP), dimension(3, 2) :: arrayA
    real(WP), dimension(2, 3) :: arrayB
    integer :: i
    ! arrayA = reshape((/ 7, 3, 5, 2, 4, 3 /), shape(arrayA))
    arrayA = reshape((/ 1, 2, 1, 2, 1, -1 /), shape(arrayA))

    call Matrix_PseudoInvert(arrayA, arrayB)

    print *, "printing solution: "
    do i = 1, size(arrayB, 1)
        print *, arrayB(i, :)
    end do

end program test