module subroutines_matrix
    implicit none
    
    public print_matrix, filling_matrix, trace

    

contains

subroutine print_matrix (A)
    real(8), intent(in) :: A(:,:)

    integer:: i

    do i =1, size(A,1)
        print *, A(i,:)
    end do

end subroutine print_matrix

subroutine filling_matrix(A)
    real(8), intent(inout) :: A (:,:)

    integer :: i, j 

    do i = 1, size(A,1)
        do j = 1 , size(A,2)
            A(i,j) = (i + j)/3.0
        end do
    end do

end subroutine filling_matrix

function trace (A) result(tr)
    implicit none 

    real(8), intent(in) :: A(:,:)
    real(8) :: tr
    integer :: i

    tr = 0.0

    if (size(A,1) < size(A,2)) then
        do i = 1 , size(A,1)
            tr = tr + A(i,i)
        end do
    else 
        do i = 1 , size(A,2)
            tr = tr + A(i,i)
        end do     
    end if

end function trace
    
end module subroutines_matrix

