module integration
    use constants
    implicit none
    
contains

pure function integrate(y,dx) result(sum)
    real(8), intent(in) :: y(:), dx
    real(8) :: sum
    integer:: ysize,i
    
    sum = 0.0 
    ysize = size(y)

    do concurrent (i = 1:ysize-1)
        sum = sum + 0.5*(y(i)+y(i+1))*dx
    end do
    
end function integrate

end module integration