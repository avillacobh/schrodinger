module hamiltonian_elements
    use constants
    use integration
    use har_oscillator
    implicit none

    private

    public H

contains

subroutine H(A, rb, l)
    real(8) :: A(:,:)
    real(8),intent(in) :: rb
    integer, intent(in) :: l 
    integer :: n,i,j 
    real(8):: y(nsteps)
    
    n = size(A, 1)


    do i = 1,n
        do j = i,n
            call faux(y,i,j,l,rb)
            A(i,j) = integrate(y,dr)  
            A(j,i) = A(i,j)      
        end do
    end do

    do concurrent (i=1:n)
        A(i,i) = A(i,i) + (2.0*l + 4.0*i + 3.0)*(hbar**2)/(rb**2)/(2.0*m)
    end do
end subroutine H


subroutine faux(y,n,np,l,rb)
    real(8) :: y(nsteps), x(nsteps)
    real(8), intent(in) :: rb
    integer, intent(in):: n,np,l
    integer:: i,j
    
    do j = 1,nsteps
        x(j) = (j-1)*dr
    end do

    do i=1,nsteps
        y(i) = (V_Po(x(i))-(hbar**2)*(x(i)**4)/(2.0*m*(rb**4)))*Rh(n,l,rb,x(i))*Rh(np,l,rb,x(i))
    end do
end subroutine faux


function V_Po(r) result(retval)
    real(8), intent(in) :: r
    real(8) :: retval

    retval =   sigma*(r**3) -4.0*alpha*r*hbar/(3.0)
    
end function V_Po
    
end module hamiltonian_elements