module variational
    use hamiltonian_elements
    use constants
    use eigenvalues
    implicit none
    
contains

function rb_function (l,Msize,rb, EV) result(E)
    integer::l, Msize, rot_num, it_num
    real(8) :: rb, E, EV(Msize), A(Msize,Msize), v(Msize,Msize)

    call H(A,rb,l)

    call jacobi_eigenvalue ( Msize, A, it_max,v, EV, it_num, rot_num )

    !print*, it_num
    !print*, EV

    E=EV(1)

end function rb_function

function rb_dev(l, Msize,rb) result(val)
    integer::l, Msize
    real(8) :: rb, val, EV(Msize)

    val = (rb_function(l,Msize,rb+epsilon,EV)-rb_function(l,Msize,rb-epsilon,EV))/(2.0*epsilon)
    
end function rb_dev

function find_minimum(l,Msize,xmin,xmax) result(EV)
    integer, intent(in)::l, Msize
    integer:: i
    real(8):: EV(Msize),a,b,c,fa,fb,fc
    real(8), intent(in)::xmin,xmax

    a=xmin
    b=xmax
    fa = rb_dev(l, MSize,a)
    fb = rb_dev(l, MSize,b)
    
    if (fa*fb>0) then
        print*, 'Bad guess. Choose other range'
    else 
        do i=1,it_max
            c=(a+b)/2.0
            fc= rb_dev(l, MSize,c)
            if (abs(fc) <= epsilon) then
                exit
            else if (fa*fc < 0) then
                fb = fc
                b=c
            else 
                fa = fc
                a = c
            end if

            if(i==it_max) then
                print*, 'Did not converge'
            end if
        end do

        print*, 'the root is rb=', c
        fc=rb_function (l,Msize,c, EV)
        do i =1,Msize
            EV(i) = EV(i) + 4.0*m
        end do
    end if
end function find_minimum

end module variational