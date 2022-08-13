module har_oscillator
    use laguerrepol
    use, intrinsic :: iso_fortran_env, only:  dp=>real64
    implicit none
    
contains

function Rh(n,l,rb,r) result(val)
    integer,intent(in) :: l,n
    real(dp), intent (in) :: rb, r
    real(dp)::lg, val, s, c
    s = r/rb 
    call lf(n,real(l+0.5,8),s**2,lg)


    c = sqrt(2.0*gamma(1.0*n+1)/gamma(n+l+1.5))/(rb**1.5)
    val = c*(s**l)*lg*EXP(-0.5*s*s)
end function Rh
    
end module har_oscillator