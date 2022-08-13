module laguerrepol
    implicit none
    
contains
subroutine lf (n, alpha, x, fx)

    !*****************************************************************************80
    !
    !! LF_FUNCTION evaluates the Laguerre function Lf(n,alpha,x).
    !
    !  Recursion:
    !
    !    Lf(0,ALPHA,X) = 1
    !    Lf(1,ALPHA,X) = 1+ALPHA-X
    !
    !    Lf(N,ALPHA,X) = (2*N-1+ALPHA-X)/N * Lf(N-1,ALPHA,X) 
    !                      - (N-1+ALPHA)/N * Lf(N-2,ALPHA,X)
    !
    !  Restrictions:
    !
    !    -1 < ALPHA
    !
    !  Special values:
    !
    !    Lf(N,0,X) = L(N,X).
    !    Lf(N,ALPHA,X) = LM(N,ALPHA,X) for ALPHA integral.
    !
    !  Norm:
    !
    !    Integral ( 0 <= X < +oo ) exp ( - X ) * Lf(N,ALPHA,X)^2 dX
    !    = Gamma ( N + ALPHA + 1 ) / N!
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    10 March 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Milton Abramowitz, Irene Stegun,
    !    Handbook of Mathematical Functions,
    !    National Bureau of Standards, 1964,
    !    ISBN: 0-486-61272-4,
    !    LC: QA47.A34.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, the number of evaluation points.
    !
    !    Input, integer ( kind = 4 ) N, the highest order function to compute.
    !
    !    Input, real ( kind = rk ) ALPHA, the parameter.  -1 < ALPHA is required.
    !
    !    Input, real ( kind = rk ) X(M), the evaluation points.
    !
    !    Output, real ( kind = rk ) CX(1:M,0:N), the functions of 
    !    degrees 0 through N evaluated at the points X.
    !
      implicit none
    
      integer, parameter :: rk = kind ( 1.0D+00 )
    
      integer ( kind = 4 ) n
    
      real ( kind = rk ) alpha
      real ( kind = rk ) fx
      real ( kind = rk ) cx(0:n)
      integer ( kind = 4 ) i
      real ( kind = rk ) x
    
      if ( alpha <= -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LF_FUNCTION - Fatal error!'
        write ( *, '(a,g14.6)' ) '  The input value of ALPHA is ', alpha
        write ( *, '(a)' ) '  but ALPHA must be greater than -1.'
        stop
      end if
     
      if ( n < 0 ) then
        return
      end if
    
      cx(0) = 1.0D+00
    
      if ( n == 0 ) then
        return
      end if
    
      cx(1) = 1.0D+00 + alpha - x

    
      do i = 2, n
        cx(i) = ( &
          ( real ( 2 * i - 1, kind = rk ) + alpha - x ) * cx(i-1)   &
        + ( real (   - i + 1, kind = rk ) - alpha          ) * cx(i-2) ) &
          / real (     i,     kind = rk )
      end do
      
      fx = cx(n)

      return
    end
end module laguerrepol