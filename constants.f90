module constants
    implicit none
    real(8),parameter :: m = 1.2185/2.0
    real(8),parameter :: hbar = 0.1973269804
    real(8), parameter :: rmax = 10000
    real(8), parameter :: dr = 0.001
    real(8), parameter ::alpha = 0.29 
    real(8), parameter ::sigma = 1.306
    integer :: nsteps = idint(rmax/dr)
    integer,parameter:: it_max = 500
    real(8), parameter :: epsilon = 10e-8
    real(8), parameter :: epsilon0 = 55.26349406
    real(8), parameter :: pi = 3.14159265359
    
end module constants