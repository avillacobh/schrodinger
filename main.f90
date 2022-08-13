program main
  use variational
  use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
  implicit none

  integer :: i,l
  integer, parameter :: n=15
  real(8) ::E0(n) = 0.0, E1(n) = 0.0, E2(n) = 0.0 ,xmin=0.1, xmax=5.0

  

  l=0
  E0=find_minimum(l,n,xmin,xmax)

  print*, 'l = ', l
  do i = 1, n
    print*,i , E0(i) 
  end do

  l=1
  E1=find_minimum(l,n,xmin,xmax)
  print*, 'l = ', l
  do i = 1, n
    print*,i , E1(i) 
  end do

  l=2
  E2=find_minimum(l,n,xmin,xmax)
  print*, 'l = ', l
  do i = 1, n
    print*,i , E2(i) 
  end do

  open(1, file = 'table.txt')

  write(1,*) 'State ', ' Mass' 
  write(1,*) '1S ', E0(1)
  write(1,*) '1P ', E1(1)
  write(1,*) '2S ', E0(2)
  write(1,*) '1D ', E2(1)
  write(1,*) '2P ', E1(2)
  write(1,*) '3S ', E0(3)

  close(1)

end program main
