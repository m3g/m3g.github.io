program randomsearch
  implicit none
  integer :: i, ntrial
  double precision :: random, x(2), f, fbest, xbest(2)
  ! Test 10000 points
  ntrial = 10000
  fbest = 1.d30
  do i = 1, ntrial
    call random_number(random)
    x(1) = -10.d0 + 20.d0*random
    call random_number(random)
    x(2) = -10.d0 + 20.d0*random
    call computef(x,f)
    if ( f < fbest ) then
      fbest = f
      xbest(1) = x(1)
      xbest(2) = x(2)
      write(*,*) i, ' New best point: ', x(1), x(2), ' f = ', f
    end if
  end do
  write(*,*) ' Best point found: ', xbest(1), xbest(2), ' f = ', fbest
end program randomsearch
!
! Subroutine that computes the function value
!
subroutine computef(x,f)
  implicit none
  double precision :: x(2), f
  f = x(1)**2 + x(2)**2
end subroutine computef
