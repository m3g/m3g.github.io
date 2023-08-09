program min2
  implicit none
  double precision :: x, dfdx, deltax, f
  double precision :: xtrial, ftrial
  deltax = 0.1d0 ! Step size
  x = 10.d0 ! Initial guess
  f = x**2
  write(*,*) ' Initial point: '
  write(*,*) x, f
  write(*,*) ' x, f, deltax, dfdx : '
  do
    ! Move x in the descent direction, with step deltax
    dfdx = 2.d0*x ! Computing the derivative
    ! If the derivative is very small, stop
    if ( dabs(dfdx) < 1.d-10 ) then
      write(*,*) ' Critical point found. '
      write(*,*) ' x = ', x, ' f = ', f, ' dfdx = ', dfdx
      stop
    end if
    ! Computing trial point
    xtrial = x - deltax * dfdx ! Move x in the -f' direction
    ! Compute function value at trial point
    ftrial = xtrial**2
    ! If the function decreased, accept trial point and increase step
    if ( ftrial < f ) then
      x = xtrial
      f = ftrial
      deltax = deltax * 2.d0
      write(*,*) ' Accepted: ', x, f, deltax, dfdx
    else
      deltax = deltax / 2.d0
      write(*,*) ' Not accepted: ', xtrial, ftrial, deltax, dfdx
    end if
  end do
end program min2
