program min1
  implicit none
  double precision :: x, dfdx, deltax, deltaf, xbest, fbest
  deltax = 0.1d0 ! Step size
  x = 14.1357d0 ! Initial guess
  fbest = x**2
  xbest = x ! Save best point
  write(*,*) ' Initial point: '
  write(*,*) x, x**2
  do
    ! Move x in the descent direction, with step deltax
    dfdx = 2.d0*x ! Computing the derivative
    x = x - deltax * dfdx/dabs(dfdx) ! Move x in the -f' direction
    ! Test new point
    deltaf = x**2 - fbest
    ! Write current point
    write(*,*) x, x**2, deltaf
    ! If the function decreased, save best point
    if ( deltaf < 0 ) then
      xbest = x
      fbest = x**2
    else
      write(*,*) ' Function is increasing. '
      write(*,*) ' Best solution found: x = ', xbest
      exit
    end if
  end do
end program min1
