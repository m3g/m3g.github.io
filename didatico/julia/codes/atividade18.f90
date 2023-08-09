program atividade18
  implicit none
  integer :: i, ntrial
  double precision :: x, dfdx, deltax, deltaf, xbest, fbest
  deltax = 1.1d0 ! Step size
  x = 10.d0 ! Initial guess
  fbest = x**2
  xbest = x ! Save best point
  write(*,*) ' Initial point: '
  write(*,*) x, x**2
  ntrial = 1000
  do i = 1, ntrial
    ! Move x in the descent direction, with step deltax
    dfdx = 2.d0*x ! Computing the derivative
    x = x - dfdx * deltax ! Move x in the -f' direction
    ! Test new point
    deltaf = x**2 - fbest
    ! Write current point
    write(*,*) x, x**2, deltaf
    ! If the function decreased, save best point
    if ( deltaf < 0 ) then
      xbest = x
      fbest = x**2
    end if
  end do
end program atividade18
