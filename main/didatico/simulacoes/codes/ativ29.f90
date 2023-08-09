program ativ29
  implicit none
  integer :: i, n
  double precision :: x(1000), g(1000), deltax, f
  double precision :: xtrial(1000), ftrial, gnorm
  deltax = 0.1d0 ! Step size
  ! Number of variables
  n = 1000
  ! Initial guess
  do i = 1, n
    x = 10.d0 
  end do
  call computef(n,x,f)
  write(*,*) ' f at initial point: ', f
  do
    ! Compute the gradient
    call computeg(n,x,g)
    ! If the derivative is very small, stop
    gnorm = 0.d0
    do i = 1, n
      gnorm = gnorm + g(i)**2
    end do
    if ( gnorm < 1.d-10 ) then
      write(*,*) ' Critical point found. '
      write(*,*) ' f = ', f, ' gnorm = ', gnorm
      stop
    end if
    ! Computing trial point
    do i = 1, n
      xtrial(i) = x(i) - deltax * g(i) ! Move x in the -f' direction
    end do
    ! Compute function value at trial point
    call computef(n,xtrial,ftrial)
    ! If the function decreased, accept trial point and increase step
    if ( ftrial < f ) then
      do i = 1, n
        x(i) = xtrial(i)
      end do
      f = ftrial
      deltax = deltax * 2.d0
      write(*,*) ' Accepted: ', f, deltax
    else
      deltax = deltax / 2.d0
      write(*,*) ' Not accepted: ', ftrial, deltax
    end if
  end do
end program ativ29
!
! Subroutine that computes the function value
!
subroutine computef(n,x,f)
  implicit none
  integer :: i, n
  double precision :: x(n), f
  f = 0.d0
  do i = 1, n
    f = f + x(i)**2
  end do
end subroutine computef
!
! Subroutine that computes the gradient
!
subroutine computeg(n,x,g)
  implicit none
  integer :: i, n
  double precision :: x(n), g(n)
  do i = 1 , n
    g(i) = 2.d0*x(i)
  end do
end subroutine computeg







