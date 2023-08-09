program simplex
  implicit none
  integer :: i, j, iter, niter
  double precision :: random
  double precision :: x(3,2), f(3), ftemp, xtemp(2)
  double precision :: xav(2)
  double precision :: xtrial(2), ftrial
  double precision :: convcrit
  ! Generate initial points
  do i = 1, 3
    call random_number(random)
    x(i,1) = -10.d0 + 20.d0*random
    call random_number(random)
    x(i,2) = -10.d0 + 20.d0*random
    ! Necessary because x in computef is of the form x(2):
    xtemp(1) = x(i,1)
    xtemp(2) = x(i,2)
    call computef(xtemp,f(i))
  end do
  write(*,*) ' Initial points: '
  do i = 1, 3
    write(*,*) x(i,1), x(i,2), f(i)
  end do
  ! Convergence criterium desired
  convcrit = 1.d-10
  ! Maximum number of iterations
  niter = 10000
  do iter = 1, niter
    write(*,*) ' ------ ITERATION: ', iter
    ! Order the points from best to worst
    do i = 1, 3
      j = i
      do while( j > 1 .and. f(j-1) > f(j) )
        ftemp = f(j-1)
        f(j-1) = f(j)
        f(j) = ftemp
        xtemp(1) = x(j-1,1)
        xtemp(2) = x(j-1,2)
        x(j-1,1) = x(j,1)
        x(j-1,2) = x(j,2)
        x(j,1) = xtemp(1)
        x(j,2) = xtemp(2)
        j = j - 1
      end do
    end do
    ! Check convergence
    if ( f(3) - f(2) < convcrit .and. f(3) - f(1) < convcrit ) then
      write(*,*) ' Precision reached. '
      write(*,*) ' Best point found: ', x(1,1), x(1,2), ' f = ', f(1)
      stop
    end if
    ! Compute averge of best points
    xav(1) = 0.5d0*(x(1,1) + x(2,1))
    xav(2) = 0.5d0*(x(1,2) + x(2,2))
    ! Compute trial point
    xtrial(1) = x(3,1) + 2.d0*(xav(1)-x(3,1))
    xtrial(2) = x(3,2) + 2.d0*(xav(2)-x(3,2))
    call computef(xtrial,ftrial)
    ! If ftrial is better than f(3), replace point 3 with trial point
    if ( ftrial < f(3) ) then
      f(3) = ftrial
      x(3,1) = xtrial(1)
      x(3,2) = xtrial(2)
      write(*,*) ' Accepted point: ', x(3,1), x(3,2), ' f = ', f(3)
    else
      write(*,*) ' Function increased. Trying line search. '
      ! Try up to 10 different points in the 
      ! direction x(3)+gamma*(xtrial-x(3))
      do j = 1, 10
        call random_number(random)
        xtemp(1) = x(3,1) + random*(xtrial(1)-x(3,1))
        xtemp(2) = x(3,2) + random*(xtrial(2)-x(3,2))
        call computef(xtemp,ftemp)
        if ( ftemp < f(3) ) then
          f(3) = ftemp
          x(3,1) = xtemp(1)
          x(3,2) = xtemp(2)
          write(*,*) '   Line search succeeded at trial ', j
          write(*,*) '   New point: ', x(3,1), x(3,2), ' f = ', f(3)
          exit
        end if
      end do
      ! If the line search didn't find a better point, stop
      if ( ftemp > f(3) ) then
        write(*,*) ' End of search. '
        write(*,*) ' Best point found: ', x(1,1), x(1,2), ' f = ', f(1)
        stop
      end if
    end if
  end do
  write(*,*) ' Maximum number of trials reached. '
  write(*,*) ' Best point found: ', x(1,1), x(1,2), ' f = ', f(1)
end program simplex
!
! Compute the function value
!
subroutine computef(x,f)
  implicit none
  double precision :: x(2), f 
  f = x(1)**2 + x(2)**2
end subroutine computef

