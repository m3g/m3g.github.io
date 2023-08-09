module reaction
  ! Initial concentrations
  double precision :: a0, b0
  ! Time step between experimental measurements
  double precision :: dt
  ! Experimental concentration of the substrate as a function of time
  double precision :: aexp(1000)
end module reaction
!
! Program kineticmodel
!
program kineticmodel
  use reaction
  implicit none
  integer :: i, j, iter, niter
  double precision :: random
  double precision :: x(3,2), f(3), ftemp, xtemp(2)
  double precision :: xav(2)
  double precision :: xtrial(2), ftrial
  double precision :: convcrit
  ! Initial concentrations
  a0 = 10.d0
  b0 = 0.d0
  ! Time-step
  dt = 1.d-1
  ! Read 'experimental' data
  open(10,file='kineticmodel.dat',action='read') 
  do i = 1, 1000
    read(10,*) aexp(i)
  end do
  close(10)
  niter = 1000
  ! Initial guess for variables (the rate constants)
  do i = 1, 3
    call random_number(random)
    x(i,1) = random ! k1
    call random_number(random)
    x(i,2) = random ! k-1
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
  iteration : do iter = 1, niter
    write(*,*) ' ------ ITERATION: ', iter
    ! Order the points from best to worst
    call sort(x,f)
    ! Check convergence
    if ( f(3) - f(2) < convcrit .and. &
         f(3) - f(1) < convcrit ) then
      write(*,*) ' Precision reached. '
      exit iteration
    end if
    ! Compute average of best points
    xav(1) = 0.5d0*(x(1,1)+x(2,1))
    xav(2) = 0.5d0*(x(1,2)+x(2,2))
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
        exit iteration
      end if
    end if
  end do iteration
  write(*,*) ' Best point found: ', x(1,1), x(1,2), ' f = ', f(1)
end program kineticmodel
!
! Subroutine to sort vectors according to their function values
!
subroutine sort(x,f)
  implicit none
  integer :: i, j
  double precision :: x(3,2), f(3)
  double precision :: xtemp(2), ftemp
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
end subroutine sort
!
! Compute the function value (does a simulation)
!
subroutine computef(x,f)
  use reaction
  implicit none
  integer :: i
  double precision :: x(2), f, asim(1000), a, b
  !
  ! Simulate the reaction using the parameters given
  !
  a = a0
  b = b0
  do i = 1, 1000
    a = a + (-1.d0*x(1)*a + x(2)*b)*dt
    b = a0 + b0 - a
    asim(i) = a
  end do
  !
  ! Compute the deviation relative to experimental data
  !
  f = 0.d0
  do i = 1, 1000
    f = f + ( aexp(i) - asim(i) )**2
  end do
end subroutine computef
