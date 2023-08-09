module reaction
  ! Initial concentrations
  double precision :: a0, b0
  ! Time step between experimental measurements
  double precision :: dt
  ! Experimental concentration of the substrate as a function of time
  double precision :: aexp(1000)
end module reaction
!
! Program modelAB
!
program modelAB
  use reaction
  implicit none
  integer :: i, info
  double precision :: x0(2), x(2), f
  double precision :: random
  external computef

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

  call random_number(random)
  x0(1) = random
  call random_number(random)
  x0(2) = random
  write(*,*) ' Initial guess: ', x0(1), x0(2)
  call minimize_2d(x0,computef,x,f,info)
  write(*,*) ' Solution: ', x(1), x(2), f

end program modelAB
!
! Compute the function value (does a MK simulation and computes the variance
! of kinetic data relative to "experimental" kinetic data).
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
    a = dmax1(a,0.d0)
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





