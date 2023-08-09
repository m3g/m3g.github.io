program sim1
  integer :: i, nsteps
  double precision :: CA, dt, time, k1
  open(10,file='sim1.dat')
  nsteps = 1000 ! Number of steps
  dt = 1.d-1 ! Time-step
  k1 = 0.1d0 ! Velocity constant
  CA = 10.d0 ! Initial concentration
  time = 0.d0
  do i = 1, nsteps
    CA = CA - k1*CA*dt
    time = time + dt
    write(10,*) time, CA
  end do
end program sim1
