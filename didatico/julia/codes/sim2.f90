program sim2
  implicit none
  integer :: i, nsteps
  double precision :: CA, CB, CA0, CB0
  double precision :: dt, time, k1, km1
  open(10,file='sim2.dat')
  nsteps = 1000 ! Number of steps
  dt = 1.d-1 ! Time-step
  k1 = 0.1d0 ! Velocity constant
  km1 = 0.05d0 ! Velocity constant
  CA0 = 10.d0 ! Initial concentration of A
  CB0 =  0.d0 ! Initial concentration of B
  time = 0.d0
  CA = CA0
  CB = CB0
  do i = 1, nsteps
    CA = CA - k1*CA*dt + km1*CB*dt
    CB = CA0 + CB0 - CA
    time = time + dt
    write(10,*) time, CA, CB
  end do
end program sim2
