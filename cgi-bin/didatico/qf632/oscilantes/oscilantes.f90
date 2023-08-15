
program oscilantes

  implicit none
  integer :: i, nsteps
  double precision :: A0, B0, C0, D0
  double precision :: k1, km1, k2, km2, k3, km3
  double precision :: time, dt, dAdt, dBdt, dCdt
  double precision, allocatable :: CA(:), CB(:), CC(:), CD(:)

  A0 = 10.
  B0 = 1.
  C0 = 0.5
  D0 = 0.

  k1 = 1.     ! A + B -> 2B  (rápida)      k2
  km1 = 1.d-4 ! 2B -> A + B (lenta)        km2
  k2 = 10.    ! B + C -> 2C (rápida)       k3
  km2 = 1.d-4 ! 2C -> B + C (lenta)        km3
  k3 = 1.     ! C -> D (rápida)            k4
  km3 = 1.d-4 ! D -> C (lenta)             km4
  
  time = 50.0
  dt = 0.0001
  
  nsteps=int(time/dt)

  allocate( CA(nsteps), CB(nsteps), CC(nsteps), CD(nsteps) )

  CA = 0.d0
  CA(1) = A0
  CB = 0.d0
  CB(1) = B0
  CC = 0.d0
  CC(1) = C0
  CD = 0 
  CD(1) = D0
  
  do i = 1, nsteps-1
      
    dAdt = -k1*CA(i)*CB(i) + km1*CB(i)**2

    dBdt = k1*CA(i)*CB(i) - km1*CB(i)**2 - k2*CB(i)*CC(i) + km2*CC(i)**2

    dCdt = k2*CB(i)*CC(i) - km2*CC(i)**2 - k3*CC(i) + km3*CD(i)

    CA(i+1) = CA(i) + dAdt * dt
    CB(i+1) = CB(i) + dBdt * dt
    CC(i+1) = CC(i) + dCdt * dt 
    CD(i+1) = A0 + B0 + C0 + D0 - CA(i+1) - CB(i+1) - CC(i+1)
      
  end do

  do i = 1, nsteps
    write(*,*) (i-1)*dt, CA(i), CB(i), CC(i), CD(i)
  end do

end

