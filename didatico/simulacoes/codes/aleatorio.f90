program aleatorio

  integer :: seed(1)
  double precision :: random

  seed = 194851
  call random_seed(seed(1))

  call random_number(random)
  write(*,*) random
  call random_number(random)
  write(*,*) random

end program aleatorio
