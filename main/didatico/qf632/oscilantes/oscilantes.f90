!
! Programa Oscilantes: Material de Apoio para QF632
!
! L. Martinez. IQ-UNICAMP, Ago 2015
!
! https://m3g.github.io/
!

program oscilantes

  implicit none
  integer :: i, nsteps, iprint
  double precision :: A0, B0, C0, D0
  double precision :: k1, km1, k2, km2, k3, km3
  double precision :: time, dt, dAdt, dBdt, dCdt, dDdt
  double precision :: CA, CB, CC, CD, mass_balance, error
  character(len=200) :: record

  ! Lendo as condições da simulação:

  call getarg(1,record) ; read(record,*) A0       ! Concentração inicial de A
  call getarg(2,record) ; read(record,*) B0       ! Concentração inicial de B 
  call getarg(3,record) ; read(record,*) C0       ! Concentração inicial de C 
  call getarg(4,record) ; read(record,*) D0       ! Concentração inicial de D 
  call getarg(5,record) ; read(record,*) k1       ! Constante de velocidade
  call getarg(6,record) ; read(record,*) km1      ! Constante de velocidade 
  call getarg(7,record) ; read(record,*) k2       ! Constante de velocidade 
  call getarg(8,record) ; read(record,*) km2      ! Constante de velocidade 
  call getarg(9,record) ; read(record,*) k3       ! Constante de velocidade 
  call getarg(10,record) ; read(record,*) km3     ! Constante de velocidade 
  call getarg(11,record) ; read(record,*) time    ! Tempo total da simulação

  ! Número total de passos do processo iterativo
  
  dt = 0.01
  nsteps=int(time/dt)
  if ( nsteps > 10000 ) stop

  ! As concentrações serão impressas 200 vezes apenas

  iprint = nsteps / 200

  ! Concentrações iniciais

  CA = A0
  CB = B0
  CC = C0
  CD = D0

  ! A soma das concentrações deve ser constante, se a aproximação numérica for boa

  mass_balance = A0 + B0 + C0 + D0

  ! Titulo da tabela

  write(*,"( t1,'#',t10,'Tempo',t37,'CA',t63,'CB',t89,'CC',t115,'CD' )")

  ! Processo iterativo - simulação

  do i = 1, nsteps-1

    ! Equações diferenciais
      
    dAdt = -k1*CA*CB + km1*CB**2
    dBdt = k1*CA*CB - km1*CB**2 - k2*CB*CC + km2*CC**2
    dCdt = k2*CB*CC - km2*CC**2 - k3*CC + km3*CD
    dDdt = k3*CC - Km3*CD

    ! Propagação das concentrações no tempo

    CA = CA + dAdt * dt
    CB = CB + dBdt * dt
    CC = CC + dCdt * dt 
    CD = CD + dDdt * dt

    ! Testar se as aproximações numéricas são razoáveis

    error = abs(CA+CB+CC+CD-mass_balance)/mass_balance
    if ( error > 1.d-3 ) then
      write(*,"( '# ERRO: Simulação instável. Erro na soma das concentrações = ', e12.7 )") error
      stop
    end if
    if ( CA > mass_balance .or. CA < 0 .or. &
         CB > mass_balance .or. CB < 0 .or. &
         CC > mass_balance .or. CC < 0 .or. &
         CD > mass_balance .or. CD < 0 ) then 
      write(*,"( '# ERRO: Simulação instável. Concentrações ficaram sem sentido.' )")
      stop
    end if

    ! Se o resto da divisão de i por iprint for 0, escrever as concentrações
      
    if ( mod(i,iprint) == 0 ) write(*,*) (i-1)*dt, CA, CB, CC, CD

  end do

end program oscilantes

