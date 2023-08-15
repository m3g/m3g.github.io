
program apos

  implicit none
  integer :: i, j, nanos, narg, iargc
  double precision :: juros, inflacao, salario, contribuicao
  double precision :: juros_mes
  double precision :: arrecadacao
  character(len=20) :: record

  nanos = 49
  juros = 6.
  inflacao = 5.
  salario = 880.
  contribuicao = 176.

  if ( nanos > 100 ) then
    write(*,*) ' anos demais '
    stop
  end if

  narg = iargc()
  if ( narg > 0 ) then
    call getarg(1,record)
    read(record,*) nanos
    call getarg(2,record)
    read(record,*) juros
    call getarg(3,record)
    read(record,*) inflacao
    call getarg(4,record)
    read(record,*) salario
    call getarg(5,record)
    read(record,*) contribuicao
  end if
  juros = 1. + juros / 100.
  inflacao = 1. + inflacao / 100.

  write(*,*)
  write(*,*) '<b> Dados: </b>'
  write(*,*)

  juros_mes = juros**(1./12.)
  write(*,"(a,f12.2,a)") ' Salário inicial = ', salario
  write(*,"(a,f12.2,a)") ' Contribuição inicial = ', contribuicao
  write(*,"(a,i3)") ' Anos de contribuição = ', nanos
  write(*,"(a,f12.2,a)") ' Juros por ano = ', (juros-1.)*100, '%'
  write(*,"(a,f12.2,a)") ' Inflação por ano = ', (inflacao-1.)*100, '%'

  arrecadacao = 0.
  do i = 1, nanos
    do j = 1, 12
      arrecadacao = arrecadacao*juros_mes + contribuicao
    end do
    arrecadacao = arrecadacao + contribuicao ! 13o
    contribuicao = contribuicao*inflacao
    salario = salario*inflacao
  end do

  write(*,*)
  write(*,*) '<b> Resultados: </b>'
  write(*,*)

  write(*,"(a,f12.2)") ' Contribuição no último mês = ', contribuicao
  write(*,"(a,f12.2)") ' Salário no último mês de contribuição = ', salario
  write(*,"(a,f12.2)") ' Total arrecadado mais rendimentos = ', arrecadacao

  i = 0
  do while( arrecadacao*juros_mes-salario > 0. ) 
    i = i + 1
    arrecadacao = arrecadacao*juros_mes - salario
    if ( mod(i-1,12) == 0 ) then
      salario = salario*inflacao
    end if
    if ( i > 50*12 ) then
      write(*,*) ' Vai durar mais que 50 anos '
      stop
    end if
  end do
  write(*,"(a,i3)") ' Número de meses que a arrecadação sustenta: = ', i
  write(*,"(a,f12.2)") ' Anos até o fim (com 13o) = ', i/13.

  ! repetindo para reportar detalhes

  nanos = 49
  juros = 8.
  inflacao = 5.
  salario = 880.
  contribuicao = 176.

  if ( nanos > 100 ) then
    write(*,*) ' anos demais '
    stop
  end if

  if ( narg > 0 ) then
    call getarg(1,record)
    read(record,*) nanos
    call getarg(2,record)
    read(record,*) juros
    call getarg(3,record)
    read(record,*) inflacao
    call getarg(4,record)
    read(record,*) salario
    call getarg(5,record)
    read(record,*) contribuicao
  end if
  juros = 1. + juros / 100.
  inflacao = 1. + inflacao / 100.


  write(*,*)
  write(*,*) ' Detalhes: Contribuicoes '
  write(*,*)
  arrecadacao = 0.
  do i = 1, nanos
    do j = 1, 12
      arrecadacao = arrecadacao*juros_mes + contribuicao
    end do
    arrecadacao = arrecadacao + contribuicao ! 13o
    contribuicao = contribuicao*inflacao
    salario = salario*inflacao
    write(*,"(a,i3,a,f12.2,a,f12.2,a,f12.2)") &
               ' Ano ', i, ': Arrecadacao = ', arrecadacao, &
                            ' Salario = ', salario, &
                            ' Contribuicao = ', contribuicao
  end do
  write(*,*)
  write(*,*) ' Detalhes: Beneficios '
  write(*,*)
  i = 0
  do while( arrecadacao*juros_mes-salario > 0. ) 
    i = i + 1
    arrecadacao = arrecadacao*juros_mes - salario
    if ( mod(i-1,12) == 0 ) then
      salario = salario*inflacao
    end if
    write(*,"(a,i3,a,f12.2,a,f12.2)") ' Mes ', i,': Reservas = ', arrecadacao, ' Salario = ', salario  
    if ( i > 50*12 ) then
      write(*,*) ' Vai durar mais que 50 anos '
      stop
    end if
  end do
  write(*,*)
  write(*,"(a,f12.2)") ' anos até o fim (com 13o) = ', i/13.

end program apos

