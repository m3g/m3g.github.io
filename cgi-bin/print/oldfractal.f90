!
! Program fractal
!
! Computes the fractal dimension of a curve y(x), by measuring
! the length of the curve with different rule sizes
!
! L. Martinez, UNICAMP, Sep 25, 2008, leandromartinez98@gmail.com.
!

program fractal

implicit none

! Variables

integer :: ixdata, iydata, status, ndata, i, j, ix, iy, nrules,&
           nsteps, narg
real :: xread, yread, stepfactor, rulemin, rulemax,&
        dmin, dmax, dij, xref, yref, xint,&
        yint, dnext, rulesize, step, drefj, a, b,&
        xmed, ymed, sumxy, sumx2, sumy2, r
real, allocatable :: x(:), y(:), logn(:), logl(:)
character(len=200) :: file, record
logical exit, check, direct

narg = iargc()
direct = .false.
if( narg > 0 ) then
  direct = .true.
  call getarg(1,file)
  call getarg(2,record)
  read(record,*) ix 
  call getarg(3,record)
  read(record,*) iy
  call getarg(4,record)
  read(record,*) rulemin 
  call getarg(5,record)
  read(record,*) rulemax 
  call getarg(6,record)
  read(record,*) nsteps
end if

! Reading the name of the input file

if ( .not. direct ) then
  write(*,*) ' ----------------------------------------------------- '
  write(*,*) ' Programa que calcula a dimensao fractal de uma curva. '
  write(*,*) ' ** O resultado so tem sentido se o coeficiente de **  '
  write(*,*) ' ** correlacao for proximo de 1.                   **  '
  write(*,*) ' ----------------------------------------------------- '
end if

exit = .false.
do while( .not. exit )

  if( .not. direct ) then
    write(*,"('  Nome do arquivo de dados (ASCII): ')",advance="no")
    read(*,*) file
    open(10,file=file,status='old',iostat=status)
    if(status /= 0) then
      write(*,*) ' Arquivo nao encontrado. '
      write(*,*) ' ---------------------------------------------------- ' 
      cycle
    end if
    write(*,"('  Coluna contendo os dados X: ')",advance="no")
    read(*,*) ix
    write(*,"('  Coluna contendo os dados Y: ')",advance="no")
    read(*,*) iy
  else
    open(10,file=file,status='old',iostat=status)
    if(status /= 0) stop
  end if
    

! Reading the number of data points in the line

  ndata = 0
  do 
    read(10,"( a200 )",iostat=status) record
    if(status /= 0) exit
    read(record,*,iostat=ixdata) (xread, i=1,ix)
    read(record,*,iostat=iydata) (yread, i=1,iy)
    if(ixdata == 0 .and. iydata == 0) then
      ndata = ndata + 1
    end if
  end do
  close(10)
  if( .not. direct ) write(*,*) ' Numero de dados: ', ndata

  if(ndata < 3) then
    write(*,*) ' Foram encontrados menos de tres pontos. '
    write(*,*) ' Impossivel continuar para este conjunto. '
    write(*,*) ' ---------------------------------------------------- '   
    if( direct ) stop
    cycle 
  end if

! Allocating arrays containing data points, and actually
! reading their values

  allocate( x(ndata), y(ndata) )
  ndata = 0
  open(10,file=file,status='old')
  do 
    read(10,"( a200 )",iostat=status) record
    if(status /= 0) exit
    read(record,*,iostat=ixdata) (xread, i=1,ix)
    read(record,*,iostat=iydata) (yread, i=1,iy)
    if(ixdata == 0 .and. iydata == 0) then
      ndata = ndata + 1 
      x(ndata) = xread
      y(ndata) = yread
    end if
  end do
  close(10)

! Obtaining the minimum and maximum distances between points

  if( .not. direct ) then
    dmin = dij(x,y,1,2)
    dmax = 0.
    do i = 1, ndata - 1
      dmin = amin1(dmin,dij(x,y,i,i+1))
      do j = i + 1, ndata
        dmax = amax1(dmax,dij(x,y,i,j))
      end do
    end do
    write(*,*) ' Minima distancia entre pontos: ', dmin
    write(*,*) ' Maxima distancia entre pontos: ', dmax

! Checking the range of rule sizes desired

    write(*,"('  Tamanho de regua minimo: ')",advance="no")
    read(*,*) rulemin
    write(*,"('  Tamanho de regua maximo: ')",advance="no")
    read(*,*) rulemax
    write(*,"('  Numero de reguas diferentes: ')",advance="no")
    read(*,*) nsteps
  end if
  allocate( logn(nsteps), logl(nsteps) )

! Checking if the user wants an output file
  
  check = .false.
  if( .not. direct ) then
    write(*,"( '  Deseja um arquivo de saida (s/n)?' )",advance="no")
    read(*,*) file
    if( file == "s" ) check = .true.
    if(check) then  
      write(*,"( '  Nome do arquivo de saida: ' )",advance="no")
      read(*,*) file
    end if
  else
    check = .true.
    file = './files/check.dat'
  end if

! Rule size parameters

  step = rulemin
  rulesize = step
  stepfactor = ( rulemax /  rulemin )**(1./real(nsteps))

! Computing the size of the line with different rule sizes

  if(check) open(10,file=file)
  do i = 1, nsteps

! Write progress

    nrules = 0
    xref = x(1)
    yref = y(1)
    j = 2 
    do while ( j <= ndata )

! Loop while the distance is smaller than the distance to next point

      dnext = drefj(xref,yref,x,y,j)
      do while (dnext < rulesize .and. j <= ndata ) 
        j = j + 1
        if(j <= ndata) dnext = drefj(xref,yref,x,y,j)   
      end do

! Sum current rulesize and compute next interesection

      if(j <= ndata) then
        nrules = nrules + 1
        call intersect(xref,yref,x,y,j,rulesize,xint,yint)
        xref = xint
        yref = yint
      end if

    end do
    logl(i) = log10(rulesize)
    logn(i) = -log10(real(nrules))
    if(check) &
    write(10,"( 2(tr2,e17.10) )") log10(rulesize), -log10(real(nrules))

! If already only one rule measure, stop

    if(nrules == 1) then
      write(*,*) ' Ultima medida com regua (n=1) = ', rulesize
      nsteps = i
      exit
    end if

    rulesize = rulesize * stepfactor
  end do

! Compute the correlation coefficient of the data, must
! be close to 1.0 in order that this makes any sense
  
  sumxy = 0.
  sumx2 =0.
  sumy2 = 0.
  xmed = 0.
  ymed = 0.
  do i = 1, nsteps
    sumxy = sumxy + logl(i) * logn(i)
    sumx2 = sumx2 + logl(i)**2
    sumy2 = sumy2 + logn(i)**2
    xmed = xmed + logl(i)
    ymed = ymed + logn(i)
  end do
  xmed = xmed / real(nsteps)
  ymed = ymed / real(nsteps)
  r = ( sumxy - nsteps * xmed * ymed )**2 / & 
      ( ( sumx2 - nsteps * xmed**2 )*( sumy2 - nsteps * ymed**2 ) )
  write(*,*) ' Correlation coefficient, R = ', r

! Computing the linear fit

  call least_square(nsteps,logl,logn,a,b)
  write(*,*) ' Fractal dimension, D = ', b

  close(10)
  if( .not. direct ) then
    write(*,*) ' ---------------------------------------------------- '
    write(*,"( '  Calcular para outra curva (s/n)?' )",advance="no")
    read(*,*) file
    if( file == 'n' ) then
      exit = .true.
    else
      write(*,*) ' ---------------------------------------------------- ' 
      deallocate ( x, y, logn, logl )
    end if
  else 
   exit = .true. 
  end if
end do

end

!
! Function to compute the distance a point and a reference
!

real function drefj(xref,yref,x,y,j)
  
real :: xref, yref, x(*), y(*)
integer :: j

drefj = sqrt( ( x(j) - xref )**2 + ( y(j) - yref )**2 )

return
end

!
! Function to compute the distance between two points
!

real function dij(x,y,i,j)
  
real :: x(*), y(*)
integer :: i, j

dij = sqrt( ( x(j) - x(i) )**2 + ( y(j) - y(i) )**2 )

return
end

!
! Subroutine intersect
!
! Finds the intersection between the circle of radius
! stepsize centered in current reference and the line
! connecting the next consecutive points
!

subroutine intersect(xref,yref,x,y,j,rulesize,xint,yint)

real :: a, b, c, delta, x(*), y(*), xref, yref, rulesize, xint,& 
        yint 
integer :: j

! Parameters of the parabola

a = ( x(j) - x(j-1) )**2 + ( y(j) - y(j-1) )**2
b = 2. * ( x(j) - x(j-1) ) * ( x(j-1) - xref ) + &
    2. * ( y(j) - y(j-1) ) * ( y(j-1) - yref )
c = ( x(j-1) - xref )**2 + ( y(j-1) - yref )**2 - rulesize**2
delta = b**2 - 4*a*c

! Check if there are solutions (there must be). If delta < 0,
! probably it is a problem of numerical accuracy when the minimum
! of the parabola is the only root. Therefore, in this case we
! check the value of the function at this minimum, before declaring
! that there is an error.

if ( delta < 0 ) then
  if ( c - 0.5 * b**2 / ( 2.*a ) > 1.d-6 ) then
    write(*,*) 
    write(*,*) ' ERROR: Delta < 0, no intersection. There is an error. '
    write(*,*) ' Delta = ', delta 
    write(*,*) ' Minimum at t = ', -b / (2.*a)
    write(*,*) ' Value at minimum, y(t) = ', c - 0.5 * b**2 / ( 2.*a )
    stop
  else
    t = -b / (2.*a) 
  end if
else 
  delta = sqrt(delta)
  t = amax1( ( -b + delta ) / ( 2.*a ), &
             ( -b - delta ) / ( 2.*a ) )
end if
xint = x(j-1) + t * ( x(j) - x(j-1) )  
yint = y(j-1) + t * ( y(j) - y(j-1) )  

return
end

!
! Subroutine that computes the least squares to provide
! directly the fractal dimension
!
! y = a + bx
!

subroutine least_square(n,x,y,a,b)
integer :: n, i
real :: a, b, x(*), y(*), a1, a2, b0, b1

a1 = 0
a2 = 0
b0 = 0
b1 = 0
do i = 1, n
  a1 = a1 + x(i)
  a2 = a2 + x(i) * x(i)
  b0 = b0 + y(i)
  b1 = b1 + y(i) * x(i)
end do
a1 = a1 / n
a2 = a2 / n
b0 = b0 / n
b1 = b1 / n
d = a1 * a1 - a2
a = a1 * b1 - a2 * b0
a = a / d
b = a1 * b0 - b1
b = b / d

return
end



