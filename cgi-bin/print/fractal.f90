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

integer :: ixdata, iydata, status, ndata, i, ix, iy,&
           nsteps, method
real :: xread, yread, stepfactor, rulemin, rulemax,&
        rulesize, step, a, b, r
real, allocatable :: x(:), y(:), logx(:), logy(:)
character(len=200) :: file, record

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

! Reading the name of the input file

open(10,file=file,status='old',iostat=status)
if(status /= 0) stop

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

if(ndata < 3) then
  write(*,*) ' Foram encontrados menos de tres pontos. '
  write(*,*) ' Impossivel continuar para este conjunto. '
  write(*,*) ' ---------------------------------------------------- '   
  stop
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

allocate( logx(nsteps), logy(nsteps) )

! Rule size parameters

step = rulemin + rulemin * 1.d-3
rulesize = step
stepfactor = ( rulemax /  rulemin )**(1./real(nsteps))

! Define method

method = 2

! Computing the fractal dimension by rule counting

if(method == 1) then 
  call rulecount(nsteps,ndata,x,y,step,stepfactor,logx,logy)
  call correlation(logx,logy,nsteps,r)
  write(*,*) ' Correlation coefficient, R = ', r
  call least_square(nsteps,logx,logy,a,b)
  write(*,*) ' Fractal dimension, D = ', b
end if

! Computing the fractal dimension by correlation matrix

if(method == 2) then
  call dimmatrix(nsteps,ndata,x,y,step,stepfactor,logx,logy)
  call correlation(logx,logy,nsteps,r)
  write(*,*) ' Correlation coefficient, R = ', r
  call least_square(nsteps,logx,logy,a,b)
  write(*,*) ' Fractal dimension, D = ', b
end if

! Write linear fit data

open(10,file='/home/leandro/public_html/fractal/files/check.dat')
do i = 1, nsteps
  write(10,*) logx(i), logy(i) 
end do
close(10)
open(10,file="/home/leandro/public_html/fractal/files/regr.dat")
write(10,*) logx(1), a + b * logx(1)
write(10,*) logx(nsteps), a + b * logx(nsteps)
close(10)

end

!
! Subroutine that computes the least squares to provide
! directly the fractal dimension
!
! y = a + bx
!

subroutine least_square(n,x,y,a,b)
implicit none
integer :: n, i
real :: a, b, x(*), y(*), a1, a2, b0, b1, d

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

!
! Function to compute the distance between two points
!

real function dij(x,y,i,j)

implicit none
real :: x(*), y(*)
integer :: i, j

dij = sqrt( ( x(j) - x(i) )**2 + ( y(j) - y(i) )**2 )

return
end

!
! Compute the fractal dimension using multiple
! rule sizes
!

subroutine rulecount(nsteps,ndata,x,y,step,stepfactor,logl,lognl)

implicit none
integer :: nsteps, i, j, ndata, nrules
real :: x(*), y(*), dnext, drefj, xref, yref, rulesize, xint, yint,&
        logl(*), lognl(*), step, stepfactor

! Computing the size of the line with different rule sizes

rulesize = step
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
  lognl(i) = -log10(rulesize*real(nrules))

  rulesize = rulesize * stepfactor
end do

return
end

!
! Function to compute the distance a point and a reference
!

real function drefj(xref,yref,x,y,j)

implicit none  
real :: xref, yref, x(*), y(*)
integer :: j

drefj = sqrt( ( x(j) - xref )**2 + ( y(j) - yref )**2 )

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

implicit none
real :: a, b, c, delta, x(*), y(*), xref, yref, rulesize, xint,& 
        yint, t
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

! Compute the correlation coefficient of the data, must
! be close to 1.0 in order that this makes any sense (and
! even if it does, the result may not make sense)

subroutine correlation(x,y,n,r)

implicit none
real :: x(*), y(*), sumxy, sumx2, sumy2, xmed, ymed, r
integer :: n, i

sumxy = 0.
sumx2 = 0.
sumy2 = 0.
xmed = 0.
ymed = 0.
do i = 1, n
  sumxy = sumxy + x(i) * y(i)
  sumx2 = sumx2 + x(i)**2
  sumy2 = sumy2 + y(i)**2
  xmed = xmed + x(i)
  ymed = ymed + y(i)
end do
xmed = xmed / real(n)
ymed = ymed / real(n)
r = ( sumxy - n * xmed * ymed )**2 / & 
    ( ( sumx2 - n * xmed**2 )*( sumy2 - n * ymed**2 ) )

return 
end

!
! Correlation matrix computation
!

subroutine dimmatrix(nsteps,ndata,x,y,step,stepfactor,logr,logc)

implicit none
integer :: j, k, nsteps, ndata, i
real :: step, x(*), y(*), rulesize, stepfactor, dij, logr(*), logc(*)
integer, allocatable :: ndist(:)

allocate ( ndist(nsteps) )
do k = 1, nsteps
  ndist(k) = 0
end do
do i = 1, ndata - 1
  do j = i + 1, ndata
    rulesize = step
    do k = 1, nsteps
      if ( dij(x,y,i,j) < rulesize ) ndist(k) = ndist(k) + 1
      rulesize = rulesize * stepfactor
    end do
  end do
end do
rulesize = step
do k = 1, nsteps
  logr(k) = log10(rulesize)
  logc(k) = log10(real(ndist(k)) / real(ndata*(ndata-1))) 
  rulesize = rulesize * stepfactor
end do

return
end


