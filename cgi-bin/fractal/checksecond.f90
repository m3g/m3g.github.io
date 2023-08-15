!
! Program checkfirst
!
! Computes the fractal dimension of a curve y(x), by measuring
! the length of the curve with different rule sizes
!
! L. Martinez, UNICAMP, Sep 25, 2008, leandromartinez98@gmail.com.
!

program fractal

implicit none

! Variables

integer :: status, ndata, i, xcol, ycol, j
real :: dmin, dmax, dij, xread, xmin, xmax, ymin, ymax
real, allocatable :: x(:), y(:)
character(len=500) :: file, record

call getarg(1,file)
call getarg(2,record)
read(record,*) xcol
call getarg(3,record)
read(record,*) ycol

open(10,file=file,status='old')
ndata = 1
do
  read(10,"( a500 )",iostat=status) record
  if(status /= 0) exit
  read(record,*,iostat=status) (xread,i=1,xcol)
  if(status /= 0) cycle
  read(record,*,iostat=status) (xread,i=1,ycol)
  if(status /= 0) cycle
  ndata = ndata + 1
end do 
close(10)
ndata = ndata - 1

allocate ( x(ndata), y(ndata) )

open(10,file=file,status='old')
ndata = 1
do
  read(10,"( a500 )",iostat=status) record
  if(status /= 0) exit
  read(record,*,iostat=status) (x(ndata),i=1,xcol)
  if(status /= 0) cycle
  read(record,*,iostat=status) (y(ndata),i=1,ycol)
  if(status /= 0) cycle
  ndata = ndata + 1
end do 
ndata = ndata - 1
close(10)

! Computing maximum intervals to normalize the data

xmin = x(1)
xmax = x(1)
ymin = y(1)
ymax = y(1)
do i = 2, ndata
  xmin = amin1(x(i),xmin)
  xmax = amax1(x(i),xmax)
  ymin = amin1(y(i),ymin)
  ymax = amax1(y(i),ymax)
end do
if ( xmax - xmin > 1.d-10 ) then
  do i = 1, ndata
    x(i) = ( x(i) - xmin ) / ( xmax - xmin )
  end do
end if
if ( ymax - ymin > 1.d-10 ) then
  do i = 1, ndata     
    y(i) = ( y(i) - ymin ) / ( ymax - ymin )
  end do
end if

ndata = ndata - 1
dmin = dij(x,y,1,2)
dmax = 0.
do i = 1, ndata - 1
  dmin = amin1(dmin,dij(x,y,i,i+1))
  do j = i + 1, ndata
    dmax = amax1(dmax,dij(x,y,i,j))
  end do
end do      

write(*,*) dmin
write(*,*) dmax

! Writing new data file with normalized data

open(10,file='/home/leandro/public_html/fractal/files/data2.dat')
write(10,"( '# xmax, xmin = ', 2(tr2,f12.5) )") xmax, xmin
write(10,"( '# ymax, ymin = ', 2(tr2,f12.5) )") ymax, ymin
do i = 1, ndata
  write(10,*) x(i), y(i)
end do
close(10)    

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





















