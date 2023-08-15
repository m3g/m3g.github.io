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

integer :: status, ndata, ncol, i
real :: x 
character(len=500) :: file, record

call getarg(1,file)

open(10,file=file,status='old')
ndata = 0
do
  read(10,"( a500 )",iostat=status) record
  if(status /= 0) exit
  read(record,*,iostat=status) x, x
  if(status /= 0) cycle
  ncol = 20
  do
    read(record,*,iostat=status) (x,i=1,ncol)
    if( status /= 0 ) then
      ncol = ncol - 1
      cycle
    else
      exit
    end if
  end do
  ndata = ndata + 1
end do 
close(10)

write(*,*) ndata
write(*,*) ncol

end























