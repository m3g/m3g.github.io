!
! Compile with: gfortran numberofstates.f90 -o numberofstates
!
! Run with: numberofstates data.dat [xmin] [xmax] [ymin] [ymax] > output.dat
!
! This program generates the data for construction the plot of 
! Supplementary Figure S8A of 
!
! L. Martínez, A. S. Nascimento, F. M. Nunes, K. Phillips, R. Aparício, S.
! M. G. Dias, A. C. M. Figueira, J. H. Lin, P. Nguyen, J. W. Apriletti, F.
! A. R. Neves, J. D. Baxter, P. Webb, M. S. Skaf, I. Polikarpov, Gaining
! Ligand Selectivity in Thyroid Hormone Receptors via Entropy, Proc. Natl.
! Acad. Sci USA, 106, 20717-20722, 2009.
!

implicit none
integer :: ndata, i, j, icount, jcount, nstates,&
           status, le, gridsize, irnd, nrandom
real :: xread, yread, size(2)
real, allocatable :: x(:), y(:)
real :: xmin, xmax, ymin, ymax, random
character(len=200) :: file, record
logical, allocatable :: count(:,:), randomcount(:,:)

call getarg(1,file)
call getarg(2,record) ; read(record,*) xmin
call getarg(3,record) ; read(record,*) xmax
call getarg(4,record) ; read(record,*) ymin
call getarg(5,record) ; read(record,*) ymax

open(10,file=file,action='read')
ndata = 0
do 
  read(10,"( a200 )",iostat=status) record 
  if( status /= 0 ) exit
  read(record,*,iostat=status) xread, yread
  if( status /= 0 ) cycle
  ndata = ndata + 1
end do
close(10)
allocate( x(ndata), y(ndata) )

open(10,file=file,action='read')
ndata = 0
do 
  read(10,"( a200 )",iostat=status) record 
  if( status /= 0 ) exit
  read(record,*,iostat=status) xread, yread
  if( status /= 0 ) cycle
  ndata = ndata + 1
  x(ndata) = xread
  y(ndata) = yread
end do
close(10)

write(*,"( '# File: ',a )") file(1:le(file))
write(*,"( '# Number of data points: ',i8 )") ndata
write(*,"( '# Maximum and minimum grid x:',2(tr1,e12.5))") xmin, xmax
write(*,"( '# Maximum and minimum grid y:',2(tr1,e12.5))") ymin, ymax
write(*,"( '#' )")

allocate( count(1000,1000) )
allocate( randomcount(1000,1000) )

write(*,"('# Column 1: Number of bins (total number of states)')")
write(*,"('# Column 2: Number of occupied bins (number of occupied states)')")
write(*,"('# Column 3: Fraction of the states occupied')")
write(*,"('# Column 4: Number of states occupied relative to random distribution')")
write(*,"('# Column 5: Entropy relative to random distribution for T=298.15K in cal / mol K')")

do gridsize = 1, 1000, 2

  size(1) = ( xmax - xmin ) / real(gridsize)
  size(2) = ( ymax - ymin ) / real(gridsize)

  nrandom = 0
  do irnd = 1, 100
    do i = 1, gridsize
      do j = 1, gridsize
        randomcount(i,j) = .false.
      end do
    end do 
    do i = 1, ndata
      call random_number(random)
      icount = int( (random*(xmax-xmin))/size(1) ) + 1
      call random_number(random)
      jcount = int( (random*(ymax-ymin))/size(2) ) + 1
      randomcount(icount,jcount) = .true.
    end do
    do i = 1, gridsize
      do j = 1, gridsize
        if(randomcount(i,j)) nrandom = nrandom + 1
      end do
    end do
  end do
  nrandom = int(real(nrandom)/100)

  do i = 1, gridsize
    do j = 1, gridsize
      count(i,j) = .false.
    end do
  end do 

  do i = 1, ndata
    icount = int( (x(i)-xmin)/size(1) ) + 1
    jcount = int( (y(i)-ymin)/size(2) ) + 1
    count(icount,jcount) = .true.
  end do

  nstates = 0
  do i = 1, gridsize
    do j = 1, gridsize
      if(count(i,j)) nstates = nstates + 1
    end do
  end do

  write(*,"( 2(tr2,i10),4(tr2,f12.5) )") &
          gridsize, nstates,&
          real(nstates)/real(gridsize*2),&
          real(nstates)/nrandom,&
          1.987*log(real(nstates)/nrandom)

end do

deallocate( count )

end

!
! Function that sets the length of a record
!

function le(record)

integer :: le
character(len=200) :: record

le = 1
do while(record(le:le) > ' ' .and. le <= 200)
  le = le + 1
end do
le = le - 1

return
end




