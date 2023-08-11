!
! Computes the covariance between two vectors, as a function of "time",
!          that is, as function of the shift between the two data-sets.
!
! Usage: ./covariance firstset.dat [col1] seconset.dat [col2]
!
! Where: firstset.dat is the file containing the first set of data
!        [col1] is a integer number indicating the column of 
!               the first file to be read
!        secondset.dat is the file containing the second set of data
!        [col2] is a integer number indicating the column of 
!               the second file to be read
!
! Compile with: gfortran -o covariance covariance.f90
!
! L. Martinez, Apr 2, 2008
!

program covariance

implicit none
integer :: narg, n1, n2, status, i, col1, col2, ndata, step, iargc
character :: firstchar
character(len=200) :: file1, file2, record
double precision :: avg1, avg2, norm1, norm2, cov
double precision, allocatable :: x1(:), x2(:), xcol(:)

! Read command line arguments

call getarg(1,file1)
call getarg(2,record)
read(record,*,iostat=status) col1
call getarg(3,file2)
call getarg(4,record)
read(record,*,iostat=status) col2
narg = iargc()
if( narg < 4 .or. status /= 0 ) then
  write(*,*) ' Run with: ./covariance data1.dat 1 data2.dat 3 '
  write(*,*) ' where 1 and 3 are the columns to be analysed in each set '
  stop
end if

! Write file data

write(*,"( a,a,/, a,a,/, a,i6,/, a,i6 )") &
        '# File of first data set: ', trim(file1),&
        '# File of second data set: ', trim(file2),&
        '# Column of 1 to be considered: ', col1,&
        '# Column of 2 to be considered: ', col2

! Alocate memory for xcol vector

allocate ( xcol(max0(col1,col2)) )

! Read the number of points of the first file

open(10,file=file1,action='read')
n1 = 0
do
  read(10,"( a200 )",iostat=status) record
  if ( status == 0 ) then
    if ( firstchar(record) /= '#' ) then
      read(record,*,iostat=status) (xcol(i),i=1,col1)
      if ( status == 0 ) then 
        n1 = n1 + 1
      end if
    end if
  else
    exit
  end if
end do
close(10)
write(*,"( a, i8 )") '# Number of point in first set: ', n1

! Read the number of points of the second

open(10,file=file2,action='read')
n2 = 0
do
  read(10,"( a200 )",iostat=status) record
  if ( status == 0 ) then
    if ( firstchar(record) /= '#' ) then
      read(record,*,iostat=status) (xcol(i),i=1,col2)
      if ( status == 0 ) then 
        n2 = n2 + 1
      end if
    end if
  else
    exit
  end if
end do
close(10)
write(*,"( a, i8 )") '# Number of point in second set: ', n2

! Allocate memory for data points

allocate ( x1(n1), x2(n2) )

! Reading data of set 1

open(10,file=file1,action='read')
n1 = 0
do
  read(10,"( a200 )",iostat=status) record
  if ( status == 0 ) then
    if ( firstchar(record) /= '#' ) then
      read(record,*,iostat=status) (xcol(i),i=1,col1)
      if ( status == 0 ) then    
        n1 = n1 + 1
        x1(n1) = xcol(col1)
      end if
    end if
  else
    exit
  end if
end do
close(10)

! Reading data of set 2

open(10,file=file2,action='read')
n2 = 0
do
  read(10,"( a200 )",iostat=status) record
  if ( status == 0 ) then
    if ( firstchar(record) /= '#' ) then
      read(record,*,iostat=status) (xcol(i),i=1,col2)
      if ( status == 0 ) then    
        n2 = n2 + 1
        x2(n2) = xcol(col2)
      end if
    end if
  else
    exit
  end if
end do
close(10)
 
! Actual number of points to be considered (minimum between two sets)

ndata = min0(n1,n2)
                   
! Computing the average of set 1

avg1 = 0.d0
do i = 1, ndata
  avg1 = avg1 + x1(i)
end do
avg1 = avg1 / dble(ndata)
  
! Computing the average of set 2

avg2 = 0.d0
do i = 1, ndata
  avg2 = avg2 + x2(i)
end do
avg2 = avg2 / float(ndata)

write(*,"( a,f12.5,/,a,f12.5 )") &
        '# Average of first set: ',avg1,&
        '# Average of second set: ',avg2

! Computing the time-dependent covariance between the two sets

write(*,"( a )") '# STEP  COVARIANCE*N/NORMS    COVARIANCE '
do step = 0, ndata - 1

! Computing averages in for this intervals

  avg1 = 0.d0
  do i = 1, ndata - step
    avg1 = avg1 + x1(i)
  end do
  avg1 = avg1 / float(ndata - step)
  avg2 = 0.d0
  do i = step+1, ndata
    avg2 = avg2 + x2(i)
  end do
  avg2 = avg2 / dble(ndata - step)

  norm1 = 0.d0
  norm2 = 0.d0
  cov = 0.d0
  do i = 1, ndata - step 
    norm1 = norm1 + ( x1(i) - avg1 )**2
    norm2 = norm2 + ( x2(i+step) - avg2 )**2
    cov = cov + ( x1(i) - avg1 )*( x2(i+step) - avg2 )
  end do
  norm1 = dsqrt(norm1)
  norm2 = dsqrt(norm2)
    
! Normalize covariance and output data

  if(norm1 < 1.d-8 .and. norm2 < 1.d-8) then
    write(*,"(i6,t15,f12.5,t29,f12.5)")&
            step, 1.d0, cov/dble(ndata-step) 
  else
    write(*,"(i6,t15,f12.5,t29,f12.5)")&
            step, cov/(norm1*norm2), cov/dble(ndata-step) 
  end if

end do

end

!
! Get the first non-blanck character from a line
!
 
function firstchar(record)
 
character(len=200) :: record
character :: firstchar
integer :: i
 
i = 1
do while(record(i:i) <= ' ')
  i = i + 1
end do
firstchar = record(i:i)
                                                                                                
return
end





