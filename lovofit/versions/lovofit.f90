!
! Program LovoFit (Linear)
!
! Find the parameters c(1)...c(n) of a general linear equation of 
! the form y = c(1)*x(1) + ... + c(n-1)*x(n-1) + c(n) that
! fits the x(1), ..., x(n) data provided. The fitting
! is performed iteratively for a percentage of the data,
! chosen by the user, in order to eliminate outliers, according
! to the Low Order Value Optimization approach. The outliers
! are automatically eliminated from the fit.
!
! Uses lapack dgesv subroutine, so compile with -llapack
!
! Run with:  ./lovofit data.dat [print]
!
! [print] may be 1 or 2. If "1" a data file with the solution is provided.
!         If "2", an grace  file is provided, which can be directly visualized
!         in xmgrace. 
!
! The data.dat input file has the form:
!
! 9  135 0.8
! Y_exp   x1  x2 ... x8  COMMENT
! ...
!
! Where "9" is the number of parameters to be computed (n). 
!       "135" is the number of data points provided, 
!       '0.8" is the fraction of the data to be considered (in this example,
!             20% of the points will be considered as outliers). 
! The data rows contain the 'experimental'
! Y value and corresponding x_1 ... x_n-1 (the last constant is
! independent). "COMMENT" is an optional string that can be added to each
! row to identify the data.
!
! L. Martinez, IFSC-USP, Nov 1, 2011.
!
program lovofit

implicit none
integer :: nvar, ndata, i, j, k, info, narg, ioerror, ndata_fraction, print,&
           iter, maxiter
double precision :: prediction, residue, ss_tot, average, fraction, tolerance,&
                    delta_residue, residue_last, xmin, xmax, rand, best_residue,&
                    best_ss_tot
integer, allocatable :: ipiv(:), best_index(:)
double precision, allocatable :: experimental(:), data(:,:), a(:,:), b(:), x(:),&
                                 experimental_ordered(:), data_ordered(:,:),&
                                 error_data(:), best_b(:)
character(len=10), allocatable :: comment(:)
character(len=200) :: data_file, record
logical :: has_comments

! For flashsort
integer :: mflash
integer, allocatable :: indflash(:), lflash(:)

print = 1
narg = iargc()
if ( narg < 1 ) call error(1,0)

if ( narg > 1 ) then
  call getarg(2,data_file)
  read(data_file,*) print
end if
call getarg(1,data_file)
open(10,file=data_file,status='old',action='read',iostat=ioerror)
if ( ioerror /= 0 ) call error(2,0)

read(10,*,iostat=ioerror) nvar, ndata, fraction
if(ioerror /= 0) call error(3,0)

if(nvar < 2) then
  write(*,*) ' Error: The number of parameters to be computed cannot be smaller than 2. '
  stop
end if

! Allocate arrays

allocate( a(nvar,nvar), b(nvar), x(nvar), ipiv(nvar),& 
          experimental(ndata), data(ndata,nvar-1),&
          comment(ndata),&
          data_ordered(ndata,nvar-1), experimental_ordered(ndata),&
          error_data(ndata),&
          best_index(ndata), best_b(nvar) )

! For flashsort

allocate( indflash(ndata), lflash(ndata) )

has_comments = .false.
do i = 1, ndata
  read(10,"( a200 )") record
  read(record,*,iostat=ioerror) experimental(i), (data(i,j), j = 1, nvar-1), comment(i)
  if(ioerror == 0 ) then
    has_comments = .true.
  else
    comment(i) = ' '
    read(record,*,iostat=ioerror) experimental(i), (data(i,j), j = 1, nvar-1)
    if(ioerror /= 0) call error(4,i)
  end if
end do
close(10)
  
maxiter = 1000
best_residue = 1.d10
do iter = 1, maxiter

! Generating random initial point

  do i = 1, ndata
    call random_number(rand)
    error_data(i) = rand
  end do

  tolerance = 1.d-5
  residue = 1.d99
  delta_residue = -1.d0
  ndata_fraction = min0(ndata,int(fraction*dfloat(ndata)))
  mflash = 1 + ndata / 10
  do while( delta_residue < -1.d0*tolerance )
  
    ! Save the residue from the last iteration
  
    residue_last = residue
  
    ! Order errors from lower to higher (on first call, all errors are zero)
  
    call flashsort(error_data, ndata, lflash, mflash, indflash)
    do i = 1, ndata
      experimental_ordered(i) = experimental(indflash(i))
      do j = 1, nvar - 1
        data_ordered(i,j) = data(indflash(i),j)
      end do
    end do
  
    ! Construction of matrices A and B of linear system Ax = B
  
    do j = 1, nvar - 1
      do k = j, nvar - 1
        a(j,k) = 0.d0
        do i = 1, ndata_fraction
          a(j,k) = a(j,k) + data_ordered(i,j)*data_ordered(i,k)
        end do
      end do
    end do
    do j = 2, nvar - 1
      do k = 1, j - 1
        a(j,k) = a(k,j)
      end do
    end do
    do j = 1, nvar - 1
      a(j,nvar) = 0.d0
      do i = 1, ndata_fraction
        a(j,nvar) = a(j,nvar) + data_ordered(i,j) 
      end do
      a(nvar,j) = a(j,nvar)
    end do
    a(nvar,nvar) = dfloat(ndata_fraction)
    
    do j = 1, nvar - 1
      b(j) = 0.d0
      do i = 1, ndata_fraction
        b(j) = b(j) + experimental_ordered(i)*data_ordered(i,j)
      end do
    end do
    b(nvar) = 0.d0
    do i = 1, ndata_fraction
      b(nvar) = b(nvar) + experimental_ordered(i)
    end do
  
    ! Solving the linear system with Lapack's dgesv
    
    call dgesv(nvar,1,a,nvar,ipiv,b,nvar,info)
    if( info /= 0 ) call error(5,0)
    
    ! Computing the residue to check convergence, for the best fraction
    
    residue = 0.d0
    do i = 1, ndata_fraction
      prediction = 0.d0
      do j = 1, nvar - 1
        prediction = prediction + b(j)*data_ordered(i,j)
      end do
      prediction = prediction + b(nvar)
      error_data(indflash(i)) = ( prediction - experimental_ordered(i) )**2
      residue = residue + error_data(indflash(i))
    end do
    delta_residue = residue - residue_last
  
    ! The other residues are computed because are needed for next sorting
  
    do i = ndata_fraction + 1, ndata
      prediction = 0.d0
      do j = 1, nvar - 1
        prediction = prediction + b(j)*data_ordered(i,j)
      end do
      prediction = prediction + b(nvar)
      error_data(indflash(i)) = ( prediction - experimental_ordered(i) )**2
    end do
  
  end do

  ! Computing the R^2 and the residue for the considered fraction (lovo score)
  
  average = 0.d0
  do i = 1, ndata_fraction
    average = average + experimental_ordered(i)
  end do
  average = average / dfloat(ndata_fraction)
  residue = 0.d0
  ss_tot = 0.d0
  do i = 1, ndata_fraction 
    ss_tot = ss_tot + (experimental_ordered(i) - average)**2
    prediction = 0.d0
    do j = 1, nvar - 1
      prediction = prediction + b(j)*data_ordered(i,j)
    end do
    prediction = prediction + b(nvar)
    residue = residue + (prediction -  experimental_ordered(i))**2
  end do

  if ( ( residue - best_residue ) < -1.d-12  ) then 
    write(*,"( '# Trial, Residue, R^2 = ', i8, e17.10, e17.10 )") iter, residue, 1 - residue/ss_tot
    best_residue = residue
    best_ss_tot = ss_tot
    do i = 1, nvar
      best_b(i) = b(i)
    end do
    do i = 1, ndata
      best_index(i) = indflash(i)
    end do
  end if

end do

! Restoring best order to print solution

residue = best_residue
ss_tot = best_ss_tot
do i = 1, nvar
  b(i) = best_b(i)
end do
do i = 1, ndata
  experimental_ordered(i) = experimental(best_index(i))
  do j = 1, nvar - 1
    data_ordered(i,j) = data(best_index(i),j)
  end do
end do

! Printing the solution (predicted values) as data file

if ( print == 1 ) then
  write(*, "('# Solution: y = c(1)*x(1) + ... + c(n-1)*x(n-1) + c(n)' )")
  do i = 1, nvar
    write(*,"( '# c(',i2,') = ',f12.7 )") i, b(i)
  end do
  write(*,"( '# LOVO Fit: Fraction of data considered = ', f8.2 )") fraction
  write(*,"( '#' )")
  write(*,"( '# Solution disconsidering outliers: ' )")
  write(*,"( '# R^2 = ', f12.7, '; SD = ', f12.7 )")& 
          1 - residue/ss_tot, dsqrt(residue/dfloat(ndata-1))
  write(*,"( '#' )")
  write(*,"( '#', t11,'Data',t19,'Prediction',t36,'Residue',t48,'COMMENT',' FIT' )")
  do i = 1, ndata
    prediction = 0.d0
    do j = 1, nvar - 1
      prediction = prediction + b(j)*data_ordered(i,j)
    end do
    prediction = prediction + b(nvar)
    if ( i <= ndata_fraction ) then
      k = 1
    else
      k = 0
    end if
    write(*,"( 3(tr2,f12.7),tr2,a10,tr2,i2 )") experimental_ordered(i), prediction,&
                                               abs(prediction - experimental_ordered(i)),& 
                                               trim(comment(best_index(i))), k
  end do
end if

! Printing the solution as xmgrace plot

if ( print == 2 ) then
  write(*,"( '# Standard Deviation = ', f12.5 )") dsqrt(residue/dfloat(ndata-1)) 
  call printxmgr(1)
  do i = 1, ndata 
    prediction = 0.d0
    do j = 1, nvar - 1
      prediction = prediction + b(j)*data_ordered(i,j)
    end do
    prediction = prediction + b(nvar)
    if ( i <= ndata_fraction ) then
      k = 1
    else
      k = 0
    end if
    write(*,*) experimental_ordered(i), prediction, k
  end do
  call printxmgr(2)
  xmin = 1.d10
  xmax = -1.d10
  do i = 1, ndata
    xmin = dmin1(xmin,experimental_ordered(i))
    xmax = dmax1(xmax,experimental_ordered(i))
  end do
  write(*,*) xmin, xmin
  write(*,*) xmax, xmax
  call printxmgr(3)
end if

end

!
! Subroutine for printing error messages
!

subroutine error(i,j)
integer :: i, j
select case(i) 
  case(1)
    write(*,*) ' Run with: ./lovofit data.dat '
  case(2) 
    write(*,*) ' ERROR: Could not open data file. '
  case(3) 
    write(*,*) ' ERROR: Could not read nvar, ndata and frac from first line of data file. '
               '        Should be something like: 2 300 0.7 '
  case(4)
    write(*,*) ' ERROR: Could not read data file, at data point: ', j
  case(5) 
    write(*,*) ' ERROR: Gaussian elimination failed. '
end select
stop

end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                             c
!     Subroutine Flashsort
!     SORTS ARRAY A WITH N ELEMENTS BY USE OF INDEX VECTOR L  c
!     OF DIMENSION M WITH M ABOUT 0.1 N.                      c
!     Karl-Dietrich Neubert, FlashSort1 Algorithm             c
!     in  Dr. Dobb's Journal Feb.1998,p.123                   c
!                                                             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine flashsort(A, N, L, M, ind)

      double precision a(1), anmin, c1, hold, flash
      integer L(1), ind(1), i, n, nmax, m, k, ihold, nmove, j, iflash
!     ============================ CLASS FORMATION ===== 


      do i = 1, n
      ind(i) = i
      end do

      ANMIN=A(1)
      NMAX=1 
      DO I=1,N
         IF( A(I).LT.ANMIN) ANMIN=A(I)
         IF( A(I).GT.A(NMAX)) NMAX=I
      END DO

      IF (ANMIN.EQ.A(NMAX)) RETURN
      C1=(M - 1) / (A(NMAX) - ANMIN)
      DO K=1,M  
         L(K)=0
      END DO 
      DO I=1,N
         K=1 + INT(C1 * (A(I) - ANMIN))
         L(K)=L(K) + 1
      END DO
      DO K=2,M
         L(K)=L(K) + L(K - 1)
      END DO
      HOLD=A(NMAX)
      A(NMAX)=A(1) 
      A(1)=HOLD

      ihold = ind(nmax)
      ind(nmax) = ind(1)
      ind(1) = ihold


!     =============================== PERMUTATION ===== 
      NMOVE=0 
      J=1
      K=M
      DO WHILE (NMOVE.LT.N - 1)
         DO WHILE (J.GT.L(K)) 
            J=J + 1 
            K=1 + INT(C1 * (A(J) - ANMIN)) 
         END DO  
         FLASH=A(J)
         iflash=ind(j)

         DO WHILE (.NOT.(J.EQ.L(K) + 1)) 
            K=1 + INT(C1 * (FLASH - ANMIN))
            HOLD=A(L(K)) 
            ihold = ind(L(k))
            A(L(K))=FLASH
            ind(L(k)) = iflash
            iflash = ihold
            FLASH=HOLD
            L(K)=L(K) - 1
            NMOVE=NMOVE + 1 
         END DO
      END DO

!     ========================= STRAIGHT INSERTION =====
      DO I=N-2,1,-1
         IF  (A(I + 1).LT.A(I)) THEN
            HOLD=A(I)
            ihold = ind(i)
            J=I
            DO WHILE  (A(J + 1).LT.HOLD)
               A(J)=A(J + 1)
               ind(j) = ind(j+1) 
               J=J + 1 
            END DO
            A(J)=HOLD 
            ind(j) = ihold
         ENDIF
      END DO

!     =========================== RETURN,END FLASH1 =====
      RETURN
      END                               

! Subroutine that prints the xmgrace plot information

subroutine printxmgr(i)

integer :: i

if ( i == 1 ) then 
  write(*,*) '#                                                                                                '
  write(*,*) '# Grace project file                                                                             '
  write(*,*) '#                                                                                                '
  write(*,*) '@version 50122                                                                                   '
  write(*,*) '@page size 594, 457                                                                              '
  write(*,*) '@page scroll 5%                                                                                  '
  write(*,*) '@page inout 5%                                                                                   '
  write(*,*) '@link page off                                                                                   '
  write(*,*) '@map font 8 to "Courier", "Courier"                                                          '
  write(*,*) '@map font 10 to "Courier-Bold", "Courier-Bold"                                               '
  write(*,*) '@map font 11 to "Courier-BoldOblique", "Courier-BoldOblique"                                 '
  write(*,*) '@map font 9 to "Courier-Oblique", "Courier-Oblique"                                          '
  write(*,*) '@map font 4 to "Helvetica", "Helvetica"                                                      '
  write(*,*) '@map font 6 to "Helvetica-Bold", "Helvetica-Bold"                                            '
  write(*,*) '@map font 7 to "Helvetica-BoldOblique", "Helvetica-BoldOblique"                              '
  write(*,*) '@map font 15 to "Helvetica-Narrow", "Helvetica-Narrow"                                       '
  write(*,*) '@map font 16 to "Helvetica-Narrow-Bold", "Helvetica-Narrow-Bold"                             '
  write(*,*) '@map font 17 to "Helvetica-Narrow-BoldOblique", "Helvetica-Narrow-BoldOblique"               '
  write(*,*) '@map font 18 to "Helvetica-Narrow-Oblique", "Helvetica-Narrow-Oblique"                       '
  write(*,*) '@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique"                                      '
  write(*,*) '@map font 20 to "NewCenturySchlbk-Bold", "NewCenturySchlbk-Bold"                             '
  write(*,*) '@map font 21 to "NewCenturySchlbk-BoldItalic", "NewCenturySchlbk-BoldItalic"                 '
  write(*,*) '@map font 22 to "NewCenturySchlbk-Italic", "NewCenturySchlbk-Italic"                         '
  write(*,*) '@map font 23 to "NewCenturySchlbk-Roman", "NewCenturySchlbk-Roman"                           '
  write(*,*) '@map font 24 to "Palatino-Bold", "Palatino-Bold"                                             '
  write(*,*) '@map font 25 to "Palatino-BoldItalic", "Palatino-BoldItalic"                                 '
  write(*,*) '@map font 26 to "Palatino-Italic", "Palatino-Italic"                                         '
  write(*,*) '@map font 27 to "Palatino-Roman", "Palatino-Roman"                                           '
  write(*,*) '@map font 12 to "Symbol", "Symbol"                                                           '
  write(*,*) '@map font 2 to "Times-Bold", "Times-Bold"                                                    '
  write(*,*) '@map font 3 to "Times-BoldItalic", "Times-BoldItalic"                                        '
  write(*,*) '@map font 1 to "Times-Italic", "Times-Italic"                                                '
  write(*,*) '@map font 0 to "Times-Roman", "Times-Roman"                                                  '
  write(*,*) '@map font 33 to "ZapfChancery-MediumItalic", "ZapfChancery-MediumItalic"                     '
  write(*,*) '@map font 13 to "ZapfDingbats", "ZapfDingbats"                                               '
  write(*,*) '@map font 35 to "LMRoman10-Bold", "LMRoman10-Bold"                                           '
  write(*,*) '@map font 36 to "LMRoman10-BoldItalic", "LMRoman10-BoldItalic"                               '
  write(*,*) '@map font 37 to "LMRoman10-BoldOblique", "LMRoman10-BoldOblique"                             '
  write(*,*) '@map font 38 to "LMRoman10-CapsOblique", "LMRoman10-CapsOblique"                             '
  write(*,*) '@map font 39 to "LMRoman10-CapsRegular", "LMRoman10-CapsRegular"                             '
  write(*,*) '@map font 40 to "LMRoman10-Demi", "LMRoman10-Demi"                                           '
  write(*,*) '@map font 41 to "LMRoman10-DemiOblique", "LMRoman10-DemiOblique"                             '
  write(*,*) '@map font 42 to "LMRoman10-Dunhill", "LMRoman10-Dunhill"                                     '
  write(*,*) '@map font 43 to "LMRoman10-DunhillOblique", "LMRoman10-DunhillOblique"                       '
  write(*,*) '@map font 44 to "LMRoman10-Italic", "LMRoman10-Italic"                                       '
  write(*,*) '@map font 45 to "LMRoman10-Oblique", "LMRoman10-Oblique"                                     '
  write(*,*) '@map font 46 to "LMRoman10-Regular", "LMRoman10-Regular"                                     '
  write(*,*) '@map font 47 to "LMRoman10-Unslanted", "LMRoman10-Unslanted"                                 '
  write(*,*) '@map font 48 to "LMSans10-Bold", "LMSans10-Bold"                                             '
  write(*,*) '@map font 49 to "LMSans10-BoldOblique", "LMSans10-BoldOblique"                               '
  write(*,*) '@map font 50 to "LMSans10-DemiCondensed", "LMSans10-DemiCondensed"                           '
  write(*,*) '@map font 51 to "LMSans10-DemiCondensedOblique", "LMSans10-DemiCondensedOblique"             '
  write(*,*) '@map font 52 to "LMSans10-Oblique", "LMSans10-Oblique"                                       '
  write(*,*) '@map font 53 to "LMSans10-Regular", "LMSans10-Regular"                                       '
  write(*,*) '@map font 54 to "LMSansQuotation8-Bold", "LMSansQuotation8-Bold"                             '
  write(*,*) '@map font 55 to "LMSansQuotation8-BoldOblique", "LMSansQuotation8-BoldOblique"               '
  write(*,*) '@map font 56 to "LMSansQuotation8-Oblique", "LMSansQuotation8-Oblique"                       '
  write(*,*) '@map font 57 to "LMSansQuotation8-Regular", "LMSansQuotation8-Regular"                       '
  write(*,*) '@map font 58 to "LMTypewriter10-CapsOblique", "LMTypewriter10-CapsOblique"                   '
  write(*,*) '@map font 59 to "LMTypewriter10-CapsRegular", "LMTypewriter10-CapsRegular"                   '
  write(*,*) '@map font 60 to "LMTypewriter10-Dark", "LMTypewriter10-Dark"                                 '
  write(*,*) '@map font 61 to "LMTypewriter10-DarkOblique", "LMTypewriter10-DarkOblique"                   '
  write(*,*) '@map font 62 to "LMTypewriter10-Italic", "LMTypewriter10-Italic"                             '
  write(*,*) '@map font 63 to "LMTypewriter10-Light", "LMTypewriter10-Light"                               '
  write(*,*) '@map font 64 to "LMTypewriter10-LightCondensed", "LMTypewriter10-LightCondensed"             '
  write(*,*) '@map font 65 to "LMTypewriter10-LightCondensedOblique", "LMTypewriter10-LightCondensedOblique"'
  write(*,*) '@map font 66 to "LMTypewriter10-LightOblique", "LMTypewriter10-LightOblique"                 '
  write(*,*) '@map font 67 to "LMTypewriter10-Oblique", "LMTypewriter10-Oblique"                           '
  write(*,*) '@map font 68 to "LMTypewriter10-Regular", "LMTypewriter10-Regular"                           '
  write(*,*) '@map font 69 to "LMTypewriterVarWd10-Dark", "LMTypewriterVarWd10-Dark"                       '
  write(*,*) '@map font 70 to "LMTypewriterVarWd10-DarkOblique", "LMTypewriterVarWd10-DarkOblique"         '
  write(*,*) '@map font 71 to "LMTypewriterVarWd10-Light", "LMTypewriterVarWd10-Light"                     '
  write(*,*) '@map font 72 to "LMTypewriterVarWd10-LightOblique", "LMTypewriterVarWd10-LightOblique"       '
  write(*,*) '@map font 73 to "LMTypewriterVarWd10-Oblique", "LMTypewriterVarWd10-Oblique"                 '
  write(*,*) '@map font 74 to "LMTypewriterVarWd10-Regular", "LMTypewriterVarWd10-Regular"                 '
  write(*,*) '@map color 0 to (255, 255, 255), "white"                                                       '
  write(*,*) '@map color 1 to (0, 0, 0), "black"                                                             '
  write(*,*) '@map color 2 to (255, 0, 0), "red"                                                             '
  write(*,*) '@map color 3 to (0, 255, 0), "green"                                                           '
  write(*,*) '@map color 4 to (0, 0, 255), "blue"                                                            '
  write(*,*) '@map color 5 to (255, 255, 0), "yellow"                                                        '
  write(*,*) '@map color 6 to (188, 143, 143), "brown"                                                       '
  write(*,*) '@map color 7 to (220, 220, 220), "grey"                                                        '
  write(*,*) '@map color 8 to (148, 0, 211), "violet"                                                        '
  write(*,*) '@map color 9 to (0, 255, 255), "cyan"                                                          '
  write(*,*) '@map color 10 to (255, 0, 255), "magenta"                                                      '
  write(*,*) '@map color 11 to (255, 165, 0), "orange"                                                       '
  write(*,*) '@map color 12 to (114, 33, 188), "indigo"                                                      '
  write(*,*) '@map color 13 to (103, 7, 72), "maroon"                                                        '
  write(*,*) '@map color 14 to (64, 224, 208), "turquoise"                                                   '
  write(*,*) '@map color 15 to (0, 139, 0), "green4"                                                         '
  write(*,*) '@map color 50 to (0, 0, 0), "grey"                                                             '
  write(*,*) '@map color 51 to (10, 10, 10), "grey"                                                          '
  write(*,*) '@map color 52 to (20, 20, 20), "grey"                                                          '
  write(*,*) '@map color 53 to (30, 30, 30), "grey"                                                          '
  write(*,*) '@map color 54 to (40, 40, 40), "grey"                                                          '
  write(*,*) '@map color 55 to (50, 50, 50), "grey"                                                          '
  write(*,*) '@map color 56 to (60, 60, 60), "grey"                                                          '
  write(*,*) '@map color 57 to (70, 70, 70), "grey"                                                          '
  write(*,*) '@map color 58 to (80, 80, 80), "grey"                                                          '
  write(*,*) '@map color 59 to (90, 90, 90), "grey"                                                          '
  write(*,*) '@map color 60 to (100, 100, 100), "grey"                                                       '
  write(*,*) '@map color 61 to (110, 110, 110), "grey"                                                       '
  write(*,*) '@map color 62 to (120, 120, 120), "grey"                                                       '
  write(*,*) '@map color 63 to (130, 130, 130), "grey"                                                       '
  write(*,*) '@map color 64 to (140, 140, 140), "grey"                                                       '
  write(*,*) '@map color 65 to (150, 150, 150), "grey"                                                       '
  write(*,*) '@map color 66 to (160, 160, 160), "grey"                                                       '
  write(*,*) '@map color 67 to (170, 170, 170), "grey"                                                       '
  write(*,*) '@map color 68 to (180, 180, 180), "grey"                                                       '
  write(*,*) '@map color 69 to (190, 190, 190), "grey"                                                       '
  write(*,*) '@map color 70 to (200, 200, 200), "grey"                                                       '
  write(*,*) '@map color 71 to (210, 210, 210), "grey"                                                       '
  write(*,*) '@map color 72 to (220, 220, 220), "grey"                                                       '
  write(*,*) '@map color 73 to (230, 230, 230), "grey"                                                       '
  write(*,*) '@map color 74 to (240, 240, 240), "grey"                                                       '
  write(*,*) '@map color 75 to (250, 250, 250), "grey"                                                       '
  write(*,*) '@map color 76 to (255, 255, 255), "grey"                                                       '
  write(*,*) '@map color 100 to (0, 0, 255), "heat"                                                          '
  write(*,*) '@map color 101 to (10, 3, 234), "heat"                                                         '
  write(*,*) '@map color 102 to (20, 11, 213), "heat"                                                        '
  write(*,*) '@map color 103 to (30, 24, 192), "heat"                                                        '
  write(*,*) '@map color 104 to (40, 66, 181), "heat"                                                        '
  write(*,*) '@map color 105 to (50, 108, 170), "heat"                                                       '
  write(*,*) '@map color 106 to (60, 150, 160), "heat"                                                       '
  write(*,*) '@map color 107 to (70, 192, 150), "heat"                                                       '
  write(*,*) '@map color 108 to (80, 210, 140), "heat"                                                       '
  write(*,*) '@map color 109 to (90, 220, 130), "heat"                                                       '
  write(*,*) '@map color 110 to (100, 230, 120), "heat"                                                      '
  write(*,*) '@map color 111 to (110, 240, 110), "heat"                                                      '
  write(*,*) '@map color 112 to (120, 255, 100), "heat"                                                      '
  write(*,*) '@map color 113 to (130, 255, 90), "heat"                                                       '
  write(*,*) '@map color 114 to (140, 255, 80), "heat"                                                       '
  write(*,*) '@map color 115 to (150, 255, 70), "heat"                                                       '
  write(*,*) '@map color 116 to (170, 255, 60), "heat"                                                       '
  write(*,*) '@map color 117 to (190, 255, 50), "heat"                                                       '
  write(*,*) '@map color 118 to (210, 255, 40), "heat"                                                       '
  write(*,*) '@map color 119 to (230, 240, 30), "heat"                                                       '
  write(*,*) '@map color 120 to (255, 220, 20), "heat"                                                       '
  write(*,*) '@map color 121 to (255, 180, 10), "heat"                                                       '
  write(*,*) '@map color 122 to (255, 140, 5), "heat"                                                        '
  write(*,*) '@map color 123 to (255, 100, 0), "heat"                                                        '
  write(*,*) '@map color 124 to (255, 60, 0), "heat"                                                         '
  write(*,*) '@map color 125 to (255, 20, 0), "heat"                                                         '
  write(*,*) '@map color 126 to (255, 0, 0), "heat"                                                          '
  write(*,*) '@map color 200 to (0, 0, 255), "h&c"                                                           '
  write(*,*) '@map color 201 to (10, 25, 255), "h&c"                                                         '
  write(*,*) '@map color 202 to (20, 50, 255), "h&c"                                                         '
  write(*,*) '@map color 203 to (30, 75, 255), "h&c"                                                         '
  write(*,*) '@map color 204 to (50, 100, 255), "h&c"                                                        '
  write(*,*) '@map color 205 to (70, 125, 255), "h&c"                                                        '
  write(*,*) '@map color 206 to (90, 150, 255), "h&c"                                                        '
  write(*,*) '@map color 207 to (110, 175, 255), "h&c"                                                       '
  write(*,*) '@map color 208 to (130, 200, 255), "h&c"                                                       '
  write(*,*) '@map color 209 to (160, 225, 255), "h&c"                                                       '
  write(*,*) '@map color 210 to (190, 250, 255), "h&c"                                                       '
  write(*,*) '@map color 211 to (220, 255, 255), "h&c"                                                       '
  write(*,*) '@map color 212 to (250, 255, 255), "h&c"                                                       '
  write(*,*) '@map color 213 to (255, 255, 255), "h&c"                                                       '
  write(*,*) '@map color 214 to (255, 255, 220), "h&c"                                                       '
  write(*,*) '@map color 215 to (255, 255, 190), "h&c"                                                       '
  write(*,*) '@map color 216 to (255, 250, 160), "h&c"                                                       '
  write(*,*) '@map color 217 to (255, 225, 130), "h&c"                                                       '
  write(*,*) '@map color 218 to (255, 200, 110), "h&c"                                                       '
  write(*,*) '@map color 219 to (255, 175, 90), "h&c"                                                        '
  write(*,*) '@map color 220 to (255, 150, 70), "h&c"                                                        '
  write(*,*) '@map color 221 to (255, 125, 50), "h&c"                                                        '
  write(*,*) '@map color 222 to (255, 100, 40), "h&c"                                                        '
  write(*,*) '@map color 223 to (255, 75, 30), "h&c"                                                         '
  write(*,*) '@map color 224 to (255, 50, 20), "h&c"                                                         '
  write(*,*) '@map color 225 to (255, 25, 10), "h&c"                                                         '
  write(*,*) '@map color 226 to (255, 0, 0), "h&c"                                                           '
  write(*,*) '@reference date 0                                                                                '
  write(*,*) '@date wrap off                                                                                   '
  write(*,*) '@date wrap year 1950                                                                             '
  write(*,*) '@default linewidth 1.0                                                                           '
  write(*,*) '@default linestyle 1                                                                             '
  write(*,*) '@default color 1                                                                                 '
  write(*,*) '@default pattern 1                                                                               '
  write(*,*) '@default font 0                                                                                  '
  write(*,*) '@default char size 1.000000                                                                      '
  write(*,*) '@default symbol size 1.000000                                                                    '
  write(*,*) '@default sformat "%.8g"                                                                        '
  write(*,*) '@background color 0                                                                              '
  write(*,*) '@page background fill on                                                                         '
  write(*,*) '@timestamp off                                                                                   '
  write(*,*) '@timestamp 0.03, 0.03                                                                            '
  write(*,*) '@timestamp color 1                                                                               '
  write(*,*) '@timestamp rot 0                                                                                 '
  write(*,*) '@timestamp font 0                                                                                '
  write(*,*) '@timestamp char size 1.000000                                                                    '
  write(*,*) '@timestamp def "Mon Dec 12 15:14:09 2011"                                                      '
  write(*,*) '@r0 off                                                                                          '
  write(*,*) '@link r0 to g0                                                                                   '
  write(*,*) '@r0 type above                                                                                   '
  write(*,*) '@r0 linestyle 1                                                                                  '
  write(*,*) '@r0 linewidth 1.0                                                                                '
  write(*,*) '@r0 color 1                                                                                      '
  write(*,*) '@r0 line 0, 0, 0, 0                                                                              '
  write(*,*) '@r1 off                                                                                          '
  write(*,*) '@link r1 to g0                                                                                   '
  write(*,*) '@r1 type above                                                                                   '
  write(*,*) '@r1 linestyle 1                                                                                  '
  write(*,*) '@r1 linewidth 1.0                                                                                '
  write(*,*) '@r1 color 1                                                                                      '
  write(*,*) '@r1 line 0, 0, 0, 0                                                                              '
  write(*,*) '@r2 off                                                                                          '
  write(*,*) '@link r2 to g0                                                                                   '
  write(*,*) '@r2 type above                                                                                   '
  write(*,*) '@r2 linestyle 1                                                                                  '
  write(*,*) '@r2 linewidth 1.0                                                                                '
  write(*,*) '@r2 color 1                                                                                      '
  write(*,*) '@r2 line 0, 0, 0, 0                                                                              '
  write(*,*) '@r3 off                                                                                          '
  write(*,*) '@link r3 to g0                                                                                   '
  write(*,*) '@r3 type above                                                                                   '
  write(*,*) '@r3 linestyle 1                                                                                  '
  write(*,*) '@r3 linewidth 1.0                                                                                '
  write(*,*) '@r3 color 1                                                                                      '
  write(*,*) '@r3 line 0, 0, 0, 0                                                                              '
  write(*,*) '@r4 off                                                                                          '
  write(*,*) '@link r4 to g0                                                                                   '
  write(*,*) '@r4 type above                                                                                   '
  write(*,*) '@r4 linestyle 1                                                                                  '
  write(*,*) '@r4 linewidth 1.0                                                                                '
  write(*,*) '@r4 color 1                                                                                      '
  write(*,*) '@r4 line 0, 0, 0, 0                                                                              '
  write(*,*) '@g0 on                                                                                           '
  write(*,*) '@g0 hidden false                                                                                 '
  write(*,*) '@g0 type XY                                                                                      '
  write(*,*) '@g0 stacked false                                                                                '
  write(*,*) '@g0 bar hgap 0.000000                                                                            '
  write(*,*) '@g0 fixedpoint off                                                                               '
  write(*,*) '@g0 fixedpoint type 0                                                                            '
  write(*,*) '@g0 fixedpoint xy 0.000000, 0.000000                                                             '
  write(*,*) '@g0 fixedpoint format general general                                                            '
  write(*,*) '@g0 fixedpoint prec 6, 6                                                                         '
  write(*,*) '@with g0                                                                                         '
  write(*,*) '@    world -40, -40, 7, 7                                                                        '
  write(*,*) '@    stack world 0, 0, 0, 0                                                                      '
  write(*,*) '@    znorm 1                                                                                     '
  write(*,*) '@    view 0.150656, 0.150000, 1.155033, 0.850000                                                 '
  write(*,*) '@    title ""                                                                                  '
  write(*,*) '@    title font 0                                                                                '
  write(*,*) '@    title size 1.500000                                                                         '
  write(*,*) '@    title color 1                                                                               '
  write(*,*) '@    subtitle ""                                                                               '
  write(*,*) '@    subtitle font 0                                                                             '
  write(*,*) '@    subtitle size 1.000000                                                                      '
  write(*,*) '@    subtitle color 1                                                                            '
  write(*,*) '@    xaxes scale Normal                                                                          '
  write(*,*) '@    yaxes scale Normal                                                                          '
  write(*,*) '@    xaxes invert off                                                                            '
  write(*,*) '@    yaxes invert off                                                                            '
  write(*,*) '@    xaxis  on                                                                                   '
  write(*,*) '@    xaxis  type zero false                                                                      '
  write(*,*) '@    xaxis  offset 0.000000 , 0.000000                                                           '
  write(*,*) '@    xaxis  bar on                                                                               '
  write(*,*) '@    xaxis  bar color 1                                                                          '
  write(*,*) '@    xaxis  bar linestyle 1                                                                      '
  write(*,*) '@    xaxis  bar linewidth 1.0                                                                    '
  write(*,*) '@    xaxis  label ""                                                                           '
  write(*,*) '@    xaxis  label layout para                                                                    '
  write(*,*) '@    xaxis  label place auto                                                                     '
  write(*,*) '@    xaxis  label char size 1.000000                                                             '
  write(*,*) '@    xaxis  label font 0                                                                         '
  write(*,*) '@    xaxis  label color 1                                                                        '
  write(*,*) '@    xaxis  label place normal                                                                   '
  write(*,*) '@    xaxis  tick on                                                                              '
  write(*,*) '@    xaxis  tick major 10                                                                        '
  write(*,*) '@    xaxis  tick minor ticks 1                                                                   '
  write(*,*) '@    xaxis  tick default 6                                                                       '
  write(*,*) '@    xaxis  tick place rounded true                                                              '
  write(*,*) '@    xaxis  tick in                                                                              '
  write(*,*) '@    xaxis  tick major size 1.000000                                                             '
  write(*,*) '@    xaxis  tick major color 1                                                                   '
  write(*,*) '@    xaxis  tick major linewidth 1.0                                                             '
  write(*,*) '@    xaxis  tick major linestyle 1                                                               '
  write(*,*) '@    xaxis  tick major grid off                                                                  '
  write(*,*) '@    xaxis  tick minor color 1                                                                   '
  write(*,*) '@    xaxis  tick minor linewidth 1.0                                                             '
  write(*,*) '@    xaxis  tick minor linestyle 1                                                               '
  write(*,*) '@    xaxis  tick minor grid off                                                                  '
  write(*,*) '@    xaxis  tick minor size 0.500000                                                             '
  write(*,*) '@    xaxis  ticklabel on                                                                         '
  write(*,*) '@    xaxis  ticklabel format general                                                             '
  write(*,*) '@    xaxis  ticklabel prec 5                                                                     '
  write(*,*) '@    xaxis  ticklabel formula ""                                                               '
  write(*,*) '@    xaxis  ticklabel append ""                                                                '
  write(*,*) '@    xaxis  ticklabel prepend ""                                                               '
  write(*,*) '@    xaxis  ticklabel angle 0                                                                    '
  write(*,*) '@    xaxis  ticklabel skip 0                                                                     '
  write(*,*) '@    xaxis  ticklabel stagger 0                                                                  '
  write(*,*) '@    xaxis  ticklabel place normal                                                               '
  write(*,*) '@    xaxis  ticklabel offset auto                                                                '
  write(*,*) '@    xaxis  ticklabel offset 0.000000 , 0.010000                                                 '
  write(*,*) '@    xaxis  ticklabel start type auto                                                            '
  write(*,*) '@    xaxis  ticklabel start 0.000000                                                             '
  write(*,*) '@    xaxis  ticklabel stop type auto                                                             '
  write(*,*) '@    xaxis  ticklabel stop 0.000000                                                              '
  write(*,*) '@    xaxis  ticklabel char size 1.000000                                                         '
  write(*,*) '@    xaxis  ticklabel font 0                                                                     '
  write(*,*) '@    xaxis  ticklabel color 1                                                                    '
  write(*,*) '@    xaxis  tick place both                                                                      '
  write(*,*) '@    xaxis  tick spec type none                                                                  '
  write(*,*) '@    yaxis  on                                                                                   '
  write(*,*) '@    yaxis  type zero false                                                                      '
  write(*,*) '@    yaxis  offset 0.000000 , 0.000000                                                           '
  write(*,*) '@    yaxis  bar on                                                                               '
  write(*,*) '@    yaxis  bar color 1                                                                          '
  write(*,*) '@    yaxis  bar linestyle 1                                                                      '
  write(*,*) '@    yaxis  bar linewidth 1.0                                                                    '
  write(*,*) '@    yaxis  label ""                                                                           '
  write(*,*) '@    yaxis  label layout para                                                                    '
  write(*,*) '@    yaxis  label place auto                                                                     '
  write(*,*) '@    yaxis  label char size 1.000000                                                             '
  write(*,*) '@    yaxis  label font 0                                                                         '
  write(*,*) '@    yaxis  label color 1                                                                        '
  write(*,*) '@    yaxis  label place normal                                                                   '
  write(*,*) '@    yaxis  tick on                                                                              '
  write(*,*) '@    yaxis  tick major 10                                                                        '
  write(*,*) '@    yaxis  tick minor ticks 1                                                                   '
  write(*,*) '@    yaxis  tick default 6                                                                       '
  write(*,*) '@    yaxis  tick place rounded true                                                              '
  write(*,*) '@    yaxis  tick in                                                                              '
  write(*,*) '@    yaxis  tick major size 1.000000                                                             '
  write(*,*) '@    yaxis  tick major color 1                                                                   '
  write(*,*) '@    yaxis  tick major linewidth 1.0                                                             '
  write(*,*) '@    yaxis  tick major linestyle 1                                                               '
  write(*,*) '@    yaxis  tick major grid off                                                                  '
  write(*,*) '@    yaxis  tick minor color 1                                                                   '
  write(*,*) '@    yaxis  tick minor linewidth 1.0                                                             '
  write(*,*) '@    yaxis  tick minor linestyle 1                                                               '
  write(*,*) '@    yaxis  tick minor grid off                                                                  '
  write(*,*) '@    yaxis  tick minor size 0.500000                                                             '
  write(*,*) '@    yaxis  ticklabel on                                                                         '
  write(*,*) '@    yaxis  ticklabel format general                                                             '
  write(*,*) '@    yaxis  ticklabel prec 5                                                                     '
  write(*,*) '@    yaxis  ticklabel formula ""                                                               '
  write(*,*) '@    yaxis  ticklabel append ""                                                                '
  write(*,*) '@    yaxis  ticklabel prepend ""                                                               '
  write(*,*) '@    yaxis  ticklabel angle 0                                                                    '
  write(*,*) '@    yaxis  ticklabel skip 0                                                                     '
  write(*,*) '@    yaxis  ticklabel stagger 0                                                                  '
  write(*,*) '@    yaxis  ticklabel place normal                                                               '
  write(*,*) '@    yaxis  ticklabel offset auto                                                                '
  write(*,*) '@    yaxis  ticklabel offset 0.000000 , 0.010000                                                 '
  write(*,*) '@    yaxis  ticklabel start type auto                                                            '
  write(*,*) '@    yaxis  ticklabel start 0.000000                                                             '
  write(*,*) '@    yaxis  ticklabel stop type auto                                                             '
  write(*,*) '@    yaxis  ticklabel stop 0.000000                                                              '
  write(*,*) '@    yaxis  ticklabel char size 1.000000                                                         '
  write(*,*) '@    yaxis  ticklabel font 0                                                                     '
  write(*,*) '@    yaxis  ticklabel color 1                                                                    '
  write(*,*) '@    yaxis  tick place both                                                                      '
  write(*,*) '@    yaxis  tick spec type none                                                                  '
  write(*,*) '@    altxaxis  off                                                                               '
  write(*,*) '@    altyaxis  off                                                                               '
  write(*,*) '@    legend on                                                                                   '
  write(*,*) '@    legend loctype view                                                                         '
  write(*,*) '@    legend 0.853719912473, 0.8                                                                  '
  write(*,*) '@    legend box color 1                                                                          '
  write(*,*) '@    legend box pattern 1                                                                        '
  write(*,*) '@    legend box linewidth 1.0                                                                    '
  write(*,*) '@    legend box linestyle 1                                                                      '
  write(*,*) '@    legend box fill color 0                                                                     '
  write(*,*) '@    legend box fill pattern 1                                                                   '
  write(*,*) '@    legend font 0                                                                               '
  write(*,*) '@    legend char size 1.000000                                                                   '
  write(*,*) '@    legend color 1                                                                              '
  write(*,*) '@    legend length 4                                                                             '
  write(*,*) '@    legend vgap 1                                                                               '
  write(*,*) '@    legend hgap 1                                                                               '
  write(*,*) '@    legend invert false                                                                         '
  write(*,*) '@    frame type 0                                                                                '
  write(*,*) '@    frame linestyle 1                                                                           '
  write(*,*) '@    frame linewidth 1.0                                                                         '
  write(*,*) '@    frame color 1                                                                               '
  write(*,*) '@    frame pattern 1                                                                             '
  write(*,*) '@    frame background color 0                                                                    '
  write(*,*) '@    frame background pattern 0                                                                  '
  write(*,*) '@    s0 hidden false                                                                             '
  write(*,*) '@    s0 type xycolor                                                                             '
  write(*,*) '@    s0 symbol 1                                                                                 '
  write(*,*) '@    s0 symbol size 0.410000                                                                     '
  write(*,*) '@    s0 symbol color 1                                                                           '
  write(*,*) '@    s0 symbol pattern 1                                                                         '
  write(*,*) '@    s0 symbol fill color 1                                                                      '
  write(*,*) '@    s0 symbol fill pattern 1                                                                    '
  write(*,*) '@    s0 symbol linewidth 1.0                                                                     '
  write(*,*) '@    s0 symbol linestyle 1                                                                       '
  write(*,*) '@    s0 symbol char 65                                                                           '
  write(*,*) '@    s0 symbol char font 0                                                                       '
  write(*,*) '@    s0 symbol skip 0                                                                            '
  write(*,*) '@    s0 line type 0                                                                              '
  write(*,*) '@    s0 line linestyle 1                                                                         '
  write(*,*) '@    s0 line linewidth 1.0                                                                       '
  write(*,*) '@    s0 line color 1                                                                             '
  write(*,*) '@    s0 line pattern 1                                                                           '
  write(*,*) '@    s0 baseline type 0                                                                          '
  write(*,*) '@    s0 baseline off                                                                             '
  write(*,*) '@    s0 dropline off                                                                             '
  write(*,*) '@    s0 fill type 0                                                                              '
  write(*,*) '@    s0 fill rule 0                                                                              '
  write(*,*) '@    s0 fill color 1                                                                             '
  write(*,*) '@    s0 fill pattern 1                                                                           '
  write(*,*) '@    s0 avalue off                                                                               '
  write(*,*) '@    s0 avalue type 2                                                                            '
  write(*,*) '@    s0 avalue char size 1.000000                                                                '
  write(*,*) '@    s0 avalue font 0                                                                            '
  write(*,*) '@    s0 avalue color 1                                                                           '
  write(*,*) '@    s0 avalue rot 0                                                                             '
  write(*,*) '@    s0 avalue format general                                                                    '
  write(*,*) '@    s0 avalue prec 3                                                                            '
  write(*,*) '@    s0 avalue prepend ""                                                                      '
  write(*,*) '@    s0 avalue append ""                                                                       '
  write(*,*) '@    s0 avalue offset 0.000000 , 0.000000                                                        '
  write(*,*) '@    s0 errorbar on                                                                              '
  write(*,*) '@    s0 errorbar place both                                                                      '
  write(*,*) '@    s0 errorbar color 1                                                                         '
  write(*,*) '@    s0 errorbar pattern 1                                                                       '
  write(*,*) '@    s0 errorbar size 1.000000                                                                   '
  write(*,*) '@    s0 errorbar linewidth 1.0                                                                   '
  write(*,*) '@    s0 errorbar linestyle 1                                                                     '
  write(*,*) '@    s0 errorbar riser linewidth 1.0                                                             '
  write(*,*) '@    s0 errorbar riser linestyle 1                                                               '
  write(*,*) '@    s0 errorbar riser clip off                                                                  '
  write(*,*) '@    s0 errorbar riser clip length 0.100000                                                      '
  write(*,*) '@    s0 comment "Cols 1:2:5"                                                                   '
  write(*,*) '@    s0 legend  ""                                                                             '
  write(*,*) '@    s1 hidden false                                                                             '
  write(*,*) '@    s1 type xycolor                                                                             '
  write(*,*) '@    s1 symbol 0                                                                                 '
  write(*,*) '@    s1 symbol size 1.000000                                                                     '
  write(*,*) '@    s1 symbol color 1                                                                           '
  write(*,*) '@    s1 symbol pattern 1                                                                         '
  write(*,*) '@    s1 symbol fill color 1                                                                      '
  write(*,*) '@    s1 symbol fill pattern 0                                                                    '
  write(*,*) '@    s1 symbol linewidth 1.0                                                                     '
  write(*,*) '@    s1 symbol linestyle 1                                                                       '
  write(*,*) '@    s1 symbol char 65                                                                           '
  write(*,*) '@    s1 symbol char font 0                                                                       '
  write(*,*) '@    s1 symbol skip 0                                                                            '
  write(*,*) '@    s1 line type 1                                                                              '
  write(*,*) '@    s1 line linestyle 1                                                                         '
  write(*,*) '@    s1 line linewidth 1.0                                                                       '
  write(*,*) '@    s1 line color 1                                                                             '
  write(*,*) '@    s1 line pattern 1                                                                           '
  write(*,*) '@    s1 baseline type 0                                                                          '
  write(*,*) '@    s1 baseline off                                                                             '
  write(*,*) '@    s1 dropline off                                                                             '
  write(*,*) '@    s1 fill type 0                                                                              '
  write(*,*) '@    s1 fill rule 0                                                                              '
  write(*,*) '@    s1 fill color 1                                                                             '
  write(*,*) '@    s1 fill pattern 1                                                                           '
  write(*,*) '@    s1 avalue off                                                                               '
  write(*,*) '@    s1 avalue type 2                                                                            '
  write(*,*) '@    s1 avalue char size 1.000000                                                                '
  write(*,*) '@    s1 avalue font 0                                                                            '
  write(*,*) '@    s1 avalue color 1                                                                           '
  write(*,*) '@    s1 avalue rot 0                                                                             '
  write(*,*) '@    s1 avalue format general                                                                    '
  write(*,*) '@    s1 avalue prec 3                                                                            '
  write(*,*) '@    s1 avalue prepend ""                                                                      '
  write(*,*) '@    s1 avalue append ""                                                                       '
  write(*,*) '@    s1 avalue offset 0.000000 , 0.000000                                                        '
  write(*,*) '@    s1 errorbar on                                                                              '
  write(*,*) '@    s1 errorbar place both                                                                      '
  write(*,*) '@    s1 errorbar color 1                                                                         '
  write(*,*) '@    s1 errorbar pattern 1                                                                       '
  write(*,*) '@    s1 errorbar size 1.000000                                                                   '
  write(*,*) '@    s1 errorbar linewidth 1.0                                                                   '
  write(*,*) '@    s1 errorbar linestyle 1                                                                     '
  write(*,*) '@    s1 errorbar riser linewidth 1.0                                                             '
  write(*,*) '@    s1 errorbar riser linestyle 1                                                               '
  write(*,*) '@    s1 errorbar riser clip off                                                                  '
  write(*,*) '@    s1 errorbar riser clip length 0.100000                                                      '
  write(*,*) '@    s1 comment "copy of setdata G0.S0"                                                        '
  write(*,*) '@    s1 legend  ""                                                                             '
  write(*,*) '@target G0.S0                                                                                    '
  write(*,*) '@type xycolor                                                                                    '
end if
if ( i == 2 ) then
  write(*,*) '&'
  write(*,*) '@target G0.S1'
  write(*,*) '@type xy'
end if
if ( i == 3 ) then
  write(*,*) '&'
end if

return
end





