c
c program lovoresult: Utility for extracting results
c                     from the lovoalign output produced
c                     in database comparisons
c
c Usage: Compile: f77 lovoresult.f -o lovoresult
c Run:   lovoresult 1npdb.pdb lovoalign.log > results.dat
c
c Where: 1npdb.pdb is the name of the file that contains
c        the protein to be analised 
c        lovoalign.log is the log file of the lovoalign
c        comparison
c        results.dat is the ordered output file
c
c Output: A list of the comparisons of protein 1npbd.pdb
c         ordered from best match (greatest score) to
c         worst match.
c
c  For being used with lovoalign:
c  http://www.ime.unicamp.br/~martinez/lovoalign
c
c  L. Martinez, 2007
c

      integer maxcomp
      parameter(maxcomp=60000)

      integer i, bije(maxcomp), gaps(maxcomp),
     +        lflash(maxcomp), mflash, index(maxcomp),
     +        ncomp, n1, n2
      real score(maxcomp), rmsd(maxcomp)

      character*30 protein, p1, p2, p2list(maxcomp)
      character*60 file
      character*200 record

c Getting command line arguments

      call getarg(1,protein)
      call getarg(2,file)

c Read file

      open(10,file=file,status='old')
      ncomp = 0
      do while(.true.)
        read(10,100,end=30,err=30) record
        if(record(1:6).eq.'PROTS:') then
          read(record(7:200),*) p1, p2, n1, n2
          if(p1.eq.protein.or.p2.eq.protein) then
            read(10,100,end=30,err=30) record
            read(10,100,end=30,err=30) record
            ncomp = ncomp + 1
            if(p1.eq.protein) p2list(ncomp) = p2 
            if(p2.eq.protein) p2list(ncomp) = p1 
            read(record(8:200),*) score(ncomp), bije(ncomp),
     +                            rmsd(ncomp), gaps(ncomp)  
          end if
        end if
      end do
      close(10)
30    continue
100   format(a200)

c Output title
      
      write(*,*)
      write(*,140)
140   format(tr2,71('-'))
      write(*,*) ' LOVOALIGN ordered results for ', protein
      write(*,*) ' Log file: ', file
      write(*,140)
      write(*,*) ' SCORE: Structal score of the alignment. '
      write(*,*) ' NBIJE: Number of atoms of the bijection. '
      write(*,*) ' GAPS: Number of gaps of the bijection. '
      write(*,140)
      write(*,150) 
150   format('  STRUCTURE ',t38,'SCORE',t48,'NBIJE',t56,'GAPS',
     +       t66,'RMSD')
      write(*,140)

c Order the results from best to worst score 

      mflash = ncomp
      call flash1(score,ncomp,lflash,mflash,index)

c Output ordered results

      do i = ncomp, 1, -1
        write(*,200) p2list(index(i)), score(i), bije(index(i)), 
     +               gaps(index(i)), rmsd(index(i))
      end do
200   format(tr2,a30,tr1,f12.6,tr1,i6,tr1,i6,tr1,f12.6)
      write(*,140)

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                             c
c     Subroutine Flash1                                       c
c     SORTS ARRAY A WITH N ELEMENTS BY USE OF INDEX VECTOR L  c
c     OF DIMENSION M WITH M ABOUT 0.1 N.                      c
c     Karl-Dietrich Neubert, FlashSort1 Algorithm             c
c     in  Dr. Dobb's Journal Feb.1998,p.123                   c
c                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine flash1 (A, N, L, M, ind)

      real a(1), anmin, c1, hold, flash
      integer L(1), ind(1), i, n, nmax, m, k, ihold, nmove, j, iflash
C     ============================ CLASS FORMATION =====


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

C     =============================== PERMUTATION =====
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

C     ========================= STRAIGHT INSERTION =====
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

C     =========================== RETURN,END FLASH1 =====
      RETURN
      END




