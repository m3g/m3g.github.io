c
c Program cylinder: Fits cylinders inside proteins
c
c How to use it: 
c
c             Compile it with:
c
c             f77 -O3 cylinder.f -o cylinder
c
c First part: Obtain one long cylinder that is fitted inside
c             all atoms. This is used to obtain the main
c             axis of the cavity, which will be the reference
c             axis for the slabing.
c
c             The input file in this first run should contain only
c             the name of the pdb file, for example (# are commentaries):
c    
c             input.inp:
c
c             #-------------------
c             structure.pdb
c             #-------------------
c
c             Run it with: 
c
c             cylinder input.inp
c
c             This run will return the file to be opened in VMD:
c
c             structure.vmd
c
c             open it with: vmd -e structure.vmd
c
c             The final lines of the output to the sceen will 
c             contain the information for running the slabing
c             part.
c 
c             If the cylinder is what you expect (if the main axis
c             of the cylinder is the main axis of you cavity, for
c             example), than you can go to the slabing part. If it
c             is not, for now you write to leandromartinez98@gmail.com
c             for help.
c
c  Second part: Obtain cylinders that fit slabs of the protein, made
c             orthogonally to a user-defined axis, usually the axis
c             obtained in the previous run. The input file must be:
c
c             input.inp:
c            
c             #----------------
c             slab 2.0
c             step 0.5
c             axis -0.289  -0.324   1.313  -0.549   0.738  -0.391
c             structure.pdb
c             #----------------
c
c             That means: 
c             The slabs of the protein will 2.0 units long. Slabs that
c             are too short contain no atoms, so are meaningless. Slabs
c             as long as the protein will return the same result as the
c             first run.
c
c             The slab will move in steps of 0.5 units. That means that
c             some atoms will be included and other atoms will be
c             excluded, but the slabs contain atoms in common. If you
c             want all slabs to be mutually exclusive, set step=slab.
c
c             The axis defines the main axis orthogonal to which the 
c             slabs will be made. The first three numbers are a point
c             belonging to the axis and the other three are a direction,
c             usually (but not necessarily) a unitary vector.
c
c             structure.pdb is the pdb file containing the structure.
c
c             The slabing includes all main-chain atoms, that means, 
c             atoms named C, N and CA of the pdb file.
c
c             Run with:
c
c             cylinder input.inp
c
c             Three files will be produced if everything goes fine:
c
c             structure_slabs.dat
c             Contains a table in which the first column is the center
c             of the slab and the second column is the largest radius
c             obtained.
c
c             structure_slabs.vmd
c             VMD control file for visualizing the results. Open it
c             with vmd -e structure.vmd
c
c             structure_slabs.atoms
c             A file containing the information per slab: First column
c             contains the center of the slab. Then it contains the
c             actual interval of atoms considered, the radius obtained
c             and finally a vmd-type selection of the CA atoms that
c             belong to the slab, to help visualization. 
c
c This code uses ALGENCAN as the optimization subroutine. ALGENCAN
c is part of the TANGO project. 
c http://www.ime.usp.br/~egbirgin/tango
c
c L. Martinez, Institut Pasteur, Nov 2007.
c leandromartinez98@gmail.com
c
c

C     =================================================================
C     File: algencanma.f
C     =================================================================

C     =================================================================
C     Module: Main program
C     =================================================================

C     Last update of any of the component of this module: 
C
C     January 30, 2007.

C     Users are encouraged to download periodically updated versions of 
C     this code at the TANGO home page:
C
C     www.ime.usp.br/~egbirgin/tango/

C     *****************************************************************
C     *****************************************************************

C     ALGENCAN solves problems of the form
C     ------------------------------------
C
C     min f(x)
C
C     subject to
C
C             c_j(x)  = 0, j in E,
C             c_j(x) <= 0, j in I,
C             l <= x <= u,
C
C     where E is the set of indices of the equality constraints, I is
C     the set of indices of the inequality constraints, and there are
C     n variables and m constraints.
C
C     ALGENCAN is an Augmented Lagrangian method that uses GENCAN to
C     solve the bound-constrained problems.
C
C     ALGENCAN is part of the TANGO Project.
C
C     Visit the TANGO home page in the web:
C
C     www.ime.usp.br/~egbirgin/tango/

C     *****************************************************************

C     TANGO LICENSE:
C     --------------
C
C     TANGO is free software; you can redistribute it and/or modify it 
C     under the terms of the GNU General Public License as published by 
C     the Free Software Foundation. Non-free versions of TANGO are 
C     available under terms different from those of the General Public 
C     License. Professors J. M. Martínez (martinez@ime.unicamp.br, 
C     martinezimecc@gmail.com) or E. G. Birgin (egbirgin@ime.usp.br, 
C     egbirgin@gmail.com) should be contacted for more information 
C     related to such a license, future developments and/or technical 
C     support.
C
C     Every published work that uses ALGENCAN should cite:
C
C     R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt, 
C     "On Augmented Lagrangian methods with general lower-level 
C     constraints", to appear in SIAM Journal on Optimization.
C
C     and
C
C     R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt, 
C     "Augmented Lagrangian methods under the Constant Positive Linear 
C     Dependence constraint qualification", to appear in Mathematical
C     Programming.
C
C     Every published work that uses GENCAN should cite:
C
C     E. G. Birgin and J. M. Martínez, "Large-scale active-set 
C     box-constrained optimization method with spectral projected 
C     gradients", Computational Optimization and Applications 23, pp. 
C     101-125, 2002.
C
C     (See other related works at the TANGO home page.)

C     *****************************************************************

C     HOW TO USE ALGENCAN TO SOLVE YOUR OWN PROBLEM
C     --------------------------------------------
C  
C     You will find below 10 subroutines and that YOU SHOULD MODIFY to 
C     solve your own problem. The modifications that you must do are:
C  
C     1) In the SUBROUTINE INIP you must set the number of variables
C     of your problem (n), the initial point (x), the lower and upper
C     bounds on the variables (l and u), the number of constraints (m)
C     (set m equal to zero if your problems has no constraints), a
C     vector (equatn) that indicates, for each constraint, whether it 
C     is and equality constraint or an inequality (least-or-equal),
C     a vector (linear) that indicates, for each constraint, whether it 
C     is a linear constraint or not, and the initial estimation of the 
C     Lagrange multipliers (lambda) (you can set it equal to zero if 
C     you do not have a better estimative).
C  
C     2) SUBROUTINES EVALF and EVALC MUST also be MODIFIED to evaluate 
C     the objective function and the constraints, respectively.
C
C     3) SUBROUTINES EVALG and EVALJAC, to compute the gradient vector
C     of the objective function and the gradients of the constraints,
C     respectively, are NOT MANDATORY but HIGHLY RECOMMENDED. They
C     must be provided by the user if he/she set gtype = 0. If the 
C     user set gtype = 1 then finite differences approximations will
C     be used. (See the description of gtype parameter in subroutine 
C     param.) 
C
C     4) OPTIONALY, you can modify SUBROUTINES EVALH and EVALHC to
C     evaluate the Hessian matrices of the objective function and the
C     constraints. Both subroutines are just optional subroutines and
C     are required when hptype = 0. (See the description of hptype
C     parameter in subroutine param.)
C
C     5) OPTIONALY, another way of using second-order information is
C     to provide, not individual subroutines to compute the Hessians
C     of the objective function and the constraints, but to provide
C     a unique SUBROUTINE EVALHLP to compute the product of an
C     arbitrary vector times the Lagrangian Hessian matrix. There is
C     not a clear advantage in coding EVALHLP instead of coding
C     EVALH and EVALHC. This alternative was incorporated to be 
C     primarily used in the AMPL and CUTEr interfaces. This 
C     subroutine is optional and is required when hptype = 1. (See 
C     the description of hptype parameter in subroutine param.)
C  
C     6) OPTIONALLY, in the SUBROUTINE PARAM you may modify some 
C     arguments like feasibility and optimality tolerances and maximum 
C     number of iterations and functional evaluations. See the detailed
C     description of each argument in subroutine PARAM. If you prefer, 
C     leave all these parameters with the suggested values. 
C  
C     7) Alternatively, if you are interested in modifying any of the
C     parameters described in (5) but you prefer not to modify subroutine
C     PARAM, you can use the algencan.dat SPECIFICATION FILE. The
C     advantage of using the specification file is that you can test
C     the method with different parameters without compiling it every
C     time. This file IS NOT MANDATORY. So, if you will not modify the
C     default parameters or if you fill more comfortable modifying
C     subroutine PARAM, you do not need to worry about algencan.dat
C     specification file.
C  
C     The specification file can have any number of lines, with at most
C     80 characters per line. Blank lines are admitted and lines starting
C     with '*' or '#' are considered comment lines and will be ignored. 
C     Each line must have just a "keyword" or a keyword followed by an 
C     integer or a real value when required. As the interpreter is not 
C     case-sensitive, you can write the keywords in lower case, upper 
C     case or any combination of them. You can also insert blanks in any 
C     place before or after the keywords. Moreover, just the first 10 
C     characters of the keyword are relevant and the rest will be ignored. 
C     Available keywords are:
C  
C     ANALYTIC-GRADIENT                     
C     FINITE-DIFFERENCES-GRADIENT           
C     HESSIANS-PROVIDED                     
C     LAGRHESS-PRODUCT-PROVIDED             
C     INCREMENTAL-QUOTIENTS                 
C     BFGS-QN-APPROXIMATION                 
C     ADAPTIVE-HESSIAN                      
C     AUTO-BDSOLVER                         
C     GENCAN-BDSOLVER                       
C     BETRA-BDSOLVER                        
C     SPARSE-BETRA-BDSOLVER                 
C     UNPRECONDITIONED-CG                   
C     BFGS-QN-PRECONDITIONER                
C     AUTO-INITIAL-PENALTY-PARAMETERS       
C     MANUAL-INITIAL-PENALTY-PARAMETERS      <real-value> 
C     COLLECTIVE-PENALTY-PARAMETERS-UPDATE  
C     INDEPENDENT-PENALTY-PARAMETERS-UPDATE 
C     DESIRED-INFEASIBILITY-FRACTION         <real-value> 
C     PENALTY-PARAMETER-MULTIPLIER-INCREMENT <real-value> 
C     CHECK-DERIVATIVES                     
C     FEASIBILITY-TOLERANCE                  <real-value> 
C     OPTIMALITY-TOLERANCE                   <real-value> 
C     MAX-OUTER-ITERATIONS                   <integer-value>
C     MAX-INNER-ITERATIONS                   <integer-value>
C     MAX-FUNCTION-EVALUATIONS               <integer-value>
C     OUTPUT-DETAIL                          <integer-value>
C     NCOMP-ARRAY                            <integer-value>
C               
C     By default, ALGENCAN uses:
C  
C     ANALYTIC-GRADIENT               
C     ADAPTIVE-HESSIAN                
C     BFGS-QN-PRECONDITIONER          
C     AUTO-BDSOLVER                    
C     BFGS-QN-PRECONDITIONER           
C     AUTO-INITIAL-PENALTY-PARAMETERS  
C     COLLECTIVE-PENALTY-PARAMETERS-UPDATE  
C     DESIRED-INFEASIBILITY-FRACTION           0.5d0 
C     PENALTY-PARAMETER-MULTIPLIER-INCREMENT 1.0d+01 
C     FEASIBILITY-TOLERANCE                  1.0d-04  
C     OPTIMALITY-TOLERANCE                   1.0d-04
C     MAX-OUTER-ITERATIONS                        50
C     MAX-INNER-ITERATIONS                   1000000
C     MAX-FUNCTION-EVALUATIONS               5000000
C     OUTPUT-DETAIL                                1
C     NCOMP-ARRAY                                  5
C
C     and derivatives are not tested.
C  
C     8) Finally, SUBROUTINE ENDP can be modified by the user to write, 
C     save or compute any type of extra information related to the 
C     solution of the problem. Subroutine ENDP will be called just once 
C     after algencan found a solution of the problem. You can modify 
C     this subroutine for, for example, save the solution in a file 
C     chosed by you and in a special format to be processed later by 
C     another software.
C   
C     Besides the modifications that you must do to solve your own 
C     problem using ALGENCAN, ALGENCAN has many other arguments to
C     control a large variety of issues of the algorithmic behaviour. 
C     To see how to control these other arguments, see their detailed 
C     description at the begining of subroutine algencan (for arguments 
C     related to Augmented Lagrangians) and at the begining of 
C     subroutine gencan (for arguments related to the bound-constraints 
C     internal solver).

C     ******************************************************************
C     ******************************************************************

      program algencanma

      implicit none

C     PARAMETERS
      integer mmax,nmax
      parameter ( mmax      =  500000 )
      parameter ( nmax      =  500000 )

C     LOCAL SCALARS
      logical checkder,rhoauto
      character * 6 precond
      integer gtype,hptype,inform,intype,iprint,m,maxoutit,maxtotfc,
     +        maxtotit,n,ncomp,outiter,rhotype,totcgcnt,totfcnt,totgcnt,
     +        totiter
      double precision epsfeas,epsopt,f,nalpsupn,rhofrac,rhomult,snorm
      real time

C     LOCAL ARRAYS
      integer wi1(nmax)
      logical equatn(mmax),linear(mmax)
      double precision l(nmax),lambda(mmax),rho(mmax),u(nmax),wd1(mmax),
     +        wd2(mmax),wd3(nmax),wd4(mmax),wd5(mmax),wd6(nmax),
     +        wd7(nmax),wd8(nmax),wd9(nmax),wd10(nmax),wd11(nmax),
     +        wd12(nmax),wd13(nmax),wd14(nmax),wd15(nmax),wd16(nmax),
     +        wd17(nmax),wd18(nmax),wd19(nmax),x(nmax)

      integer maxatom
      parameter(maxatom=5000)
      integer i, j, natom, narg, length, nref, ncall, icall
      double precision coor(maxatom,3), cref(maxatom,3), step, slab,
     +                 xa(7), xbest(7), bnorm(maxatom),
     +                 bmax, bmin, bstep, xanorm
      character*4 name
      character*200 record, keyword, keyvalue, value, file,
     +              lines(maxatom)
      logical findaxis
      common/atoms/coor,natom

C     EXTERNAL SUBROUTINES
      external solver

c Reading parameters and coordinates of the proteins for the cylinder
c problem
 
c Reading input file from command line

      narg = iargc()
      if(narg.eq.0) then
        write(*,100)
100     format('Run with: ./cylinder list.dat',/,
     +         'where list.dat is a list of pdb files or',
     +         ' a single pdb file.')
        stop
      end if
      call getarg(1,record)

c Some default parameters

      step = 1.d0
      slab = 3.d0
      findaxis = .true.

c Number of different initial points of each problem

      ncall = 100

c Reading data from file

      open(10,file=record,status='old')
      do while(.true.)

        read(10,110,end=40,err=40) record
        if(keyword(record).eq.'step') then
          record = keyvalue(record,1)
          read(record,*) step
        else if(keyword(record).eq.'slab') then
          record = keyvalue(record,1)
          read(record,*) slab
        else if(keyword(record).eq.'axis') then
          findaxis = .false.
          value = keyvalue(record,1)      
          read(value,*) xa(1) 
          value = keyvalue(record,2)      
          read(value,*) xa(2) 
          value = keyvalue(record,3)      
          read(value,*) xa(3) 
          value = keyvalue(record,4)      
          read(value,*) xa(4) 
          value = keyvalue(record,5)      
          read(value,*) xa(5) 
          value = keyvalue(record,6)      
          read(value,*) xa(6) 
          write(*,*) ' Provided main axis: ', (xa(i),i=1,6)
          xanorm = dsqrt(xa(4)**2+xa(5)**2+xa(6)**2)
          xa(4) = xa(4) / xanorm
          xa(5) = xa(5) / xanorm
          xa(6) = xa(6) / xanorm

c Reading coordinates from pdb file

        else if(record(1:1).ne.'#') then
          file = keyword(record)
          open(20,file=file,status='old')
          write(*,*) ' Reading file: ', record(1:length(record))
          nref = 0
          do while(.true.)
            read(20,110,end=20,err=20) record
            if(record(1:4).eq.'ATOM') then
              read(record(13:16),*,err=15,end=15) name
              if(name.eq.'CA'.or.
     +           name.eq.'C'.or.
     +           name.eq.'N') then
                nref = nref + 1
                read(record(31:38),*,err=15,end=15) cref(nref,1)  
                read(record(39:46),*,err=15,end=15) cref(nref,2)  
                read(record(47:54),*,err=15,end=15) cref(nref,3) 
                lines(nref) = record 
              end if
            end if
          end do
15        continue
          write(*,*) ' ERROR READING COORDINATES FOR ATOM ', nref 
          stop
20        continue
110       format(a200)
          write(*,*) ' Number of atoms: ', nref
          close(20)
        end if
      end do
40    continue

      if(nref.lt.0) then
        write(*,*) ' ERROR: Number of atoms found < 3 ' 
        stop
      end if

c Printing run data:

      write(*,*) ' Step of coordinates: ', step
      write(*,*) ' Slab of coordinates: ', slab

c If the main axis was not provided will try to obtain it by fitting a
c cylinder to all atoms

      if(findaxis) then
        write(*,*) ' Fitting cylinder to all atoms, trying to find',
     +             ' main axis '

        natom = nref
        do i = 1, natom
          coor(i,1) = cref(i,1) 
          coor(i,2) = cref(i,2) 
          coor(i,3) = cref(i,3)
        end do 

C     SET UP PROBLEM DATA

        call inip(n,x,l,u,m,lambda,equatn,linear,1)

        xbest(7) = 0.d0
        write(*,*) ' Searching for cylinder: ', ncall,' trials.' 
        do icall = 1, ncall

c Initial point

          call xini(n,x,l,u,icall,1,xa)

          call param(rhoauto,rhotype,rhomult,rhofrac,m,rho,gtype,hptype,
     +intype,precond,checkder,epsfeas,epsopt,maxoutit,maxtotit,maxtotfc,
     +iprint,ncomp)
          call solver(n,x,l,u,m,lambda,equatn,linear,rhoauto,rhotype,
     +rhomult,rhofrac,rho,gtype,hptype,intype,precond,checkder,epsfeas,
     +epsopt,maxoutit,maxtotit,maxtotfc,iprint,ncomp,f,snorm,nalpsupn,
     +outiter,totiter,totfcnt,totgcnt,totcgcnt,time,inform,wi1,wd1,wd2,
     +wd3,wd4,wd5,wd6,wd7,wd8,wd9,wd10,wd11,wd12,wd13,wd14,wd15,wd16,
     +wd17,wd18,wd19)
 
c Save the best solution

          if(x(7).gt.xbest(7).and.inform.eq.0) then
            do j = 1, 7
              xbest(j) = x(j)
            end do 
            write(*,*) ' Found valid radius of: ', x(7),
     +                 ' on trial ', icall
          end if
        end do

c Write output file for viewing with vmd

        call tovmd(file,xbest,step)
        write(*,*) ' ---------------------------------------'
        write(*,*) ' RESULTS: ' 
        write(*,*) ' ---------------------------------------'
        write(*,*) ' Best radius obtained: ', xbest(7)
        write(*,*) ' Reference point: ', (xbest(i),i=1,3)
        write(*,*) ' Axis vector: ', (xbest(i),i=4,6)
        write(*,*) ' ---------------------------------------'
        write(*,*) ' Wrote VMD input file: ',
     +             file(1:length(file)-3)//'.vmd'
        write(*,*) ' Check if the cylinder is what you expect.'  
        write(*,*) ' ---------------------------------------'
        write(*,207) (xbest(i),i=1,6),file(1:length(file))
207     format('  The following lines should be added to',
     +         ' the input file for slabing: ',/,
     +         '    step 0.5 ',/, 
     +         '    slab 2.0 ',/,
     +         '    axis ',6(f8.3),/, 
     +         t5,a,/, 
     +         ' ---------------------------------------')
        stop               
      end if 

c
c Now will slab the cylinder in pieces and obtain the cylinder that fits
c best to each slab
c

c Obtaining the projection of each point into the axis of the main
c cylinder

      do i = 1, nref
        bnorm(i) = (cref(i,1)-xa(1))*xa(4) 
     +           + (cref(i,2)-xa(2))*xa(5) 
     +           + (cref(i,3)-xa(3))*xa(6)
      end do

c Find maximum and minimum bnorms and move them to the zero
      
      bmin = bnorm(1)
      bmax = bnorm(1)
      do i = 2, nref
        bmin = dmin1(bmin,bnorm(i))
        bmax = dmax1(bmax,bnorm(i))
      end do
 
c Reference position

      xa(1) = xa(1) + bmin*xa(4)
      xa(2) = xa(2) + bmin*xa(5)
      xa(3) = xa(3) + bmin*xa(6)

c Now distances are from zero (current xa position)

      do i = 1, nref
        bnorm(i) = bnorm(i) - bmin
      end do   
      bmax = bmax - bmin   
      write(*,*) ' Perpendicular length of the protein: ', bmax

c Fractionate protein and optimize cylinder if there at least 3 atoms in
c the current fraction

      open(97,file=file(1:length(file)-4)//'_slabs.atoms')
      open(98,file=file(1:length(file)-4)//'_slabs.dat')
      write(98,195)
195   format('# PR: Position relative to reference point',/,
     +       '# RADIUS: Radius of the corresponding cylinder',/,
     +       '# PZ: Position relative to minimum coordinate',/,
     +       '# PR    RADIUS    PZ ') 
      open(99,file=file(1:length(file)-4)//'_slabs.vmd')
      write(99,200) file(1:length(file)),
     +              file(1:length(file)),
     +              file(1:length(file)) 
200   format('#!/usr/local/bin/vmd',/,
     +       'proc vmd_draw_arrow {mol start end} {',/,
     +       'set middle [vecadd $start', 
     +       ' [vecscale 0.95 [vecsub $end $start]]]',/,
     +       'graphics $mol cylinder $start $middle', 
     +       ' radius 0.15 resolution 30',/,
     +       'graphics $mol cone $middle $end radius 0.5 ',
     +       ' resolution 30 }',/,
     +       'mol new ',a,' type pdb first 0 last -1 step 1',
     +       ' filebonds 1 autobonds 1 waitfor all',/,
     +       'mol new ',a,' type pdb first 0 last -1 step 1',
     +       ' filebonds 1 autobonds 1 waitfor all',/,
     +       'mol new ',a,' type pdb first 0 last -1 step 1',
     +       ' filebonds 1 autobonds 1 waitfor all',/,
     +       'mol top 1',/,
     +       'mol rename top axis',/,
     +       'mol delrep 0 top')
      write(99,202) xa(1)-xa(4)*bmax/5.d0,
     +              xa(2)-xa(5)*bmax/5.d0,
     +              xa(3)-xa(6)*bmax/5.d0,
     +              xa(1)+xa(4)*(bmax + bmax/5.d0),
     +              xa(2)+xa(5)*(bmax + bmax/5.d0),
     +              xa(3)+xa(6)*(bmax + bmax/5.d0),
     +              xa(1),xa(2),xa(3)
202   format('graphics top color orange',/,
     +       'draw arrow {',3(f8.3),'} {',3(f8.3),'}',/,
     +       'graphics top sphere {',3(f8.3),'} ',
     +       ' radius 0.4 resolution 20',/,
     +       'mol top 2',/,
     +       'mol rename top cylinders',/,
     +       'mol delrep 0 top',/,
     +       'graphics top color blue')
      write(97,204)
204   format('#   step  ',
     +       t16,'interval',
     +       t32,'radius',
     +       t39,'VMD-selection for CA atoms')

      ncall = 20
      bstep = -step
      do while(bstep.lt.(bmax-step))
        bstep = bstep + step
        natom = 0
        write(97,205) bstep+slab/2.d0, bstep, bstep+slab
205     format(f8.3,' [',f8.3,',',f8.3,'] ',$)
        value='name CA and resid'
        do i = 1, nref
          if(bnorm(i).ge.bstep.and.bnorm(i).le.(bstep+slab)) then
            natom = natom + 1
            coor(natom,1) = cref(i,1)
            coor(natom,2) = cref(i,2)
            coor(natom,3) = cref(i,3)
            record = lines(i)
            read(record(13:16),*) name
            if(name.eq.'CA') then         
              value=value(1:length(value))//' '//record(23:26)
            end if
          end if
        end do 
        if(natom.lt.3) write(97,*)
        write(*,*) ' ------------------------------------- ' 
        write(*,*) ' Searching within ', bstep, bstep+slab
        write(*,*) ' Number of atoms: ', natom
        write(*,*) ' ------------------------------------- '
        if(natom.ge.3) then

          x(1) = xa(1) + xa(4)*(bstep + slab/2.d0)
          x(2) = xa(2) + xa(5)*(bstep + slab/2.d0)
          x(3) = xa(3) + xa(6)*(bstep + slab/2.d0)
          x(4) = xa(4)
          x(5) = xa(5)
          x(6) = xa(6)
          call inip(n,x,l,u,m,lambda,equatn,linear,2)

          xbest(7) = 0.d0
          do icall = 1, ncall

c Initial point

            call xini(n,x,l,u,icall,2,xa)

            call param(rhoauto,rhotype,rhomult,rhofrac,m,rho,gtype,
     +hptype,intype,precond,checkder,epsfeas,epsopt,maxoutit,maxtotit,
     +maxtotfc,iprint,ncomp)
            call solver(n,x,l,u,m,lambda,equatn,linear,rhoauto,rhotype,
     +rhomult,rhofrac,rho,gtype,hptype,intype,precond,checkder,epsfeas,
     +epsopt,maxoutit,maxtotit,maxtotfc,iprint,ncomp,f,snorm,nalpsupn,
     +outiter,totiter,totfcnt,totgcnt,totcgcnt,time,inform,wi1,wd1,wd2,
     +wd3,wd4,wd5,wd6,wd7,wd8,wd9,wd10,wd11,wd12,wd13,wd14,wd15,wd16,
     +wd17,wd18,wd19)
 
c Save the best solution

            if(x(7).gt.xbest(7).and.inform.eq.0) then
              do j = 1, 7
                xbest(j) = x(j)
              end do
              write(*,*) ' Best radius up to now: ', xbest(7) 
            end if    
          end do

c Write vmd line for drawing cylinder

        write(*,*) ' Best radius found: ', xbest(7)
        write(98,*) bmin + bstep + slab/2.d0, xbest(7), 
     +              bstep + slab/2.d0
        write(99,30) xbest(1)-xbest(4)*step/2.d0,
     +               xbest(2)-xbest(5)*step/2.d0,
     +               xbest(3)-xbest(6)*step/2.d0,
     +               xbest(1)+xbest(4)*step/2.d0,
     +               xbest(2)+xbest(5)*step/2.d0,
     +               xbest(3)+xbest(6)*step/2.d0,
     +               xbest(7)
30      format('graphics top cylinder', 
     +         ' {',3(tr1,f8.3),' }', 
     +         ' {',3(tr1,f8.3),' }',         
     +         ' radius ',f8.3,
     +         ' resolution 30',
     +         ' filled 1')
        write(97,206) xbest(7),' '//value(1:length(value))
206     format(f8.3,a)

        end if 
      end do
      close(98)

      write(99,201)
201   format('mol top 0',/,
     +       'mol delrep 0 top',/,
     +       'mol representation Trace 0.300000 6.000000',/,
     +       'mol color Name',/,
     +       'mol selection {all}',/,
     +       'mol material Opaque',/,
     +       'mol addrep top',/,
     +       'mol selupdate 0 top 0',/,
     +       'mol colupdate 0 top 0',/,
     +       'mol scaleminmax top 0 0.000000 0.000000',/,
     +       'mol smoothrep top 0 0')

      close(99)

      write(*,*) '--------------------------------------'
      write(*,*) ' See results with: xmgrace ',
     +           file(1:length(file)-4)//'_slabs.dat'
      write(*,*) ' See results with: vmd -e ',
     +           file(1:length(file)-4)//'_slabs.vmd'
      write(*,*) ' CA atoms belonging to each coordinate are in file:'
      write(*,*) ' ',file(1:length(file)-4)//'_slabs.atoms'
      write(*,*) '--------------------------------------'
 
      stop
      end


C     ******************************************************************
C     ******************************************************************

      subroutine param(rhoauto,rhotype,rhomult,rhofrac,m,rho,gtype,
     +hptype,intype,precond,checkder,epsfeas,epsopt,maxoutit,maxtotit,
     +maxtotfc,iprint,ncomp)

C     SCALAR ARGUMENTS
      character * 6 precond
      logical checkder,rhoauto
      integer gtype,hptype,iprint,intype,m,maxoutit,maxtotfc,maxtotit,
     +        ncomp,rhotype
      double precision epsfeas,epsopt,rhofrac,rhomult

C     ARRAY ARGUMENTS
      double precision rho(m)

C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     M integer
C     ---------
C
C     Number of constraints. This parameter can be used by the user in 
C     case he/she would like to initialize the penalty parameters (one
C     per constraint).
C
C
C     On Return:
C     ==========
C
C     RHOAUTO logical
C     ---------------
C
C     indicates whether the initial penalty parameters will be the
C     ones given by the user (rhoauto = .false.) or will be computed 
C     automatically by the Augmented Lagrangian solver (rhoauto = 
C     .true.).
C
C     RHOTYPE integer
C     ---------------
C
C     indicates whether the penalty parameters will be updated all
C     toghether (rhotype=1) or individualy (rhotype=2)
C
C     RHOMULT double precision
C     ------------------------
C
C     when a penalty parameter is updated, it is multiplied by rhofrac.
C     So, rhofrac must be greater than 1. Suggested value rhofrac=10.0d0.
C
C     RHOFRAC double precision
C     ------------------------
C
C     this paramater between 0 and 1 is the improvement required in the 
C     infeasibility for not to update the penalty parameters. 
C
C     When rhotype=1, if the sup-norm of the constraints is smaller than 
C     or equal to rhofrac times the sup-norm of the constraints of the 
C     previous iteration then the penalty parameters are not updated. 
C     Otherwise, they are multiplied by rhomult.
C
C     When rhotype=2, the improvement in feasibility of each constraint 
C     is evaluated in separate. If the feasibility of a particular 
C     constraint is smaller than or equal to rhofrac times the
C     feasibility of the previous iteration then the penalty parameter
C     is not modified. Otherwise, it is multiplied by rhomult.
C
C     The suggested value for rhofrac is 0.5d0.
C
C     RHO double precision rho(m)
C     ---------------------------
C
C     This the vector of penalty parameters. The user can leave to the
C     method the task to find adequate initial values for the penalty 
C     parameters setting rhoauto = .true. Otherwise, he/she can set 
C     rhoauto = .false. and atribute a value for each one of the penalty 
C     parameters rho(i), i = 1, ..., m,
C
C     GTYPE integer
C     -------------
C
C     Type of first derivatives calculation according to the following 
C     convention:
C
C     0 means true first derivatives. In this case, subroutines evalg
C       and evaljac must be modified by the user to compute the 
C       gradient of the objective function and the gradients of the 
C       constraints, respectively.
C
C     1 means that a finite difference approximation will be used. In 
C       this case, subroutines evalg and evaljac may have an empty body 
C       but must be present. It is also recommended that those 
C       empty-body subroutines set flag = - 1. Last but not least, the 
C       option gtype = 1 is not cheap neither safe.
C
C     HPTYPE integer
C     --------------
C
C     Type of Hessian-vector product according to the following 
C     convention:
C
C     We will first describe the meaning if the choices in the Augmented 
C     Lagrangian framework. The way in which the product of the Hessian 
C     matrix of the Augmented Lagrangian by a vector will be done 
C     depends on the value of the parameter hptype in the following way:
C
C     9 means that an incremental quotients approximation without any 
C       extra consideration will be used. This option requires the 
C       evaluation of an extra gradient at each Conjugate Gradient 
C       iteration. If gtype = 0 then this gradient evaluation will be 
C       done using the user supplied subroutines evalg and evaljac 
C       (evalc will also be used). On the other hand, if gtype = 1, the 
C       gradient calculation will be done using just calls to the user 
C       provided subroutines evalf and evalc. nind calls will be done, 
C       where nind is the dimension of the current face of the 
C       active-set method. This option is not cheap neither safe.
C
C       If you did not code subroutines evalg and evaljac, to compute 
C       the gradient of the objective function and the Jacobian of the 
C       constraints then your options finished here.
C
C     0 means that subroutines to compute the Hessian of the objective 
C       function (evalh) and the Hessians of the constraints (evalhc) 
C       were provided by the user. So, the product of the Hessian of the 
C       Augmented Lagrangian times a vector will be computed using the
C       Hessians provided by these subroutines and then adding the first 
C       order term (for the first order term the user-supplied 
C       subroutine to compute the Jacobian of the constraints (evaljac) 
C       is also used).
C
C     1 means that, instead of providing individual subroutines to 
C       compute the Hessians of the objective function and the 
C       constraints, the user provided a subroutine to compute the 
C       product of an arbitrary vector times the Hessian of the 
C       Lagrangian.
C
C     2 means that incremental quotients will be used. The difference 
C       between hptype = 9 and hptype = 2 is that, in the latter case, 
C       the non-differentiability of the Hessian of the Augmented 
C       Lagrangian will be taken into account. In particular, when 
C       computing the gradient of the Augmented Lagrangian at 
C       (x + step p), the constraints that will be considered will be 
C       the same constraints that contributed to the computation of the
C       gradient of the Augmented Lagrangian at x.
C
C       If GENCAN is been used to solve a bound-constrained problem 
C       which is not the subproblem of an Augmented Lagrangian method, 
C       but an isolated bound-constrained problem, then there is no 
C       difference between this option and hptype = 9.
C
C       This option also requires the evaluation of an extra gradient at 
C       each Conjugate Gradient iteration. 
C
C     3 is similar to hptype = 2. The difference is that the 
C       contribution of the linear constraints in the Hessian matrix 
C       will be computed explicitly. If the problem has not linear 
C       constraints then this option is identical to hptype = 2. 
C       Moreover, if the problem has no constraints then this option is 
C       equal to hptype = 9.
C
C       This option also requires the evaluation of an extra gradient at 
C       each Conjugate Gradient iteration. 
C
C     4 means that the Hessian matrix will be approximated and then the 
C       product of the Hessian approximation by the vector will be 
C       computed exactly. In particular, the Hessian matrix will be 
C       approximated doing a BFGS correction to the Gauss-Newton
C       approximation of the Hessian. Before the BFGS correction, a 
C       structured spectral correction is done to force the Gauss-Newton 
C       approximation to be positive definite.
C
C       If the problem has not constraints then the approximation 
C       reduces to a BFGS approximation of the Hessian (without memory) 
C       and using the spectral approximation (instead of the identity) 
C       as initial approximation.
C
C       Numerical experiments suggested that this option is convenient 
C       just for constrained problems. This motivated the introduction 
C       of the next option.
C
C       This option does NOT require an extra gradient evaluation per 
C       iteration and, in this sense, each CG iteration is 
C       computationally cheaper than a CG iteration of the previous 
C       choices. However, the approximation of the Hessian matrix 
C       requires some information (mainly the Jacobian of the 
C       constraints) that must be saved during the gradient evaluation. 
C       To save this information requires an amount of memory 
C       proportional to the number of non-null elements of the Jacobian 
C       matrix.
C 
C       Quadratic subproblems are convex with this choice.
C
C     5 is an adaptive strategy that choose, at every iteration, between
C       2 and 4. When the gradient of the Augmented Lagrangian is 
C       computed, it is verified if at least a constraint contributes to 
C       the calculation. If this is the case, 4 is used. Otherwise, 2 is 
C       used.
C
C       For problems with equality constraints (that always contributes 
C       to the Augmented Lagrangian function) this option is identical 
C       to 4.
C
C       For problems without constraints this option is identical to 2.
C
C     6 is identical to 5 but the choice is made between 3 and 4 instead 
C       of between 2 and 4. 
C
C       For problems with equality constraints (that always contributes 
C       to the Augmented Lagrangian function) this option is identical 
C       to 4. 
C
C       For problems without constraints this option is identical to 3.
C
C     We will now describe the meaning if the choices for unconstrained
C     and bound-constrained problems. In this context the way in which 
C     the product of the Hessian matrix by a vector will be done depends 
C     on the value of the parameter hptype in the following way:
C
C     0 means that the subroutine to compute the Hessian of the 
C       objective function (evalh) was provided by the user. So, the 
C       product of the Hessian times a vector will be computed using the 
C       Hessian provided by this subroutine.
C
C     1 means that a subroutine (evalhlp) to compute the product of the 
C       Hessian of the objective function times an arbitrary vector is
C       being provided by the user.
C
C     9 means that an incremental quotients approximation will be used. 
C       This option requires the evaluation of an extra gradient at each
C       Conjugate Gradient iteration. If gtype = 0 then this gradient 
C       evaluation will be done using the user supplied subroutine 
C       evalg. On the other hand, if gtype = 1, the gradient calculation 
C       will be done using just calls to the user provided subroutine 
C       evalf. nind calls will be done, where nind is the dimension of 
C       the current face of the active-set method.
C
C       If you did not code subroutine evalg to compute the gradient of 
C       the objective function then your options finished here.
C
C     4 means that the Hessian matrix will be approximated and then the 
C       product of the Hessian approximation by the vector will be 
C       computed exactly. The approximation is a BFGS approximation of 
C       the Hessian (without memory) and using the spectral 
C       approximation (instead of the identity) as initial 
C       approximation.
C
C       Numerical experiments suggested that this option is not 
C       convenient for unconstrained or just bound-constrained problems. 
C       (Note that this option was developed to be used in the Augmented 
C       Lagrangian framework.)
C
C       This option does NOT require an extra gradient evaluation per 
C       iteration and, in this sense, each CG iteration is 
C       computationally cheaper than a CG iteration of the previous 
C       choices.
C
C       Quadratic subproblems are convex with this choice.
C
C     In the bound-constrained context, options hptype = 2,3,5 and 6
C     are all identical to hptype = 9.
C
C     INTYPE integer
C     --------------
C
C     This parameter can be used to select the algorithm that will be
C     used to solve the Augmented Lagrangian subproblems, or the problem
C     itself if it an unconstrained problem. In the present implementation
C     just one algorithm is available. So, the value of this variable is
C     ignored.
C
C     PRECOND character * 6
C     ---------------------
C
C     Indicates the type of preconditioning that will be used for 
C     Conjugates Gradients according to the following convention:
C
C     'NONE'   means no preconditioner at all.
C
C     'QNAGNC' means Quasi-Newton Correction of the Gauss-Newton 
C              approximation of the Hessian. The exact form is this 
C              preconditioner is described in:
C 
C              E. G. Birgin and J. M. Martínez, "Structured minimal-
C              memory inexact quasi-Newton method and secant 
C              preconditioners for Augmented Lagrangian Optimization", 
C              submitted, 2005.
C
C     CHECKDER logical
C     ----------------
C
C     If you are using finite differences aproximations for the 
C     derivatives (gtype = 1) then checkder must be set to FALSE. 
C     On the other hand, if you are using your own coded derivatives 
C     you may would like to test them against a finite differences 
C     approximation to verify if they are o.k. In this case, set 
C     checkder = TRUE.
C
C     EPSFEAS double precision
C     ------------------------
C
C     Feasibility tolerance for the sup-norm of the constraints.
C     (Ignored in the unconstrained and bound-constrained cases.)
C
C     EPSOPT double precision
C     -----------------------
C
C     Optimality tolerance for the sup-norm of the projected gradient of 
C     the Augmented Lagrangian in the constrained case and the sup-norm
C     of the projected gradient of the objective function in the 
C     unconstrained and the bound-constrained cases.
C
C     MAXOUTIT integer
C     ----------------
C
C     Maximum number of Augmented Lagrangian (outer) iterations.
C     (Ignored in the unconstrained and bound-constrained cases.)
C
C     MAXTOTIT integer
C     ----------------
C
C     Maximum total number of inner iterations in the Augmented 
C     Lagrangian context (total means summing up the inner iterations of 
C     each outer iteration). In the unconstrained and bound-constrained
C     cases it means just the maximum number of iterations.
C
C     MAXTOTFC integer
C     ----------------
C
C     Idem MAXTOTIT but for number of functional evaluations.
C
C     IPRINT integer
C     --------------
C                
C     Controls the ammount of information of the output according to the 
C     following convention:
C
C     0 means no output at all.
C
C     1 means information at each outer iteration but without any 
C       information of how the subproblems are being solved.
C
C     2 means the same as 1 plus information of each inner iteration.
C
C     3 means the same as 2 plus information of the line searches and 
C       the calculation of the truncated Newton direction (using CG) of 
C       each inner iteration.
C
C     In all cases, an output file named solution.txt with the final 
C     point, Lagrange mutipliers and penalty parameters will be 
C     generated. Moreover, the same output of the screen will be saved 
C     in a file named algencan.out.
C
C     NCOMP integer
C     -------------
C
C     Every time a vector is printed, just its first ncomp component 
C     will be displayed.

      gtype    = 0
      hptype   = 6

      intype   = 0

      precond  = 'QNCGNA'

      rhoauto  = .true.

      rhotype  = 1
      rhomult  = 1.0d+01
      rhofrac  = 0.5d0

c voltar
      checkder = .false.

      epsfeas  = 1.0d-04
      epsopt   = 1.0d-04

      maxoutit = 50
      maxtotit = 1000000
      maxtotfc = 5 * maxtotit

c voltar
      iprint   = 0
      ncomp    = 5

      end

       


C     =================================================================
C     File: toyprob.f
C     =================================================================

C     =================================================================
C     Module: Subroutines that define the problem
C     =================================================================

C     Last update of any of the component of this module: 
 
C     February 7, 2006.

C     Users are encouraged to download periodically updated versions of 
C     this code at the TANGO home page:
 
C     www.ime.usp.br/~egbirgin/tango/ 

C     ******************************************************************
C     ******************************************************************

      subroutine inip(n,x,l,u,m,lambda,equatn,linear,itype)

      implicit none

C     This subroutine must set some problem data. For achieving this 
C     objective YOU MUST MODIFY it according to your problem. See below 
C     where your modifications must be inserted.
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     This subroutine has no input parameters.
C
C     On Return
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              initial point,
C
C     l        double precision l(n),
C              lower bounds on x,
C
C     u        double precision u(n),
C              upper bounds on x,
C
C     m        integer,
C              number of constraints (excluding the bounds),
C
C     lambda   double precision lambda(m),
C              initial estimation of the Lagrange multipliers,
C
C     equatn   logical equatn(m)
C              for each constraint j, set equatn(j) = .true. if it is an 
C              equality constraint of the form c_j(x) = 0, and set 
C              equatn(j) = .false. if it is an inequality constraint of 
C              the form c_j(x) <= 0,
C
C     linear   logical linear(m)
C              for each constraint j, set linear(j) = .true. if it is a 
C              linear constraint, and set linear(j) = .false. if it is a
C              nonlinear constraint.

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(*),linear(*)
      double precision l(*),lambda(*),u(*),x(*)

C     LOCAL SCALARS
      integer i

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                              c
c MODIFICATIONS THAT ARE DEPENDENT ON THE PARTICULAR PROBLEM   c
c BEING SOLVED (I. E. THE CYLINDER PROBLEM) START HERE         c
c                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer maxatom
      parameter(maxatom=5000)

      integer natom, itype
      double precision coor(maxatom,3), cmx, cmy, cmz,
     +                 xmax, ymax, zmax, xmin, ymin, zmin

      common/atoms/coor,natom
      common/bounds/cmx,cmy,cmz,xmin,ymin,zmin,xmax,ymax,zmax

c Computing the center of mass of the CA atoms, to use as an
c initial point for the axis of the cylinder

      cmx = 0.d0
      cmy = 0.d0
      cmz = 0.d0
      do i = 1, natom
        cmx = cmx + coor(i,1)
        cmy = cmy + coor(i,2)
        cmz = cmz + coor(i,3)
      end do
      cmx = cmx / dfloat(natom)
      cmy = cmy / dfloat(natom)
      cmz = cmz / dfloat(natom)

c
c Specifications for the optimization problem
c
      
c Number of variables

      n = 7
c
c Lower and upper bounds for the variables
c

c Reference point of the cylinder

      xmax = coor(1,1)
      ymax = coor(1,2)
      zmax = coor(1,3)
      xmin = coor(1,1)
      ymin = coor(1,2)
      zmin = coor(1,3)
      do i = 2, natom
        xmax = dmax1(xmax,coor(i,1)) 
        ymax = dmax1(ymax,coor(i,2)) 
        zmax = dmax1(zmax,coor(i,3)) 
        xmin = dmin1(xmin,coor(i,1)) 
        ymin = dmin1(ymin,coor(i,2)) 
        zmin = dmin1(zmin,coor(i,3))
      end do
      if(itype.eq.1) then
        l(1) = xmin 
        l(2) = ymin 
        l(3) = zmin 
        u(1) = xmax 
        u(2) = ymax 
        u(3) = zmax
      else if(itype.eq.2) then
        l(1) = x(1) - 1.d-3 
        l(2) = x(2) - 1.d-3 
        l(3) = x(3) - 1.d-3 
        u(1) = x(1) + 1.d-3 
        u(2) = x(2) + 1.d-3 
        u(3) = x(3) + 1.d-3
      end if
 
      if(itype.eq.1) then
        write(*,*) ' Maximum and minimum coordinates: ' 
        write(*,*) xmin, ymin, zmin, xmax, ymax, zmax
      end if

      if(itype.eq.1) then
        xmin = cmx - 5.d0 
        ymin = cmy - 5.d0 
        zmin = cmz - 5.d0 
        xmax = cmx + 5.d0
        ymax = cmy + 5.d0
        zmax = cmz + 5.d0
      else if(itype.eq.2) then
        xmin = -1.d5
        ymin = -1.d5
        zmin = -1.d5
        xmax = 1.d5
        ymax = 1.d5
        zmax = 1.d5
      end if

      if(itype.eq.1) then
        write(*,*) ' Number of atoms in this group: ', natom

        write(*,*) ' Center of mass of this group of atoms: ' 
        write(*,*) cmx, cmy, cmz
      
        write(*,*) ' Bounds for the center of the cylinder: ' 
        write(*,*) l(1),l(2),l(3),u(1),u(2),u(3)

        write(*,*) ' Restraints for the center of the cylinder: ' 
        write(*,*) xmin, ymin, zmin, xmax, ymax, zmax
      end if

c Direction

      if(itype.eq.1) then
        do i = 4, 6
          l(i) = -1.d5
          u(i) = 1.d5
        end do
      else if(itype.eq.2) then
        l(4) = x(4) - 1.d-3
        l(5) = x(5) - 1.d-3
        l(6) = x(6) - 1.d-3
        u(4) = x(4) + 1.d-3
        u(5) = x(5) + 1.d-3
        u(6) = x(6) + 1.d-3
      end if

      if(itype.eq.1) then
        write(*,*) ' Bounds for the possible axis: '
        write(*,*) l(4),l(5),l(6),u(4),u(5),u(6) 
      end if

c Radius

      l(7) = 0.d0
      u(7) = (xmax - xmin)**2 +
     +       (ymax - ymin)**2 +
     +       (zmax - zmin)**2  

      if(itype.eq.1) then
        write(*,*) ' Bounds for the values of the radius: '
        write(*,*) l(7), u(7)
      end if

c Number of constraints (the number of CA atoms, which cannot be inside
c of the cylinder at the solution)

      m = natom + 7

C     Lagrange multipliers approximation. Most users prefer to use the 
C     null initial Lagrange multipliers estimates. However, if the 
C     problem that you are solving is "slightly different" from a 
C     previously solved problem of which you know the correct Lagrange 
C     multipliers, we encourage you to set these multipliers as initial 
C     estimates. Of course, in this case you are also encouraged to use 
C     the solution of the previous problem as initial estimate of the 
C     solution. Similarly, most users prefer to use rho = 10 as initial 
C     penalty parameters. But in the case mentioned above (good 
C     estimates of solution and Lagrange multipliers) larger values of 
C     the penalty parameters (say, rho = 1000) may be more useful. More 
C     warm-start procedures are being elaborated.

      do i = 1,m
        lambda(i) =  0.0d0
      end do

C     For each constraint i, set equatn(i) = .true. if it is an equality
C     constraint of the form c_i(x) = 0, and set equatn(i) = .false. if 
C     it is an inequality constraint of the form c_i(x) <= 0.

      equatn(1) = .true.      
      do i = 2, m
        equatn(i) = .false.
      end do

C     For each constraint i, set linear(i) = .true. if it is a linear
C     constraint, otherwise set linear(i) = .false.

      linear(1) = .false.
      do i = 2, 7
        linear(i) = .true.
      end do
      do i = 8, m
        linear(i) = .false.
      end do

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE INIP.
C     ******************************************************************

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalf(n,x,f,flag)

      implicit none

C     This subroutine must compute the objective function. For achieving 
C     this objective YOU MUST MODIFY it according to your problem. See 
C     below where your modifications must be inserted.
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     On Return
C
C     f        double precision,
C              objective function value at x,
C
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the objective 
C              function. (For example, trying to compute the square root 
C              of a negative number, dividing by zero or a very small 
C              number, etc.) If everything was o.k. you must set it 
C              equal to zero.

C     SCALAR ARGUMENTS
      integer flag,n
      double precision f

C     ARRAY ARGUMENTS
      double precision x(n)

C     Objective function

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR OBJECTIVE
C     FUNCTION:
C     ******************************************************************

      flag = 0

c For each atom, computes the distance of the atom to the axis of the
c cylinder. The objective function is the radius of the cylinder, which
c is also the variable number 7.      

      f = -x(7)

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALF.
C     ******************************************************************

      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine evalg(n,x,g,flag)

      implicit none

C     This subroutine must compute the gradient vector of the objective 
C     function. For achieving these objective YOU MUST MODIFY it in the 
C     way specified below. However, if you decide to use numerical 
C     derivatives (we dont encourage this option at all!) you dont need
C     to modify evalg.
C
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     On Return
C
C     g        double precision g(n),
C              gradient vector of the objective function evaluated at x,
C
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of any component 
C              of the gradient vector. (For example, trying to compute 
C              the square root of a negative number, dividing by zero or 
C              a very small number, etc.) If everything was o.k. you 
C              must set it equal to zero.

C     SCALAR ARGUMENTS
      integer flag,n

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     LOCAL SCALARS
      integer i

C     Gradient vector

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET THE GRADIENT
C     VECTOR OF YOUR OBJECTIVE FUNCTION: 
C     ******************************************************************

      flag = 0

      do i = 1, 6
        g(i) = 0.d0
      end do
      g(7) = -1.d0

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALG. 
C     ******************************************************************
 
      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine evalh(n,x,hlin,hcol,hval,nnzh,flag)

      implicit none

C     This subroutine might compute the Hessian matrix of the objective 
C     function. For achieving this objective YOU MAY MODIFY it according 
C     to your problem. To modify this subroutine IS NOT MANDATORY. See 
C     below where your modifications must be inserted.
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     On Return
C
C     nnzh     integer,
C              number of perhaps-non-null elements of the computed 
C              Hessian,
C
C     hlin     integer hlin(nnzh),
C              see below,
C
C     hcol     integer hcol(nnzh),
C              see below,
C
C     hval     double precision hval(nnzh),
C              the non-null value of the (hlin(k),hcol(k)) position 
C              of the Hessian matrix of the objective function must 
C              be saved at hval(k). Just the lower triangular part of
C              Hessian matrix must be computed,
C
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the Hessian
C              matrix of the objective funtion. (For example, trying 
C              to compute the square root of a negative number, 
C              dividing by zero or a very small number, etc.) If 
C              everything was o.k. you must set it equal to zero.

C     SCALAR ARGUMENTS
      integer flag,n,nnzh

C     ARRAY ARGUMENTS
      integer hcol(*),hlin(*)
      double precision hval(*),x(n)

C     ******************************************************************
C     FROM HERE ON YOU MAY (OPTIONALY) MODIFY THE SUBROUTINE TO SET THE 
C     HESSIAN MATRIX OF YOUR OBJECTIVE FUNCTION: 
C     ******************************************************************

      flag = 0

      nnzh = 0 

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALH. 
C     ******************************************************************
 
      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine evalc(n,x,ind,c,flag)

      implicit none

C     This subroutine must compute the ind-th constraint of your 
C     problem. For achieving this objective YOU MUST MOFIFY it 
C     according to your problem. See below the places where your 
C     modifications must be inserted.
C
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     ind      integer,
C              index of the constraint to be computed,
C
C     On Return
C
C     c        double precision,
C              ind-th constraint evaluated at x,
C
C     flag     integer
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the 
C              constraint. (For example, trying to compute the square 
C              root of a negative number, dividing by zero or a very 
C              small number, etc.) If everything was o.k. you must set 
C              it equal to zero.
 
C     SCALAR ARGUMENTS

      integer maxatom
      parameter(maxatom=5000)
      integer ind,flag,n,natom
      double precision c,coor(maxatom,3), aa1, aa2, aa3, bb1, 
     +                 a1, a2, a3, sx7, d2,
     +                 xmin,ymin,zmin,xmax,ymax,zmax,
     +                 cmx, cmy, cmz

C     ARRAY ARGUMENTS
      double precision x(n)
      common/atoms/coor,natom
      common/bounds/cmx,cmy,cmz,xmin,ymin,zmin,xmax,ymax,zmax

C     Constraints

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET YOUR
C     CONSTRAINTS: 
C     ******************************************************************

      flag = 0

c First restraint: the direction vector is unitary

      if(ind.eq.1) then
        c = x(4)**2 + x(5)**2 + x(6)**2 - 1.d0

c Restraints from 2 to 7: The cylinder must be within
c the atoms (approximatelly)

      else if(ind.eq.2) then
        c = xmin - ( x(1) - x(7) )

      else if(ind.eq.3) then
        c = ymin - ( x(2) - x(7) )

      else if(ind.eq.4) then
        c = zmin - ( x(3) - x(7) )

      else if(ind.eq.5) then
        c = x(1) + x(7) - xmax

      else if(ind.eq.6) then
        c = x(2) + x(7) - ymax

      else if(ind.eq.7) then
        c = x(3) + x(7) - zmax

      else

c Restraints from 7 to natom + 7: atoms must not be inside the cylinder
c This computes the distance of a point to the axis of the cylinder
c defined at the current point, d2 is the squared distance

        a1 = x(1) - coor(ind-7,1) 
        a2 = x(2) - coor(ind-7,2)
        a3 = x(3) - coor(ind-7,3)
        aa1 = x(5)*a3 - x(6)*a2
        aa2 = x(6)*a1 - x(4)*a3
        aa3 = x(4)*a2 - x(5)*a1
        bb1 = 1.d0 / ( x(4)**2 + x(5)**2 + x(6)**2 )
        d2 = ( aa1**2 + aa2**2 + aa3**2 ) * bb1

        sx7 = x(7)**2
        c = sx7 - d2

      end if

      return

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALC. 
C     ******************************************************************
 
      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine evaljac(n,x,ind,indjac,valjac,nnzjac,flag)

      implicit none

C     This subroutine must compute the gradient of the ind-th constraint.
C     For achieving these objective YOU MUST MODIFY it in the way 
C     specified below.
C
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     ind      integer,
C              index of the constraint whose gradient will be computed,
C
C     On Return
C
C     nnzjac   integer,
C              number of perhaps-non-null elements of the computed 
C              gradient,
C
C     indjac   integer indjac(nnzjac),
C              see below,
C
C     valjac   double precision valjac(nnzjac),
C              the non-null value of the partial derivative of the 
C              ind-th constraint with respect to the indjac(k)-th 
C              variable must be saved at valjac(k).
C
C     flag     integer
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the 
C              constraint. (For example, trying to compute the square 
C              root of a negative number, dividing by zero or a very 
C              small number, etc.) If everything was o.k. you must set 
C              it equal to zero.

C     SCALAR ARGUMENTS
      integer maxatom
      parameter(maxatom=5000)
      integer flag,ind,n,nnzjac

C     ARRAY ARGUMENTS
      integer indjac(n), natom, i
      double precision x(n),valjac(n),coor(maxatom,3),
     +                 a1, a2, a3, aa1, aa2, aa3, bb1, bb2, 
     +                 dfdx(7), d2, sx7, aa12, aa22, aa32,
     +                 xmin,ymin,zmin,xmax,ymax,zmax,
     +                 cmx,cmy,cmz
      common/atoms/coor,natom
      common/bounds/cmx,cmy,cmz,xmin,ymin,zmin,xmax,ymax,zmax

C     Sparse gradient vector of the ind-th constraint

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET THE GRADIENTS  
C     OF YOUR CONSTRAINTS: 
C     ******************************************************************

      flag = 0

      if(ind.eq.1) then
        nnzjac = 3
        indjac(1) = 4
        indjac(2) = 5
        indjac(3) = 6
        valjac(1) = 2.d0*x(4)
        valjac(2) = 2.d0*x(5)
        valjac(3) = 2.d0*x(6)

      else if(ind.eq.2) then
        nnzjac = 2
        indjac(1) = 1
        indjac(2) = 7 
        valjac(1) = -1.d0
        valjac(2) = 1.d0

      else if(ind.eq.3) then
        nnzjac = 2
        indjac(1) = 2
        indjac(2) = 7 
        valjac(1) = -1.d0
        valjac(2) = 1.d0

      else if(ind.eq.4) then
        nnzjac = 2
        indjac(1) = 3
        indjac(2) = 7 
        valjac(1) = -1.d0
        valjac(2) = 1.d0

      else if(ind.eq.5) then
        nnzjac = 2
        indjac(1) = 1
        indjac(2) = 7 
        valjac(1) = 1.d0
        valjac(2) = 1.d0

      else if(ind.eq.6) then
        nnzjac = 2
        indjac(1) = 2
        indjac(2) = 7 
        valjac(1) = 1.d0
        valjac(2) = 1.d0

      else if(ind.eq.7) then
        nnzjac = 2
        indjac(1) = 3
        indjac(2) = 7 
        valjac(1) = 1.d0
        valjac(2) = 1.d0

      else
               
        a1 = x(1) - coor(ind-7,1) 
        a2 = x(2) - coor(ind-7,2)
        a3 = x(3) - coor(ind-7,3)
        aa1 = x(5)*a3 - x(6)*a2
        aa2 = x(6)*a1 - x(4)*a3
        aa3 = x(4)*a2 - x(5)*a1
        bb1 = 1.d0 / ( x(4)**2 + x(5)**2 + x(6)**2 )
        d2 = ( aa1**2 + aa2**2 + aa3**2 ) * bb1

        sx7 = x(7)**2

        aa12 = aa1**2
        aa22 = aa2**2
        aa32 = aa3**2
        bb2 = bb1**2
                   
        dfdx(1) = - ( 2.d0*aa2*x(6) + 2.d0*aa3*(-x(5)) ) * bb1
        dfdx(2) = - ( 2.d0*aa3*x(4) + 2.d0*aa1*(-x(6)) ) * bb1
        dfdx(3) = - ( 2.d0*aa1*x(5) + 2.d0*aa2*(-x(4)) ) * bb1
        dfdx(4) = - ( 2.d0*aa2*(-a3)*bb1 - bb2*2.d0*x(4)*aa22 )   
     +            - ( 2.d0*aa3*a2*bb1 - bb2*2.d0*x(4)*aa32 )
     +            - ( -2.d0*bb2*x(4)*aa12 ) 
        dfdx(5) = - ( 2.d0*aa1*a3*bb1 - bb2*2.d0*x(5)*aa12 ) 
     +            - ( 2.d0*aa3*(-a1)*bb1 - bb2*2.d0*x(5)*aa32 )
     +            - ( -2.d0*bb2*x(5)*aa22 )
        dfdx(6) = - ( 2.d0*aa1*(-a2)*bb1 - bb2*2.d0*x(6)*aa12 )
     +            - ( 2.d0*aa2*a1*bb1 - bb2*2.d0*x(6)*aa22 )   
     +            - ( -2.d0*x(6)*bb2*aa32 )
        dfdx(7) = 2.d0*x(7)

        nnzjac = 0 
        do i = 1, 7
          if(abs(dfdx(i)).gt.1.d-10) then
            nnzjac = nnzjac + 1
            indjac(nnzjac) = i
            valjac(nnzjac) = dfdx(i)
          end if
        end do
      end if

      return

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALJAC. 
C     ******************************************************************
 
      end

C     ******************************************************************
C     ******************************************************************
 
      subroutine evalhc(n,x,ind,hclin,hccol,hcval,nnzhc,flag)

      implicit none

C     This subroutine might compute the Hessian matrix of the ind-th
C     constraint. For achieving this objective YOU MAY MODIFY it 
C     according to your problem. To modify this subroutine IS NOT 
C     MANDATORY. See below where your modifications must be inserted.
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     ind      integer,
C              index of the constraint whose Hessian will be computed,
C
C     On Return
C
C     nnzhc    integer,
C              number of perhaps-non-null elements of the computed 
C              Hessian,
C
C     hclin    integer hclin(nnzhc),
C              see below,
C
C     hccol    integer hccol(nnzhc),
C              see below,
C
C     hcval    double precision hcval(nnzhc),
C              the non-null value of the (hclin(k),hccol(k)) position 
C              of the Hessian matrix of the ind-th constraint must 
C              be saved at hcval(k). Just the lower triangular part of
C              Hessian matrix must be computed,
C
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the Hessian
C              matrix of the ind-th constraint. (For example, trying 
C              to compute the square root of a negative number, 
C              dividing by zero or a very small number, etc.) If 
C              everything was o.k. you must set it equal to zero.

C     SCALAR ARGUMENTS
      integer flag,ind,n,nnzhc

C     ARRAY ARGUMENTS
      integer hccol(*),hclin(*)
      double precision hcval(*),x(n)

C     ******************************************************************
C     FROM HERE ON YOU MAY (OPTIONALY) MODIFY THE SUBROUTINE TO SET THE 
C     HESSIANS OF YOUR CONSTRAINTS: 
C     ******************************************************************

      flag = -1

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE EVALHC. 
C     ******************************************************************
 
      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalhlp(n,x,m,lambda,p,hp,goth,flag)

      implicit none

C     This subroutine might compute the product of the Hessian of the
C     Lagrangian times vector p (just the Hessian of the objective 
C     function in the unconstrained or bound-constrained case). 
C     
C     Parameters of the subroutine:
C
C     On Entry:
C
C     n        integer,
C              number of variables,
C
C     x        double precision x(n),
C              current point,
C
C     m        integer,
C              number of constraints,
C
C     lambda   double precision lambda(m),
C              vector of Lagrange multipliers,
C
C     p        double precision p(n),
C              vector of the matrix-vector product,
C
C     goth     logical,
C              can be used to indicate if the Hessian matrices were
C              computed at the current point. It is set to .false.
C              by the optimization method every time the current
C              point is modified. Sugestion: if its value is .false. 
C              then compute the Hessians, save them in a common 
C              structure and set goth to .true.. Otherwise, just use 
C              the Hessians saved in the common block structure,
C
C     On Return
C
C     hp       double precision hp(n),
C              Hessian-vector product,
C
C     goth     logical,
C              see above,
C              
C     flag     integer,
C              You must set it to any number different of 0 (zero) if 
C              some error ocurred during the evaluation of the 
C              Hessian-vector product. (For example, trying to compute 
C              the square root of a negative number, dividing by zero 
C              or a very small number, etc.) If everything was o.k. you 
C              must set it equal to zero.

C     SCALAR ARGUMENTS
      logical goth
      integer flag,m,n

C     ARRAY ARGUMENTS
      double precision hp(n),lambda(m),p(n),x(n)

      flag = - 1

      end

C     ******************************************************************
C     ******************************************************************

      subroutine endp(n,x,l,u,m,lambda,equatn,linear)

      implicit none

C     This subroutine can be used to do some extra job after the solver
C     has found the solution,like some extra statistics, or to save the
C     solution in some special format or to draw some graphical
C     representation of the solution. If the information given by the
C     solver is enough for you then leave the body of this subroutine
C     empty.
C     
C     Parameters of the subroutine:
C
C     The paraemters of this subroutine are the same parameters of
C     subroutine inip. But in this subroutine there are not output
C     parameter. All the parameters are input parameters.

C     SCALAR ARGUMENTS
      integer m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),u(n),x(n)

C     ******************************************************************
C     FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO COMPLEMENT THE
C     INFORMATION RELATED TO THE SOLUTION GIVEN BY THE SOLVER
C     ******************************************************************

C     ******************************************************************
C     STOP HERE YOUR MODIFICATIONS OF SUBROUTINE ENDP
C     ******************************************************************

      end


C     =================================================================
C     File: algencan.f
C     =================================================================

C     =================================================================
C     Module: Auxiliary subroutines
C     =================================================================

C     Last update of any of the component of this module: 
C
C     January 30, 2007.

C     ******************************************************************
C     ******************************************************************

      subroutine solver(n,x,l,u,m,lambda,equatn,linear,rhoauto,rhotype,
     +rhomult,rhofrac,rho,gtype,hptype,intype,precond,checkder,epsfeas,
     +epsopt,maxoutit,maxtotit,maxtotfc,iprint,ncomp,f,snorm,nalpsupn,
     +outiter,totiter,totfcnt,totgcnt,totcgcnt,time,inform,wi1,wd1,wd2,
     +wd3,wd4,wd5,wd6,wd7,wd8,wd9,wd10,wd11,wd12,wd13,wd14,wd15,wd16,
     +wd17,wd18,wd19)

      implicit none

C     This subroutine: 
C
C     (1) Computes some machine-dependent constants. 
C
C     (2) Process the optionally input specification file algencan.dat. 
C
C     (3) Check analytic derivatives if desired by the user.
C
C     (4) Open and close the output file algencan.out. 
C
C     (5) Computes the CPU elapsed time used by the solver.
C
C     (6) Calls the solver (GENCAN or ALGENCAN). 
C
C     (7) Write an output file called solution.txt with the solution.

C     SCALAR ARGUMENTS
      logical checkder,rhoauto
      character * 6 precond
      integer gtype,hptype,inform,intype,iprint,m,maxoutit,maxtotfc,
     +        maxtotit,n,ncomp,outiter,rhotype,totcgcnt,totfcnt,totgcnt,
     +        totiter
      double precision epsfeas,epsopt,f,nalpsupn,rhomult,rhofrac,snorm
      real time

C     ARRAY ARGUMENTS
      integer wi1(n)
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),rho(m),u(n),wd1(m),wd2(m),wd3(n),
     +        wd4(m),wd5(m),wd6(n),wd7(n),wd8(n),wd9(n),wd10(n),wd11(n),
     +        wd12(n),wd13(n),wd14(n),wd15(n),wd16(n),wd17(n),wd18(n),
     +        wd19(n),x(n)

C     LOCAL SCALARS
      integer i
      double precision bignum,macheps
C     double precision d1mach

      integer maxatom
      parameter(maxatom=5000)
      integer natom
      double precision coor(maxatom,3)
      common/atoms/coor,natom

C     DECLARATIONS RELATED TO THE TIME MEASUREMENT USING DTIME
      real dtime
      external dtime

      real dum(2)
      data dum/0.0,0.0/

C     COMPUTE MACHINE-DEPENDENT CONSTANTS

      bignum = 1.0d+99

C     macheps = d1mach(4)
      macheps = 1.0d-16

C     OPEN THE OUTPUT FILE

      open(10,file='algencan.out')

C     SET SOME SOLVER ARGUMENTS USING THE SPECIFICATION FILE

      call fparam(rhoauto,rhotype,rhomult,rhofrac,m,rho,gtype,hptype,
     +intype,precond,checkder,epsfeas,epsopt,maxoutit,maxtotit,maxtotfc,
     +iprint,ncomp)

C     TEST DERIVATIVES

c#ifndef CUTEr

      if ( checkder ) then
          call checkd(m,n,l,u,wi1,wd3,wd6,wd7,wd8,wd9,gtype,hptype,
     +    macheps,inform)

          if ( inform .lt. 0 ) then
              if ( iprint .ge. 1 ) then
                  write(*, 2000) inform
                  write(10,2000) inform
              end if

              go to 500
          end if
      end if

c#endif

C     CALL ALGENCAN IF THE PROBLEM HAS CONSTRAINTS. OTHERWISE, CALL GENCAN

      time = dtime(dum)

      if ( m .gt. 0 ) then

C         CALL THE AUGMENTED LAGRANGIAN SOLVER ALGENCAN

          call easyalgencan(n,x,l,u,m,lambda,equatn,linear,rhoauto,
     +    rhotype,rhomult,rhofrac,rho,gtype,hptype,intype,precond,
     +    epsfeas,epsopt,maxoutit,maxtotit,maxtotfc,iprint,ncomp,
     +    macheps,bignum,f,snorm,nalpsupn,outiter,totiter,totfcnt,
     +    totgcnt,totcgcnt,inform,wi1,wd1,wd2,wd3,wd4,wd5,wd6,wd7,wd8,
     +    wd9,wd10,wd11,wd12,wd13,wd14,wd15,wd16,wd17,wd18,wd19)

      else

C         CALL THE BOX-CONSTRAINTS SOLVER GENCAN

          call easygencan(n,x,l,u,m,lambda,equatn,linear,rho,gtype,
     +    hptype,intype,precond,epsopt,maxtotit,maxtotfc,iprint,ncomp,
     +    macheps,bignum,f,wd3,nalpsupn,totiter,totfcnt,totgcnt,
     +    totcgcnt,inform,wi1,wd6,wd7,wd8,wd9,wd10,wd11,wd12,wd13,wd14,
     +    wd15,wd16,wd17,wd18,wd19)

          outiter =     0
          snorm   = 0.0d0

      end if

      time = dtime(dum)
      time = dum(1)

C     CLOSE THE OUTPUT FILE

      close(10)

c 
c Saving solution for cylinder problem
c

c Project each atom into the cylinder line to get
c first and last cylinder coordinates for VMD

      open(10,file='vmd.txt')
      write(10,*) 'graphics 0 cylinder ',
     +           ' "',x(1)-10*x(4),x(2)-10*x(5),x(3)-10*x(6),'"',
     +           ' "',x(1)+10*x(4),x(2)+10*x(5),x(3)+10*x(6),'"',
     +           ' radius ',x(7),
     +           ' resolution 30 '
      close(10)

C     SAVE THE SOLUTION

      open(10,file='solution.txt')

C     Solution point
      write(10,9000)
      do i = 1,n
          write(10,9010) i,x(i)
      end do

C     Lagrange multipliers and penalty parameters
      if ( m .gt. 0 ) then
          write(10,9020)
          do i = 1,m
              write(10,9030) i,lambda(i),rho(i)
          end do
      end if

      close(10)

C     WRITE STATISTICS

      open(10,file='algencan-tabline.out')

      write(10,9050) time,inform,n,m,outiter,totiter,totfcnt,totgcnt,
     +totcgcnt,f,snorm,nalpsupn

      close(10)

C     RETURN

 500  continue

C     NON-EXECUTABLE STATEMENTS

 2000 format(/,' Flag of ALGENCAN = ',I3,' Fatal Error',/,
     +       /,' The following codes means: ',/,
     +       /,' -90 : error in evalf   subroutine', 
     +       /,' -91 : error in evalc   subroutine', 
     +       /,' -92 : error in evalg   subroutine', 
     +       /,' -93 : error in evaljac subroutine',
     +       /,' -94 : error in evalh   subroutine', 
     +       /,' -95 : error in evalhc  subroutine', 
     +       /,' -96 : error in evalhlp subroutine',/) 
 9000 format(/,'FINAL POINT:',//,2X,'INDEX',16X,'X(INDEX)')
 9010 format(  I7,1P,D24.16)
 9020 format(/,'FINAL ESTIMATION OF THE LAGRANGE MULTIPLIERS AND ',
     +         'PENALTY PARAMETERS: ',//,2X,'INDEX',11X,
     +         'LAMBDA(INDEX)',15X,'RHO(INDEX)')
 9030 format(I7,1P,D24.16,1X,1P,D24.16)
 9050 format(F8.2,1X,I3,1X,I6,1X,I6,1X,I3,1X,I7,1X,I7,1X,I7,1X,I7,1X,1P,
     +       D24.16,1X,1P,D7.1,1X,1P,D7.1)

      end

C     ******************************************************************
C     ******************************************************************
      subroutine fparam(rhoauto,rhotype,rhomult,rhofrac,m,rho,gtype,
     +hptype,intype,precond,checkder,epsfeas,epsopt,maxoutit,maxtotit,
     +maxtotfc,iprint,ncomp)

      implicit none

C     This subroutine set some algencan arguments related to stopping 
C     criteria and output. The setting values are taken fron file 
C     algencan.dat. Nothing is done if the file does not exist. File 
C     algencan.dat is an alternative to run ALGENCAN several times 
C     with different arguments without having to compile it again.

C     SCALAR ARGUMENTS
      logical checkder,rhoauto
      character * 6 precond
      integer gtype,hptype,intype,iprint,m,maxoutit,maxtotfc,maxtotit,
     +        ncomp,rhotype
      double precision epsfeas,epsopt,rhofrac,rhomult
 
C     ARRAY ARGUMENTS
      double precision rho(m)

C     Parameters of the subroutine:
C     =============================
C
C     The parameters of this subroutine are the same parameters of
C     subroutine param. See the description of them in subroutine 
C     param.

C     PARAMETERS
      integer nwords
      parameter ( nwords = 27 )

C     DATA BLOCKS
      character * 1 lower(26),upper(26)
      data lower /'a','b','c','d','e','f','g','h','i','j','k','l','m',
     +            'n','o','p','q','r','s','t','u','v','w','x','y','z'/
      data upper /'A','B','C','D','E','F','G','H','I','J','K','L','M',
     +            'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/

      character * 10 dictionary(nwords)
      data dictionary( 1) /'ANALYTIC-G'/
      data dictionary( 2) /'FINITE-DIF'/
      data dictionary( 3) /'HESSIANS-P'/
      data dictionary( 4) /'LAGRHESS-P'/
      data dictionary( 5) /'INCREMENTA'/
      data dictionary( 6) /'BFGS-QN-AP'/
      data dictionary( 7) /'ADAPTIVE-H'/
      data dictionary( 8) /'AUTO-BDSOL'/
      data dictionary( 9) /'GENCAN-BDS'/
      data dictionary(10) /'BETRA-BDSO'/
      data dictionary(11) /'SPARSE-BET'/
      data dictionary(12) /'UNPRECONDI'/
      data dictionary(13) /'BFGS-QN-PR'/
      data dictionary(14) /'AUTO-INITI'/
      data dictionary(15) /'MANUAL-INI'/
      data dictionary(16) /'COLLECTIVE'/
      data dictionary(17) /'INDEPENDEN'/
      data dictionary(18) /'DESIRED-IN'/
      data dictionary(19) /'PENALTY-PA'/
      data dictionary(20) /'CHECK-DERI'/
      data dictionary(21) /'FEASIBILIT'/
      data dictionary(22) /'OPTIMALITY'/
      data dictionary(23) /'MAX-OUTER-'/
      data dictionary(24) /'MAX-INNER-'/
      data dictionary(25) /'MAX-FUNCTI'/
      data dictionary(26) /'OUTPUT-DET'/
      data dictionary(27) /'NCOMP-ARRA'/

      character * 38 description(nwords)
      data description( 1) /'ANALYTIC-GRADIENT                     '/
      data description( 2) /'FINITE-DIFFERENCES-GRADIENT           '/
      data description( 3) /'HESSIANS-PROVIDED                     '/
      data description( 4) /'LAGRHESS-PRODUCT-PROVIDED             '/
      data description( 5) /'INCREMENTAL-QUOTIENTS                 '/
      data description( 6) /'BFGS-QN-APPROXIMATION                 '/
      data description( 7) /'ADAPTIVE-HESSIAN                      '/
      data description( 8) /'AUTO-BDSOLVER                         '/
      data description( 9) /'GENCAN-BDSOLVER                       '/
      data description(10) /'BETRA-BDSOLVER                        '/
      data description(11) /'SPARSE-BETRA-BDSOLVER                 '/
      data description(12) /'UNPRECONDITIONED-CG                   '/
      data description(13) /'BFGS-QN-PRECONDITIONER                '/
      data description(14) /'AUTO-INITIAL-PENALTY-PARAMETERS       '/
      data description(15) /'MANUAL-INITIAL-PENALTY-PARAMETERS     '/
      data description(16) /'COLLECTIVE-PENALTY-PARAMETERS-UPDATE  '/
      data description(17) /'INDEPENDENT-PENALTY-PARAMETERS-UPDATE '/
      data description(18) /'DESIRED-INFEASIBILITY-FRACTION        '/
      data description(19) /'PENALTY-PARAMETER-MULTIPLIER-INCREMENT'/
      data description(20) /'CHECK-DERIVATIVES                     '/
      data description(21) /'FEASIBILITY-TOLERANCE                 '/
      data description(22) /'OPTIMALITY-TOLERANCE                  '/
      data description(23) /'MAX-OUTER-ITERATIONS                  '/
      data description(24) /'MAX-INNER-ITERATIONS                  '/
      data description(25) /'MAX-FUNCTION-EVALUATIONS              '/
      data description(26) /'OUTPUT-DETAIL                         '/
      data description(27) /'NCOMP-ARRAY                           '/

      character * 1 addinfo(nwords)
      data addinfo( 1) /' '/
      data addinfo( 2) /' '/
      data addinfo( 3) /' '/
      data addinfo( 4) /' '/
      data addinfo( 5) /' '/
      data addinfo( 6) /' '/
      data addinfo( 7) /' '/
      data addinfo( 8) /' '/
      data addinfo( 9) /' '/
      data addinfo(10) /' '/
      data addinfo(11) /' '/
      data addinfo(12) /' '/
      data addinfo(13) /' '/
      data addinfo(14) /' '/
      data addinfo(15) /'D'/
      data addinfo(16) /' '/
      data addinfo(17) /' '/
      data addinfo(18) /'D'/
      data addinfo(19) /'D'/
      data addinfo(20) /' '/
      data addinfo(21) /'D'/
      data addinfo(22) /'D'/
      data addinfo(23) /'I'/
      data addinfo(24) /'I'/
      data addinfo(25) /'I'/
      data addinfo(26) /'I'/
      data addinfo(27) /'I'/

C     LOCAL SCALARS
      integer i,ifirst,ikey,ilast,inum,j
      double precision dnum

C     LOCAL ARRAYS
      character * 80 line
      character * 10 keyword

C     OPENING THE SPECIFICATION FILE
      open(20,err=300,file='algencan.dat',status='old')

C     MAIN LOOP
      write(*, 9005)
      write(10,9005)

 100  continue

C     READING LINES
      read(20,fmt=1000,err=400,end=200) line

C     PROCESS LINES

C     Find first character
      i = 1
 110  if ( i .le. 80 .and. line(i:i) .eq. ' ' ) then
          i = i + 1
          go to 110
      end if

C     Skip blank lines
      if ( i .gt. 80 ) then
          go to 100
      end if

      ifirst = i

C     Skip comments
      if ( line(ifirst:ifirst) .eq. '*' .or. 
     +     line(ifirst:ifirst) .eq. '#' ) then
          go to 100
      end if

C     Find the end of the keyword
      i = ifirst + 1
 120  if ( i .le. 80 .and. line(i:i) .ne. ' ' ) then
          i = i + 1
          go to 120
      end if

      ilast = i - 1

C     Obtain the first 10 characters and convert to upper-case
      keyword = '          '
      do i = 1,min( 10, ilast - ifirst + 1 )
          keyword(i:i) = line(ifirst+i-1:ifirst+i-1)
          do j = 1,26
              if ( keyword(i:i) .eq. lower(j) ) then
                  keyword(i:i) = upper(j)
              end if
          end do
      end do

C     Look up the keyword in the dictionary
      i = 1
 130  if ( i .le. nwords .and. keyword .ne. dictionary(i) ) then
          i = i + 1
          go to 130
      end if

C     Ignore unknown keywords
      if ( i .gt. nwords ) then
          write(*, 9020) line(ifirst:ilast)
          write(10,9020) line(ifirst:ilast)
          go to 100
      end if

      ikey = i

C     Read additional information if needed
      if ( addinfo(ikey) .ne. ' ' ) then

C         Skip blanks
          i = ilast + 1
 140      if ( i .le. 80 .and. line(i:i) .eq. ' ' ) then
              i = i + 1
              go to 140
          end if

C         Ignore keywords without the required information
          if ( i .gt. 80 ) then
              write(*, 9030) description(ikey)
              write(10,9030) description(ikey)
              go to 100
          end if 

C         Read additional information
          if ( addinfo(ikey) .eq. 'I' ) then
              read(unit=line(i:80),fmt=2000) inum

          else if ( addinfo(ikey) .eq. 'D' ) then
              read(unit=line(i:80),fmt=3000) dnum
          end if

      end if

C     Process keyword
      if ( addinfo(ikey) .eq. ' ' ) then
          write(*, 9040) description(ikey)
          write(10,9040) description(ikey)
      else if ( addinfo(ikey) .eq. 'I' ) then
          write(*, 9041) description(ikey),inum
          write(10,9041) description(ikey),inum
      else if ( addinfo(ikey) .eq. 'D' ) then
          write(*, 9042) description(ikey),dnum
          write(10,9042) description(ikey),dnum
      end if

C     Set the corresponding algencan argument
      if ( ikey .eq. 1 ) then
          gtype = 0
      else if ( ikey .eq.  2 ) then
          gtype = 1
      else if ( ikey .eq.  3 ) then
          hptype = 0
      else if ( ikey .eq.  4 ) then
          hptype = 1
      else if ( ikey .eq.  5 ) then
          hptype = 3
      else if ( ikey .eq.  6 ) then
          hptype = 4
      else if ( ikey .eq.  7 ) then
          hptype = 6
      else if ( ikey .eq.  8 ) then
          intype  = 0
      else if ( ikey .eq.  9 ) then
          intype  = 1
      else if ( ikey .eq. 10 ) then
          intype  = 2
      else if ( ikey .eq. 11 ) then
          intype  = 3
      else if ( ikey .eq. 12 ) then
          precond = 'NONE'
      else if ( ikey .eq. 13 ) then
          precond = 'QNCGNA'
      else if ( ikey .eq. 14 ) then
          rhoauto = .true.
      else if ( ikey .eq. 15 ) then
          rhoauto = .false.
          do i = 1,m
              rho(i) = dnum
          end do
      else if ( ikey .eq. 16 ) then
          rhotype  = 1
      else if ( ikey .eq. 17 ) then
          rhotype  = 2
      else if ( ikey .eq. 18 ) then
          rhofrac  = dnum
      else if ( ikey .eq. 19 ) then
          rhomult  = dnum
      else if ( ikey .eq. 20 ) then
          checkder = .true.
      else if ( ikey .eq. 21 ) then
          epsfeas = dnum
      else if ( ikey .eq. 22 ) then
          epsopt = dnum
      else if ( ikey .eq. 23 ) then
          maxoutit = inum
      else if ( ikey .eq. 24 ) then
          maxtotit = inum
      else if ( ikey .eq. 25 ) then
          maxtotfc = inum
      else if ( ikey .eq. 26 ) then
          iprint = inum
      else if ( ikey .eq. 27 ) then
          ncomp = inum
      end if

C     IIERATE
      go to 100

C     END OF LOOP

C     TERMINATIONS

C     CLOSING SPECIFICATION FILE
 200  continue
      close(20)
      go to 500

C     NO SPECIFICATION FILE
 300  continue
c      write(*, 9000)
      write(10,9000)
      go to 500

C     ERROR READING THE SPECIFICATION FILE
 400  continue
      write(*, 9010)
      write(10,9010)
      go to 500

C     PRINTING PARAMETERS VALUES
 500  continue

c      write(*, 4000) rhoauto,rhotype,rhomult,rhofrac,gtype,hptype,
c     +intype,precond,checkder,epsfeas,epsopt,maxoutit,maxtotit,maxtotfc,
c     +iprint,ncomp

c      write(10,4000) rhoauto,rhotype,rhomult,rhofrac,gtype,hptype,
c     +intype,precond,checkder,epsfeas,epsopt,maxoutit,maxtotit,maxtotfc,
c     +iprint,ncomp

C     NON-EXECUTABLE STATEMENTS
 1000 format(A80)
 2000 format(BN,I20)
 3000 format(BN,F24.0)
 4000 format(/,' ALGENCAN PARAMETERS:',
     +       /,' rhoauto  = ',     19X,L1,
     +       /,' rhotype  = ',        I20,
     +       /,' rhomult  = ',8X,1P,D12.4,
     +       /,' rhofrac  = ',8X,1P,D12.4,
     +       /,' gtype    = ',        I20,
     +       /,' hptype   = ',        I20,
     +       /,' intype   = ',        I20,
     +       /,' precond  = ',     14X,A6,
     +       /,' checkder = ',     19X,L1,
     +       /,' epsfeas  = ',8X,1P,D12.4,
     +       /,' epsopt   = ',8X,1P,D12.4,
     +       /,' maxoutit = ',        I20,
     +       /,' maxtotit = ',        I20,
     +       /,' maxtotfc = ',        I20,
     +       /,' iprint   = ',        I20,
     +       /,' ncomp    = ',        I20)
 9000 format(/,' The optional specification file algencan.dat was not',
     +         ' found in the current',/,' directory (this is not a',
     +         ' problem nor an error). The default values for the',/,
     +         ' ALGENCAN parameters will be used.')
 9005 format(/,' Specification file algencan.dat is being used.')
 9010 format(/,' Error reading specification file algencan.dat.')
 9020 format(  ' Ignoring unknown keyword ',A32)
 9030 format(  ' Ignoring incomplete keyword ',A32)
 9040 format(1X,A32)
 9041 format(1X,A32,5X,I20)
 9042 format(1X,A32,1X,1P,D24.8)

      end

c#ifndef CUTEr

C     *****************************************************************
C     *****************************************************************

      subroutine checkd(m,n,l,u,wi,wd1,wd2,wd3,wd4,x,gtype,hptype,
     +macheps,inform)

      implicit none

C     This subrotutine checks the user supplied first and second
C     derivatives subroutines (evalg, evalh, evaljac and evalhc) for
C     computing the objective function gradient and Hessian and the
C     constraints gradients and Hessians, respectively.

C     SCALAR ARGUMENTS
      integer gtype,hptype,inform,m,n
      double precision macheps

C     ARRAY ARGUMENTS
      integer wi(n)
      double precision l(n),u(n),wd1(n),wd2(n),wd3(n),wd4(n),x(n)

C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     M integer: Number of constraints.
C     ---------
C
C     N integer: Number of variables.
C     ---------
C
C     L double precision l(n): Lower bounds.
C     -----------------------
C
C     U double precision u(n): Upper bounds.
C     -----------------------
C
C     WI double precision wi(n): N-dimensional integer working vector.
C     -----------------------
C
C     WD1 ... WD4, X double precision wd1(n) ... wd4(n), x(n) 
C     -------------------------------------------------------
C
C     N-dimensional double precision working vectors.
C
C     GTYPE integer: Type of first derivatives calculation.
C     -------------
C
C     HPTYPE integer: Type of Hessian-vector product.
C     --------------
C
C     MACHEPS double precision
C     ------------------------
C
C     Smallest positive number such that 1 + macheps is not equal to 1.
C
C     On Return:
C     ==========
C
C     INFORM integer: Output flag.
C     --------------

C     LOCAL SCALARS
      character answer
      integer i,j
      double precision drand,seed,smalll,smallu

      inform = 0

C     SET A RANDOM POINT 

      seed = 17325.0d0
      do i = 1,n
          smalll = max( l(i), - 1.0d4 )
          smallu = min( u(i),   1.0d4 )
          if ( .not. smalll .lt. smallu ) then
              smalll = l(i)
              smallu = u(i)
          end if
          x(i) = smalll + ( smallu - smalll ) * drand(seed)
      end do
 
      write(*,900)

      x(1) = 8.
      x(2) = 18.
      x(3) = 4.
      x(4) = 8.
      x(5) = 7.
      x(6) = 9.
      x(7) = 1.d5

      do i = 1,n
          write(*,910) i,x(i)
      end do

C     CHECK OBJECTIVE FUNCTION GRADIENT

      write(*,920)

      read(*,*) answer

      if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
          go to 500

      else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
          call checkg(n,x,wd1,macheps,inform)
          if ( inform .lt. 0 ) then
              return
          end if 
      end if

C     CHECK CONSTRAINTS GRADIENTS

      j = 1

 110  if ( j .le. m ) then

          write(*,930) j

          read(*,*) answer

          if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
              go to 500

          else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
              go to 200

          else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
              call checkjac(n,x,j,wi,wd1,wd2,macheps,inform)
              if ( inform .lt. 0 ) then
                  return
              end if 
          end if

          j = j + 1

          go to 110

      end if

C     CHECK HESSIAN OF THE OBJECTIVE FUNCTION

 200  continue

      write(*,940)

      read(*,*) answer

      if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
          go to 500

      else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
          call checkh(n,x,wd1,wd2,wd3,macheps,inform)
          if ( inform .lt. 0 ) then
              return
          end if 
      end if

C     CHECK HESSIANS OF THE CONSTRAINTS

      j = 1

 310  if ( j .le. m ) then

          write(*,950) j

          read(*,*) answer

          if ( answer .eq. 'A' .or. answer .eq. 'a' ) then
              go to 500

          else if ( answer .eq. 'S' .or. answer .eq. 's' ) then
              go to 500

          else if ( answer .eq. 'Y' .or. answer .eq. 'y' ) then
              call checkhc(n,x,j,wi,wd1,wd2,wd3,wd4,macheps,inform)
              if ( inform .lt. 0 ) then
                  return
              end if 
          end if

          j = j + 1

          go to 310

      end if

C     RETURN

 500  continue

      write(*,*) 'Hit any letter to continue.'
      read(*,*) answer

 900  format(/,1X,'Derivatives will be tested at the random point: ')
 910  format(  1X,'x(',I6,') = ',1P,D15.8)
 920  format(/,1X,'Check the gradient of the objective function?',
     +       /,1X,'Type Y(es), N(o) or A(bort checking): ')
 930  format(/,1X,'Check the gradient of constraint ',I5,'?',
     +       /,1X,'Type Y(es), N(o), A(bort checking) or ',
     +            'S(kip constraints gradients): ')
 940  format(/,1X,'Check the Hessian matrix of the objective function?',
     +       /,1X,'Type Y(es), N(o) or A(bort checking): ')
 950  format(/,1X,'Check the Hessian of constraint ',I5,'?',
     +       /,1X,'Type Y(es), N(o), A(bort checking) or ',
     +            'S(kip constraints gradients): ')

      end

C     *****************************************************************
C     *****************************************************************

      subroutine checkg(n,x,g,macheps,inform)

      implicit none

C     This subrotutine checks the user supplied subroutine evalg for 
C     computing the gradient of the objective function using central
C     finite differences with two different discretization steps.

C     SCALAR ARGUMENTS
      integer inform,n
      double precision macheps

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     N integer: Number of variables.
C     ---------
C
C     X double precision x(n): Point at which the gradient will be tested.
C     -----------------------
C
C     G double precision g(n) 
C     -----------------------
C
C     N-dimensional double precision working vector.
C
C     MACHEPS double precision
C     ------------------------
C
C     Smallest positive number such that 1 + macheps is not equal to 1.
C
C     On Return:
C     ==========
C
C     INFORM integer: Output flag.
C     --------------

C     LOCAL SCALARS
      integer flag,i
      double precision fminus,fplus,gdiff1,gdiff2,maxerr,step1,step2,tmp

      inform = 0

      call evalg(n,x,g,flag)
      if ( flag .ne. 0 ) then
          inform = - 92
          return
      end if

      write(*,100)

      maxerr = 0.0d0

      do i = 1,n
          tmp  = x(i)

          step1 = macheps ** (1.0d0/3.0d0) * max( abs( tmp ), 1.0d0 )

          x(i) = tmp + step1
          call evalf(n,x,fplus,flag)
          if ( flag .ne. 0 ) then
              inform = - 90
              return
          end if

          x(i) = tmp - step1
          call evalf(n,x,fminus,flag)
          if ( flag .ne. 0 ) then
              inform = - 90
              return
          end if

          gdiff1 = ( fplus - fminus ) / ( 2.0d0 * step1 )

          step2 = macheps ** (1.0d0/3.0d0) * max( abs( tmp ), 1.0d-03 )

          x(i) = tmp + step2
          call evalf(n,x,fplus,flag)
          if ( flag .ne. 0 ) then
              inform = - 90
              return
          end if

          x(i) = tmp - step2
          call evalf(n,x,fminus,flag)
          if ( flag .ne. 0 ) then
              inform = - 90
              return
          end if

          x(i) = tmp

          gdiff2 = ( fplus - fminus ) / ( 2.0d0 * step2 )

          tmp = min( abs( g(i) - gdiff1 ), abs( g(i) - gdiff2 ) )
          write(*,110) i,g(i),gdiff1,gdiff2,tmp
          maxerr = max( maxerr, tmp )

      end do

      write(*,120) maxerr

      return

 100  format(/,1X,'Gradient vector of the objective function.',
     +       /,1X,'Index',13X,'evalg',2X,'Central diff (two different ',
     +            'steps)',4X,'Absolute error')
 110  format(  1X,I5,4(3X,1P,D15.8))
 120  format(  1X,'Maximum absolute error = ',1P,D15.8)

      end

C     *****************************************************************
C     *****************************************************************

      subroutine checkh(n,x,g,gplus1,gplus2,macheps,inform)

      implicit none

C     This subrotutine checks the user supplied subroutine evalh for 
C     computing the Hessian of the objective function using central 
C     finite differences with two different discretization steps.

C     SCALAR ARGUMENTS
      integer inform,n
      double precision macheps

C     ARRAY ARGUMENTS
      double precision g(n),gplus1(n),gplus2(n),x(n)

C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     N integer: Number of variables.
C     ---------
C
C     X double precision x(n): Point at which the Hessian will be tested.
C     -----------------------
C
C     G, GPLUS1, GPLUS2 double precision g(n), gplus1(n), gplus2(n) 
C     -------------------------------------------------------------
C
C     N-dimensional double precision working vectors.
C
C     MACHEPS double precision
C     ------------------------
C
C     Smallest positive number such that 1 + macheps is not equal to 1.
C
C     On Return:
C     ==========
C
C     INFORM integer: Output flag.
C     --------------

C     PARAMETERS
      integer nsmax
      parameter ( nsmax = 1000 )

C     LOCAL SCALARS
      logical nullcol
      integer flag,i,j,nnzh
      double precision elem,hdiff1,hdiff2,step1,step2,tmp

C     LOCAL ARRAYS
      integer hlin(nsmax**2),hcol(nsmax**2)
      double precision H(nsmax,nsmax),hval(nsmax**2),maxerr(nsmax)

      inform = 0

C     Test viability

      if ( n .gt. nsmax ) then
          write(*,*) 'Subroutine CheckH uses dense matrices up to ',
     +               'dimension ',nsmax,' times ',nsmax,'. The ',
     +               'Hessian checking will be skipped.'
          go to 500
      end if

C     Compute the gradient of the objective function at x

      call evalg(n,x,g,flag)
      if ( flag .ne. 0 ) then
          inform = - 92
          return
      end if

C     Compute the Hessian of the objective function at x and save in a
C     dense matrix

      call evalh(n,x,hlin,hcol,hval,nnzh,flag)
      if ( flag .ne. 0 ) then
          inform = - 94
          return
      end if

      do j = 1,n
          do i = 1,n
              H(i,j) = 0.0d0
          end do
      end do

      do i = 1,nnzh
          H(hlin(i),hcol(i)) = H(hlin(i),hcol(i)) + hval(i)
      end do

C     Test column by column

      write(*,100)

      do j = 1,n

          tmp  = x(j)

          step1 = sqrt( macheps ) * max( abs( tmp ), 1.0d0 )

          x(j) = tmp + step1
          call evalg(n,x,gplus1,flag)
          if ( flag .ne. 0 ) then
              inform = - 92
              return
          end if

          step2 = sqrt( macheps ) * max( abs( tmp ), 1.0d-03 )

          x(j) = tmp + step2
          call evalg(n,x,gplus2,flag)
          if ( flag .ne. 0 ) then
              inform = - 92
              return
          end if

          x(j) = tmp

          write(*,105) j

          maxerr(j) = 0.0d0

          nullcol = .true.

          do i = 1,n
              if ( i .ge. j ) then
                  elem = H(i,j)
              else
                  elem = H(j,i)
              end if
              hdiff1 = ( gplus1(i) - g(i) ) / step1
              hdiff2 = ( gplus2(i) - g(i) ) / step2
              tmp = min( abs( elem - hdiff1 ), abs( elem - hdiff2 ) )
              if ( elem .ne. 0.0d0 .or. 
     +             hdiff1 .ne. 0.0d0 .or. hdiff2 .ne. 0.0d0 ) then
                  if ( nullcol ) then
                      nullcol = .false.
                      write(*,106)
                  end if
                  write(*,110) i,elem,hdiff1,hdiff2,tmp
              end if
              maxerr(j) = max( maxerr(j), tmp )
          end do

          if ( nullcol ) then
              write(*,115)
          else
              write(*,120) maxerr(j)
          end if

      end do

      write(*,*)
      do j = 1,n
          write(*,130) j,maxerr(j)
      end do

 500  continue

      return

 100  format(/,1X,'Hessian matrix of the objective function column by ',
     +            'column.')
 105  format(/,1X,'Column:  ',I6)
 106  format(/,1X,'Index',13X,'evalh',3X,'Incr. Quoc. (two different ',
     +            'steps)',4X,'Absolute error')
 110  format(  1X,I5,4(3X,1P,D15.8))
 115  format(  1X,'All the elements of this column are null.')
 120  format(  1X,'Maximum absolute error = ',1P,D15.8)
 130  format(  1X,'Column ',I6,' Maximum absolute error = ',1P,D15.8)

      end

C     *****************************************************************
C     *****************************************************************

      subroutine checkjac(n,x,ind,indjac,valjac,g,macheps,inform)

      implicit none

C     This subrotutine checks the user supplied subroutine evaljac for 
C     computing the gradients of the constraints using central finite 
C     differences with two different discretization steps.

C     SCALAR ARGUMENTS
      integer ind,inform,n
      double precision macheps

C     ARRAY ARGUMENTS
      integer indjac(n)
      double precision g(n),valjac(n),x(n)

C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     N integer: Number of variables.
C     ---------
C
C     X double precision x(n) 
C     -----------------------
C
C     Point at which the gradient of the ind-th constraint will be 
C     tested.
C
C     IND integer
C     -----------
C
C     Index of the constraint whose gradient will be tested.
C
C     INDJAC integer indjac(n): N-dimensional integer working vector.
C     ------------------------
C
C     VALJAC, G double precision valjac(n), g(n) 
C     ------------------------------------------
C
C     N-dimensional double precision working vectors.
C
C     MACHEPS double precision
C     ------------------------
C
C     Smallest positive number such that 1 + macheps is not equal to 1.
C
C     On Return:
C     ==========
C
C     INFORM integer: Output flag.
C     --------------

C     LOCAL SCALARS
      logical nullcol
      integer flag,i,nnzjac
      double precision cminus,cplus,jacdiff1,jacdiff2,maxerr,step1,
     +        step2,tmp

      inform = 0

C     COMPUTE THE GRADIENT OF THE CONSTRAINT AND SAVE IT INTO A DENSE 
C     VECTOR

      call evaljac(n,x,ind,indjac,valjac,nnzjac,flag)
      if ( flag .ne. 0 ) then
          inform = - 93
          return
      end if

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,nnzjac
          g(indjac(i)) = g(indjac(i)) + valjac(i)
      end do

C     COMPARE WITH CENTRAL FINITE DIFFERENCES

      write(*,100) ind

      maxerr = 0.0d0

      nullcol = .true.

      do i = 1,n
          tmp  = x(i)

          step1 = macheps ** (1.0d0/3.0d0) * max( abs( tmp ), 1.0d0 )
          step1 = 1.d-8

          x(i) = tmp + step1
          call evalc(n,x,ind,cplus,flag)
          if ( flag .ne. 0 ) then
              inform = - 91
              return
          end if

          x(i) = tmp - step1
          call evalc(n,x,ind,cminus,flag)
          if ( flag .ne. 0 ) then
              inform = - 91
              return
          end if

          jacdiff1 = ( cplus - cminus ) / ( 2.0d0 * step1 )

          step2 = macheps * (1.0d0/3.0d0) * max( abs( tmp ), 1.0d-03 )
          step2 = 1.d-5

          x(i) = tmp + step2
          call evalc(n,x,ind,cplus,flag)
          if ( flag .ne. 0 ) then
              inform = - 91
              return
          end if

          x(i) = tmp - step2
          call evalc(n,x,ind,cminus,flag)
          if ( flag .ne. 0 ) then
              inform = - 91
              return
          end if

          x(i) = tmp

          jacdiff2 = ( cplus - cminus ) / ( 2.0d0 * step2 )

          tmp = min( abs( g(i) - jacdiff1 ), abs( g(i) - jacdiff2 ) )
              if ( nullcol ) then
                  nullcol = .false.
                  write(*,105)
              end if
              write(*,110) i,g(i),jacdiff1,jacdiff2,tmp
          maxerr = max( maxerr, tmp )
      end do

      if ( nullcol ) then
          write(*,115)
      else
          write(*,120) maxerr
      end if

      return

 100  format(/,1X,'Gradient vector of constraints ',I5,'.')
 105  format(/,1X,'Index',11X,'evaljac',2X,'Central diff (two ',
     +            'different steps)',4X,'Absolute error')
 110  format(  1X,I5,4(3X,1P,D15.8))
 115  format(  1X,'All the elements of this gradient are null.')
 120  format(  1X,'Maximum absolute error = ',1P,D15.8)

      end

C     *****************************************************************
C     *****************************************************************
      subroutine checkhc(n,x,ind,indjac,g,gplus1,gplus2,valjac,macheps,
     +inform)

      implicit none

C     This subrotutine checks the user supplied subroutine evalhc for 
C     computing the Hessians of the constraints using finite
C     differences.

C     SCALAR ARGUMENTS
      integer ind,inform,n
      double precision macheps

C     ARRAY ARGUMENTS
      integer indjac(n)
      double precision g(n),gplus1(n),gplus2(n),valjac(n),x(n)

C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     N integer: Number of variables.
C     ---------
C
C     X double precision x(n) 
C     -----------------------
C
C     Point at which the Hessian of the ind-th constraint will be 
C     tested.
C
C     IND integer
C     -----------
C
C     Index of the constraint whose gradient will be tested.
C
C     INDJAC integer indjac(n): N-dimensional integer working vector.
C     ------------------------
C
C     G      double precision g(n)
C     GPLUS1 double precision gplus1(n)
C     GPLUS2 double precision gplus2(n)
C     VALJAC double precision valjac(n)
C     ---------------------------------
C
C     N-dimensional double precision working vectors.
C
C     MACHEPS double precision
C     ------------------------
C
C     Smallest positive number such that 1 + macheps is not equal to 1.
C
C     On Return:
C     ==========
C
C     INFORM integer: Output flag.
C     --------------

C     PARAMETERS
      integer nsmax
      parameter ( nsmax = 1000 )

C     LOCAL SCALARS
      logical nullcol
      integer flag,i,j,nnzh,nnzjac
      double precision elem,hdiff1,hdiff2,step1,step2,tmp

C     LOCAL ARRAYS
      integer hlin(nsmax**2),hcol(nsmax**2)
      double precision H(nsmax,nsmax),hval(nsmax**2),maxerr(nsmax)

      inform = 0

C     Test viability

      if ( n .gt. nsmax ) then
          write(*,*) 'Subroutine CheckHc uses dense matrices up to ',
     +               'dimension ',nsmax,' times ',nsmax,'. The ',
     +               'Hessian checking will be skipped.'
          go to 500
      end if

C     Compute the constraint gradient at x and save in a dense vector

      call evaljac(n,x,ind,indjac,valjac,nnzjac,flag)
      if ( flag .ne. 0 ) then
          inform = - 93
          return
      end if

      do i = 1,n
          g(i) = 0.0d0
      end do

      do i = 1,nnzjac
          g(indjac(i)) = g(indjac(i)) + valjac(i)
      end do

C     Compute the Hessian of the constraints at x and save in 
C     dense matrix

      call evalhc(n,x,ind,hlin,hcol,hval,nnzh,flag)
      if ( flag .ne. 0 ) then
          inform = - 95
          return
      end if

      do j = 1,n
          do i = 1,n
              H(i,j) = 0.0d0
          end do
      end do

      do i = 1,nnzh
          H(hlin(i),hcol(i)) = H(hlin(i),hcol(i)) + hval(i)
      end do

      write(*,100) ind

      do j = 1,n

          tmp  = x(j)

C         Compute the constraint gradient at xplus1 and save in a 
C         dense vector

          step1 = sqrt( macheps ) * max( abs( tmp ), 1.0d0 )

          x(j) = tmp + step1

          call evaljac(n,x,ind,indjac,valjac,nnzjac,flag)
          if ( flag .ne. 0 ) then
              inform = - 93
              return
          end if

          do i = 1,n
              gplus1(i) = 0.0d0
          end do

          do i = 1,nnzjac
              gplus1(indjac(i)) = valjac(i)
          end do

C         Compute the constraint gradient at xplus2 and save in a 
C         dense vector

          step2 = sqrt( macheps ) * max( abs( tmp ), 1.0d-03 )

          x(j) = tmp + step2

          call evaljac(n,x,ind,indjac,valjac,nnzjac,flag)
          if ( flag .ne. 0 ) then
              inform = - 93
              return
          end if

          do i = 1,n
              gplus2(i) = 0.0d0
          end do

          do i = 1,nnzjac
              gplus2(indjac(i)) = valjac(i)
          end do

          x(j) = tmp

          write(*,105) j

          maxerr(j) = 0.0d0

          nullcol = .true.

          do i = 1,n
              if ( i .ge. j ) then
                  elem = H(i,j)
              else
                  elem = H(j,i)
              end if
              hdiff1 = ( gplus1(i) - g(i) ) / step1
              hdiff2 = ( gplus2(i) - g(i) ) / step2
              tmp = min( abs( elem - hdiff1 ), abs( elem - hdiff2 ) )
              if ( elem .ne. 0.0d0 .or.
     +             hdiff1 .ne. 0.0d0 .or. hdiff2 .ne. 0.0d0 ) then
                  if ( nullcol ) then
                      nullcol = .false.
                      write(*,106)
                  end if
                  write(*,110) i,elem,hdiff1,hdiff2,tmp
              end if
              maxerr(j) = max( maxerr(j), tmp )
          end do

          if ( nullcol ) then
              write(*,115)
          else
              write(*,120) maxerr(j)
          end if

      end do

      write(*,*)
      do j = 1,n
          write(*,130) j,maxerr(j)
      end do

 500  continue

      return

 100  format(/,1X,'Hessian matrix of constraint ',I5,' column by ',
     +            'column.')
 105  format(/,1X,'Column:  ',I6)
 106  format(/,1X,'Index',12X,'evalhc',3X,'Incr. Quoc. (two different ',
     +            'steps)',4X,'Absolute error')
 110  format(  1X,I5,4(3X,1P,D15.8))
 115  format(  1X,'All the elements of this column are null.')
 120  format(  1X,'Maximum absolute error = ',1P,D15.8)
 130  format(  1X,'Column ',I6,' Maximum absolute error = ',1P,D15.8)

      end

c#endif

C     =================================================================
C     Module: Augmented Lagrangian subroutines
C     =================================================================

C     Last update of any of the component of this module: 
C
C     January 30, 2007.

C     ******************************************************************
C     ******************************************************************

      subroutine easyalgencan(n,x,l,u,m,lambda,equatn,linear,rhoauto,
     +rhotype,rhomult,rhofrac,rho,gtype,hptype,intype,precond,epsfeas,
     +epsopt,maxoutit,maxtotit,maxtotfc,iprint,ncomp,macheps,bignum,f,
     +snorm,nalpsupn,outiter,totiter,totfcnt,totgcnt,totcgcnt,inform,
     +wi1,wd1,wd2,wd3,wd4,wd5,wd6,wd7,wd8,wd9,wd10,wd11,wd12,wd13,wd14,
     +wd15,wd16,wd17,wd18,wd19)

      implicit none

C     This subroutine aims to simplify the usage of algencan.

C     SCALAR ARGUMENTS
      logical rhoauto
      character * 6 precond
      integer gtype,hptype,inform,intype,iprint,m,maxoutit,maxtotfc,
     +        maxtotit,n,ncomp,outiter,rhotype,totcgcnt,totfcnt,totgcnt,
     +        totiter
      double precision bignum,epsfeas,epsopt,f,macheps,nalpsupn,rhomult,
     +        rhofrac,snorm

C     ARRAY ARGUMENTS
      integer wi1(n)
      logical equatn(m),linear(m)
      double precision l(n),lambda(m),rho(m),u(n),wd1(m),wd2(m),wd3(n),
     +        wd4(m),wd5(m),wd6(n),wd7(n),wd8(n),wd9(n),wd10(n),wd11(n),
     +        wd12(n),wd13(n),wd14(n),wd15(n),wd16(n),wd17(n),wd18(n),
     +        wd19(n),x(n)

C     LOCAL SCALARS
      integer maxitncp
      double precision epsopfin,epsopini,epsopcon,lammin,lammax

C     SET AUGMENTED-LAGRANGIAN PARAMETERS

      lammin    = - 1.0d+20         
      lammax    =   1.0d+20         

      epsopini  =    epsopt
      epsopfin  =    epsopt
      epsopcon  =   1.0d-01

      maxitncp  =         9

C     CALL ALGENCAN

      call algencan(n,x,l,u,m,lambda,equatn,linear,rhoauto,rhotype,
     +rhomult,rhofrac,rho,lammin,lammax,gtype,hptype,intype,precond,
     +epsfeas,epsopini,epsopfin,epsopcon,maxitncp,maxoutit,maxtotfc,
     +maxtotit,iprint,ncomp,macheps,bignum,f,wd1,wd2,snorm,wd3,nalpsupn,
     +outiter,totiter,totfcnt,totgcnt,totcgcnt,inform,wd4,wd5,wi1,wd6,
     +wd7,wd8,wd9,wd10,wd11,wd12,wd13,wd14,wd15,wd16,wd17,wd18,wd19)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine algencan(n,x,l,u,m,lambda,equatn,linear,rhoauto,
     +rhotype,rhomult,rhofrac,rho,lammin,lammax,gtype,hptype,intype,
     +precond,epsfeas,epsopini,epsopfin,epsopcon,maxitncp,maxoutit,
     +maxtotfc,maxtotit,iprint,ncomp,macheps,bignum,f,c,sigma,snorm,nal,
     +nalpsupn,outiter,totiter,totfcnt,totgcnt,totcgcnt,inform,lambar,
     +sigpre,wi1,wd1,wd2,wd3,wd4,wd5,wd6,wd7,wd8,wd9,wd10,wd11,wd12,
     +wd13,wd14)

      implicit none

C     Solves the general-constrained minimization problem
C
C     min f(x)
C
C     subject to
C
C             c_j(x)  = 0, j in E,
C             c_j(x) <= 0, j in I,
C             l <= x <= u,
C
C     where E is the set of indices of the equality constraints, I is
C     the set of indices of the inequality constraints, and there are
C     n variables and m constraints, using the method of multipliers 
C     described in:
C 
C     R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt, 
C     "On Augmented Lagrangian methods with general lower-level 
C     constraints", Technical Report MCDO-050304, Department of Applied 
C     Mathematics, UNICAMP, Brazil, 2005.

C     SCALAR ARGUMENTS
      logical rhoauto
      character * 6 precond
      integer gtype,hptype,inform,intype,iprint,m,maxitncp,maxoutit,
     +        maxtotfc,maxtotit,n,ncomp,outiter,rhotype,totcgcnt,
     +        totfcnt,totgcnt,totiter
      double precision bignum,epsfeas,epsopcon,epsopfin,epsopini,f,
     +        lammax,lammin,macheps,nalpsupn,rhofrac,rhomult,snorm

C     ARRAY ARGUMENTS
      integer wi1(n)
      logical equatn(m),linear(m)
      double precision c(m),l(n),lambar(m),lambda(m),nal(n),rho(m),
     +        sigma(m),sigpre(m),u(n),wd1(n),wd2(n),wd3(n),wd4(n),
     +        wd5(n),wd6(n),wd7(n),wd8(n),wd9(n),wd10(n),wd11(n),
     +        wd12(n),wd13(n),wd14(n),x(n)

C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     N integer: Number of variables.
C     ---------
C
C     X double precision x(n): Initial estimation of the solution.
C     -----------------------
C
C     L double precision l(n): Lower bounds for x.
C     -----------------------
C
C     U double precision u(n): Upper bounds for x.
C     -----------------------
C
C     M integer
C     ---------
C
C     Number of constraints (equalities plus inequalities) without 
C     considering the bound constraints.
C
C     LAMBDA double precision lambda(m)
C     ---------------------------------
C
C     Initial estimation of the Lagrange multipliers.
C
C     RHO double precision rho(m): Initial penalty parameters.
C     ---------------------------
C
C     EQUATN logical equatn(m)
C     ------------------------
C
C     Logical array that, for each constraint, indicates whether the 
C     constraint is an equality constraint (TRUE) or an inequality 
C     constraint (FALSE).
C
C     LINEAR logical linear(m)
C     ------------------------
C
C     Logical array that, for each constraint, indicates whether the 
C     constraint is a linear constraint (TRUE) or a nonlinear constraint 
C     (FALSE).
C
C     GTYPE integer
C     -------------
C
C     Type of derivatives calculation according to the following 
C     convention:
C
C     0 means true first derivatives. In this case, subroutines evalg
C       and evaljac must be modified by the user to compute the 
C       gradient of the objective function and the gradients of the 
C       constraints, respectively.
C
C     1 means that a finite difference approximation will be used. In 
C       this case, subroutines evalg and evaljac may have an empty body 
C       but must be present. It is also recommended that those 
C       empty-body subroutines set flag = - 1. Last but not least, the 
C       option gtype = 1 is not cheap neither safe.
C
C     HPTYPE integer
C     --------------
C
C     Type of Hessian-vector product according to the following 
C     convention:
C
C     We will first describe the meaning if the choices in the Augmented 
C     Lagrangian framework. The way in which the product of the Hessian 
C     matrix of the Augmented Lagrangian by a vector will be done 
C     depends on the value of the parameter hptype in the following way:
C
C     9 means that an incremental quotients approximation without any 
C       extra consideration will be used. This option requires the 
C       evaluation of an extra gradient at each Conjugate Gradient 
C       iteration. If gtype = 0 then this gradient evaluation will be 
C       done using the user supplied subroutines evalg and evaljac 
C       (evalc will also be used). On the other hand, if gtype = 1, the 
C       gradient calculation will be done using just calls to the user 
C       provided subroutines evalf and evalc. nind calls will be done, 
C       where nind is the dimension of the current face of the 
C       active-set method. This option is not cheap neither safe.
C
C       If you did not code subroutines evalg and evaljac, to compute 
C       the gradient of the objective function and the Jacobian of the 
C       constraints then your options finished here.
C
C     0 means that subroutines to compute the Hessian of the objective 
C       function (evalh) and the Hessians of the constraints (evalhc) 
C       were provided by the user. So, the product of the Hessian of the 
C       Augmented Lagrangian times a vector will be computed using the
C       Hessians provided by these subroutines and then adding the first 
C       order term (for the first order term the user-supplied 
C       subroutine to compute the Jacobian of the constraints (evaljac) 
C       is also used).
C
C     1 means that, instead of providing individual subroutines to 
C       compute the Hessians of the objective function and the 
C       constraints, the user provided a subroutine to compute the 
C       product of an arbitrary vector times the Hessian of the 
C       Lagrangian.
C
C     2 means that incremental quotients will be used. The difference 
C       between hptype = 9 and hptype = 2 is that, in the latter case, 
C       the non-differentiability of the Hessian of the Augmented 
C       Lagrangian will be taken into account. In particular, when 
C       computing the gradient of the Augmented Lagrangian at 
C       (x + step p), the constraints that will be considered will be 
C       the same constraints that contributed to the computation of the
C       gradient of the Augmented Lagrangian at x.
C
C       If GENCAN is been used to solve a bound-constrained problem 
C       which is not the subproblem of an Augmented Lagrangian method, 
C       but an isolated bound-constrained problem, then there is no 
C       difference between this option and hptype = 9.
C
C       This option also requires the evaluation of an extra gradient at 
C       each Conjugate Gradient iteration. 
C
C     3 is similar to hptype = 2. The difference is that the 
C       contribution of the linear constraints in the Hessian matrix 
C       will be computed explicitly. If the problem has not linear 
C       constraints then this option is identical to hptype = 2. 
C       Moreover, if the problem has no constraints then this option is 
C       equal to hptype = 9.
C
C       This option also requires the evaluation of an extra gradient at 
C       each Conjugate Gradient iteration. 
C
C     4 means that the Hessian matrix will be approximated and then the 
C       product of the Hessian approximation by the vector will be 
C       computed exactly. In particular, the Hessian matrix will be 
C       approximated doing a BFGS correction to the Gauss-Newton
C       approximation of the Hessian. Before the BFGS correction, a 
C       structured spectral correction is done to force the Gauss-Newton 
C       approximation to be positive definite.
C
C       If the problem has not constraints then the approximation 
C       reduces to a BFGS approximation of the Hessian (without memory) 
C       and using the spectral approximation (instead of the identity) 
C       as initial approximation.
C
C       Numerical experiments suggested that this option is convenient 
C       just for constrained problems. This motivated the introduction 
C       of the next option.
C
C       This option does NOT require an extra gradient evaluation per 
C       iteration and, in this sense, each CG iteration is 
C       computationally cheaper than a CG iteration of the previous 
C       choices. However, the approximation of the Hessian matrix 
C       requires some information (mainly the Jacobian of the 
C       constraints) that must be saved during the gradient evaluation. 
C       To save this information requires an amount of memory 
C       proportional to the number of non-null elements of the Jacobian 
C       matrix.
C 
C       Quadratic subproblems are convex with this choice.
C
C     5 is an adaptive strategy that choose, at every iteration, between
C       2 and 4. When the gradient of the Augmented Lagrangian is 
C       computed, it is verified if at least a constraint contributes to 
C       the calculation. If this is the case, 4 is used. Otherwise, 2 is 
C       used.
C
C       For problems with equality constraints (that always contributes 
C       to the Augmented Lagrangian function) this option is identical 
C       to 4.
C
C       For problems without constraints this option is identical to 2.
C
C     6 is identical to 5 but the choice is made between 3 and 4 instead 
C       of between 2 and 4. 
C
C       For problems with equality constraints (that always contributes 
C       to the Augmented Lagrangian function) this option is identical 
C       to 4. 
C
C       For problems without constraints this option is identical to 3.
C
C     We will now describe the meaning if the choices for unconstrained
C     and bound-constrained problems. In this context the way in which 
C     the product of the Hessian matrix by a vector will be done depends 
C     on the value of the parameter hptype in the following way:
C
C     0 means that the subroutine to compute the Hessian of the 
C       objective function (evalh) was provided by the user. So, the 
C       product of the Hessian times a vector will be computed using the 
C       Hessian provided by this subroutine.
C
C     1 means that a subroutine (evalhlp) to compute the product of the 
C       Hessian of the objective function times an arbitrary vector is
C       being provided by the user.
C
C     9 means that an incremental quotients approximation will be used. 
C       This option requires the evaluation of an extra gradient at each
C       Conjugate Gradient iteration. If gtype = 0 then this gradient 
C       evaluation will be done using the user supplied subroutine 
C       evalg. On the other hand, if gtype = 1, the gradient calculation 
C       will be done using just calls to the user provided subroutine 
C       evalf. nind calls will be done, where nind is the dimension of 
C       the current face of the active-set method.
C
C       If you did not code subroutine evalg to compute the gradient of 
C       the objective function then your options finished here.
C
C     4 means that the Hessian matrix will be approximated and then the 
C       product of the Hessian approximation by the vector will be 
C       computed exactly. The approximation is a BFGS approximation of 
C       the Hessian (without memory) and using the spectral 
C       approximation (instead of the identity) as initial 
C       approximation.
C
C       Numerical experiments suggested that this option is not 
C       convenient for unconstrained or just bound-constrained problems. 
C       (Note that this option was developed to be used in the Augmented 
C       Lagrangian framework.)
C
C       This option does NOT require an extra gradient evaluation per 
C       iteration and, in this sense, each CG iteration is 
C       computationally cheaper than a CG iteration of the previous 
C       choices.
C
C       Quadratic subproblems are convex with this choice.
C
C     In the bound-constrained context, options hptype = 2,3,5 and 6
C     are all identical to hptype = 9.
C
C     PRECOND character * 6
C     ---------------------
C
C     Indicates the type of preconditioning that will be used for 
C     Conjugates Gradients accordig to the followng convention:
C
C     'NONE'   means no preconditioner at all.
C
C     'QNAGNC' means Quasi-Newton Correction of the Gauss-Newton 
C              approximation of the Hessian. The exact form is this 
C              preconditioner is described in:
C 
C              E. G. Birgin and J. M. Martínez, "Structured minimal-
C              memory inexact quasi-Newton method and secant 
C              preconditioners for Augmented Lagrangian Optimization", 
C              submitted, 2005.
C
C     LAMMIN double precision: See below.
C     -----------------------
C
C     LAMMAX double precision
C     -----------------------
C
C     The Lagrange multipliers approximations used for the definition of 
C     the subproblems must be inside the interval [lammin,lammax].
C
C     RHOMULT double precision
C     ------------------------
C
C     When increased, the penalty parameters will be multiplied by 
C     rhomult.
C
C     RHOFRAC double precision
C     ------------------------
C
C     The penalty parameters will be, basically, increased when the 
C     infeasibility is larger than a fraction of the previous 
C     infeasibility. This fraction is rhofrac.
C
C     RHOTYPE integer
C     ---------------
C
C     Indicates the type of penalty parameters that will be used. 
C     rhotype = 1 means just a unique penalty parameter for all the 
C     constraints while rhotype = 2 means a penalty parameter per 
C     constraint.
C
C     EPSFEAS double precision 
C     ------------------------
C
C     Required precision for feasibility and complementarity.
C
C     EPSOPINI double precision
C     ------------------------- 
C
C     Required precision for the first subproblem optimality condition 
C     (sup-norm).
C
C     EPSOPFIN double precision: Required precision for optimality.
C     ------------------------- 
C
C     EPSOPCON double precision 
C     -------------------------
C
C     The precision used for each subproblem is allways computed as 
C     epsopk = max( epsopfin, epsopcon * epsop_k-1 ).
C
C     MAXITNCP integer 
C     ----------------
C
C     Maximum number of outer iterations withoutprogress in feasibility.
C
C     MAXOUTIT integer
C     ----------------
C
C     Maximum number of Augmented Lagrangian (outer) iterations.
C     (Ignored in the unconstrained and bound-constrained cases.)
C
C     MAXTOTIT integer
C     ----------------
C
C     Maximum total number of inner iterations in the Augmented 
C     Lagrangian context (total means summing up the inner iterations of 
C     each outer iteration). In the unconstrained and bound-constrained
C     cases it means just the maximum number of iterations.
C
C     MAXTOTFC integer
C     ----------------
C
C     Idem MAXTOTIT but for number of functional evaluations.
C
C     IPRINT integer
C     --------------
C                
C     Controls the ammount of information of the output according to the 
C     following convention:
C
C     0 means no output at all.
C
C     1 means information at each outer iteration but without any 
C       information of how the subproblems are being solved.
C
C     2 means the same as 1 plus information of each inner iteration.
C
C     3 means the same as 2 plus information of the line searches and 
C       the calculation of the truncated Newton direction (using CG) of 
C       each inner iteration.
C
C     In all cases, an output file named solution.txt with the final 
C     point, Lagrange mutipliers and penalty parameters will be 
C     generated. Moreover, the same output of the screen will be saved 
C     in a file named algencan.out.
C
C     NCOMP integer
C     -------------
C
C     Every time a vector is printed, just its first ncomp component 
C     will be displayed.
C
C     MACHEPS double precision
C     ------------------------
C
C     Smallest positive number such that 1 + macheps is not equal to 1.
C
C     BIGNUM double precision: A big number like 1.0d+99.
C     -----------------------
C
C     LAMBAR and SIGPRE double precision lambar(m),sigpre(m)
C     ------------------------------------------------------
C
C     M-dimensional double precision working vectors.
C
C     WI1 integer wi1(n): N-dimensional integer working vector.
C     ------------------
C
C     WD1 ... WD14 double precision wd1(n) ... wd14(n)
C     ------------------------------------------------
C
C     N-dimensional double precision working vectors.
C              
C     On Return:
C     ==========
C
C     X double precision x(n): Final estimation of the solution.
C     -----------------------
C
C     LAMBDA double precision lambda(m)
C     ---------------------------------
C
C     Final estimation of the Lagrange multipliers.
C
C     RHO double precision rho(m): Final penalty parameters.
C     ---------------------------
C
C     F double precision: Objective functional value at the solution.
C     ------------------
C
C     C double precision c(m): Constraints at the solution.
C     -----------------------
C
C     SIGMA double precision sigma(m)
C     -------------------------------
C
C     Complementarity measurement at the solution: c(i) for equality 
C     constraints and max(c(i), - lambar(i) / rho(i)) for inequalities.
C
C     SNORM double precision: Sup-norm of the constraints.
C     ----------------------
C
C     NAL double precision nal(n)
C     ---------------------------
C
C     Gradient of the Augmented Lagrangian at the solution.
C
C     NALPSUPN double precision
C     -------------------------
C 
C     Sup-norm of the continuous projected gradient of the Augmented 
C     Lagrangian at the solution
C
C     OUTITER integer: Number of outer iterations.
C     ---------------
C 
C     TOTITER integer: Total number of inner iterations.
C     ---------------
C 
C     TOTFCNT integer
C     ---------------
C 
C     Total number of Augmented Lagrangian function evaluations.
C
C     TOTGCNT integer
C     ---------------
C 
C     Total number of Augmented Lagrangian gradient evaluations.
C
C     TOTCGCNT integer
C     ----------------
C 
C     Total number of Conjugate Gradients iterations.
C
C     INFORM integer
C     --------------
C
C     This output parameter tells what happened in this subroutine, 
C     according to the following conventions:
C 
C       0 means convergence with feasibility, optimality and 
C         complementarity.
C
C       1 means that it was achieved the maximum number of outer 
C         iterations (maxoutit).
C
C       2 means that it was achieved the maximum total number of inner 
C         iterations (maxtotit).
C
C       3 means that it was achieved the maximum total number of 
C         functional evaluations (maxtotfc).
C
C       4 means that the algorithm stopped by ``lack of feasibility 
C         progress'', i.e., the current point is infeasible (the 
C         constraints violation is larger than the tolerance epsfeas) 
C         and the constraint violation has not even simple decrease 
C         during maxitncp consecutive iterations. In this case, the 
C         problem may be infeasible.
C
C     -90 means that subroutine evalf   retuned an error flag.
C
C     -91 means that subroutine evalc   retuned an error flag.
C
C     -92 means that subroutine evalg   retuned an error flag.
C
C     -93 means that subroutine evaljac retuned an error flag.
C
C     -94 means that subroutine evalh   retuned an error flag.
C 
C     -95 means that subroutine evalhc  retuned an error flag.
C 
C     -96 means that subroutine evalhlp retuned an error flag.

C     LOCAL SCALARS
      integer cgcnt,fcnt,gcnt,i,icrit,iter,mprint,nprint,maxit,maxfc,
     +        rhoincr
      double precision al,deafea,epsopk,nalpi,rhoini,rhomax,snormprev,
     +        sumc

C     SET UP SOME CONSTANTS

C     just for printing
      nprint   = min0(n,ncomp)
      mprint   = min0(m,ncomp)

C     INTIALIZATION

      icrit    = 0
      outiter  = 0
      totiter  = 0
      totfcnt  = 0
      totgcnt  = 0
      totcgcnt = 0

C     COMPUTE OBJECTIVE FUNCTION AND CONSTRAINTS

      call evalobjc(n,x,f,m,c,inform)
      totfcnt = totfcnt + 1

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 1 ) then
              write(*, 2000) inform
              write(10,2000) inform
          end if

          go to 500
      end if

C     SET PENALTY PARAMETERS:

C     AUGMENTED LAGRANGIAN FUNCTION IS OF THE FORM
C
C         AL = F(X) + \SUM_i P(LAMBDA_i,C_i(X),RHO),
C
C     WHERE
C
C         P(LAMBDA,Y,RHO) = Y ( LAMBDA + 0.5 RHO Y ), 
C
C     IF C_i(X) IS AN EQUALITY CONSTRAINT OR LAMBDA + RHO Y > 0, AND
C
C         P(LAMBDA,Y,RHO) = - 0.5 LAMBDA^2 / RHO,
C
C     OTHERWISE.
C
C     ASSUMING THAT LAMBDA_i = 0 FOR ALL i, IT IS CLEAR THAT 
C
C         P(LAMBDA_i,C_i(X),RHO) = 0.5 RHO C_i(X)^2
C
C     AND THAT THE VALUE OF RHO THAT BALANCES F(X) AND 
C     \SUM_i P(LAMBDA_i,C_i(X),RHO) IS GIVEN BY
C
C     RHO = F(X) / ( 0.5 \SUM_i C_i(X)^2 ),
C
C     WHERE THE SUM IS MADE OVER THE VIOLATED CONSTRAINTS.

      if ( rhoauto ) then

          sumc = 0.0d0
          do i = 1,m
              if ( equatn(i) .or. c(i) .gt. 0.0d0 ) then
                  sumc = sumc + 0.5d0 * c(i) ** 2
              end if
          end do

          if ( sumc .eq. 0.0d0 ) then
              rhoini = 10.0d0
          else
              rhoini = 10.0d0 * max( 1.0d0, abs( f ) ) / sumc
              rhoini = max( 1.0d-06, min( rhoini, 10.0d0 ) )
          end if 

          do i = 1,m
              rho(i) = rhoini
          end do

      end if

C     COMPUTE MAXIMUM RHO

      rhomax = 0.0d0
      do i = 1,m
          rhomax = max( rhomax, rho(i) )
      end do

C     COMPUTE LAGRANGE MULTIPLIERS

      do i = 1,m
          lambar(i) = max( lammin, min( lammax, lambda(i) ) )
      end do

C     COMPUTE COMPLEMENTARITY AND FEASIBILITY VIOLATIONS
 
      snormprev = bignum

      snorm = 0.0d0
      do i = 1,m
          if ( equatn(i) ) then
              sigma(i) = c(i)
          else
              sigma(i) = max( c(i), - lambar(i) / rho(i) )
          end if
          snorm = max( snorm, abs( sigma(i) ) )
      end do

C     COMPUTE CONTINUOUS PROJECTED GRADIENT NORM

      call setpoint(x)

      call evalnl(n,x,m,lambar,equatn,linear,nal,gtype,macheps,inform)
      totgcnt = totgcnt + 1

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 1 ) then
              write(*, 2000) inform
              write(10,2000) inform
          end if

          go to 500
      end if

      nalpsupn = 0.0d0
      do i = 1,n
          nalpi    = max( l(i), min( x(i) - nal(i), u(i) ) ) - x(i)
          nalpsupn = max( nalpsupn, abs( nalpi ) )
      end do

C     PRINT INITIAL INFORMATION

      if  ( iprint .ge. 1 ) then
          write(*,  900) n,m
          write(*,  910) nprint,(l(i),i=1,nprint)
          write(*,  920) nprint,(u(i),i=1,nprint)

          write(10, 900) n,m
          write(10, 910) nprint,(l(i),i=1,nprint)
          write(10, 920) nprint,(u(i),i=1,nprint)
      end if

C     MAIN LOOP

 100  continue

C     PRINT ITERATION INFORMATION

      if ( iprint .ge. 1 ) then
          write(*,  970) outiter
          write(*,  980) nprint,(x(i),i=1,nprint)
          write(*,  990) mprint,(lambda(i),i=1,mprint)
          write(*, 1000) mprint,(rho(i),i=1,mprint)
          write(*, 1020) f,snorm,nalpsupn,totiter,totfcnt,totgcnt

          write(10, 970) outiter
          write(10, 980) nprint,(x(i),i=1,nprint)
          write(10, 990) mprint,(lambda(i),i=1,mprint)
          write(10,1000) mprint,(rho(i),i=1,mprint)
          write(10,1020) f,snorm,nalpsupn,totiter,totfcnt,totgcnt
      end if

C     SAVING INTERMEDIATE DATA FOR CRASH REPORT

      open(20,file='algencan-tabline.out')
      write(20,3000) 0.0d0,-1,n,m,outiter,totiter,totfcnt,totgcnt,
     +totcgcnt,f,snorm,nalpsupn
      close(20)

C     TEST STOPPING CRITERIA

      if ( snorm .le. epsfeas .and. nalpsupn .le. epsopfin ) then

        ! THE POINT IS FEASIBLE AND OPTIMAL

          inform = 0

          if ( iprint .ge. 1 ) then
              write(*, 1060) inform,epsopfin,epsfeas,epsfeas
              write(10,1060) inform,epsopfin,epsfeas,epsfeas
          end if

          go to 500
      end if

      if ( outiter .ge. maxoutit ) then
        ! MAXIMUM NUMBER OF OUTER ITERATIONS EXCEEDED, STOP
          inform = 1

          if ( iprint .ge. 1 ) then
              write(*, 1070) inform,maxoutit
              write(10,1070) inform,maxoutit
          end if

          go to 500
      end if
 
      if ( totiter .gt. maxtotit ) then
        ! MAXIMUM NUMBER OF INNER ITERATIONS EXCEEDED, STOP
          inform = 2

          if ( iprint .ge. 1 ) then
              write(*, 1080) inform,maxtotit
              write(10,1080) inform,maxtotit
          end if

          go to 500
      end if
 
      if ( totfcnt .gt. maxtotfc ) then
        ! MAXIMUM NUMBER OF FUNCTION EVALUATIONS EXCEEDED, STOP
          inform = 3

          if ( iprint .ge. 1 ) then
              write(*, 1090) inform,maxtotfc
              write(10,1090) inform,maxtotfc
          end if

          go to 500
      end if

C     TEST IF THE INFEASIBILITY DECREASED SUFFICIENTLY

      if ( snorm .gt. epsfeas .and. snorm .ge. snormprev ) then

          icrit = icrit + 1

          if ( icrit .ge. maxitncp ) then
            ! LACK OF FEASIBILITY PROGRESS, STOP
              inform = 4

              if ( iprint .ge. 1 ) then
                  write(*, 1100) inform,maxitncp
                  write(10,1100) inform,maxitncp
              end if

              go to 500
          end if

      else
          icrit = 0
      end if

C     TEST IF THE PENALTY PARAMETERS ARE TOO LARGE

      if ( rhomax .gt. 1.0d+20 ) then
        ! TOO LARGE PENALTY PARAMETERS
          inform = 5

          if ( iprint .ge. 1 ) then
              write(*, 1110) inform,rhomax
              write(10,1110) inform,rhomax
          end if

          go to 500

      end if

C     DO AN OUTER ITERATION

      outiter = outiter + 1      

C     OPTIMALITY REQUERIMENT FOR THE SUBPROBLEM

      if ( outiter .eq. 1 ) then
          epsopk = epsopini
      else
          epsopk = max( epsopfin, epsopk * epsopcon )
      end if

      maxit =        100
      maxfc = 10 * maxit

C     CALL THE INNER-SOLVER

      call easygencan(n,x,l,u,m,lambar,equatn,linear,rho,gtype,hptype,
     +intype,precond,epsopk,maxit,maxfc,iprint,ncomp,macheps,bignum,al,
     +nal,nalpsupn,iter,fcnt,gcnt,cgcnt,inform,wi1,wd1,wd2,wd3,wd4,wd5,
     +wd6,wd7,wd8,wd9,wd10,wd11,wd12,wd13,wd14)

      totiter  = totiter  + iter
      totfcnt  = totfcnt  + fcnt
      totgcnt  = totgcnt  + gcnt
      totcgcnt = totcgcnt + cgcnt

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 1 ) then
              write(*, 2000) inform
              write(10,2000) inform
          end if

          go to 500
      end if

C     COMPUTE OBJECTIVE FUNCTION AND CONSTRAINTS

      call evalobjc(n,x,f,m,c,inform)
      totfcnt = totfcnt + 1

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 1 ) then
              write(*, 2000) inform
              write(10,2000) inform
          end if

          go to 500
      end if

C     UPDATE LAGRANGE MULTIPLIERS APPROXIMATION 

      do i = 1,m
          call evaldpdy(c(i),rho(i),lambar(i),equatn(i),lambda(i))
          lambar(i) = max( lammin, min( lammax, lambda(i) ) )
      end do

C     COMPUTE COMPLEMENTARITY AND FEASIBILITY VIOLATIONS
 
      snormprev = snorm

      snorm = 0.0d0
      do i = 1,m
          sigpre(i) = sigma(i)

          if ( equatn(i) ) then
              sigma(i) = c(i)
          else
              sigma(i) = max( c(i), - lambar(i) / rho(i) )
          end if

          snorm = max( snorm, abs( sigma(i) ) )
      end do

C     UPDATE PENALTY PARAMETERS
 
      if ( outiter .eq. 1 ) then

          if ( iprint .ge. 1 ) then
              write(*, 1025)
              write(10,1025)
          end if

      else

          if ( rhotype .eq. 1 ) then

              deafea = max( epsfeas, rhofrac * snormprev )

              if ( snorm .gt. deafea ) then

                  do i = 1,m
                      rho(i) = rhomult * rho(i)
                  end do

                  rhomax = rho(1)

                  if ( iprint .ge. 1 ) then
                      write(*, 1030) deafea,rhomult
                      write(10,1030) deafea,rhomult
                  end if

              else

                  if ( iprint .ge. 1 ) then
                      write(*, 1050)
                      write(10,1050)
                  end if

              end if

          else ! if ( rhotype .eq. 2 ) then

              rhoincr = 0
              do i = 1,m
                  if ( abs(sigma(i)) .gt.
     +            max ( epsfeas, rhofrac * abs(sigpre(i)) ) ) then
                      rho(i)  = rhomult * rho(i)
                      rhomax  = max( rhomax, rho(i) )
                      rhoincr = rhoincr + 1
                  end if
              end do

              if ( iprint .ge. 1 ) then

                  if ( rhoincr .ne. 0 ) then
                      write(*, 1040) rhoincr,m,rhomult
                      write(10,1040) rhoincr,m,rhomult
                  else
                      write(*, 1050)
                      write(10,1050)
                  end if

              end if

          end if

      end if

C     ITERATE

      go to 100

C     STOP

 500  continue

C     NON-EXECUTABLE STATEMENTS 

  900 format( /,' Entry to ALGENCAN.',/,
     +        /,' Number of variables: ',I7,
     +          ' Number of constraints: ',I7)
  910 format( /,' Lower bounds (first ',I7, ' components): ',
     +        /,6(1X,1PD11.4))
  920 format( /,' Upper bounds (first ',I7, ' components): ',
     +        /,6(1X,1PD11.4))
  970 format( /,' ALGENCAN outer iteration: ',I7)
  980 format( /,' Current point (first ',I7, ' components): ',
     +        /,6(1X,1PD11.4))
  990 format( /,' Updated Lagrange multipliers (first ',I7,
     +          ' components): ',
     +        /,6(1X,1PD11.4))
 1000 format( /,' Updated penalty parameters (first ',I7,
     +          ' components): ',
     +        /,6(1X,1PD11.4))
 1020 format( /,' Objective function value',38X,' = ',1PD11.4,
     +        /,' Sup-norm of the constraints',35X,' = ',1PD11.4,
     +        /,' Sup-norm of the projected gradient of the',
     +          ' Lagrangian',10X,' = ',1PD11.4,/,
     +        /,' Up-to-now number of inner iterations',26X,' = ',4X,I7,
     +        /,' Up-to-now number of Augmented Lagrangian function',
     +          ' evaluations',1X,' = ',4X,I7,
     +        /,' Up-to-now number of Augmented Lagrangian gradient',
     +          ' evaluations',1X,' = ',4X,I7)
 1025 format( /,' The penalty parameters are not modified',
     +        /,' after the resolution of the first subproblem.')
 1030 format( /,' The desired infeasibility was not achieved (',
     +            1PD11.4,'). The penalty',
     +        /,' parameter will be increased multiplying by rhofrac',
     +          ' = ',1PD11.4,'.')
 1040 format( /,' The desired infeasibility was not achieved in ',I7,
     +          ' over ',I7,
     +        /,' constraints.',
     +        /,' Penalty parameters will be increased multiplying',
     +          ' by rhofrac = ',1PD11.4,'.')
 1050 format( /,' Desired feasibility improvement was achieved,',
     +        /,' so, penalty parameters will not be modified.')
 1060 format( /,' Flag of ALGENCAN = ',I2,
     +          ' (Convergence with Sup-norm of the Lagrangian',
     +          ' projected',/,' gradient smaller than ',1PD11.4,',',
     +          ' constraints Sup-norm smaller than',
     +        /,' ',1PD11.4,', and largest complementarity violation',
     +          ' smaller than ',1PD11.4,')')
 1070 format( /,' Flag of ALGENCAN = ',I2,
     +          ' (It was exceeded the maximum allowed number of outer',
     +        /,' iterations (maxoutit =',I7,').)')
 1080 format( /,' Flag of ALGENCAN = ',I2,
     +          ' (It was exceeded the maximum allowed number of inner',
     +        /,' iterations (maxinnit =',I7,').)')
 1090 format( /,' Flag of ALGENCAN = ',I2,
     +          ' (It was exceeded the maximum allowed number of',
     +        /,' functional evaluations (maxfc =',I7,').)')
 1100 format( /,' Flag of ALGENCAN = ',I2,
     +          ' (Constraint violations have not decreased',
     +          ' substantially',/,' over',I2,' outer iterations.',
     +          ' Problem possibly infeasible.)')
 1110 format( /,' Flag of ALGENCAN = ',I2,
     +          ' (Too large penalty paremeters: ',1PD11.4,'.)')
 2000 format( /,' Flag of ALGENCAN = ',I3,' Fatal Error',/,
     +        /,' The following codes means: ',/,
     +        /,' -90 : error in evalf   subroutine', 
     +        /,' -91 : error in evalc   subroutine', 
     +        /,' -92 : error in evalg   subroutine', 
     +        /,' -93 : error in evaljac subroutine', 
     +        /,' -94 : error in evalh   subroutine', 
     +        /,' -95 : error in evalhc  subroutine', 
     +        /,' -96 : error in evalhlp subroutine',/) 
 3000 format(F8.2,1X,I3,1X,I6,1X,I6,1X,I3,1X,I7,1X,I7,1X,I7,1X,I7,1X,1P
     +       D24.16,1X,1P,D7.1,1X,1P,D7.1)

      end

c#ifndef CUTEr

C     *****************************************************************
C     *****************************************************************

      subroutine evalobjc(n,x,f,m,c,inform)

      implicit none

C     This subroutine computes the objective function and the 
C     constraints.

C     SCALAR ARGUMENTS
      double precision f
      integer inform,m,n

C     ARRAY ARGUMENTS
      double precision c(m),x(n)

C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     N integer: Number of variables.
C     ---------
C
C     X double precision x(n): Current point.
C     -----------------------
C
C     M integer: Number of constraints.
C     ---------
C
C     On Return:
C     ==========
C
C     F double precision: Objective function value at x.
C     ------------------
C
C     C double precision c(m): Constraints at x.
C     -----------------------
C
C     INFORM integer
C     --------------
C
C     0 means that the evaluation was successfuly done.
C
C     Any other value means that some error occured during the 
C     evaluation.

C     LOCAL SCALARS
      integer flag,i

      inform = 0

      call setpoint(x)

C     COMPUTE OBJECTIVE FUNTION

      call evalf(n,x,f,flag)
      if ( flag .ne. 0 ) then
          inform = - 90
          return  
      end if

C     COMPUTE CONSTRAINTS

      do i = 1,m

C         COMPUTE THE i-TH CONSTRAINT
          call evalc(n,x,i,c(i),flag)
          if ( flag .ne. 0 ) then
              inform = - 91
              return
          end if

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine evalal(n,x,m,lambda,rho,equatn,linear,al,inform)

      implicit none

C     This subroutine computes the Augmented Lagrangian function.

C     SCALAR ARGUMENTS
      double precision al
      integer inform,m,n

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),rho(m),x(n)

C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     N integer: Number of variables.
C     ---------
C
C     X double precision x(n): Current point.
C     -----------------------
C
C     M integer: Number of constraints.
C     ---------
C
C     LAMBDA double precision lambda(m)
C     ---------------------------------
C
C     Initial estimation of the Lagrange multipliers.
C
C     RHO double precision rho(m): Initial penalty parameters.
C     ---------------------------
C
C     EQUATN logical equatn(m)
C     ------------------------
C
C     Logical array that, for each constraint, indicates whether the 
C     constraint is an equality constraint (TRUE) or an inequality 
C     constraint (FALSE).
C
C     LINEAR logical linear(m)
C     ------------------------
C
C     Logical array that, for each constraint, indicates whether the 
C     constraint is a linear constraint (TRUE) or a nonlinear constraint 
C     (FALSE).
C
C     On Return:
C     ==========
C
C     AL double precision al: Value of the Augmented Lagrangian function.
C     ----------------------
C
C     INFORM integer
C     --------------
C
C     0 means that the evaluation was successfuly done.
C
C     Any other value means that some error occured during the 
C     evaluation.

C     LOCAL SCALARS
      integer flag,i
      double precision c,f,p

      inform = 0

      call setpoint(x)

C     COMPUTE OBJECTIVE FUNTION

      call evalf(n,x,f,flag)
      if ( flag .ne. 0 ) then
          inform = - 90
          return
      end if

C     COMPUTES AL(x) = f(x) + \sum_j=1^m P(c(j),rho(j),lambda(j))

      al = f

      do i = 1,m

C         COMPUTE i-TH CONSTRAINT
          call evalc(n,x,i,c,flag)
          if ( flag .ne. 0 ) then
              inform = - 91
              return
          end if

C         ADD P(c(i),rho(i),lambda(i))
          call evalp(c,rho(i),lambda(i),equatn(i),p)
          al = al + p

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine evalnl(n,x,m,lambda,equatn,linear,nl,gtype,macheps,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      integer gtype,inform,m,n
      double precision macheps

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nl(n),x(n)

C     This subroutine computes the gradient of the Lagrangian function.


C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     N integer: Number of variables.
C     ---------
C
C     X double precision x(n): Current point.
C     -----------------------
C
C     M integer: Number of constraints.
C     ---------
C
C     LAMBDA double precision lambda(m)
C     ---------------------------------
C
C     Lagrange multipliers.
C
C     EQUATN logical equatn(m)
C     ------------------------
C
C     Logical array that, for each constraint, indicates whether the 
C     constraint is an equality constraint (TRUE) or an inequality 
C     constraint (FALSE).
C
C     LINEAR logical linear(m)
C     ------------------------
C
C     Logical array that, for each constraint, indicates whether the 
C     constraint is a linear constraint (TRUE) or a nonlinear constraint 
C     (FALSE).
C
C     On Return:
C     ==========
C
C     NL double precision nl(n): 
C     --------------------------
C
C     Gradient of the Lagrangian function.
C
C     INFORM integer
C     --------------
C
C     0 means that the evaluation was successfuly done.
C
C     Any other value means that some error occured during the 
C     evaluation.

C     PARAMETERS
      integer mmax,nmax,jcnnzmax
      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( jcnnzmax  =  5000000 )

C     COMMON SCALARS
      logical constrc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     LOCAL SCALARS
      integer flag,i,ind,j

C     COMMON BLOCKS
      common /graddata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc

      inform = 0

C     COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION

      if ( gtype .eq. 0 ) then
          call evalg(n,x,nl,flag)
          if ( flag .ne. 0 ) then
              inform = - 92
              return
          end if
      else
          call evalgdiff(n,x,nl,macheps,inform)
          if ( inform .lt. 0 ) then
              return
          end if
      end if

      do i = 1,n
          g(i) = nl(i)
      end do

C     COMPUTE \nabla L(x) = \nabla f(x) + 
C                           \sum_{j=1}^m lambda_j * \nabla c_j(x)

      ind = 0

      constrc = .false.

      do j = 1,m

          if ( equatn(j) .or. lambda(j) .gt. 0.0d0 ) then

              jcsta(j) = ind + 1

C             COMPUTE THE GRADIENT OF THE j-TH CONSTRAINT
              if ( gtype .eq. 0 ) then
                  call evaljac(n,x,j,jcvar(ind+1),jcval(ind+1),jclen(j),
     +            flag)
                  if ( flag .ne. 0 ) then
                      inform = - 93
                      return
                  end if
              else
                  call evaljacdiff(n,x,j,jcvar(ind+1),jcval(ind+1),
     +            jclen(j),macheps,inform)
                  if ( inform .lt. 0 ) then
                      return
                  end if
              end if

              ind = ind + jclen(j)

C             ADD lambda_j * \nabla c_j(x)
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  nl(jcvar(i)) = nl(jcvar(i)) + lambda(j) * jcval(i)
              end do

              if ( .not. linear(j) ) then
                  do i = jcsta(j),jcsta(j) + jclen(j) - 1
                      g(jcvar(i)) = g(jcvar(i)) + lambda(j) * jcval(i)
                  end do
              end if

              constrc = .true.

          end if

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine evalnal(n,x,m,lambda,rho,equatn,linear,nal,gtype,
     +macheps,inform)

      implicit none

C     This subroutine computes the gradient of the Augmented Lagrangian 
C     function.

C     SCALAR ARGUMENTS
      integer gtype,inform,m,n
      double precision macheps

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nal(n),rho(m),x(n)

C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     N integer: Number of variables.
C     ---------
C
C     X double precision x(n): Current point.
C     -----------------------
C
C     M integer: Number of constraints.
C     ---------
C
C     LAMBDA double precision lambda(m)
C     ---------------------------------
C
C     Lagrange multipliers.
C
C     RHO double precision rho(m): Initial penalty parameters.
C     ---------------------------
C
C     EQUATN logical equatn(m)
C     ------------------------
C
C     Logical array that, for each constraint, indicates whether the 
C     constraint is an equality constraint (TRUE) or an inequality 
C     constraint (FALSE).
C
C     LINEAR logical linear(m)
C     ------------------------
C
C     Logical array that, for each constraint, indicates whether the 
C     constraint is a linear constraint (TRUE) or a nonlinear constraint 
C     (FALSE).
C
C     On Return:
C     ==========
C
C     NAL double precision nal(n): 
C     ---------------------------
C
C     Gradient of the Augmented Lagrangian function.
C
C     INFORM integer
C     --------------
C
C     0 means that the evaluation was successfuly done.
C
C     Any other value means that some error occured during the 
C     evaluation.

C     PARAMETERS
      integer mmax,nmax,jcnnzmax
      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( jcnnzmax  =  5000000 )

C     COMMON SCALARS
      logical constrc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     LOCAL SCALARS
      integer flag,j
      double precision c

C     COMMON BLOCKS
      common /graddata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc

      inform = 0

      do j = 1,m

C         COMPUTE THE j-TH CONSTRAINT
          call evalc(n,x,j,c,flag)
          if ( flag .ne. 0 ) then
              inform = - 91
              return
          end if

C         COMPUTE dP/dc
          call evaldpdy(c,rho(j),lambda(j),equatn(j),dpdc(j))

      end do

C     COMPUTE GRADIENT OF THE LAGRANGIAN WITH DPDC INSTEAD OF LAMBDA
      call evalnl(n,x,m,dpdc,equatn,linear,nal,gtype,macheps,inform)

      end

C     *****************************************************************
C     *****************************************************************

      subroutine evalnalu(n,xp,m,lambda,rho,equatn,linear,nalp,gtype,
     +macheps,inform)

      implicit none

C     This subroutine computes the gradient of the Augmented Lagrangian 
C     function at a point xp, which is near to x, taking care of the
C     non-differentiability. The Augmented Lagrangian gradient must be
C     previously computed at x.

C     SCALAR ARGUMENTS
      integer gtype,inform,m,n
      double precision macheps

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),nalp(n),rho(m),xp(n)

C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     N integer: Number of variables.
C     ---------
C
C     X double precision x(n): Current point.
C     -----------------------
C
C     M integer: Number of constraints.
C     ---------
C
C     LAMBDA double precision lambda(m)
C     ---------------------------------
C
C     Initial estimation of the Lagrange multipliers.
C
C     RHO double precision rho(m): Initial penalty parameters.
C     ---------------------------
C
C     EQUATN logical equatn(m)
C     ------------------------
C
C     Logical array that, for each constraint, indicates whether the 
C     constraint is an equality constraint (TRUE) or an inequality 
C     constraint (FALSE).
C
C     LINEAR logical linear(m)
C     ------------------------
C
C     Logical array that, for each constraint, indicates whether the 
C     constraint is a linear constraint (TRUE) or a nonlinear constraint 
C     (FALSE).
C
C     On Return:
C     ==========
C
C     AL double precision al: 
C     ----------------------
C
C     Gradient of the Augmented Lagrangian function.
C
C     INFORM integer
C     --------------
C
C     0 means that the evaluation was successfuly done.
C
C     Any other value means that some error occured during the 
C     evaluation.

C     PARAMETERS
      integer mmax,nmax,jcnnzmax
      parameter ( mmax      =   500000 )
      parameter ( nmax      =   500000 )
      parameter ( jcnnzmax  =  5000000 )

C     COMMON SCALARS
      logical constrc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     LOCAL SCALARS
      integer flag,i,j,jcpnnz
      double precision cp,dpdcp

C     LOCAL ARRAYS
      integer jcpvar(nmax)
      double precision jcpval(nmax)

C     COMMON BLOCKS
      common /graddata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc

      inform = 0

C     COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION AT xp

      if ( gtype .eq. 0 ) then
          call evalg(n,xp,nalp,flag)
          if ( flag .ne. 0 ) then
              inform = - 92
              return
          end if
      else
          call evalgdiff(n,xp,nalp,macheps,inform)
          if ( inform .lt. 0 ) then
              return
          end if
      end if

C     COMPUTE \nabla L(x) = \nabla f(x) + \sum_{i=1}^m dPdc * dcdx

      do j = 1,m

          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

C             COMPUTE THE i-TH CONSTRAINT
              call evalc(n,xp,j,cp,flag)
              if ( flag .ne. 0 ) then
                  inform = - 91
                  return
              end if

C             COMPUTE dP/dc
              dpdcp = lambda(j) + rho(j) * cp

              if ( dpdcp .ne. 0.0d0 ) then

C                 COMPUTE THE GRADIENT OF THE j-TH CONSTRAINT
                  if ( gtype .eq. 0 ) then
                      call evaljac(n,xp,j,jcpvar,jcpval,jcpnnz,flag)
                      if ( flag .ne. 0 ) then
                          inform = - 93
                          return
                      end if
                  else
                      call evaljacdiff(n,xp,j,jcpvar,jcpval,jcpnnz,
     +                macheps,inform)
                      if ( inform .lt. 0 ) then
                          return
                      end if
                  end if

C                 ADD dPdc * dcdx
                  do i = 1,jcpnnz
                      nalp(jcpvar(i)) = 
     +                nalp(jcpvar(i)) + dpdcp * jcpval(i)
                  end do

              end if

          end if

      end do

      end

C     *****************************************************************
C     *****************************************************************

      subroutine evalhalp(n,x,p,nal,m,lambda,rho,equatn,linear,s,y,
     +ssupn,seucn,yeucn,sts,sty,lspgmi,lspgma,samefa,gtype,hptype,
     +aptype,hp,xp,tmp,macheps,inform,goth,hlspg,hds,hstds)

      implicit none

C     SCALAR ARGUMENTS
      logical goth,samefa
      character * 6 aptype
      integer gtype,hptype,inform,m,n
      double precision hlspg,hstds,lspgma,lspgmi,macheps,seucn,ssupn,
     +        sts,sty,yeucn
 
C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision hds(n),hp(n),lambda(m),nal(n),p(n),rho(m),s(n),
     +        tmp(n),x(n),xp(n),y(n)

C     PARAMETERS
      integer mmax,nmax,hnnzmax,jcnnzmax
      parameter ( mmax      =  500000 )
      parameter ( nmax      =  500000 )
      parameter ( hnnzmax   = 5000000 )
      parameter ( jcnnzmax  = 5000000 )

C     COMMON SCALARS
      logical constrc

C     COMMON ARRAYS
      integer hcol(hnnzmax),hlen(0:mmax),hlin(hnnzmax),hsta(0:hnnzmax),
     +        jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision dpdc(mmax),g(nmax),hval(hnnzmax),jcval(jcnnzmax)

C     LOCAL SCALARS
      logical exalin
      integer i,ind,flag,j,jcpnnz
      double precision atp,c1,c2,cp,dpdcp,psupn,ptds,pty,step,xsupn

C     LOCAL ARRAYS
      integer jcpvar(nmax)
      double precision jcpval(nmax)

C     COMMON BLOCKS
      common /graddata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc
      common /hessdata/ hval,hlin,hcol,hsta,hlen

      inform = 0

C     ==================================================================
C     HESSIAN APPROXIMATION TYPE
C     ==================================================================

C     True hessian-vector product using the user-provided second 
C     derivatives of the objective function and the constraints
      if ( hptype .eq. 0 ) then
          aptype = 'TRUEHE'

C     For user user-provided Hessian-of-the-Lagrangian times vector 
C     subroutine
      else if ( hptype .eq. 1 ) then
          aptype = 'HLPROD'

C     Incremental quotients considering the non-differentiability of 
C     the Hessian matrix of the Augmented Lagrangian
      else if ( hptype .eq. 2 ) then
          aptype = 'INCQUO'
          exalin = .false.

C     Idem 2 but considering explicitly the contribution of the 
C     linear constraints in the Hessian matrix
      else if ( hptype .eq. 3 ) then
          aptype = 'INCQUO'
          exalin = .true.

C     Quasi-Newton correction of a Gauss-Newton approximation of the 
C     Hessian matrix
      else if ( hptype .eq. 4 ) then
          aptype = 'QNCGNA'

C     Adaptative choice: idem 4 when there is some contribution of 
C     the constraints in the Augmented Lagrangian function and idem 
C     2 when there is no contribution at all.
      else if ( hptype .eq. 5 ) then
          if ( constrc ) then
              aptype = 'QNCGNA'
          else
              aptype = 'INCQUO'
              exalin = .false.
          end if

C     Idem 5 but using 3 when there is no contribution of the
C     constraints in the Augmented Lagrangian function
      else if ( hptype .eq. 6 ) then
          if ( constrc ) then
              aptype = 'QNCGNA'
          else
              aptype = 'INCQUO'
              exalin = .true.
          end if

C     Incremental quotients without considering the non-differentiability 
C     of the Hessian matrix of the Augmented Lagrangian
      else ! if ( hptype .eq. 9 ) then
          aptype = 'PUREIQ'
      end if

      if ( aptype .eq. 'TRUEHE' ) then
          go to 100
      else if ( aptype .eq. 'HLPROD' ) then
          go to 200
      else if ( aptype .eq. 'QNCGNA' ) then
          go to 300
      else if ( aptype .eq. 'INCQUO' ) then
          go to 400
      else ! if ( aptype .eq. 'PUREIQ' ) then
          go to 500
      end if

C     ==================================================================
C     EXACT SECOND-DERIVATIVES OF THE OBJECTIVE FUNCTION AND THE 
C     CONSTRAINTS WHERE PROVIDED BY THE USER. THE TRUE HESSIAN-VECTOR 
C     PRODUCT WILL BE COMPUTED
C     ==================================================================

 100  continue

C     ------------------------------------------------------------------
C     COMPUTE THE HESSIAN IF THIS IS THE FIRST TIME THIS SUBBROUTINE IS 
C     BEING CALLED
C     ------------------------------------------------------------------

      if ( .not. goth ) then
          goth = .true.

C         COMPUTE THE HESSIAN OF THE OBJECTIVE FUNCTION
          hsta(0) = 1

          call evalh(n,x,hlin,hcol,hval,hlen(0),flag)
          if ( flag .ne. 0 ) then
              inform = - 94
              return
          end if

          ind = hlen(0)

C         COMPUTE THE HESSIANS OF THE CONSTRAINTS
          do j = 1,m
              if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then
                  hsta(j) = ind + 1

                  call evalhc(n,x,j,hlin(ind+1),hcol(ind+1),hval(ind+1),
     +            hlen(j),flag)
                  if ( flag .ne. 0 ) then
                      inform = - 95
                      return
                  end if

                  ind = ind + hlen(j)
              end if
          end do
      end if

C     ------------------------------------------------------------------
C     COMPUTE hp = \nabla^2 f(x) p
C     ------------------------------------------------------------------

      do i = 1,n
          hp(i) = 0.0d0
      end do

      do i = hsta(0),hsta(0) + hlen(0) - 1
          if ( hlin(i) .ge. hcol(i) ) then
              hp(hlin(i)) = hp(hlin(i)) + hval(i) * p(hcol(i))
              if ( hlin(i) .ne. hcol(i) ) then
                  hp(hcol(i)) = hp(hcol(i)) + hval(i) * p(hlin(i))
              end if
          end if
      end do

C     ------------------------------------------------------------------
C     ADD ( dpdcj \nabla^2 cj(x) ) p
C     ------------------------------------------------------------------

      do j = 1,m

          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

C             ADD \nabla^2 cj(x) p
              do i = hsta(j),hsta(j) + hlen(j) - 1
                  if ( hlin(i) .ge. hcol(i) ) then
                      hp(hlin(i)) = 
     +                hp(hlin(i)) + dpdc(j) * hval(i) * p(hcol(i))
                      if ( hlin(i) .ne. hcol(i) ) then
                          hp(hcol(i)) = 
     +                    hp(hcol(i)) + dpdc(j) * hval(i) * p(hlin(i))
                      end if
                  end if
              end do

          end if

      end do

C     ------------------------------------------------------------------
C     ADD ( rhoj \nabla cj(x) \nabla cj(x)^t ) p
C     ------------------------------------------------------------------

      do j = 1,m

          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

C             COMPUTE THE INNER PRODUCT <a,p>
              atp = 0.0d0
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  atp = atp + jcval(i) * p(jcvar(i))
              end do

              atp = atp * rho(j)

C             ADD rhoj * atp * a
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  hp(jcvar(i)) = hp(jcvar(i)) + atp * jcval(i)
              end do

          end if

      end do

      go to 900

C     ==================================================================
C     END OF USER-PROVIDED SECOND DERIVATIVES
C     ==================================================================

C     ==================================================================
C     A SUBROUTINE TO COMPUTE THE PRODUCT OF A VECTOR TIMES THE HESSIAN
C     OF THE LAGRANGIAN WAS PROVIDED BY THE USER. THIS SUBROUTINE WILL
C     BE USED AND THEN THE FIRST-ORDER TERM WILL BE ADDED.
C     ==================================================================

 200  continue

C     ------------------------------------------------------------------
C     COMPUTE hp = ( \nabla^2 f(x) + dpdcj \nabla^2 cj(x) ) p
C     ------------------------------------------------------------------

      call evalhlp(n,x,m,dpdc,p,hp,goth,flag)
      if ( flag .ne. 0 ) then
          inform = - 96
          return
      end if

C     ------------------------------------------------------------------
C     ADD ( rhoj \nabla cj(x) \nabla cj(x)^t ) p
C     ------------------------------------------------------------------

      do j = 1,m

          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

C             COMPUTE THE INNER PRODUCT <a,p>
              atp = 0.0d0
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  atp = atp + jcval(i) * p(jcvar(i))
              end do

              atp = atp * rho(j)

C             ADD rhoj * atp * a
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  hp(jcvar(i)) = hp(jcvar(i)) + atp * jcval(i)
              end do

          end if

      end do

      go to 900

C     ==================================================================
C     END OF USER-PROVIDED HESSIAN-OF-THE-LAGRANGIAN X VECTOR PRODUCT
C     ==================================================================

C     ==================================================================
C     QUASI-NEWTON CORRECTION OF THE GAUSS-NEWTON APPROXIMATION OF THE 
C     HESSIAN MATRIX
C     ==================================================================

 300  continue

C     ------------------------------------------------------------------
C     COMPUTE THE QUASI-NEWTON CORRECTION OF THE GAUSS-NEWTON 
C     APPROXIMATION OF H IF THIS IS THE FIRST TIME THIS SUBBROUTINE IS 
C     BEING CALLED
C     ------------------------------------------------------------------

      if ( .not. goth ) then
          goth = .true.
          call comph(n,m,lambda,rho,equatn,linear,s,y,ssupn,seucn,yeucn,
     +    sts,sty,lspgmi,lspgma,hlspg,hds,hstds)
      end if

C     ------------------------------------------------------------------
C     COMPUTE THE STRUCTURED SPECTRAL CORRECTION S d, WHERE 
C     S = ( lamspg I ) AND lamspg = s^t (y - rho A^T A s) / (s^t s)
C     ------------------------------------------------------------------

C     hp = hlspg p

      do i = 1,n
          hp(i) = hlspg * p(i)
      end do

C     ------------------------------------------------------------------
C     ADD THE BFGS CORRECTION B p, WHERE
C     B = [ y y ^t / ( y^t s ) ] - [ D s ( D s )^t / ( s^t D s ) ], AND
C     D = hlspg I + rho A^T A
C     ------------------------------------------------------------------

C     hp = hp + B p = ( hlspg I + B ) p

      if ( samefa .and. sty .gt. 1.0d-08 * seucn * yeucn ) then

          pty = 0.0d0
          ptds = 0.0d0
          do i = 1,n
              pty = pty + p(i) * y(i)
              ptds = ptds + p(i) * hds(i)
          end do

          c1 = pty / sty
          c2 = ptds / hstds
          do i = 1,n
              hp(i) = hp(i) + c1 * y(i) - c2 * hds(i)
          end do

      end if

C     ------------------------------------------------------------------
C     ADD ( rhoj \nabla cj(x) \nabla cj(x)^t ) p
C     ------------------------------------------------------------------

C     hp = hp + B p = ( hlspg I + B + rho A^T A ) p

      do j = 1,m

          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

C             COMPUTE THE INNER PRODUCT <a,p>
              atp = 0.0d0
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  atp = atp + jcval(i) * p(jcvar(i))
              end do

              atp = atp * rho(j)

C             ADD rho * atp * a
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  hp(jcvar(i)) = hp(jcvar(i)) + atp * jcval(i)
              end do

          end if

      end do

      go to 900

C     ==================================================================
C     END OF GAUSS-NEWTON APPROXIMATION OF THE HESSIAN MATRIX
C     ==================================================================

C     ==================================================================
C     INCREMENTAL QUOTIENTS APPROXIMATION OF THE HESSIAN-VECTOR PRODUCT
C     TAKING CARE OF THE NON-DIFFERENTIABILITY
C     ==================================================================

 400  continue

C     COMPUTE INCREMENTAL QUOTIENTS STEP
      xsupn = 0.0d0
      psupn = 0.0d0
      do i = 1,n
          xsupn = max( xsupn, abs( x(i) ) )
          psupn = max( psupn, abs( p(i) ) )
      end do

      step = sqrt( macheps ) * max( xsupn / psupn, 1.0d0 )

C     SET THE POINT AT WHICH THE GRADIENT OF THE AUGMENTED LAGRANGIAN
C     WILL BE COMPUTED
      do i = 1,n
          xp(i) = x(i) + step * p(i)
      end do

      call setpoint(xp)

C     COMPUTE THE GRADIENT OF THE OBJECTIVE FUNCTION AT xp
      if ( gtype .eq. 0 ) then
          call evalg(n,xp,hp,flag)
          if ( flag .ne. 0 ) then
              inform = - 92
              return
          end if
      else
          call evalgdiff(n,xp,hp,macheps,inform)
          if ( inform .lt. 0 ) then
              return
          end if
      end if

C     COMPUTE \nabla L(x) = \nabla f(x) + \sum_{i=1}^m dPdc * dcdx

      do i = 1,n
          tmp(i) = 0.0d0
      end do

      do j = 1,m

          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

C             LINEAR CONSTRAINTS
              if ( linear(j) .and. exalin ) then

C                 COMPUTE THE GRADIENT OF THE i-TH CONSTRAINT
                  if ( gtype .eq. 0 ) then
                      call evaljac(n,xp,j,jcpvar,jcpval,jcpnnz,flag)
                      if ( flag .ne. 0 ) then
                          inform = - 93
                          return
                      end if
                  else
                      call evaljacdiff(n,xp,j,jcpvar,jcpval,jcpnnz,
     +                macheps,inform)
                      if ( inform .lt. 0 ) then
                          return
                      end if
                  end if

C                 COMPUTE THE INNER PRODUCT <dcdx,p>
                  atp = 0.0d0
                  do i = 1,jcpnnz
                      atp = atp + jcpval(i) * p(jcpvar(i))
                  end do

                  atp = atp * rho(j)

C                 ADD rho * atp * dcdx
                  if ( atp .ne. 0.0d0 ) then
                      do i = 1,jcpnnz
                          tmp(jcpvar(i)) = 
     +                    tmp(jcpvar(i)) + atp * jcpval(i)
                      end do
                  end if

C             NONLINEAR CONSTRAINTS
              else

C                 COMPUTE THE i-TH CONSTRAINT
                  call evalc(n,xp,j,cp,flag)
                  if ( flag .ne. 0 ) then
                      inform = - 91
                      return
                  end if

C                 COMPUTE dP/dc
                  dpdcp = lambda(j) + rho(j) * cp

                  if ( dpdcp .ne. 0.0d0 ) then

C                     COMPUTE THE GRADIENT OF THE i-TH CONSTRAINT
                      if ( gtype .eq. 0 ) then
                          call evaljac(n,xp,j,jcpvar,jcpval,jcpnnz,flag)
                          if ( flag .ne. 0 ) then
                              inform = - 93
                              return
                          end if
                      else
                          call evaljacdiff(n,xp,j,jcpvar,jcpval,jcpnnz,
     +                    macheps,inform)
                          if ( inform .lt. 0 ) then
                              return
                          end if
                      end if

C                     ADD dPdc * dcdx
                      do i = 1,jcpnnz
                          hp(jcpvar(i)) = 
     +                    hp(jcpvar(i)) + dpdcp * jcpval(i)
                      end do

                  end if

              end if

          end if

      end do

C     COMPUTE INCREMENTAL QUOTIENTS
      if ( exalin ) then
          do i = 1,n
              hp(i) = ( hp(i) - g(i) ) / step + tmp(i)
          end do 
      else
          do i = 1,n
              hp(i) = ( hp(i) - nal(i) ) / step
          end do 
      end if

      go to 900

C     ==================================================================
C     END OF INCREMENTAL QUOTIENTS APPROXIMATION OF THE HESSIAN-VECTOR 
C     PRODUCT TAKING CARE OF THE NON-DIFFERENTIABILITY
C     ==================================================================

C     ==================================================================
C     PURE INCREMENTAL QUOTIENTS APPROXIMATION OF THE HESSIAN-VECTOR 
C     PRODUCT (WITHOUT TAKING CARE OF THE NON-DIFFERENTIABILITY)
C     ==================================================================

 500  continue

C     COMPUTE INCREMENTAL QUOTIENTS STEP
      xsupn = 0.0d0
      psupn = 0.0d0
      do i = 1,n
          xsupn = max( xsupn, abs( x(i) ) )
          psupn = max( psupn, abs( p(i) ) )
      end do

      step = sqrt( macheps ) * max( xsupn / psupn, 1.0d0 )

C     SET THE POINT AT WHICH THE GRADIENT OF THE AUGMENTED LAGRANGIAN
C     WILL BE COMPUTED
      do i = 1,n
          xp(i) = x(i) + step * p(i)
      end do

      call setpoint(xp)

C     COMPUTE THE GRADIENT OF THE AUGMENTED LAGRANGIAN AT xp
      call evalnalu(n,xp,m,lambda,rho,equatn,linear,hp,gtype,macheps,
     +inform)
      if ( inform .lt. 0 ) then
          return
      end if

C     COMPUTE INCREMENTAL QUOTIENTS
      do i = 1,n
          hp(i) = ( hp(i) - nal(i) ) / step
      end do 

      go to 900

C     ==================================================================
C     END OF PURE INCREMENTAL QUOTIENTS APPROXIMATION OF THE HESSIAN-
C     VECTOR PRODUCT (WITHOUT TAKING CARE OF THE NON-DIFFERENTIABILITY)
C     ==================================================================

 900  continue

      end

C     *****************************************************************
C     *****************************************************************

      subroutine comph(n,m,lambda,rho,equatn,linear,s,y,ssupn,seucn,
     +yeucn,sts,sty,lspgmi,lspgma,hlspg,hds,hstds)

      implicit none

C     SCALAR ARGUMENTS
      integer m,n
      double precision lspgma,lspgmi,hlspg,hstds,seucn,ssupn,sts,sty,
     +        yeucn

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision hds(n),lambda(m),rho(m),s(n),y(n)

C     Consider the Hessian matrix of the Augmented Lagrangian function
C
C     M = \nabla^2 f(x) + \sum_i \rho_i \nabla^2 c(x) + 
C         \sum_i \rho_i \nabla c(x) \nabla c(x)^T,
C
C     where the indices i are such that constraint c_i(x) is an 
C     equality constraint or \lambda_i + \rho_i c_i(x) \geq 0.
C
C     This subroutine computes an approximation H of matrix M following
C     a very simple idea: "discard the second order terms and then 
C     correct the remaining matrix in order to satisfy a secant 
C     equation".
C
C     Hence, H takes the form
C
C         H = B + S + rho A^T A, 
C
C     where S is the spectral correction of (rho A^T A) and B is 
C     the BFGS correction of (S + diag(rho A^T A)). More specifically,
C
C     S = hlspg I,
C
C     where
C
C     hlspg = max(lspgmi, min(lspgma, s^T (y - rho A^T A s) / s^T s))
C
C     and
C
C     B = [ y y ^t / ( y^t s ) ] - [ B s ( B s )^t / ( s^t B s ) ].
C
C     Note that this subroutine does not compute matrix H explicitly,
C     but computes some quantities that will be used latter, by 
C     subroutine evalhd, to compute the product of H by a vector d.
C
C     The quantities computed by this subroutine are:
C
C     (a) hlspg = s^T (y - rho A^T A s) / (s^T s)
C
C     (b) hds = ( hlspg I + rho A^T A ) s, and
C
C     (c) hstds = <s,ds>.

C     PARAMETERS
      integer mmax,nmax,jcnnzmax
      parameter ( mmax      =  500000 )
      parameter ( nmax      =  500000 )
      parameter ( jcnnzmax  = 5000000 )

C     COMMON SCALARS
      logical constrc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     LOCAL SCALARS
      integer i,j
      double precision ats

C     COMMON BLOCKS
      common /graddata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc

C     ==================================================================
C     GAUSS-NEWTON CORRECTION OF THE HESSIAN MATRIX
C     ==================================================================

C     COMPUTE hds = rho A^t A s

      do i = 1,n
          hds(i) = 0.0d0
      end do

      do j = 1,m

          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then

C             COMPUTE THE INNER PRODUCT <a,s>
              ats = 0.0d0
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  ats = ats + jcval(i) * s(jcvar(i))
              end do

              ats = ats * rho(j)

C             ADD rho * ats * a
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  hds(jcvar(i)) = hds(jcvar(i)) + ats * jcval(i) 
              end do

          end if

      end do

      hstds = 0.0d0
      do i = 1,n
          hstds = hstds + s(i) * hds(i)
      end do

C     COMPUTE THE STRUCTURED SPECTRAL CORRECTION

C     ------------------------------------------------------------------
C     hlspg = s^t (y - rho A^T A s) / (s^t s)
C     ------------------------------------------------------------------

      if ( sty - hstds .le. 0.0d0 ) then
          hlspg = lspgmi
      else
          hlspg = max( lspgmi, min( (sty - hstds) / sts, lspgma ) )
      end if

      do i = 1,n
          hds(i) = hds(i) + hlspg * s(i)
      end do

      hstds = hstds + hlspg * sts

C     ==================================================================
C     END OF GAUSS-NEWTON APPROXIMATION OF THE HESSIAN MATRIX
C     ==================================================================

      end

C     *****************************************************************
C     *****************************************************************
      subroutine compp(n,m,lambda,rho,equatn,linear,s,y,ssupn,seucn,
     +yeucn,sts,sty,lspgmi,lspgma,samefa,pdiag,plspg,psmdy,psmdyty)

      implicit none

C     SCALAR ARGUMENTS
      logical samefa
      integer m,n
      double precision lspgma,lspgmi,plspg,psmdyty,seucn,ssupn,sts,sty,
     +        yeucn

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),rho(m),pdiag(n),psmdy(n),s(n),y(n)

C     Consider the preconditioner 
C
C         P = Q + E + diag(rho A^T A) 
C
C     for matrix 
C
C         H = B + S + rho A^T A, 
C
C     where E is the spectral correction of diag(rho A^T A) and Q is 
C     the BFGS correction of (E + diag(rho A^T A)), while S and B are
C     the spectral and BFGS corrections of matrix (rho A^T A), 
C     respectively. 
C
C     This subroutine computes:
C
C     (a) pdiag = diag(rho A^T A),
C
C     (b) plspg such that E = plspg I, 
C
C     (c) psmdy = s - D^-1 y, where D = E + diag(rho A^T A), and
C
C     (d) the inner product psmdty = <psmdy,y>.
C
C     These quantities will be used latter, in subroutine applyp, to
C     compute z = P^{-1} r.

C     PARAMETERS
      integer mmax,nmax,jcnnzmax
      parameter ( mmax      =  500000 )
      parameter ( nmax      =  500000 )
      parameter ( jcnnzmax  = 5000000 )

C     COMMON SCALARS
      logical constrc

C     COMMON ARRAYS
      integer jcvar(jcnnzmax),jcsta(mmax),jclen(mmax)
      double precision dpdc(mmax),g(nmax),jcval(jcnnzmax)

C     LOCAL SCALARS
      integer i,j
      double precision sttmp

C     COMMON BLOCKS
      common /graddata/ g,dpdc,jcval,jcvar,jcsta,jclen,constrc

C     ------------------------------------------------------------------
C     COMPUTE THE DIAGONAL OF (rho A^T A)
C     ------------------------------------------------------------------

      do i = 1,n
          pdiag(i) = 0.0d0
      end do

      do j = 1,m
          if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then
              do i = jcsta(j),jcsta(j) + jclen(j) - 1
                  pdiag(jcvar(i)) = 
     +            pdiag(jcvar(i)) + rho(j) * jcval(i) ** 2
              end do
          end if
      end do

      sttmp = 0.0d0
      do i = 1,n
          sttmp = sttmp + pdiag(i) * s(i) ** 2
      end do

C     ------------------------------------------------------------------
C     COMPUTE THE SPECTRAL CORRECTION E OF diag(rho A^T A)
C     ------------------------------------------------------------------

      if ( sty - sttmp .le. 0.0d0 ) then
          plspg = lspgmi
      else
          plspg = max( lspgmi, min( ( sty - sttmp ) / sts, lspgma ) )
      end if

C     ------------------------------------------------------------------
C     COMPUTE THE BFGS CORRECTION Q OF (E + diag(rho A^T A))
C
C     Q = [ (s - D^-1 y) s^t + s (s - D^-1 y)^t ] / s^t y -
C         [ <s - D^-1 y, y> s s^t ] / (s^t y)^2,
C
C     WHERE D = (E + diag(rho A^T A))
C     ------------------------------------------------------------------
 
      if ( samefa .and. sty .gt. 1.0d-08 * seucn * yeucn ) then

          psmdyty = 0.0d0
          do i = 1,n
              psmdy(i) = s(i) - y(i) / ( plspg + pdiag(i) )
              psmdyty = psmdyty + psmdy(i) * y(i)
          end do

      end if

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evalgdiff(n,x,g,macheps,inform)

      implicit none

C     This subroutine approximates, by central finite differences, the 
C     gradient of the objective function.

C     SCALAR ARGUMENTS
      integer inform,n
      double precision macheps

C     ARRAY ARGUMENTS
      double precision g(n),x(n)

C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     N integer: Number of variables.
C     ---------
C
C     X double precision x(n): Current point.
C     -----------------------
C
C     On Return:
C     ==========
C
C     G double precision g(n): Gradient of the objective function.
C     -----------------------
C
C     INFORM integer
C     --------------
C
C     0 means that the evaluation was successfuly done.
C
C     Any other value means that some error occured during the 
C     evaluation.

C     LOCAL SCALARS
      integer flag,j
      double precision fminus,fplus,step,tmp

      do j = 1,n
          tmp  = x(j)

          step = macheps ** (1.0d0/3.0d0) * max( abs( tmp ), 1.0d0 )

          x(j) = tmp + step
          call evalf(n,x,fplus,flag)
          if ( flag .ne. 0 ) then
              inform = - 90
              return
          end if

          x(j) = tmp - step
          call evalf(n,x,fminus,flag)
          if ( flag .ne. 0 ) then
              inform = - 90
              return
          end if

          g(j) = ( fplus - fminus ) / ( 2.0d0 * step )
          x(j) = tmp
      end do

      end

C     ******************************************************************
C     ******************************************************************

      subroutine evaljacdiff(n,x,ind,indjac,valjac,nnzjac,macheps,
     +inform)

      implicit none

C     This subroutine approximates, by central finite differences, the 
C     gradient of the ind-th constraint.

C     SCALAR ARGUMENTS
      integer ind,inform,n,nnzjac
      double precision macheps

C     ARRAY ARGUMENTS
      integer indjac(n)
      double precision valjac(n),x(n)

C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     N integer: Number of variables.
C     ---------
C
C     X double precision x(n): Current point.
C     -----------------------
C
C     IND integer: index of the constraint.
C     -----------
C
C     On Return:
C     ==========
C
C     NNZJAC integer: number of non-null elements in the gradient vector.
C     --------------
C
C     INDJAC integer indjac(nnzjac): see below.
C     -----------------------------
C
C     VALJAC double precision valjac(nnzjac):
C     --------------------------------------
C
C     The value of the indjac(i)-th non-null element of the gradient 
C     vector is valjac(i), for i = 1, ..., nnzjac.
C
C     INFORM integer
C     --------------
C
C     0 means that the evaluation was successfuly done.
C
C     Any other value means that some error occured during the 
C     evaluation.

C     LOCAL SCALARS
      integer flag,j
      double precision cminus,cplus,step,tmp

      nnzjac = 0

      do j = 1,n
          tmp  = x(j)

          step = macheps ** (1.0d0/3.0d0) * max( abs( tmp ), 1.0d0 )

          x(j) = tmp + step
          call evalc(n,x,ind,cplus,flag)
          if ( flag .ne. 0 ) then
              inform = - 91
              return
          end if

          x(j) = tmp - step
          call evalc(n,x,ind,cminus,flag)
          if ( flag .ne. 0 ) then
              inform = - 91
              return
          end if

          indjac(nnzjac + 1) = j
          valjac(nnzjac + 1) = ( cplus - cminus ) / ( 2.0d0 * step )

          if ( abs( valjac(nnzjac + 1) ) .gt. 0.0d0 ) then
              nnzjac = nnzjac + 1
          end if

          x(j) = tmp
      end do

      end

c#endif

C     *****************************************************************
C     *****************************************************************

      subroutine evalp(y,rho,lambda,equatn,p)

      implicit none

C     SCALAR ARGUMENTS
      logical equatn
      double precision lambda,p,rho,y

      if ( equatn ) then
          p  = y * ( lambda + 0.5d0 * rho * y )
      else
          if ( lambda + rho * y .ge. 0.0d0 ) then
              p = y * ( lambda + 0.5d0 * rho * y )
          else
              p = - 0.5d0 * lambda ** 2 / rho
          end if
      end if

      return

      end

C     *****************************************************************
C     *****************************************************************

      subroutine evaldpdy(y,rho,lambda,equatn,dpdy)

      implicit none

C     SCALAR ARGUMENTS
      logical equatn
      double precision y,rho,lambda,dpdy

      if ( equatn ) then
          dpdy = lambda + rho * y
      else
          dpdy = max( 0.0d0, lambda + rho * y )
      end if

      return

      end

C     *****************************************************************
C     *****************************************************************

      subroutine applyp(n,r,m,lambda,rho,equatn,linear,s,y,ssupn,seucn,
     +yeucn,sts,sty,lspgmi,lspgma,samefa,gotp,pdiag,plspg,psmdy,psmdyty,
     +z)

      implicit none

C     SCALAR ARGUMENTS
      logical gotp,samefa
      integer m,n
      double precision lspgma,lspgmi,plspg,psmdyty,seucn,ssupn,sts,sty,
     +        yeucn

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision lambda(m),pdiag(n),psmdy(n),r(n),rho(m),s(n),
     +        y(n),z(n)

C     Consider the preconditioner 
C
C         P = Q + E + diag(rho A^T A) 
C
C     for matrix 
C
C         H = B + S + rho A^T A, 
C
C     where E is the spectral correction of diag(rho A^T A) and Q is 
C     the BFGS correction of (E + diag(rho A^T A)), while S and B are
C     the spectral and BFGS corrections of matrix (rho A^T A), 
C     respectively. 
C
C     Given the quantities computed in advance in the subroutine compp,
C     this subroutine computes z = P^{-1} r.

C     LOCAL SCALARS
      integer i
      double precision c1,c2,psmdytr,str

C     ------------------------------------------------------------------
C     COMPUTE P IF THIS IS THE FIRST TIME THIS SUBBROUTINE IS BEING 
C     CALLED
C     ------------------------------------------------------------------

      if ( .not. gotp ) then
          gotp = .true.
          call compp(n,m,lambda,rho,equatn,linear,s,y,ssupn,seucn,
     +    yeucn,sts,sty,lspgmi,lspgma,samefa,pdiag,plspg,psmdy,psmdyty)
      end if

C     ------------------------------------------------------------------
C     COMPUTE z = (E + diag(rho A^T A))^{-1} r
C     ------------------------------------------------------------------

      do i = 1,n
          z(i) = r(i) / ( plspg + pdiag(i) )
      end do

C     ------------------------------------------------------------------
C     ADD Q^{-1} r, WHERE
C
C     Q^{-1} = [ (s - D^-1 y) s^t + s (s - D^-1 y)^t ] / s^t y -
C              [ <s - D^-1 y, y> s s^t ] / (s^t y)^2
C
C     AND
C
C     D = (E + diag(rho A^T A))
C     ------------------------------------------------------------------

      if ( samefa .and. sty .gt. 1.0d-08 * seucn * yeucn ) then

          str = 0.0d0
          psmdytr = 0.0d0
          do i = 1,n
              str = str + s(i) * r(i)
              psmdytr = psmdytr + psmdy(i) * r(i)
          end do

          c1 = str / sty
          c2 = psmdytr / sty - psmdyty * str / sty ** 2

          do i = 1,n
              z(i) = z(i) + c1 * psmdy(i) + c2 * s(i) 
          end do

      end if

      end

C     =================================================================
C     Module: Bound-constraints solver GENCAN
C     =================================================================

C     Last update of any of the component of this module: 
C
C     May 11, 2007.

      subroutine easygencan(n,x,l,u,m,lambda,equatn,linear,rho,gtype,
     +hptype,intype,precond,epsgpsn,maxit,maxfc,iprint,ncomp,macheps,
     +bignum,f,g,gpsupn,iter,fcnt,gcnt,cgcnt,inform,wi1,wd1,wd2,wd3,wd4,
     +wd5,wd6,wd7,wd8,wd9,wd10,wd11,wd12,wd13,wd14)

      implicit none

C     SCALAR ARGUMENTS
      character * 6 precond
      integer cgcnt,fcnt,gcnt,gtype,hptype,intype,m,maxfc,maxit,n,ncomp,
     +        inform,iprint,iter
      double precision bignum,epsgpsn,f,gpsupn,macheps

C     ARRAY ARGUMENTS
      integer wi1(n)
      logical equatn(m),linear(m)
      double precision g(n),l(n),lambda(m),rho(m),u(n),wd1(n),wd2(n),
     +        wd3(n),wd4(n),wd5(n),wd6(n),wd7(n),wd8(n),wd9(n),wd10(n),
     +        wd11(n),wd12(n),wd13(n),wd14(n),x(n)

C     This subroutine aims to simplify the use of GENCAN. For this 
C     purpose it gives values to most of the GENCAN arguments and 
C     leaves to the user those arguments which he/she may would like to 
C     set by him/herself.
C
C     The arguments of EASYGENCAN are the input and output arguments of 
C     GENCAN that are supposed to be useful for a common user. The input 
C     arguments are mostly related to basic problem information, like 
C     dimension and bounds, and the initial point. There are also input 
C     arguments related to simple stopping criteria (like norm of the 
C     projected gradient, and maximum number of iterations and 
C     functional evaluations). There are also two input arguments 
C     related to control the amount of information written into the 
C     screen. The output arguments are related to information of the 
C     solution and some few performance measurements. Basically, on 
C     return, EASYGENCAN gives to the user the solution, the objective 
C     functional value and its gradient at the solution, Euclidian and 
C     sup-norm of the projected gradient at the solution, the number of 
C     iterations, functional and gradient evaluations, and Conjugate 
C     Gradient iterations used to reach the solution, and, finally, a 
C     flag that indicates the stopping criterion that was satisfied.
C
C     All the other arguments of GENCAN are setted with its default 
C     values by EASYGENCAN. EASYGENCAN divides the arguments of GENCAN 
C     in two sets. Those that are related to the behaviour of GENCAN are 
C     declared as Fortran parameters (constants). The other arguments of 
C     GENCAN, most of them related to alternative stopping criteria, and 
C     that may depend of, for example, maxit, are declared as local 
C     variables of EASYGENCAN.
C
C     GENCAN arguments that are defined as Fortran parameters in this 
C     subroutine are GENCAN arguments that should not be modified by a 
C     common user. They are arguments that modify the behaviour of 
C     GENCAN and whos values were selected because they are classical 
C     values in some cases or because some numerical experiments seemed 
C     to indicate that they are the best choices.
C
C     GENCAN arguments that are declared as local variables in this 
C     subroutine are GENCAN arguments that may be modified if, with 
C     their suggested values, GENCAN does not give the desired result. 
C     Most of them are related to Conjugate Gradients or to disabled 
C     stopping criteria that may be useful in bad-scaled problems or 
C     problems with not trustable derivatives.
C
C     Finally, this subroutine declares as local variables some 
C     arguments of GENCAN which in fact are output arguments. Most of 
C     them are related to quantities that can be used for statistics 
C     related to the GENCAN performance, like number Spectral Projected 
C     Gradient iterations, Truncated Newton iterations, Conjugate 
C     Gradient iterations, etc. As we assume that this values are not 
C     useful for the common user, this subroutine throw all of them 
C     away.
C
C     We describe below the meaning of the arguments of the EASYGENCAN
C     subroutine. More detailed descriptions as well as the descriptions 
C     of all the other GENCAN arguments that are not arguments of 
C     EASYGENCAN are also described at the begining of the GENCAN 
C     subroutine.
C     
C     Parameters of the subroutine:
C
C     On entry:
C
C     n        integer 
C              number of variables
C
C     x        double precision x(n)
C              initial estimation of the solution
C
C     l        double precision l(n)
C              lower bounds on the variables
C
C     u        double precision u(n)
C              upper bounds on the variables
C
C     m        integer
C     lambda   double precision lambda(m)
C     rho      double precision rho(m)
C     equatn   logical equatn(m)
C     linear   logical linear(m)
C              These five parameters are not used nor modified by 
C              GENCAN and they are passed as arguments to the user-
C              defined subroutines evalal and evalnal to compute the 
C              objective function and its gradient, respectively. 
C              Clearly, in an Augmented Lagrangian context, if GENCAN is 
C              being used to solve the bound-constrained subproblems, m 
C              would be the number of constraints, lambda the Lagrange 
C              multipliers approximation and rho the penalty parameters.
C              equatn is logical array that, for each constraint, 
C              indicates whether the constraint is an equality constraint
C              (.true.) or an inequality constraint (.false.). Finally,
C              linear is logical array that, for each constraint, 
C              indicates whether the constraint is a linear constraint
C              (.true.) or a nonlinear constraint (.false.)
C
C     epsgpsn  double precision
C              GENCAN stops declaring convergence if it finds a point 
C              whos projected gradient sup-norm is smaller than or equal 
C              to epsgpsn
C
C     maxit    integer
C              GENCAN stops declaring ''maximum number of iteration 
C              achieved'' if the number of iterations exceeds maxit
C
C     maxfc    integer
C              the same as before but with the number of functional 
C              evaluations
C
C     iprint   integer
C              indicates the degree of details of the output generated 
C              by GENCAN. Setting iprint to a value smaller than 2 will 
C              make GENCAN to generate no output at all. An iprint value 
C              greater than or equal to 2 will generate information of 
C              every GENCAN iteration. An iprint value greater than or 
C              equal to 3 will also show information of the Conjugate 
C              Gradient iterations (used to compute the Truncated Newton 
C              direction) and also information related to the line 
C              search procedures in the Spectral Projected Gradient 
C              direction and the Truncated Newton direction.
C
C     ncomp    integer,         
C              every time a vector is printed, just its first ncomp
C              component will be displayed.
C
C     gtype    integer,
C              type of derivatives calculation according to the 
C              following convention:
C
C              = 0, true derivatives. In this case, subroutines evalg
C                   and evaljac must be modified by the user to compute
C                   the derivatives of the objective function and the
C                   constraints.
C              = 1, finite differences approximation provided by 
C                   ALGENCAN. In this case, subroutines evalg and 
C                   evaljac may have an empty body but must be present.
C                   It is also recommended that those empty-body
C                   subroutines set flag = - 1. Last but not least,
C                   the option gtype = 1 is not cheap neither safe.
C
C     hptype   integer
C              The way in which the product of the Hessian matrix (of 
C              the Augmented Lagrangian) by a vector will be done 
C              depends on the value of the parameter hptype in the 
C              following way:
C
C              Read this explanation if your are solving a problem
C              with constraints (other than the bound constraints).
C              Otherwise, see the explanation below.
C
C              = 9, means that an incremental quotients approximation 
C                   without any extra consideration will be used. This 
C                   option requires the evaluation of an extra gradient 
C                   at each Conjugate Gradient iteration. If gtype = 0 
C                   then this gradient evaluation will be done using the 
C                   user supplied subroutines evalg and evaljac (evalc 
C                   will also be used). On the other hand, if gtype = 1, 
C                   the gradient calculation will be done using just 
C                   calls to the user provided subroutines evalf and 
C                   evalc. nind calls will be done, where nind is the 
C                   dimension of the current face of the active-set 
C                   method.
C
C              If you did not code subroutines evalg and evaljac, to 
C              compute the gradient of the objective function and the 
C              Jacobian of the constraints, respectively, then your 
C              options finished here.
C
C              = 0, means that subroutines to compute the Hessian of the
C                   objective function (evalh) and the Hessians of the 
C                   constraints (evalhc) were provided by the user. So, 
C                   the product of the Hessian of the Augmented 
C                   Lagrangian times a vector will be computed using the
C                   Hessians provided by these subroutines and then 
C                   adding the first order term (for the first order term 
C                   the user-supplied subroutine to compute the Jacobian 
C                   of the constraints (evaljac) is also used).
C
C              = 1, means that, instead of providing individual 
C                   subroutines to compute the Hessians of the objective
C                   function and the constraints, the user provided a
C                   subroutine to compute the product of an arbitrary 
C                   vector times the Hessian of the Lagrangian. This
C                   option was basically designed for the AMPL and CUTEr
C                   interfaces.
C
C              = 2, means that incremental quotients will be used. The 
C                   difference between hptype = 9 and hptype = 2 is that, 
C                   in the latter case, the non-differentiability of the 
C                   Hessian of the Augmented Lagrangian will be taken 
C                   into account. In particular, when computing the 
C                   gradient of the Augmented Lagrangian at (x + step p), 
C                   the constraints that will be considered will be the 
C                   same constraints that contributed to the computation 
C                   of the gradient of the Augmented Lagrangian at x.
C
C                   If GENCAN is been used to solve a bound-constrained 
C                   problem which is not the subproblem of an Augmented 
C                   Lagrangian method, but an isolated bound-constrained 
C                   problem, then there is no difference between this 
C                   option and hptype = 9.
C
C                   This option also requires the evaluation of an extra 
C                   gradient at each Conjugate Gradient iteration. 
C
C              = 3, is similar to hptype = 2. The difference is that
C                   the contribution of the linear constraints in the 
C                   Hessian matrix will be computed explicitly. If the 
C                   problem has not linear constraints then this option 
C                   is identical to hptype = 2. Moreover, if the problem
C                   has no constraints then this option is equal to
C                   hptype = 9.
C
C                   This option also requires the evaluation of an extra 
C                   gradient at each Conjugate Gradient iteration. 
C
C              = 4, means that the Hessian matrix will be approximated
C                   and then the product of the Hessian approximation
C                   by the vector will be computed exactly. In 
C                   particular, the Hessian matrix will be approximated
C                   doing a BFGS correction to the Gauss-Newton
C                   approximation of the Hessian. Before the BFGS
C                   correction, a structured spectral correction is done
C                   to force the Gauss-Newton approximation to be
C                   positive definite.
C
C                   If the problem has not constraints then the 
C                   approximation reduces to a BFGS approximation of
C                   the Hessian (without memory) and using the spectral
C                   approximation (instead of the identity) as initial
C                   approximation.
C
C                   Numerical experiments suggested that this option is
C                   convenient just for constrained problems. This 
C                   motivated the introduction of the next option.
C
C                   This option does NOT require an extra gradient
C                   evaluation per iteration and, in this sense, each
C                   CG iteration is computationally cheaper than a CG
C                   iteration of the previous choices. However, the 
C                   approximation of the Hessian matrix requires some 
C                   information (mainly the Jacobian of the constraints) 
C                   that must be saved during the gradient evaluation. 
C                   To save this information requires an amount of memory 
C                   proportional to the number of non-null elements of 
C                   the Jacobian matrix.
C 
C                   Quadratic subproblems are convex with this choice.
C
C              = 5, is an adaptive strategy that choose, at every 
C                   iteration, between 2 or 4. When the gradient of
C                   the Augmented Lagrangian is computed, it is verified
C                   if at least a constraint contributes to the
C                   calculation. If this is the case, 4 is used. 
C                   Otherwise, 2 is used.
C
C                   For problems with equality constraints (that always
C                   contributes to the Augmented Lagrangian function) this
C                   option is identical to 4.
C
C                   For problems without constraints this option is
C                   identical to 2.
C
C              = 6, is identical to 5 but the choice is made between 3
C                   and 4 instead of between 2 and 4. 
C
C                   For problems with equality constraints (that always
C                   contributes to the Augmented Lagrangian function) this
C                   option is identical to 4.
C
C                   For problems without constraints this option is
C                   identical to 3.
C
C              Read this explanation if your are solving an unconstrained
C              or bound-constrained problem:
C
C              = 0, means that a subroutine provided by the user will be
C                   used to compute the Hessian-vector product.
C
C              = 0, means that the subroutine to compute the Hessian of 
C                   the objective function (evalh) was provided by the 
C                   user. So, the product of the Hessian times a vector 
C                   will be computed using the Hessian provided by this 
C                   subroutine.
C
C              = 9, means that an incremental quotients approximation 
C                   will be used. This option requires the evaluation of 
C                   an extra gradient at each Conjugate Gradient 
C                   iteration. If gtype = 0 then this gradient evaluation 
C                   will be done using the user supplied subroutine 
C                   evalg. On the other hand, if gtype = 1, the gradient 
C                   calculation will be done using just calls to the user 
C                   provided subroutine evalf. nind calls will be done, 
C                   where nind is the dimension of the current face of 
C                   the active-set method.
C
C              In the unconstrained or bound-constrained context, options 
C              hptype = 2 and hptype = 3 are both identical to hptype = 9.
C
C              If you did not code subroutine evalg to compute the 
C              gradient of the objective function then your options 
C              finished here.
C
C              = 4, means that the Hessian matrix will be approximated
C                   and then the product of the Hessian approximation
C                   by the vector will be computed exactly. The 
C                   approximation is a BFGS approximation of the 
C                   Hessian (without memory) and using the spectral
C                   approximation (instead of the identity) as initial
C                   approximation.
C
C                   Numerical experiments suggested that this option is
C                   not convenient for unconstrained or just bound-
C                   constrained problems. (Note that this option was
C                   developed to be used in the Augmented Lagrangian
C                   framework.)
C
C                   This option does NOT require an extra gradient
C                   evaluation per iteration and, in this sense, each
C                   CG iteration is computationally cheaper than a CG
C                   iteration of the previous choices.
C
C                   Quadratic subproblems are convex with this choice.
C
C              In the unconstrained or bound-constrained context, 
C              options hptype = 5 and hptype = 6 are both identical to 
C              hptype = 9.
C
C     precond  character * 6
C              indicates the type of preconditioning that will be used
C              for Conjugates Gradients.
C
C              'NONE'   means no preconditioner at all,
C
C              'QNAGNC' means Quasi-Newton Correction of the Gauss-
C                       Newton approximation of the Hessian. The exact
C                       form is this preconditioner is described in:
C 
C                       E. G. Birgin and J. M. Martínez, "Structured 
C                       minimal-memory inexact quasi-Newton method and 
C                       secant preconditioners for Augmented Lagrangian 
C                       Optimization", submitted, 2005.
C
C     macheps  double precision
C              macheps is the smallest positive number such that
C              1 + macheps is not equal to 1
C
C     bignum   double precision
C              a big number like 1.0d+99
C
C     wi1      integer wi1(n)
C              n-dimensional working vector of integers
C
C     wd1, ..., wd14 double precision wd1(n), ..., wd14(n)
C              n-dimensional working vectors of doubles
C              
C     On return:
C
C     x        double precision x(n)
C              estimation of the solution
C
C     f        double precision
C              objective function value at the solution
C
C     g        double precision g(n)
C              gradient of the objective function at the solution
C
C     gpsupn   double precision
C              sup-norm of the continuous projected gradient
C
C     iter     integer
C              number of iterations used to reach the solution
C
C     fcnt     integer
C              number of functional evaluations
C
C     gcnt     integer
C              number of gradient evaluations
C
C     cgcnt    integer
C              number of Conjugate Gradient iterations
C
C     inform   integer
C              termination criteria. inform equal to 1 means that 
C              GENCAN converged with the sup-norm of the continuous 
C              projected gradient stopping criterion (inform equal to 0
C              means the same but with the Euclidian norm). Other 
C              positive values means that GENCAN stopped by a may be not  
C              successful stopping criteria. A negative value means that 
C              there was an error in the user-defined subroutines that 
C              computes the objective function (subroutine evalal), the 
C              gradient (subroutine evalnal), or the Hessian-vector
C              product (subroutine evalhd). See the GENCAN description 
C              for more details.

C     HERE STARTS THE DESCRIPTION OF SOME GENCAN ARGUMENTS THAT ARE 
C     BEING SETTED INSIDE EASYGENCAN. THE FIRST SET OF ARGUMENTS ARE 
C     THOSE ARGUMENTS THAT WE WILL CALL ''CONSTANTS'' AND THAT, AS THEIR 
C     VALUES ALTER THE BEHAVIOUR OF GENCAN, SHOULD NOT BE MODIFIED BY A 
C     COMMON USER.

C     CONSTANTS FOR CLASSICAL LINE-SEARCH CONDITIONS

C     beta is the constant for the ''beta condition''. We use this 
C     condition to test whether is promising to extrapolate or not.

C     gamma is the constant for the sufficient decrease ''Armijo 
C     condition''.

C     theta is the constant for the ''angle condition''.

C     sigma1 and sigma2 are the constants for the safeguarding quadratic 
C     interpolations. We use them in a rather unusual way. Instead of 
C     discarding a new step anew if it does not belong to the interval 
C     [ sigma1 * aprev, sigma2 * aprev ], we discard it if it does not 
C     belong to the interval [ sigma1, sigma2 * aprev ]. In such a case 
C     we take something similar to ''anew = aprev / 2''.

      double precision beta,gamma,theta,sigma1,sigma2
      parameter ( beta   =   0.5d0 )
      parameter ( gamma  = 1.0d-04 )
      parameter ( theta  = 1.0d-06 )
      parameter ( sigma1 =   0.1d0 )
      parameter ( sigma2 =   0.9d0 )

C     CONSTANTS FOR SPECIFIC PROCEDURES (NOT SO CLASSICAL)

C     In line searches, when interpolating, the step may become so 
C     small that we should declare a line search failure indicating that 
C     direction may not be a descent direction. This decision is never 
C     take before doing at least mininterp interpolations.

C     In line searches, the beta condition (see above) may recommend to
C     extrapolate. We never do more than maxextrap extrapolations.

C     In the line searches, when we need to interpolate and the result 
C     of the quadratic interpolation is rejected, the new step is 
C     computed as anew = aprev / etaint. When the beta condition 
C     recommends to extrapolate, we compute anew = aprev * etaext.

C     When computing the Newton direction by Conjugate Gradients we 
C     never go further an artificial ''trust region''. This ''trust 
C     radius'' is never smaller than delmin.

C     In active set strategies, constants eta is used to decide whether 
C     the current face should be abandoned or not. In particular, the 
C     current face is abandoned when the norm of the internal to face 
C     component of the continuous projected gradient is smaller than 
C     ( 1 - eta ) times the norm of the continuous projected gradient. 
C     In this way, values of eta near 1 makes the method to work hard 
C     inside the faces and values of eta near 0 makes the method to 
C     abandon the faces very quickly.

C     We always use as a first step in a line search procedure along a
C     first order direction the spectral steplength. This steplength 
C     must belong to the interval [lspgmi,lspgma].

      integer maxextrap,mininterp
      parameter ( maxextrap = 100 )
      parameter ( mininterp =   4 )

      double precision etaint,etaext,delmin,eta,lspgma,lspgmi
      parameter ( etaint  =   2.0d0 )
      parameter ( etaext  =   2.0d0 )
      parameter ( delmin  = 1.0d+04 )
      parameter ( eta     =   0.9d0 )
      parameter ( lspgma  = 1.0d+10 )
      parameter ( lspgmi  = 1.0d-10 )

C     HERE STARTS THE DESCRIPTION OF THE OTHER ARGUMENTS OF GENCAN BEING 
C     SETTED BY EASYGENCAN. THESE ARGUMENTS MAY BE MODIFIED BY A COMMON 
C     USER IF, WITH THEIR SUGGESTED VALUES, GENCAN DOES NOT GIVE THE 
C     EXPECTED RESULT.

C     GENCAN INPUT ARGUMENTS THAT WILL BE SETTED BELOW

      integer cgscre,maxitnfp,maxitngp,maxitnqmp,trtype
      double precision cgepsf,cgepsi,cggpnf,cgmia,cgmib,delta0,epsgpen,
     +        epsnfp,epsnqmp,fmin

C     GENCAN OUTPUT ARGUMENTS THAT WILL BE DISCARDED

      integer spgfcnt,spgiter,tnexbcnt,tnexgcnt,tnexbfe,tnexgfe,tnfcnt,
     +        tnintcnt,tnintfe,tniter,tnstpcnt
      double precision gpeucn2

C     ARGUMENTS RELATED TO STOPPING CRITERIA

C     Besides the stopping criterion related to the sup-norm of the 
C     continuous projected gradient, there is another stopping criterion 
C     related to its Euclidian norm. So, GENCAN stops the process if it 
C     finds a point at which the Euclidian norm of the continuous 
C     projected gradient is smaller than epsgpen.

      epsgpen   =    0.0d0

C     Sometimes, is the problem is bad scaled, to request a small 
C     gradient norm at the solution may be inadequate. For this reason, 
C     if this norm is not decreasing during maxitngp (MAXimum of 
C     ITerations with No Gradient Progress) consecutive iterations then 
C     we stop the method with a warning. 
     
      maxitngp  =    maxit

C     maxitnfp means MAXimum of allowed number of iterations with No 
C     Progress in the objective functional value. ''Progress'' from one 
C     iteration to the next one refers to ( fnew - fprev ). Since the 
C     begining of the algorithm we save the ''best progress'' and 
C     consider that there was no progress in an iteration if the 
C     progress of this iterations was smaller than epsnfp times the best 
C     progress. Finally, the algorithm stops if there was no progress 
C     during maxitnfp consecutive iterations.

      maxitnfp  =    maxit
      epsnfp    =    0.0d0

C     There is a stopping criterion that stops the method if a point 
C     with a functional value smaller than fmin is found. The idea 
C     behind this stopping criterion is to stop the method if the 
C     objective function is not bounded from bellow.

      fmin      = - bignum

C     ARGUMENTS RELATED TO CONJUGATE GRADIENTS

C     When computing the Truncated Newton direction by Conjugate 
C     Gradients there is something similar to a ''trust-region radius''. 
C     This trust radius is updated from iteration to iteration depending 
C     on the agreement of the objective function and its quadratic 
C     model. But an initial value for the trust radius is required. If 
C     the user has a good guess for this initial value then it should be 
C     passed to GENCAN using the delta0 arguments. On the other hand, if 
C     delta0 is set to -1, a default value depending on the norm of the 
C     current point will be used.

      delta0    =  - 1.0d0

C     The ''trust-region'' can be like a ball (using Euclidian norm) or 
C     like a box (using sup-norm). This choice can be made using trtype 
C     (TRust region TYPE) argument. trtype equal to 0 means Euclidian 
C     norm and trtype equal to 1 means sup-norm.

      trtype    =        0

C     When the method is far from the solution, it may be not useful to 
C     do a very large effort in computing the Truncated Newton direction 
C     precisely. To avoid it, a fixed maximum number of iterations for 
C     Conjugate Gradients can be given to GENCAN. If the user would like 
C     to choose this maximum number of iterations for Conjugate 
C     Gradient then it should use cgmia and cgmib arguments. On the other 
C     hand, if he/she prefers to leave this task to GENCAN then he/she 
C     should set cgmia = -1.0d0 and cgmib = -1.0d0.
 
      cgmia     =  - 1.0d0
      cgmib     =  - 1.0d0

C     If the task of deciding the accuracy for computing the Truncated 
C     Newton direction is leaved to GENCAN then a default strategy based 
C     on increasing accuracies will be used. The proximity to the 
C     solution is estimated observing the norm of the projected gradient 
C     at the current point and locating it between that norm at the 
C     initial point and the expected value of that norm at the solution. 
C     Then the accuracy for the Truncated Newton direction of the 
C     current iteration will be computed taking a precision located in 
C     the same relative position with respect to two given values for 
C     the accuracies for the first and the last Truncated Newton 
C     direction calculations. These two accuracies (cgepsi and cgepsf, 
C     respectively) must be given by the user. Moreover, the expected 
C     value of the projected gradient norm at the solution (cggpnf) must
C     also be given by the user who must indicate setting argument 
C     cgscre to 1 or 2 if that norm is the Euclidian or the sup-norm.
      
      cggpnf    =  max( 1.0d-04, max( epsgpen, epsgpsn ) ) 
      cgscre    =        2
      cgepsi    =  1.0d-01
      cgepsf    =  1.0d-08

C     The next two arguments are used for an alternative stopping 
C     criterion for Conjugate Gradients. Conjugate Gradients method is 
C     stopped if the quadratic model makes no progress during maxitnqmp 
C     (MAXimum of ITerations with No Quadratic Model Progress) 
C     consecutive iterations. In this context, ''no progress'' means 
C     that the progress is smaller than epsnqmp (EPSilon to measure the 
C     No Quadratic Model Progress) times the best progress obtained 
C     during the previous iterations.

      epsnqmp   =  1.0d-08
      maxitnqmp =        5

C     FINALLY, CALL GENCAN

      call gencan(n,x,l,u,m,lambda,equatn,linear,rho,gtype,hptype,
     +intype,precond,epsgpen,epsgpsn,maxitnfp,epsnfp,maxitngp,fmin,
     +maxit,maxfc,delta0,cgmia,cgmib,cgscre,cggpnf,cgepsi,cgepsf,
     +epsnqmp,maxitnqmp,etaint,etaext,mininterp,maxextrap,trtype,iprint,
     +ncomp,macheps,bignum,f,g,gpeucn2,gpsupn,iter,fcnt,gcnt,cgcnt,
     +spgiter,spgfcnt,tniter,tnfcnt,tnstpcnt,tnintcnt,tnexgcnt,tnexbcnt,
     +tnintfe,tnexgfe,tnexbfe,inform,wd1,wd2,wd3,wd4,wd5,wi1,wd6,wd7,
     +wd8,wd9,wd10,wd11,wd12,wd13,wd14,eta,delmin,lspgma,lspgmi,theta,
     +gamma,beta,sigma1,sigma2)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine gencan(n,x,l,u,m,lambda,equatn,linear,rho,gtype,hptype,
     +intype,precond,epsgpen,epsgpsn,maxitnfp,epsnfp,maxitngp,fmin,
     +maxit,maxfc,udelta0,ucgmia,ucgmib,cgscre,cggpnf,cgepsi,cgepsf,
     +epsnqmp,maxitnqmp,etaint,etaext,mininterp,maxextrap,trtype,iprint,
     +ncomp,macheps,bignum,f,g,gpeucn2,gpsupn,iter,fcnt,gcnt,cgcnt,
     +spgiter,spgfcnt,tniter,tnfcnt,tnstpcnt,tnintcnt,tnexgcnt,tnexbcnt,
     +tnintfe,tnexgfe,tnexbfe,inform,s,y,d,xprev,gprev,ind,wd1,wd2,wd3,
     +wd4,wd5,wd6,wd7,wd8,wd9,eta,delmin,lspgma,lspgmi,theta,gamma,beta,
     +sigma1,sigma2)

      implicit none

C     SCALAR ARGUMENTS
      character * 6 precond
      integer cgcnt,cgscre,fcnt,gcnt,gtype,hptype,inform,intype,iprint,
     +        iter,m,maxextrap,maxfc,maxit,maxitnfp,maxitngp,maxitnqmp,
     +        mininterp,n,ncomp,spgfcnt,spgiter,tnexbcnt,tnexbfe,
     +        tnexgcnt,tnexgfe,tnfcnt,tnintcnt,tnintfe,tniter,tnstpcnt,
     +        trtype
      double precision beta,bignum,cgepsf,cgepsi,cggpnf,delmin,epsgpen,
     +        epsgpsn,epsnfp,epsnqmp,eta,etaext,etaint,f,fmin,gamma,
     +        gpeucn2,gpsupn,lspgma,lspgmi,macheps,sigma1,sigma2,theta,
     +        ucgmia,ucgmib,udelta0

C     ARRAY ARGUMENTS
      integer ind(n)
      logical equatn(m),linear(m)
      double precision d(n),g(n),gprev(n),l(n),lambda(m),rho(m),s(n),
     +        u(n),wd1(n),wd2(n),wd3(n),wd4(n),wd5(n),wd6(n),wd7(n),
     +        wd8(n),wd9(n),x(n),xprev(n),y(n)

C     Solves the box-constrained minimization problem
C
C                         Minimize f(x)
C
C                         subject to 
C 
C                                  l <= x <= u
C     
C     using a method described in 
C
C     E. G. Birgin and J. M. Martinez, ''Large-scale active-set box-
C     constrained optimization method with spectral projected 
C     gradients'', Computational Optimization and Applications 23, pp. 
C     101-125, 2002.  
C
C     Description of GENCAN arguments:
C
C     On Entry
C
C     n        integer 
C              number of variables
C
C     x        double precision x(n)
C              initial estimation of the solution
C
C     l        double precision l(n)
C              lower bounds on the variables
C
C     u        double precision u(n)
C              upper bounds on the variables
C
C     m        integer
C     lambda   double precision lambda(m)
C     rho      double precision rho(m)
C     equatn   logical equatn(m)
C     linear   logical linear(m)
C              These five parameters are not used nor modified by 
C              GENCAN and they are passed as arguments to the user-
C              defined subroutines evalal and evalnal to compute the 
C              objective function and its gradient, respectively. 
C              Clearly, in an Augmented Lagrangian context, if GENCAN is 
C              being used to solve the bound-constrained subproblems, m 
C              would be the number of constraints, lambda the Lagrange 
C              multipliers approximation and rho the penalty parameters.
C              equatn is logical array that, for each constraint, 
C              indicates whether the constraint is an equality constraint
C              (.true.) or an inequality constraint (.false.). Finally,
C              linear is logical array that, for each constraint, 
C              indicates whether the constraint is a linear constraint
C              (.true.) or a nonlinear constraint (.false.)
C
C     epsgpen  double precision
C              epsgpen means EPSilon for the Projected Gradient Euclidian
C              Norm. It is a small positive number for declaring 
C              convergence when the Euclidian norm of the continuous 
C              projected gradient is less than or equal to epsgpen
C
C              RECOMMENDED: epsgpen = 1.0d-05
C
C              CONSTRAINTS: epsgpen >= 0.0
C
C     epsgpsn  double precision
C              epsgpsn means EPSilon for the Projected Gradient Sup Norm.
C              It is a small positive number for declaring convergence 
C              when the sup norm of the continuous projected gradient is 
C              less than or equal to epsgpsn
C
C              RECOMMENDED: epsgpsn = 1.0d-05
C
C              CONSTRAINTS: epsgpsn >= 0.0
C
C     maxitnfp integer
C              maxitnfp means MAXimum of ITerations with No Function 
C              Progress. See below for more details.
C
C     epsnfp   double precision
C              epsnfp means EPSilon for No Function Progress. It is a
C              small positive number for declaring ''lack of progress in 
C              the objective function value'' if f(x_k) - f(x_{k+1}) <= 
C              epsnfp * max{ f(x_j) - f(x_{j+1}, j < k } during maxitnfp 
C              consecutive iterations. This stopping criterion may be 
C              inhibited setting maxitnfp equal to maxit.
C
C              RECOMMENDED: maxitnfp = 5 and epsnfp = 1.0d-02
C
C              CONSTRAINTS: maxitnfp >= 1 and epsnfp >= 0.0
C
C     maxitngp integer
C              maxitngp means MAXimum of ITerations with No Gradient
C              Progress. If the order of the Euclidian norm of the 
C              continuous projected gradient did not change during 
C              maxitngp consecutive iterations then the execution stops. 
C
C              RECOMMENDED: maxitngp = 10
C
C              CONSTRAINTS: maxitngp >= 1
C
C     fmin     double precision
C              function value for the stopping criteria f <= fmin
C
C              There is a stopping criterion that stops GENCAN if a 
C              point with a functional value smaller than fmin is found. 
C              The idea behind this stopping criterion is to stop the 
C              method if the objective function is not bounded from 
C              bellow.
C
C              RECOMMENDED: fmin = - bignum
C
C              CONSTRAINTS: there are no constraints for this argument
C
C     maxit    integer
C              maximum number of allowed iterations
C
C              RECOMMENDED: maxit = 1000
C
C              CONSTRAINTS: maxit >= 0
C
C     maxfc    integer
C              maximum allowed number of functional evaluations
C
C              RECOMMENDED: maxfc = 5 * maxit
C
C              CONSTRAINTS: maxfc >= 1
C
C     udelta0  double precision
C              initial ''trust-radius'' for Conjugate Gradients. The 
C              default value max( delmin, 0.1 * max( 1, ||x|| ) ) is 
C              used if the user sets udelta0 <= 0. 
C
C              RECOMMENDED: udelta0 = - 1.0
C
C              CONSTRAINTS: there are no constraints for this argument
C
C     ucgmia   double precision
C              see below
C
C     ucgmib   double precision
C              the maximum allowed number of iterations for each run of 
C              the Conjugate Gradient subalgorithm will be 
C
C              max( 1, int( ucgmia * nind + ucgmib ) ),
C
C              where nind is the number of variables of the subproblem.
C
C              The default value for this maximum number of CG iterations
C              is a linear function of the projected gradient sup-norm.
C              It goes from max( 1, 10 * log( nind ) ) when the method is
C              far from the solution to nind when the method is near to 
C              the solution, where nind is the number of variables of the 
C              subproblem (equal to the number of free variables). 
C
C              The default value will be used if the user sets ucgmia or 
C              ucgmib to any non-positive value. 
C
C              RECOMMENDED: ucgmia = - 1 and ucgmib = - 1
C
C              CONSTRAINTS: there are no constraints for this argument
C
C     cgscre   integer
C              See below
C
C     cggpnf   double precision
C              cgscre means conjugate gradient stopping criterion 
C              relation, and cggpnf means Conjugate Gradients projected 
C              gradient final norm. Both are related to a stopping 
C              criterion of Conjugate Gradients. This stopping criterion 
C              depends on the norm of the residual of the linear system. 
C              The norm of the residual should be less or equal than a 
C              ''small'' quantity which decreases as we are 
C              approximating the solution of the minimization problem 
C              (near the solution, better the truncated-Newton direction 
C              we aim). Then, the log of the required accuracy requested 
C              to Conjugate Gradient has a linear dependence on the log 
C              of the norm of the continuous projected gradient. This 
C              linear relation uses the squared Euclidian norm of the 
C              projected gradient if cgscre is equal to 1 and uses the 
C              sup-norm if cgscre is equal to 2. In addition, the 
C              precision required to CG is equal to cgepsi (conjugate 
C              gradient initial epsilon) at x0 and cgepsf (conjugate 
C              gradient final epsilon) when the Euclidian- or sup-norm 
C              of the projected gradient is equal to cggpnf (conjugate 
C              gradients projected gradient final norm) which is an 
C              estimation of the value of the Euclidian- or sup-norm of 
C              the projected gradient at the solution.
C
C              RECOMMENDED: cgscre = 1, cggpnf = epsgpen; or
C                           cgscre = 2, cggpnf = epsgpsn.
C
C              CONSTRAINTS:  allowed values for cgscre are just 1 or 2
C                            cggpnf >= 0.0
C
C     cgepsi   double precision
C              See below
C
C     cgepsf   double precision
C              small positive numbers for declaring convergence of the 
C              Conjugate Gradients subalgorithm when ||r||_2 < cgeps * 
C              ||rhs||_2, where r is the residual and rhs is the right 
C              hand side of the linear system, i.e., CG stops when the 
C              relative error of the solution is smaller than cgeps. 
C
C              cgeps varies from cgepsi to cgepsf in a way that depends 
C              on cgscre as follows:
C
C              i) CASE cgscre = 1: log10(cgeps^2) depends linearly on 
C              log10(||g_P(x)||_2^2) which varies from ||g_P(x_0)||_2^2 
C              to epsgpen^2
C
C              ii)  CASE cgscre = 2: log10(cgeps) depends linearly on 
C              log10(||g_P(x)||_inf) which varies from ||g_P(x_0)||_inf 
C              to epsgpsn
C
C              RECOMMENDED: cgepsi = 1.0d-01, cgepsf = 1.0d-05
C
C              CONSTRAINTS: cgepsi >= cgepsf >= 0.0
C
C     epsnqmp  double precision
C              See below
C
C     maxitnqmp integer
C              This and the previous argument are used for a stopping 
C              criterion of the Conjugate Gradients subalgorithm. If the 
C              progress in the quadratic model is smaller than fraction 
C              of the best progress ( epsnqmp * bestprog ) during 
C              maxitnqmp consecutive iterations then CG is stopped 
C              declaring ''not enough progress of the quadratic model''.
C
C              RECOMMENDED: epsnqmp = 1.0d-04, maxitnqmp = 5
C
C              CONSTRAINTS: epsnqmp >= 0.0, maxitnqmp >= 1.
C
C     etaint   double precision
C              Constant for the interpolation. See the description of 
C              sigma1 and sigma2 above. Sometimes, in a line search, we 
C              take the new trial step as the previous one divided by 
C              etaint
C
C              RECOMMENDED: etaint = 2.0
C
C              CONSTRAINTS: etaint > 1.0.
C
C     etaext   double precision
C              Constant for the extrapolation. When extrapolating we 
C              try alpha_new = alpha * etaext
C
C              RECOMMENDED: etaext = 2.0
C
C              CONSTRAINTS: etaext > 1.0
C
C     mininterp integer
C              Constant for testing if, after having made at least 
C              mininterp interpolations, the steplength is too small. In
C              that case, failure of the line search is declared (may be 
C              the direction is not a descent direction due to an error 
C              in the gradient calculations). Use mininterp greater 
C              than or equal to maxfc for inhibit this stopping 
C              criterion
C
C              RECOMMENDED: mininterp = 4 
C
C              CONSTRAINTS: mininterp >= 1
C
C     maxextrap integer
C              Constant to limit the number of extrapolations in the 
C              Truncated Newton direction.
C
C              RECOMMENDED: maxextrap = 100 
C
C              CONSTRAINTS: maxextrap >= 0
C
C     gtype    integer
C              gtype indicates in which way the gradient of the 
C              objective function will be computed. See the detailed
C              description of this parameter in subroutine param.
C
C              RECOMMENDED: gtype = 0
C
C              CONSTRAINTS: allowed values are just 0 or 1.
C
C     hptype   integer
C              hptype indicates in which way the Hessian (or the matrix-
C              vector product) will be approximated. See the detailed
C              description of this parameter in subroutine param.
C 
C              RECOMMENDED: hptype = 6
C
C              CONSTRAINTS: any value from 0 to 6 and 9.
C
C     trtype   integer
C              Type of Conjugate Gradients ''trust-radius''. trtype 
C              equal to 0 means Euclidian-norm trust-radius and trtype 
C              equal to 1 means sup-norm trust radius
C
C              RECOMMENDED: trtype = 0
C
C              CONSTRAINTS: allowed values are just 0 or 1.
C
C     iprint   integer
C              Commands printing. Nothing is printed if iprint is 
C              smaller than 2. If iprint is greater than or equal to 
C              2, GENCAN iterations information is printed. If iprint 
C              is greater than or equal to 3, line searches and 
C              Conjugate Gradients information is printed.
C
C              RECOMMENDED: iprint = 2
C
C              CONSTRAINTS: allowed values are just 2 or 3.
C
C     ncomp    integer
C              This constant is just for printing. In a detailed 
C              printing option, ncomp component of some vectors will be 
C              printed
C
C              RECOMMENDED: ncomp = 5
C
C              CONSTRAINTS: ncomp >= 0
C
C     s,y,d,xprev,xgprev double precision s(n),y(n),d(n),xprev(n),gprev(n)
C              n-dimensional working vectors of doubles
C
C     ind      integer ind(n)
C              n-dimensional working vectors of integers
C
C     wd1, ..., wd9 double precision wd1(n), ..., wd9(n)
C              n-dimensional working vectors of doubles
C
C     eta      double precision
C              Constant for deciding abandon the current face or not. We 
C              abandon the current face if the norm of the internal 
C              gradient (here, internal components of the continuous 
C              projected gradient) is smaller than ( 1 - eta ) times the 
C              norm of the continuous projected gradient. Using eta = 
C              0.9 is a rather conservative strategy in the sense that 
C              internal iterations are preferred over SPG iterations. 
C
C              RECOMMENDED: eta = 0.9
C
C              CONSTRAINTS: 0.0 < eta < 1.0
C
C     delmin   double precision
C              Smaller Conjugate Gradients ''trust radius'' to compute 
C              the Truncated Newton direction
C
C              RECOMMENDED: delmin = 0.1
C
C              CONSTRAINTS: delmin > 0.0
C
C     lspgmi   double precision
C              See below
C
C     lspgma   double precision
C              The spectral steplength, called lamspg, is projected onto 
C              the box [lspgmi,lspgma] 
C
C              RECOMMENDED: lspgmi = 1.0d-10 and lspgma = 1.0d+10
C 
C              CONSTRAINTS: lspgma >= lspgmi > 0.0
C
C     theta    double precision
C              Constant for the angle condition, i.e., at iteration k we 
C              need a direction dk such that <gk,dk> <= - theta 
C              ||gk||_2 ||dk||_2, where gk is \nabla f(xk)
C
C              RECOMMENDED: theta = 10^{-6}
C
C              CONSTRAINTS: 0.0 < theta < 1.0
C
C     gamma    double precision
C              Constant for the Armijo criterion
C              f(x + alpha d) <= f(x) + gamma * alpha * <g,d>
C
C              RECOMMENDED: gamma = 1.0d-04
C
C              CONSTRAINTS: 0.0 < gamma < 0.5.
C
C     beta     double precision
C              Constant for the beta condition <dk, g(xk + dk)>  < beta 
C              * <dk,gk>. If (xk + dk) satisfies the Armijo condition 
C              but does not satisfy the beta condition then the point is 
C              accepted, but if it satisfied the Armijo condition and 
C              also satisfies the beta condition then we know that there 
C              is the possibility for a successful extrapolation
C
C              RECOMMENDED: beta = 0.5
C
C              CONSTRAINTS: 0.0 < beta < 1.0.
C
C     sigma1   double precision
C              See below
C
C     sigma2   double precision
C              Constant for the safeguarded interpolation. If alpha_new 
C              is not inside the interval [sigma1, sigma * alpha] then 
C              we take alpha_new = alpha / etaint
C
C              RECOMMENDED: sigma1 = 0.1 and sigma2 = 0.9
C
C              CONSTRAINTS: 0 < sigma1 < sigma2 < 1.
C
C     On Return
C
C     x        double precision x(n)
C              final estimation to the solution
C
C     f        double precision
C              function value at the final estimation 
C
C     g        double precision g(n)
C              gradient at the final estimation
C
C     gpeucn2  double precision
C              squared euclidian-norm of the continuous projected 
C              gradient at the final estimation
C
C     gpsupn   double precision
C              sup-norm of the continuous projected gradient at the 
C              final estimation
C
C     iter     integer
C              number of iterations
C
C     fcnt     integer
C              number of function evaluations   
C
C     gcnt     integer
C              number of gradient evaluations   
C
C     cgcnt    integer
C              number of Conjugate Gradients iterations   
C
C     spgiter  integer
C              number of Spectral Projected Gradient iterations
C
C     spgfcnt  integer
C              number of functional evaluations along Spectral Projected
C              Gradient directions
C
C     tniter   integer
C              number of Truncated-Newton iterations
C
C     tnfcnt   integer
C              number of functional evaluations along Truncated-Newton
C              directions
C
C     tnintcnt integer
C              number of times a backtracking in a Truncated-Newton
C              direction was needed
C
C     tnexgcnt integer
C              number of times an extrapolation in a Truncated-Newton
C              direction successfully decreased the objective funtional
C              value
C
C     tnexbcnt integer
C              number of times an extrapolation was aborted in the first
C              extrapolated point by an increase in the objective 
C              functional value
C
C     tnstpcnt integer
C              number of times the Newton point was accepted (without
C              interpolations nor extrapolations)
C
C     tnintfe  integer
C              number of functional evaluations used in interpolations 
C              along Truncated-Newton directions
C
C     tnexgfe  integer
C              number of functional evaluations used in successful 
C              extrapolations along Truncated-Newton directions
C
C     tnexbfe  integer
C              number of functional evaluations used in unsuccessful 
C              extrapolations along Truncated-Newton directions
C
C     inform   integer
C              This output parameter tells what happened in this 
C              subroutine, according to the following conventions:
C 
C              0 = convergence with small Euclidian norm of the 
C                  continuous projected gradient (smaller than epsgpen);
C
C              1 = convergence with small sup-norm of the continuous 
C                  projected gradient (smaller than epsgpsn);
C
C              2 = the algorithm stopped by ''lack of progress'', that 
C                  means that f(xk) - f(x_{k+1}) <= epsnfp * 
C                  max{ f(x_j) - f(x_{j+1}, j < k } during maxitnfp 
C                  consecutive iterations. If desired, set maxitnfp 
C                  equal to maxit to inhibit this stopping criterion.
C
C              3 = the algorithm stopped because the order of the 
C                  Euclidian norm of the continuous projected gradient 
C                  did not change during maxitngp consecutive 
C                  iterations. Probably, we are asking for an 
C                  exaggerated small norm of continuous projected 
C                  gradient for declaring convergence. If desired, set
C                  maxitngp equal to maxit to inhibit this stopping 
C                  criterion.
C
C              4 = the algorithm stopped because the functional value 
C                  is very small (smaller than fmin). If desired, set 
C                  fmin equal to minus bignum to inhibit this stopping 
C                  criterion.
C
C              6 = too small step in a line search. After having made at 
C                  least mininterp interpolations, the steplength 
C                  becames small. ''small steplength'' means that we are 
C                  at point x with direction d and step alpha, such that
C
C                  | alpha * d(i) | <= macheps * max( | x(i) |, 1 ) 
C
C                  for all i.
C
C                  In that case failure of the line search is declared 
C                  (may be the direction is not a descent direction due 
C                  to an error in the gradient calculations). If 
C                  desired, set mininterp equal to maxfc to inhibit this 
C                  stopping criterion.
C
C              7 = it was achieved the maximum allowed number of 
C                  iterations (maxit);
C
C              8 = it was achieved the maximum allowed number of 
C                  function evaluations (maxfc);
C
C            < 0 = error in evalal, evalnal or evalhd subroutines.

C     DATA BLOCKS
      character * 41 ittext(6)
      data ittext(1) /'(Used an SPG iteration)                  '/
      data ittext(2) /'(Used TN with the ANALYTIC HESSIAN)      '/
      data ittext(3) /'(Used TN with the USER HP PRODUCT)       '/
      data ittext(4) /'(Used TN with a QN HESSIAN APPROX)       '/
      data ittext(5) /'(Used TN with INCREMENTAL QUOTIENTS)     '/
      data ittext(6) /'(Used TN with PURE INCREMENTAL QUOTIENTS)'/

C     LOCAL SCALARS
      logical samefa
      character * 6 aptype
      integer cgiter,cgmaxit,fcntprev,i,infotmp,ittype,itnfp,itngp,nind,
     +        nprint,rbdind,rbdtype,tnexbprev,tnexgprev,tnintprev
      double precision acgeps,amax,amaxx,bestprog,bcgeps,cgeps,currprog,
     +        delta,epsgpen2,fprev,gieucn2,gpi,gpsupnprev,lamspg,ometa2,
     +        seucn,ssupn,sts,sty,yeucn,xeucn2,xsupn,tsmall
C             gpeucn20,gpsupn0


C     ==================================================================
C     Initialization
C     ==================================================================

C     Set some initial values:

C     just for printing,
      nprint   = min0( n, ncomp )

C     for testing convergence,
      epsgpen2 = epsgpen ** 2

C     for testing whether to abandon the current face or not
C     (ometa2 means '(one minus eta) squared'),
      ometa2   = ( 1.0d0 - eta ) ** 2

C     for testing progress in f and decrease of the sup-norm of the 
C     projected gradient,
      fprev      = bignum
      gpsupnprev = bignum
      bestprog   =  0.0d0
      itnfp      =      0
      itngp      =      0

C     counters.
      iter       =  0
      fcnt       =  0
      gcnt       =  0
      cgcnt      =  0

      spgiter    =  0
      spgfcnt    =  0

      tniter     =  0
      tnfcnt     =  0

      tnstpcnt   =  0
      tnintcnt   =  0
      tnexgcnt   =  0
      tnexbcnt   =  0

      tnintfe    =  0
      tnexgfe    =  0
      tnexbfe    =  0

C     Print problem information

      if( iprint .ge. 2 ) then
          write(*, 977) n
          write(*, 978) nprint,(l(i),i=1,nprint)
          write(*, 979) nprint,(u(i),i=1,nprint)
          write(*, 980) nprint,(x(i),i=1,nprint)

          write(10,977) n
          write(10,978) nprint,(l(i),i=1,nprint)
          write(10,979) nprint,(u(i),i=1,nprint)
          write(10,980) nprint,(x(i),i=1,nprint)
      end if

C     Project initial guess. If the initial guess is infeasible, 
C     projection puts it into the box.

      do i = 1,n
          x(i) = max( l(i), min( x(i), u(i) ) )
      end do

C     Compute x norms

      xeucn2 = 0.0d0
      xsupn = 0.0d0
      do i = 1,n
          xeucn2 = xeucn2 + x(i) ** 2
          xsupn = max( xsupn, abs( x(i) ) )
      end do

C     Compute function and gradient at the initial point

      call evalal(n,x,m,lambda,rho,equatn,linear,f,inform)
      fcnt = fcnt + 1

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 2 ) then
              write(*, 1000) inform
              write(10,1000) inform
          end if

          return
      end if

      call evalnal(n,x,m,lambda,rho,equatn,linear,g,gtype,macheps,
     +inform)
      gcnt = gcnt + 1

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 2 ) then
              write(*, 1000) inform
              write(10,1000) inform
          end if

          return
      end if

C     Compute continuous-project-gradient Euclidian and Sup norms,
C     internal gradient Euclidian norm, and store in nind the number of
C     free variables and in array ind their identifiers. Also set
C     variable samefa (same face) indicating that the "previous iterate"
C     does not belong to the current face.

      nind    = 0
      gpsupn  = 0.0d0
      gpeucn2 = 0.0d0
      gieucn2 = 0.0d0
      do i = 1,n
          gpi     = min( u(i), max( l(i), x(i) - g(i) ) ) - x(i)
          gpsupn  = max( gpsupn, abs( gpi ) )
          gpeucn2 = gpeucn2 + gpi ** 2
          if ( x(i) .gt. l(i) .and. x(i) .lt. u(i) ) then
              gieucn2   = gieucn2 + gpi ** 2
              nind      = nind + 1
              ind(nind) = i
          end if
      end do

      samefa = .false.

C     Initial spectral steplength

C     Compute a small step and set the point at which the auxiliary
C     gradient will be computed

      tsmall = sqrt( macheps ) * max( xsupn / gpsupn, 1.0d0 )

      do i = 1,n
          gpi  = min( u(i), max( l(i), x(i) - g(i) ) ) - x(i)
          s(i) = x(i) + tsmall * gpi
      end do

      call setpoint(s)

C     Compute the gradient at the auxiliary point

      call evalnalu(n,s,m,lambda,rho,equatn,linear,y,gtype,macheps,
     +inform)
      gcnt = gcnt + 1

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 2 ) then
              write(*, 1000) inform
              write(10,1000) inform
          end if

          return
      end if

C     Compute s = x_aux - x_0, y = g_aux - g_0, <s,s> and <s,y>

      sts = 0.0d0
      sty = 0.0d0
      ssupn = 0.0d0
      yeucn = 0.0d0
      do i = 1,n
          s(i) = s(i) - x(i)
          y(i) = y(i) - g(i)
          sts  = sts + s(i) ** 2
          sty  = sty + s(i) * y(i)
          ssupn = max( ssupn, abs( s(i) ) )
          yeucn = yeucn + y(i) ** 2
      end do
      seucn = sqrt( sts )
      yeucn = sqrt( yeucn )

C     Compute a linear relation between gpeucn2 and cgeps2, i.e.,
C     scalars a and b such that 
c
C         a * log10(||g_P(x_0)||_2^2) + b = log10(cgeps_0^2) and
c
C         a * log10(||g_P(x_f)||_2^2) + b = log10(cgeps_f^2),
c
C     where cgeps_0 and cgeps_f are provided. Note that if 
C     cgeps_0 is equal to cgeps_f then cgeps will be always 
C     equal to cgeps_0 and cgeps_f.

C     We introduce now a linear relation between gpsupn and cgeps also.

      if ( cgscre .eq. 1 ) then
          acgeps = 2.0d0 * log10( cgepsf / cgepsi ) / 
     +                     log10( cggpnf ** 2 / gpeucn2 )
          bcgeps = 2.0d0 * log10( cgepsi ) - acgeps * log10( gpeucn2 )
      else ! if ( cgscre .eq. 2 ) then
          acgeps = log10( cgepsf / cgepsi ) / log10( cggpnf / gpsupn )
          bcgeps = log10( cgepsi ) - acgeps * log10( gpsupn )
      end if 

C     And it will be used for the linear relation of cgmaxit

C     gpsupn0  = gpsupn
C     gpeucn20 = gpeucn2

C     Print initial information

      if( iprint .ge. 2 ) then
          write(*, 981) iter
          write(*, 985) nprint,(x(i),i=1,nprint)
          write(*, 986) nprint,(g(i),i=1,nprint)
          write(*, 987) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),i=1,
     +    nprint)
          write(*, 988) min0(nprint,nind),nind,(ind(i),i=1,min0(nprint,
     +    nind))
          write(*, 1002) f,sqrt(gpeucn2),sqrt(gieucn2),gpsupn,nind,n,
     +    spgiter,tniter,fcnt,gcnt,cgcnt

          write(10,981) iter
          write(10,985) nprint,(x(i),i=1,nprint)
          write(10,986) nprint,(g(i),i=1,nprint)
          write(10,987) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),i=1,
     +    nprint)
          write(10,988) min0(nprint,nind),nind,(ind(i),i=1,min0(nprint,
     +    nind))
          write(10,1002) f,sqrt(gpeucn2),sqrt(gieucn2),gpsupn,nind,n,
     +    spgiter,tniter,fcnt,gcnt,cgcnt
      end if

C     SAVING INTERMEDIATE DATA FOR CRASH REPORT

      open(20,file='gencan-tabline.out')
      write(20,3000) 0.0d0,-1,n,0,0,iter,fcnt,gcnt,cgcnt,f,0.0d0,gpsupn
      close(20)

C     ==================================================================
C     Main loop
C     ==================================================================

 100  continue

C     ==================================================================
C     Test stopping criteria
C     ==================================================================

C     Test whether the continuous-projected-gradient Euclidian norm
C     is small enough to declare convergence

      if ( gpeucn2 .le. epsgpen2 ) then
          inform = 0

          if ( iprint .ge. 2 ) then
              write(*, 990) inform,epsgpen
              write(10,990) inform,epsgpen
          end if

          go to 500
      end if

C     Test whether the continuous-projected-gradient Sup norm
C     is small enough to declare convergence

      if ( gpsupn .le. epsgpsn ) then
          inform = 1

          if ( iprint .ge. 2 ) then
              write(*, 991) inform,epsgpsn
              write(10,991) inform,epsgpsn
          end if

          go to 500
      end if

C     Test whether we performed many iterations without good progress
C     of the functional value

      currprog = fprev - f
      bestprog = max( currprog, bestprog )

      if ( currprog .le. epsnfp * bestprog ) then

          itnfp = itnfp + 1

          if ( itnfp .ge. maxitnfp ) then
              inform = 2

              if ( iprint .ge. 2 ) then
                  write(*, 992) inform,epsnfp,maxitnfp
                  write(10,992) inform,epsnfp,maxitnfp
              end if

              go to 500
          end if

      else
          itnfp = 0
      end if

C     Test whether we have performed many iterations without good 
C     reduction of the sup-norm of the projected gradient

      if ( gpsupn .ge. gpsupnprev ) then

          itngp = itngp + 1

          if ( itngp .ge. maxitngp ) then
              inform = 3

              if ( iprint .ge. 2 ) then
                  write(*, 993) inform,maxitngp
                  write(10,993) inform,maxitngp
              end if

              go to 500
          end if

      else
          itngp = 0
      end if

C     Test whether the functional value is very small

      if ( f .le. fmin ) then

          inform = 4

          if ( iprint .ge. 2 ) then
              write(*, 994) inform,fmin
              write(10,994) inform,fmin
          end if

          go to 500

      end if

C     Test whether the number of iterations is exhausted

      if ( iter .ge. maxit ) then

          inform = 7

          if ( iprint .ge. 2 ) then
              write(*, 997) inform,maxit
              write(10,997) inform,maxit
          end if

          go to 500

      end if

C     Test whether the number of functional evaluations is exhausted

      if ( fcnt .ge. maxfc ) then

          inform = 8

          if ( iprint .ge. 2 ) then
              write(*, 998) inform,maxfc
              write(10,998) inform,maxfc
          end if

          go to 500

      end if

C     ==================================================================
C     The stopping criteria were not satisfied, a new iteration will be 
C     made
C     ==================================================================

      iter = iter + 1

C     ==================================================================
C     Save current values, f, x and g
C     ==================================================================

      fprev = f
      gpsupnprev = gpsupn

      do i = 1,n
          xprev(i) = x(i)
          gprev(i) = g(i)
      end do

C     ==================================================================
C     Compute new iterate
C     ==================================================================

C     We abandon the current face if the norm of the internal gradient
C     (here, internal components of the continuous projected gradient)
C     is smaller than (1-eta) times the norm of the continuous 
C     projected gradient. Using eta=0.9 is a rather conservative 
C     strategy in the sense that internal iterations are preferred over 
C     SPG iterations. Replace eta = 0.9 by other tolerance in (0,1) if 
C     you find it convenient. 

      if ( gieucn2 .le. ometa2 * gpeucn2 ) then

C         ==============================================================
C         Some constraints should be abandoned. Compute the new iterate 
C         using an SPG iteration
C         ==============================================================

          spgiter = spgiter + 1

C         Set iteration type

          ittype = 1

C         Compute spectral steplength

          if ( sty .le. 0.0d0 ) then
              lamspg = max( 1.0d0, sqrt( xeucn2 ) ) / sqrt( gpeucn2 )
          else
              lamspg = sts / sty
          end if
          lamspg = min( lspgma, max( lspgmi, lamspg ) )

C         Perform a line search with safeguarded quadratic interpolation 
C         along the direction of the spectral continuous projected 
C         gradient

          fcntprev = fcnt

          call spgls(n,x,m,lambda,rho,equatn,linear,f,g,l,u,lamspg,
     +    etaint,mininterp,fmin,maxfc,iprint,fcnt,inform,wd1,wd2,gamma,
     +    sigma1,sigma2,macheps,bignum) 

          spgfcnt = spgfcnt + ( fcnt - fcntprev ) 

          if ( inform .lt. 0 ) then

              if ( iprint .ge. 2 ) then
                  write(*, 1000) inform
                  write(10,1000) inform
              end if

              return
          end if

C         Compute the gradient at the new iterate

          call evalnal(n,x,m,lambda,rho,equatn,linear,g,gtype,macheps,
     +    inform)
          gcnt = gcnt + 1
 
          if ( inform .lt. 0 ) then

              if ( iprint .ge. 2 ) then
                  write(*, 1000) inform
                  write(10,1000) inform
              end if

              return
          end if

      else

C         ==============================================================
C         The new iterate will belong to the closure of the current face
C         ==============================================================

          tniter = tniter + 1

C         Compute trust-region radius

          if ( iter .eq. 1 ) then
              if( udelta0 .le. 0.0d0 ) then
                  if ( trtype .eq. 0 ) then
                      delta = max( delmin, 
     +                        0.1d0 * max( 1.0d0, sqrt( xeucn2 ) ) )
                  else ! if ( trtype .eq. 1 ) then
                      delta = max( delmin, 0.1d0 * max( 1.0d0, xsupn ) )
                  end if
              else
                  delta = max( delmin, udelta0 )
              end if
          else
              if ( trtype .eq. 0 ) then
                  delta = max( delmin, 10.0d0 * sqrt( sts ) )
              else ! if ( trtype .eq. 1 ) then
                  delta = max( delmin, 10.0d0 * ssupn )
              end if
          end if

C         Shrink the point, its gradient and the bounds

          call shrink(nind,ind,n,x)
          call shrink(nind,ind,n,g)
          call shrink(nind,ind,n,l)
          call shrink(nind,ind,n,u)

C         Compute the descent direction solving the newtonian system by 
C         conjugate gradients

C         Set conjugate gradient stopping criteria. Default values are 
C         taken if you set ucgeps < 0 or ucgmia < 0 or ucgmib < 0. 
C         Otherwise, the parameters cgeps and cgmaxit will be the ones 
C         set by the user.

          if( ucgmia .lt. 0.0d0 .or. ucgmib .lt. 0.0d0 ) then
c             if ( cgscre .eq. 1 ) then
c                 kappa = log10( gpeucn2 / gpeucn20 )/
c    +                    log10( epsgpen2 / gpeucn20 )
c             else ! if ( cgscre .eq. 2 ) then
c                 kappa= log10( gpsupn / gpsupn0 ) / 
c    +                   log10( epsgpsn / gpsupn0 )
c             end if
c             kappa = max( 0.0d0, min( 1.0d0, kappa ) )
c             cgmaxit = int( ( 1.0d0 - kappa ) * max( 1.0d0, 
c    +        min( dfloat( nind ), 10.0d0 * 
c    +        log10( dfloat( nind ) ) ) ) + kappa * dfloat( nind ) )
              cgmaxit = max( 1, min( 2 * nind, 10000 ) )
          else
              cgmaxit = max( 1, int( ucgmia * nind + ucgmib ) )
          end if

          if ( cgscre .eq. 1 ) then
              cgeps = sqrt( 10.0d0 ** ( acgeps * log10( gpeucn2 ) + 
     +                bcgeps ) )
          else ! if ( cgscre .eq. 2 ) then
              cgeps = 10.0d0 ** ( acgeps * log10( gpsupn ) + bcgeps )
          end if
          cgeps = max( cgepsf, min( cgepsi, cgeps ) )

C         Call conjugate gradients

          call cgm(nind,ind,n,x,m,lambda,rho,equatn,linear,g,delta,l,u,
     +    cgeps,epsnqmp,maxitnqmp,cgmaxit,gtype,hptype,trtype,precond,
     +    samefa,aptype,s,y,ssupn,seucn,yeucn,sts,sty,lspgmi,lspgma,
     +    iprint,ncomp,d,cgiter,rbdtype,rbdind,inform,wd1,wd2,wd3,wd4,
     +    wd5,wd6,wd7,wd8,wd9,theta,macheps,bignum)

	  cgcnt = cgcnt + cgiter

          if ( inform .lt. 0 ) then

              if ( iprint .ge. 2 ) then
                  write(*, 1000) inform
                  write(10,1000) inform
              end if

              return

          end if

C         Set iteration type

          if ( aptype .eq. 'TRUEHE' ) then
              ittype = 2
          else if ( aptype .eq. 'HLPROD' ) then
              ittype = 3
          else if ( aptype .eq. 'QNCGNA' ) then
              ittype = 4
          else if ( aptype .eq. 'INCQUO' ) then
              ittype = 5
          else if ( aptype .eq. 'PUREIQ' ) then
              ittype = 6
          end if

C         Compute maximum step

          if ( inform .eq. 2 ) then
              amax = 1.0d0
          else
              amax = bignum
              rbdtype = 0
              do i = 1,nind
                  if ( d(i) .gt. 0.0d0 ) then
                      amaxx = ( u(i) - x(i) ) / d(i)
                      if ( amaxx .lt. amax ) then
                          amax    = amaxx
                          rbdind  = i
                          rbdtype = 2
                      end if
                  else if ( d(i) .lt. 0.0d0 ) then
                      amaxx = ( l(i) - x(i) ) / d(i)
                      if ( amaxx .lt. amax ) then
                          amax    = amaxx
                          rbdind  = i
                          rbdtype = 1
                      end if
                  end if
               end do
          end if

C         Perform the line search

          tnintprev = tnintcnt
          tnexgprev = tnexgcnt
          tnexbprev = tnexbcnt

          fcntprev = fcnt

          call tnls(nind,ind,n,x,m,lambda,rho,equatn,linear,l,u,f,g,d,
     +    amax,rbdtype,rbdind,etaint,etaext,mininterp,maxextrap,fmin,
     +    maxfc,gtype,iprint,fcnt,gcnt,tnintcnt,tnexgcnt,tnexbcnt,
     +    inform,wd1,wd2,wd3,gamma,beta,sigma1,sigma2,macheps,bignum)

          if ( inform .lt. 0 ) then

              if ( iprint .ge. 2 ) then
                  write(*, 1000) inform
                  write(10,1000) inform
              end if

              return

          end if

          if ( tnintcnt .gt. tnintprev ) then
              tnintfe = tnintfe + ( fcnt - fcntprev )
          else if ( tnexgcnt .gt. tnexgprev ) then
              tnexgfe = tnexgfe + ( fcnt - fcntprev )
          else if ( tnexbcnt .gt. tnexbprev ) then
              tnexbfe = tnexbfe + ( fcnt - fcntprev )
          else
              tnstpcnt = tnstpcnt + 1
          end if

          tnfcnt = tnfcnt + ( fcnt - fcntprev )

C         Expand the point, its gradient and the bounds

          call expand(nind,ind,n,x)
          call expand(nind,ind,n,g)
          call expand(nind,ind,n,l)
          call expand(nind,ind,n,u)

C         If the line search (interpolation) in the Truncated Newton
C         direction stopped due to a very small step (inform = 6), we 
C         will discard this iteration and force a SPG iteration

C         Note that tnls subroutine was coded in such a way that in case
C         of inform = 6 termination the subroutine discards all what was 
C         done and returns with the same point it started

          if ( inform .eq. 6 ) then

              if ( iprint .ge. 2 ) then
                  write(*,*)  
                  write(*,*)  
     +            '     The previous Newtonian iteration was discarded',
     +            '     due to a termination for very small step in   ',
     +            '     the line search. A SPG iteration will be      ',
     +            '     forced now.                                   '

                  write(10,*)  
                  write(10,*)  
     +            '     The previous Newtonian iteration was discarded',
     +            '     due to a termination for very small step in   ',
     +            '     the line search. A SPG iteration will be      ',
     +            '     forced now.                                   '
              end if

              spgiter = spgiter + 1

C             Set iteration type

              ittype = 1

C             Compute spectral steplength

              if ( sty .le. 0.0d0 ) then
                  lamspg = max( 1.0d0, sqrt( xeucn2 ) ) / sqrt(gpeucn2)
              else
                  lamspg = sts / sty
              end if
              lamspg = min( lspgma, max( lspgmi, lamspg ) )

C             Perform a line search with safeguarded quadratic 
C             interpolation along the direction of the spectral 
C             continuous projected gradient

              fcntprev = fcnt

              call spgls(n,x,m,lambda,rho,equatn,linear,f,g,l,u,lamspg,
     +        etaint,mininterp,fmin,maxfc,iprint,fcnt,inform,wd1,wd2,
     +        gamma,sigma1,sigma2,macheps,bignum) 

              spgfcnt = spgfcnt + ( fcnt - fcntprev )

              if ( inform .lt. 0 ) then

                  if ( iprint .ge. 2 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return
              end if

C             Compute the gradient at the new iterate

              infotmp = inform

              call evalnal(n,x,m,lambda,rho,equatn,linear,g,gtype,
     +        macheps,inform)
              gcnt = gcnt + 1
 
              if ( inform .lt. 0 ) then

                  if ( iprint .ge. 2 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return
              end if

              inform = infotmp

          end if

      end if

C     ==================================================================
C     Prepare for the next iteration 
C     ==================================================================

C     This adjustment/projection is ''por lo que las putas pudiera''

      do i = 1,n
          if ( x(i) .le. l(i) + 
     +    macheps * max( abs( l(i) ), 1.0d0 ) ) then
              x(i) = l(i)
          else if (x(i). ge. u(i) - 
     +    macheps * max( abs( u(i) ), 1.0d0 ) ) then  
              x(i) = u(i)
          end if
      end do

C     Compute x Euclidian norm

      xeucn2 = 0.0d0
      xsupn = 0.0d0
      do i = 1,n
          xeucn2 = xeucn2 + x(i) ** 2
          xsupn = max( xsupn, abs( x(i) ) )
      end do

C     Compute continuous-project-gradient Euclidian and Sup norms,
C     internal gradient Euclidian norm, and store in nind the number of
C     free variables and in array ind their identifiers.

      nind    = 0
      gpsupn  = 0.0d0
      gpeucn2 = 0.0d0
      gieucn2 = 0.0d0
      do i = 1,n
          gpi     = min( u(i), max( l(i), x(i) - g(i) ) ) - x(i)
          gpsupn  = max( gpsupn, abs( gpi ) )
          gpeucn2 = gpeucn2 + gpi ** 2
          if ( x(i) .gt. l(i) .and. x(i) .lt. u(i) ) then
              gieucn2   = gieucn2 + gpi ** 2
              nind      = nind + 1
              ind(nind) = i
          end if
      end do

C     Verify whether the previous point xprev belongs to the current 
C     face or not

      samefa = .true.
      do i = 1,n
          if ( ( ( x(i) .eq. l(i) .or. xprev(i) .eq. l(i) ) .and. 
     +             x(i) .ne. xprev(i) ) .or.
     +         ( ( x(i) .eq. u(i) .or. xprev(i) .eq. u(i) ) .and. 
     +             x(i) .ne. xprev(i) ) ) then
              samefa = .false.
          end if
      end do

C     Compute s = x_{k+1} - x_k, y = g_{k+1} - g_k, <s,s> and <s,y>

      sts = 0.0d0
      sty = 0.0d0
      ssupn = 0.0d0
      yeucn = 0.0d0
      do i = 1,n
          s(i) = x(i) - xprev(i)
          y(i) = g(i) - gprev(i)
          sts  = sts + s(i) ** 2
          sty  = sty + s(i) * y(i)
          ssupn = max( ssupn, abs( s(i) ) )
          yeucn = yeucn + y(i) ** 2
      end do
      seucn = sqrt( sts )
      yeucn = sqrt( yeucn )

C     Print information of this iteration

      if ( iprint .ge. 2 ) then 
          write(*, 983) iter,ittext(ittype)
          write(*, 985) nprint,(x(i),i=1,nprint)
          write(*, 986) nprint,(g(i),i=1,nprint)
          write(*, 987) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),i=1,
     +    nprint)
          write(*, 988) min0(nprint,nind),nind,(ind(i),i=1,min0(nprint,
     +    nind))
          write(*, 1002) f,sqrt(gpeucn2),sqrt(gieucn2),gpsupn,nind,n,
     +    spgiter,tniter,fcnt,gcnt,cgcnt

          write(10,983) iter,ittext(ittype)
          write(10,985) nprint,(x(i),i=1,nprint)
          write(10,986) nprint,(g(i),i=1,nprint)
          write(10,987) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),i=1,
     +    nprint)
          write(10,988) min0(nprint,nind),nind,(ind(i),i=1,min0(nprint,
     +    nind))
          write(10,1002) f,sqrt(gpeucn2),sqrt(gieucn2),gpsupn,nind,n,
     +    spgiter,tniter,fcnt,gcnt,cgcnt
      end if

C     SAVING INTERMEDIATE DATA FOR CRASH REPORT

      open(20,file='gencan-tabline.out')
      write(20,3000) 0.0d0,-1,n,0,0,iter,fcnt,gcnt,cgcnt,f,0.0d0,gpsupn
      close(20)

C     ==================================================================
C     Test some stopping criteria that may occur inside the line 
C     searches 
C     ==================================================================

      if ( inform .eq. 6 ) then

          if ( iprint .ge. 2 ) then
              write(*, 996) inform,mininterp
              write(10,996) inform,mininterp
          end if

          go to 500

      end if

C     ==================================================================
C     Iterate 
C     ==================================================================

      go to 100

C     ==================================================================
C     End of main loop
C     ==================================================================

C     ==================================================================
C     Report output status and return
C     ==================================================================

 500  continue

C     Print final information

      if ( iprint .ge. 2 ) then
          write(*, 982) iter
          write(*, 985) nprint,(x(i),i=1,nprint)
          write(*, 986) nprint,(g(i),i=1,nprint)
          write(*, 987) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),i=1,
     +    nprint)
          write(*, 988) min0(nprint,nind),nind,(ind(i),i=1,min0(nprint,
     +    nind))
          write(*, 1002) f,sqrt(gpeucn2),sqrt(gieucn2),gpsupn,nind,n,
     +    spgiter,tniter,fcnt,gcnt,cgcnt

          write(10,982) iter
          write(10,985) nprint,(x(i),i=1,nprint)
          write(10,986) nprint,(g(i),i=1,nprint)
          write(10,987) nprint,(min(u(i),max(l(i),x(i)-g(i)))-x(i),i=1,
     +    nprint)
          write(10,988) min0(nprint,nind),nind,(ind(i),i=1,min0(nprint,
     +    nind))
          write(10,1002) f,sqrt(gpeucn2),sqrt(gieucn2),gpsupn,nind,n,
     +    spgiter,tniter,fcnt,gcnt,cgcnt
      end if

      return 

C     Non-executable statements

 977  format(/1X, 'Entry to GENCAN. Number of variables: ',I7)
 978  format(/1X,'Lower bounds (first ',I6, ' components): ',
     */,6(1X,1PD11.4))
 979  format(/1X,'Upper bounds (first ',I6, ' components): ',
     */,6(1X,1PD11.4))
 980  format(/1X,'Initial point (first ',I6, ' components): ',
     */,6(1X,1PD11.4))
 981  format(/1X,'GENCAN iteration: ',I6,' (Initial point)')
 982  format(/1X,'GENCAN iteration: ',I6,' (Final point)')
 983  format(/1X,'GENCAN iteration: ',I6,1X,A41)
 985  format(1X,'Current point (first ',I6, ' components): ',
     */,6(1X,1PD11.4))
 986  format(1X,'Current gradient (first ',I6, ' components): ',
     */,6(1X,1PD11.4))
 987  format(1X,'Current continuous projected gradient (first ',I6, 
     *' components): ',/,6(1X,1PD11.4))
 988  format(1X,'Current free variables (first ',I6,
     *', total number ',I6,'): ',/,10(1X,I6))
 990  format(/1X,'Flag of GENCAN = ',I3,
     *' (convergence with Euclidian-norm of the projected',
     */,1X,'gradient smaller than ',1PD11.4,')')
 991  format(/1X,'Flag of GENCAN = ',I3,
     *' (convergence with sup-norm of the projected gradient',
     */,1X,'smaller than ',1PD11.4,')')
 992  format(/1X,'Flag of GENCAN= ',I3,
     *' (The algorithm stopped by lack of enough progress. This means',
     */,1X,'that  f(x_k) - f(x_{k+1}) .le. ',1PD11.4,
     *' * max [ f(x_j)-f(x_{j+1}, j < k ]',/,1X,'during ',I7,
     *' consecutive iterations')
 993  format(/1X,'Flag of GENCAN = ',I3,
     *' (The algorithm stopped because the order of the',
     */,1X,'Euclidian-norm of the continuous projected gradient did',
     *' not change during ',/,1X,I7,' consecutive iterations.',
     *' Probably, an exaggerated small norm of the',/,1X,'continuous',
     *' projected gradient is required for declaring convergence')
 994  format(/1X,'Flag of GENCAN = ',I3,
     *' (The algorithm stopped because the functional value is',
     */,1X,'smaller than ',1PD11.4)
 996  format(/1X,'Flag of GENCAN = ',I3,
     *' (After having made at least ',I7,' interpolations, the',/,
     *' line search step became very small.')
 997  format(/1X,'Flag of GENCAN = ',I3,
     *' (It was exceeded the maximum allowed number of iterations',
     */,1X,'(maxit=',I7,')')
 998  format(/1X,'Flag of GENCAN = ',I3,
     *' (It was exceeded the maximum allowed number of functional',
     */,1X,'evaluations (maxfc=',I7,')')
 1002 format(1X,'Functional value: ', 1PD11.4,
     */,1X,'Euclidian-norm of the continuous projected gradient: ',
     *1PD11.4,
     */,1X,'Euclidian-norm of the internal projection of gp: ',1PD11.4,
     */,1X,'Sup-norm of the continuous projected gradient: ',1PD11.4,
     */,1X,'Free variables at this point: ',I7,
     *' (over a total of ',I7,')',
     */,1X,'SPG iterations: ',I7,
     */,1X,'TN iterations: ',I7,
     */,1X,'Functional evaluations: ',I7,
     */,1X,'Gradient evaluations: ',I7,
     */,1X,'Conjugate gradient iterations: ',I7)
 1000 format(/1X,'Flag of GENCAN = ',I3,' Fatal Error')
 3000 format(F8.2,1X,I3,1X,I6,1X,I6,1X,I3,1X,I7,1X,I7,1X,I7,1X,I7,1X,1P
     +       D24.16,1X,1P,D7.1,1X,1P,D7.1)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine spgls(n,x,m,lambda,rho,equatn,linear,f,g,l,u,lamspg,
     +etaint,mininterp,fmin,maxfc,iprint,fcnt,inform,xtrial,d,gamma,
     +sigma1,sigma2,macheps,bignum) 

      implicit none

C     SCALAR ARGUMENTS
      integer fcnt,m,maxfc,mininterp,n,inform,iprint
      double precision bignum,f,fmin,gamma,lamspg,macheps,etaint,sigma1,
     +        sigma2

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      double precision d(n),g(n),l(n),lambda(m),rho(m),u(n),x(n),
     +        xtrial(n)
 
C     Safeguarded quadratic interpolation, used in the Spectral 
C     Projected Gradient directions.
C
C     On Entry
C
C     n        integer
C              the order of the x
C
C     x        double precision x(n)
C              current point
C
C     m        integer
C     lambda   double precision lambda(m)
C     rho      double precision rho(m)
C     equatn   logical equatn(m)
C     linear   logical linear(m)
C              These five parameters are not used nor modified by 
C              GENCAN and they are passed as arguments to the user-
C              defined subroutines evalal and evalnal to compute the 
C              objective function and its gradient, respectively. 
C              Clearly, in an Augmented Lagrangian context, if GENCAN is 
C              being used to solve the bound-constrained subproblems, m 
C              would be the number of constraints, lambda the Lagrange 
C              multipliers approximation and rho the penalty parameters.
C              equatn is logical array that, for each constraint, 
C              indicates whether the constraint is an equality constraint
C              (.true.) or an inequality constraint (.false.)
C
C     f        double precision
C              function value at the current point
C
C     g        double precision g(n)
C              gradient vector at the current point
C
C     l        double precision l(n)
C              lower bounds
C
C     u        double precision u(n)
C              upper bounds
C
C     lamspg   double precision
C              spectral steplength
C
C     etaint   double precision
C              constant for the interpolation. See the description of
C              sigma1 and sigma2 above. Sometimes we take as a new 
C              trial step the previous one divided by etaint
C
C              RECOMMENDED: etaint = 2.0
C
C     mininterp integer
C              constant for testing if, after having made at least 
C              mininterp interpolations, the steplength is so small. In 
C              that case failure of the line search is declared (may be 
C              the direction is not a descent direction due to an error 
C              in the gradient calculations) 
C
C              RECOMMENDED: mininterp = 4
C
C     fmin     double precision
C              functional value for the stopping criterion f <= fmin
C
C     maxfc    integer
C              maximum number of functional evaluations
C
C     iprint   integer
C              Commands printing. Nothing is printed if iprint is 
C              smaller than 2. If iprint is greater than or equal to 
C              2, GENCAN iterations information is printed. If iprint 
C              is greater than or equal to 3, line searches and 
C              Conjugate Gradients information is printed.
C
C              RECOMMENDED: iprint = 2
C
C              CONSTRAINTS: allowed values are just 2 or 3.
C
C     xtrial,d double precision xtrial(n),d(n)
C              n-dimension working vectors of doubles
C
C     gamma    double precision
C              constant for the Armijo criterion
C              f(x + alpha d) <= f(x) + gamma * alpha * <\nabla f(x),d>
C
C              RECOMMENDED: gamma = 10^{-4}
C
C     sigma1   double precision
C     sigma2   double precision
C              constant for the safeguarded interpolation
C              if alpha_new \notin [sigma1, sigma*alpha] then we take
C              alpha_new = alpha / etaint
C
C              RECOMMENDED: sigma1 = 0.1 and sigma2 = 0.9
C
C     On Return
C
C     x        double precision
C              final estimation of the solution
C
C     f        double precision
C              functional value at the final estimation 
C
C     fcnt     integer
C              number of functional evaluations used in the line search   
C
C     inform   integer
C              This output parameter tells what happened in this 
C              subroutine, according to the following conventions:
C
C              0 = convergence with an Armijo-like criterion
C                  (f(xnew) <= f(x) + gamma * alpha * <g,d>);
C
C              4 = the algorithm stopped because the functional value
C                  is smaller than fmin;
C
C              6 = too small step in the line search. After having made 
C                  at least mininterp interpolations, the steplength 
C                  becames small. ''small steplength'' means that we are 
C                  at point x with direction d and step alpha, such that
C
C                  | alpha * d(i) | <= macheps * max( | x(i) |, 1 ) 
C
C                  for all i.
C 
C                  In that case failure of the line search is declared 
C                  (maybe the direction is not a descent direction due
C                  to an error in the gradient calculations). Use
C                  mininterp > maxfc to inhibit this criterion;
C
C              8 = it was achieved the maximum allowed number of
C                  function evaluations (maxfc);
C
C            < 0 = error in evalf subroutine.

C     LOCAL SCALARS
      logical samep
      integer i,interp
      double precision alpha,atmp,ftrial,gtd

C     Print presentation information

      if ( iprint .ge. 3 ) then
          write(*, 980) lamspg
          write(10,980) lamspg
      end if

C     Initialization

      interp = 0

C     Compute first trial point, spectral projected gradient direction, 
C     and directional derivative <g,d>.

      alpha = 1.0d0

      gtd = 0.0d0
      do i = 1,n
          xtrial(i) = min( u(i), max( l(i), x(i) - lamspg * g(i) ) )
          d(i)      = xtrial(i) - x(i)
          gtd       = gtd + g(i) * d(i)
      end do

      call evalal(n,xtrial,m,lambda,rho,equatn,linear,ftrial,inform)
      fcnt = fcnt + 1

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 3 ) then
              write(*, 1000) inform
              write(10,1000) inform
          end if

          return

      end if

C     Print information of the first trial

      if ( iprint .ge. 3 ) then
          write(*, 999) alpha,ftrial,fcnt
          write(10,999) alpha,ftrial,fcnt
      end if

C     Main loop

 100  continue

C     Test Armijo stopping criterion

      if ( ftrial .le. f + gamma * alpha * gtd ) then

          f = ftrial

          do i = 1,n
              x(i) = xtrial(i)
          end do

          inform = 0

          if ( iprint .ge. 3 ) then
              write(*, 990) inform
              write(10,990) inform
          end if

          go to 500

      end if

C     Test whether f is very small

      if ( ftrial .le. fmin ) then

          f = ftrial

          do i = 1,n
              x(i) = xtrial(i)
          end do

          inform = 4

          if ( iprint .ge. 3 ) then
              write(*, 994) inform
              write(10,994) inform
          end if

          go to 500

      end if

C     Test whether the number of functional evaluations is exhausted

      if ( fcnt .ge. maxfc ) then

          if ( ftrial .lt. f ) then

              f = ftrial

              do i = 1,n
                  x(i) = xtrial(i)
              end do

          end if

          inform = 8

          if ( iprint .ge. 3 ) then
              write(*, 998) inform
              write(10,998) inform
          end if

          go to 500

      end if

C     Compute new step (safeguarded quadratic interpolation)

      interp = interp + 1

      if ( alpha .lt. sigma1 ) then
          alpha = alpha / etaint      

      else
          atmp = ( - gtd * alpha ** 2 ) / 
     +           ( 2.0d0 * ( ftrial - f - alpha * gtd ) )

          if ( atmp .lt. sigma1 .or. atmp .gt. sigma2 * alpha ) then
              alpha = alpha / etaint

          else
              alpha = atmp
          end if
      end if

C     Compute new trial point

      do i = 1,n
          xtrial(i) = x(i) + alpha * d(i)
      end do

      call evalal(n,xtrial,m,lambda,rho,equatn,linear,ftrial,inform)
      fcnt = fcnt + 1

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 3 ) then
              write(*, 1000) inform
              write(10,1000) inform
          end if

          return

      end if

C     Print information of the current trial

      if ( iprint .ge. 3 ) then
          write(*, 999) alpha,ftrial,fcnt
          write(10,999) alpha,ftrial,fcnt
      end if

C     Test whether at least mininterp interpolations were made and two 
C     consecutive iterates are close enough

      samep = .true.
      do i = 1,n
         if ( abs( alpha * d(i) ) .gt. 
     +   macheps * max( abs( x(i) ), 1.0d0 ) ) then
             samep = .false.
         end if
      end do

      if ( interp .ge. mininterp .and. samep ) then

          if ( ftrial .lt. f ) then

              f = ftrial

              do i = 1,n
                  x(i) = xtrial(i)
              end do

          end if

          inform = 6

          if ( iprint .ge. 3 ) then
              write(*, 996) inform
              write(10,996) inform
          end if
  
          go to 500

      end if

C     Iterate

      go to 100

C     Return

 500  continue

      return

C     Non-executable statements

 980  format(/,6x,'SPG (spectral steplength ',1PD11.4,')',/,/,
     *         6x,'SPG Line search')
 999  format(6x,'Alpha= ',1PD11.4,' F= ',1PD11.4,' FE= ',I5)
 990  format(6x,'Flag of SPG Line search = ',I3,
     *          ' (Armijo-like criterion satisfied)')
 994  format(6x,'Flag of SPG Line search = ',I3,
     *          ' (Small functional value, smaller than ',/,
     *       6X,'parameter fmin)')
 996  format(6x,'Flag of SPG Line search = ',I3,
     *          ' (Too small step in the interpolation)')
 998  format(6x,'Flag of SPG Line search = ',I3,
     *          ' (Too many functional evaluations)')
 1000 format(6x,'Flag of SPG Line search = ',I3,' Fatal Error')

      end

C     ******************************************************************
C     ******************************************************************

      subroutine cgm(nind,ind,n,x,m,lambda,rho,equatn,linear,g,delta,l,
     +u,eps,epsnqmp,maxitnqmp,maxit,gtype,hptype,trtype,precond,samefa,
     +aptype,s,y,ssupn,seucn,yeucn,sts,sty,lspgmi,lspgma,iprint,ncomp,d,
     +iter,rbdtype,rbdind,inform,p,hp,r,z,pdiag,psmdy,hds,wdn1,wdn2,
     +theta,macheps,bignum)

      implicit none

C     SCALAR ARGUMENTS
      logical samefa
      character * 6 aptype,precond
      integer gtype,hptype,inform,iprint,iter,m,maxit,maxitnqmp,n,
     +        ncomp,nind,trtype,rbdind,rbdtype
      double precision bignum,delta,eps,epsnqmp,lspgma,lspgmi,macheps,
     +        seucn,ssupn,sts,sty,theta,yeucn

C     ARRAY ARGUMENTS
      integer ind(nind)
      logical equatn(m),linear(m)
      double precision d(n),g(n),hds(n),hp(n),l(n),lambda(m),p(n),
     +        pdiag(n),psmdy(n),r(n),rho(m),s(n),u(n),wdn1(n),wdn2(n),
     +        x(n),y(n),z(n)

C     This subroutine implements the Conjugate Gradients method for 
C     minimizing the quadratic approximation q(d) of L(x,lambda,rho) 
C     at x
C
C     q(d) = 1/2 d^T H d + g^T d,
C
C     where H is an approximation of the Hessian matrix of the 
C     Augmented Lagrangian and g is its gradient vector,
C
C     subject to || d || <= delta and l <= x + d <= u.
C
C     In the constraint ''|| d || <= delta'', the norm will be the
C     Euclidian-norm if the input parameter trtype is equal to 0, and
C     it will be the sup-norm if trtype is equal to 1.
C
C     The method returns an approximation d of the solution such that
C 
C     (a) ||H d + g||_2 <= eps * ||g||_2, 
C
C     (b) ||d|| = delta or x + d is in the boundary of the box, or
C
C     (c) ( p such that p^t H p = 0 ) and ( d = - amax g if such p was 
C         found during the first CG iteration or the current point d 
C         of CG if such p was found in any other iteration ).
C
C     On Entry
C
C     nind     integer
C              number of free variables (this is thee dimension in 
C              which this subroutine will work)
C
C     ind      integer ind(n)
C              array which contains, in the first nind positions, the
C              identifiers of the free variables
C
C     n        integer
C              dimension of the full space
C
C     x        double precision x(n)
C              point at which function L is being approximated by the
C              quadratic model
C
C              The first nind positions of x contains the free variables 
C              x_ind(1), x_ind(2), ..., x_ind(nind).
C
C     m        integer
C     lambda   double precision lambda(m)
C     rho      double precision rho(m)
C     equatn   logical equatn(m)
C     linear   logical linear(m)
C              These five parameters are not used nor modified by 
C              GENCAN and they are passed as arguments to the user-
C              defined subroutines evalal and evalnal to compute the 
C              objective function and its gradient, respectively. 
C              Clearly, in an Augmented Lagrangian context, if GENCAN is 
C              being used to solve the bound-constrained subproblems, m 
C              would be the number of constraints, lambda the Lagrange 
C              multipliers approximation and rho the penalty parameters.
C              equatn is logical array that, for each constraint, 
C              indicates whether the constraint is an equality constraint
C              (.true.) or an inequality constraint (.false.). Finally,
C              linear is logical array that, for each constraint, 
C              indicates whether the constraint is a linear constraint
C              (.true.) or a nonlinear constraint (.false.)
C
C     g        double precision g(n)
C              linear coefficient of the quadratic function
C
C              This is \nabla L(x) and it also contains in the first 
C              nind positions the components g_ind(1), g_ind(2), ..., 
C              g_ind(nind).
C    
C              IMPORTANT: the linear algebra of this subroutine lies in 
C              a space of dimension nind. The value of the full 
C              dimension n, the non-free variables (which are at the end 
C              of array x) and its gradient components (which are at the 
C              and of array g) are, at this moment, being used to 
C              approximate the Hessian times vector products by 
C              incremental quotients.
C
C     delta    double precision
C              "trust region radius" ( ||d|| <= delta ) 
C    
C     l        double precision l(n)
C              lower bounds on x + d. Its components are ordered in the 
C              same way as x and g.
C
C     u        double precision u(n)
C              upper bounds on x + d. Its components are ordered in the 
C              same way as x, g and l.
C
C     eps      double precision
C              tolerance for the stopping criterion ||H d + g||_2 < eps 
C              * ||g||_2
C
C     epsnqmp  double precision
C              See below
C
C     maxitnqmp integer
C              This and the previous parameter are used for a stopping 
C              criterion of the conjugate gradient subalgorithm. If the 
C              progress in the quadratic model is less than or equal to 
C              a fraction of the best progress ( epsnqmp * bestprog ) 
C              during maxitnqmp consecutive iterations then CG stops by 
C              not enough progress of the quadratic model.
C
C              RECOMMENDED: epsnqmp = 1.0d-4, maxitnqmp = 5
C
C     maxit    integer
C              maximum number of iterations.
C
C              RECOMMENDED: maxit = nind
C
C     gtype    integer
C              type of gradient calculation
C              gtype = 0 means user suplied evalng subroutine,
C              gtype = 1 means central difference approximation.
C
C              RECOMMENDED: gtype = 0
C
C              (provided you have the evalg subroutine)
C
C     hptype   integer
C              Type of Hessian-vector product calculation. See the 
C              detailed explanation in the algencan parameters
C              description.
C
C              RECOMMENDED: hptype = 6
C
C     trtype   integer
C              type of "trust region"
C              trtype = 0 means Euclidian-norm trust-region
C              trtype = 1 means sup-norm trust-region
C
C              RECOMMENDED: trtype = 0
C
C     iprint   integer
C              Commands printing. Nothing is printed if iprint is 
C              smaller than 2. If iprint is greater than or equal to 
C              2, GENCAN iterations information is printed. If iprint 
C              is greater than or equal to 3, line searches and 
C              Conjugate Gradients information is printed.
C
C              RECOMMENDED: iprint = 2
C
C              CONSTRAINTS: allowed values are just 2 or 3.
C
C     ncomp    integer
C              This constant is just for printing. In a detailed 
C              printing option, ncomp component of some vectors will be 
C              printed
C
C              RECOMMENDED: ncomp = 5
C
C              CONSTRAINTS: ncomp >= 0
C
C     hp,xtmp,r,p double precision hp(n),xtmp(n),r(n),p(n)
C              n-dimensional working vectors of doubles
C
C     theta    double precision
C              constant for the angle condition, i.e., at iteration k we 
C              need a direction pk such that 
C
C              <gk,pk> <= - theta ||gk||_2 ||pk||_2, 
C
C              where gk is \nabla L(xk)
C
C              RECOMMENDED: theta = 10^{-6}
C
C     On Return
C
C     d        double precision d(n)
C              final estimation of the solution
C
C     iter     integer
C              number of Conjugate Gradient iterations performed
C
C     inform   integer
C              termination parameter:
C
C              0 = convergence with ||H d + g||_2 <= eps * ||g||_2;
C
C              1 = convergence to the boundary of ||d|| <= delta;
C
C              2 = convergence to the boundary of l <= x + d <= u;
C
C              3 = stopping with d = dk  such that <gk,dk> <= - theta 
C                  ||gk||_2 ||dk||_2 and <gk,d_{k+1}> > - theta 
C                  ||gk||_2 ||d_{k+1}||_2;
C
C              4 = not enough progress of the quadratic model during
C                  maxitnqmp iterations, i.e., during maxitnqmp 
C                  iterations | q - qprev | <= macheps * max( | q |, 1 )
C
C              6 = very similar consecutive iterates, for two 
C                  consecutive iterates x1 and x2 we have that
C
C                  | x2(i) - x1(i) | <= macheps * max ( | x1(i) |, 1 )
C
C                  for all i.
C
C              7 = stopping with p such that p^T H p = 0 and g^T p = 0;
C
C              8 = too many iterations;
C
C            < 0 = error in evalhalp subroutine.

C     LOCAL SCALARS
      character * 4 cgtype
      character * 5 rbdtypea
      character * 6 prectmp
      logical goth,gotp,negcur,restarted,samep
      integer i,itertmp,itnqmp,rbdposaind,rbdposatype
      double precision aa,alpha,amax,amax1,amax2,amax2x,bb,bestprog,
     +        beta,cc,currprog,dd,dnorm2,gnorm2,gtd,gtp,hlspg,hstds,
     +        norm2s,pnorm2,plspg,psmdyty,ptd,pthp,ptr,q,qprev,rnorm2,
     +        ztrprev,ztr,znorm2

C     ==================================================================
C     Initialization
C     ==================================================================

      restarted = .false.

 001  continue

      goth = .false.
      gotp = .false.

      gnorm2   = norm2s(nind,g)

      iter     =      0
      itnqmp   =      0
      qprev    = bignum
      bestprog =  0.0d0

      do i = 1,nind
          d(i) = 0.0d0
          r(i) =  g(i)
      end do

      q        =  0.0d0
      gtd      =  0.0d0
      dnorm2   =  0.0d0
      rnorm2   = gnorm2

      ztr      =  0.0d0

C     ==================================================================
C     Print initial information
C     ==================================================================

      if ( iprint .ge. 3 ) then

          if ( precond .eq. 'NONE' ) then
              cgtype = '    '
          else
              cgtype = 'PREC'
          end if

          write(*, 980) cgtype,maxit,eps
          if ( trtype .eq. 0 ) then
              write(*, 981) delta
          else if ( trtype .eq. 1 ) then
              write(*, 982) delta
          else
              write(*, 983)
          end if
          write(*, 984) iter,sqrt(rnorm2),sqrt(dnorm2),q

          write(10,980) cgtype,maxit,eps
          if ( trtype .eq. 0 ) then
              write(10,981) delta
          else if ( trtype .eq. 1 ) then
              write(10,982) delta
          else
              write(10,983)
          end if
          write(10,984) iter,sqrt(rnorm2),sqrt(dnorm2),q

      end if

C     ==================================================================
C     Main loop
C     ==================================================================

 100  continue

C     ==================================================================
C     Test stopping criteria
C     ==================================================================

C     if ||r||_2 = ||H d + g||_2 <= eps * ||g||_2 then stop

      if ( iter .ne. 0 .and. 
     +   ( ( rnorm2 .le. eps ** 2 * gnorm2 .and. iter .ge. 4 ) .or. 
     +     ( rnorm2 .le. 1.0d-16 ) ) ) then 

          inform = 0

          if ( iprint .ge. 3 ) then
              write(*, 990) inform
              write(10,990) inform
          end if
  
          go to 500

      end if

C     if the maximum number of iterations was achieved then stop

      if ( iter .ge. max(4, maxit) ) then

          inform = 8

          if ( iprint .ge. 3 ) then
              write(*, 998) inform
              write(10,998) inform
          end if
  
          go to 500

      end if

C     ==================================================================
C     Preconditioner
C     ==================================================================

      if ( precond .eq. 'NONE' ) then

          do i = 1,nind
              z(i) = r(i)
          end do

          ztrprev = ztr
          ztr     = rnorm2
          znorm2  = rnorm2

      else if ( precond .eq. 'QNCGNA' ) then

          call calcpz(nind,ind,n,r,m,lambda,rho,equatn,linear,s,y,ssupn,
     +    seucn,yeucn,sts,sty,lspgmi,lspgma,samefa,gotp,pdiag,plspg,
     +    psmdy,psmdyty,z)

          ztrprev = ztr

          ztr = 0.0d0
          do i = 1,nind
              ztr = ztr + z(i) * r(i)
          end do

          znorm2 = norm2s(nind,z)

      end if

C     ==================================================================
C     Compute direction
C     ==================================================================

      if ( iter .eq. 0 ) then

          do i = 1,nind
              p(i) = - z(i)
          end do

          ptr    = - ztr
          pnorm2 =   znorm2

      else

          beta = ztr / ztrprev

          do i = 1,nind
              p(i) = - z(i) + beta * p(i)
          end do

          if ( precond .eq. 'NONE' ) then

              pnorm2 = rnorm2 - 2.0d0 * beta * ( ptr + alpha * pthp ) 
     +               + beta ** 2 * pnorm2
              ptr = - rnorm2 + beta * ( ptr + alpha * pthp )

          else if ( precond .eq. 'QNCGNA' ) then

              ptr = 0.0d0
              pnorm2 = 0.0d0
              do i = 1,nind
                  ptr = ptr + p(i) * r(i)
                  pnorm2 = pnorm2 + p(i) ** 2
              end do

          end if

      end if

C     Force p to be a descent direction of q(d), i.e.,
C     <\nabla q(d), p> = <H d + g, p> = <r, p> \le 0.

      if ( ptr .gt. 0.0d0 ) then

          do i = 1,nind
              p(i) = - p(i)
          end do

          ptr = - ptr

      end if

C     ==================================================================
C     Compute p^T H p
C     ==================================================================

C     hp = H p

      call calchalp(nind,ind,x,p,g,n,x,m,lambda,rho,equatn,linear,s,y,
     +ssupn,seucn,yeucn,sts,sty,lspgmi,lspgma,samefa,gtype,hptype,
     +aptype,hp,wdn1,wdn2,macheps,inform,goth,hlspg,hds,hstds)

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 3 ) then
              write(*, 1000) inform
              write(10,1000) inform
          end if

          return
      end if

C     Compute p^T hp

      pthp = 0.0d0
      do i = 1,nind
          pthp = pthp + p(i) * hp(i)
      end do 

C     ==================================================================
C     Compute maximum steps
C     ==================================================================

C     amax1 is the value of alpha such that ||d + alpha * p||_2 or 
C     ||d + alpha * p||_\infty = delta 

      ptd = 0.0d0
      do i = 1,nind
          ptd = ptd + p(i) * d(i)
      end do

C     Euclidian-norm trust radius

      if ( trtype .eq. 0 ) then

          aa = pnorm2
          bb = 2.0d0 * ptd
          cc = dnorm2 - delta ** 2
          dd = sqrt( bb ** 2 - 4.0d0 * aa * cc )

          amax1 = ( - bb + dd ) / ( 2.0d0 * aa )

C     Sup-norm trust radius

      else if ( trtype .eq. 1 ) then

          amax1 = bignum

          do i = 1,nind
              if ( p(i) .gt. 0.0d0 ) then
                  amax1 = min( amax1,  (   delta - d(i) ) / p(i) )
              else if ( p(i) .lt. 0.0d0 ) then
                  amax1 = min( amax1,  ( - delta - d(i) ) / p(i) )
              end if
          end do

      end if

C     amax2 is the maximum values of alpha such that 
C     l <= x + d + alpha * p <= u

      amax2 = bignum
      rbdposatype = 0
      do i = 1,nind
          if ( p(i) .gt. 0.0d0 ) then
                  amax2x = ( u(i) - x(i) - d(i) ) / p(i)
                  if ( amax2x .lt. amax2 ) then
                      amax2       = amax2x
                      rbdposaind  = i
                      rbdposatype = 2
                  end if
          else if ( p(i) .lt. 0.0d0 ) then
                  amax2x = ( l(i) - x(i) - d(i) ) / p(i)
                  if ( amax2x .lt. amax2 ) then
                      amax2       = amax2x
                      rbdposaind  = i
                      rbdposatype = 1
                  end if
          end if
      end do

C     Compute amax as the minimum among amax1 and amax2 

      amax  = min( amax1 , amax2  )

C     ==================================================================
C     Compute the step
C     ==================================================================

      negcur = .false.

C     If p^T H p > 0 then take the conjugate gradients step

      if ( pthp .gt. 0.0d0 ) then

          alpha = min( amax, ztr / pthp )

C     Else, if we are at iteration zero then take the maximum 
C     positive step in the minus gradient direction

      else if ( iter .eq. 0 ) then

          alpha = amax

          negcur = .true.

C     Otherwise, stop at the current iterate

      else

          inform = 7

          if ( iprint .ge. 3 ) then
              write(*, 997) inform
              write(10,997) inform
          end if
  
          go to 500

      end if

C     ==================================================================
C     Test the angle condition
C     ==================================================================

      gtp = 0.0d0
      do i = 1,nind
          gtp = gtp + g(i) * p(i)
      end do

      gtd = gtd + alpha * gtp
      dnorm2 = dnorm2 + alpha ** 2 * pnorm2 + 2.0d0 * alpha * ptd


      if ( gtd .gt. 0.0d0 .or. 
     +gtd ** 2 .lt. theta ** 2 * gnorm2 * dnorm2 ) then

          if ( precond .ne. 'NONE' .and. iter .eq. 0 ) then

              if ( iprint .ge. 3 ) then
                  write(*, 986)
                  write(10,986)
              end if

              restarted = .true.
              itertmp   = iter
              prectmp   = precond
              precond   = 'NONE'
              go to 001

          end if

          inform = 3

          if ( iprint .ge. 3 ) then
              write(*, 993) inform
              write(10,993) inform
          end if

          go to 500

      end if

C     ==================================================================
C     Compute the quadratic model functional value at the new point
C     ==================================================================

      qprev = q

      q = q + 0.5d0 * alpha ** 2 * pthp + alpha * ptr

C     ==================================================================
C     Compute new d
C     ==================================================================

      do i = 1,nind
          d(i) = d(i) + alpha * p(i)
      end do

C     ==================================================================
C     Compute the residual r = H d + g
C     ==================================================================

      do i = 1,nind
          r(i) = r(i) + alpha * hp(i)
      end do

      rnorm2 = norm2s(nind,r)

C     ==================================================================
C     Increment number of iterations
C     ==================================================================

      iter = iter + 1

C     ==================================================================
C     Print information of this iteration
C     ==================================================================

      if ( iprint .ge. 3 ) then
          write(*, 984) iter,sqrt(rnorm2),sqrt(dnorm2),q
          write(10,984) iter,sqrt(rnorm2),sqrt(dnorm2),q
      end if

C     ==================================================================
C     Test other stopping criteria
C     ==================================================================

C     Boundary of the box constraints

      if ( alpha .eq. amax2 ) then

          rbdind  = rbdposaind
          rbdtype = rbdposatype

          if ( rbdtype .eq. 1 ) then
              rbdtypea = 'lower'
          else if (rbdtype.eq.2) then
              rbdtypea = 'upper'
          end if
 
          inform = 2

          if ( iprint .ge. 3 ) then
              if ( negcur ) then
                  write(*, 999)
                  write(10,999)
              end if

              if ( rbdtype .eq. 0 ) then
                  write(*, 987) inform
                  write(10,987) inform
              else
                  write(*, 992) inform,ind(rbdind),rbdtypea
                  write(10,992) inform,ind(rbdind),rbdtypea
              end if
          end if
  
          go to 500

      end if

C     Boundary of the "trust region"

      if ( alpha .eq. amax1 ) then

          inform = 1

          if ( iprint .ge. 3 ) then
              if ( negcur ) then
                  write(*, 999)
                  write(10,999)
              end if

              write(*, 991) inform
              write(10,991) inform
          end if
  
          go to 500

      end if

C     Two consecutive iterates are too much close

      samep = .true.
      do i = 1,nind
         if ( abs( alpha * p(i) ) .gt. 
     +   macheps * max( abs( d(i) ) , 1.0d0 ) ) then 
              samep = .false.
          end if
      end do

      if ( samep ) then

          inform = 6

          if ( iprint .ge. 3 ) then
              write(*, 996) inform
              write(10,996) inform
          end if
  
          go to 500

      end if

C     Many iterations without good progress of the quadratic model 

      currprog = qprev - q
      bestprog = max( currprog, bestprog )

      if ( currprog .le. epsnqmp * bestprog ) then

          itnqmp = itnqmp + 1

          if ( itnqmp .ge. maxitnqmp ) then
              inform = 4

              if ( iprint .ge. 3 ) then
                  write(*, 994) inform,itnqmp,epsnqmp,bestprog
                  write(10,994) inform,itnqmp,epsnqmp,bestprog
              end if

              go to 500
          endif

      else
          itnqmp = 0
      endif

C     ==================================================================
C     Iterate
C     ==================================================================

      go to 100

C     ==================================================================
C     End of main loop
C     ==================================================================

C     ==================================================================
C     Return
C     ==================================================================

 500  continue

C     Print final information

      if ( iprint .ge. 3 ) then
          write(*, 985) min0(nind,ncomp),(d(i),i=1,min0(nind,ncomp))
          write(10,985) min0(nind,ncomp),(d(i),i=1,min0(nind,ncomp))
      end if

      if ( restarted ) then
          iter = iter + itertmp
          precond = prectmp
      end if

      return

C     Non-executable statements

 980  format(/,6x,'Conjugate gradients ',A4,' (maxit= ',I7,' acc= ',
     *1PD11.4,')')
 981  format(6x,'Using Euclidian trust region (delta= ',1PD11.4,
     *')')
 982  format(6x,'Using sup-norm trust region (delta= ',1PD11.4,')')
 983  format(6x,'Unknown trust-region type')
 984  format(6x,'CG iter= ',I5,' rnorm: ',1PD11.4,' dnorm= ',1PD11.4,
     *' q= ',1PD11.4)
 985  format(/,6x,'Truncated Newton direction (first ',I6, 
     *' components): ',/,1(6x,6(1PD11.4,1x)))
 986  format(6x,'The first CG-PREC iterate did not satisfy the angle ',
     *' condition. CG will be restarted without preconditioner)')
 990  format(6x,'Flag of CG = ',I3,' (Convergence with small residual)')
 991  format(6x,'Flag of CG = ',I3,
     *' (Convergence to the trust region boundary)')
 992  format(6x,'Flag of CG = ',I3,
     *' (Convergence to the boundary of the box constraints,',/,6x,
     *'taking step >= 1, variable ',I6,' will reaches its ',A5,
     *' bound)')
 987  format(6x,'Flag of CG = ',I3,
     *' (Taking a too large step, no variable will reach its bound)')
 993  format(6x,'Flag of CG = ',I3,
     *' (The next CG iterate will not satisfy the angle condition)')
 994  format(6x,'Flag of CG = ',I3,
     *' (Not enough progress in the quadratic model. This means',/,6x,
     *'that the progress of the last ',I7,' iterations was smaller ', 
     *'than ',/,6x,1PD11.4,' times the best progress (',1PD11.4,')')
 996  format(6x,'Flag of CG = ',I3,
     *' (Very near consecutive iterates)')
 997  format(6x,'Flag of CG= ',I3,
     *' (p such that p^T H p = 0 was found)')
 998  format(6x,'Flag of CG = ',I3,' (Too many GC iterations)')
 999  format(6x,'p such that p^T H p = 0 was found. ',
     *       'Maximum step was taken.')
 1000 format(6x,'Flag of CG = ',I3,' Fatal Error')

      end

C     *****************************************************************
C     *****************************************************************
      subroutine tnls(nind,ind,n,x,m,lambda,rho,equatn,linear,l,u,f,g,d,
     +amax,rbdtype,rbdind,etaint,etaext,mininterp,maxextrap,fmin,maxfc,
     +gtype,iprint,fcnt,gcnt,intcnt,exgcnt,exbcnt,inform,xplus,xtmp,
     +xbext,gamma,beta,sigma1,sigma2,macheps,bignum)

      implicit none

C     SCALAR ARGUMENTS
      integer exbcnt,exgcnt,fcnt,gcnt,gtype,inform,intcnt,iprint,m,
     +        maxextrap,maxfc,mininterp,n,nind,rbdind,rbdtype
      double precision amax,beta,bignum,etaext,etaint,f,fmin,gamma,
     +        macheps,sigma1,sigma2

C     ARRAY ARGUMENTS
      integer ind(nind)
      logical equatn(m),linear(m)
      double precision d(n),g(n),l(n),lambda(m),rho(m),u(n),x(n),
     +        xbext(n),xplus(n),xtmp(n)

C     This subroutine implements the line search used in the Truncated
C     Newton direction.
C
C     On Entry
C
C     nind     integer
C              number of free variables (this is thee dimension in 
C              which this subroutine will work)
C
C     ind      integer ind(n)
C              array which contains, in the first nind positions, the
C              identifiers of the free variables
C
C     n        integer
C              dimension of the full space
C
C     x        double precision x(n)
C              current point
C
C              The first nind positions of x contains the free variables 
C              x_ind(1), x_ind(2), ..., x_ind(nind).
C
C     m        integer
C     lambda   double precision lambda(m)
C     rho      double precision rho(m)
C     equatn   logical equatn(m)
C     linear   logical linear(m)
C              These five parameters are not used nor modified by 
C              GENCAN and they are passed as arguments to the user-
C              defined subroutines evalal and evalnal to compute the 
C              objective function and its gradient, respectively. 
C              Clearly, in an Augmented Lagrangian context, if GENCAN is 
C              being used to solve the bound-constrained subproblems, m 
C              would be the number of constraints, lambda the Lagrange 
C              multipliers approximation and rho the penalty parameters.
C              equatn is logical array that, for each constraint, 
C              indicates whether the constraint is an equality constraint
C              (.true.) or an inequality constraint (.false.). Finally,
C              linear is logical array that, for each constraint, 
C              indicates whether the constraint is a linear constraint
C              (.true.) or a nonlinear constraint (.false.)
C
C     l        double precision l(nind)
C              lower bounds on x. It components are ordered in the
C              same way as x and g.
C
C     u        double precision u(nind)
C              upper bounds on x. It components are ordered in the
C              same way as x, g and l.
C
C     f        double precision
C              functional value at x
C
C     g        double precision g(n)
C              gradient vector at x
C
C              It also contains in the first nind positions the 
C              components g_ind(1), g_ind(2), ..., g_ind(nind).
C    
C              IMPORTANT: the linear algebra of this subroutine lies in 
C              a space of dimension nind. The value of the full 
C              dimension n, the non-free variables (which are at the end 
C              of array x) and its gradient components (which are at the 
C              end of array g) are also used and updated any time the 
C              gradient is being computed.
C
C     d        double precision d(nind)
C              descent direction 
C    
C     amax     double precision
C
C     rbdtype  integer
C
C     rbdind   integer
C
C     etaint   double precision
C              constant for the interpolation. See the description of
C              sigma1 and sigma2 above. Sometimes we take as a new 
C              trial step the previous one divided by etaint
C
C              RECOMMENDED: etaint = 2.0
C
C     etaext   double precision
C              constant for the extrapolation
C              when extrapolating we try alpha_new = alpha * etaext
C
C              RECOMMENDED: etaext = 2.0
C
C     mininterp integer
C              constant for testing if, after having made at least 
C              mininterp interpolations, the steplength is so small. 
C              In that case failure of the line search is declared (may 
C              be the direction is not a descent direction due to an 
C              error in the gradient calculations) 
C
C              RECOMMENDED: mininterp = 4
C
C     maxextrap integer
C              constant to limit the number of extrapolations
C
C              RECOMMENDED: maxextrap = 1000 (a big number)
C
C     fmin     double precision
C              functional value for the stopping criteria f <= fmin
C
C     maxfc    integer
C              maximum number of functional evaluations
C
C     gtype    integer
C              type of gradient calculation
C              gtype = 0 means user suplied evalg subroutine,
C              gtype = 1 means central difference approximation.
C
C              RECOMMENDED: gtype = 0
C
C              (provided you have the evalg subroutine)
C
C     iprint   integer
C              Commands printing. Nothing is printed if iprint is 
C              smaller than 2. If iprint is greater than or equal to 
C              2, GENCAN iterations information is printed. If iprint 
C              is greater than or equal to 3, line searches and 
C              Conjugate Gradients information is printed.
C
C              RECOMMENDED: iprint = 2
C
C              CONSTRAINTS: allowed values are just 2 or 3.
C
C     xplus,xtmp,xbext double precision xplus(nind),xtmp(nind),xbext(nind) 
C              n-dimensional working vectors of doubles
C
C     gamma    double precision
C              constant for the Armijo criterion
C              f(x + alpha d) <= f(x) + gamma * alpha * <\nabla f(x),d>
C
C              RECOMMENDED: gamma = 10^{-4}
C
C     beta     double precision
C              constant for the beta condition <dk, g(xk + dk)>  <  beta 
C              * <dk,gk>. If (xk + dk) satisfies the Armijo condition 
C              but does not satisfy the beta condition then the point is 
C              accepted, but if it satisfied the Armijo condition and 
C              also satisfies the beta condition then we know that there 
C              is the possibility for a successful extrapolation
C
C              RECOMMENDED: beta = 0.5
C
C     sigma1   double precision
C     sigma2   double precision
C              constant for the safeguarded interpolation
C              if alpha_new \notin [sigma1, sigma*alpha] then we take
C              alpha_new = alpha / etaint
C
C              RECOMMENDED: sigma1 = 0.1 and sigma2 = 0.9
C
C     On Return
C
C     x        double precision x(n)
C              new current point
C
C     f        double precision
C              functional value at x
C
C     g        double precision g(n)
C              gradient vector at x
C
C     fcnt     integer
C              number of functional evaluations used in this line search
C
C     gcnt     integer
C              number of gradient evaluations used in this line search
C
C     intcnt   integer
C              number of interpolations
C
C     exgcnt   integer
C              number of good extrapolations
C
C     exbcnt   integer
C              number of bad extrapolations
C
C     inform   integer
C              This output parameter tells what happened in this 
C              subroutine, according to the following conventions:
C
C              0 = convergence with an Armijo-like criterion
C                  (f(xnew) <= f(x) + 1.0d-4 * alpha * <g,d>);
C
C              4 = the algorithm stopped because the functional value
C                  is very small (f <= fmin);
C
C              6 = so small step in the line search. After having made 
C                  at least mininterp interpolations, the steplength 
C                  becames small. ``small steplength'' means that we are 
C                  at point x with direction d and step alpha such that
C
C                  |alpha * d(i)| <= macheps * max ( |x(i)|, 1 )
C
C                  for all i. 
C 
C                  In that case failure of the line search is declared 
C                  (may be the direction is not a descent direction 
C                  due to an error in the gradient calculations). Use
C                  mininterp > maxfc for inhibit this criterion;
C
C              8 = it was achieved the maximum allowed number of
C                  function evaluations (maxfc);
C
C            < 0 = error in evalf or evalg subroutines.

C     LOCAL SCALARS
      logical samep
      integer extrap,i,interp
      double precision alpha,atmp,fbext,fplus,ftmp,gptd,gtd

C     ==================================================================
C     Initialization 
C     ==================================================================

C     ==================================================================
C     Compute directional derivative
C     ==================================================================

      gtd = 0.0d0
      do i = 1,nind
          gtd = gtd + g(i) * d(i)
      end do

C     ==================================================================
C     Compute first trial
C     ==================================================================

      alpha = min( 1.0d0, amax )

      do i = 1,nind
          xplus(i) = x(i) + alpha * d(i)
      end do

      if ( alpha .eq. amax ) then
          if ( rbdtype .eq. 1 ) then
              xplus(rbdind) = l(rbdind)
          else if (rbdtype.eq.2) then
              xplus(rbdind) = u(rbdind)
          end if
      end if

      call calcal(nind,ind,xplus,n,x,m,lambda,rho,equatn,linear,fplus,
     +inform)
      fcnt = fcnt + 1

      if ( inform .lt. 0 ) then
 
          if ( iprint .ge. 3 ) then
              write(*, 1000) inform
              write(10,1000) inform
          end if

          return

      end if

C     Print initial information

      if ( iprint .ge. 3 ) then
          write(*, 980) amax
          write(*, 999) alpha,fplus,fcnt

          write(10,980) amax
          write(10,999) alpha,fplus,fcnt
      end if

C     ==================================================================
C     Test Armijo and beta-condition and decide for accepting the trial 
C     point, interpolate or extrapolate.
C     ==================================================================

      if ( amax .gt. 1.0d0 ) then

C         x + d belongs to the interior of the feasible set
          if ( iprint .ge. 3 ) then
              write(*, *) '     x+d belongs to int of the feasible set'
              write(10,*) '     x+d belongs to int of the feasible set'
          end if

C         Verify Armijo

          if ( fplus .le. f + gamma * alpha * gtd ) then

C             Armijo condition holds  
              if ( iprint .ge. 3 ) then
                  write(*, *) '     Armijo condition holds' 
                  write(10,*) '     Armijo condition holds' 
              end if

              call calcnal(nind,ind,xplus,n,x,m,lambda,rho,equatn,
     +        linear,g,gtype,macheps,inform)
              gcnt = gcnt + 1

              if ( inform .lt. 0 ) then

                  if ( iprint .ge. 3 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return

              end if

              gptd = 0.0d0
              do i = 1,nind
                  gptd = gptd + g(i) * d(i)
              end do

C             Verify directional derivative (beta condition)

              if ( gptd .lt. beta * gtd ) then

C                 Extrapolate
                  if ( iprint .ge. 3 ) then
                      write(*, *)'     The beta-condition does not hold'
                      write(*, *)'     We will extrapolate'
                      write(10,*)'     The beta-condition does not hold'
                      write(10,*)'     We will extrapolate'
                  end if

C                 f and x before extrapolation
                  fbext = fplus

                  do i = 1,nind
                      xbext(i) = xplus(i)
                  end do

                  go to 100

              else

C                 Step = 1 was ok, finish the line search
                  if ( iprint .ge. 3 ) then
                      write(*, *) '     The beta condition is also true'
                      write(*, *) '     Line search is over'
                      write(10,*) '     The beta condition is also true'
                      write(10,*) '     Line search is over'
                  end if

                  f = fplus

                  do i = 1,nind
                      x(i) = xplus(i)
                  end do

                  inform = 0

                  if ( iprint .ge. 3 ) then
                      write(*, 990) inform
                      write(10,990) inform
                  end if

                  go to 500

              end if

          else 

C             Interpolate
              if ( iprint .ge. 3 ) then
                  write(*, *) '     Armijo does not hold'
                  write(*, *) '     We will interpolate'
                  write(10,*) '     Armijo does not hold'
                  write(10,*) '     We will interpolate'
              end if

              go to 200

          end if

      else

C         x + d does not belong to the feasible set (amax <= 1)
          if ( iprint .ge. 3 ) then
              write(*, *) '     x+d does not belong to box-interior'
              write(10,*) '     x+d does not belong to box-interior'
          end if

          if ( fplus .lt. f ) then

C             Extrapolate
              if ( iprint .ge. 3 ) then
                  write(*, *) '     f(x+d) < f(x)'
                  write(*, *) '     We will extrapolate'
                  write(10,*) '     f(x+d) < f(x)'
                  write(10,*) '     We will extrapolate'
              end if

C             f and x before extrapolation
              fbext = fplus

              do i = 1,nind
                  xbext(i) = xplus(i)
              end do

              go to 100

          else 

C             Interpolate
              if ( iprint .ge. 3 ) then
                  write(*, *) '     f(x+d) >= f(x)'
                  write(*, *) '     We will interpolate'
                  write(10,*) '     f(x+d) >= f(x)'
                  write(10,*) '     We will interpolate'
              end if

              go to 200

          end if

      end if


C     ==================================================================
C     Extrapolation
C     ==================================================================

 100  continue

      extrap = 0

C     Test f going to -inf

 120  if ( fplus .le. fmin ) then

C         Finish the extrapolation with the current point

          f = fplus

          do i = 1,nind
              x(i) = xplus(i)
          end do

          if ( extrap .ne. 0 .or. amax .le. 1.0d0 ) then

              call calcnal(nind,ind,x,n,x,m,lambda,rho,equatn,linear,g,
     +        gtype,macheps,inform)
              gcnt = gcnt + 1

              if ( inform .lt. 0 ) then
 
                  if ( iprint .ge. 3 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return

              end if

              if ( f .lt. fbext ) then
                  exgcnt = exgcnt + 1
              else
                  exbcnt = exbcnt + 1
              end if

          end if

          inform = 4

          if ( iprint .ge.3 ) then
              write(*, 994) inform
              write(10,994) inform
          end if

          go to 500

      end if

C     Test maximum number of functional evaluations

      if ( fcnt .ge. maxfc ) then

C         Finish the extrapolation with the current point

          f = fplus

          do i = 1,nind
              x(i) = xplus(i)
          end do

C         If extrap=0 and amax>1 the gradient was computed for testing 
C         the beta condition and it is not necessary to compute it again
          if ( extrap .ne. 0 .or. amax .le. 1.0d0 ) then

              call calcnal(nind,ind,x,n,x,m,lambda,rho,equatn,linear,g,
     +        gtype,macheps,inform)
              gcnt = gcnt + 1

              if ( inform .lt. 0 ) then

                  if ( iprint .ge. 3 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return

              end if

              if ( f .lt. fbext ) then
                  exgcnt = exgcnt + 1
              else
                  exbcnt = exbcnt + 1
              end if

          end if

          inform = 8

          if ( iprint .ge. 3 ) then
              write(*, 998) inform
              write(10,998) inform
          end if

          go to 500

      end if

C     Test if the maximum number of extrapolations was exceeded

      if ( extrap .ge. maxextrap ) then

C         Finish the extrapolation with the current point

          f = fplus

          do i = 1,nind
              x(i) = xplus(i)
          end do

C         If extrap=0 and amax>1 the gradient was computed for testing 
C         the beta condition and it is not necessary to compute it again
          if ( extrap .ne. 0 .or. amax .le. 1.0d0 ) then

              call calcnal(nind,ind,x,n,x,m,lambda,rho,equatn,linear,g,
     +        gtype,macheps,inform)
              gcnt = gcnt + 1

              if ( inform .lt. 0 ) then

                  if ( iprint .ge. 3 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return

              end if

              if ( f .lt. fbext ) then
                  exgcnt = exgcnt + 1
              else
                  exbcnt = exbcnt + 1
              end if

          end if

          inform = 7

          if ( iprint .ge. 3 ) then
              write(*, 997) inform
              write(10,997) inform
          end if

          go to 500

      end if

C     Chose new step 

      if ( alpha .lt. amax .and. etaext * alpha .gt. amax ) then
          atmp = amax
      else
          atmp = etaext * alpha
      end if

C     Compute new trial point

      do i = 1,nind
          xtmp(i) = x(i) + atmp * d(i)
      end do

      if ( atmp .eq. amax ) then
          if ( rbdtype .eq. 1 ) then
              xtmp(rbdind) = l(rbdind)
          else if ( rbdtype .eq. 2 ) then
              xtmp(rbdind) = u(rbdind)
          end if
      end if

C     Project

      if ( atmp .gt. amax ) then
          do i = 1,nind
              xtmp(i) = max( l(i), min( xtmp(i), u(i) ) )
          end do
      end if

C     Test if this is not the same point as the previous one.
C     This test is performed only when alpha > amax.

      if( alpha .gt. amax ) then

          samep = .true.
          do i = 1,nind
              if ( abs( xtmp(i) - xplus(i) ) .gt. 
     +        macheps * max( abs( xplus(i) ), 1.0d0 ) ) then
                  samep = .false.
              end if
          end do

          if ( samep ) then

C             Finish the extrapolation with the current point

              f = fplus

              do i = 1,nind
                  x(i) = xplus(i)
              end do

C             If extrap=0 and amax>1 the gradient was computed for 
C             testing the beta condition and it is not necessary to 
C             compute it again
              if ( extrap .ne. 0 .or. amax .le. 1.0d0 ) then

                  call calcnal(nind,ind,x,n,x,m,lambda,rho,equatn,
     +            linear,g,gtype,macheps,inform)
                  gcnt = gcnt + 1

                  if ( inform .lt. 0 ) then

                      if ( iprint .ge. 3 ) then
                          write(*, 1000) inform
                          write(10,1000) inform
                      end if

                      return

                  end if

                  if ( f .lt. fbext ) then
                      exgcnt = exgcnt + 1
                  else
                      exbcnt = exbcnt + 1
                  end if

              end if

              inform = 0

              if ( iprint .ge. 3 ) then
                  write(*, 990) inform
                  write(10,990) inform
              end if

              go to 500

          end if

      end if

C     Evaluate function

      call calcal(nind,ind,xtmp,n,x,m,lambda,rho,equatn,linear,ftmp,
     +inform)
      fcnt = fcnt + 1

      if ( inform .lt. 0 ) then

C         if ( iprint .ge. 3 ) then
C             write(*, 1000) inform
C             write(10,1000) inform
C         end if

C         return

C         If the objective function is not well defined in an 
C         extrapolated point, we discard all the extrapolated points
C         and return to a safe region (where the point before
C         starting the extrapolations is)

          f = fbext

          do i = 1,nind
              x(i) = xbext(i)
          end do

C         If extrap=0 and amax>1 the gradient was computed for testing 
C         the beta condition and it is not necessary to compute it again
          if ( extrap .ne. 0 .or. amax .le. 1.0d0 ) then

              call csetpoint(nind,ind,x,n,x)

              call calcnal(nind,ind,x,n,x,m,lambda,rho,equatn,linear,g,
     +        gtype,macheps,inform)
              gcnt = gcnt + 1

              if ( inform .lt. 0 ) then

                  if ( iprint .ge. 3 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return

              end if

              exbcnt = exbcnt + 1

          end if

          inform = 0

          if ( iprint .ge. 3 ) then
              write(*, 1010) inform
              write(10,1010) inform
          end if

          go to 500

      end if

C     Print information of this iteration

      if ( iprint .ge. 3 ) then
          write(*, 999) atmp,ftmp,fcnt
          write(10,999) atmp,ftmp,fcnt
      end if

C     If the functional value decreases then set the current point and 
C     continue the extrapolation

      if ( ftmp .lt. fplus ) then

          alpha = atmp

          fplus = ftmp

          do i = 1,nind
              xplus(i) = xtmp(i)
          end do

          extrap = extrap + 1

          go to 120

C     If the functional value does not decrease then discard the last 
C     trial and finish the extrapolation with the previous point

      else

          f = fplus

          do i = 1,nind
              x(i) = xplus(i)
          end do

C         If extrap=0 and amax>1 the gradient was computed for testing 
C         the beta condition and it is not necessary to compute it again
          if ( extrap .ne. 0 .or. amax .le. 1.0d0 ) then

              call csetpoint(nind,ind,x,n,x)

              call calcnal(nind,ind,x,n,x,m,lambda,rho,equatn,linear,g,
     +        gtype,macheps,inform)
              gcnt = gcnt + 1

              if ( inform .lt. 0 ) then

                  if ( iprint .ge. 3 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return

              end if

              if ( f .lt. fbext ) then
                  exgcnt = exgcnt + 1
              else
                  exbcnt = exbcnt + 1
              end if

          end if

          inform = 0

          if ( iprint .ge.3 ) then
              write(*, 990) inform
              write(10,990) inform
          end if

          go to 500

      end if
C     ==================================================================
C     End of extrapolation
C     ==================================================================

C     ==================================================================
C     Interpolation
C     ==================================================================

 200  continue

      intcnt = intcnt + 1

      interp = 0

 210  continue

C     Test f going to -inf

      if ( fplus .le. fmin ) then

C         Finish the interpolation with the current point

          f = fplus

          do i = 1,nind
              x(i) = xplus(i)
          end do

          call calcnal(nind,ind,x,n,x,m,lambda,rho,equatn,linear,g,
     +    gtype,macheps,inform)
          gcnt = gcnt + 1

          if ( inform .lt. 0 ) then

              if ( iprint .ge. 3 ) then
                  write(*, 1000) inform
                  write(10,1000) inform
              end if

              return

          end if

          inform = 4

          if ( iprint .ge. 3 ) then
              write(*, 994) inform
              write(10,994) inform
          end if

          go to 500

      end if

C     Test maximum number of functional evaluations

      if ( fcnt .ge. maxfc ) then

C         As this is an abrupt termination then the current point of the 
C         interpolation may be worst than the initial one

C         If the current point is better than the initial one then
C         finish the interpolation with the current point else discard
C         all we did inside this line search and finish with the initial
C         point

          if ( fplus .lt. f ) then

              f = fplus

              do i = 1,nind
                  x(i) = xplus(i)
              end do

              call calcnal(nind,ind,x,n,x,m,lambda,rho,equatn,linear,g,
     +        gtype,macheps,inform)
              gcnt = gcnt + 1

              if ( inform .lt. 0 ) then

                  if ( iprint .ge. 3 ) then
                      write(*, 1000) inform
                      write(10,1000) inform
                  end if

                  return

              end if

          end if

          inform = 8

          if ( iprint .ge. 3 ) then
              write(*, 998) inform
              write(10,998) inform
          end if

          go to 500

      end if

C     Test Armijo condition

      if ( fplus .le. f + gamma * alpha * gtd ) then

C         Finish the line search
      
          f = fplus

          do i = 1,nind
              x(i) = xplus(i)
          end do

          call calcnal(nind,ind,x,n,x,m,lambda,rho,equatn,linear,g,
     +    gtype,macheps,inform)
          gcnt = gcnt + 1

          if ( inform .lt. 0 ) then

              if ( iprint .ge. 3 ) then
                  write(*, 1000) inform
                  write(10,1000) inform
              end if

              return

          end if

          inform = 0

          if ( iprint .ge. 3 ) then
              write(*, 990) inform
              write(10,990) inform
          end if

          go to 500

      end if

C     Compute new step

      interp = interp + 1

      if ( alpha .lt. sigma1 ) then
          alpha = alpha / etaint      

      else
          atmp = ( - gtd * alpha **2 ) / 
     +           (2.0d0 * ( fplus - f - alpha * gtd ) )

          if ( atmp .lt. sigma1 .or. atmp .gt. sigma2 * alpha ) then
              alpha = alpha / etaint

          else
              alpha = atmp
          end if
      end if

C     Compute new trial point

      do i = 1,nind
          xplus(i) = x(i) + alpha * d(i)
      end do

      call calcal(nind,ind,xplus,n,x,m,lambda,rho,equatn,linear,fplus,
     +inform)
      fcnt = fcnt + 1

      if ( inform .lt. 0 ) then

          if ( iprint .ge. 3 ) then
              write(*, 1000) inform
              write(10,1000) inform
          end if

          return

      end if

C     Print information of this iteration

      if ( iprint .ge. 3 ) then
          write(*, 999) alpha,fplus,fcnt
          write(10,999) alpha,fplus,fcnt
      end if

C     Test whether at least mininterp interpolations were made and two 
C     consecutive iterates are much close

      samep = .true.
      do i = 1,nind
         if ( abs( alpha * d(i) ) .gt. 
     +   macheps * max( abs( x(i) ), 1.0d0 ) ) then
             samep = .false.
         end if
      end do

      if ( interp .ge. mininterp .and. samep ) then

C         As this is an abrupt termination then the current point of the 
C         interpolation may be worst than the initial one

C         If the current point is better than the initial one then
C         finish the interpolation with the current point else discard 
C         all we did inside this line search and finish with the initial 
C         point

C         if ( fplus .lt. f ) then

C             f = fplus

C             do i = 1,nind
C                 x(i) = xplus(i)
C             end do

C             call calcnal(nind,ind,x,n,x,m,lambda,rho,equatn,linear,g,
C    +        gtype,macheps,inform)
C             gcnt = gcnt + 1

C             if ( inform .lt. 0 ) then 

C                 if ( iprint .ge. 3 ) then
C                     write(*, 1000) inform
C                     write(10,1000) inform
C                 end if

C                 return

C             end if

C         end if

C         The previous lines were commented because, as it is been used, 
C         this subroutine must return with the initial point in case of 
C         finding a very small interpolation step. From that initial 
C         point, something different will be tried.

          inform = 6

          if ( iprint .ge. 3 ) then
              write(*, 996) inform
              write(10,996) inform
          end if
  
          go to 500

      end if

C     Else, iterate

      go to 210
C     ==================================================================
C     End of interpolation
C     ==================================================================

 500  continue

C     ==================================================================
C     Return
C     ==================================================================

      return

C     Non-executable statements

 980  format(/,6X,'TN Line search (alphamax= ',1PD11.4,')')
 999  format(6X,'Alpha= ',1PD11.4,' F= ',1PD11.4,' FE= ',I5)
 990  format(6X,'Flag of TN Line search= ',I3,
     +          ' (Armijo-like criterion satisfied)')
 994  format(6X,'Flag of TN Line search= ',I3,
     +          ' (Small functional value, smaller than ',/,
     +       6X,'parameter fmin)')
 996  format(6X,'Flag of TN Line search= ',I3,
     +          ' (Too small step in the interpolation)')
 997  format(6X,'Flag of TN Line search= ',I3,
     +          ' (Too many extrapolations)')
 998  format(6X,'Flag of TN Line search= ',I3,
     +          ' (Too many functional evaluations)')
 1000 format(6X,'Flag of TN Line search = ',I3,' Fatal Error')
 1010 format(6X,'Flag of TN Line search= ',I3,
     +          ' (Fatal Error in an extrapolated point)')

      end

C     ******************************************************************
C     ******************************************************************

      subroutine calcal(nind,ind,x,n,xc,m,lambda,rho,equatn,linear,f,
     +inform)

      implicit none

C     SCALAR ARGUMENTS
      integer nind,n,m,inform
      double precision f

C     ARRAY ARGUMENTS
      integer ind(nind)
      logical equatn(m),linear(m)
      double precision x(n),xc(n),lambda(m),rho(m)

C     This subroutines computes the objective function. 
C
C     It is called from the reduced space (dimension nind), expands the
C     point x where the function will be evaluated and call the 
C     subroutine evalf to compute the objective function Finally, 
C     shrinks vector x to the reduced space. 
C
C     About subroutines named calc[something]. The subroutines whos 
C     names start with ``calc'' work in (are called from) the reduced 
C     space. Their tasks are (i) expand the arguments to the full space, 
C     (ii) call the corresponding ``eval'' subroutine (which works in 
C     the full space), and (iii) shrink the parameters again and also 
C     shrink a possible output of the ``eval'' subroutine. Subroutines
C     of this type are: calcal, calcnal and calchd. 
C     The corresponding subroutines in the full space are the user 
C     defined subroutines evalf, evalg and evalhd.

C     LOCAL SCALARS
      integer i

C     Complete x

      do i = nind + 1,n
          x(i) = xc(i)
      end do

C     Expand x to the full space

      call expand(nind,ind,n,x)

C     Compute f calling the user supplied subroutine evalf

      call evalal(n,x,m,lambda,rho,equatn,linear,f,inform)

C     Shrink x to the reduced space

      call shrink(nind,ind,n,x)

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine calcnal(nind,ind,x,n,xc,m,lambda,rho,equatn,linear,g,
     +gtype,macheps,inform)

      implicit none

C     SCALAR ARGUMENTS
      integer gtype,inform,m,n,nind
      double precision macheps

C     ARRAY ARGUMENTS
      integer ind(nind)
      logical equatn(m),linear(m)
      double precision x(n),xc(n),lambda(m),rho(m),g(n)

C     This subroutine computes the gradient vector g of the objective 
C     function. 
C
C     It is called from the reduced space (dimension nind), expands the
C     point x where the gradient will be evaluated and calls the user 
C     supplied subroutine evalg to compute the gradient vector. Finally, 
C     shrinks vectors x and g to the reduced space. 
C
C     About subroutines named calc[something]. The subroutines whos 
C     names start with ``calc'' work in (are called from) the reduced 
C     space. Their tasks are (i) expand the arguments to the full space, 
C     (ii) call the corresponding ``eval'' subroutine (which works in 
C     the full space), and (iii) shrink the parameters again and also 
C     shrink a possible output of the ``eval'' subroutine. Subroutines
C     of this type are: calcal, calcnal and calchd
C     The corresponding subroutines in the full space are the user 
C     defined subroutines evalf, evalg and evalhd.

C     LOCAL SCALARS
      integer i

C     Complete x

      do i = nind + 1,n
          x(i) = xc(i)
      end do

C     Expand x to the full space

      call expand(nind,ind,n,x)

C     Compute the gradient vector calling the user supplied subroutine 
C     evalg

      call evalnal(n,x,m,lambda,rho,equatn,linear,g,gtype,macheps,
     +inform)

C     Shrink x and g to the reduced space

      call shrink(nind,ind,n,x)
      call shrink(nind,ind,n,g)

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine calcpz(nind,ind,n,r,m,lambda,rho,equatn,linear,s,y,
     +ssupn,seucn,yeucn,sts,sty,lspgmi,lspgma,samefa,gotp,pdiag,plspg,
     +psmdy,psmdyty,z)

      implicit none

C     SCALAR ARGUMENTS
      logical gotp,samefa
      integer m,n,nind
      double precision lspgma,lspgmi,plspg,psmdyty,seucn,ssupn,sts,sty,
     +        yeucn

C     ARRAY ARGUMENTS
      logical equatn(m),linear(m)
      integer ind(nind)
      double precision lambda(m),r(n),rho(m),s(n),pdiag(n),psmdy(n),
     +        y(n),z(n)

C     LOCAL SCALARS
      integer i

C     Complete r with zeroes

      do i = nind + 1,n
          r(i) = 0.0d0
      end do

C     Expand r to the full space

      call expand(nind,ind,n,r)

C     Solve P z = r

      call applyp(n,r,m,lambda,rho,equatn,linear,s,y,ssupn,seucn,yeucn,
     +sts,sty,lspgmi,lspgma,samefa,gotp,pdiag,plspg,psmdy,psmdyty,z)

C     Shrink r and z to the reduced space

      call shrink(nind,ind,n,r)
      call shrink(nind,ind,n,z)

      end

C     ******************************************************************
C     ******************************************************************

      subroutine calchalp(nind,ind,x,p,g,n,xc,m,lambda,rho,equatn,
     +linear,s,y,ssupn,seucn,yeucn,sts,sty,lspgmi,lspgma,samefa,gtype,
     +hptype,aptype,hp,wdn1,wdn2,macheps,inform,goth,hlspg,hds,hstds)

      implicit none

C     This subroutine computes the product Hessian times vector p. As it
C     is called from the reduced space, it expands vectors x and p,  
C     calls subroutine evalhalp to compute the Hessian times vector p 
C     product, and shrinks vectors x, p and hp. 

C     SCALAR ARGUMENTS
      logical goth,samefa
      character * 6 aptype
      integer gtype,hptype,inform,m,n,nind
      double precision hlspg,hstds,lspgma,lspgmi,macheps,seucn,ssupn,
     +        sts,sty,yeucn

C     ARRAY ARGUMENTS
      integer ind(nind)
      logical equatn(m),linear(m)
      double precision g(n),hds(n),hp(n),lambda(m),p(n),rho(m),wdn1(n),
     +        wdn2(n),x(n),xc(n),s(n),y(n)

C     LOCAL SCALARS
      integer i

C     Complete p with zeroes

      do i = nind + 1,n
          p(i) = 0.0d0
      end do

C     Complete x

      do i = nind + 1,n
          x(i) = xc(i)
      end do

C     Expand x and p to the full space

      call expand(nind,ind,n,x)
      call expand(nind,ind,n,p)
      call expand(nind,ind,n,g)

C     Compute the Hessian-vector product

      call evalhalp(n,x,p,g,m,lambda,rho,equatn,linear,s,y,ssupn,seucn,
     +yeucn,sts,sty,lspgmi,lspgma,samefa,gtype,hptype,aptype,hp,wdn1,
     +wdn2,macheps,inform,goth,hlspg,hds,hstds)

C     Shrink x, p and hp to the reduced space

      call shrink(nind,ind,n,x)
      call shrink(nind,ind,n,p)
      call shrink(nind,ind,n,g)
      call shrink(nind,ind,n,hp)
      
      end

C     ******************************************************************
C     ******************************************************************

      subroutine shrink(nind,ind,n,v)

      implicit none

C     This subroutine shrinks vector v from the full dimension space 
C     (dimension n) to the reduced space (dimension nind).

C     SCALAR ARGUMENTS
      integer n,nind

C     ARRAY ARGUMENTS
      integer ind(nind)
      double precision v(n)

C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     NIND integer: Dimension of the reduced space.
C     ------------
C
C     IND integer ind(nind)
C     ---------------------
C
C     Components ind(1)-th, ..., ind(nind)-th are the components that 
C     belong to the reduced space.
C
C     N integer: Dimension of the full space.
C     ---------
C
C     V double precision v(n): Vector to be shrinked.
C     -----------------------
C
C     On Return:
C     ==========
C
C     V double precision v(n): Shrinked vector.
C     -----------------------

C     LOCAL SCALARS
      integer i,indi
      double precision tmp

      do i = 1,nind
           indi = ind(i)
           if ( i .ne. indi ) then
               tmp     = v(indi)
               v(indi) = v(i)
               v(i)    = tmp
          end if
      end do

      return

      end

C     ******************************************************************
C     ******************************************************************

      subroutine expand(nind,ind,n,v)

      implicit none

C     This subroutine expands vector v from the reduced space 
C     (dimension nind) to the full space (dimension n).

C     SCALAR ARGUMENTS
      integer n, nind

C     ARRAY ARGUMENTS
      integer ind(nind)
      double precision v(n)

C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     NIND integer: Dimension of the reduced space.
C     ------------
C
C     IND integer ind(nind)
C     ---------------------
C
C     Components ind(1)-th, ..., ind(nind)-th are the components that 
C     belong to the reduced space.
C
C     N integer: Dimension of the full space.
C     ---------
C
C     V double precision v(n): Vector to be expanded.
C     -----------------------
C
C     On Return:
C     ==========
C
C     V double precision v(n): Expanded vector.
C     -----------------------

C     LOCAL SCALARS
      integer i,indi
      double precision tmp

      do i = nind,1,- 1
          indi = ind(i)
          if ( i .ne. indi ) then
              tmp     = v(indi)
              v(indi) = v(i)
              v(i)    = tmp
          end if
      end do
     
      return

      end
 
C     ******************************************************************
C     ******************************************************************

      subroutine csetpoint(nind,ind,x,n,xc)

      implicit none

C     SCALAR ARGUMENTS
      integer nind,n

C     ARRAY ARGUMENTS
      integer ind(nind)
      double precision x(n),xc(n)

C     LOCAL SCALARS
      integer i

C     Complete x

      do i = nind + 1,n
          x(i) = xc(i)
      end do

C     Expand x to the full space

      call expand(nind,ind,n,x)

C     Call setpoint

      call setpoint(x)

C     Shrink x

      call shrink(nind,ind,n,x)

      return

      end

C     *****************************************************************
C     *****************************************************************

      double precision function norm2s(n,x)

C     This subroutine computes the squared Euclidian norm of an 
C     n-dimensional vector.

C     SCALAR ARGUMENTS
      integer n

C     ARRAY ARGUMENTS
      double precision x(n)

C     Parameters of the subroutine:
C     =============================
C
C     On Entry:
C     =========
C
C     N integer: Dimension.
C     ---------
C
C     X double precision x(n): Vector.
C     -----------------------
C
C     On Return:
C     ==========
C
C     The function return the squared Euclidian norm of the 
C     n-dimensional vector x.

      double precision hsldnrm2

      norm2s = hsldnrm2(n,x,1) ** 2

      return 

      end

C     ******************************************************************
C     ******************************************************************

      DOUBLE PRECISION FUNCTION HSLDNRM2(N,DX,INCX)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION CUTLO,CUTHI
      PARAMETER (CUTLO=8.232D-11,CUTHI=1.304D19)
      INTEGER INCX,N
      DOUBLE PRECISION DX(*)
      DOUBLE PRECISION HITEST,SUM,XMAX
      INTEGER I,J,NEXT,NN
      INTRINSIC DABS,DSQRT,FLOAT
      IF (N.GT.0) GO TO 10
      HSLDNRM2 = ZERO
      GO TO 300
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N*INCX
      I = 1
   20 GO TO NEXT
   30 IF (DABS(DX(I)).GT.CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
   50 IF (DX(I).EQ.ZERO) GO TO 200
      IF (DABS(DX(I)).GT.CUTLO) GO TO 85
      ASSIGN 70 TO NEXT
      GO TO 105
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM/DX(I))/DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
   70 IF (DABS(DX(I)).GT.CUTLO) GO TO 75
  110 IF (DABS(DX(I)).LE.XMAX) GO TO 115
      SUM = ONE + SUM* (XMAX/DX(I))**2
      XMAX = DABS(DX(I))
      GO TO 200
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
   75 SUM = (SUM*XMAX)*XMAX
   85 HITEST = CUTHI/DFLOAT(N)
      DO 95 J = I,NN,INCX
        IF (DABS(DX(J)).GE.HITEST) GO TO 100
   95 SUM = SUM + DX(J)**2
      HSLDNRM2 = DSQRT(SUM)
      GO TO 300
  200 CONTINUE
      I = I + INCX
      IF (I.LE.NN) GO TO 20
      HSLDNRM2 = XMAX*DSQRT(SUM)
  300 CONTINUE
      RETURN
      END

C     ******************************************************************
C     ******************************************************************

      double precision function drand(ix)

C     This is the random number generator of Schrage:
C
C     L. Schrage, A more portable Fortran random number generator, ACM
C     Transactions on Mathematical Software 5 (1979), 132-138.

      double precision ix

      double precision a,p,b15,b16,xhi,xalo,leftlo,fhi,k
      data a/16807.d0/,b15/32768.d0/,b16/65536.d0/,p/2147483647.d0/

      xhi= ix/b16
      xhi= xhi - dmod(xhi,1.d0)
      xalo= (ix-xhi*b16)*a
      leftlo= xalo/b16
      leftlo= leftlo - dmod(leftlo,1.d0)
      fhi= xhi*a + leftlo
      k= fhi/b15
      k= k - dmod(k,1.d0)
      ix= (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
      if (ix.lt.0) ix= ix + p
      drand= ix*4.656612875d-10

      return

      end

c#ifndef AMPL

C     ******************************************************************
C     ******************************************************************

      subroutine setpoint(x)

      implicit none

C     ARRAY ARGUMENTS
      double precision x(*)

      end

c#endif

c#ifndef BLAS

C     ******************************************************************
C     ******************************************************************

      double precision function d1mach(dum)

C     This function computes the machine epsilon. The BLAS subroutine
C     is prefered and the present one is used when BLAS is not available.

C     SCALAR ARGUMENTS
      integer dum

C     LOCAL SCALARS
      double precision macheps

      macheps = 1.0d0
 10   if ( 1.0d0 + macheps .ne. 1.0d0 ) then
          macheps = macheps / 2.0d0
          go to 10
      end if
      macheps = 2.0d0 * macheps

      d1mach = macheps

      return

      end

c#endif

C     ******************************************************************
C     ******************************************************************
C
C     Report of modifications of algencan.
C
C     May 22, 2007.
C
C     1) A few changes were made related to variable rbdtype (which 
C     indicates whether a lower or an upper bound is reached when a step
C     to the boundary of the box constraints is taken). In abnormal 
C     cases in which the gradient evaluation returns NaN, rbdtype
C     remained unset. Then, it may cause a segmentation fault. The bug
C     was fixed.
C
C     May 11, 2007.
C
C     1) A small bug was fixed in the CG subroutine. The Euclidian norm
C     of the new direction (used to test the angle condition) was being
C     computed in the wrong place.
C
C     February 27, 2007.
C
C     1) A messagge was added to distinguished the cases in which a
C     maximum step is taken by chance and when it is taken by purpouse
C     because a negative curvature direction was found.
C
C     2) Messagge "p^t H p = 0 and g^t p = 0" was wrong. The correct
C     one is just "p^t H p = 0". It was corrected. The description of the
C     Conjugate Gradient subrotuine was corrected acordingly.
C
C     January 30, 2007.
C
C     1) In option TRUEHE of evalhalp subroutine, we changed
C
C         COMPUTE THE HESSIANS OF THE CONSTRAINTS
C         do j = 1,m
C             if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then
C                 hsta(j) = ind + 1
C
C                 call evalhc(n,x,j,hlin(ind+1),hcol(ind+1),hval(ind+1),
C    +            hlen(j),flag)
C                 if ( flag .ne. 0 ) then
C                     inform = - 95
C                     return
C                 end if
C
C             end if
C
C             ind = ind + hlen(j) <-----------------------------
C         end do
C
C     by
C
C         COMPUTE THE HESSIANS OF THE CONSTRAINTS
C         do j = 1,m
C             if ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) then
C                 hsta(j) = ind + 1
C
C                 call evalhc(n,x,j,hlin(ind+1),hcol(ind+1),hval(ind+1),
C    +            hlen(j),flag)
C                 if ( flag .ne. 0 ) then
C                     inform = - 95
C                     return
C                 end if
C
C                 ind = ind + hlen(j) <-----------------------------
C             end if
C         end do
C
C     The previous version was wrong because hlen(j) is not computed
C     when ( equatn(j) .or. dpdc(j) .gt. 0.0d0 ) is false.
C
C     2) Two constants that affect the stopping criteria of CG were 
C     changed in order to achieve a more strict precision. Namely:
C     cgepsf change from 1.0d-05 to 1.0d-08 and a constant used to
C     stop CG by absolute error of the residual was also changed from
C     1.0d-12 to 1.0d-16.
C
C     3) Penalty parameter rho was moved from subrotuine inip to
C     subroutine param. In addition, a new parameter rhoauto was also
C     added to subroutine param. The user can choose rhoauto=.true.
C     to let the task to initialize the penalty parameters to algencan.
C     Otherwise, the user should choose rhoauto=.false. and set the 
C     values of vector rho. Both parameters were also added to
C     subroutine fparam.
C
C     4) The maximum number of iterations and functional evaluations
C     of gencan were changed from maxit=5000 and maxfc=10*maxit to
C     maxit=1000 and maxfc=2*maxit.
C
C     5) A subroutine (evalnl) to compute the gradient of the Lagrangian 
C     was coded; and the subroutine (evalnal) to compute the gradient
C     of the Augmented Lagrangian was modified to use evalnl. Finally,
C     at the initial point of algencan, the call to evalnal was replaced
C     by a call to evalnl. The intention of this modification was to
C     report correctly the sup-norm of the projected gradient of the
C     Lagrangian at the initial point (the norm of the gradient of the
C     AUGMENTED Lagrangian was being reported). This mistake was realised
C     when starting algencan from its final point. In this cases, the final
C     and the initial norms should coincide. 
C
C     6) A modification was made to do not modify the penalty parameters
C     after the first outer iteration. this is because the feasibility
C     at the initial point (against which we should compare to decide
C     whether to increase the penalty parameters or not) is a ver arbitray
C     quantity.
C
C     7) Parameter intype was added, to choose the type of algorithm that
C     will be used inside the faces of the active-set method used to solve
C     the bound-constrained subproblems.
C
C     8) nmax was replaced by nsmax in subroutines checkh and checkhc.
C
C     9) Some comments were modified.
C 
C     November 9, 2006.
C
C     1) Inform parameter was being left undefined within checkd 
C     subroutine (causing an abnormal termination) when the derivatives
C     checking were aborted. inform was initialited with zero within
C     checkd.
C
C     2) Even within checkd, the random point used to check derivatives
C     was being randomly selected within a reduced box that might be
C     empty. That was also corrected.
C
C     3) Finally, within checkh, checkjac and checkhc, we now print the 
C     non-null elements and the overall error per column at the end.
C
C     4) For the specification file algencan.dat, character '#' will
C     also be considered (thogether with '*') as the begining of a
C     comment line.
C
C     5) The message saying that the specification file was not found
C     was rephrase to avoid the wrong impression that it could be an
C     error message.
C
C     October 11, 2006:
C
C     1) A new stopping criterion to test wether the penalty parameters
C     are too large was added.
C
C     2) A modification to update the penalty parameters just if the
C     current point does not satisfy the feasibility tolerance was
C     added.
C
C     March 27, 2006:
C
C     1) Just "#ifndef" is used with the following objective: if the 
C     compiler does not recognize precompiling options then everything
C     is included and at least the fortran stand-alone version can be
C     used.
C
C     March 22, 2006:
C
C     1) A few simple modifications where done:
C
C     1.a) To avoid a warning given by FTN95 compiler, the name of 
C     variable nint was changed to etaint and, to mantain the coherence
C     we also changed next by etaext. The problem is that there exist
C     an intrinsic Fortran functiona named nint.
C
C     1.b) The real vector dum(2) was, although unnecessary, initialized
C     with zero using a data statement, just to avoid a warning given
C     by ftnchek-3.0.4:
C
C     Variables used before set DUM used at line 121 file algencan.f; 
C     never set
C
C     2) Due to the modifications of March 21, 2006, a (forgotten) new 
C     ittext ("Used QN with PURE INCREMENTAL QUOTIENTS") was defined.
C
C     3) Moreover, a forgotten statement that fix hptype=9 whenever
C     the user select gtype=1 was deleted from the begining of GENCAN
C     subroutine. 
C
C     March 21, 2006:
C
C     1) Many of the last updates of algencan related to the Hessian-
C     vector product and the computation of preconditioners were 
C     implemented without considering the possibility of being used in 
C     combination wit finite differences approximations of first 
C     derivatives. Due to a suggestion of Prof. N. Shamsundar 
C     (University of Houston) we completely rearranged the subroutines 
C     to compute first-derivatives approximations. In this rearrengment:
C
C     1.a) Subroutine calcnaldiff, evalnaldiff and calchalpdiff were 
C          deleted. 
C
C     1.b) Subroutines evalgdiff, evaljcadiff and evalgjcdiff were coded.
C
C     1.c) Parameter m was added to evalgjc.
C
C     February 7, 2006:
C
C     1) Several suggestions given by Prof. N. Shamsundar (University of 
C     Houston), to whom we are very grateful, were incorporated. They are 
C     described in the following two e-mails:
C
C     e-mail 1:
C
C     1. If Algencan.dat is used, and the user types in, for example, 
C     1D-4 for one of the tolerances, it gets read as 1D-12 because of the 
C     missing decimal point and the rules of Fortran, and the program 
C     fails to find a solution. The fix is to change the Format 3000 in 
C     subroutine FPARAM to something such as BN,F24.0 instead of 1P,D24.8. 
C     Another solution would be to do a list directed read, but this is 
C     not required to work for internal reads by the rules of Fortran-77.
C
C     2. Emphasize that all expressions in subroutine EVALF, especially 
C     constants, should be double precision. (IMPLICIT NONE in the user 
C     supplied routines would be beneficial.) Typing in 0.2+x(1)*x(3), for 
C     example, may cause checking of derivatives to fail; the perturbation 
C     of x is in the lower half of the significant digits, since you use 
C     (quite properly so) sqrt(macheps). If the expression for the 
C     analytical derivative is dominated by real*4 constants, numerical 
C     derivatives come out as zero!
C
C     3. ANALITICAL -> ANALYTICAL (spelling) in dictionary and description 
C     of parameters (at least in English!).
C
C     e-mail 2:
C
C     The suggestion is based on well known facts concerning the 
C     precision of numerically computed derivatives, which you probably 
C     know but I will repeat. I assume that the objective function and 
C     constraint functions are calculated to full machine precision.
C     
C     The error in the numerical derivative of f(x), computed using a 
C     one sided difference, is 2 eps_mach f(x)/h + h/2 f'(x). This is 
C     minimized by using a step size h ~ sqrt(eps_mach); this is what 
C     you have used throughout your program. The corresponding error is 
C     also ~ sqrt(eps_mach).     
C     
C     However, in EVALNALDIFF (and, less importantly, in CHECKG and 
C     CHECKJAC), you use a central difference. Here, the error in the 
C     f'(x) estimate is eps_mach f/h + h^2/6 f'''(x). This error is 
C     minimized by h ~ cube_root(eps_mach). The corresponding error ~ 
C     (eps_mach)^(2/3), which is orders of magnitude better (lower) 
C     than with the first order difference.
C
C     January 12, 2006:
C
C     1) A new option for the Hessian-of-the-Augmented-Lagrangian times
C     a vector product was added. The user can now provide a subroutine
C     to compute the product of the Hessian-of-the Lagrangian times an
C     arbitrary vector. (Instead of providing independet subroutines to
C     compute the Hessians of the objective function and the 
C     constraints.) This option was developed to be able to use second 
C     derivatives with AMPL and CUTEr.
C
C     2) Some input parameter explanations were re-written.
C
C     3) Subroutines calcf and calcg were renamed as calcal and calcnal.
C
C     4) Subroutine evalhp was renamed as evalhalp.
C
C     November 3, 2005:
C
C     1) Update of checkd subroutine. The Hessians checking was 
C     incorporated.
C
C     2) The Augmented Lagrangian hessian is automatically computed 
C     using the user provided Hessians. So evalh and evalhc subroutines
C     were incorporated.
C
C     3) Each gencan iteration discriminates in the output between TN
C     and QN iterations.
C
C     4) If the problem has just bound constraints, GENCAN is called
C     instead of ALGENCAN. This task is made by the new subroutine called
C     solver.
C
C     September 26, 2005.
C
C     1) Some modifications in the subroutine to check derivatives:
C     less precision in the write sentence; to different steps in the
C     finite differences calculation.
C
C     September 20, 2005.
C
C     1) For compatibility with the R interface, we change the name of 
C     the conjugate gradients subroutine from cg to cgm.
C
C     2) The specification fiel functionality was added.
C
C     3) The CPU time is measured by algencan and it is now an output
C     parameter.
C
C     4) The call to subroutine param, the opening of the output file
C     algencan.out and the creation of the output file solution.txt
C     were changed from the main program algencama to to subroutine
C     algencan. It was made to simplify the interfaces which, basically,
C     needs the main file to be recoded in a different language.
C
C     June 7th, 2005.
C
C     1) Parameter jcnnzmax was increased.
C
C     2) A break line was added in the message mentioned below.
C
C     April 29th, 2005.
C
C     1) Phrase "infeasibility improvement" was replaced by just 
C     "infeasibility" in the messages related to the increase of the
C     penalty parameters.
C
C     2) A small error that had been forced the method to do an extra 
C     outer iteration was corrected. The sup-norm of sigma that is 
C     being used as the stopping criterion related to feasibility and
C     complementarity of the inequality constraints must be computed
C     _after_ the update of the Lagrange multipliers, and it had been
C     being test before. Toghether with this modification, the 
C     secondary loop related to the subproblem precision necessary
C     to prove the limitation of the penalty parameters was also 
C     deleted. In fact, it was the secondary loop, that required to
C     continue solving the same subproblem (without updating the
C     Lagrange multipliers) that leaved us to compute the sup-norm of
C     sigma before updating the Lagrange multipliers. It was because
C     the test to verify if the subproblem was abandoned prematurely
C     used the sup-norm of sigma in the comparison.
C
C     April 25th, 2005.
C
C     1) After exaustive tests, we discover that using the new version 
C     of subroutine evalhd the method performs bad in a few problems.
C     With this hint, a bug was fixed.
C      
C     April 20th, 2005.
C
C     In the algencanma-tabline.out and algencan-tabline.out one-line
C     files, the number of active bounds at the final point was added.
C     This quantity is being computed in the main program algencanma
C     and in subroutine algencan just for informative purposes. 
C
C     April 18th, 2005.
C
C     1) evalal and evalnal subroutines were rewritten using different
C     CUTE subroutines.
C
C     2) The contribution of the linear constraints to the matrix-vector
C     product computed in evalhd subroutine was modified. It was changed
C     from incremental quotients to exact calculation. 
C
C     April 13th, 2005.
C
C     1) Vector linear was added to indicate which constraints are
C     linear constraints.
C
C     2) In the output files algencan-tabline.out and 
C     algencanma-tabline.out the number of active bounds at the final
C     point was added.
C
C     March 11th, 2005
C
C     1) equatn is not any more in a common block. It is now a parameter
C     of almost all the subroutines.
C 
C     March 10th, 2005.
C
C     1) Subroutines evalnal and evalhd were modified to avoid the
C     unnecessary calculation of the gradients of the constraints that
C     would appear multiplied by a null coefficient in the Augmented
C     Lagrangian gradient formula. 
C
C     2) Output files with the output in the screen and with the final 
C     estimation of the solution, the Lagrange multipliers and the 
C     penalty parameters are now being generated by algencan.
C
C     3) The initial penalty parameters are now being seted by the user
C     in the subroutine inip.
C
C     March 9th, 2005.
C
C     1) Subroutines evalc and evaljac were re-written to compute just a 
C     constraint and the sparse gradient of a unique constraint 
C     (instead of all the constraints and the Jacobian matrix), 
C     respectively.
C
C     2) Due to the modifications related above, subroutine evalgjc was 
C     deleted and subroutine evalfc was rewritten.
C
C     3) Three local arrays of algencan were replaced by three working
C     vectors passed as arguments of algencan. easyalgencan was also
C     modified to call algencan with its new prototype.
C
C     4) The position of lammin and lammax in the calling sequence of
C     algencan was altered.
C
C     5) Explanations of lammin and lammax were added to algencan. 
C
C     March 7th, 2005.
C
C     1) The reference to the submitted report was modified.
C
C     February 18th, 2005.
C
C     1) The maximum number of iterations and functional evaluations that
C     ALGENCAN gives to GENCAN to solve a subproblem were increased from
C     1000 both to 5000 iterations and 10 times the maximum number of
C     iterations for the maximum number of functional evaluations. This 
C     change was in fact needed to solve some packing problems.
C
C     2) The format statements of the messages related to infeasibility
C     improvements were rewritten. The final ALGENCAN statement was also
C     corrected (it was written larger instead of largest). 
C
C     3) An unused format statement related the the desired and the 
C     obtained infeasibility reduction and the required and obtained
C     projected gradient norms was deleted.
C
C     February 16th, 2005.
C
C     1) The implementation of the matrix-vetor product subroutine that
C     takes care of the discontinuity of the second derivatives had
C     been made in GENCAN, inside calchddiff subroutine. The problem
C     with this choice was that calchddiff was calling evalc and evalgjc
C     subroutines to do this task. It was innapropriate as those
C     subroutines are subroutines of the augmented lagrangian framework
C     that are not present when GENCAN is being used stand-alone to
C     solve a bound-constrained problem. So, calchddiff subroutine of
C     GENCAN implements now the classical incremental quotients version
C     of the matrix-vector product approximation while the specific
C     implementation that takes care of the second-derivatives 
C     discontinuity of the augmented lagrangian function is now 
C     implemented in the evalhd subroutine that comes toghether with 
C     ALGENCAN.
C
C     IMPORTANT: as a consequence of this modification, argument hptype
C     of GENCAN must be set equal to 0 (user defined evalhd subroutine)
C     instead of 1 (calchddiff subroutine is used by GENCAN in this 
C     case) when GENCAN is called by ALGENCAN.
C
C     2) The common block that now contains array equatn was splitted
C     into probdata and cutedata. The information related to the 
C     constraints that is just needed by evalfc, evalgjc and evalc 
C     subroutines (that belong to the main algencanma program and not to 
C     the algencan subroutine) is now in the common block cutedata.

C     ******************************************************************
C     ******************************************************************
C
C     Report of modifications of gencan.
C
C     September 18th, 2006.
C
C     1) A bug in the order of the variables when writing format 1060
C     in algencan was fixed.
C
C     August 17th, 2005.
C
C     1) Once again, the stopping criterion for not enough progress in
C     the projected gradient norm was reformulated. Now, the method 
C     stops if there is not simple decrease of the sup-norm of the 
C     projected gradient over maxitngp consecutive iterations. It means
C     that the method stops if it is stagned at a stationary point of a
C     bad scaled problem. 
C
C     April 21th, 2005.
C
C     1) An option for computing the exact contribution of the linear
C     constraints in the Hessian-vector product was finally concluded.
C     The most difficult part was to obtain a version of the algorithm
C     using the same CUTE subroutines we were using in the past. It was
C     done the objective of being able to compare the effect of the 
C     exact calculation when compared with pure incremental quotients 
C     version without the interference of other factor (as, for example, 
C     to call different CUTE subroutines that use different times and
C     compute the values in different orders).
C
C     2) User subroutines evalf, evalg, evalc and evaljac were deleted
C     from the cute algencan version. Moreover, the way in which 
C     inequality constraints of the form cl <= c(x) <= cu were being
C     handle was also modified.
C
C     April 20th, 2005.
C
C     1) We add the array linear as a parameter of almost all the 
C     subroutines. The array indicates, for each constraint wether it is
C     a linear constraint (.true.) or not (.false.). This information
C     is currently being used for the computation of Hessian-vector
C     product in the Conjugate Gradients subalgorithm.
C
C     2) The contribution of the linear constraints in the 
C     Hessian-vector product is now being computed exactly instead of
C     approximating it by incremental quotients. Assuming that
C     computing the gradient of a linear constraint is (computationally)
C     cheap we gave preference to storage space savings instead of
C     time savings. So, the gradients of the linear constraints are
C     being computed every time they are necessary instead of being
C     saved.
C
C     April 12th, 2005.
C
C     1) When the initial "trust region radius" of CG is given by the
C     user, instead of using it directly, we use thea maximum between
C     the initial trust regions radius given by the user and delmin.
C
C     March 22th, 2005
C
C     1) The working vector wd of dimension ( 8 * n ) was splitted into
C     8 working vectors of dimension n. It was done in order to avoid
C     an artificial limit in the dimension of the biggest problem that
C     the method can solve. Now, the limitation is given by the maximum
C     integer number that the Fortran can represent, and not by this
C     number divided by 8.
C
C     March 18th, 2005
C
C     1) Some small changes that affect the computation of delta in the
C     cases in which the Sup-norm trust region is used in CG. These
C     changes implied in the calculation of the sup-norm of the step s
C     to be used instead of the so far used Euclidian norm.
C
C     2) In addtion to the previous change, another small change that
C     should not have any effect was made. We change the name of the
C     variable xnorm to xeucn (Eunclidian norm of x) and then, we change
C     it again to xeucn2 (squared Euclidian norm). So, in the places 
C     where xeucn were used now it is being used sqrt( xeucn2 ).
C
C     3) The change in (1) was motivated for a poor gencan performance
C     in very big problems. It was observed that the default value of
C     delmin (0.1) is very small for big problems when the 
C     Euclidian-norm trust-region is being used. So, in this cases, 
C     delmin should be increased or the Sup-norm trust-region should be 
C     used.
C
C     March 11th, 2005
C
C     1) equatn is not any more in a common block. It is now a parameter
C     of almost all the subrotines.
C
C     2) We modified the calling sequence moving gtype and hptype near
C     to the relative and absolute steps, eps and infs.
C 
C     March 10th, 2005.
C
C     1) We replaced "Convergence with an Armijo-like stopping criterion"
C     by "Armijo-like criterion satisfied" in the explanation of the line
C     serch flags.
C
C     March 9th, 2005.
C
C     1) The computation of the maximum number of conjugate gradient
C     iterations considers a linear relation between the norm (Euclidian
C     or Sup) of the projected gradient and a quantity which goes from
C     10 * log10( nind ) to nind. In the past, it was detected that it
C     is necessary to consider the maximum between 10 * log10( nind ) and
C     1. Now, we realise that we also need to consider the minimum 
C     between this quantity and nind. Moreover, we do not want (based
C     on Packmol) that the number of conjugate gradient iterations be 
C     greater than 10,000. So, we also add this upper bound. 
C
C     March 7th, 2005.
C
C     1) The integer ucgmaxit argument was replaced by two double 
C     precision arguments ucgmia and ucgmib. Now, if both, ucgmia and
C     ucgmib are setted by the user with strictly positive values, the
C     maximum allowed number of CG iterations for a subproblem of 
C     dimension nind is max( 1, int( ucgmia * nind + ucgmib ) ).
C     In the previous version this number was ucgmaxit, without any
C     dependency of the subproblem dimension, which, clearly, was not
C     useful at all. 
C
C     2) Doing this modification, we realize that, as it is, for solving
C     each subproblem CG is forced to, at least, 4 iterations. This fact
C     needs to be studied.
C
C     February 18th, 2005.
C
C     1) An unsed format statement, previously used to automaticaly
C     generates some tables, was deleted.
C
C     2) An unmateched parenthesis was corrected in the format
C     statement used to stop GENCAN due to a small step in a line search.
C
C     February 16th, 2005.
C
C     1) The evalhd subroutine used by default in GENCAN is now the one
C     implemented in calchddiff, which approximates the Hessian-vector
C     product by incremental quotients. The implementation used to 
C     overcome the non twice continuously differentiability of the 
C     classical (PHR) Augmented Lagrangian function is now part of 
C     ALGENCAN (and not GENCAN). So, to use GENCAN inside ALGENCAN, 
C     hptype argument must be set equal to 0 (ZERO).
C
C     2) The commented version of the empty function evalhd that must
C     be added when GENCAN is beinf used stand-alone was wrong. The
C     arguments declarations had been copied from evalnal. It was 
C     corrected.
C
C     November 10th, 2004.
C
C     1) After several test, all references to nonmontone line search
C     schemes were deleted.
C
C     September 28th, 2004.
C
C     1) Subroutines were checked an some absent arguments explanations
C     were added
C
C     2) Some calling sequences were modified to group related arguments
C
C     3) Arguments and local variables declarations were reordered in
C     alphabetical order.
C
C     3) Shrink and expand subroutines were modified to deal with just
C     one vector at a time. In this way, they are now being called from
C     calc* subroutines.
C
C     September 27th, 2004.
C
C     1) All comments were arranged to fit into the 72-columns format
C
C     2) Unused variable goth, which was prepared to indicate whether 
C     the Hessian matrix have been evaluated at the current point, was 
C     deleted from CG subroutine.
C
C     3) A spell check was used to correct the comments
C
C     September 21th, 2004.
C
C     1) In the stopping criterion where the progress in the objective 
C     function is verified, ''itnfp .ge. maxitnfp'' was changed for 
C     ''itnfp .gt. maxitnfp'', to make the choice maxitnfp equal to 1 
C     sounds reasonable.
C
C     2) Moreover, the previous chance came from the addition in the 
C     comments of GENCAN of the ''constraints'' information which makes 
C     clear to the user the values each argument may assume.
C
C     3) In the calculations of the first ''trust-radius'' for Conjugate 
C     Gradients, ''if( udelta0 .lt. 0.d0 ) then'' was changed by ''if 
C     ( udelta0 .le. 0.0d0 ) then'' to also make the default GENCAN 
C     choice of this initial trust-radius in the case of the user have 
C     been setted udelta = 0 by mistake.
C
C     4) The same for ucgmaxit.
C
C     5) In the line search subroutines spgls and tnls, ''if ( interp 
C     .gt. mininterp .and. samep ) then'' was changes by ''.ge.''.
C
C     6) Some comments of GENCAN arguments were re-written.
C
C     September 16th, 2004.
C
C     1) With the reconfiguration of the calc* subroutines (see (1) 
C     below) there were a number of redundant parameters in calchd and 
C     evalhd subroutines. These parameters were eliminated.
C
C     September 13th, 2004.
C
C     1) Subroutines named calc* that work in the reduced space always
C     call the corresponding eval* subroutine. As it was, calcnal (that
C     computes the gradient in the reduced space) called evalg or 
C     evalgdiff depending on gtype parameter. The same was for calchd. 
C     Now, calcnal calls evalg, calchd calls evalhd, and calchddiff (new) 
C     approximates the Hessian times vector product by incremental 
C     quotients calling calcnal or calcnaldiff depending on gtype parameter.
C     An improvement of this modification is that calcnal does not call 
C     evalg or evalgdiff (both work in the full space) any more but it 
C     approximates the gradient vector in the reduced space (by central 
C     finite differences) calling 2 * nind times evalf subroutine.
C
C     2) Some comments were added inside evalg and evalhd user supplied
C     subroutines alerting about the relation of these subroutines and
C     the parameters gtype and hptype, respectively.
C
C     3) Description of tnls subroutine was slightly modified.
C
C     4) The description of hptype parameter in gencan was again 
C     slightly modified.
C
C     5) With the introduction of the parameter lambda (that in the
C     context of Augmented Lagrangians is used to store the 
C     approximation of the Lagrange multipliers) the name of the 
C     variable used for spectral steplength was changed from lambda to 
C     lamspg. In addition, lammax was changed to lspgma and lammin to 
C     lspgmi.
C
C     6) Modifications introduced in June 15th, 2004 and May 5th, 2004
C     were, in fact, made in this version on September 13th, 2004.
C
C     June 15th, 2004.
C
C     1) The fmin stopping criterion and the maximum number of
C     functional evaluation stopping criterion were erroneously being 
C     tested before the main loop. It was just redundant and, for this 
C     reason, deleted.
C
C     May 5th, 2004.
C
C     1) Incorporated into an Augmented Lagrangian framework.
C
C     a) evalf and evalg were renamed as evalal and evalnal, 
C        respectively.
C
C     b) m,lambda,rho were added as parameters of the subroutines evalal 
C        and evalnal, and, as a consequence, as parameters of almost all 
C        the other subroutines.
C
C     2) The comment of hptype parameter of gencan was in portuguese
C     and it was translated into english.
C
C     3) A nonmonotone version of gencan is starting to be studied.
C     Parameters p and lastfv(0:p-1) were added to gencan, spgls, and
C     tnls to allow a nonmonotone line search. Array lastfv is now 
C     been updated for saving the last p functional values and the 
C     nonmonotone line searches are been done in a SPG or a 
C     Truncated Newton direction. p = 1 means monotone line search 
C     and is recommended until this study finish.
C
C     April 13th, 2004.
C
C     1) The modifications introduced in the occasion of the IRLOC 
C     development and re-development (October 21th, 2003 and February 
C     19th, 2003, respectively) were in fact made in this version on 
C     April 13th, 2004. The motivation to do this was to unify two 
C     parallel and different version of GENCAN (created, obviously, by 
C     mistake).
C
C     2) The complete reference of the GENCAN paper was finally added.
C
C     May 14th, 2003.
c
C     1) The way amax2 and amax2n were being computing may caused a 
C     segmentation fault. Its initialization was changed from infty and
C     -infty to 1.0d+99 and -1.0d+99, respectively. Using infty, when
C     combined with a big trust region radius, the final value of amax2
C     or amax2n may cause the impression that a bound is being attained, 
C     when it is not. "Redundant" ifs inside the amax2 and anax2n 
C     calculation were deleted. It should considered the possibility of 
C     using two constants, namely, bignum = 1.0d+20 and infty = 1.0d+99, 
C     instead of just infty. 
C
C     Modification introduced in October 21, 2003 in occasion of the
C     IRLOC re-development:
C
C     1) The stooping criteria related to functional value smaller than
C     fmin and exhaustion of maximum allowed number of functional 
C     evaluations have been done after the line search. And the 
C     questions were done as "if line search flag is equal to 4" or "if 
C     line search flag is equal to 8". But it was wrong in the case, for 
C     example, inside the line search, a functional value such that f <= 
C     fmin and the Armijo criterion was satisfied. In such case, the 
C     line search flag was being setted to 0 and not to 4. And gencan 
C     did not stop by the fmin criterion. Now, both stooping criteria 
C     are tested at the begining of the main gencan loop and just the 
C     stooping criteria by small line search step is tested after the 
C     line search.
C
C     Modification introduced in February 19, 2003 in occasion of the
C     IRLOC development:
C
C     1) The description of epsnfp parameter of GENCAN was modified. It
C     was written that to inhibit the related stopping criterion (lack
C     of function progress) it was necessary just set epsnfp = 0 when
C     it is also necessary to set maxitnfp = maxit. it was added in the
C     explanation.
C
C     2) In the explanation at the beginning of GENCAN it was written 
C     that cgscre parameter should be double precision. This comment was 
C     wrong. The correct type for cgscre parameter is integer.
C
C     Modifications introduced near April 1st 2003 in occasion of the 
C     PHR and inequality-constraints Augmented Lagrangian methods 
C     development:
C
C     1) The use of iprint was redefined and iprint2 was deleted.
C
C     2) The way to detect no progress in the log of the projected 
C     gradient norm was changed. As it was, ''no progress'' means no
C     reduction in the projected gradient norm over M iterations.
C     But this criterion implicitly assumed that the projected
C     gradient norm must decrease monotonously. Is it is clearly not
C     true, the criterion was changed by a non-monotone decrease
C     criterion. Now, progress means that the projected gradient
C     norm is, at each iteration, smaller than the maximum over the
C     last M iterations. And "no progress" means the it does not 
C     occurs during  not smaller than the 
C
C     3 ) The computation of qamaxn inside cg subroutine was in the 
C     wrong place (it was being used before computed) and it may was 
C     the reason for which the option nearlyq = .true. never worked 
C     properly. With this correction this option should be tested again.
C
C     On September 29th, 2004, we did a new test using the 41 bound
C     constrained problems with quadratic objective function from the
C     CUTE collection. The behaviour of GENCAN setting nearly equal
C     to true or false was indistinguishable. The test did not
C     include the different choices for the maximum number of CG
C     iterations being restricted to evaluate the different
C     alternatives for the case of finding a direction d such that
C     d^t H d <= 0. As a conclusion of this experiment we continue
C     recommending as a default choice to set nearlyq equal to false.
C
C     Modifications introduced from March 1st to March 21th of 2002
C     in occasion of the ISPG development:
C
C     1) Comments of some new parameters introduced in the previous
C     modification
C
C     2) As it was, in the first iteration of GENCAN (when kappa takes
C     value equal 1) and for one-dimensional faces, cgmaxit(the maximum 
C     number of Conjugate Gradient iterations to compute the internal to
C     the face truncated-Newton direction) was being 0. As it is 
C     obviously wrong, we add a max between what was being computed and 
C     one to allow at least one CG iteration.
C
C     3) Parameter inform in subroutines evalf, evalg and evalhd 
C     supplied by the user was added
C
C     Modifications introduced from May 31th to November 2nd of 2001
C     in occasion of the ALGENCAN development:
C
C     Fixed bugs:
C
C     1) The first spectral steplength was not been projected in the
C     [lspgmi,lspgma] interval.
C
C     2) The conjugate gradients accuracy (cgeps) which is linearly
C     dependent of the Euclidian norm of the projected gradient, was
C     also not been projected in the interval [cgepsi,cgepsf].
C
C     3) Conjugate gradients said that it was being used an Euclidian
C     norm trust region when it has really being used an infinite norm
C     trust region and viceversa.
C
C     4) Sometimes, the analytic gradient has been used although the
C     user choose the finite differences option.
C
C     Modifications:
C
C     1) To avoid roundoff errors, an explicit detection of at least one
C     variable reaching its bound when a maximum step is being made was
C     added.
C
C     2) The way in which two points were considered very similar in, 
C     for example, the interpolations and the extrapolations (which was 
C     dependent of the infinity norm of the points) showed to be very
C     scale dependent. A new version which test the difference 
C     coordinate to coordinate was done. In this was the calculus of the 
C     current point x and the descent direction sup-norm is not done any
C     more.
C
C     3) The same constants epsrel and epsabs were used as small 
C     relative and absolute values for, for example, detecting similar
C     points and for finite differences. Now, epsrel and epsabs are used 
C     for detecting similar points (and the recommended values are 
C     10^{-10} and 10^{-20}, respectively) and new constants sterel and 
C     steabs were introduced for finite differences (and the recommended 
C     values are 10^{-7} and 10^{-10}, respectively).
C
C     4) Two new stopping criteria for CG were added: (i) we stop if
C     two consecutive iterates are too  close; and (ii) we also
C     stop if there is no enough quadratic model progress during
C     maxitnqmp iterations.
C
C     5) The linear relation between the conjugate gradient accuracy
C     and the norm of the projected gradient can be computed using
C     the Euclidian- and the sup-norm of the projected gradient (only
C     Euclidian norm version was present in the previous version. The
C     linear relation is such that the CG accuracy is cgepsi when the
C     projected gradient norm value is equal to the value corresponding
C     to the initial guess and the CG accuracy is cgepsf when the
C     projected gradient norm value is cgrelf).
C
C     6) Inside Conjugate Gradients, the Euclidian-norm is been computed 
C     using an algorithm developed by C.L.LAWSON, 1978 JAN 08. Numerical 
C     experiments showed that the performance of GENCAN depends 
C     basically on the conjugate gradients performance and stopping
C     criteria and that the conjugate gradients depends on the way the
C     Euclidian-norm is been computed. These things deserve further 
C     research.
C
C     7) In the Augmented Lagrangian algorithm ALGENCAN, which uses
C     GENCAN to solve the bounded constrained subproblems, the maximum
C     number of Conjugate Gradients iterations (cgmaxit), which in this
C     version is linearly dependent of the projected gradient norm, was 
C     set to 2 * (# of free variables). As CG is not using restarts we 
C     do not know very well what this means. On the other hand, the 
C     accuracy (given by cgeps) continues being more strict when we are 
C     near to the solution and less strict when we ar far from the 
C     solution. 
C
C     8) Many things in the output were changed.


c
c AUXILIAR SUBROUTINES FOR THE CYLINDER PROBLEM
c

c
c Subroutine that computes de lengthgth of a string
c

      integer function length(string)

      character*200 string

      length = 200
      do while(string(length:length).le.' '.and.length.gt.0)
        length = length - 1
      end do

      return
      end

c
c Subroutine that gets the keyword from a line of the input file
c

      
      function keyword(string)
      character*200 keyword, string
      integer i, if, il

      i = 1      
      do while(string(i:i).le.' ') 
        i = i + 1
      end do
      if = i
      do while(string(i:i).gt.' ')
        i = i + 1
      end do
      il = i - 1
      keyword = string(if:il)

      return
      end

c
c Subroutine that gets the keyword from a line of the input file
c

      
      function keyvalue(string,n)
      character*200 keyvalue, string
      integer i, j, if, il, n

      i = 1      
      do while(string(i:i).le.' ') 
        i = i + 1
      end do
      do while(string(i:i).gt.' ')
        i = i + 1
      end do
      do j = 1, n
        do while(string(i:i).le.' ') 
          i = i + 1
        end do
        if = i
        do while(string(i:i).gt.' ')
          i = i + 1
        end do
        il = i - 1
      end do
      keyvalue = string(if:il)

      return
      end

c
c Subroutine that sets the bounds for the variables
c

      subroutine xini(n,x,l,u,icall,itype,xa)

      integer n, icall, itype
      double precision x(*), l(*), u(*), cmx, cmy, cmz,
     +                 xmax, ymax, zmax, xmin, ymin, zmin,
     +                 seed, drand, xnorm, xa(*)
      common/bounds/cmx,cmy,cmz,xmin,ymin,zmin,xmax,ymax,zmax
             
c The initial point for the reference point of the cylinder axis is the
c center of mass

      seed = 11891911.d0 + dfloat(icall)
      if(itype.eq.1) then
        x(1) = l(1) + drand(seed)*( u(1) - l(1) )
        x(2) = l(2) + drand(seed)*( u(2) - l(2) )
        x(3) = l(3) + drand(seed)*( u(3) - l(3) )

c Initial approximation for the axis direction is the z-direction

        x(4) = drand(seed)
        x(5) = drand(seed)
        x(6) = drand(seed)
        xnorm = dsqrt(x(4)**2 + x(5)**2 + x(6)**2)
        x(4) = x(4) / xnorm
        x(5) = x(5) / xnorm
        x(6) = x(6) / xnorm

      else if(itype.eq.2) then

        x(1) = xa(1) - 1.d3 + 2.d-3*drand(seed)
        x(2) = xa(2) - 1.d3 + 2.d-3*drand(seed)
        x(3) = xa(3) - 1.d3 + 2.d-3*drand(seed)
        x(4) = xa(4)
        x(5) = xa(5)
        x(6) = xa(6)

      end if

c Initial radius
     
      x(7) = 1.d-2*(xmax - xmin)+drand(seed)*1.d-2*(ymax-ymin)
                              
      return
      end

c
c Subroutine that writes the output in vmd format
c

      subroutine tovmd(file,x,step)

      integer length
      double precision x(*), step
      character*200 file

      step = 30.
      open(10,file=file(1:length(file)-4)//'.vmd')
      write(10,100) file(1:length(file)),
     +              x(1) - x(4)*step/2.d0, 
     +              x(2) - x(5)*step/2.d0, 
     +              x(3) - x(6)*step/2.d0, 
     +              x(1) + x(4)*step/2.d0, 
     +              x(2) + x(5)*step/2.d0, 
     +              x(3) + x(6)*step/2.d0,
     +              x(7) 
      close(10)

100   format(
     +'#!/usr/local/bin/vmd',/,
     +'mol new ',a, 
     +' type pdb first 0 last -1 step 1 filebonds 1', 
     +' autobonds 1 waitfor all',/,
     +'graphics top cylinder', 
     +' {',3(tr1,f8.3),' }', 
     +' {',3(tr1,f8.3),' }', 
     +' radius ',f8.3,
     +' resolution 30 filled 0',/,
     +'mol delrep 0 top',/,
     +'mol representation NewCartoon 0.300000 6.000000 4.100000 0',/,
     +'mol color Structure',/,
     +'mol selection {all}',/,
     +'mol material Opaque',/,
     +'mol addrep top',/,
     +'mol selupdate 0 top 0',/,
     +'mol colupdate 0 top 0',/,
     +'mol scaleminmax top 0 0.000000 0.000000',/,
     +'mol smoothrep top 0 0',/,
     +'mol drawframes top 0 {now}',/,
     +'set topmol [molinfo top]')    

      return
      end

















