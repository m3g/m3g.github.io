module topolink_data

  use ioformat, only : max_string_length
  implicit none

  !
  ! Types of data
  ! 
  
  ! PDB residue data

  type pdbresidue

    integer :: index, firstatom, lastatom
    character(len=4) :: name, chain = "    "
    character(len=1) :: conformation
    logical :: accessible
  
  end type pdbresidue

  ! PDB atom data

  type pdbatom

    integer :: index 
    type(pdbresidue) :: residue
    character(len=4) :: name 
    double precision :: x, y, z
    real :: b = 0., occup = 0.
    logical :: accessible

  end type pdbatom

  ! Specific link data

  type experiment_in_link

    ! observed: The atom pair was observed to link in this experiment
    ! type_reactive: The atom pair is reactive, according to the link types of this experiment
    ! obs_reactive: The atom pair is reactive, according to the observed reactivity in this experiment
    ! type_consistent: The topological distance of this atom pair is consistent with 
    !                  the observed reactivity of this atom pair, according to atom types
    ! obs_consistent: The topological distance is consistent with the observed reactivity of this
    !                 atom pair, acoording to the observed atom reactivities
     
    logical :: observed, type_reactive, obs_reactive, type_consistent, obs_consistent

  end type experiment_in_link

  type specific_link

    ! atom1: First atom of this pair
    ! atom2: Second atom of this pair
    ! nbeads: Number of beads of the linker, for optimization
    ! status: Final classification of the linker
    ! n_type_expected: Number of expected links, according to atom types
    ! n_obs_expected: Number of expected links, according to observed reactivity
    ! n_type_consistent: Number of observations consistent with type-reactivity
    ! n_obs_consistent: Number of observations consistent with observed reactivity
    ! euclidean: euclidean distance
    ! topodist: topological distance, if computed 
    ! dmaxlink: Maximum length of the linker, according to experimental linkers used
    ! dmax: Maximum distance consistent with observations
    ! dmin: Minimum distance consistent with observations
    ! dsearch: Maximum distance to search for a topological distance 
    ! observed: Linker was observed to form in some experiment
    ! type_reactive: Both atoms are reactive according to some linker types
    ! obs_reactive: Both atoms are reactive according to some observed reactivity
    ! found: A topological distance was found
    ! exp: Data for the experiments, concerning this atom pair

    type(pdbatom) :: atom1, atom2
    integer :: nbeads, status
    integer :: n_type_expected, n_obs_expected, n_type_consistent, n_obs_consistent
    double precision :: euclidean, topodist, dmaxlink, dmax, dmin, dsearch
    logical :: observed, type_reactive, obs_reactive, found
    type(experiment_in_link), allocatable :: exp(:)

  end type specific_link

  ! The data type for a link type

  type link_type
   
    type(pdbatom) :: atom1, atom2
    double precision :: dist

  end type link_type
 
  ! Observed link type

  type observed_link
      
    integer :: type
    double precision :: score
    type(pdbresidue) :: residue1, residue2

  end type observed_link

  ! Dead Ends

  type observed_deadend

    type(pdbresidue) :: residue   

  end type observed_deadend

  ! Experiment data

  type experiment_data

    integer :: nobs, ntypes, ndeadends, ngood, nbad
    integer :: nreactive_type, nreactive_obs
    integer :: nreach_type, nreach_obs
    integer :: noutreach_type, noutreach_obs
    integer :: nmiss_type, nmiss_obs
    character(len=max_string_length) :: name
    type(observed_deadend), allocatable :: deadend(:)
    type(observed_link), allocatable :: observed(:)
    type(link_type), allocatable :: linktype(:)
    double precision :: likelihood, userlikelihood, pgood, pbad, score

  end type experiment_data

  ! Model Data, for result analysis

  type modeldata

    ! Model log file name
    character(len=max_string_length) :: name
    ! Number of links in this file
    integer :: nlinks
    ! Score (rosetta?, TM-score?, GDT?, G-score) 
    double precision :: score
    ! Degree (P-value of G-score)
    double precision :: degree
    ! Link results for this model
    type(specific_link), allocatable :: link(:)
    ! Results
    integer :: nobscons  ! RESULT0: Number of observations that are consistent with the structure.
    integer :: ntopcons  ! RESULT1: Number of topological distances consistent with all observations.
    integer :: ntopnot   ! RESULT2: Number of topological distances NOT consistent with observations.
    integer :: nmiss     ! RESULT3: Number of missing links.
    integer :: nminmax   ! Number of links with min and max data that are consistent
    integer :: nobsgood  ! Number of observed links that are consistent with observations
    double precision :: sumscores  ! RESULT4: Sum of scores of observed links of all experiments.
    double precision :: likeli     ! RESULT5: Likelihood of the set of experimental results.
    double precision :: loglikeli  ! RESULT6: Log-likelihood of the set of experimental results.
    double precision :: usrlike    ! RESULT7: Likelihood of the set of experimental results. (with user pgood and pbad)
    double precision :: usrloglike ! RESULT8: Log-likelihood of the set of experimental results. (with user pgood and pbad)
    ! Index of link in overall model links lists
    integer, allocatable :: linkindex(:)

  end type modeldata

  ! 
  !  Functions 
  ! 

  contains

     ! subroutine that reads atom information from a PDB line

     function read_atom(record,error)

       use ioformat, only : max_string_length
       integer :: ioerr
       logical :: error
       type(pdbatom) :: read_atom
       character(len=max_string_length) :: record
       
       error = .false.
       if ( record(1:4) == "ATOM" .or. record(1:6) == "HETATM" ) then

         read(record(13:16),*,iostat=ioerr) read_atom%name
         read_atom%name = trim(adjustl(read_atom%name))
         if ( ioerr /= 0 ) error = .true.

         read(record(17:21),*,iostat=ioerr) read_atom%residue%name
         read_atom%residue%name = trim(adjustl(read_atom%residue%name))
         if ( ioerr /= 0 ) error = .true.
         call alternate_conformation(read_atom%residue)

         if ( record(22:22) /= " " ) then
           read(record(22:22),*,iostat=ioerr) read_atom%residue%chain
           if ( ioerr /= 0 ) error = .true.
           read_atom%residue%chain = trim(adjustl(read_atom%residue%chain))
         else
           read_atom%residue%chain = "0"
         end if

         read(record(23:26),*,iostat=ioerr) read_atom%residue%index
         if ( ioerr /= 0 ) error = .true.

         read(record(31:38),*,iostat=ioerr) read_atom%x
         if ( ioerr /= 0 ) error = .true.

         read(record(39:46),*,iostat=ioerr) read_atom%y
         if ( ioerr /= 0 ) error = .true.

         read(record(47:54),*,iostat=ioerr) read_atom%z
         if ( ioerr /= 0 ) error = .true.

       else
         error = .true.
       end if
       
     end function read_atom

     ! Function that checks if an atom is a Hydrogen atom

     function ishydrogen(atom)
      
       integer :: i, j, ioerr
       character(len=4) :: namestring
       logical :: ishydrogen
       type(pdbatom) :: atom

       ishydrogen = .false.
       namestring = adjustl(atom%name)
       i = 0
       do while( i < 4 )
         i = i + 1
         read(namestring(i:i),*,iostat=ioerr) j
         if ( ioerr == 0 ) cycle
         if ( namestring(i:i) == "H" ) ishydrogen = .true.
         exit
       end do

     end function ishydrogen

     ! Prints the data of an atom

     function print_atom(atom)

        type(pdbatom) :: atom
        character(len=25) print_atom

        write(print_atom,"( 2(tr2,a4),tr2,i5,tr2,a4 )") &
                                atom%residue%name, &
                                atom%residue%chain, &
                                atom%residue%index, &
                                atom%name


     end function print_atom

     ! Reads the data of a computed link from a previous log file

     function read_link(record)

        use ioformat, only : max_string_length
        implicit none
        integer :: i
        character(len=max_string_length) :: record
        character(len=13) :: charstat
        character(len=3) :: charobs
        character(len=9) :: chardmax
        type(specific_link) :: read_link
      
        read(record(9:12),*) read_link%atom1%residue%name
        read(record(14:14),*) read_link%atom1%residue%chain
        read(record(16:19),*) read_link%atom1%residue%index
        read(record(21:24),*) read_link%atom1%name
        
        read(record(26:29),*) read_link%atom2%residue%name
        read(record(31:31),*) read_link%atom2%residue%chain
        read(record(33:36),*) read_link%atom2%residue%index
        read(record(38:41),*) read_link%atom2%name
      
        read(record(43:50),*) read_link%euclidean

        read(record(64:66),*) charobs
        if ( charobs == "YES" ) then
          read_link%observed = .true.
        else
          read_link%observed = .false.
        end if

        read(record(69:77),*) read_link%dmin

        read(record(79:87),*) chardmax
        do i = 1, 9
          if ( chardmax(i:i) == ">" ) chardmax(i:i) = " "
        end do
        read(chardmax,*) read_link%dmax

        read(record(89:101),"( a13 )") charstat

        if ( charstat == "    OK: FOUND" ) then
          read_link%status = 0
          read_link%found = .true.
          read(record(52:60),*) read_link%topodist
        end if
        if ( charstat == "   BAD: SHORT" ) then
          read_link%status = 1 
          read_link%found = .true.
          read(record(52:60),*) read_link%topodist
        end if
        if ( charstat == "    BAD: LONG" ) then
          read_link%status = 2
          read_link%found = .true.
          read(record(52:60),*) read_link%topodist
        end if
        if ( charstat == "    BAD: EUCL" ) then
          read_link%status = 3
          read_link%found = .false.
          read_link%topodist = -1.d0
        end if
        if ( charstat == "BAD: NOTFOUND" ) then
          read_link%status = 4
          read_link%found = .false.
          read_link%topodist = -1.d0
        end if
        if ( charstat == " BAD: MISSING" ) then
          read_link%status = 5
          read_link%found = .true.
          read(record(52:60),*) read_link%topodist
        end if
        if ( charstat == "     OK: LONG" ) then
          read_link%status = 6
          read_link%found = .true.
          read(record(52:60),*) read_link%topodist
        end if
        if ( charstat == "     OK: EUCL" ) then
          read_link%status = 7
          read_link%found = .false.
          read_link%topodist = -1.d0
        end if
        if ( charstat == " OK: NOTFOUND" ) then
          read_link%status = 8
          read_link%found = .false.
          read_link%topodist = -1.d0
        end if

     end function read_link

     ! Prints the data of a general link

     function print_link(link)

       type(specific_link) :: link
       character(len=90) :: print_link

       write(print_link,"(2(tr2,a4),tr2,i5,2(tr2,a4),tr2,i5)") &
                          link%atom1%residue%name, &
                          link%atom1%residue%chain, &
                          link%atom1%residue%index, &
                          link%atom2%residue%name, &
                          link%atom2%residue%chain, &
                          link%atom2%residue%index

     end function print_link

     ! Prints the data of an observed link

     function print_obs(observed)

       type(observed_link) :: observed
       character(len=90) :: print_obs

       write(print_obs,"(2(tr2,a4),tr2,i5,2(tr2,a4),tr2,i5)") &
                          observed%residue1%name, &
                          observed%residue1%chain, &
                          observed%residue1%index, &
                          observed%residue2%name, &
                          observed%residue2%chain, &
                          observed%residue2%index

     end function print_obs

     ! Prints the data of an observed dead end

     function print_deadend(deadend)

       type(observed_deadend) :: deadend
       character(len=90) :: print_deadend

       write(print_deadend,"( 2(tr2,a4),tr2,i5 )") &
                              deadend%residue%name, &
                              deadend%residue%chain, &
                              deadend%residue%index

     end function print_deadend

     ! Prints the data of a linktype

     function print_linktype(link)

       type(link_type) :: link
       character(len=90) :: print_linktype
        
       write(print_linktype,"( 6(tr2,a4),tr2,f8.3)") &
                               link%atom1%residue%name,&
                               link%atom1%residue%chain,&
                               link%atom1%name,&
                               link%atom2%residue%name,&
                               link%atom2%residue%chain,&
                               link%atom2%name,&
                               link%dist

     end function print_linktype

     ! Prints a PDB atom line with coordinates
     
     function print_pdbatom(atom)

       use ioformat, only : max_string_length
       type(pdbatom) :: atom
       character(len=max_string_length) :: pdbformat, print_pdbatom

       pdbformat = "('ATOM',t7,i5,t13,a4,t18,a4,t22,a1,t23,i4,t31,f8.3,t39,f8.3,t47,f8.3,t55,f6.2,t61,f6.2)"
       write(print_pdbatom,pdbformat) &
                           atom%index, adjustl(atom%name), adjustl(atom%residue%name), atom%residue%chain, &
                           atom%residue%index, atom%x, atom%y, atom%z, atom%occup, atom%b

     end function print_pdbatom

     ! Prints a PDB hetero-atom line with coordinates
     
     function print_pdbhetatm(atom)

       use ioformat, only : max_string_length
       type(pdbatom) :: atom
       character(len=max_string_length) :: pdbformat, print_pdbhetatm

       pdbformat = "('HETATM',t7,i5,t13,a4,t18,a4,t22,a1,t23,i4,t31,f8.3,t39,f8.3,t47,f8.3,t55,f6.2,t61,f6.2)"
       write(print_pdbhetatm,pdbformat) &
                             atom%index, adjustl(atom%name), adjustl(atom%residue%name), atom%residue%chain, &
                             atom%residue%index, atom%x, atom%y, atom%z

     end function print_pdbhetatm

     !
     ! function isprotein: check if an atom is a protein atom
     !
     logical function isprotein(atom)
     
       implicit none
       type(pdbatom) :: atom
     
       isprotein = .false.
       select case ( atom%residue%name ) 
         case ( "ALA" ) ; isprotein = .true. ; return
         case ( "ARG" ) ; isprotein = .true. ; return
         case ( "ASN" ) ; isprotein = .true. ; return
         case ( "ASP" ) ; isprotein = .true. ; return
         case ( "ASX" ) ; isprotein = .true. ; return
         case ( "CYS" ) ; isprotein = .true. ; return
         case ( "GLU" ) ; isprotein = .true. ; return
         case ( "GLN" ) ; isprotein = .true. ; return
         case ( "GLX" ) ; isprotein = .true. ; return
         case ( "GLY" ) ; isprotein = .true. ; return
         case ( "HIS" ) ; isprotein = .true. ; return
         case ( "HSE" ) ; isprotein = .true. ; return
         case ( "HSD" ) ; isprotein = .true. ; return
         case ( "ILE" ) ; isprotein = .true. ; return
         case ( "LEU" ) ; isprotein = .true. ; return
         case ( "LYS" ) ; isprotein = .true. ; return
         case ( "MET" ) ; isprotein = .true. ; return
         case ( "PHE" ) ; isprotein = .true. ; return
         case ( "PRO" ) ; isprotein = .true. ; return
         case ( "SER" ) ; isprotein = .true. ; return
         case ( "THR" ) ; isprotein = .true. ; return
         case ( "TRP" ) ; isprotein = .true. ; return
         case ( "TYR" ) ; isprotein = .true. ; return
         case ( "VAL" ) ; isprotein = .true. ; return
         case default ; return
       end select
     
     end function isprotein

     !
     ! subroutine alternate_conformation: check if an atom is from an alternate
     !                                    conformation, and pick only the "A" conformation
     !
     !
     subroutine alternate_conformation(residue)
     
       implicit none
       type(pdbresidue) :: residue
     
       select case ( residue%name(2:4) ) 
         case ( "ALA" , &
                "ARG" , & 
                "ASN" , & 
                "ASP" , & 
                "ASX" , & 
                "CYS" , & 
                "GLU" , & 
                "GLN" , & 
                "GLX" , & 
                "GLY" , & 
                "HIS" , & 
                "HSE" , & 
                "HSD" , & 
                "ILE" , & 
                "LEU" , & 
                "LYS" , & 
                "MET" , & 
                "PHE" , & 
                "PRO" , & 
                "SER" , & 
                "THR" , & 
                "TRP" , & 
                "TYR" , & 
                "VAL" ) 
           residue%conformation = residue%name(1:1)
           residue%name = residue%name(2:4)
         case default 
           residue%conformation = " "
       end select
     
     end subroutine alternate_conformation

end module topolink_data

