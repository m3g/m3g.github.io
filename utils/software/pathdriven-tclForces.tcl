#
# pathdriven-tclForces.tcl
#
# L. Martinez - Institute of Chemistry - State University of Campinas
# leandro@iqm.unicamp.br
#
# Main script for the pathdriven method. Calls other scripts to
# do the work. This script is to be called from the NAMD input file,
# by using the TCL Boundary forces interface. 
# 
# The method is described in this reference, which must be cited if
# the script or the algorithm is used:
#
# L. Martínez, I. Polikarpov, M. S. Skaf, Only Subtle Protein
# Conformational Adaptations Are Required for Ligand Binding to Thyroid
# Hormone Receptors: Simulations Using a Novel Multipoint Steered
# Molecular Dynamics Approach, J. Phys. Chem. B, 112, 10741, 2008.
#
# More information: http://leandro.iqm.unicamp.br
#
# -----------------------------------------------------------------------------------------
#
# How to use this script:
#
# INPUT:
# 
# A script like the one below must be included in a NAMD simulation input
# file, which will run the production run in which the forces will be
# applied:
#
#    |#
#    |# Script to apply a smd-like force with variable direction
#    |#
#    |tclForces on
#    |tclForcesScript {
#    |#
#    |# Input parameters
#    |#
#    |set pdpath /home/Bis/leandro/path/pathdriv
#    |set coordinates ../psfgen/1K8Tcam.ready.pdb
#    |set patomsfile ./pulled.pdb
#    |set pathfile ./path.pdb
#    |set tol 0.
#    |set SMDk 4.00
#    |set harmk 0.25
#    |set trajsteps 2500000
#    |set OutputFreq 5000
#    |set OutputName ./1K8Tcam.open
#    |set restartpdfreq 1000
#    |#
#    |# For reestarting from crashed calculation
#    |#
#    |#set restartfromfile ./run2/open1.restart
#    |#
#    |# Start script
#    |#
#    |source $pdpath/pathdriven-tclForces.tcl
#    |}
#
#  The lines above must be added BEFORE the "run" directive of the
#  simulation. The parameters are:
#  
#  pdpath : path to the pathdriven-tclForces.tcl script (this one)
#  coordinates : the PDB coordinates file of the simulation (the actual
#                coordinates are not used, this file must be provided for
#                the definition of atom names, residues, etc.)
#  patomsfile : A PDB file (which can be similar to the coordinate file) in
#               which the B-factor column of the defines to which atoms
#               the pulling force will be applied. Atoms with B-factors
#               different from 0.00 will be considered as atoms to which
#               the forces will be applied. 
#  pathfile : A PDB file containing the desired trajectory for the atoms
#             to which the forces will be applied. Therefore, this file
#             contains a SEQUENCE of PDB coordinates for each atom,
#             which define which is the desired trajectory. For example,
#             if the forces will be applied to a single atom, which is
#             a CA carbon atom from a ARG282 residue:
#  
#             ATOM      5  CA  ARG X 292
#  
#             The pathfile must contain a sequence of coordinates for this
#             atom, in PDB format (with the same nomenclature as in the
#             coordinates file). For example:
#  
#             ATOM      5  CA  ARG X 292      11.034 189.455   5.417  1.00  0.00      1K8T
#             ATOM      5  CA  ARG X 292      11.178 189.770   5.499  1.00  0.00      1K8T    
#             ATOM      5  CA  ARG X 292      11.486 190.286   5.626  1.00  0.00      1K8T    
#             ATOM      5  CA  ARG X 292      11.928 190.847   6.373  1.00  0.00      1K8T    
#             ATOM      5  CA  ARG X 292      12.527 191.459   7.799  1.00  0.00      1K8T    
#  
#             This means that the forces will be applied such that this atom
#             is induced to move towards the first coordinates. When it
#             reaches (in the present method sense) this point, it starts
#             to be pulled in such a way to reach the second point, etc.
#             
#             There are not many restrictions to the format of the path
#             file. If different trajectories are to be used for different
#             atoms in the simulation, you need to set the sequential
#             coordinates for each atom. The atoms can be intercalated, or
#             not, it doesn't matter. It is important only that the format
#             of the atom, residue, number, segment, etc. is identical
#             to the ones in the coordinate file. The path file does not
#             need to contain data for all atoms of the simulation, only
#             of those to which forces will be applied. 
#  tol : The precision to which it is considered that a final point in 
#        the trajectory step was reached. Since the criterion is that
#        the atom crossed a plane containing the final point and
#        perpendicular to the current direction, this precision can
#        be set to zero. If greater than zero, it will mean that the
#        direction will be changed before the atoms cross that plane.
#  SMDk : The harmonic force constant of the Steered MD force, with its
#         usual meaning.
#  harmk : When an atom reaches the final end point of its trajectory,
#          you can harmonically restrain it to its final position. harmk
#          defines the force constant of this harmonic potential. 
#  trajsteps : The total length of the SMD simulation. Note that it may
#              be different from the "run" directive, because you may
#              be restarting your simulation from a previous run. This
#              value is used to compute the velocity of the SMD force
#              for each atom, in such a way that, roughly, each atom 
#              will reach the end point of the desired trajectory at the final
#              step. 
#  OutputFreq : The frequency of output of the data on the SMD pulling.
#  OutputName : The base name of the output files. 
#  restartpdfreq : The frequency in which the restart files are written.
#  
#  If a simulation crashed, you can restart it as described in the 
#  example, where:
#  
#  restartfromfile : the .restart file that was last produced by the 
#                    pathdriven script.
#  
#  OUTPUT:
#  
#  The script will output three files:
#  
#  1: .restart : file containing the data of the last restart output,
#                to be used for restarting the simulation from a 
#                previous one.
#  
#  2: .total.dat : A file containing the total force applied to all
#                  atoms, the total work, the total force modulus, 
#                  for each time-step.
#  
#  3: .atoms.dat : A file containing, for each time-step, for each
#                  atom, the path step it is, the force applied, 
#                  the displacement from reference coordinates, 
#                  the modulus of the force and the work. 
#  
# -----------------------------------------------------------------------------------------
#
# Should not be modified from here on.
#
# L. Martinez, October 2007. Institut Pasteur.
# leandromartinez98@gmail.com
#
# This script is used to apply the
# generalized steered molecular dynamics approach to a potentially large
# number of atoms. The following script is called only once during the
# simulation.
#
# It reads from a PDB file, in column B (b-factors) the
# atoms for which trajectory data will be provided. This number of each
# atom for which an external force will be applied is stored in the
# patoms array.
#
# It also reads the expected trajectory for each atom. This information
# will be stored in the traj array. The patoms array and the traj array
# contain the information for each atoms consecutivelly, in the same
# order.
#
# Finally, the script reads the simulation time (in nano-seconds),
# and the force constant for the pulling. The pulling velocity is
# atom-dependent, and will be automatically set according to the length
# of the trajectory (that will be computed later) and the simulation
# time.  
#
# Data comming out from this script:
#
# Scallars:
#   SMDk: The force constant
#   trajsteps: Total pulling time, in number of steps.
#   OutputFreq: Frequency of output.
#   npull: Number of pulled atoms.
# Names:
#   patomsfile: Name of the file containing the pulled atoms.
#   pathfile: Name of the path file.
#   OutputName: Base name for output files.
# Arrays:
#   pulled: For i=1 to npull, contains the pdb-file index for the atom.
#   iatom: For each value in pulled, contains the index of pulled containing
#          that value (the inverse vector of pulled, very useful).
#   nptraj: For i=1 to npull, contains the number of trajectory steps of
#           each pulled atom.
#   step(pulled counter, step counter, 3): Contains the coordinates of
#           each trajectory step, for each step, for each atom, this
#           ones labeled by their index in the pdb file. 
#
# L. Martinez, October 2007. Institut Pasteur.
# leandromartinez98@gmail.com
#
# Checking the existance of necessary data in the input file
#
if { [ info exists coordinates ] == 0 } {
  print "Error: coordinates file not provided."
  exit
}
if { [ info exists patomsfile ] == 0 } {
  print "Error: patomsfile not provided."
  exit
}
if { [ info exists pathfile ] == 0 } {
  print "Error: pathfile not provided."
  exit
}
if { [ info exists OutputName ] == 0 } {
  print "Error: OutputName not provided."
  exit
}
if { [ info exists SMDk ] == 0 } {
  print "Error: SMDk not provided."
  exit
}                   
if { [ info exists trajsteps ] == 0 } {
  print "Error: trajsteps not provided."
  exit
}                   
if { [ info exists OutputFreq ] == 0 } {
  print "Error: OutputFreq not provided."
  exit
}
if { [ info exists harmk ] == 0 } {
  print "Error: harmk not provided."
  exit
}
if { [ info exists restartpdfreq ] == 0 } {
  print "Error: restartpdfreq frequency not provided."
  exit
}
# 
# Reading file names and steered MD parameters (trajsteps, SMDk,
# OutputFreq, patomsfile and pathfile)
print "Will apply the generalized SMD approach to induce a path."
print "patomsfile: $patomsfile"
print "pathfile: $pathfile"
print "Base name for output files: $OutputName"
print "Force constant (SMDk): $SMDk"
print "Total pulling time: $trajsteps steps"
print "Output frequency: $OutputFreq"
print "Force constant for harmonic potential: $harmk "
print "Printing restart data at every $restartpdfreq steps"
#
# Reading the PDB file containing information on the atoms that will be
# pulled.
#
puts -nonewline "TCL: Reading $patomsfile ... "; flush stdout
set file $patomsfile
set file [ open $file r ]
set file [ read $file ]
set file [ split $file "\n" ]
set npull 0
foreach line $file {
  if { [ string range $line 0 3 ] == "ATOM" |
       [ string range $line 0 5 ] == "HETATM" } {
    if { [ expr [ string range $line 60 66 ] ] != "0" } {
      set atom [ string trim [ string range $line 12 15 ] ]
      set resid [ string trim [ string range $line 17 20 ] ]
      set resn [ string trim [ string range $line 23 26 ] ]
      set seg [ string trim [ string range $line 72 76 ] ]
      set pulled($atom$resid$resn$seg) 1
      incr npull
      set chars($npull) $atom$resid$resn$seg
    }
  }
}
#
# Checking if all selected atoms exist in the coordinate file of the
# simulation

set pdb [ open $coordinates r ]
set pdb [ read $pdb ]
set pdb [ split $pdb "\n" ]
for { set i 1 } { $i <= $npull } { incr i } {
  set exist 0
  foreach line $pdb {
    if { [ string range $line 0 3 ] == "ATOM" |
         [ string range $line 0 5 ] == "HETATM" } {
      set atom [ string trim [ string range $line 12 15 ] ]
      set resid [ string trim [ string range $line 17 20 ] ]
      set resn [ string trim [ string range $line 23 26 ] ]
      set seg [ string trim [ string range $line 72 76 ] ]
      if { $chars($i) == "$atom$resid$resn$seg" } {
        set exist 1; break 
      }
    }
  }
  if { $exist == 0 } {
    print "failed." 
    print "Error: atom to be pulled in patomsfile not found "
    print "       in coordinate file: $chars($i)"
    exit
  }
}
puts "done."
print "Number of atoms to be pulled: $npull"
if { $npull < 1 } { 
  print "Error: No atom to be pulled. "
  exit
}
#
# Now, from the coordinate file, get the number of the atoms and the
# order in which they appear in the actual simulation files, order that
# will be fundamental within the calcforces procedure. The include
# vector will be only used at the first time step to drop out atoms for
# which no external forces will be applied. The pulled array will be
# used only for checking purposes.
#
set pdb [ open $coordinates r ]
set pdb [ read $pdb ]
set pdb [ split $pdb "\n" ]
set i 0
foreach line $pdb { 
  if { [ string range $line 0 3 ] == "ATOM" |
       [ string range $line 0 5 ] == "HETATM" } {
    set atom [ string trim [ string range $line 12 15 ] ]
    set resid [ string trim [ string range $line 17 20 ] ]
    set resn [ string trim [ string range $line 23 26 ] ]
    set seg [ string trim [ string range $line 72 76 ] ]
    if { [ info exists pulled($atom$resid$resn$seg) ] == 1 } {
      incr i
      set ipull($atom$resid$resn$seg) $i 
      set iatom [ string trim [ string range $line 6 10 ] ]
      set include($iatom) 1
      set pdbind($i) $iatom 
      set chars($iatom) $atom$resid$resn$seg
    }
  }
}
#
# Reading the PDB file that contain a sequence of structures defining
# the trajectory of each atom. It is not required that all atoms appear
# in all structures, in such a way that each atom may have a different
# number of trajectory points. The script will therefore read all atom
# coordinates and whenever it finds atom coordinates for an atom listed
# in the pulled array, it will add a trajectory point to this atom. On
# the other side, the pulling velocities will be set according the
# total trajectory length and simulation time, therefore if the point
# has two or ten points, but the same trajectory length, the pulling
# velocity will be the same. This is to remark that the pulling velocity
# will be set for each atom in order that all atoms reach the final
# point at the same time (this may not necessarily happen though). This
# part of the script is probably slow, but fortunately it must only be
# called once. It also sets the arrays with larger dimensions, since it
# will store the trajectories for each pulled atom in memory.
#

# First set the number of trajectory points per atom as zero

for { set i 1 } { $i <= $npull } { incr i } {
  set nptraj($pdbind($i)) 0
}

# Reading the file

puts -nonewline \
     "TCL: Reading $pathfile ... This may take a while ... " 
      flush stdout
set file [ open $pathfile r ]
set file [ read $file ]
set file [ split $file "\n" ]
foreach line $file {
  if { [ string range $line 0 3 ] == "ATOM" |
       [ string range $line 0 5 ] == "HETATM" } {
    set atom [ string trim [ string range $line 12 15 ] ]
    set resid [ string trim [ string range $line 17 20 ] ]
    set resn [ string trim [ string range $line 23 26 ] ]
    set seg [ string trim [ string range $line 72 76 ] ]
    if { [ info exists pulled($atom$resid$resn$seg) ] } {
      set ip $ipull($atom$resid$resn$seg)
      set ip $pdbind($ip)
      incr nptraj($ip)
      set step($ip,$nptraj($ip),1) \
          [ string trim [ string range $line 30 37 ] ]
      set step($ip,$nptraj($ip),2) \
          [ string trim [ string range $line 38 45 ] ]
      set step($ip,$nptraj($ip),3) \
          [ string trim [ string range $line 46 53 ] ]
    }
  }
}
puts "done."
for { set i 1 } { $i <= $npull } { incr i } {
  if { $nptraj($pdbind($i)) < 2 } {
    print "Error: The trajectory of all atoms must have at least \
           two points. Atom $pdbind($i) doesn't."
    exit
  }
}
#
# Script that performs initial computations that are done
# only once (the path length, for instance, and the pulling
# velocities of each atom)
#
#
# initial.tcl
#
# This script will receive the data from readdata.tcl and compute the
# necessary things to start a path driven molecular dynamics. It is
# called only once and should not be editted.
#
# L. Martinez, October 2007, Institut Pasteur.
# leandromartinez98@gmail.com
#
# Procedure to do a power 2
#
proc p2 { val } {
  set p2 [ expr $val*$val ]
  return $p2
}
#
# Computing the length of the trajectory for each atom and setting the
# pulling velocity.
#
puts -nonewline "TCL: Computing the trajectory length of each atom ... "
flush stdout
for { set ii 1 } { $ii <= $npull } { incr ii } {
  set i $pdbind($ii)
  set tlength($i) 0.
  for { set j 2 } { $j <= $nptraj($i) } { incr j } {
    set dist [ expr \
        [ p2 [ expr $step($i,$j,1) - $step($i,[ expr $j - 1 ],1) ] ] + \
        [ p2 [ expr $step($i,$j,2) - $step($i,[ expr $j - 1 ],2) ] ] + \
        [ p2 [ expr $step($i,$j,3) - $step($i,[ expr $j - 1 ],3) ] ] ]
    set dist [ expr sqrt($dist) ]
    set tlength($i) [ expr $tlength($i) + $dist ]
  }
  set vnorm($i) [ expr $tlength($i) / $trajsteps ]
}
puts "done."
#                                                                               
# Writting some data of the path to file (can be removed)
#       
set file [ open $OutputName.pathinfo.dat w ]
puts $file "Information about the path:"
for { set ii 1 } { $ii <= $npull } { incr ii } {
  set i $pdbind($ii)
  set dist [ expr \
      [ p2 [ expr $step($i,$nptraj($i),1) - $step($i,1,1) ] ] + \
      [ p2 [ expr $step($i,$nptraj($i),2) - $step($i,1,2) ] ] + \
      [ p2 [ expr $step($i,$nptraj($i),3) - $step($i,1,3) ] ] ]
  set dist [ expr sqrt($dist) ]
  puts $file "$ii atom: $pdbind($ii) n_points: $nptraj($i)\
              length: $tlength($i) vnorm: $vnorm($i)\
              linear: $dist char: $chars($ii)"
}
close $file
#
# Reseting arrays that will be used in the simulation (are in separate
# loops just for clarity)
#
# Work performed by the force applied on each atom (must be computed
# every step)
#
for { set ii 1 } { $ii <= $npull } { incr ii } {
  set i $pdbind($ii)
  set work($i) 0.
}

#
# Path counter, sets in which trajectory step the atom currently is
#
for { set ii 1 } { $ii <= $npull } { incr ii } {
  set i $pdbind($ii)
  set pstep($i) 1
}
#
# First reference point (the position in the first step), and also the
# first position in the xlast vector that is used to compute the work
#
for { set ii 1 } { $ii <= $npull } { incr ii } {
  set i $pdbind($ii)
  set ref($i,1) $step($i,1,1) 
  set ref($i,2) $step($i,1,2) 
  set ref($i,3) $step($i,1,3) 
  set xlast($i,1) $step($i,1,1) 
  set xlast($i,2) $step($i,1,2) 
  set xlast($i,3) $step($i,1,3) 
}
#
# First pulling direction
#
for { set ii 1 } { $ii <= $npull } { incr ii } {
  set i $pdbind($ii)
  set vel($i,1) [ expr $step($i,2,1) - $step($i,1,1) ]
  set vel($i,2) [ expr $step($i,2,2) - $step($i,1,2) ]
  set vel($i,3) [ expr $step($i,2,3) - $step($i,1,3) ]
  set veln [ expr sqrt( [ p2 $vel($i,1) ] + \
                        [ p2 $vel($i,2) ] + \
                        [ p2 $vel($i,3) ] ) ]
  set vel($i,1) [ expr $vel($i,1) / $veln ] 
  set vel($i,2) [ expr $vel($i,2) / $veln ] 
  set vel($i,3) [ expr $vel($i,3) / $veln ] 
}
#
# Other parameters required for the control of the simulation
#
# Time step counter

set ts -1

# Force at timestep -1, the force at the last step is used to compute
# the work by w = flast * ( x - xlast )

for { set ii 1 } { $ii <= $npull } { incr ii } {
  set i $pdbind($ii)
  set flast($i,1) 0.
  set flast($i,2) 0.
  set flast($i,3) 0.
}

#
# Atoms that will be used in the calculations
#
for { set ii 1 } { $ii <= $npull } { incr ii } {
  addatom $pdbind($ii)
}

# If this is not a restart run, move old output files and start new ones

if { [ info exists restartfromfile ] == 0 } {

  if { [ file exists "$OutputName.atoms.dat" ] } {
    exec mv $OutputName.atoms.dat $OutputName.atoms.BAK
  } 
  if { [ file exists "$OutputName.total.dat" ] } {
    exec mv $OutputName.total.dat $OutputName.total.BAK
  }

  set file [ open $OutputName.atoms.dat w ]
  puts $file "# atom  path-step fx fy fz disp fmod work"
  close $file  

  set file [ open $OutputName.total.dat w ]
  puts $file "# time-step fx fy fz fmod work"
  close $file  
}

# If this is a restart run, continue existing output files and read
# reestart data to reestart the pulling where it was

if { [ info exists restartfromfile ] == 1 } {

  print " Restarting pulling from previous run: $restartfromfile "

  set file [ open $OutputName.atoms.dat a+ ]
  puts $file "# atom  path-step fx fy fz disp fmod work"
  close $file  

  set file [ open $OutputName.total.dat a+ ]
  puts $file "# time-step fx fy fz fmod work"
  close $file  

  set file [ open $restartfromfile r ]
  set restartfile [ read $file ]
  close $file
  set restartfile [ split $restartfile "\n" ]
  set i 0
  foreach line $restartfile {
    incr i
    set restartline($i) $line
  }
  set ts [ string trim $restartline(1) ]
  for { set ii 1 } { $ii <= $npull } { incr ii } {
    set i $pdbind($ii)
    set line $restartline([ expr $ii + 1 ])  
    set line [ split $line " " ]
    set j 0
    foreach value $line {
      if { [ string trim $value ] > " " } { incr j; set data($j) $value }
    }
    set pstep($i) $data(1)
    set work($i) $data(2)
    set flast($i,1) $data(3)
    set flast($i,2) $data(4)
    set flast($i,3) $data(5)
    set xlast($i,1) $data(6)
    set xlast($i,2) $data(7)
    set xlast($i,3) $data(8)
    set ref($i,1) $data(9)
    set ref($i,2) $data(10)
    set ref($i,3) $data(11)
    set vel($i,1) $data(12)
    set vel($i,2) $data(13)
    set vel($i,3) $data(14)
  }
}
#
# The script that contains the calcforces routine, that will be
# called at every step. Everything here must be optimized at maximum, so
# if you find anything to improove, do it (and tell me!). I'm not
# sure how to optimize calculations done with TCL scripting...
#
#
# allstep.tcl
#
# This scripts performs the computations at every step of the simulation
# to add forces to the atoms in order that they follow the trajectory.
# Everything must be optimized here to have good performance, but
# probably, if the pulling is applied to a lot of atoms, the simulation
# will be much slower, since computations in scripting languague are not
# efficient at all. It would be better to just include this in the code.
#
# L. Martinez, October 2007, Institut Pasteur
# leandromartinez98@gmail.com
#
proc calcforces { } {
  global OutputName OutputFreq vnorm vel step nptraj npull xlast\
         ref work pstep ts pdbind tol harmk SMDk tforce flast\
         restartpdfreq

# Increase time-step counter

  incr ts

# Print time step to restart file

  set restart 0 
  if { [ expr $ts % $restartpdfreq ] == 0 } { 
    set restart 1 
  }
  if { $restart == 1 } {
    catch { exec mv $OutputName.restart $OutputName.restart.BAK } msg
    if { $msg > " " } { 
      print "Created pathdriven restart file." 
    } else {
      print "Saving old pathdriven restart file..."
    }
    set restartfile [ open $OutputName.restart w ]
    puts $restartfile $ts
  }

# Check if there will be printing in this step

  set print 0
  if { [ expr $ts % $OutputFreq ] == 0 } { set print 1 }

# Reset total work

  set twork 0.

# Set total forces to zero (two values, a total vector force and a
# scalar sum of the moduli of the forces (which is probably more useful)

  if { $print == 1 } {
    set tforce(1) 0.
    set tforce(2) 0.
    set tforce(3) 0.
    set fsum 0.
    set file [ open $OutputName.atoms.dat a+ ]
    puts $file "TIME STEP $ts"
  }

# Iterate over all atoms using TclBC commands
 
  loadcoords coords
  for { set ii 1 } { $ii <= $npull } { incr ii } {

    set i $pdbind($ii)
    set xcoords $coords($i)
    set j 1; foreach val $xcoords { set x($j) $val; incr j }

# Print data to restart file 

    if { $restart == 1 } {
      puts $restartfile "$pstep($i) $work($i)\
                         $flast($i,1) $flast($i,2) $flast($i,3)\
                         $xlast($i,1) $xlast($i,2) $xlast($i,3)\
                         $ref($i,1) $ref($i,2) $ref($i,3)\
                         $vel($i,1) $vel($i,2) $vel($i,3)"
    }
 
#
# If the path has ended for this atom, simply apply an harmonic
# constranint to its position relative to the last target point
#

    if { $pstep($i) == $nptraj($i) } { 
      set dx [ expr $step($i,$nptraj($i),1) - $x(1) ]
      set dy [ expr $step($i,$nptraj($i),2) - $x(2) ]
      set dz [ expr $step($i,$nptraj($i),3) - $x(3) ]
      set fx [ expr 2. * $harmk * $dx ]
      set fy [ expr 2. * $harmk * $dy ]
      set fz [ expr 2. * $harmk * $dz ]
      addforce $i " $fx $fy $fz " 
      set work($i) [ expr $work($i) + $flast($i,1)*($x(1)-$xlast($i,1)) + \
                                      $flast($i,2)*($x(2)-$xlast($i,2)) + \
                                      $flast($i,3)*($x(3)-$xlast($i,3)) ]
      set twork [ expr $twork + $work($i) ]            
      if { $print == 1 } {
        set tforce(1) [ expr $tforce(1) + $fx ] 
        set tforce(2) [ expr $tforce(2) + $fy ] 
        set tforce(3) [ expr $tforce(3) + $fz ]
        set fmod [ expr sqrt( $fx*$fx + $fy*$fy + $fz*$fz ) ]
        set fsum [ expr $fsum + $fmod ]  
        set disp [ expr sqrt ( $dx*$dx + $dy*$dy + $dz*$dz ) ]
        puts $file [ format "%8i %8i %14.3f %14.3f %14.3f %14.3f %14.3f %14.3f" \
                     $i $pstep($i) $disp $fx $fy $fz $fmod $work($i) ]
      }
      set flast($i,1) $fx
      set flast($i,2) $fy
      set flast($i,3) $fz
      set xlast($i,1) $x(1)
      set xlast($i,2) $x(2)
      set xlast($i,3) $x(3)
      continue
    }

#
# Computing the displacement vector at this point
#

    set dx [ expr $x(1) - $ref($i,1) ]
    set dy [ expr $x(2) - $ref($i,2) ]
    set dz [ expr $x(3) - $ref($i,3) ]

#
# Computing the external force to be applied to this atom
#

    set fx [ expr $SMDk * ( $vnorm($i)*$vel($i,1) * $ts - $dx ) ]
    set fy [ expr $SMDk * ( $vnorm($i)*$vel($i,2) * $ts - $dy ) ]
    set fz [ expr $SMDk * ( $vnorm($i)*$vel($i,3) * $ts - $dz ) ]

# Will only apply the force if it is in the desired direction (meaning
# the the steered force will not help to resist to the movement)
                            
    set direc [ expr $fx*$vel($i,1) + \
                     $fy*$vel($i,2) + \
                     $fz*$vel($i,3) ]  
    if { $direc > 0. } {
      addforce $i " $fx $fy $fz "  
    } else {
      set fx 0.
      set fy 0.
      set fz 0.
    }
 
# Computing the work performed on this atom

    set work($i) [ expr $work($i) + $flast($i,1)*($x(1)-$xlast($i,1)) + \
                                    $flast($i,2)*($x(2)-$xlast($i,2)) + \
                                    $flast($i,3)*($x(3)-$xlast($i,3)) ]
    set twork [ expr $twork + $work($i) ]             

# Adding the force to the total force array and printing force
# information at this step, with displacement and work

    if { $print == 1 } {
      set tforce(1) [ expr $tforce(1) + $fx ] 
      set tforce(2) [ expr $tforce(2) + $fy ] 
      set tforce(3) [ expr $tforce(3) + $fz ]
      set fmod [ expr sqrt( $fx*$fx + $fy*$fy + $fz*$fz ) ]
      set fsum [ expr $fsum + $fmod ] 
      set disp [ expr sqrt ( $dx*$dx + $dy*$dy + $dz*$dz ) ]
      puts $file [ format "%8i %8i %14.3f %14.3f %14.3f %14.3f %14.3f %14.3f" \
                   $i $pstep($i) $disp $fx $fy $fz $fmod $work($i) ]
    }
 
# Updating the xlast and flast vectors for this atom

    set xlast($i,1) $x(1)
    set xlast($i,2) $x(2)
    set xlast($i,3) $x(3)
    set flast($i,1) $fx
    set flast($i,2) $fy
    set flast($i,3) $fz
                                            
#
# Now we check if the current path step has ended and we will perform
# all necessary operations to change to the next path step if that is
# the case, that means, we will update the reference position and the
# velocity vector
#
 
# Compute the distance of the pulled atoms to the plane that
# contains the next target point and is perpendicular to the current
# pulling direction. 

    set target [ expr $pstep($i) + 1 ]

    set dplan [ expr ( $step($i,$target,1) - $x(1) ) * $vel($i,1) + \
                     ( $step($i,$target,2) - $x(2) ) * $vel($i,2) + \
                     ( $step($i,$target,3) - $x(3) ) * $vel($i,3) ]

#
# If the path step has not ended, just go to the computations for the
# next atom
#

    if { $dplan > $tol } { continue }

#
# If the path has ended, then will update references and directions
#
 
# Update target point counter 

    incr pstep($i) 

# Update the velocity vector

    set vel($i,1) [ expr $step($i,$pstep($i),1) - $x(1) ]
    set vel($i,2) [ expr $step($i,$pstep($i),2) - $x(2) ]
    set vel($i,3) [ expr $step($i,$pstep($i),3) - $x(3) ]
    set veln [ expr sqrt( [ p2 $vel($i,1) ] + \
                          [ p2 $vel($i,2) ] + \
                          [ p2 $vel($i,3) ] ) ]
    set vel($i,1) [ expr $vel($i,1) / $veln ]
    set vel($i,2) [ expr $vel($i,2) / $veln ]
    set vel($i,3) [ expr $vel($i,3) / $veln ]

#
# Update the reference position, this is the crucial step in applying
# this technique.
#

# Obtain the projection of the current force along the new pulling
# direction

    set fv [ expr $fx*$vel($i,1) + $fy*$vel($i,2) + $fz*$vel($i,3) ]

# The new displacement vector is set in order to preserve this component
# of the force

    set xpar(1) [ expr $vel($i,1) * ( $ts * $vnorm($i) - $fv / $SMDk ) ]
    set xpar(2) [ expr $vel($i,2) * ( $ts * $vnorm($i) - $fv / $SMDk ) ]
    set xpar(3) [ expr $vel($i,3) * ( $ts * $vnorm($i) - $fv / $SMDk ) ]

# The new reference position is the current coordinate less than the
# desired displacement vector (this is a simplification that considers
# zero the orthogonal displacement, since the current position is used
# to compute the next direction, in order that it is always nule). 

    set ref($i,1) [ expr $x(1) - $xpar(1) ]
    set ref($i,2) [ expr $x(2) - $xpar(2) ]
    set ref($i,3) [ expr $x(3) - $xpar(3) ]
  }

# Printing total forces

  if { $print == 1 } {

# Closing the atom data file (opened before atom loop) 

    puts $file "END TIME STEP"
    close $file

# Opening and printing to the overall force file

    set file [ open $OutputName.total.dat a+ ]
    puts $file [ format "%8i %14.3f %14.3f %14.3f %14.3f %14.3f" \
                $ts $tforce(1) $tforce(2) $tforce(3) $fsum $twork ] 
    close $file
  }
  if { $restart == 1 } {
    print "Wrote restart information to $OutputName.restart at step $ts"
    close $restartfile
  }
}

#
# The end.
#
