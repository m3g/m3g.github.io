#!/usr/bin/tclsh
#
# This is a script that writes SCORE and RMSD matrices from
# the output of a lovoalign all-on-all aligment
# 
# Must be run with:
#
# ./matrix.tcl result.dat
#
# where "result.dat" is an output file of a lovoalign all-on-all
# alignment (produced with a ./lovoalign -pdblist list.txt run).
#
# L. Martinez, Aug 25, 2009.
#

set i 0
foreach data [ split $argv " " ] {
  if { $data > " " } { incr i; set file($i) $data }
}
if { $i != 1 } { 
  puts " Run with: ./matrix.tcl input.dat > output.dat"
  puts "    where input.dat is an output file of a "
  puts "    \"lovoalign -pdblist\" type run. "
  exit
}

set score_type "TM-Score"

if { $score_type == "TM-Score" } {
  set maxscore 1.0
}
if { $score_type == "Structal" } {
  set maxscore 20.0
}

puts "# Input file: $file(1)"

set input [ open $file(1) r ]
set input [ read $input ]
set input [ split $input "\n" ]

set i 0
foreach line $input {
  if { [ string range $line 0 0 ] == "#" | $line <= " " } { continue }
  
  set iarg 0
  set line [ split $line " " ]
  foreach data $line { 
    if { $data <= " " } { continue }
    incr iarg
    set arg($iarg) $data
  }
  if { [ info exists found($arg(1)) ] == 0 } {
    incr i
    set pdb($i) $arg(1) 
    set found($pdb($i)) 1
    set index($pdb($i)) $i
  }
  if { [ info exists found($arg(2)) ] == 0 } {
    incr i
    set pdb($i) $arg(2) 
    set found($pdb($i)) 1
    set index($pdb($i)) $i
  }

  set score($index($arg(1)),$index($arg(2))) $arg(3)
  set rmsd($index($arg(1)),$index($arg(2))) $arg(5)

}
set nfiles $i

for { set i 1 } { $i <= $nfiles } { incr i } {
  puts "# File $i : $pdb($i) "
}

# Printing score matrix

puts "# SCORE MATRIX: "
for { set i 1 } { $i <= $nfiles } { incr i } {
  for { set j 1 } { $j <= $nfiles } { incr j } {
    if { $i == $j } {  
      puts -nonewline " [format %12.6f $maxscore ]" 
    } else {
      if { [ info exists score($i,$j) ] != 0 } {
        puts -nonewline " [format %12.6f $score($i,$j)]" 
      } else {
        puts -nonewline " [format %12.6f $score($j,$i)]" 
      }
    }
  }
  puts ""
}

# Printing rmsd matrix

puts "# RMSD MATRIX: "
for { set i 1 } { $i <= $nfiles } { incr i } {
  for { set j 1 } { $j <= $nfiles } { incr j } {
    if { $i == $j } {  
      puts -nonewline " [format %12.6f 0.0 ]" 
    } else {
      if { [ info exists rmsd($i,$j) ] != 0 } {
        puts -nonewline " [format %12.6f $rmsd($i,$j)]" 
      } else {
        puts -nonewline " [format %12.6f $rmsd($j,$i)]" 
      }
    }
  }
  puts ""
}






