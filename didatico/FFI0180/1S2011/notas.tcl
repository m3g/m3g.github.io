#!/usr/bin/tclsh

proc parse_line { line char } {
  set line [ split $line $char ]
  set i 0
  foreach word $line { 
    if { $word > " " & $word != $char } { incr i; set words($i) $word }
  }
  array get words 
}

set file [ open notas.dat r ]
set notas [ read $file ]
set notas [ split $notas "\n" ]

set iline 0
foreach line $notas {

  if { $line > " " } {
 
  incr iline

  array set data [ parse_line $line " " ] 
  
  set n_notas 0
  set media_provinha 0.
  set media_relatorio 0.

  set pior_nota 10.
  set ipior 0
  for { set i 1 } { $i <= 11 } { set i [ expr $i + 2 ] } {
    if { $data($i) == "F" | $data($i) == "0.0" | $data($i) == "-" } {
      set pior_nota 0.0
      set npior $i
    } elseif { $data($i) > "0" & $data($i) != "-" } {
      set media_pratica [ expr 0.3*$data($i) + 0.7*$data([expr $i + 1]) ]
      if { $media_pratica < $pior_nota } { 
        set pior_nota $media_pratica 
        set npior $i
      }
    }
  } 

  set n_notas 0
  set media_provinha 0.
  set media_relatorio 0.
  for { set i 1 } { $i <= 11 } { set i [ expr $i + 2 ] } {
    if { $i != $npior } {
      if { $data($i) == "F" | $data($i) == "0.0" } {
        incr n_notas
      } elseif { $data($i) > "0" & $data($i) != "-" } {
        incr n_notas
        set media_provinha [ expr $media_provinha + $data($i) ]
        set media_relatorio [ expr $media_relatorio + $data([expr $i + 1]) ]
      }
    }
  } 
  set media_provinha [ expr $media_provinha / $n_notas ] 
  set media_relatorio [ expr $media_relatorio / $n_notas ] 
  set media [ expr 0.7*$media_relatorio + 0.3*$media_provinha ]
  puts "$iline $n_notas [ format "%2.2f %2.2f %2.2f"  $media_provinha $media_relatorio $media ]"


  }

}




