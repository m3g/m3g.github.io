#!/usr/bin/tclsh

set file [ open notas.dat r ]
set file [ read $file ]
set file [ split $file "\n" ]
set n_exame 0
set n_aprovados 0
foreach line $file {
  if { $line > " " } {
    set data [ split $line " " ]
    set i 0
    foreach dat $data {
      if { $dat > " " } {
        incr i
        set nota($i) $dat
      }
    }
    set media [ expr (2.0*$nota(1) + 2.0*$nota(2) + 3.0*$nota(3))/7.0 ]
    set media [ format %2.1f [ expr $media + 0.04999 ] ]
    if { $media >= 5.0 } { set sit "Aprovado" }
    if { $media < 5.0 & $media >= 3.0 } { set sit "Direito a exame" }
    if { $media < 3.0 } { set sit "Reprovado" }
    puts "$media  $sit"
  }
}
