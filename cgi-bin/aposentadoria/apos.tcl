#!/usr/bin/tclsh

set basehtml "http://leandro.iqm.unicamp.br/leandro/shtml/didatico/qf632/oscilantes"
puts "Content-type: text/html"
puts ""
puts "<head>"
puts "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">"
puts "<link rel=\"stylesheet\" type=\"text/css\" href=$basehtml/estilo.css>"
puts "<table width=60% align=center><tr><td>"

catch {
  set data [read stdin $env(CONTENT_LENGTH)]
  set data [ split $data "\n" ]
  set i 0
  set ivar 0
  set ival 0
  foreach arg $data {
    incr i
    set val($i) $arg
    if { [ string first "name=" $val($i) ] > 0 } {
      set ivar $i
    }
    if { $i == [ expr $ivar + 2 ] } {
      incr ival
      set value($ival) $val($i)
    }
  
  }
  
  set nanos $value(1)
  set juros $value(2)
  set inflacao $value(3)
  set salario $value(4)
  set contribuicao $value(5)
}

catch {
  exec ./apos $nanos $juros $inflacao $salario $contribuicao
} msg
set result [ split $msg "\n" ]
foreach line $result {
puts "$line <br>"
}

puts "</td></tr></table>"

