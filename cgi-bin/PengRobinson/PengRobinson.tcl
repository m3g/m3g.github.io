#!/usr/bin/tclsh

set basehtml "http://leandro.iqm.unicamp.br/leandro/html/didatico/qf632/PengRobinson"
set filehtml "$basehtml/files"
set files "/home/leandro/public_html/leandro/html/didatico/qf632/PengRobinson/files"

puts "Content-type: text/html"
puts ""
puts "<head>"
puts "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">"
puts "<link rel=\"stylesheet\" type=\"text/css\" href=$basehtml/estilo.css>"
puts "<title>Peng Robinson</title>"
puts "</head>"
puts "<body>"
puts "<table width=80% align=center><tr><td align=center>"
puts "<br><h2>Resultado:</h2><br>"

catch {
set data [read stdin $env(CONTENT_LENGTH)]
set data [ split $data "\n" ]
set i 0
set ivar 0
set ival 0
foreach arg $data {
  incr i
  set val($i) $arg
  set found_name [ string first "name=" $val($i) ]
  if { $found_name > 0 } {
    incr ival
    set ivar $i
    set iname [ expr $found_name + 6 ]
    set lname [ string first \" [ string range $val($i) $iname [ string length $val($i) ] ] ]
    set lname [ expr $lname - 1 ]
    set name($ival) [ string range $val($i) $iname [ expr $iname + $lname ] ]
  }
  if { $i == [ expr $ivar + 2 ] } {
    #incr ival
    set value($ival) $val($i)
    set $name($ival) $value($ival)
  }
}
set nval $ival


catch { 
exec -ignorestderr ./PengRobinson $tc $pc $w $acp $bcp $ccp $dcp $tr $pr $te $p $opcao
} msg
#puts "<table align=center>"
#set i 1
#foreach data $msg {
##puts "$i $data <br>"
#  if { $i == 1 } { puts "<tr><td> $data =" }
#  if { $i == 3 } { puts "$data </td></tr>" }
#  incr i
#  if { $i == 4 } { set i 1 }
#}
#puts "</table>"
puts $msg

puts "<br><br><a href=http://leandro.iqm.unicamp.br/leandro/html/didatico/qf632/PengRobinson>\[voltar\]</a>"

puts "<br><br><b>Dados fornecidos:</b><br><br>"
for { set ival 1 } { $ival <= $nval } { incr ival } {
puts "$name($ival) = $value($ival) <br>"
}

} msg
puts $msg
puts "</td></tr></table>"
puts "</body>"
