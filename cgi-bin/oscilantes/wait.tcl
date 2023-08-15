#!/usr/bin/tclsh

set basehtml "http://leandro.iqm.unicamp.br/leandro/shtml/didatico/qf632/oscilantes"

puts "Content-type: text/html"
puts ""
puts "<head>"
puts "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">"
puts "<link rel=\"stylesheet\" type=\"text/css\" href=$basehtml/estilo.css>"
puts "<title>Reações Oscilantes</title>"
puts "</head>"
puts "<body>"
puts "<table width=80% align=center><tr><td align=center>"
puts "<br><h2>Calculando .</h2><br>"

catch {
after 500
puts -nonewline "." ; flush stdout
after 500
puts -nonewline "." ; flush stdout
after 500
puts -nonewline "." ; flush stdout
after 500
puts -nonewline "." ; flush stdout
after 500
puts -nonewline "." ; flush stdout
} msg
puts $msg
puts "</td></tr></table>"
puts "</body>"
