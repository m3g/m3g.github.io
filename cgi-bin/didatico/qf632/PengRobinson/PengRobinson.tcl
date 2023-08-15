#!/usr/bin/tclsh

set pagepath="http://leandro.iqm.unicamp.br/leandro/shtml/didatico/qf632/PengRobinson"

puts "Content-type: text/html"
puts ""

puts "<head>"
puts "<meta charset=utf-8>"
puts "<link rel=\"stylesheet\" type=\"text/css\" href=$pagepath/estilo.css>"
puts "<title>Agriprec</title>"
puts "</head>"

puts "<body>"

puts "<table width=70% align=center><tr><td>"

set files "./files"

catch {
  set data [read stdin $env(CONTENT_LENGTH)]
  set data [ split $data "\n" ]
puts $data
exit
  set file [ open $files/data.dat w ]
  foreach line $data { 
    if { [ string range $line 0 3 ] == "Cont" ||
         [ string range $line 0 3 ] == "----" ||
         [ string range [ string trim $line ] 0 0 ] <= " " } { continue }
    puts $file $line 
  }
}
close $file

catch {
  exec ./run.sh
} error
puts $error

catch {
  puts "<br><h2>Resultados:</h2><br>"
  set result [ open $files/result.dat r ] 
  set result [ read $result ]
  set result [ split $result "\n" ]
  foreach line $result { 
    if { [ string range [ string trim $line ] 0 0 ] != "#"  } { 
      if { [ string length [ string trim $line ] ] == 0 } { 
        puts "<br>" 
      } else {
        puts [ string range $line 1 [ string length $line ] ] 
      }
    }
  } 
} error
puts $error

puts "<br><br><br><br>"
puts "<br><br><br><br>"
puts "<br><br><br><br>"
puts "<br><br><br><br>"
puts "<br><br><br><br>"

puts "</td></tr></table>"
puts "</body>"
