#!/usr/bin/tclsh

set files "/home/leandro/public_html/fractal/files"

puts "Content-type: text/html"
puts ""
puts "<table width=80% align=center><tr><td>"
puts "<br><h2>Results:</h2><br></tr></td>"
puts "<tr><td>"

set data [read stdin $env(CONTENT_LENGTH)]
set data [ split $data "&" ]
foreach line $data {
  if { [ string first "rulemin" $line ] != -1 } { 
    set mark [ string first "=" $line ]
    set mark [ expr 1 + $mark ]
    set rulemin [ string trim [ string range $line $mark 100 ] ]
  }
  if { [ string first "rulemax" $line ] != -1 } { 
    set mark [ string first "=" $line ]
    set mark [ expr 1 + $mark ]
    set rulemax [ string trim [ string range $line $mark 100 ] ]
  }
  if { [ string first "xdata" $line ] != -1 } { 
    set mark [ string first "=" $line ]
    set mark [ expr 1 + $mark ]
    set xcol [ string trim [ string range $line $mark 100 ] ]
  }
  if { [ string first "ydata" $line ] != -1 } { 
    set mark [ string first "=" $line ]
    set mark [ expr 1 + $mark ]
    set ycol [ string trim [ string range $line $mark 100 ] ]
  }
}
catch {
set check [ exec ./fractal $files/data2.dat 1 2 $rulemin $rulemax 30 ]
set check [ split $check "\n" ]
} msg

set i 0
foreach line $check {
  incr i 
  set n($i) [ string trim $line ]
}

puts " <h4>$n(2)</h4> <br><br>"

puts " <b>Method used:</b><br> The <i>correlation method</i> reported in <br>"
puts " J. Theiler, \"Estimating the Fractal Dimension\", <br>"
puts " <i> J. Opt. Soc. Am. A</i> Vol. 7, No. 6, pp. 1055-1073, 1990."

puts "<br><br>"

puts "<a href=http://leandro.iqm.unicamp.br/fractal> Click here start again.</a>"

puts "<br><br>"
puts "<h3> Important: </h3>"
puts "The fitting was done with all data. Please check if this is reasonable.<br>"
puts "If not, change the rulesize range or fit yourself the data below. <br><br><br>"
puts "Fitted data: <a href=http://leandro.iqm.unicamp.br/fractal/files/check.dat>check.dat</a><br><br>"
puts "If the graphs do not correspond to your data, reload the page.<br><br>"
puts "<b> $n(1) </b> <br>"
catch { exec gnuplot ./plot.gnp } msg
puts "<br></i>$msg</i><br>"
set date [ exec date ]
puts "<img src=\"http://leandro.iqm.unicamp.br/fractal/files/check.jpg?$date\" width=70% align=center>"

puts "<br><br>"
puts "This is the corresponding data: <br>"
puts "<tt>The values in x and y are normalized by the ranges"
puts "(x<sub>max</sub>-x<sub>min</sub> and y<sub>max</sub>-y<sub>min</sub>)"
puts "in order to minimize dimensionality problems.</tt>"
catch { exec gnuplot ./plot2.gnp } msg
puts "<br><i>$msg</i><br>"
puts "<img src=\"http://leandro.iqm.unicamp.br/fractal/files/data.jpg?$date\" width=70% align=center>"
puts "</td></tr></table>"
puts "</body>"
