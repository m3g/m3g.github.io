#!/usr/bin/tclsh

set files "/home/leandro/public_html/fractal/files"

puts "Content-type: text/html"
puts ""
puts "<table width=80% align=center><tr><td>"
puts "<br><h2>Properties of the measures:</h2><br>"

set data [read stdin $env(CONTENT_LENGTH)]
set data [ split $data "&" ]
foreach line $data {
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
set check [ exec ./checksecond $files/data.dat $xcol $ycol ]
set check [ split $check "\n" ]
} msg

set i 0
foreach line $check {
  incr i 
  set n($i) [ string trim $line ]
}

puts " <h4>Choose range of ruler sizes:</h4> <br>"
puts " Suggested values are 5 times minimum distance between points) and  "
puts " 1/5 of the maximum distance. <br><br> "

catch {
puts "
<form action=\"/cgi-bin/fractal/fractal3.cgi\" method=\"post\" ENCTYPE=\"multipar
t/form-data\">
<!-- 
<form action=\"/cgi-bin/fractal/upload.php\" method=\"post\" ENCTYPE=\"multipart/form-data\">
-->
<table>
<tr><td>Minimum ruler size: </td><td><input type=\"text\" name=\"rulemin\" size=\"20\" value=\"[ expr $n(1) * 5 ]\"></td>
<tr><td>Maximum ruler size: </td><td><input type=\"text\" name=\"rulemax\" size=\"20\" value=\"[ expr $n(2) / 5. ]\"></td><br>
<td align=right>&nbsp;&nbsp;<input type=\"submit\" value=\"Continue\"></td></tr>
</table><br><br>
Minimum distance between points: $n(1) <br>
Maximum distance between points: $n(2) <br>
<br><br>
<table>
<tr><td><b>Previous options (do not change here): </b></td></tr>
<tr><td>X column: </td><td><input type=\"text\" name=\"xdata\" size=\"5\" value=\"$xcol\"></td></tr>
<tr><td>Y column: </td><td><input type=\"text\" name=\"ydata\" size=\"5\" value=\"$ycol\"></td></tr>
</table>
</form>"

puts "</td></tr></table>"
puts "</body>"
} msg
puts $msg
