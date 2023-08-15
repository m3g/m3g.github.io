#!/usr/bin/tclsh

set files "/home/leandro/public_html/fractal/files"

puts "Content-type: text/html"
puts ""
puts "<table width=80% align=center><tr><td>"
puts "<br><h2>Properties of the data file:</h2><br>"
puts "<h3>Check before proceeding.</h3><br><br>"

catch {
set data [read stdin $env(CONTENT_LENGTH)]
set data [ split $data "\n" ]
set file [ open $files/data.dat w ]
foreach line $data { puts $file $line }
close $file
set check [ exec ./checkfirst $files/data.dat ]
set check [ split $check "\n" ]

set i 0
foreach line $check {
  incr i 
  set n($i) [ string trim $line ]
} 

puts " Number of data lines found: <b>$n(1)</b> <br><br>"
puts " Number of data columns found: <b>$n(2)</b> <br>"

if { $n(1) < 4 } { puts " ERROR: Found less than four lines of data. "; exit }
if { $n(2) < 2 } { puts " ERROR: Found less than two columns of data. "; exit }

puts "
<form action=\"/cgi-bin/fractal/fractal2.cgi\" method=\"post\" ENCTYPE=\"multipar
t/form-data\">
<!-- 
<form action=\"/cgi-bin/fractal/upload.php\" method=\"post\" ENCTYPE=\"multipart/form-data\">
-->
<table>
<tr><td>Column containing X data: </td><td><input type=\"text\" name=\"xdata\" size=\"5\" value=\"1\"></td>
<tr><td>Column containing Y data: </td><td><input type=\"text\" name=\"ydata\" size=\"5\" value=\"2\"></td><br>
<td align=right>&nbsp;&nbsp;<input type=\"submit\" value=\"Continue\"></td></tr>
</form>"

} msg
puts $msg

puts "</td></tr></table>"
puts "</body>"
