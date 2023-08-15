#!/usr/bin/tclsh

source ./procs.tcl
puts "[ header ]"

package require ncgi
::ncgi::parse

set project [ ::ncgi::exists project ]
set myemail [ ::ncgi::value myemail ]

# The title

# Read user file
eval [ read_users $archives ]

puts "
<html>
<table align=center width=90%>
<tr><td align=center colspan=3> <br> - Candela - <br><br> </td></tr>
<tr><td align=right><a href=http://leandro.iqm.unicamp.br/candela>\[Logout\]</a></td></tr>
<tr><td align=center colspan=3 bgcolor=lightblue> User: $name($myemail) <br> </td></tr>
</table>"

# This project

set project [ ::ncgi::value project ]
puts "
<table align=center width=90%>
<tr><td><br><br></td></tr>
<tr><td align=center colspan=3 bgcolor=lightblue> Editing authors of project: <b>$project</b><br> </td><tr>
</table>"   
  
# Read project file to read which are the authors of this project

eval [ read_project_authors $project $archives ]

# Show current list of authors 

puts "
<table align=center width=90%>
<tr><td width=20% valign=top> Current authors: </td>
    <td valign=top>"
for { set i 1 } { $i <= $n_authors } { incr i } {
  if { $i < $n_authors } {
    puts "$name($author($i)),"
  } else {
    puts "$name($author($i))"
  }
}

# Display the list of all users with checkboxes: 

puts "</td></tr>"
puts "<tr><td width=20% valign=top> &nbsp; </td><td>"
puts "<tr><td width=20% valign=top> Select authors: </td><td>"

puts "<form action=./edit_users_in_project_apply.tcl method=\"get\"><ul class=\"checklist\">"
puts "[ hidden myemail $myemail ]
      [ hidden project $project ]"
set checked " "
for { set i 1 } { $i <= $n_users } { incr i } {
  set checked " "
  for { set j 1 } { $j <= $n_authors } { incr j } {
    if { $email($i) == $author($j) } { set checked "checked" } 
  }
  puts "<li><label for=\"teste:$email($i)\"><input id=\"$email($i)\" name=\"author=$email($i)\" 
                   type=\"checkbox\" $checked/>"
  puts "$name($email($i)) ($email($i))"
  puts "</label></li>"
}

# Warn the user about unchecking himself:

puts "<tr><td bgcolor=white colspan=3> &nbsp; </td></tr>"  
puts "<tr><td bgcolor=red colspan=3> Important: </td></tr>"  
puts "<tr><td bgcolor=white colspan=3> If you uncheck yourself, you will be removed from the project."
puts " To be able again to edit the project, other user will have to add you. </td></tr>"  
puts "</tr></td>"
puts "<tr><td bgcolor=white colspan=3 align=right>"

# Submit changes:

puts "[ button update_authors {Update Authors} ]"
puts "</ul></form>"                  

puts "</tr></td>"
puts "</table>"

puts "</html>"   
