#!/usr/bin/tclsh

source ./procs.tcl
puts "[ header ]"

source ./errors.tcl ; catch { 

package require ncgi
::ncgi::parse
set try_again [ ::ncgi::exists try_again ]    

puts "
<html>
<table align=center width=70%>
<tr><td align=center colspan=3> <br> - Candela - <br><br> </td></tr>
<tr><td align=center colspan=3 bgcolor=lightblue> Version system for papers<br> </td><tr> 
</table>
<br><br>
<table align=center>
<tr><td bgcolor=lightgreen colspan=2> Registered user </td></tr>
[ form login.tcl ]
<tr><td> e-mail: </td><td><input type=\"text\" name=\"myemail\" size=30><br> </td></tr>
<tr><td> Password: </td><td> <input type=\"Password\" name=\"password\" size=30></td></tr>
<tr><td></td><td align=right>[ button login Login ]</td></tr>
</form>
"

if { $try_again } {
  set try_again [ ::ncgi::value try_again ]
  if { $try_again == 3 } {
    puts "<tr><td colspan=2 bgcolor=red> Invalid username or password. </td></tr>" 
  }
}

puts "</table>"

puts "
<table align=center>
<tr><td><br><br><br></td></tr>
<tr><td bgcolor=lightgrey colspan=2> New user </td></tr>
[ form add_user.tcl ]
<tr><td> Full name: (ex. Lionel A. Messi): </td><td><input type=\"text\" name=\"new_name\" size=30><br> </td></tr>
<tr><td> e-mail (ex. messi@barca.es): </td><td><input type=\"text\" name=\"new_email\" size=30><br> </td></tr>
<tr><td> Password: </td><td> <input type=\"Password\" name=\"new_password\" size=30></td></tr>
<tr><td></td><td align=right>[ button create Create ]</td></tr>
</form>
"

package require ncgi
::ncgi::parse
set try_again [ ::ncgi::exists try_again ]    

if { $try_again } {
  set try_again [ ::ncgi::value try_again ]
  if { $try_again == 1 } {
    puts "<tr><td colspan=2 bgcolor=red> Username already in use. </td></tr>" 
  }
  if { $try_again == 2 } {
    puts "<tr><td colspan=2 bgcolor=red> Use only letters and numbers. </td></tr>" 
  }
  if { $try_again == 4 } {
    puts "<tr><td colspan=2 bgcolor=red> Invalid username or password. </td></tr>" 
  }
}

puts " </table> </html> "

} status ; puts "[ errors $status ]"

