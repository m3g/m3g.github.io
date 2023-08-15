#!/usr/bin/tclsh

source ./procs.tcl
puts "[ header ]"

source ./errors.tcl ; catch { 

package require ncgi
::ncgi::parse

set myemail [ ::ncgi::value myemail ]
set password [ ::ncgi::value password ]

if { [ file isfile $archives/users/users.txt ] != 1 } {

# If the users file does not exist (there are no users registered), return

  puts "[ reload_form start.tcl begin ]
        [ hidden try_again 3 ]
        [ reload_form start.tcl end ]"   

# Checking the username and password
  
} else {

  set users [ open $archives/users/users.txt r ]
  set users [ read $users ]
  set users [ split $users "\n" ]
  foreach line $users {
    array set data [ parse_line $line ";" ]
    if { $myemail == $data(2) & [ encryptPassword $password ] == $data(3) } {
      puts "[ reload_form edit_project.tcl begin ]
            [ hidden myemail $myemail ]
            [ reload_form edit_project.tcl end ]"
      exit
    } 
  }
  
# If the username or password are incorrect, return
  
  puts "[ reload_form start.tcl begin ]
        [ hidden try_again 3 ]
        [ reload_form start.tcl end ]"   
}
  
} status ; puts "[ errors $status ]"      

