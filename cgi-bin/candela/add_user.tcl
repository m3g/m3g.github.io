#!/usr/bin/tclsh

source ./procs.tcl
puts "[ header ]"

package require ncgi
::ncgi::parse

set new_name [ ::ncgi::value new_name ]
set new_email [ ::ncgi::value new_email ]
set new_password [ ::ncgi::value new_password ]

# Check if the users file exists, if not, create one

if { [ file exists $archives/users ] != 1 } {
  exec mkdir $archives/users
}
if { [ file isfile $archives/users/users.txt ] != 1 } {
  exec touch $archives/users/users.txt
}

# Error for empty username or password

if { [ string trim $new_email ] <= " " | 
     [ string trim $new_name ] <= " "  | 
     [ string trim $new_password ] <= " " } {
    puts "[ reload_form start.tcl begin ]
          [ hidden try_again 4 ]
          [ reload_form start.tcl end ]"
    exit
}

# Error for spaces in the username or password

if { [ string first " " [ string trim $new_email ] ] != -1 | 
     [ string first " " [ string trim $new_password ] ] != -1 } {
    puts "[ reload_form start.tcl begin ]
          [ hidden try_again 2 ]
          [ reload_form start.tcl end ]"
    exit
}

# Error for email that is already registered

set temp [ open $archives/users/users.txt r ]
set user_file [ read $temp ]
close $temp
set user_file [ split $user_file "\n" ]
foreach line $user_file {
  array set data [ parse_line $line ";" ]
  if { $new_email  == $data(2) } {
    puts "[ reload_form start.tcl begin ]
          [ hidden try_again 1 ]
          [ reload_form start.tcl end ]"
    exit
  }
}

# No errors found, create user

set user_file [ open $archives/users/users.txt a ]
puts $user_file "$new_name;$new_email;[ encryptPassword $new_password ]"
close $user_file
puts "[ reload_form login.tcl begin ]
      [ hidden myemail $new_email ]
      [ hidden password $new_password ]
      [ reload_form login.tcl end ]"


