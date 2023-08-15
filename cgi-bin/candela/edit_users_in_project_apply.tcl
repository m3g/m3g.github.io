#!/usr/bin/tclsh

source ./procs.tcl
puts "[ header ]"

package require ncgi
::ncgi::parse


set myemail [ ::ncgi::value myemail ]
set project [ ::ncgi::value project ]

# Read users file

eval [ read_users $archives ]

# Checking which authors were checked:

set new_authors " "
for { set i 1 } { $i <= $n_users } { incr i } { 
  set variable "author=$email($i)"
  if { [ ::ncgi::empty $variable ] != 1 } { 
    set new_authors "$new_authors $email($i)"
  }
}

# Updating the projects file in order to put the new author list

set temp [ open $archives/projects/projects.txt r ]
set projects [ read $temp ]
close $temp
set output [ open $archives/projects/projects.txt w ]
set projects [ split $projects "\n" ]
foreach line $projects {
  if { $line > " " } {
    set line [ split $line ";" ]
    set i 0
    foreach word $line {
      if { $word > " " & $word != ";" } { incr i; set data($i) $word }
    }
    if { $data(2) == "$project" } {
      set data(3) [ string trim $new_authors ]
      puts $output "$data(1);$data(2);$data(3)"
    } else {
      puts $output "$data(1);$data(2);$data(3)"
    }
  }
}

# Reloading project page

puts "
      [ reload_form edit_project.tcl begin ]
      [ hidden myemail $myemail ]
      [ hidden project $project ]
      [ reload_form edit_project.tcl end ]
     "
exit

