#!/usr/bin/tclsh
#
source ./procs.tcl
puts "[ header ]"

package require ncgi
::ncgi::parse

set project [ ::ncgi::value project ]
set myemail [ ::ncgi::value myemail ]
set add_user [ ::ncgi::value add_user ]

set temp [ open $archives/projects/projects.txt r ]
set projects [ read $temp ]
close $temp
set projects [ split $projects "\n" ]

set temp [ open $archives/projects/projects.txt w ]
foreach line $projects {
  if { [ string trim $line ] > " " } {
    array set data [ parse_line $line " " ] 
    if { $data(1) == $project } { 

# If the user is not in the project, add it

      if { [ string first " $add_user" $line ] == -1 } {
        set line "$line $add_user"

# If the user is in the project, remove it

      } else {
        set line [ string map [ list "$add_user" "" ] $line ]
      }
    }
    puts $temp [ string trim $line ] 
  }
}
close $temp

puts "
[ reload_form edit_project.tcl begin ]
[ hidden myemail $myemail ]
[ hidden project $project ]
[ reload_form edit_project.tcl end ]
"   

