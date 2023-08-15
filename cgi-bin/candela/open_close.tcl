#!/usr/bin/tclsh

source ./procs.tcl
puts "[ header ]"

package require ncgi
::ncgi::parse

set project [ ::ncgi::value project ]
set myemail [ ::ncgi::value myemail ]

# Close project

set close [ ::ncgi::exists close ]
if { $close } {
  set temp [ open $archives/projects/projects.txt r ]
  set project_file [ read $temp ]
  close $temp
  set project_file [ split $project_file "\n" ]
  set write_projects [ open $archives/projects/projects.txt w ]
  foreach line $project_file {
    if { [ string trim $line ] > " " } {
      array set data [ parse_line $line ";" ]
      if { $project == $data(2) } {
        puts $write_projects "CLOSED;$data(2);$data(3)"
      } else {
        puts $write_projects [ string trim $line ]
      }
    }   
  }
  close $write_projects
  puts "[ reload_form edit_project.tcl begin ]
        [ hidden myemail $myemail ]
        [ reload_form edit_project.tcl end ]"
  exit
}

# Reopen

set open_project [ ::ncgi::value open_project ]
set temp [ open $archives/projects/projects.txt r ]
set project_file [ read $temp ]
close $temp
set project_file [ split $project_file "\n" ]
set write_projects [ open $archives/projects/projects.txt w ]
foreach line $project_file {
  if { [ string trim $line ] > " " } {
    array set data [ parse_line $line ";" ]
    if { $open_project == $data(2) } {
      puts $write_projects "OPEN;$data(2);$data(3)"
    } else {
      puts $write_projects $line
    }
  }
}   
close $write_projects
puts "[ reload_form edit_project.tcl begin ]
      [ hidden myemail $myemail ]
      [ hidden project $open_project ]
      [ reload_form edit_project.tcl end ]"
exit

