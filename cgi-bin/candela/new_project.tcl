#!/usr/bin/tclsh

source ./procs.tcl
puts "[ header ]"

package require ncgi
::ncgi::parse

set new_project [ ::ncgi::value new_project ]
set myemail [ ::ncgi::value myemail ]

catch {

# Error for empty project name

if { [ string trim $new_project ] <= " " } {
    puts "[ reload_form edit_project.tcl begin ]
          [ hidden myemail $myemail ]
          [ hidden new_project_error {Invalid project name} ]
          [ reload_form edit_project.tcl end ]"
    exit
}

# Error for spaces or semi-colons in project name

if { [ string first " " [ string trim $new_project ] ] != -1 |
     [ string first "@" [ string trim $new_project ] ] != -1 |
     [ string first "," [ string trim $new_project ] ] != -1 |
     [ string first ";" [ string trim $new_project ] ] != -1 } {
    puts "[ reload_form edit_project.tcl begin ]
          [ hidden myemail $myemail ]
          [ hidden new_project_error {Use only standard characters.} ]
          [ reload_form edit_project.tcl end ]"
    exit
}

# Create project file and directory if they don't exist already

if { [ file exists $archives/projects ] != 1 } {
  exec mkdir $archives/projects
}
if { [ file isfile $archives/projects/projects.txt ] != 1 } {
  exec touch $archives/projects/projects.txt
}

# Error for project name that already exists

set temp [ open $archives/projects/projects.txt r ]
set project_file [ read $temp ]
close $temp
set project_file [ split $project_file "\n" ]
foreach line $project_file {
  array set name [ parse_line $line ";" ]
  if { $new_project == $name(2) } {
    puts "[ reload_form edit_project.tcl begin ]
          [ hidden myemail $myemail ]
          [ hidden new_project_error {Project name already exists.} ]
          [ reload_form edit_project.tcl end ]"
    exit
  }
}

# No errors found, create new project

set project_file [ open $archives/projects/projects.txt a ]
exec mkdir $archives/projects/$new_project
exec mkdir $archives/projects/$new_project/files
exec mkdir $archives/projects/$new_project/deprecated
exec touch $archives/projects/$new_project/files.txt

puts $project_file "OPEN;$new_project;$myemail"
close $project_file
puts "[ reload_form edit_project.tcl begin ]
      [ hidden myemail $myemail ]
      [ hidden project $new_project ]
      [ reload_form edit_project.tcl end ]"

} msg
puts $msg

