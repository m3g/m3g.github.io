#!/usr/bin/tclsh

source ./procs.tcl
puts "[ header ]"

source ./errors.tcl ; catch { 

package require ncgi
::ncgi::parse
set project [ ::ncgi::value project ] 
set myemail [ ::ncgi::value myemail ]
set file_data [ ::ncgi::value file_data]
set file_desc [ ::ncgi::value file_desc ]
set file_name [ ::ncgi::importFile -client file_data ]

# Removing path from file name

if { [ string last "\\" $file_name ] != -1 } {
  set file_name [ string range $file_name [ expr [ string last "\\" $file_name ] + 1 ] [ string length $file_name ] ]
} else {
  set file_name [ file tail $file_name ]
}

# If the new file is empty, report and error and reload

if { [ string bytelength $file_data ] < 1 } {
puts "
  [ reload_form edit_project.tcl begin ]
  [ hidden myemail $myemail ]
  [ hidden project $project ]
  [ hidden upload_error {Empty file or empty file name. Try again.} ]
  [ reload_form edit_project.tcl end ]
"
exit
}

# Check if file exists already

set new_file $archives/projects/$project/files/$file_name
if { [ file exists $new_file ] } {
puts "
  [ reload_form edit_project.tcl begin ]
  [ hidden myemail $myemail ]
  [ hidden project $project ]
  [ hidden upload_error {File exists. Upload a new version instead or deprecate current before.} ]
  [ reload_form edit_project.tcl end ]
"
exit
}

# Check if a description was provided

if { [ string trim $file_desc ] <= " " } {
puts "
  [ reload_form edit_project.tcl begin ]
  [ hidden myemail $myemail ]
  [ hidden project $project ]
  [ hidden upload_error {Please add a brief description.} ]
  [ reload_form edit_project.tcl end ]
"
exit
}

# Set the version of this file

set date [ exec date +%d%h%y ] 
set i_version 0
set version "$date-$i_version"
while { [ file exists $archives/projects/$project/deprecated/$version.$file_name ] } {
  incr i_version
  set version "$date-$i_version"
}

# If everything is ok, upload the new file

set writefile [ open $new_file w ]
puts -nonewline $writefile $file_data
close $writefile

# Update file list

set filelist [ open $archives/projects/$project/files.txt a ]
puts $filelist "--file--"
puts $filelist "$file_name"
puts $filelist "$file_desc"
puts $filelist "$version"
puts $filelist "Available to modify"

puts "
[ reload_form edit_project.tcl begin ]
[ hidden myemail $myemail ]
[ hidden project $project ]
[ reload_form edit_project.tcl end ]
" 

} status ; puts "[ errors $status ]"






