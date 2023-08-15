#!/usr/bin/tclsh

source ./procs.tcl
puts "[ header ]"

package require ncgi
::ncgi::parse
set project [ ::ncgi::value project ] 
set myemail [ ::ncgi::value myemail ]
set new_file_data [ ::ncgi::value new_file_data]
set new_file_desc [ ::ncgi::value new_file_desc ]
set old_file_name [ ::ncgi::value old_file_name ]
set old_version [ ::ncgi::value old_version ]
set old_description [ ::ncgi::value old_description ]
set new_file_name [ ::ncgi::importFile -client new_file_data ]

# Removing path from file name

if { [ string last "\\" $new_file_name ] != -1 } {
  set new_file_name [ string range $new_file_name [ expr [ string last "\\" $new_file_name] + 1 ] [ string length $new_file_name ] ]
} else {
  set new_file_name [ file tail $new_file_name ]
}

# If the new file is empty, report and error and reload

if { [ string bytelength $new_file_data ] < 1 } {
puts "
[ reload_form file_operations.tcl begin ]
[ hidden myemail $myemail ]
[ hidden project $project ]
[ hidden file_action Upload ]
[ hidden file_name $old_file_name ]
[ hidden file_version $old_version ]
[ hidden file_description $old_description ]
[ hidden file_error "empty_file" ]
[ reload_form file_operations.tcl end ]
"
exit
}

# If the new file name exists but is not the current file, report an error and reload

if { $new_file_name != $old_file_name } {
  set test_file $archives/projects/$project/files/$new_file_name
  if { [ file exists $test_file ] } {
  puts "
  [ reload_form file_operations.tcl begin ]
  [ hidden myemail $myemail ]
  [ hidden project $project ]
  [ hidden file_action Upload ]
  [ hidden file_name $old_file_name ]
  [ hidden file_version $old_version ]
  [ hidden file_description $old_description ]
  [ hidden file_error "file_exists" ]
  [ reload_form file_operations.tcl end ]
  "
  exit
  }
}

# Deprecate current file

eval [ deprecate $project $old_file_name $old_version $archives ]

# Set the version of this file

set version [ version $new_file_name $project $archives ]

# If everything is ok, upload the new file

catch {
set writefile [ open $archives/projects/$project/files/$new_file_name w ]
puts -nonewline $writefile $new_file_data
close $writefile
} msg
puts $msg

# If the new description is empty, keep the old one

if { [ string trim $new_file_desc ] < " " } {
  set new_file_desc $old_description
}

# Update file list

eval [ read_project_files $project $archives ] 
for { set i 1 } { $i <= $n_files } { incr i } {
  if { $file_name($i) == $old_file_name } {
    set file_name($i) $new_file_name
    set file_data($i,description) $new_file_desc
    set file_data($i,version) $version
    set file_data($i,status) "Available to modify"
  }
}
eval [ write_project_files $project $archives ]

puts "
[ reload_form edit_project.tcl begin ]
[ hidden myemail $myemail ]
[ hidden project $project ]
[ reload_form edit_project.tcl end ]
" 




