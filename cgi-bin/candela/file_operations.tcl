#!/usr/bin/tclsh

source ./procs.tcl
puts "[ header ]"

package require ncgi
::ncgi::parse

set project [ ::ncgi::value project ]
set myemail [ ::ncgi::value myemail ]
set modify_file_name [ ::ncgi::value file_name ]
set action [ ::ncgi::value file_action ]
set modify_file_version [ ::ncgi::value file_version ]
set modify_file_description [ ::ncgi::value file_description ]
set file_error [ ::ncgi::value file_error ]
eval [ read_project_files $project $archives ]  

# Do nothing

if { $action == "none" } {
  puts "
        [ reload_form edit_project.tcl begin ]
        [ hidden myemail $myemail ]
        [ hidden project $project ]
        [ reload_form edit_project.tcl end ]
       "
}

# Mark file as locked

eval [ read_users $archives ]
if { $action == "Lock" } {
  for { set i_file 1 } { $i_file <= $n_files } { incr i_file } {
    if { $file_name($i_file) == $modify_file_name } {
      set file_data($i_file,status) "Locked by $name($myemail) <!-- =$myemail= -->"
    }
  }
  eval [ write_project_files $project $archives ]
  puts "
        [ reload_form edit_project.tcl begin ]
        [ hidden myemail $myemail ]
        [ hidden project $project ]
        [ reload_form edit_project.tcl end ]
       "
}

# Mark file as available to modify

if { $action == "Unlock" } {
  for { set i_file 1 } { $i_file <= $n_files } { incr i_file } {
    if { $file_name($i_file) == $modify_file_name } {
      set file_data($i_file,status) "Available to modify"
    }
  }
  eval [ write_project_files $project $archives ]
  puts "[ reload_form edit_project.tcl begin ]
        [ hidden myemail $myemail ]
        [ hidden project $project ]
        [ reload_form edit_project.tcl end ]
       "
}

# Deprecate file

if { $action == "Deprecate" } {
  for { set i_file 1 } { $i_file <= $n_files } { incr i_file } {
    if { $file_name($i_file) == $modify_file_name } {
      set mv_file $i_file
      break
    }
  }
  set new_name "$file_data($mv_file,version).$file_name($mv_file)"
  exec mv "$archives/projects/$project/files/$file_name($i_file)" $archives/projects/$project/deprecated/$new_name 
  for { set i $mv_file } { $i <= [ expr $n_files - 1 ] } { incr i } { 
    set file_name($i)             $file_name([ expr $i + 1 ])             
    set file_data($i,description) $file_data([ expr $i + 1 ],description) 
    set file_data($i,version)     $file_data([ expr $i + 1 ],version)     
    set file_data($i,status)      $file_data([ expr $i + 1 ],status)      
  }
  set n_files [ expr $n_files - 1 ] 
  eval [ write_project_files $project $archives ]
  puts "
        [ reload_form edit_project.tcl begin ]
        [ hidden myemail $myemail ]
        [ hidden project $project ]
        [ reload_form edit_project.tcl end ]
       "
}

# Upload new file version

if { $action == "Upload" } {

  puts "
  <html>
  <table align=center width=70%>
  <tr><td align=center colspan=3> <br> - Candela - <br><br> </td></tr>
  </table>"  

  puts "<table align=center width=70%>"
  puts "<tr><td bgcolor=lightblue> Project $project </td></tr>"
  puts "<tr><td></td></tr>"
  puts "<tr><td> Upload new version of file <b>$modify_file_name</b> </tr></td>"
 
  if { $file_error == "empty_file" } {
    puts "<tr><td bgcolor=red> <b>Error:</b> File empty or no file name. Try again. </tr></td>"
  }
  if { $file_error == "file_exists" } {
    puts "<tr><td bgcolor=red> <b>Error:</b> Uploaded file exists in project, but is not current file </tr></td>"
    puts "<tr><td bgcolor=red> You need to deprecate it before uploading another file with same name. </tr></td>"
  }

  puts "
  <tr><td colspan=2>
  [ form upload_new_version.tcl ]
  [ hidden myemail $myemail ]
  [ hidden project $project ]
  [ hidden old_file_name $modify_file_name ]
  [ hidden old_version $modify_file_version ]
  [ hidden old_description $modify_file_description ]
  File: <input type=\"file\" name=\"new_file_data\" class=\"browse\"><br>
  Description: 
  <input type=\"text\" name=\"new_file_desc\" size=40 
         value=\"$modify_file_description\">
  [ button Upload Upload ]
  </form>
  </td></tr>"

  set upload_error [ ::ncgi::exists upload_error ]
  if { $upload_error } {
    set upload_error [ ::ncgi::value upload_error ]  
    puts "<tr><td colspan=2 bgcolor=red> $upload_error </td></tr>"
  }                        
   
  puts "</table>"

  exit

}

# Put file on top of the list

if { $action == "Put_on_top" } {
  set i 2 
  for { set i_file 1 } { $i_file <= $n_files } { incr i_file } {
    if { $file_name($i_file) == $modify_file_name } {
      set temp_file_name(1)             $file_name($i_file)
      set temp_file_data(1,description) $file_data($i_file,description) 
      set temp_file_data(1,version)     $file_data($i_file,version)     
      set temp_file_data(1,status)      $file_data($i_file,status)      
    } else {
      set temp_file_name($i)             $file_name($i_file)
      set temp_file_data($i,description) $file_data($i_file,description) 
      set temp_file_data($i,version)     $file_data($i_file,version)     
      set temp_file_data($i,status)      $file_data($i_file,status)      
      incr i
    }
  }
  for { set i 1 } { $i <= $n_files } { incr i } { 
    set file_name($i)              $temp_file_name($i)            
    set file_data($i,description)  $temp_file_data($i,description)
    set file_data($i,version)      $temp_file_data($i,version)    
    set file_data($i,status)       $temp_file_data($i,status)     
  }
  eval [ write_project_files $project $archives ]
  puts "
        [ reload_form edit_project.tcl begin ]
        [ hidden myemail $myemail ]
        [ hidden project $project ]
        [ reload_form edit_project.tcl end ]
       "
  exit
}

# Move file up one position

if { $action == "Move_up" } { 
  for { set i_file 1 } { $i_file <= $n_files } { incr i_file } {
    if { $file_name($i_file) == $modify_file_name } {
      set i_swap $i_file
    }
  }
  if { $i_swap > 1 } {
    set i_up [ expr $i_swap - 1 ]
    set tmp_file_name  $file_name($i_up)
    set tmp_file_desc  $file_data($i_up,description)
    set tmp_file_vers  $file_data($i_up,version)
    set tmp_file_stat  $file_data($i_up,status)
    set file_name($i_up)             $file_name($i_swap)
    set file_data($i_up,description) $file_data($i_swap,description)
    set file_data($i_up,version)     $file_data($i_swap,version)
    set file_data($i_up,status)      $file_data($i_swap,status)
    set file_name($i_swap)             $tmp_file_name
    set file_data($i_swap,description) $tmp_file_desc
    set file_data($i_swap,version)     $tmp_file_vers
    set file_data($i_swap,status)      $tmp_file_stat
  }
  eval [ write_project_files $project $archives ]
  puts "
        [ reload_form edit_project.tcl begin ]
        [ hidden myemail $myemail ]
        [ hidden project $project ]
        [ reload_form edit_project.tcl end ]
       "
  exit
}

# Move file down one position

if { $action == "Move_down" } { 
  for { set i_file 1 } { $i_file <= $n_files } { incr i_file } {
    if { $file_name($i_file) == $modify_file_name } {
      set i_swap $i_file
    }
  }
  if { $i_swap < $n_files } {
    set i_down [ expr $i_swap + 1 ]
    set tmp_file_name  $file_name($i_down)
    set tmp_file_desc  $file_data($i_down,description)
    set tmp_file_vers  $file_data($i_down,version)
    set tmp_file_stat  $file_data($i_down,status)
    set file_name($i_down)             $file_name($i_swap)
    set file_data($i_down,description) $file_data($i_swap,description)
    set file_data($i_down,version)     $file_data($i_swap,version)
    set file_data($i_down,status)      $file_data($i_swap,status)
    set file_name($i_swap)             $tmp_file_name
    set file_data($i_swap,description) $tmp_file_desc
    set file_data($i_swap,version)     $tmp_file_vers
    set file_data($i_swap,status)      $tmp_file_stat
  }
  eval [ write_project_files $project $archives ]
  puts "
        [ reload_form edit_project.tcl begin ]
        [ hidden myemail $myemail ]
        [ hidden project $project ]
        [ reload_form edit_project.tcl end ]
       "
  exit
}



