#!/usr/bin/tclsh

source ./procs.tcl
puts "[ header ]"

package require ncgi
::ncgi::parse

set project [ ::ncgi::exists project ]
set myemail [ ::ncgi::value myemail ]

# Obtain the name of the user

eval [ read_users $archives ]
set myname $name($myemail)


# The title

puts "
<html>
<table align=center width=90%>
<tr><td align=center colspan=3> <br> - Candela - <br><br> </td></tr>
<tr><td align=right><a href=http://leandro.iqm.unicamp.br/candela>\[Logout\]</a></td></tr>
<tr><td align=center colspan=3 bgcolor=lightblue> Projects of $myname <br> </td></tr>
</table>"

# Check if any project is loaded

if { $project > 0 } { 
  set project [ ::ncgi::value project ]
  if { $project == "no_project" } { set project 0 }
}

# List of projects

puts "<table align=center width=90%><tr><td>"
if { [ file isfile $archives/projects/projects.txt ] != 1 } {
  set projects " "
} else {
  set temp [ open $archives/projects/projects.txt r ]
  set projects [ read $temp ]
  close $temp
}
set projects_list [ split $projects "\n" ] 
puts "[ form edit_user_project.tcl ]
      [ hidden myemail $myemail ]
      <select name=\"project\" class=\"sel\">"
      puts "<option bg=grey value=\"no_project\"> Select a project"

foreach line $projects_list { 
  if { [ string first $myemail $line ] != -1 } {
    array set data [ parse_line $line ";" ]
    if { $data(1) != "CLOSED" } {
      puts "<option bg=grey value=\"$data(2)\"> $data(2)"
    }
  }
}
puts "</select>
      [ hidden myemail $myemail ]
      [ button Load Load ]
      </form></td>"   

# Create new project 

puts "</td><td valign=top align=right>
   [ form new_project.tcl ]
   New project:
   <input type=\"text\" name=\"new_project\" size=10>
   [ hidden myemail $myemail ]
   [ button newproject {Create} ]
   </form>"
set new_project_error [ ::ncgi::exists new_project_error ] 
if { $new_project_error } {
  set new_project_error [ ::ncgi::value new_project_error ]
  puts "<tr><td colspan=2 bgcolor=red align=right> $new_project_error </td></tr>" 
}
puts "</td></tr></table>"


# This project

if { $project > 0 } {
  eval [ read_project_authors $project $archives ]  
  eval [ read_users $archives ]
  puts "
  <table align=center width=90%>
  <tr><td><br><br></td></tr>
  <tr><td align=center colspan=3 bgcolor=lightblue> Project <b>$project</b><br> </td><tr>
  </table>"   

# Output current list of authors
  
  puts "
  <table align=center width=90%>
  <tr><td width=10% valign=top> Authors: </td>
      <td valign=top>"
      for { set i 1 } { $i <= $n_authors } { incr i } {
        if { $i < $n_authors } { 
          puts "$name($author($i)), "
        } else {
          puts "$name($author($i))"
        }
      }

# Button for the edit authors form

puts "<td valign=top align=right>"
puts "[ form edit_users_in_project.tcl ]
      [ hidden myemail $myemail ]
      [ hidden project $project ]
      [ button Add "Edit authors" ]
      </form></td></tr></table>"   

# File list

  puts "
  <table align=center width=90%>
  <tr><td colspan=5><br><br></td></tr>
  <tr><td colspan=5 align=center><font color=red>Warning: always reload the project just before downloading 
         a file for modification to be sure it is unlocked! And then lock it!</font>
      [ form edit_user_project.tcl ] 
      [ hidden project $project ]
      [ hidden myemail $myemail ]
      [ button Load Reload\ Now! ]
      </form></td></tr>    
  <tr><td bgcolor=lightgreen width=10%> File </td>
  <td bgcolor=lightgreen width=30%> Name </td> 
  <td bgcolor=lightgreen width=10%> Version </td> 
  <td bgcolor=lightgreen width=25%> Status </td> 
  <td bgcolor=lightblue align=right width=35%> Actions </td>
  </tr>
  "
  
  set temp [ open $archives/projects/$project/files.txt ]
  set file [ read $temp ]
  close $temp
  set file [ split $file "\n" ]
  
  set i_file 0
  set color "white"
  
  eval [ read_project_files $project $archives ]

  for { set i_file 1 } { $i_file <= $n_files } { incr i_file } { 
  
    if { $color == "#E6E6E6" } { set color "white" } else { set color "#E6E6E6" } 
  
    puts "<tr>"

# A random number after the link url is necessary so that the browser always reload the file

   set random [ expr rand() ]

# Puts file_name as hiperlink to download file

    puts "<td bgcolor=$color valign=top> 
          <a target=\"_blank\" href=\"$archives/projects/$project/files/$file_name($i_file)\?$random\">
          $file_name($i_file)</td>
          </a></td>"

# File description, version and status

    puts "<td bgcolor=$color valign=top> $file_data($i_file,description)</td>"
    puts "<td bgcolor=$color valign=top> $file_data($i_file,version)</td> "
    puts "<td bgcolor=$color valign=top>"
    if { $file_data($i_file,status) == "Available to modify" } {
      puts "$file_data($i_file,status)"
    } else {
      puts "<font color=red>$file_data($i_file,status)</font>"
    }

# Check if file is available to be modified or locked

    set file_status [ string trim $file_data($i_file,status) ] 
    if { $file_status == "Available to modify" } {
      set lock_action "Lock"
      set locked_by_user "NoneXXX"
    } else {
      set lock_action "Unlock"
      set file_status [ split $file_status " " ]
      array set word [ parse_line $file_status "=" ]
      set locked_by_user $word(2)
    }

# Output file action options

    puts "</td>
          <td bgcolor=$color align=right>
          [ form file_operations.tcl ]
          [ hidden myemail $myemail ]
          [ hidden project $project ]
          [ hidden file_name $file_name($i_file) ]
          [ hidden file_version $file_data($i_file,version) ] 
          [ hidden file_description $file_data($i_file,description) ] 
          [ hidden file_error "none" ]
         <select name=\"file_action\" class=\"sel\" 
                 style=\"width:150px;background-color:white;\">
         <option value=\"none\"> Choose an action...
         <option value=\"$lock_action\"> $lock_action"
    if { $lock_action == "Lock" | $locked_by_user == $myemail } {
         puts "<option value=\"Upload\"> Upload new version"
         puts "<option value=\"Deprecate\"> Deprecate"
    }
    puts "<option value=\"Put_on_top\"> Put on top
          <option value=\"Move_up\"> Move one up
          <option value=\"Move_down\"> Move one down
          </select>
          [ button Go Go ]
          </form>
          </td></tr>"
  }


# Link to the comments page
#  puts "<tr><td colspan=5>"
#  puts "<tr><td><br><br><br></td><td></td></tr>"
#  puts "<tr><td colspan=5 bgcolor=lightgreen align=center> 
#<a href=$archives/projects/$project/comments.html onClick=\"return popup(this,'Help')\">
#Add or view comments.
#</a></td></tr>"


# Add new file to project
  
  puts "<tr><td colspan=5>"
  puts "<tr><td><br><br><br></td><td></td></tr>"
  puts "<tr><td colspan=5 bgcolor=lightgreen> Add new file to the project
        (to upload a new version of an existing file, choose the appropriate action above): </td></tr>"
  puts "
  <tr><td colspan=5>
  [ form upload.tcl ]
  [ hidden myemail $myemail ]
  [ hidden project $project ]
  File: <input type=\"file\" name=\"file_data\" class=\"browse\"><br>
  Description: <input type=\"text\" name=\"file_desc\" size=40 value=\"Figure 1\">
  [ button Upload Upload ]
  </form>
  </td></tr>"
  
  set upload_error [ ::ncgi::exists upload_error ]
  if { $upload_error } {
    set upload_error [ ::ncgi::value upload_error ]  
    puts "<tr><td colspan=5 bgcolor=red> $upload_error </td></tr>"
  } 

# List of deprecated project files
  
  puts "<tr><td><br><br><br></td><td></td></tr>"
  puts "<tr><td colspan=5 bgcolor=lightgrey> 
        <font size=-1>Deprecated project files:</font> </td></tr>"
  puts " "
  
  set deprecated [ exec ls -t $archives/projects/$project/deprecated/ ]
  set deprecated [ split $deprecated "\n" ]
  
  puts "<tr><td colspan=5 width=90%><font size=-1>"
  foreach file $deprecated {
    if { [ string trim $file ] > " " } { 
      puts "<a target=\"_blank\" href=\"$archives/projects/$project/deprecated/$file\"> $file </a> &nbsp; "
    }
  }
  puts "</td></tr></table>"

# If no project is selected

} else {

  puts "<table align=center width=90%><tr><td bgcolor=lightgreen> No project is selected </td></tr></table>"

}

# Project wide control

puts "<br><br><br><br>
      <table width=90% align=center><tr><td colspan=3><br></td></tr>
      <br><td colspan=3 bgcolor=lightgrey> Additional project-wide options </td></tr>"

# Mark project as closed

if { $project > 0 } {
  puts "<tr><td>
        [ form open_close.tcl ]
        [ hidden myemail $myemail ]
        [ hidden project $project ]
        [ hidden close close ]
        [ button open_close Mark\ project\ $project\ as\ CLOSED ]
        </form>
        </td>"

} else {
  puts "<tr><td></td>"
}

# Re-open closed project

puts "<td valign=top align=right> Re-open closed project: </td><td align=right width=30%> 
      [ form open_close.tcl ]
      [ hidden myemail $myemail ]
      [ hidden project $project ]
      <select name=\"open_project\" class=\"sel\">"
foreach line $projects_list {
  if { [ string first $myemail $line ] != -1 } {
    array set data [ parse_line $line ";" ]
    if { $data(1) == "CLOSED" } {
      puts "<option value=\"$data(2)\"> $data(2) "
    }
  }
}
puts "</select>
      [ button re_open Open ]
      </form></td></tr>"
puts "</table>"
puts "</html>"   


# User options

#
# Change password
#
#puts "<table width=90% align=center><tr><td colspan=3><br></td></tr>
#      <br><td colspan=4 bgcolor=lightgrey> User account options </td></tr>"
#
#puts "<tr>
#      [ form change_password ]
#      [ hidden myemail $myemail ]
#      [ hidden project $project ]
#      <td>Change password:</td>
#      <td>Current password:</td>
#      <td>New password:</td>
#      <td align=right>[ button change_password Change ]</td>
#      </tr>
#      "






