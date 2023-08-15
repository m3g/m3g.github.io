#!/usr/bin/tclsh
#
# procs.tcl
#

# Packare required for encryption

package require des
package require base64

# The path to archives

set archives "../../candela"

# The header for all pages

proc header { } {
  return {Content-type: text/html

<style type="text/css"><!--
a:visited {
	text-decoration: none;
	font-size: 12;
	font-family: arial, helvetica, sans-serif;
}

a:link {
	text-decoration: none;
	font-size: 12;
	font-family: arial, helvetica, sans-serif;
}
table {
	border: flat;
}

td {
	font-size: 12;
	font-family: arial, helvetica, sans-serif;
	border: none;
}
body {
	font-size: 12;
	font-family: arial, helvetica, sans-serif;
}
input.btn { 
      color:black; 
      font-size: 12;
      font-family: arial, helvetica, sans-serif;
      background-color: white;
      border: 1px solid; 
      border-color: black black black black; 
} 
select.sel { 
      color:black; 
      font-size: 12;
      font-family: arial, helvetica, sans-serif;
      background-color: transparent;
      border: 1px solid; 
      border-color: black black black black; 
} 
file.browse { 
      color:black; 
      font-size: 12;
      font-family: arial, helvetica, sans-serif;
      background-color: transparent;
      border: 1px solid; 
      border-color: black black black black; 
} 
.checklist {
    border: 1px solid #ccc;
    list-style: none;
    height: 20em;
    overflow: auto;
    width: 400px;
    align: left;
}
.checklist, .checklist li { margin: 0; padding: 0; }
.checklist label {
    display: lock;
    padding-left: 25px;
    text-indent: -25px;
}
.checklist label:hover { background: #777; color: #fff; }
--></style>

<SCRIPT TYPE="text/javascript">
<!--
function popup(mylink, windowname)
{
if (! window.focus)return true;
var href;
if (typeof(mylink) == 'string')
   href=mylink;
else
   href=mylink.href;
window.open(href, windowname, 'width=650,height=600,scrollbars=yes');
return false;
}
function initChecklist() {
    if (document.all && document.getElementById) {
        // Get all unordered lists
         var lists = document.getElementsByTagName("ul");
        
        for (i = 0; i < lists.length; i++) {
            var theList = lists[i];
            
              // Only work with those having the class "checklist"
            if (theList.className.indexOf("checklist") > -1) {
                var labels = theList.getElementsByTagName("label");
                
                   // Assign event handlers to labels within
                for (var j = 0; j < labels.length; j++) {
                    var theLabel = labels[j];
                    theLabel.onmouseover = function() { this.className += " hover"; };
                    theLabel.onmouseout = function() { this.className = this.className.replace(" hover", ""); };
                }
            }
        }
    }
}
//-->
</SCRIPT>


  }
}

# Separate a line in words(1) and words(2) ...

proc parse_line { line char } {
  set line [ split $line $char ]
  set i 0
  foreach word $line { 
    if { $word > " " & $word != $char } { incr i; set words($i) $word }
  }
  array get words 
}

# Count the number of words

proc count_words { line char } {
  set line [ split $line $char ]
  set i 0
  foreach word $line { 
    if { $word > " " & $word != $char } { incr i }
  }
  return $i
}

# form title and end

proc form { script } {
  return "<form action=\"/cgi-bin/candela/$script\" \
          method=\"POST\" enctype=\"multipart/form-data\">"
}

# Hidden forms with problem data

proc hidden { name value } {
  return "<input type=\"hidden\" name=\"$name\" value=\"$value\">"
}

# Simplify button syntax

proc button { name value } {
  return "<input name=\"$name\" type=\"submit\" value=\"$value\" class=\"btn\">" 
}

# Procedure that reads the file containg files of project

proc read_project_files { project archives } {
  global n_files file_name file_data

  set temp [ open $archives/projects/$project/files.txt r ] 
  set files [ read $temp ]
  close $temp
  set files [ split $files "\n" ]

  set n_files 0
  foreach line $files {
    if { [ string trim $line ] == "--file--" } { incr n_files ; set i 0 }
    set data($n_files,$i) $line
    incr i
  }

  for { set i_file 1 } { $i_file <= $n_files } { incr i_file } {
    set file_name($i_file) $data($i_file,1)
    set file_data($i_file,description) $data($i_file,2)
    set file_data($i_file,version) $data($i_file,3)
    set file_data($i_file,status) $data($i_file,4)
  }
}

# Procedure that writes the file containing files of project

proc write_project_files { project archives } {
  global n_files file_name file_data
 
  set output [ open $archives/projects/$project/files.txt w ] 

  for { set i_file 1 } { $i_file <= $n_files } { incr i_file } {
    puts $output "--file--"  
    puts $output "$file_name($i_file)"
    puts $output "$file_data($i_file,description)"
    puts $output "$file_data($i_file,version)"
    puts $output "$file_data($i_file,status)"    
  }

  close $output

}

# Reload some page sending some information in a form

proc reload_form { page action } {
  if { $action == "begin" } {
    return "
            <head></head><body>
            <form name=\"reload_page\" method=\"POST\" action=\"/cgi-bin/candela/$page\">
            "
  }
  
  if { $action == "end" } {
    return "</form>
            <script type=\"text/javascript\" language=\"JavaScript\"><!-- 
            document.reload_page.submit(); 
            //--></script>
            </body>
            "              
  }
}

# The current files directory of a project

proc files_dir { project archives } {
  return "$archives/projects/$project/files"
}

# The deprecated files directory of a project

proc deprecated_dir { project archives } {
  return "$archives/projects/$project/deprecated"
}

# Procedure to deprecate file
                                      
proc deprecate { project file version archives } {
  if { [ string trim $file ] > " " } {
    catch { 
    exec mv [ files_dir $project $archives ]/$file [ deprecated_dir $project $archives ]/$version.$file
    } msg
    puts $msg
  }
  return
}

# Procedure that sets the version of a newly uploaded file

proc version { file project archives } { 
  set date [ exec date +%d%h%y ] 
  set i_version 0
  set version "$date-$i_version"
  while { [ file exists [ deprecated_dir $project $archives ]/$version.$file ] } {
    incr i_version
    set version "$date-$i_version"
  }
  return $version
}
      
# Procedure to read the users file

proc read_users { archives } {
global n_users name email
  set temp [ open $archives/users/users.txt r ]
  set users [ read $temp ]
  close $temp
  set users [ split $users "\n" ]
  set users [ lsort -dictionary $users ]
  set i_user 0 
  foreach line $users {
    if { $line > " " } {
      incr i_user
      set line [ split $line ";" ]
      set i 0
      foreach word $line { 
        if { $word > " " & $word != ";" } { incr i; set data($i) $word }
      }
      set email($i_user) $data(2)
      set name($data(2)) $data(1) 
    }
  }
  set n_users $i_user
  return
}

# Procedure to read the authors of a project

proc read_project_authors { project archives } {
  global n_authors author
  set temp [ open $archives/projects/projects.txt r ] 
  set projects [ read $temp ]
  close $temp
  set projects [ split $projects "\n" ]
  foreach line $projects {
    set line [ split $line ";" ]
    set i 0
    foreach word $line { 
      if { $word > " " & $word != ";" } { incr i; set data($i) $word }
    }
    if { $data(2) == "$project" } { 
      set authors_temp [ split $data(3) " " ]
      set n_authors 0
      foreach email $authors_temp {
        incr n_authors
        set author($n_authors) $email
      }
      return
    }
  }
}

# Procedure to "encrypt" passwords
#
proc encryptPassword { password } {
  set encryptedPassword "00000000$password"
  set encryptedPassword [ string map { - % } "$encryptedPassword" ]
  set encryptedPassword \
    [ string range $encryptedPassword \
    [ expr [ string length $encryptedPassword ] - 8 ] \
    [ string length $encryptedPassword ] ]
  set encryptedPassword "[ base64::encode -wrapchar {} [DES::des -dir encrypt -key "abcdefgh" "$encryptedPassword" ] ]"
  return $encryptedPassword
}









