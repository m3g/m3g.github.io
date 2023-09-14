#!/usr/bin/tclsh

puts "Content-type: text/html"
puts ""

# Start page

set url "http://leandro.iqm.unicamp.br/topolink"
set topolink_dir "../../../topolink"

# html header

set head [ open $topolink_dir/head.html r ] 
set head [ read $head ] 
set head [ split $head "\n" ]
foreach line $head { 
puts $line
}

# Procedure that reads an argument from a line of data
# at a specified position

proc read_data { line position } {
  set line [ split $line " " ]
  set i 0
  foreach string $line {
    if { [ string trim $string ] > " " } {
      incr i
      set data($i) $string
    }
  }
  return $data($position)
}

# Create temporary folder for this execution

set random [ expr rand() ]
set random [ split $random "." ]
foreach a $random { set tmp_folder $a }
set tmp_folder "$topolink_dir/serverfiles/TEMP$tmp_folder"
exec -ignorestderr mkdir $tmp_folder

# Read form data

set data [read stdin $env(CONTENT_LENGTH)]
set data [ split "$data" "\n" ]
set data [ split "$data" "&" ]
set i 0
foreach value $data { 
  incr i
  set value [ split "$value" "=" ]
  set j 0
  foreach val $value {
    incr j
    if { $j == 2 } { set option_val $val }
  }
  if { $i == 1 } { set protein_type $option_val } 
  if { $i == 2 } { set protein_size $option_val } 
  if { $i == 3 } { set linker_length $option_val } 
}

if { $protein_type == "alpha-beta" } { set type_name "All" }
if { $protein_type == "alpha" } { set type_name "Mostly alpha" }
if { $protein_type == "beta" } { set type_name "Mostly beta" }

#puts " Protein type: $type_name <br>"
#puts " Protein size: $protein_size <br>"
#puts " Linker length: $linker_length <br>"

catch {        
exec ./plot_by_linklength.py $protein_type $protein_size $linker_length $tmp_folder
} msg
puts $msg


catch {
  set returnpage [ open ./stat1-return.html "r" ] 
  set returnpage [ read $returnpage ]
  set returnpage [ split $returnpage "\n" ]
  foreach line $returnpage {
  
    if { $protein_type == "alpha-beta" } { 
      set line [ string map { "CHECKEDAB" "checked" } $line ]
      set line [ string map { "CHECKEDAONLY" " " } $line ]
      set line [ string map { "CHECKEDBONLY" " " } $line ]
    }
    if { $protein_type == "alpha" } { 
      set line [ string map { "CHECKEDAB" " " } $line ]
      set line [ string map { "CHECKEDAONLY" "checked" } $line ]
      set line [ string map { "CHECKEDBONLY" " " } $line ]
    }
    if { $protein_type == "beta" } { 
      set line [ string map { "CHECKEDAB" " " } $line ]
      set line [ string map { "CHECKEDAONLY" " " } $line ]
      set line [ string map { "CHECKEDBONLY" "checked" } $line ]
    }

    set map "PROTEIN_SIZE $protein_size"
    set line [ string map $map $line ]

    set map "LINKER_LENGTH $linker_length"
    set line [ string map $map $line ]

    set map "TMP_FOLDER $tmp_folder"
    set line [ string map $map $line ]
  
    puts $line
  }
} msg
puts $msg

set tail [ open $topolink_dir/tail.html r ] 
set tail [ read $tail ] 
set tail [ split $tail "\n" ]
foreach line $tail { 
puts $line
}


