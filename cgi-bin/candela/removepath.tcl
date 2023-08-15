#!/usr/bin/tclsh

if { [ string last "\\" $file ] != -1 } {
  set file_name [ string range $file [ expr [ string last "\\" $file ] + 1 ] [ string length $file ] ] 
} else { 
  set file_name [ file tail $file ]
}


puts $file_name
