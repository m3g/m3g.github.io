#
# Procedure to report errors to log file. The log file
# will be available at /cgi-bin/candela/log/errors.log
#

proc errors { msg } {
  if { $msg > " " } {
    set log [ open $archives/log/errors.log a ]
    puts $log "[ exec date ]"
    puts $log $msg
    return " <br><br> Errors were found. Please warn the administrator. "
  } else {
    return ""
  }
}


