#!/usr/bin/tclsh
#
# password.tcl
#

package require des
package require base64

# Procedure to "encrypt" passwords
#
proc encryptPassword { password } {
  set encryptedPassword "00000000$password"
  set encryptedPassword \
    [ string range $encryptedPassword \
    [ expr [ string length $encryptedPassword ] - 8 ] \
    [ string length $encryptedPassword ] ]
  set encryptedPassword "[ base64::encode -wrapchar {} [DES::des -dir encrypt -key "abcdefgh" "$encryptedPassword" ] ]"
  return $encryptedPassword
}

# Packare required for encryption

puts " New password: "
gets stdin password
set new_password [ encryptPassword $password ]
puts $new_password









