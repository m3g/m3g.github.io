<?php
// ==============
// Configuration
// ==============
$uploaddir = "./files"; 
// ==============
// Upload Part
// ==============
if(is_uploaded_file($_FILES['file1']['tmp_name']))
{
move_uploaded_file($_FILES['file1']['tmp_name'],$uploaddir.'/file1.pdb');
}
if(is_uploaded_file($_FILES['file2']['tmp_name']))
{
move_uploaded_file($_FILES['file2']['tmp_name'],$uploaddir.'/file2.pdb');
}
?>
<meta HTTP-EQUIV="REFRESH" content="0; url=../../cgi-bin/lovoalignonline/align.cgi"> 
