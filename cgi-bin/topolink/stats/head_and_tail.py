#!/usr/bin/python

url="http://leandro.iqm.unicamp.br/topolink" 
source_files="../../../topolink"

def htmlTop() :
  print "Content-type: text/html"
  print ""

def pageHead() :
  file=open(source_files+"/head.shtml","r")
  for line in file :
    if "title.shtml" in line :
      title_file=open(source_files+"/title.shtml","r")
      title_file=title_file.read()
      print title_file
    if "menu.shtml" in line :
      menu_file=open(source_files+"/menu.shtml","r")
      menu_file=menu_file.read()
      print menu_file
    print line

def pageTail() :
  file=open(source_files+"/tail.shtml")
  file_html=file.read()
  print file_html


