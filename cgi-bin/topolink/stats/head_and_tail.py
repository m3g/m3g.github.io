#!/usr/bin/python

url="https://m3g.github.io/topolink" 
source_files="../../../topolink"

def htmlTop() :
  print "Content-type: text/html"
  print ""

def pageHead() :
  file=open(source_files+"/head.html","r")
  for line in file :
    if "title.html" in line :
      title_file=open(source_files+"/title.html","r")
      title_file=title_file.read()
      print title_file
    if "menu.html" in line :
      menu_file=open(source_files+"/menu.html","r")
      menu_file=menu_file.read()
      print menu_file
    print line

def pageTail() :
  file=open(source_files+"/tail.html")
  file_html=file.read()
  print file_html


