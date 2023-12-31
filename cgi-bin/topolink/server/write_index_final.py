#!/usr/bin/python

import os.path
from sys import argv
from sys import path
path.insert(0,'../stats')
import head_and_tail

import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import sys, time, os

output_dir = argv[1]
pdb_name = argv[2]
inputfile = argv[3]

def LinkData(line) :

  residue1 = line[8:12].strip()
  chain1 = line[13]
  resnum1 = line[15:19].strip()

  residue2 = line[25:29].strip()
  chain2 = line[30]
  resnum2 = line[32:36].strip()

  atom1 = line[20:24].strip()
  atom2 = line[37:41].strip()

  return residue1, chain1, resnum1, residue2, chain2, resnum2, atom1, atom2

def LinkLine(line,pdb_name,found) :

  residue1, chain1, resnum1, residue2, chain2, resnum2, atom1, atom2 = LinkData(line)

  euclidean = line[42:50].strip()
  topological = line[51:60].strip()

  solvent1 = line[118]
  solvent2 = line[119]

  line = "<tr>"+\
         "<td align=center>"+residue1+"</td>"+\
         "<td align=center>"+chain1+"</td>"+\
         "<td align=center>"+resnum1+"</td>"+\
         "<td align=center>"+atom1+"</td>"+\
         "<td align=center>"+residue2+"</td>"+\
         "<td align=center>"+chain2+"</td>"+\
         "<td align=center>"+resnum2+"</td>"+\
         "<td align=center>"+atom2+"</td>"+\
         "<td align=center>"+euclidean+"</td>"+\
         "<td align=center>"+topological+"</td>"+\
         "<td align=center>"+solvent1+"</td>"+\
         "<td align=center>"+solvent2+"</td>"

  check_file = os.path.isfile(output_dir+"/links/"+pdb_name[:-4]+"_"+residue1+chain1+resnum1+atom1+"-"\
                                                                    +residue2+chain2+resnum2+atom2+".pdb")
  if check_file :
    link_file = "./links/"+pdb_name[:-4]+"_"+residue1+chain1+resnum1+atom1+"-"\
                                            +residue2+chain2+resnum2+atom2+".pdb"
    line = line+"<td align=center><a target=_blank_ href="+link_file+">[PDB]</a></td>"

  else :
    link_file = " "
    line = line+"<td align=center> </td>"
  line = line+"</tr>"

  return line, link_file

print output_dir, pdb_name

index_base=open("./index_final_base.html","r")
html_base=index_base.read()
index_base.close()
html_base=html_base.replace("PDB_NAME",pdb_name)
html_base=html_base.split("\n")

html_out=open(output_dir+"/index.html","w")

# Check if there are errors in the topolink log file

log_error = False
try : 
  log=open(output_dir+"/topolink.log","r")
  topolink_log=log.read()
  log.close()
  if "ERROR" in topolink_log : 
    log_error = True
    log_error_message = "ERROR: An execution error was reported by TopoLink. Please check the log file carefully."
  topolink_log=topolink_log.split("\n")
      
except :
  log_error = True
  log_error_message = "ERROR: Could not load TopoLink log file. Please contact: leandromartinez98@gmail.com"

for line in html_base :

  if "PDB_NAME" in line : 
    line=line.replace("PDB_NAME",pdb_name)
    html_out.write(line+"\n")

  if "ERROR" in line :  
    if log_error : 
      html_out.write("<center><font color=red><b>"+log_error_message+"</font></center>"+"\n")

      html_out.write(r'<textarea rows="50" cols="124">')
      for log_line in topolink_log :
        html_out.write(log_line)
 
      html.out.write(r'</textarea></center>')
      
      html_out.write(head_and_tail.pageTail())
      quit()

  elif line[0:15] == "LINK_LIST_FOUND" :  

    for logline in topolink_log : 
      if logline[2:7] == "LINK:" : 
        if logline[94:101] == "MISSING" or logline[92:101] == "OK: FOUND" :
          found=True
          line, link_file = LinkLine(logline,pdb_name,found)
          html_out.write(line+"\n")

  elif line[0:18] == "LINK_LIST_NOTFOUND" :  

    for logline in topolink_log : 
      if logline[2:7] == "LINK:" : 
        if ( not "MISSING" in logline ) and ( not "OK: FOUND" in logline ) : 
          found=False
          line, link_file = LinkLine(logline,pdb_name,found)
          html_out.write(line+"\n")

  else :
    html_out.write(line+"\n")

html_out.close()

#
# Create pymol session
#

import pymol

pymol.finish_launching()

for logline in topolink_log : 
  if logline[2:7] == "LINK:" : 
    residue1, chain1, resnum1, residue2, chain2, resnum2, atom1, atom2 = LinkData(logline)
    link_file_name = output_dir+"/links/"+pdb_name[:-4]+"_"+residue1+chain1+resnum1+atom1+"-"\
                                                           +residue2+chain2+resnum2+atom2+".pdb"
    link_file = os.path.isfile(link_file_name)

    if ("LONG" or "OK: FOUND" or "BAD: SHORT" or "MISSING") in logline :
      if link_file :
        selname=residue1+resnum1+chain1+"-"+residue2+resnum2+chain2
        pymol.cmd.load(link_file_name,selname)
        pymol.cmd.show_as("sticks",selname)

selname=pdb_name[:-4]
pymol.cmd.load(output_dir+"/"+pdb_name,selname)
pymol.cmd.show_as("cartoon",selname)
pymol.cmd.color("cyan",selname)
pse_file=output_dir+"/links.pse"
pymol.cmd.save(pse_file)
pymol.cmd.bg_color(color="white")

pymol.cmd.center(selname)
pymol.cmd.png(output_dir+"/links-1.png",200,200,quiet=0)
pymol.cmd.rotate("z",120)
pymol.cmd.png(output_dir+"/links-2.png",200,200,quiet=0)
pymol.cmd.rotate("z",120)
pymol.cmd.png(output_dir+"/links-3.png",200,200,quiet=0)

os.system("convert -trim "+output_dir+"/links-1.png "+output_dir+"/links-1.png")
os.system("convert -trim "+output_dir+"/links-2.png "+output_dir+"/links-2.png")
os.system("convert -trim "+output_dir+"/links-3.png "+output_dir+"/links-3.png")

pymol.cmd.quit()

