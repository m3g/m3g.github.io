#!/usr/bin/python

from sys import path
import cgi
import os.path
import comtypes.client

#
# Get form data
#

formData = cgi.FieldStorage()

file=formData.getvalue('file')
pdb_name=formData['pdb_file'].filename

email=formData.getvalue('email')

disciplina=formData.getvalue('disciplina')

#email=formData.getvalue('email')
#if ('@' not in email) or ('.' not in email) : 
#  QuitPage("ERROR: Invalid e-mail address.")


docfile=open("/tmp/"+file,'w')
docfile.write(file)

extension = os.path.splittext(file)[1]
if extension == ".docx" :
  sys.
  


if natoms > 1000 :
  QuitPage("ERROR: For server runs, only structures with less than 1000 residues are accepted.")

# Open the topolink input file

topolink_input_file=open(output_dir+'/inputfile.inp','w')

#
# If the user provided the complete input file, just change the PDB file name
#

if input_type == "input_file" : 
  topolink_input=topolink_input.split("\n")
  for line in topolink_input :
    if "linkdir" in line : 
      line="linkdir "+output_dir+"/links"
    if "printlinks" in line : 
      line="printlinks yes"
    topolink_input_file.write(line+"\n")

#
# If the user selected a specific linker, create the input file accordingly
#

if input_type == "linker_type" :
  linker_data=getLinker(linker_type)
  base_input=open('./inputfile.inp','r')
  topolink_input=base_input.read()
  base_input.close()
  topolink_input=topolink_input.split("\n")
  for line in topolink_input :
     if line[0:7] == "linkdir" :
       line=line.replace("TMPDIR",output_dir)
       topolink_input_file.write(line+"\n")
     elif line[0:7] == "pdbfile" : 
       line=line.replace("PDB_NAME",pdb_name)
       topolink_input_file.write(line+"\n")
     elif line[0:11] == "#interchain" :
       if interchain : 
         topolink_input_file.write("interchain"+"\n")
       else :
         topolink_input_file.write(line+"\n")
     elif line[0:11] == "EXPERIMENT" :
       topolink_input_file.write("experiment "+linker_type+"\n")
       topolink_input_file.write("  # Possible types of cross-links and maximum distances"+"\n")
       topolink_input_file.write("  #        ResType  Chain  ResNum   AtomType    ResType  Chain  ResNum   AtomType  MaxDist"+"\n")
       for link in linker_data : 
         topolink_input_file.write\
         ("  linktype   "+link['res1']+"     all     all      CB          "\
                         +link['res2']+"     all     all      CB        "\
                         +'{:5.2f}'.format(link['dist'])+"\n")
       topolink_input_file.write("end experiment "+linker_type+"\n")
     else : 
       topolink_input_file.write(line+"\n")


#
# If the user selected a residue pair, create the input file accordingly
#

if input_type == "residue_pair" :
  base_input=open('./inputfile.inp','r')
  topolink_input=base_input.read()
  base_input.close()
  topolink_input=topolink_input.split("\n")
  for line in topolink_input :
     if line[0:7] == "linkdir" :
       line=line.replace("TMPDIR",output_dir)
       topolink_input_file.write(line+"\n")
     elif line[0:7] == "pdbfile" : 
       line=line.replace("PDB_NAME",pdb_name)
       topolink_input_file.write(line+"\n")
     elif line[0:11] == "EXPERIMENT" :
       topolink_input_file.write("experiment "+residue1+'-'+residue2+"\n")
       topolink_input_file.write("  # Possible types of cross-links and maximum distances"+"\n")
       topolink_input_file.write("  #        ResType  Chain  ResNum   AtomType    ResType  Chain  ResNum   AtomType  MaxDist"+"\n")
       topolink_input_file.write\
         ("  linktype   "+residue1+"     all     all      CB          "\
                         +residue2+"     all     all      CB        "\
                         +'{:5.2f}'.format(linker_length)+"\n")
       topolink_input_file.write("end experiment "+residue1+'-'+residue2+"\n")
     else : 
       topolink_input_file.write(line+"\n")

#
# Close the input file
#

topolink_input_file.close()

# Run TopoLink

try : 
  #system("./topolink.sh "+output_dir+"/inputfile.inp "+output_dir+"/"+pdb_name+" "+output_dir+"/topolink.log "+email)
  system("./topolink.sh "+output_dir+"/inputfile.inp "+output_dir+"/"+pdb_name+" "+output_dir+"/topolink.log ")
except : 
  QuitPage("ERROR: Could not run TopoLink on the server. Please report the error to: leandromartinez89@gmail.com")

#
# Create temporary index.shtml file on output dir
#

index_base=open("./index_temp_base.shtml","r")
html_base=index_base.read()
index_base.close()

html_base=html_base.replace("PDB_NAME",pdb_name)

html_out=open(output_dir+"/index.shtml","w")
html_out.write(html_base)
html_out.close()

# Return page


head_and_tail.htmlTop()
print r'<html>'
print r'<meta HTTP-EQUIV="REFRESH" content="0; url='+output_dir+r'/index.shtml">'
print r'</html>'













