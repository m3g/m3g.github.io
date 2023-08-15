#!/usr/bin/python

import cgi
import head_and_tail 
from numpy import loadtxt, size, zeros, abs, arange
from random import randint
from os import makedirs

import itertools
style=itertools.cycle(["-","--",":",":",".","h","H"])

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 

def getLinker(linker_type) :
  linker_file="./linkers/"+linker_type+".linktype"
  linker_data = loadtxt(linker_file,
                        dtype={'names':('res1','res2','dist'),
                               'formats':('S3','S3','f8')})
  return linker_data

def getResidues(residue_pair) :
  residue1 = residue_pair[0:3]
  residue2 = residue_pair[4:7]
  return residue1, residue2

#def sidechainLength(residue)
#  if residue == "SER" : length=1.45
#  if residue == "ARG" : length=5.35
#  if residue == "ASP" : length=4.6
#  if residue == "GLU" : length=5.9

# Print html header and page head

head_and_tail.htmlTop()
head_and_tail.pageHead()

# Get form data

formData = cgi.FieldStorage()

protein_type=formData.getvalue('protein_type')

linker_option=formData.getvalue('linker_option')

linker_type=formData.getvalue('linker_type')

residue_pair=formData.getvalue('residue_pair')
linker_length=float(formData.getvalue('linker_length'))

residue1, residue2 = getResidues(residue_pair)

# Name of protein type

if protein_type == "alpha" : 
  name_type = r'mostly $\alpha$-helical'
  checked = "CHECKEDA"
if protein_type == "beta" : 
  name_type = r'mostly $\beta$-sheet'
  checked = "CHECKEDB"
if protein_type == "alpha-beta" : 
  name_type = r'mixed $\alpha/\beta$'
  checked = "CHECKEDALL"

# Output directory

output_dir='../../../topolink/serverfiles/TEMP'+'{}'.format(randint(0,100000000)+100000000)
makedirs(output_dir)

## Output data file
#
#output_data_file=open(output_dir+"/data.dat","w")
#output_data_file.write("# Protein type: "+protein_type+"\n")
#output_data_file.write("# Sequence length: "+'{:5}'.format(protein_size)+"\n")
#output_data_file.write("# Linker type: "+linker_type+"\n")
#output_data_file.write("# Residue1   Residue2   Number of crosslinkable pairs"+"\n")

#
# If a predefined linker was chosen
#

if linker_option == 'by_linker_type' :
  linker_data=getLinker(linker_type)

  plt.suptitle('Proteins of type: '+name_type,fontsize=14)
  plt.title('Linker type: '+linker_type,fontsize=14)
  
  plt.xlabel('Sequence length',fontsize=14)
  plt.ylabel('Number of cross-linkable'+'\n'+' residue pairs',fontsize=14)
  plt.xlim(50,300)

  ylim = 0.
  for linktype in linker_data :

    residue1=linktype['res1']
    residue2=linktype['res2']
    dist=linktype['dist']

    try :
      data_file='./data_files/'+protein_type+'_'+residue1+residue2+'.dist'
      avg_size, link_dist, link_count_avg, pdb_count, size_min, size_max = \
         loadtxt(data_file,usecols=(0,1,2,4,5,6),unpack=True,comments="#",dtype=float)
    except :
      data_file='./data_files/'+protein_type+'_'+residue2+residue1+'.dist'
      avg_size, link_dist, link_count_avg, pdb_count, size_min, size_max = \
         loadtxt(data_file,usecols=(0,1,2,4,5,6),unpack=True,comments="#",dtype=float)

    ndata = size(avg_size)
    filter=zeros(ndata,dtype=bool)

    #
    # Filter by linker length, and plot as a function of protein size
    #
    
    i=1
    for value in link_dist :
      if (dist >= link_dist[i-1]) and (dist < link_dist[i]) :
        filter[i] = True
      i=i+1
      if i == ndata : break

    x = avg_size[filter]
    y = link_count_avg[filter]
    #plt.plot(x,y,label=residue1+'-'+residue2)
    plt.plot(x,y,label=residue1+'-'+residue2,linestyle=style.next(),linewidth=3)
    i=0
    for value in x : 
      if value > 300 : break
      i=i+1
    ylim = max(ylim,y[i+1])

  plt.xticks(fontsize=14)
  plt.yticks(fontsize=14)
  plt.legend(loc=0,prop={'size':14})
  plt.ylim(0,ylim)
  plt.subplots_adjust(left=0.15,
                      bottom=0.15,
                      right=0.90, 
                      top=0.90, 
                      wspace=0.3, 
                      hspace=0.3)
  plt.gcf().set_size_inches(6,6)
  
  pngplot=output_dir+'/plot.png'
  plt.savefig(pngplot)

#
# else, if a residue pair and linker length was chosen
#

if linker_option == 'by_residue_pair' :
  
  try : 
    data_file='./data_files/'+protein_type+'_'+residue1+residue2+'.dist'
    avg_size, link_dist, link_count_avg, pdb_count, size_min, size_max = \
       loadtxt(data_file,usecols=(0,1,2,4,5,6),unpack=True,comments="#",dtype=float)
  except : 
    data_file='./data_files/'+protein_type+'_'+residue2+residue1+'.dist'
    avg_size, link_dist, link_count_avg, pdb_count, size_min, size_max = \
       loadtxt(data_file,usecols=(0,1,2,4,5,6),unpack=True,comments="#",dtype=float)

  ndata = size(avg_size) 
  filter = zeros(ndata, dtype=bool)
 
  #
  # Filter by linker length, and plot as a function of protein size
  #
  
  i=1
  for value in link_dist :
    if (linker_length >= link_dist[i-1]) and (linker_length < link_dist[i]) :
      filter[i] = True
    i=i+1
    if i == ndata : break
  
  x = avg_size[filter]
  y = link_count_avg[filter]

  # Create plot

  plt.suptitle('Proteins of type: '+name_type,fontsize=14)
  plt.title('Residue pair: '+residue1+'-'+residue2+' Linker length: '+'{:5.2f}'.format(linker_length))
  
  plt.xlabel('Sequence length',fontsize=14)
  plt.ylabel('Number of cross-linkable'+'\n'+' residue pairs')
  plt.plot(x,y)
  plt.xlim(50,300)
  i=0
  for value in x : 
    if value > 300 : break
    i=i+1
     
  plt.ylim(0,y[i+1])
  
  plt.subplots_adjust(left=0.13,
                      bottom=0.20,
                      right=0.97, 
                      top=0.90, 
                      wspace=0.3, 
                      hspace=0.3)
  plt.gcf().set_size_inches(6,6)
  
  pngplot=output_dir+'/plot.png'
  plt.savefig(pngplot)

# Write data file

#output_data_file.close()

# Return page

return_page=open("./by_protein_size-return.shtml","r")
return_page=return_page.read()

return_page=return_page.replace(checked,"checked")
return_page=return_page.replace("CHECKEDALL"," ")
return_page=return_page.replace("CHECKEDA"," ")
return_page=return_page.replace("CHECKEDB"," ")

if linker_type == "DSS" : 
  selected="DSS_SELECTED"
if linker_type == "DSG" : 
  selected="DSG_SELECTED"
if linker_type == "16Diamine" : 
  selected="16DIAMINE_SELECTED"
if linker_type == "13Diamine" : 
  selected="13DIAMINE_SELECTED"
if linker_type == "Zero-length" : 
  selected="ZEROLENGTH_SELECTED"

if linker_option == "by_linker_type" : 
  return_page=return_page.replace("LINKER_TYPE_CHECKED","checked")
  return_page=return_page.replace("RESIDUE_PAIR_CHECKED"," ")
if linker_option == "by_residue_pair" : 
  return_page=return_page.replace("LINKER_TYPE_CHECKED"," ")
  return_page=return_page.replace("RESIDUE_PAIR_CHECKED","checked")

return_page=return_page.replace(selected,"selected")
return_page=return_page.replace("DSS_SELECTED"," ")
return_page=return_page.replace("DSG_SELECTED"," ")
return_page=return_page.replace("16DIAMINE_SELECTED"," ")
return_page=return_page.replace("13DIAMINE_SELECTED"," ")
return_page=return_page.replace("ZEROLENGTH_SELECTED"," ")

return_page=return_page.replace("LINKER_LENGTH",'{:5.2f}'.format(linker_length))
return_page=return_page.replace("PNG_PLOT",pngplot)

return_page = return_page.split("\n")
for return_line in return_page : 

  if return_line[4:13] == "PAIR_LIST" :
    pair_list = open("../../../topolink/stats/pairs.shtml","r")
    pair_list = pair_list.read()
    pair_list = pair_list.split("\n")
    for line in pair_list :
      if line[15:22] == residue_pair : 
        print "<option value="+residue_pair+" selected >"+residue_pair+"</option>"
      else :
        print line

  else :

   if linker_option == "by_linker_type" : 
     if return_line[9:26] == "LINKER_REACTIVITY" : 
       print "Linker reactivity:"
       print "<table align=center>"
       print "<tr><td>Residue</td><td>Residue</td><td>Reach (99%)</td></tr>"
       for linktype in linker_data :
         residue1=linktype['res1']
         residue2=linktype['res2']
         dist=linktype['dist']
         print "<tr><td>"+residue1+"</td><td>"+residue2+"</td><td align=center>"+'{:5.1f}'.format(dist)+"</td></tr>"
       print "</table>"

  print return_line

return_page.close()
pair_list.close()

# Tail of standard page

head_and_tail.pageTail()


