#!/usr/bin/python

import cgi
import head_and_tail 
from numpy import loadtxt, size, zeros, abs, arange, max
from random import randint
from os import makedirs
from sys import exit

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 

def getLinker(linker_type) :
  linker_file="./linkers/"+linker_type+".linktype"
  linker_data = loadtxt(linker_file,
                        dtype={'names':('res1','res2','dist'),
                               'formats':('S3','S3','f8')})
  return linker_data

# topolink main page directory

topolink_dir="../../topolink"

# Get form data

formData = cgi.FieldStorage()

protein_type=formData.getvalue('protein_type')
protein_size=float(formData.getvalue('protein_size'))
linker_type=formData.getvalue('linker_type')

# Print html header and page head

head_and_tail.htmlTop()
head_and_tail.pageHead()

linker_data=getLinker(linker_type)

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

nlinks=[]
labels=[]
for linktype in linker_data :

  residue1=linktype['res1']
  residue2=linktype['res2']
  dist=linktype['dist']

  try : 
    data_file='./data_files/'+protein_type+'_'+residue1+residue2+'.dist'
    avg_size, linker_length, link_count_avg, pdb_count, size_min, size_max = \
       loadtxt(data_file,usecols=(0,1,2,4,5,6),unpack=True,comments="#",dtype=float)
  except : 
    data_file='./data_files/'+protein_type+'_'+residue2+residue1+'.dist'
    avg_size, linker_length, link_count_avg, pdb_count, size_min, size_max = \
       loadtxt(data_file,usecols=(0,1,2,4,5,6),unpack=True,comments="#",dtype=float)

  ndata = size(avg_size)
  filter=zeros(ndata,dtype=bool)

  #print (residue1+' '+residue2+' '+'{:5.1f}'.format(dist)+'<br>')

  i=0
  x = 0.
  for value in linker_length :
    if i == ndata : 
      break
    if (dist >= linker_length[i] and dist <= linker_length[i+1] ) and \
       (protein_size >= size_min[i] and protein_size <= size_max[i]) :

      x = link_count_avg[i+1]

      # Interpolate on linker length

      x = link_count_avg[i-1] + \
          (dist-linker_length[i-1])/(value-linker_length[i-1])*(link_count_avg[i]-link_count_avg[i-1])

      # Interpolate on protein size
      
      x = x * protein_size/avg_size[i]

      break
    i = i + 1

  #output_data_file.write(residue1+" "+residue2+" "+'{:5.2f}'.format(x)+"\n")
  labels.append(residue1+'-'+residue2)
  nlinks.append(x)

# Write data file

#output_data_file.close()

# Prepare plot

ind = arange(size(labels))

if protein_type == "alpha" : 
  name_type = r'mostly $\alpha$-helical'
  checked = "CHECKEDA"
if protein_type == "beta" : 
  name_type = r'mostly $\beta$-sheet'
  checked = "CHECKEDB"
if protein_type == "alpha-beta" : 
  name_type = r'mixed $\alpha/\beta$'
  checked = "CHECKEDALL"

plt.suptitle('Proteins of type: '+name_type,fontsize=14)
plt.title('Linker type:'+linker_type+' - Sequence length:'+'{:3}'.format(int(protein_size)),fontsize=14)

plt.xlabel('Residue pair',fontsize=14)
plt.ylabel('Number of cross-linkable'+'\n'+' residue pairs',fontsize=14)
plt.bar(ind+0.20,nlinks,0.5,alpha=0.5)
plt.xticks(ind+0.18,labels,rotation=90,size=14)
plt.yticks(size=14)

plt.ylim(0,max(nlinks)+1)
plt.xlim(-1,size(labels))

plt.subplots_adjust(left=0.15,
                    bottom=0.20,
                    right=0.90, 
                    top=0.90, 
                    wspace=0.3, 
                    hspace=0.3)
plt.gcf().set_size_inches(6,6)

pngplot=output_dir+'/plot.png'
plt.savefig(pngplot)


# Return page

return_page=open("./by_linker_type-return.shtml","r")
return_page=return_page.read()

return_page=return_page.replace(checked,"checked")
return_page=return_page.replace("CHECKEDALL"," ")
return_page=return_page.replace("CHECKEDA"," ")
return_page=return_page.replace("CHECKEDB"," ")

return_page=return_page.replace("PROTEIN_SIZE",'{}'.format(int(protein_size)))

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

return_page=return_page.replace(selected,"selected")
return_page=return_page.replace("DSS_SELECTED"," ")
return_page=return_page.replace("DSG_SELECTED"," ")
return_page=return_page.replace("16DIAMINE_SELECTED"," ")
return_page=return_page.replace("13DIAMINE_SELECTED"," ")
return_page=return_page.replace("ZEROLENGTH_SELECTED"," ")

return_page=return_page.replace("PNGPLOT",pngplot)

print return_page

i=0
for linktype in linker_data :
  residue1=linktype['res1']
  residue2=linktype['res2']
  dist=linktype['dist'] 
  print "<tr><td>"+residue1+"</td><td>"+residue2+\
         "</td><td align=center>"+'{:5.1f}'.format(dist)+"</td>"+\
         "<td align=center>"+'{:5.1f}'.format(nlinks[i])+"</td></tr>"
  i=i+1


print "</table>"
print "</td></tr>"
print "</table>"

# Tail of standard page

head_and_tail.pageTail()


             
