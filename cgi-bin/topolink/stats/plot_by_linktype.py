#!/usr/bin/python

import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt

#matplotlib.rc('text', usetex=True)     

protein_type = sys.argv[1]
protein_size = int(sys.argv[2])
residue1 = sys.argv[3]
residue2 = sys.argv[4]
length = float(sys.argv[5])

#print(protein_type, protein_size, residue1, residue2, length)

# Name of each type of protein:

name_type = range(3)
name_type[0] = r'mostly $\alpha$-helical'
name_type[1] = r'mostly $\beta$-sheet'
name_type[2] = r'mixed $\alpha/\beta$'

if protein_type == 'alpha' : itype = 0
if protein_type == 'beta' : itype = 1
if protein_type == 'alpha-beta' : itype = 2

try :
  data_file = open( protein_type+"_"+residue1+residue2+'.dist' )
except :
  try : 
    data_file = open( protein_type+"_"+residue2+residue1+'.dist' )
  except : 
    print(' ERROR: could not find data file. ')


avg_size, linker_length, link_count_avg, pdb_count, size_min, size_max = \
   np.loadtxt(data_file,usecols=(0,1,2,4,5,6),unpack=True,comments="#",dtype=float)

ndata = np.size(avg_size) 
filter = np.zeros(ndata, dtype=bool)

#
# Filter by linker length, and plot as a function of protein size
#

i=0
for value in linker_length :
  if np.abs(value-length) < 0.1 :
    filter[i] = True
  i=i+1

x = avg_size[filter]
y = link_count_avg[filter]


plt.subplot(211)

plt.title('Proteins of type: '+name_type[itype]+'\n'+\
          'Residue types: '+residue1+'-'+residue2)

plt.xlabel('Protein size (number of residues)')
plt.ylabel('Number of cross-linkable'+'\n'+' residue pairs')
plt.annotate('Linker length:'+'{:5.1f}'.format(length),xy=(0.02,0.9),size=12,xycoords='axes fraction')
plt.plot(x,y)

#
# Filter by protein size, and plot as a funciton of linker length
#

filter = np.zeros(ndata, dtype=bool)
i=0
for value in size_min :
  if (protein_size >= size_min[i]) and (protein_size <= size_max[i]) :
    filter[i] = True
  i=i+1

x = linker_length[filter]
y = link_count_avg[filter]

plt.subplot(212)

#plt.title('Proteins of type: '+protein_type+'\n'+\
#          'Linker length: '+'{:5.1f}'.format(length)+'\n'+\
#          'Residue types: '+residue1+' '+residue2)

plt.xlabel('Linker length')
plt.ylabel('Number of cross-linkable'+'\n'+' residue pairs')
plt.annotate('Protein size:'+'{:3}'.format(protein_size),xy=(0.02,0.9),size=12,xycoords='axes fraction')
plt.plot(x,y)


plt.subplots_adjust(left=0.14, 
                    bottom=0.10, 
                    right=0.95, 
                    top=0.90, 
                    wspace=0.3, 
                    hspace=0.3)
plt.gcf().set_size_inches(6,6)

plt.savefig('~/Dropbox/temp/plot.pdf')













