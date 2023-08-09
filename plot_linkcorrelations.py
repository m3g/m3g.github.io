#!/usr/bin/python

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from sys import exit

font = {'fontname':'FreeSans','size':'14'}

# Label the crosslinks, in this case 26 crosslinks

labels=range(26)
labels[0]="M1-K17"
labels[1]="M1-K113"
labels[2]="K6-E9"
labels[3]="K6-K17"
labels[4]="K6-K113"
labels[5]="E9-K17"
labels[6]="E13-K17"
labels[7]="E13-E20"
labels[8]="E13-D22"
labels[9]="K17-K113"
labels[10]="E28-E53"
labels[11]="E32-E53"
labels[12]="D38-E53"
labels[13]="D38-S57"
labels[14]="E62-E86"
labels[15]="E62-D135"
labels[16]="E66-E86"
labels[17]="S71-D75"
labels[18]="S71-S78"
labels[19]="D75-S78"
labels[20]="E86-K99"
labels[21]="K99-S131"
labels[22]="K99-S133"
labels[23]="E111-D115"
labels[24]="E111-S116"
labels[25]="K113-S116"

# Read output of linkcorrelations

data = np.loadtxt("./linkcorrelations.dat",unpack=True,comments="#")

# Plot correlations matrix, with labels

fig, ax = plt.subplots()

ax.imshow(data, cmap='RdBu', interpolation='none', 
          norm=colors.Normalize(vmin=-1.,vmax=1.))

for (i, j), z in np.ndenumerate(data):
    ax.text(j, i, '{:0.1f}'.format(z), ha='center', va='center')

plt.ylabel("Crosslink",**font)
plt.gca().invert_yaxis()
start, end = ax.get_ylim()
start=start+0.5
end=end+0.5
ax.yaxis.set_ticks(np.arange(start,end,1.))
ax.set_yticklabels(labels,**font)

plt.xlabel("Crosslink",**font)
start, end = ax.get_xlim()
start=start+0.5
end=end+0.5
ax.xaxis.set_ticks(np.arange(start,end,1.))
ax.set_xticklabels(labels,**font)
plt.xticks(rotation='vertical')

plt.show()

