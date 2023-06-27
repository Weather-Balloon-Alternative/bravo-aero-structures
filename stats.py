from matplotlib import pyplot as plt
import matplotlib
import json
import numpy as np

print('jydfgasjdh')
stations = ['boulder', 'broadmedows', 'debilt', 'macquarie', 'paramaribo', 'southpole']
station_name = {'boulder':'Boulder, US', 'macquarie': 'Macquarie, AU', 'debilt':'De Bilt, NL', 'paramaribo':'Paramaribo, SR', 'southpole':'South Pole, AQ', 'broadmedows':'Broadmeadows, AU'}
station_pcts = {'boulder':['1.0%','99.0%'], 'macquarie': ['36.1%','63.9%'], 'debilt':['9.1%', '90.9%'], 'paramaribo':['0.0%', '100%'], 'southpole':['1.0%', '99.0%'], 'broadmedows':['3.5%', '96.5%']}
filename = 'results_bravo_mini.json'

ax_ids = [] 

with open(filename) as f:
     data = json.load(f)

fig, ax = plt.subplots(2,3)
for i in range(2):
     for j in range(3):
          ax_ids.append([i, j])

i = 0
bins = np.arange(-300, 440, 20)
cmap = matplotlib.colors.ListedColormap(['red', 'green'])
fig.tight_layout()
for station in stations:
    x = data[station]
    edges, bin, patches = ax[ax_ids[i][0],ax_ids[i][1]].hist(x, bins, color='red')
    mmx = np.max(edges)
    
    for j, p in enumerate(patches):
          if p.get_x() >= 0:
               p.set_facecolor('green')
    ax[ax_ids[i][0],ax_ids[i][1]].set_title(station_name[station])
    ax[ax_ids[i][0],ax_ids[i][1]].set_xlabel('km')
    ax[ax_ids[i][0],ax_ids[i][1]].set_ylabel('number of launches')
    ax[ax_ids[i][0],ax_ids[i][1]].text(-320, mmx*0.98, station_pcts[station][0], size=12, color='red')
    ax[ax_ids[i][0],ax_ids[i][1]].text(345, mmx*0.98, station_pcts[station][1], size=12, color='green')
    #plt.hist(x, bins=50, ax=ax[ax_ids[i][0],ax_ids[i][1]]) 

    i += 1
print('jydfgasjdh')
plt.show()


