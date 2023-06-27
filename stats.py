from matplotlib import pyplot as plt
import matplotlib
import json
import numpy as np

print('jydfgasjdh')
stations = ['boulder', 'broadmedows', 'debilt', 'lindenberg', 'macquarie', 'paramaribo', 'southpole', 'wallops']
station_name = {'boulder':'Boulder, US', 'macquarie': 'Macquarie, AU', 'debilt':'De Bilt, NL', 'paramaribo':'Paramaribo, SR', 'southpole':'South Pole, AQ', 'broadmedows':'Broadmeadows, AU', 'wallops':'Wallops Island, US', 'lindenberg':'Lindenberg, DE'}
station_pcts = {'boulder':['1.0%','99.0%'], 'macquarie': ['36.1%','63.9%'], 'debilt':['9.1%', '90.9%'], 'paramaribo':['0.0%', '100%'], 'southpole':['1.0%', '99.0%'], 'broadmedows':['3.5%', '96.5%'], 'wallops':['13.3%','86.7%'], 'lindenberg':['8.3%','91.7%']}
filename = 'results_bravo_mini.json'

colors = {'green': (108/255, 194/255, 74/255), 'red':(237/255, 104/255, 66/255)}

ax_ids = [] 

with open(filename) as f:
     data = json.load(f)
size_x, size_y = 2, 4
fig, ax = plt.subplots(size_x,size_y)
for i in range(size_x):
     for j in range(size_y):
          ax_ids.append([i, j])

i = 0
bins = np.arange(-300, 440, 20)
cmap = matplotlib.colors.ListedColormap(['red', 'green'])
fig.tight_layout()
for station in stations:
    x = data[station]
    edges, bin, patches = ax[ax_ids[i][0],ax_ids[i][1]].hist(x, bins, color=colors['red'])
    mmx = np.max(edges)
    
    for j, p in enumerate(patches):
          if p.get_x() >= 0:
               p.set_facecolor(colors['green'])
    ax[ax_ids[i][0],ax_ids[i][1]].set_title(station_name[station])
    ax[ax_ids[i][0],ax_ids[i][1]].set_xlabel('km')
    ax[ax_ids[i][0],ax_ids[i][1]].set_ylabel('number of launches')
    ax[ax_ids[i][0],ax_ids[i][1]].text(-320, mmx*0.98, station_pcts[station][0], size=12, color=colors['red'])
    ax[ax_ids[i][0],ax_ids[i][1]].text(310, mmx*0.98, station_pcts[station][1], size=12, color=colors['green'])
    #plt.hist(x, bins=50, ax=ax[ax_ids[i][0],ax_ids[i][1]]) 

    i += 1
print('jydfgasjdh')
plt.show()


