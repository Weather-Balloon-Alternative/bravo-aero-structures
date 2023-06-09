import json
import matplotlib.pyplot as plt
import numpy as np
import ambiance

def read_json(file_name):
    with open(file_name) as json_file:
        dict = json.load(json_file)
    return dict

atmos_dict = read_json('atmospheric_characteristics.json')

t_exp = []


dh = 10
h_range = np.arange(38_000, -420, -dh)

t_mod = ambiance.Atmosphere(h_range).temperature

t_exp = []
h_exp = []
t_d = []

with open('NLM00006260-data-2.txt', 'r') as file:
    lines = file.readlines()
    for line in lines:
        line = line.replace("B-", " -").replace("A-", " -")
        if line.strip().split()[0] != "#NLM00006260" and line.strip().split()[4] != '-9999' and line.strip().split()[3] != '-9999' and line.strip().split()[4] != '-8888' and line.strip().split()[3] != '-8888' :
            t_exp.append(int(line.strip().split()[4].replace("B", "").replace("A", ""))/10 + 273.15)
            h_exp.append(int(line.strip().split()[3].replace("B", "").replace("A", '')))

for i_h, h in enumerate(h_exp):
    t_d.append(t_exp[i_h] - ambiance.Atmosphere(h).temperature)

print(len(t_d))
plt.boxplot(t_d)

plt.show()

plt.plot(t_mod, h_range, c = "r", zorder = 10)
plt.scatter(t_exp, h_exp, zorder = 0)
plt.ylim(-3000, 50000)
plt.xlim(190, 315)
plt.show()

