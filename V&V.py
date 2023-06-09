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


dh = 1000
h_range = np.arange(30_000, 1_000, -dh)
for i, h in enumerate(h_range):
    t_h = atmos_dict[str(h)]["temperature"]
    print(t_h)
    t_exp.append(t_h)

t_mod = ambiance.Atmosphere(h_range).temperature

plt.plot(t_mod, h_range)
plt.plot(t_exp, h_range)
plt.show()