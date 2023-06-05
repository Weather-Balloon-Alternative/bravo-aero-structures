import json
import numpy as np
from matplotlib import pyplot as plt
def read_json(file_name):
    with open(file_name) as json_file:
        dict = json.load(json_file)
    return dict

atmos_dict = read_json('atmospheric_characteristics.json')

def glide_ratio(atmos_dict, h_start, h_end, max_flighttime, distance, max_mach=0.3):
    dh = 10 #m
    glide_ratio = []
    v_ground_avg = distance / max_flighttime
    #print(v_ground_avg)
    for h in np.arange(h_start, h_end+dh, -dh):
        a = atmos_dict[str(h)]['a'] #m/s
        V_max = max_mach * a #m/s
        V_horizontal_allowable = atmos_dict[str(h)]['maximum_windspeed'] + v_ground_avg
        if V_max < V_horizontal_allowable:
            print(V_horizontal_allowable)
        V_vertical_allowable = np.sqrt(V_max ** 2 - V_horizontal_allowable ** 2)
        required_glide_ratio = V_horizontal_allowable / V_vertical_allowable
        glide_ratio.append(required_glide_ratio)

    return glide_ratio


def Vgnd_over_D(C_L, rho, V_wind, C_D0, AR, oss, mass, S, g):
    C_D = C_D0 + C_L ** 2 / (np.pi * AR * oss)
    return (np.sqrt((mass * g) / S * 2 / rho * 1 / C_L) - V_wind) / (C_D / C_L * (mass * g))


def Vgnd_over_sink(C_L, rho, V_wind, C_D0, AR, oss, mass, S, g):
    C_D = C_D0 + C_L ** 2 / (np.pi * AR * oss)
    V_T = np.sqrt((2 * mass * g) / (rho * S * C_L))
    return (V_T - V_wind) / (V_T * (C_D / C_L))


g = 9.81
C_D0 = 0.0223
AR = 10
oss = 0.9
mass = 4.5
S = 0.21
controlled_h = 27_000
C_L = np.arange(0.01, 1.51, 0.01)
dh = 100
h_range = np.arange(30_000, 0, -dh)


flight_dict = {"h": [], "C_L": [], "L_over_D": [], "V_T": [], "V_ground": [], "V_descent": [], "D_t": []}
distance = 0

for h in h_range:
    flight_dict["h"].append(int(h))

    Vgnd_over_sink_h = Vgnd_over_sink(C_L, atmos_dict[str(h)]['density'],
                                    atmos_dict[str(h)]['average_windspeed'], C_D0, AR, oss, mass, S, g)

    plt.plot(C_L, Vgnd_over_sink_h)

    index = np.where(Vgnd_over_sink_h == np.max(Vgnd_over_sink_h))
    #flight_dict["C_L"].append(C_L[index[0]])

    L_over_D_opt = C_L[index[0]] / (C_D0 + (C_L[index[0]] ** 2) / (np.pi * AR * oss))
    #flight_dict["L_over_D"].append(L_over_D_opt)

    V_T = np.sqrt((mass * g) / S * 2 / atmos_dict[str(h)]['density'] * 1 / C_L[index[0]])
    flight_dict["V_T"].append(V_T)

    V_ground = V_T - atmos_dict[str(h)]['average_windspeed']
    flight_dict["V_ground"].append(V_T - atmos_dict[str(h)]['average_windspeed'])

    V_descent = V_T / L_over_D_opt
    flight_dict["V_descent"].append(V_T / L_over_D_opt)

    D_t = dh / V_descent
    flight_dict["D_t"].append(D_t)

    if h <= controlled_h:
        distance = distance + dh / V_descent * V_ground

def write_JSON(dictionary, output_file):
    with open(output_file, "w") as outfile:
        json.dump(dictionary, outfile)

file_name = "flight_dict.json"
write_JSON(flight_dict, file_name)

print(distance / 1000)
print(sum(flight_dict['D_t']) / 3600)
plt.show()
plt.plot(flight_dict["L_over_D"], h_range)
plt.show()
