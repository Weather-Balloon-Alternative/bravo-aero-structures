import json
import numpy as np
from matplotlib import pyplot as plt
def read_json(file_name):
    with open(file_name) as json_file:
        dict = json.load(json_file)
    return dict

atmos_dict = read_json('atmospheric_characteristics.json')
wind_dict = read_json('windprofiles_sonde_data.json')

def Vgnd_over_D(C_L, rho, V_wind, C_D0, AR, oss, mass, S, g):
    C_D = C_D0 + C_L ** 2 / (np.pi * AR * oss)
    return (np.sqrt((mass * g) / S * 2 / rho * 1 / C_L) - V_wind) / (C_D / C_L * (mass * g))


def Vgnd_over_sink(C_L, rho, V_wind, C_D0, AR, oss, mass, S, g):
    C_D = C_D0 + C_L ** 2 / (np.pi * AR * oss)
    V_T = np.sqrt((2 * mass * g) / (rho * S * C_L))
    return (V_T - V_wind) / (V_T * (C_D / C_L))



g = 9.81
C_D0 = 0.035
AR = 15
oss = 0.9
mass = 18.671 / g
S = 0.0878
mac = 0.077
controlled_h = 27_000
C_L = np.arange(0.01, 1.51, 0.01)
dh = 10
h_range = np.arange(30_000, 1_000, -dh)


flight_dict = {"h": [], "C_L": [], "L_over_D": [], "V_T": [], "V_ground": [], "V_descent": [], "D_t": [], "distance_travelled": [], "reynolds_number": []}
distance = 0

windspeed_ssd = 2 # 0, 1, 2, 3 ssd

def windspeed(wind_dict, windspeed_ssd, h):
    return wind_dict[str(h)]['mean'] + windspeed_ssd * wind_dict[str(h)]['ssd']

def reynolds(rho, V_TAS, D, mu):
    return (rho * V_TAS * D) / mu

for h in h_range:
    rho_h = atmos_dict[str(h)]['density']
    windspeed_h = windspeed(wind_dict, windspeed_ssd, h)
    dyn_mu = atmos_dict[str(h)]['dynamic_viscosity']

    Vgnd_over_sink_h = Vgnd_over_sink(C_L, rho_h,
                                    windspeed_h, C_D0, AR, oss, mass, S, g)

    plt.plot(C_L, Vgnd_over_sink_h)

    Vgnd_over_sink_opt = np.max(Vgnd_over_sink_h)
    index = np.where(Vgnd_over_sink_h == Vgnd_over_sink_opt)

    flight_dict["h"].append(int(h))

    flight_dict["C_L"].append(C_L[index[0]])

    L_over_D_opt = C_L[index[0]] / (C_D0 + (C_L[index[0]] ** 2) / (np.pi * AR * oss))
    flight_dict["L_over_D"].append(L_over_D_opt)

    V_TAS = np.sqrt((mass * g) / S * 2 / rho_h * 1 / C_L[index[0]])
    flight_dict["V_T"].append(V_TAS)

    V_ground = V_TAS - windspeed_h
    flight_dict["V_ground"].append(V_ground)

    V_descent = V_TAS / L_over_D_opt
    flight_dict["V_descent"].append(V_TAS / L_over_D_opt)

    D_t = dh / V_descent
    flight_dict["D_t"].append(D_t)

    Re_number = reynolds(rho_h, V_TAS, mac, dyn_mu)
    flight_dict["reynolds_number"].append(Re_number)

    if h <= controlled_h:
        dist_at_alt = dh * Vgnd_over_sink_opt
        distance += dist_at_alt
        flight_dict["distance_travelled"].append(distance)
    else:
        distance += 0
        flight_dict["distance_travelled"].append(distance)


def write_JSON(dictionary, output_file):
    with open(output_file, "w") as outfile:
        json.dump(dictionary, outfile)

file_name = "flight_dict.json"
#write_JSON(flight_dict, file_name)

print(distance / 1000)
print(sum(flight_dict['D_t']) / 3600)

plt.xlabel("Lift Coefficient, C_L [-]")
plt.ylabel("Ground Speed over Descent Speed, V_gnd/h_dot [-]")
plt.show()

plt.plot(flight_dict["L_over_D"], h_range)
plt.xlabel("Lift over Drag, L/D [-]")
plt.ylabel("Altitude, h [m]")
plt.show()

plt.plot(flight_dict["C_L"], h_range)
plt.xlabel("Lift coefficient, C_L [-]")
plt.ylabel("Altitude, h [m]")
plt.show()

plt.plot(flight_dict["V_T"], h_range)
plt.xlabel("True Airspeed, V_TAS [m/s]")
plt.ylabel("Altitude, h [m]")
plt.show()

plt.plot(flight_dict["V_ground"], h_range)
plt.xlabel("Ground speed, V_ground [m/s]")
plt.ylabel("Altitude, h [m]")
plt.show()

plt.plot(flight_dict["V_descent"], h_range)
plt.xlabel("Descent Speed, V_descent [m/s]")
plt.ylabel("Altitude, h [m]")
plt.show()

plt.plot(flight_dict["distance_travelled"], h_range)
plt.xlabel("Distance, D [m]")
plt.ylabel("Altitude, h [m]")
plt.show()

plt.plot(flight_dict["reynolds_number"], h_range)
plt.xlabel("Reynolds Number, Re [...]")
plt.ylabel("Altitude, h [m]")
plt.show()
