import json
import numpy as np
from matplotlib import pyplot as plt
def read_json(file_name):
    with open(file_name) as json_file:
        dict = json.load(json_file)
    return dict

atmos_dict = read_json('atmospheric_characteristics.json')
wind_dict = read_json('windprofiles_sonde_data.json')



def Vgnd_over_sink(C_L, rho, V_wind, C_D0, AR, oss, mass, S, g):
    C_D = C_D0 + C_L ** 2 / (np.pi * AR * oss)
    V_T = np.sqrt((2 * mass * g) / (rho * S * C_L))
    return (V_T - V_wind) / (V_T * (C_D / C_L))

def windspeed(wind_dict, windspeed_ssd, h):
    return wind_dict[str(h)]['mean'] + windspeed_ssd * wind_dict[str(h)]['ssd']

def reynolds(rho, V_TAS, D, mu):
    return (rho * V_TAS * D) / mu

def Optimal_flight(C_L_lst, rho_h, windspeed_h, C_D0, AR, oss, mass, S, g, plot=False):

    Vgnd_over_sink_h = Vgnd_over_sink(C_L_lst, rho_h, windspeed_h, C_D0, AR, oss, mass, S, g)
    Vgnd_over_sink_opt = np.max(Vgnd_over_sink_h)
    index = np.where(Vgnd_over_sink_h == Vgnd_over_sink_opt)
    C_L_opt = float(C_L_lst[index[0]])
    if plot == True:
        plt.plot(C_L_lst, Vgnd_over_sink_h)


    L_over_D_opt = C_L_opt / (C_D0 + (C_L_opt ** 2) / (np.pi * AR * oss))
    V_TAS_opt = np.sqrt((mass * g) / S * 2 / rho_h * 1 / C_L_opt)
    return C_L_opt, L_over_D_opt, V_TAS_opt, Vgnd_over_sink_opt
g = 9.81
C_D0 = 0.035
AR = 14
oss = 0.9
mass = 18.671 / g
S = 0.0878
mac = 0.077
#Re_limit = 50_000
#controlled_h = 27_000
windspeed_ssd = 2 # 0, 1, 2, 3 ssd
C_L_lst = np.arange(0.01, 2, 0.01)
dh = 10
h_range = np.arange(30_000, 1_000, -dh)


flight_dict = {"h": [], "C_L": [], "L_over_D": [], "V_T": [], "mach": [], "V_ground": [], "V_descent": [], "D_t": [], "distance_travelled": [], "reynolds_number": []}
distance = 0
Re_number = 0
V_TAS = 0
trigger = False

for i_h, h in enumerate(h_range):

    rho_h = atmos_dict[str(h)]['density']
    windspeed_h = windspeed(wind_dict, windspeed_ssd, h)
    dyn_mu = atmos_dict[str(h)]['dynamic_viscosity']
    a_h = atmos_dict[str(h)]["a"]
    C_L_opt, L_over_D_opt, V_TAS_opt, Vgnd_over_sink_opt = \
        Optimal_flight(C_L_lst, rho_h, windspeed_h, C_D0, AR, oss, mass, S, g, plot=True)

    if (V_TAS < V_TAS_opt) and (trigger is False):
        C_L = 0

        L_over_D = 0

        V_TAS = np.sqrt(2 * g * (h_range[0] - h))

        V_ground = 0

        V_descent = V_TAS

        D_t = dh / V_descent

        distance = 0

        Re_number = reynolds(rho_h, V_TAS, mac, dyn_mu)



    else:
        trigger = True
        print(h)
        C_L = C_L_opt

        L_over_D = L_over_D_opt

        V_TAS = V_TAS_opt

        V_ground = V_TAS - windspeed_h

        V_descent = V_TAS / L_over_D_opt

        D_t = dh / V_descent

        Re_number = reynolds(rho_h, V_TAS, mac, dyn_mu)

        dist_at_alt = dh * Vgnd_over_sink_opt
        distance += dist_at_alt

    flight_dict["h"].append(int(h))
    flight_dict["C_L"].append(C_L)
    flight_dict["L_over_D"].append(L_over_D)
    flight_dict["V_T"].append(V_TAS)
    flight_dict["mach"].append(V_TAS / a_h)
    flight_dict["V_ground"].append(V_ground)
    flight_dict["V_descent"].append(V_descent)
    flight_dict["D_t"].append(D_t)
    flight_dict["reynolds_number"].append(Re_number)
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

plt.plot(flight_dict["mach"], h_range)
plt.xlabel("Mach Number, M [-]")
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
