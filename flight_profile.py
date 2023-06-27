import json
import numpy as np
from matplotlib import pyplot as plt
from flight_performance import FlightPerformance


def read_json(file_name):
    with open(file_name) as json_file:
        dict = json.load(json_file)
    return dict


def windspeed(wind_dict, windspeed_ssd, h):
    if windspeed_ssd != None:
        return wind_dict[str(h)]['mean'] + windspeed_ssd * wind_dict[str(h)]['ssd']
    else:
        return 0

atmos_dict = read_json('atmospheric_characteristics.json')
wind_dict = read_json('windprofiles_sonde_data.json')

wind_profile = {}
for h in wind_dict:
    wind_profile[h] = windspeed(wind_dict, 2, h)



AC_params_z = [{"AR": 12, "e": 0.9, "m": 0.75, "S": 0.05, "CD_0": None, "CD_0_base": 0.02771, "CD_0_h": 5.714285714285715e-07, "CL_max": 0.666 * 1.55,
                 "CL_alpha": 6.1859591823509295 * np.pi / 180}, {"AR": 12, "e": 0.9, "m": 2.621, "S": 0.215, "CD_0": None, "CD_0_base": 0.02561, "CD_0_h": 5.714285714285715e-07, "CL_max": 0.666 * 1.55,
                 "CL_alpha": 6.1859591823509295 * np.pi / 180}, {"AR": 12, "e": 0.9, "m": 3.021, "S": 0.215, "CD_0": None, "CD_0_base": 0.02561, "CD_0_h": 5.714285714285715e-07, "CL_max": 0.666 * 1.55,
                 "CL_alpha": 6.1859591823509295 * np.pi / 180}, {"AR": 12, "e": 0.9, "m": 3.521, "S": 0.215, "CD_0": None, "CD_0_base": 0.02561, "CD_0_h": 5.714285714285715e-07, "CL_max": 0.666 * 1.55,
                 "CL_alpha": 6.1859591823509295 * np.pi / 180}, {"AR": 12, "e": 0.9, "m": 4.021, "S": 0.215, "CD_0": None, "CD_0_base": 0.02561, "CD_0_h": 5.714285714285715e-07, "CL_max": 0.666 * 1.55,
                 "CL_alpha": 6.1859591823509295 * np.pi / 180}]

for AC_params in AC_params_z:
    FP = FlightPerformance(wind_profile, AC_params, 33_000, 1_500, 1000)
    FP.flight_sim()
    flight_dict = FP.get_flightdict()
    print(flight_dict["distance_travelled"][-1] / 1000)


def plotting(data, x_label, y_label):
    if np.shape(data)[0] > 2:
        for i in range(0, np.shape(data)[0], 2):
            plt.plot(data[i], data[i+1])
    else:
        plt.plot(data[0], data[1])
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()

plotting(np.array([flight_dict["V_descent"], flight_dict["h"]]), "Descent speed, V_descent [m/s]", "Altitude, h [m]")

'''print(flight_dict["distance_travelled"][-1] / 1000)
print(sum(flight_dict['D_t']) / 3600)'''

'''plotting(np.array([flight_dict["L_over_D"], flight_dict["h"]]), "Lift over Drag, L/D [-]", "Altitude, h [m]")

#plotting(np.array([alpha_range, CL_alpha * alpha_range]), "Angle of attack, a [deg]", "Lift coefficient, C_L [-]")

plotting(np.array([flight_dict["C_L"], flight_dict["h"]]), "Lift coefficient, C_L [-]", "Altitude, h [m]")

plotting(np.array([np.array(flight_dict["C_L"]) / AC_params["CL_alpha"], flight_dict["h"]]), "Angle of attack, a [deg]", "Altitude, h [m]")

plotting(np.array([flight_dict["V_stall"], flight_dict["h"], flight_dict["V_T"], flight_dict["h"]]),
         "Stall speed, V_stall and True Airspeed, V_TAS [m/s]", "Altitude, h [m]")

plotting(np.array([flight_dict["mach"], flight_dict["h"]]), "Mach Number, M [-]", "Altitude, h [m]")

plotting(np.array([flight_dict["V_ground"], flight_dict["h"]]), "Ground speed, V_ground [m/s]", "Altitude, h [m]")

plotting(np.array([flight_dict["V_descent"], flight_dict["h"]]), "Descent Speed, V_descent [m/s]", "Altitude, h [m]")

plotting(np.array([flight_dict["distance_travelled"], flight_dict["h"]]), "Distance, D [m]", "Altitude, h [m]")

plotting(np.array([flight_dict["reynolds_number"], flight_dict["h"]]), "Reynolds Number, Re [...]", "Altitude, h [m]")'''


