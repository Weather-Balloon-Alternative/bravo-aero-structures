import json
import numpy as np
from scipy.constants import g
from matplotlib import pyplot as plt
from ambiance import Atmosphere


class FlightPerformance:
    def __init__(self, wind_profile=None, AC_parameters=None, h_max=None, h_spiral=None, dh=None):
        self.wind_profile = wind_profile
        self.AR = AC_parameters["AR"]
        self.e = AC_parameters["e"]
        self.m = AC_parameters["m"]
        self.S = AC_parameters["S"]
        self.CD_0 = AC_parameters["CD_0"]
        self.CD_0_base = AC_parameters["CD_0_base"]
        self.CD_0_h = AC_parameters["CD_0_h"]
        self.CL_max = AC_parameters["CL_max"]
        self.CL_alpha = AC_parameters["CL_alpha"]
        self.h_range = np.arange(h_max, h_spiral, -dh)
        self.C_L_lst = np.arange(0.01, 1.5, 0.01)
        self.dh = dh
        self.flight_dict = {"h": [], "C_L": [], "L_over_D": [], "V_gnd_over_V_sink": [],
                            "V_T": [], "V_stall": [], "mach": [], "V_ground": [], "V_descent": [],
                            "D_t": [], "distance_travelled": [], "reynolds_number": []}

        if self.CD_0 == None:
            self.CD_0 = {}
            for h in self.h_range:
                self.CD_0[str(h)] = self.CD_0_base + self.CD_0_h * h


    def get_flightdict(self):
        return self.flight_dict

    def get_windprofie(self):
        return self.wind_profile

    def set_windprofile(self, dict):
        self.wind_profile = dict

    def get_hrange(self):
        return self.h_range

    def set_hrange(self, h_max, h_spiral, dh):
        self.h_range = np.arange(h_max, h_spiral, -dh)
        self.dh = dh

    def m_ac(self):
        return np.sqrt(self.AR * self.S) / self.AR

    def gamma(self, C_L, C_D):
        #print(180 / np.pi * np.arctan(C_L / C_D))
        return np.arctan(C_D / C_L)

    def Lift(self, m, gamma):
        return (m * g) / (np.sin(gamma) * np.tan(gamma) + np.cos(gamma))

    def V_TAS(self, Lift, S, rho, C_L):
        return np.sqrt(Lift / S * 2 / rho * 1 / C_L)

    def Vgnd_over_sink_variable_test(self, C_L, C_D, h):
        gamma = self.gamma(C_L, C_D(C_L, h))
        Lift = self.Lift(self.m, gamma)
        V_TAS = self.V_TAS(Lift, self.S, Atmosphere(h).density[0], C_L)
        return (V_TAS * np.cos(gamma) - self.wind_profile[str(h)]) / (V_TAS * np.sin(gamma))

    def C_D_constant(self, C_L, h):
        return self.CD_0 + C_L ** 2 / (np.pi * self.AR * self.e)

    def Vgnd_over_sink_constant(self, C_L, C_D, h):
        V_TAS = np.sqrt((2 * self.m * g) / (Atmosphere(h).density[0] * self.S * C_L))
        return (V_TAS - self.wind_profile[str(h)]) / (V_TAS * (C_D(C_L, h) / C_L))

    def C_D_variable(self, C_L, h):
        return self.CD_0[str(h)] + C_L ** 2 / (np.pi * self.AR * self.e)

    def Vgnd_over_sink_variable(self, C_L, C_D, h):
        V_TAS = np.sqrt((2 * self.m * g) / (Atmosphere(h).density[0] * self.S * C_L))
        return (V_TAS - self.wind_profile[str(h)]) / (V_TAS * (C_D(C_L, h) / C_L))

    def reynolds(self, h, V_TAS):
        return (Atmosphere(h).density[0] * V_TAS * self.m_ac()) / Atmosphere(h).dynamic_viscosity[0]

    def optimal_flight(self, Vgnd_over_sink, C_D, h, plot=False):

        Vgnd_over_sink_h = Vgnd_over_sink(self.C_L_lst, C_D, h)
        Vgnd_over_sink_opt = np.max(Vgnd_over_sink_h)
        index = np.where(Vgnd_over_sink_h == Vgnd_over_sink_opt)
        C_L_opt = float(self.C_L_lst[index[0]])
        L_over_D_opt = C_L_opt / C_D(C_L_opt, h)
        V_TAS_opt = np.sqrt((self.m * g) / self.S
                            * 2 / Atmosphere(h).density[0]
                            * 1 / C_L_opt)

        if plot == True:
            plt.plot(self.C_L_lst, Vgnd_over_sink_h)

        return C_L_opt, L_over_D_opt, V_TAS_opt, Vgnd_over_sink_opt

    def flight_sim(self, plot=False):

        if type(self.CD_0) == dict:
            Vgnd_over_sink = self.Vgnd_over_sink_variable_test
            C_D = self.C_D_variable
        else:
            Vgnd_over_sink = self.Vgnd_over_sink_constant
            C_D = self.C_D_constant

        distance = 0
        V_TAS = 0
        trigger = False

        for i_h, h in enumerate(self.h_range):

            C_L_opt, L_over_D_opt, V_TAS_opt, Vgnd_over_sink_opt = self.optimal_flight(Vgnd_over_sink, C_D, h, plot)

            if (V_TAS < V_TAS_opt) and (trigger is False):

                C_L = 0
                L_over_D = 0
                V_TAS = np.sqrt(2 * g * (self.h_range[0] - h))
                V_ground = 0
                V_descent = V_TAS
                distance = 0

                if V_descent != 0:

                    D_t = self.dh / V_descent

                else:

                    D_t = 0

            else:

                trigger = True
                C_L = C_L_opt
                L_over_D = L_over_D_opt
                V_TAS = V_TAS_opt
                V_ground = V_TAS - self.wind_profile[str(h)]
                V_descent = V_TAS / L_over_D_opt
                D_t = self.dh / V_descent
                dist_at_alt = self.dh * Vgnd_over_sink_opt
                distance += dist_at_alt

            self.flight_dict["h"].append(h)
            self.flight_dict["C_L"].append(C_L)
            self.flight_dict["L_over_D"].append(L_over_D)
            self.flight_dict["V_gnd_over_V_sink"].append(Vgnd_over_sink_opt)
            self.flight_dict["V_T"].append(V_TAS)
            self.flight_dict["mach"].append(V_TAS / Atmosphere(h).speed_of_sound[0])
            self.flight_dict["V_ground"].append(V_ground)
            self.flight_dict["V_descent"].append(V_descent)
            self.flight_dict["D_t"].append(D_t)
            self.flight_dict["reynolds_number"].append(self.reynolds(h, V_TAS))
            self.flight_dict["distance_travelled"].append(distance)
            self.flight_dict["V_stall"].append(np.sqrt((self.m * g) / self.S
                            * 2 / Atmosphere(h).density[0]
                            * 1 / self.CL_max))

        if plot == True:
            plt.xlabel("Lift Coefficient, C_L [-]")
            plt.ylabel("Ground Speed over Descent Speed, V_gnd/h_dot [-]")
            plt.show()


