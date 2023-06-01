import ambiance
import numpy as np

def V_ground(highest_altitude, lowest_altitude,max_mach_no):
    """
    ARGS: 
        highest_altitude: float,            The highest service altitude of the vehicle
        lowest_altitude: float,             Lowest service altitude-> landing height
        max_mach_no: float btwn 0 and 0.8,  maximum mach number for stall.
    OUTPUTS:
        V_ground: float,                    Speed at which the aircraft will stall at ground level
    """
    #calculate the stall speed at maximum service ceilling
    V_max = ambiance.Atmosphere(highest_altitude).speed_of_sound*max_mach_no
    #calculate the density at service ceilling
    rho_max = ambiance.Atmosphere(highest_altitude).density
    #calculate the density at lowest altitude
    rho_0 = ambiance.Atmosphere(lowest_altitude).density
    #use relation to get the stall speed at ground level
    V_ground = V_max * (rho_max/rho_0)**0.5
    return V_ground

def span_cord_area(V_ground, W, C_L, AR, lowest_altitude):
    """
    ARGS: 
           V_ground: float,                     Speed at which the aircraft will stall at ground level
           W: float,                            Weight of the airframe
           C_L: float,                          maximum lift coefficient
           AR: float,                           aspect ratio of the wing
           lowest_altitude: float,              lowest altitude in the atmosphere at which landing occurs                     
    OUTPUTS:
            c: float,                           mean aerodynamic cord of the wing
            b: float,                           span of the wing
            S: float,                           Surface area of the wing

    """
    rho_0 = ambiance.Atmosphere(lowest_altitude).density
    #calculate surface area based on inputs
    S = 2*W/(rho_0*V_ground**2*C_L)
    #use aspect ratio for span and cord
    c = (S/AR)**0.5
    b = AR*c
    return c, b, S

def wing_mass(b, S, taper_ratio, rho_material, fill_coefficient, thickness ):
    """
    ARGS: 
            b: float,                           span of the wing
            S: float,                           Surface area of the wing
            taper_ratio: float btwn 0 and 1     taper ratio of the wing
            rho_material: float,                denisty of the filling material
            fill_coefficient; float btwn 0 and 1fill coefficient of the airfoil
            thickness: float,                   thickness of the airfoil
    OUTPUTS:

    """
    #calculate rootcord
    C_r = (2*S/b)/(1+taper_ratio)
    #calculate tipcord
    C_t = C_r*taper_ratio
    #calculate volume following this
    Volume = 2*((b/2)*((C_r+C_t)/2)*fill_coefficient)*thickness
    #calculate mass of the material
    mass = Volume*rho_material
    return mass



if __name__ =="__main__":
    #input variables
    rho_epo = 25 #kg/m^3
    tr = 0.6 
    W_init = 50
    C_L = 1.1
    AR = 10
    highest_altitude = 33_000
    lowest_altitude = 0
    max_mach_no =0.8
    fill_coefficient = 0.70
    thickness = 0.12
    g_0 = 9.81

    #design specific variables
    V_fus = 0.017
    W_electrics = 1.2*g_0
    W_PL = 2.0*g_0                      #TODO: different payload options
    W_tail = 0.20 *g_0                  #TODO: actual weight
    spar_weight_per_meter = 0.4*g_0 #N

    W_fus = rho_epo*V_fus*g_0

    #calculations
    W = W_init 
    V_0 = V_ground(highest_altitude, lowest_altitude, max_mach_no)

    W_last = 0

    while abs(W_last-W)>0.01:
        W_last = W
        c, b, S = span_cord_area(V_0, W, C_L, AR, lowest_altitude)
        W_wing = wing_mass(b, S, tr, rho_epo, fill_coefficient, thickness)*g_0
        W_spar = b*spar_weight_per_meter
        W = W_wing+W_PL+W_fus+W_electrics+W_spar+W_tail
        # print(W, W_last)
    print(f"##################REPORT#################")
    print("total weight:", round(float(W),3), "[N]")
    print("wing weight:", round(float(W_wing),3), "[N]")
    print("spar weight:", round(float(W_spar),3), "[N]")
    print("payload weight: ", round(float(W_PL),3), "[N]")
    print("fuselage weight:", round(float(W_fus),3), "[N]")
    print("tail weight:", round(float(W_tail),3), "[N]")
    print("electronics weight:", round(float(W_electrics),2), "[N]")
    print("surface area: ", round(float(S),3) , "[m^2]")
    print("mean aerodynamic cord",  round(float(c),3), "[m]")
    print("span",  round(float(b),3), "[m]")
    print("landing speed:",round(float(V_0*1.15),3),"[m/s]")
    print("################END REPORT###############")

    #save results:
    # userinput= input("do you want to print the results? (y/n)")
    # print(str(userinput))
    # if "y" in str(userinput) or "Y" in str(userinput):
    #     print("print")
