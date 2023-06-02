import numpy as np
import ambiance 
import matplotlib.pyplot as plt
import pandas as pd

def flight_velocity_profile(CD0, AR, e, rho, V_wind, W, S, SOS):
    """
    ARGS:
            CD0: float,                 Zero lift drag of the complete body
            AR: float,                  Aspect ratio of the vehicle
            e: float btwn 0 and 1,      oswald efficiency factor 
            rho: array of floats,       denisty at all altitudes of length n
            V_wind: array of floats,    Wind at all altitudes of length n
            W: float,                   Weight of the vehicle
            S: float,                   Sufrace area of the wing of the vehicle
            SOS: array of floats,       Speed of sounds at all altitutes of length n
    OUTPUTS:
            V: array of floats,         optimum velocities at the altitudes length n
            C_L_array: array of floats, lift coefficients belonging to velocities of length n
    """
    #range of operable C_L coefficients:
    x=np.arange(0.01,1.1,0.01)

    #terms that are required for the following equations
    a = CD0**2
    b = 2*CD0*1/(np.pi*AR*e)
    c = 1/(np.pi**2*AR**2*e**2)
    
    #making lift coefficients for every altitutde
    C_L_array =[]
    for i in range(len(rho)):
        CD = CD0 + x**2/(np.pi * AR * e)
        y = (np.sqrt((W/S)*(2/rho[i])*(1/x)) -V_wind[i])/((CD/x)*W)
        max_x = np.argwhere(y == np.max(y))
        C_L_array.append(x[max_x][0][0])
    C_L_array = np.array(C_L_array)
    V = ((W/S)*(2/rho)*(1/C_L_array))**0.5
    
    #update velocity and CL in case of crossing 0.8 speed of sound barrier
    for i in range(len(V)):
        if V[i]>0.8*SOS[i]:
            V[i]=0.8*SOS[i]
            C_L_array[i] = float((W/S)*(2/rho[i])*(1/V[i]**2))
    C_L_array = np.array(C_L_array)
    plt.plot(V,np.arange(len(V))*10)
    plt.plot(V-V_wind, np.arange(len(V))*10)
    plt.xlabel("height [m]")
    plt.ylabel("speed [m/s]")
    return V, C_L_array

def descent_range_and_time_calculation(rho, V, V_wind, C_L_array, CD0, AR, e, W, S, max_control_alt, resolution):
    """
    ARGS:
            rho: array of floats,       denisty at all altitudes of length n
            V: array of floats,         optimum velocities at the altitudes length n
            V_wind: array of floats,    Wind at all altitudes of length n
            C_L_array: array of floats, lift coefficients belonging to velocities of length n
            CD0: float,                 Zero lift drag of the complete body
            AR: float,                  Aspect ratio of the vehicle
            e: float btwn 0 and 1,      oswald efficiency factor 
            W: float,                   Weight of the vehicle
            S: float,                   Sufrace area of the wing of the vehicle
            max_control_alt: int,       maximum control altitude
    OUTPUTS:
            Range: float,               Range of theoretical flight profile
            descent_time: float,        Expected descent time
    
    """
    #calculate CD for all altitudes
    C_D = CD0+C_L_array**2/(np.pi*AR*e)
    #calcualte the drag for all altitudes
    D = 0.5*rho*V**2*S*C_D
    #calculate the power required for all altitudes
    P_req = D*V
    #calculate the required descent velocity for that energy height.
    V_descent = P_req/W
    #calculate the time that will be traversed at these altitudes.
    T= resolution/V_descent
    #calculate the range that can be attained based on this
    R = ((V**2-V_descent**2)**0.5-V_wind)*T
    #sum the steerable range (from which one can argue that they are in control of the vehicle)
    Range = sum(R[:int(max_control_alt/resolution)+1])
    #sum the total descent time
    descent_time = sum(T)
    return Range, descent_time


if __name__ =="__main__":
    rho = ambiance.Atmosphere(np.array(range(33_000))).density
    SOS = ambiance.Atmosphere(np.array(range(33_000))).speed_of_sound

    W = 44.561
    S = 0.21
    CD0 = 0.0229#very preliminary
    AR = 10
    e = 0.9
    x=np.arange(0.01,1.1,0.01)
    a = CD0**2
    b = 2*CD0*1/(np.pi*AR*e)
    c = 1/(np.pi**2*AR**2*e**2)
    df = pd.read_json('atmospheric_characteristics.json')
    V_wind= df.loc["average_windspeed"]
    V_wind = V_wind.to_numpy()
    V_wind = V_wind[20:]        #20 because the first like 15 measurements are not a number
    rho = df.loc["density"]
    rho = rho.to_numpy()
    rho = rho[20:]
    SOS = df.loc["a"]
    SOS = SOS.to_numpy()
    SOS = SOS[20:]
    resolution = int(df.loc["h"].iloc[1]-df.loc["h"].iloc[0])
    # V_wind = np.ones(33_000)*10       #m/s headwind, negative values tailwind
    #calculations for average airspeed
    V, C_L_array = flight_velocity_profile(CD0, AR, e, rho, V_wind, W, S,SOS)
    print(descent_range_and_time_calculation(rho, V, V_wind, C_L_array, CD0, AR, e, W, S, 28_000,resolution))
    plt.title("maximum wind speed conditions")
    plt.legend(("True airspeed", "True groundspeed"))
    plt.show()
