import numpy as np
import ambiance 
import matplotlib.pyplot as plt

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
    x=np.arange(0.01,1.1,0.01)
    a = CD0**2
    b = 2*CD0*1/(np.pi*AR*e)
    c = 1/(np.pi**2*AR**2*e**2)
    C_L_array =[]
    for i in range(len(rho)):
        y = ((W/S)*(2/rho[i]))**0.5*(a-x**2*(b+3*c*x**2))/(2*x**0.5*(a+b*x**2+c*x**4)**1.5)-2*x**2/(CD0**2+2*CD0*x**2/(np.pi*AR*e)+x**4/(np.pi*AR*e)**2)*V_wind[i]

        if y[-1]>0:
            C_L =x[-1]
            C_L_array.append(float(C_L))
        else:
            C_L= x[np.argwhere(y<0)[0]]
            C_L_array.append(float(C_L))
    print(len(C_L_array))
    V = ((W/S)*(2/rho)*(1/C_L))**0.5
    for i in range(len(V)):
        if V[i]>0.8*SOS[i]:
            V[i]=0.8*SOS[i]
            C_L_array[i] = float((W/S)*(2/rho[i])*(1/V[i]**2))
    C_L_array = np.array(C_L_array)
    return V, C_L_array

def descent_range_and_time_calculation(rho, V, V_wind, C_L_array, CD0, AR, e, W, S, max_control_alt):
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
    C_D = CD0+C_L_array**2/(np.pi*AR*e)
    D = 0.5*rho*V**2*S*C_D
    P_req = D*V
    V_descent = P_req/W
    T= 1/V_descent
    R = ((V**2-V_descent**2)**0.5-V_wind)*T
    Range = sum(R[:max_control_alt+1])
    descent_time = sum(T)
    return Range, descent_time



rho = ambiance.Atmosphere(np.array(range(33_000))).density
SOS = ambiance.Atmosphere(np.array(range(33_000))).speed_of_sound

W = 39.32
S = 0.104
CD0 =0.25 #very preliminary
AR = 8
e = 0.9
x=np.arange(0.01,1.1,0.01)
a = CD0**2
b = 2*CD0*1/(np.pi*AR*e)
c = 1/(np.pi**2*AR**2*e**2)
V_wind = np.ones(33_000)*10       #m/s headwind, negative values tailwind
V, C_L_array = flight_velocity_profile(CD0, AR, e, rho, V_wind, W, S,SOS)
print(descent_range_and_time_calculation(rho, V, V_wind, C_L_array, CD0, AR, e, W, S, 28_000))