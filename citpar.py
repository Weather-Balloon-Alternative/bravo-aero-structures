# Import from Aerodynamics and Structures
import aeroloader
import ambiance
import numpy as np

# Main inputs:
coeff = aeroloader.loadaero('Aero_output/v12/v12steadywithtail_v2.xlsx')
coeff_notail = aeroloader.loadaero('Aero_output/v12/v12steadynotail_v2.xlsx')
# steady = aeroloader.loadaero('Aero_output/v10/steady_analysis_v10-2.xlsx')
Vh =  (0.2 * 30 )*( 0.71 * 5.968)  # Horizontal tail volume TODO: From Airplane Design (Sadraey)
lhc =  0.71 * 5.968     # Ratio of the tail length over the chord TODO: From Airplane Design
ShS =  0.2   # Ratio of the horizontal tail area over the main wing area
shift = 0 # CG shift
incidence_h = -2*np.pi/180 # Tail incidence angle
m = 4547.8 # mass
mainwingloc = 0.154 # Main wing location
xcg =  # CG location TODO: find what 0.30c for location is
xac_xroottip = 0.030930075 # AC location relative to the main wing root tip
b = 15.911
c = 2.0569	# MAC
deda = 1-0.725    # derivative of the downwash angle w.r.t. alpha TODO: Estimated from Sailplane Design (Thomas)
VhV = 1     # ratio of the velocity at the tail over the aircraft velocity TODO: Guesstimate, just neglect for now
th0= 0
W=m*9.81

# From requirements:
SM = 0.05 # Stability margin
h = 0 # height - at max velocity, h=29000m at min velocity h=1000m
atmos = ambiance.Atmosphere(h)
rho = float(atmos.density) # density
g0 = float(atmos.grav_accel) # gravitational acceleration
V0 = 30 # velocity - Max velocity 150m/s Min velocity 34m/s
V0c = V0/c
theta = 0 # rad
dampratio_phugoid_min = 0.04 # Minimum dampratio of the phugoid eigenmotion From: Sadraey
dampratio_shortperiod_min = 0.3 # Minimum dampratio of the short-period eigenmotion From: Sadraey
dampratio_spiral_min = 0.0 # Minimum dampratio of the spiral eigenmotion From: requirements
dampratio_dutchroll_min = 0.19 # Minimum dampratio of the Dutch roll eigenmotion From: Sadraey
dampratio_aperiodic_min = 0.0 # Minimum dampratio of the aperiodic roll eigenmotion From: requirements
t_to_60deg_bank = 1.3 # Maximum time to achieve a 60 degree bank angle
delta_a_max = 20*np.pi/180 # rad, maximum aileron deflection
delta_e_max = 20*np.pi/180 # rad, maximum elevator deflection
delta_r_max = 0*np.pi/180 # rad, maximum rudder deflection =0, no rudder

# From Aerodynamics:
## Geometry
AR = b/c     # Aspect ratio of the main wing
S = 30    # Wing surface 
ARh = 3     # Aspect ratio of the horizontal tail TODO: Very much taken from the A-10, can change
xacbar = (xac_xroottip)/c  # location of the aerodynamic centre divided by the MAC
lh = lhc*c  # Tail length

# From Structures:
deltaxcg = 0.1
xcgbar = (xcg-mainwingloc)/c # location of the centre of gravity divided by the MAC
xcg_fwBAR = xcgbar-0.5*deltaxcg  # most forward location of the centre of gravity divided by the MAC TODO: Guesstimate, get from Florian and Marten
xcg_aftBAR = xcgbar+0.5*deltaxcg  # most aft location of the centre of gravity divided by the MAC TODO: Guesstimate, get from Florian and Marten
# Ix =  # MMOI around x-axis TODO Get actual value
# Iy =  # MMOI around y-axis TODO Get actual value
# Iz =  # MMOI around z-axis TODO Get actual value
KX2    = 0.019
KZ2    = 0.042
KXZ    = 0.002
KY2    = 1.25 * 1.114
nmax = 2.5 # Maximum allowed load factor

## Stability stuff for symmetric
muc = m/(rho*S*c) #relative density [=m/(rho*S*c)]
CZadot = -0.00350 #derivative of Cz w.r.t. alpha_dot*c/V0 #TODO Find a better alternative
Cmadot = +0.17800 #derivative of Cm w.r.t. alpha_dot*c/V0 #TODO Find a better alternative
# KY2 =  #square of the nondimensional radius of gyration about the Y-axis [=I_y/(m*b)]
CXu = -0.09500 #derivative of X w.r.t. u divided by 1/(0.5*rho*V*S)
CXa = +0.47966 #derivative of CX w.r.t. alpha
CXq = -0.28170 #derivative of CX w.r.t. qc/V0
CZ0 = -W * np.cos(th0) / (0.5 * rho * V0 ** 2 * S) #CZ in steady flight
CZu = -0.37616 #derivative of Z w.r.t. u divided by 1/(0.5*rho*V*S)
CZa = -5.74340#derivative of CZ w.r.t. alpha
CX0 =  W * np.sin(th0) / (0.5 * rho * V0 ** 2 * S)#CX in steady flight
CZq =  -5.66290 #derivative of CZ w.r.t. qc/V0
Cmu =  +0.06990 #derivative of M w.r.t. u divided by 1/(0.5*rho*V*S)
Cma =  -0.4300 #derivative of Cm w.r.t. alpha
Cmq =  -8.79415 #derivative of Cm w.r.t. qc/V0
CXde = -0.03728 #derivative of CX w.r.t. delta_e # Commonly neglected, =0
CNde = 0#derivative of N w.r.t. delta_e, basically what force the elevator needs to generate
CZde = -0.69612#derivative of CZ w.r.t. delta_e
Cmde =  -1.5530#elevator efficiency, derivative of Cm w.r.t. delta_e

## Stability stuff for asymmetric
CYbdot = 0#derivative of CY w.r.t. beta_dot
mub = m/(rho*S*b) #relative density [=m/(rho*S*b)]
# KX2 = Ix/(m*b) #square of the nondimensional radius of gyration about the X-axis [=I_x/(m*b)]
# KXZ = 0 #nondimensional product of inertia [=J_xz/(m*b)] #Assume 0
# KZ2 = Iz/(m*b) #square of the nondimensional radius of gyration about the Z-axis [=I_z/(m*b)]
CYb = -0.7500 #derivative of CY w.r.t. beta
CL  =  1.1360 #lift coefficient
CYp = -0.0304 #derivative of CY w.r.t. p*b/(2*V0)
CYr = +0.8495 #derivative of CY w.r.t. r*b/(2*V0)
Clb = -0.10260 #derivative of Cl w.r.t. beta
Clp = -0.71085 #derivative of Cl w.r.t. p*b/(2*V0)
Clr = +0.23760 #derivative of Cl w.r.t. r*b/(2*V0)
Cnb = +0.1348 #static directional stability, derivative of Cn w.r.t. beta
Cnp = -0.0602 #derivative of Cn w.r.t. p*b/(2*V0)
Cnr = -0.2061 #derivative of Cn w.r.t. r*b/(2*V0)
CYda = -0.0400 #derivative of CY w.r.t. delta_a # =0
CYdr = +0.2300 #derivative of CY w.r.t. delta_r # No rudder, =0
Clda = -0.23088 #aileron efficiency, derivative of Cl w.r.t. delta_a #Should be negative
Cldr = +0.03440 #derivative of Cl w.r.t. delta_r # No rudder, =0
Cnda = -0.0120 #derivative of Cn w.r.t. delta_a # Adverse yaw, is positive
Cndr = -0.0939 #derivative of Cn w.r.t. delta_r # No rudder, =0

## Misc.
CLAh  =    # lift coefficient of the aircraft without the tail
CLa   = 5.74340   # derivative of the lift coefficient of the aircraft
CLaAh =    # derivative of the lift coefficient of the aircraft without the tail w.r.t alpha 
CLah  =     # derivative of the lift coefficient of the tail w.r.t. alpha 
CLh   =     # lift coefficient of the tail
Cmac  =     # moment coefficient of the aerodynamic centre 

# ###### TEST DATA #####
# Vh = 0.41339
# lhc = 2.066
# SM = 0.05


# # xcgbar = 0.4
# xcg_fwBAR = 0.1
# xcg_aftBAR = 0.4

# S = 30.00
# b = 15.911
# AR = b**2/S
# ARh = 5.5893
# xacbar = 0.25
# c = 2.0569
# Clah = 2 * 3.14 * ARh / (ARh + 2)
# ClaAh = 5.35 - Clah*Vh
# deda = 0.05
# VhV = 0.8
# CLh = -0.35*ARh**(1/3)
# CLAh = 0.4 - CLh
# Cmac = 0
# lh = lhc*c