# Import from Aerodynamics and Structures
import aeroloader
import ambiance
import numpy as np

coeff = aeroloader.loadaero('Aero_output/glider_v7.xlsx')
coeff_notail = aeroloader.loadaero('Aero_output/glider_v7_wingonly.xlsx')

# From requirements:
SM = 0.05 # Stability margin
h = 27000 # height
atmos = ambiance.Atmosphere(h)
rho = atmos.density # density
g0 = atmos.grav_accel # gravitational acceleration
V0 = 150 # velocity TODO: get actual figure from flight planning

# From statistics:
Vh = 0.5    # Horizontal tail volume TODO: From Airplane Design (Sadraey)
lhc = 4+0.008796802253736823     # Ratio of the tail length over the chord TODO: From Airplane Design
ShS=0.10

# From Structures:
deltaxcg = 0.1
xcgbar = 0.25+0.008796802253736823  # location of the centre of gravity divided by the MAC TODO: Guesstimate, get from Florian and Marten
xcg_fwBAR = xcgbar-0.5*deltaxcg  # most forward location of the centre of gravity divided by the MAC TODO: Guesstimate, get from Florian and Marten
xcg_aftBAR = xcgbar+0.5*deltaxcg  # most aft location of the centre of gravity divided by the MAC TODO: Guesstimate, get from Florian and Marten
m = 19 / g0 # mass
Ix = 0 # MMOI around x-axis
Iy = 0 # MMOI around y-axis
Iz = 0 # MMOI around z-axis

# From Aerodynamics:
## Geometry
b = 0.767928166   # Span TODO: Initial guess from Marten
AR = 8.532535173     # Aspect ratio of the main wing TODO: Initial guess from Marten
S = b**2/AR    # Wing surface 
ARh = 4     # Aspect ratio of the horizontal tail TODO: Very much taken from the A-10, can change
xacbar = 0.137686198  # location of the aerodynamic centre divided by the MAC TODO: From Marten
c = S/b       # MAC
lh = lhc*c  # Tail length TODO: Ideally this would be obtained from iteration

## Stability stuff for symmetric
muc = m/(rho*S*c) #relative density [=m/(rho*S*c)]
CZadot = 0 #derivative of Cz w.r.t. alpha_dot*c/V0
Cmadot = 0 #derivative of Cm w.r.t. alpha_dot*c/V0
KY2 = 0 #square of the nondimensional radius of gyration about the Y-axis [=I_y/(m*b)]
CXu = 0 #derivative of X w.r.t. u divided by 1/(0.5*rho*V*S)
CXa = 0 #derivative of CX w.r.t. alpha
CZ0 = 0 #CZ in CX w.r.t. qc/V0
CZu = 0 #derivative of Z w.r.t. u divided by 1/(0.5*rho*V*S)
CZa = 0 #derivative of CZ w.r.t. alpha
CX0 = 0 #CX in steady flight
CZq = 0 #derivative of CZ w.r.t. qc/V0
Cmu = 0 #derivative of M w.r.t. divided by 1/(0.5*rho*V*S)
Cma = 0 #derivative of Cm w.r.t. alpha
Cmq = 0 #derivative of Cm w.r.t. qc/V0
CXde = 0 #derivative of CX w.r.t. delta_e
CZde = 0 #derivative of CZ w.r.t. delta_e
Cmde = 0 #elevator efficiency, derivative of Cm w.r.t. delta_e

## Stability stuff for asymmetric
CYbdot = 0 #derivative of CY w.r.t. beta_dot
mub = 0 #relative density [=m/(rho*S*b)]
b = 0 #wing span
KX2 = 0 #square of the nondimensional radius of gyration about the X-axis [=I_x/(m*b)]
KXZ = 0 #nondimensional product of inertia [=J_xz/(m*b)]
KZ2 = 0 #square of the nondimensional radius of gyration about the Z-axis [=I_z/(m*b)]
CYb = 0 #derivative of CY w.r.t. beta
CL = coeff['CL'][-1] #lift coefficient
CYp = 0 #derivative of CY w.r.t. p*b/(2*V0)
CYr = 0 #derivative of CY w.r.t. r*b/(2*V0)
Clb = 0 #derivative of Cl w.r.t. beta
Clp = 0 #derivative of Cl w.r.t. p*b/(2*V0)
Clr = 0 #derivative of Cl w.r.t. r*b/(2*V0)
Cnb = 0 #static directional stability, derivative of Cn w.r.t. beta
Cnp = 0 #derivative of Cn w.r.t. p*b/(2*V0)
Cnr = 0 #derivative of Cn w.r.t. r*b/(2*V0)
CYda = 0 #derivative of CY w.r.t. delta_a
CYdr = 0 #derivative of CY w.r.t. delta_r
Clda = 0 #derivative of Cl w.r.t. delta_a
Cldr = 0 #derivative of Cl w.r.t. delta_r
Cnda = 0 #derivative of Cn w.r.t. delta_a
Cndr = 0 #derivative of Cn w.r.t. delta_r

## Misc.
deda = 1-0.725    # derivative of the downwash angle w.r.t. alpha TODO: Estimated from Sailplane Design (Thomas)
VhV = 1     # ratio of the velocity at the tail over the aircraft velocity TODO: Guesstimate, just neglect for now
CLAh = coeff_notail['CL'][-1]   # lift coefficient of the aircraft without the tail
CLh = (CL-CLAh)/(VhV**2*ShS)    # lift coefficient of the tail
CLa = (coeff['CL1'][-1] - CL)*180/np.pi # derivative of the lift coefficient of the aircraft
CLaAh = (coeff_notail['CL1'][-1] - CL)*180/np.pi   # derivative of the lift coefficient of the aircraft without the tail w.r.t alpha 
# CLah = (CLa - CLaAh) / (VhV**2*ShS*(1-deda))    # derivative of the lift coefficient of the tail w.r.t. alpha 
CLah = 2 * np.pi * ARh / (ARh + 2)    # derivative of the lift coefficient of the tail w.r.t. alpha 
Cmac = coeff['CMy'][-1] - CLAh*(xcgbar-xacbar)/c + (ShS*lhc*VhV**2)*CLh    # moment coefficient of the aerodynamic centre 

# ###### TEST DATA ######
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