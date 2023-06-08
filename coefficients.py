# Import from Aerodynamics and Structures
import aeroloader
import ambiance
import numpy as np

coeff = aeroloader.loadaero('Aero_output/v10/gliderv10withtail.xlsx')
coeff_notail = aeroloader.loadaero('Aero_output/v10/gliderv10notail.xlsx')

# From requirements:
SM = 0.05 # Stability margin
h = 0 # height - at max velocity, h=29000m at min velocity h=1000m
atmos = ambiance.Atmosphere(h)
rho = atmos.density # density
g0 = atmos.grav_accel # gravitational acceleration
V0 = 34 # velocity - Max velocity 150m/s Min velocity 34m/s
theta = 0 # rad
dampratio_phugoid_min = 0.04 # Minimum dampratio of the phugoid eigenmotion From: Sadraey
dampratio_shortperiod_min = 0.3 # Minimum dampratio of the short-period eigenmotion From: Sadraey
dampratio_spiral_min = 0.0 # Minimum dampratio of the spiral eigenmotion From: requirements
dampratio_dutchroll_min = 0.19 # Minimum dampratio of the Dutch roll eigenmotion From: Sadraey
dampratio_aperiodic_min = 0.0 # Minimum dampratio of the aperiodic roll eigenmotion From: requirements
t_to_60deg_bank = 1.3 # Maximum time to achieve a 60 degree bank angle
delta_a_max = 20*np.pi/180 # rad, maximum aileron deflection
delta_e_max = 20*np.pi/180 # rad, maximum elevator deflection


#Shift
shift = 0

# From statistics:
Vh = 0.5    # Horizontal tail volume TODO: From Airplane Design (Sadraey)
lhc = 4     # Ratio of the tail length over the chord TODO: From Airplane Design
ShS=0.05

# From OpenVSP:
# mainwingloc = 0.154
mainwingloc = 0.154
xcg = 0.162942294
xac_xroottip = 0.030930075
# xac_xroottip = 0.0

# Stuff
incidence_h = -2*np.pi/180

# From Aerodynamics:
## Geometry
b = 0.774596669   # Span
c = 0.064549722      # MAC
AR = b/c     # Aspect ratio of the main wing TODO: Initial guess from Marten
S = b**2/AR    # Wing surface 
ARh = 3     # Aspect ratio of the horizontal tail TODO: Very much taken from the A-10, can change
xacbar = (xac_xroottip)/c +0.03259401019870793  # location of the aerodynamic centre divided by the MAC TODO: From Marten
lh = lhc*c  # Tail length TODO: Ideally this would be obtained from iteration

# From Structures:
deltaxcg = 0.1
xcgbar = (xcg-mainwingloc)/c+shift # location of the centre of gravity divided by the MAC TODO: Guesstimate, get from Florian and Marten
xcg_fwBAR = xcgbar-0.5*deltaxcg  # most forward location of the centre of gravity divided by the MAC TODO: Guesstimate, get from Florian and Marten
xcg_aftBAR = xcgbar+0.5*deltaxcg  # most aft location of the centre of gravity divided by the MAC TODO: Guesstimate, get from Florian and Marten
m = 0.92504571 # mass
Ix = 0 # MMOI around x-axis
Iy = 0 # MMOI around y-axis
Iz = 0 # MMOI around z-axis
nmax = 2.5 # Maximum allowed load factor

## Stability stuff for symmetric
muc = m/(rho*S*c) #relative density [=m/(rho*S*c)]
CZadot = 0 #derivative of Cz w.r.t. alpha_dot*c/V0
Cmadot = 0 #derivative of Cm w.r.t. alpha_dot*c/V0
KY2 = 0 #square of the nondimensional radius of gyration about the Y-axis [=I_y/(m*b)]
CXu = 0 #derivative of X w.r.t. u divided by 1/(0.5*rho*V*S)
CXa = (coeff['CFx1'][-1] - coeff['CFx'][-1])*180/np.pi #derivative of CX w.r.t. alpha
CXq = 0 #derivative of CX w.r.t. qc/V0
CZ0 = -m*g0*np.cos(theta)/(0.5*rho*V0**2*S) #CZ in steady flight
CZu = 0 #derivative of Z w.r.t. u divided by 1/(0.5*rho*V*S)
CZa = (coeff['CFz1'][-1] - coeff['CFz'][-1])*180/np.pi #derivative of CZ w.r.t. alpha
CX0 = -m*g0*np.sin(theta)/(0.5*rho*V0**2*S) #CX in steady flight
CZq = 0 #derivative of CZ w.r.t. qc/V0
Cmu = 0 #derivative of M w.r.t. u divided by 1/(0.5*rho*V*S)
Cma = (coeff['CMx1'][-1] - coeff['CMx'][-1])*180/np.pi #derivative of Cm w.r.t. alpha
Cmq = 0 #derivative of Cm w.r.t. qc/V0
CXde = 0 #derivative of CX w.r.t. delta_e
CZde = 0 #derivative of CZ w.r.t. delta_e
Cmde = 0 #elevator efficiency, derivative of Cm w.r.t. delta_e

## Stability stuff for asymmetric
CYbdot = 0 #derivative of CY w.r.t. beta_dot
mub = m/(rho*S*b) #relative density [=m/(rho*S*b)]
KX2 = 0 #square of the nondimensional radius of gyration about the X-axis [=I_x/(m*b)]
KXZ = 0 #nondimensional product of inertia [=J_xz/(m*b)]
KZ2 = 0 #square of the nondimensional radius of gyration about the Z-axis [=I_z/(m*b)]
CYb = (coeff['CFy2'][-1] - coeff['CFy'][-1])*180/np.pi #derivative of CY w.r.t. beta
CL = coeff['CL'][-1] #lift coefficient
CYp = 0 #derivative of CY w.r.t. p*b/(2*V0)
CYr = 0 #derivative of CY w.r.t. r*b/(2*V0)
Clb = (coeff['CMx2'][-1] - coeff['CMx'][-1])*180/np.pi #derivative of Cl w.r.t. beta
Clp = 0 #derivative of Cl w.r.t. p*b/(2*V0)
Clr = 0 #derivative of Cl w.r.t. r*b/(2*V0)
Cnb = (coeff['CMz2'][-1] - coeff['CMz'][-1])*180/np.pi #static directional stability, derivative of Cn w.r.t. beta
Cnp = 0 #derivative of Cn w.r.t. p*b/(2*V0)
Cnr = 0 #derivative of Cn w.r.t. r*b/(2*V0)
CYda = 0 #derivative of CY w.r.t. delta_a
CYdr = 0 #derivative of CY w.r.t. delta_r # No rudder, =0
Clda = 0 #derivative of Cl w.r.t. delta_a
Cldr = 0 #derivative of Cl w.r.t. delta_r # No rudder, =0
Cnda = 0 #derivative of Cn w.r.t. delta_a
Cndr = 0 #derivative of Cn w.r.t. delta_r # No rudder, =0

## Misc.
deda = 1-0.725    # derivative of the downwash angle w.r.t. alpha TODO: Estimated from Sailplane Design (Thomas)
VhV = 0.8     # ratio of the velocity at the tail over the aircraft velocity TODO: Guesstimate, just neglect for now
CLAh = coeff_notail['CL'][-1]   # lift coefficient of the aircraft without the tail
# CLh = (CL-CLAh)/(VhV**2*ShS)    # lift coefficient of the tail
CLa = (coeff['CL1'][-1] - CL)*180/np.pi # derivative of the lift coefficient of the aircraft
CLaAh = (coeff_notail['CL1'][-1] - CL)*180/np.pi   # derivative of the lift coefficient of the aircraft without the tail w.r.t alpha 
# CLah = (CLa - CLaAh) / (VhV**2*ShS*(1-deda))    # derivative of the lift coefficient of the tail w.r.t. alpha 
CLah = 2 * np.pi * ARh / (ARh + 2)    # derivative of the lift coefficient of the tail w.r.t. alpha 
CLh = (CL-CLAh)/(VhV**2*ShS) + CLah*incidence_h    # lift coefficient of the tail
Cmac = coeff['CMy'][-1] - CLAh*(xcgbar-xacbar)/c + (ShS*lhc*VhV**2)*CLh    # moment coefficient of the aerodynamic centre 

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