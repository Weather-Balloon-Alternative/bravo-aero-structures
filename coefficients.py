# Import from Aerodynamics and Structures
import aeroloader
import ambiance
import numpy as np

# Main inputs:
coeff = aeroloader.loadaero         ('Aero_output/v12g/v12steadywithtail.xlsx')
coeff_notail = aeroloader.loadaero  ('Aero_output/v12g/v12steadynotail.xlsx')
# steady = aeroloader.loadaero('Aero_output/v10/steady_analysis_v10-2.xlsx')
Vh = 0.5    # Horizontal tail volume TODO: From Airplane Design (Sadraey)
lhc = 4     # Ratio of the tail length over the chord TODO: From Airplane Design
ShS = 0.10    # Ratio of the horizontal tail area over the main wing area
shift = (-0.033873816880506724-0.028995398560739183-0.024816579620861368-0.021242461015908987-0.01818186826215651-0.01555642039680949-0.01331733179777228-0.011398012948470204-0.009758521811473297-0.008348901658708052-0.0071416850441190505)*0.064549722/4 # CG shift
incidence_h = -0*np.pi/180 # Tail incidence angle
m = 0.743091757 # mass
mainwingloc = 0.121 # Main wing location
xcg = 0.161268648+shift # CG location
xac = 0.16137931
xac_xroottip = 0.030930075 # AC location relative to the main wing root tip
b = 0.774596669 # Span
c = 0.064549722 # MAC
deda = 1-0.725    # derivative of the downwash angle w.r.t. alpha TODO: Estimated from Sailplane Design (Thomas)
VhV = 0.95     # ratio of the velocity at the tail over the aircraft velocity TODO: Guesstimate, just neglect for now


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
delta_a_max = -15*np.pi/180 # rad, maximum aileron deflection
delta_e_max = -15*np.pi/180 # rad, maximum elevator deflection
delta_r_max = 0*np.pi/180 # rad, maximum rudder deflection =0, no rudder

# From Aerodynamics:
## Geometry
AR = b/c     # Aspect ratio of the main wing
S = b**2/AR    # Wing surface 
ARh = 3     # Aspect ratio of the horizontal tail TODO: Very much taken from the A-10, can change
xacbar = (xac_xroottip)/c  # location of the aerodynamic centre divided by the MAC
lh = lhc*c  # Tail length

# From Structures:
deltaxcg = 0.1
xcgbar = (xcg-mainwingloc)/c # location of the centre of gravity divided by the MAC
xcg_fwBAR = xcgbar-0.5*deltaxcg  # most forward location of the centre of gravity divided by the MAC TODO: Guesstimate, get from Florian and Marten
xcg_aftBAR = xcgbar+0.5*deltaxcg  # most aft location of the centre of gravity divided by the MAC TODO: Guesstimate, get from Florian and Marten
Ix = 0.00599954 # MMOI around x-axis
Iy = 0.01031032 # MMOI around y-axis
Iz = 0.014749995 # MMOI around z-axis
nmax = 2.5 # Maximum allowed load factor

## Stability stuff for symmetric
muc = m/(rho*S*c) #relative density [=m/(rho*S*c)]
CZadot = -coeff['CFz_q'][-1]/3 #derivative of Cz w.r.t. alpha_dot*c/V0 #TODO Find a better alternative
Cmadot = coeff['CMm_q'][-1]/3 #derivative of Cm w.r.t. alpha_dot*c/V0 #TODO Find a better alternative
KY2 = Iy/(m*b) #square of the nondimensional radius of gyration about the Y-axis [=I_y/(m*b)]
CXu = -coeff['CFx_U'][-1] #derivative of X w.r.t. u divided by 1/(0.5*rho*V*S)
CXa = -coeff['CFx_Alpha'][-1] #derivative of CX w.r.t. alpha
CXq = -coeff['CFx_q'][-1] #derivative of CX w.r.t. qc/V0
CZ0 = -coeff['CFz_Total'][-1] #CZ in steady flight
CZu = -coeff['CFz_U'][-1] #derivative of Z w.r.t. u divided by 1/(0.5*rho*V*S)
CZa = -coeff['CFz_Alpha'][-1] #derivative of CZ w.r.t. alpha
CX0 = -coeff['CFx_Total'][-1] #CX in steady flight
CZq = -coeff['CFz_q'][-1] #derivative of CZ w.r.t. qc/V0
Cmu = coeff['CMm_U'][-1] #derivative of M w.r.t. u divided by 1/(0.5*rho*V*S)
Cma = coeff['CMm_Alpha'][-1] #derivative of Cm w.r.t. alpha
Cmq = coeff['CMm_q'][-1] #derivative of Cm w.r.t. qc/V0
CXde = 0 #derivative of CX w.r.t. delta_e # Commonly neglected, =0
CNde = 0 #derivative of N w.r.t. delta_e, basically what force the elevator needs to generate
CZde = -0.06200000000021344 #derivative of CZ w.r.t. delta_e
Cmde = lhc*CZde #elevator efficiency, derivative of Cm w.r.t. delta_e

## Stability stuff for asymmetric
CYbdot = coeff['CFy_r'][-1] #derivative of CY w.r.t. beta_dot
mub = m/(rho*S*b) #relative density [=m/(rho*S*b)]
KX2 = Ix/(m*b) #square of the nondimensional radius of gyration about the X-axis [=I_x/(m*b)]
KXZ = 0 #nondimensional product of inertia [=J_xz/(m*b)] #Assume 0
KZ2 = Iz/(m*b) #square of the nondimensional radius of gyration about the Z-axis [=I_z/(m*b)]
CYb = coeff['CFy_Beta'][-1] #derivative of CY w.r.t. beta
CL = coeff['CL_Total'][-1] #lift coefficient
CYp = coeff['CFy_p'][-1] #derivative of CY w.r.t. p*b/(2*V0)
CYr = coeff['CFy_r'][-1] #derivative of CY w.r.t. r*b/(2*V0)
Clb = coeff['CMl_Beta'][-1] #derivative of Cl w.r.t. beta
Clp = coeff['CMl_p'][-1] #derivative of Cl w.r.t. p*b/(2*V0)
Clr = coeff['CMl_r'][-1] #derivative of Cl w.r.t. r*b/(2*V0)
Cnb = coeff['CMn_Beta'][-1] #static directional stability, derivative of Cn w.r.t. beta
Cnp = coeff['CMn_p'][-1] #derivative of Cn w.r.t. p*b/(2*V0)
Cnr = coeff['CMn_r'][-1] #derivative of Cn w.r.t. r*b/(2*V0)
CYda = 0 #derivative of CY w.r.t. delta_a # =0
CYdr = 0 #derivative of CY w.r.t. delta_r # No rudder, =0
Clda = 0 #aileron efficiency, derivative of Cl w.r.t. delta_a #Should be negative
Cldr = 0 #derivative of Cl w.r.t. delta_r # No rudder, =0
Cnda = 0 #derivative of Cn w.r.t. delta_a # Adverse yaw, is positive #Neglect for now
Cndr = 0 #derivative of Cn w.r.t. delta_r # No rudder, =0

## Misc.
CLAh = coeff_notail['CL_Total'][-1]   # lift coefficient of the aircraft without the tail
CLa = coeff['CL_Alpha'][-1] # derivative of the lift coefficient of the aircraft
CLaAh = coeff_notail['CL_Alpha'][-1]   # derivative of the lift coefficient of the aircraft without the tail w.r.t alpha 
CLah = 2 * np.pi * ARh / (ARh + 2)    # derivative of the lift coefficient of the tail w.r.t. alpha 
CLh = (CL-CLAh)/(VhV**2*ShS) + CLah*incidence_h    # lift coefficient of the tail
Cmac = coeff['CMy_Total'][-1] - CLAh*(xcgbar-xacbar) + (ShS*lhc*VhV**2)*CLh    # moment coefficient of the aerodynamic centre 

