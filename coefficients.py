# Import from Aerodynamics and Structures
import aeroloader

coeff_dict = aeroloader.loadaero('Aero_output/glider_v7.xlsx')
# From requirements:
SM = 0.05 # Stability margin

# From statistics:
Vh = 0.5    # Horizontal tail volume TODO: From Airplane Design (Sadraey)
lhc = 4     # Ratio of the tail length over the chord TODO: From Airplane Design
ShS=0.15

# From Structures:
xcgbar = 0.25  # location of the centre of gravity divided by the MAC TODO: Guesstimate, get from Florian and Marten
xcg_fwBAR = 0.2  # most forward location of the centre of gravity divided by the MAC TODO: Guesstimate, get from Florian and Marten
xcg_aftBAR = 0.3  # most aft location of the centre of gravity divided by the MAC TODO: Guesstimate, get from Florian and Marten

# From Aerodynamics:
## Geometry
b = 0.767928166   # Span TODO: Initial guess from Marten
AR = 8.532535173     # Aspect ratio of the main wing TODO: Initial guess from Marten
S = b**2/AR    # Wing surface 
ARh = 4     # Aspect ratio of the horizontal tail TODO: Very much taken from the A-10, can change
xacbar = 0.137686198  # location of the aerodynamic centre divided by the MAC TODO: From Marten
c = S/b       # MAC
lh = lhc*c  # Tail length TODO: Ideally this would be obtained from iteration

print(coeff_dict["Cmx"])
## Stability stuff for symmetric
muc = 0 #relative density [=m/(rho*S*c)]
V0 = 0 #aircraft speed
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
CL = 0 #lift coefficient
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
CLa = 4.946591 # derivative of the lift coefficient of the aircraft
ClaAh = 5.73   # derivative of the lift coefficient of the aircraft without the tail w.r.t alpha TODO: FILL IN
Clah = -5.73    # derivative of the lift coefficient of the tail w.r.t. alpha TODO: In Rad, from Sailplane Design (Thomas)
deda = 1-0.725    # derivative of the downwash angle w.r.t. alpha TODO: Estimated from Sailplane Design (Thomas)
VhV = 1     # ratio of the velocity at the tail over the aircraft velocity TODO: Guesstimate, just neglect for now
CLh = -0.35*ARh**(1/3)     # lift coefficient of the tail
CLAh = CL - VhV**2*ShS*CLh    # lift coefficient of the aircraft without the tail
Cmac = -0.015    # moment coefficient of the aerodynamic centre TODO: Estimate from Marten



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