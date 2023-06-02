# Import from Aerodynamics and Structures

# From requirements:
SM = 0.05 # Stability margin

# From statistics:
Vh = 0.5    # Horizontal tail volume TODO: From Airplane Design (Sadraey)
lhc = 5     # Ratio of the tail length over the chord TODO: From Airplane Design

# From Structures:
xcgbar = 0.5  # location of the centre of gravity divided by the MAC TODO: Guesstimate, get from Florian and Marten
delta_xcgbar = 0.1 # xcg excursion divided by the MAC #TODO: Again guesstimate, pretty concept dependent 

# From Aerodynamics:
AR = 8      # Aspect ratio of the main wing TODO: Initial guess from Marten
ARh = 4     # Aspect ratio of the horizontal tail TODO: Very much taken from the A-10, can change
xacbar = 0.25  # location of the aerodynamic centre divided by the MAC TODO: Like probably???
c = 0.17       # MAC TODO: Guesstimate from Menno and Marten
Clah = -5.73    # derivative of the lift coefficient of the tail w.r.t. alpha TODO: In Rad, from Sailplane Design (Thomas)
ClaAh = 5.73   # derivative of lift coefficient of the aircraft without the tail w.r.t alpha TODO: FILL IN
deda = 1-0.725    # derivative of the downwash angle w.r.t. alpha TODO: Estimated from Sailplane Design (Thomas)
VhV = 1     # ratio of the velocity at the tail over the aircraft velocity TODO: Guesstimate, just neglect for now
CLh = -0.35*ARh**(1/3)     # lift coefficient of the tail
CLAh = 0.4 - VhV    # lift coefficient of the aircraft without the tail TODO: Guesstimate
Cmac = 0    # moment coefficient of the aerodynamic centre TODO: FILL IN
lh = lhc*c  # Tail length TODO: Ideally this would be obtained from iteration