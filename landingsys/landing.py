import numpy as np
import matplotlib.pyplot as plt



class NetDesign:

    def __init__(self):
        #general stuff
        self.max_decel = 50 #m/s^2
        self.land_speed_mult = 1.15
        self.max_nav_error = 2 #m
        self.nav_error_safety = 1.3
        self.approach_angle = np.deg2rad(20) #deg->rad
        self.ground_clearance = 0.5 #m

        #define stuff small glider
        self.root_chord = 0.0807 #m
        self.wing_location = 0.25 #m from nose
        self.mass_ac = 0.8 #kg
        self.fuselage_width = 0.105 #m
        self.wing_span = 0.775 #m
        self.stall_speed = 14.6 #m/s
        self.acceptable_angle = np.deg2rad(30) #deg->rad

        #define stuff big glider
        self.root_chord_xl = 0.167 #m
        self.wing_location_xl = 0.234 #m from nose
        self.mass_ac_xl = 3.98 #kg
        self.fuselage_width_xl = 0.14 #m
        self.wing_span_xl = 1.6 #m
        self.stall_speed_xl = 14.6 #m/s
        self.acceptable_angle_xl = np.deg2rad(30) #deg->rad

        #calculate stuff
        self.fuselage_nose_length = self.wing_location - 0.5*self.root_chord #m
        self.angle_nose_fuse = np.arctan(2*self.fuselage_nose_length / self.wing_span) #rad
        self.land_speed = self.stall_speed * self.land_speed_mult

        self.fuselage_nose_length_xl = self.wing_location_xl - 0.5*self.root_chord_xl #m
        self.angle_nose_fuse_xl = np.arctan(2*self.fuselage_nose_length_xl / self.wing_span_xl) #rad
        self.land_speed_xl = self.stall_speed_xl * self.land_speed_mult

        #execute methods
        self.calculate_net_radius()


    def calculate_net_radius(self):

        #radious based on deceleration:
        r1 = (self.land_speed*self.land_speed) / 2*self.max_decel
        r1_xl = (self.land_speed_xl*self.land_speed_xl) / 2*self.max_decel

        #radius based on nav error
        r2 = self.max_nav_error / np.tan(max(self.angle_nose_fuse, self.acceptable_angle))
        r2_xl = self.max_nav_error / np.tan(max(self.angle_nose_fuse_xl, self.acceptable_angle_xl))

        #radius based on wind gusts
        #TODO: calculate here
        r3 = 0
        r3_xl = 0

        
        #calc net dimensions
        self.net_top_radius = max(r1, r2, r3, r1_xl, r2_xl, r3_xl)
        self.net_height = (self.max_nav_error * 2) * self.nav_error_safety
        self.net_btm_radius = self.net_top_radius + self.net_height / np.tan(np.deg2rad(90) - self.approach_angle)


    def simulate_impact(self):
        pass


def calc_dy(l, d, dx):
	dl = l-d
	return (l - ((d**2 - dx**2)**0.5 + (dl**2 - dx**2)**0.5))


d = np.arange(1.5,3.5,0.01)
l = np.ones_like(d)*5
dx = np.ones_like(d)
dy = calc_dy(l,d,dx)
print(d)
print(dy)
plt.plot(d,dy)
plt.show()

nd = NetDesign()
print(nd.net_btm_radius, nd.net_top_radius)