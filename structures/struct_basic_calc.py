import numpy as np
from airfoils import Airfoil
import matplotlib.pyplot as plt
import scipy.integrate
import sections





def get_aero_loads(planform, W, n, Cmac, CLmax):
	#Calculate the loads at a location on the wing based on a given weight and load factor
	#for now only at root
	L_center = planform["span"]*(4/(3*np.pi)) #elliptical lift distr
	L = W*n
	D_center = planform["span"]/4 #rough estimate
	D = (L**2)/(np.pi*planform["aspect"]*planform["efficiency"]) #estimate drag based on current lift, ignore cd0
	M_z = (Cmac/CLmax)*planform["mac"] # directions are based on structural analysis coordinate system
	M_x = L*L_center
	M_y = D*D_center
	return {"M_x":M_x, "M_y":M_y, "M_z":M_z}

def landing_loads(planform, W, net_length, landing_speed, landing_angle):
	a = (landing_speed**2)/(2*net_length)

	F_x = sin(landing_angle)*a*W/9.81
	F_y = cos(landing_angle)*a*W/9.81

	wing_loading_point = 1 # 1 is tip, 0 is root

	M_y = F_x*planform["span"]*0.5*wing_loading_point
	M_x = F_y*planform["span"]*0.5*wing_loading_point

	return {"M_y": planform_prop["span"]*0.5*78.15*0.9*(Wt/9.81)*0.7, "M_x":planform_prop["span"]*0.5*78.15*0.9*(Wt/9.81)*0.7}

	


def gen_airfoil_shape(c, tc, plot=False):
	#generate an airfoil profile based on a scaled naca 2412, just used as a generic shape for now
	points = np.genfromtxt('NACA2412.dat', delimiter=',').T
	n_half = 41
	
	points[0] *= c #scale the points by dimensions given
	points[1] *= (tc/0.12)*c # the 0.12 is the default t/c ratio of the airfoil used

	area = -scipy.integrate.simpson(points[1,:n_half], points[0,:n_half]) - scipy.integrate.simpson(points[1,n_half:], points[0,n_half:])
	#print(points[0,:n_half])
	print(area)

	if plot:
		plt.plot(points[0], points[1])
		plt.grid(visible=True)
		plt.show()
	#all points, upper points, lower points, cross-sectional area
	return points, points[:,:n_half], points[:,n_half:], area

def wing_volume(cr, tapre_rat, span, tc=0.12):
	A_root = gen_airfoil_shape(cr, tc)[3]
	A_tip = gen_airfoil_shape(cr*tapre_rat, tc)[3]
	V = span*(A_root+A_tip)*0.5
	return V


def planform(S, Ar, tapre_rat, tc, e):
	#calculate wing platform properties based on span, aspect ratio and tapre ratio
	b = np.sqrt(Ar*S)
	mac = b/Ar
	c_root = (2*S)/((1+tapre_rat)*b)
	c_tip = tapre_rat*c_root
	return {"S":S, "mac":mac, "c_root":c_root, "c_tip":c_tip, "span":b, "aspect":Ar, "tapre":tapre_rat, "tc":tc, "efficiency":e}

def wing_cs_sizing(planform, loads, f_S):

	a_ellipse = planform["c_root"]*0.5 # ellipse half width
	b_ellipse = planform["tc"]*planform["c_root"]*0.5 # ellipse half height


	mat_prop_skin = (24.1e9, 409e6)
	t_min = 1e-8
	t_max = 1e-2
	wing_skin = sections.ellipse((a_ellipse-t_min, b_ellipse-t_min, a_ellipse, b_ellipse), (0,0), mat_prop_skin)
	print("test", wing_skin.I_xx, wing_skin.I_yy)
	print(f_S)


	'''
	for i in range(20):
		t_mid = (t_min + t_max)*0.5
		
		wing_skin.set_inner_t(t_min)
		f_min = (loads["M_x"]*b_ellipse)/wing_skin.I_xx - wing_skin.sigma_y/f_S

		wing_skin.set_inner_t(t_mid)
		f_mid = (loads["M_x"]*b_ellipse)/wing_skin.I_xx - wing_skin.sigma_y/f_S

		if f_min*f_mid < 0:
			t_max = t_mid
		else:
			t_min = t_mid
	'''
	for i in range(20):
		t_mid = (t_min + t_max)*0.5
		
		wing_skin.set_inner_t(t_min)
		f_min = wing_skin.get_max_stress(loads) - wing_skin.sigma_y/f_S

		wing_skin.set_inner_t(t_mid)
		f_mid = wing_skin.get_max_stress(loads) - wing_skin.sigma_y/f_S

		if f_min*f_mid < 0:
			t_max = t_mid
		else:
			t_min = t_mid

		#print(wing_skin.I_xx, f_mid, t_mid)
	
	wing_skin.set_inner_t(t_max)
	stress_max = wing_skin.get_max_stress(loads)
	print("thickness required: {} mm, max stress: {} MPa".format(t_max*1000, stress_max/1e6))








	
if __name__=="__main__":
	#Sizing values
	Wt = 19.1 
	#WoS = 100000
	nt = 2.5
	#St = Wt/WoS
	St = 0.0787
	ARt = 14
	tct = 0.12
	tapre_ratt = 0.6
	planform_prop = planform(St, ARt, tapre_ratt, tct, 0.85) #generate planform
	print(planform_prop)
	planform_prop = {"S":St, "mac":0.07497, "c_root":0.09372, "c_tip":0.05623214, "span":1.04967, "aspect":ARt, "tapre":tapre_ratt, "tc":tct, "efficiency":0.8}
	print(planform_prop)
	#planform_prop = planf

	loads_t = get_aero_loads(planform_prop, Wt, nt, -0.04, 2) # this could also come directly from openvsb

	
	
	wing_cs_sizing(planform_prop, loads_t, 1.5)

	#gen_airfoil_shape(1, 0.16, plot=True)
	#print(wing_volume(planform_prop['c_root'], tapre_ratt, bt)*25*1000)

