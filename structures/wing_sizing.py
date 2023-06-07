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
	M_z = (Cmac/CLmax)*planform["mac"]*n*W # directions are based on structural analysis coordinate system
	M_x = L*L_center
	M_y = D*D_center
	return {"M_x":M_x, "M_y":M_y, "M_z":M_z}

def landing_loads(planform, d_cg_n, W, net_length, landing_speed, incidence_angle):
	a = (landing_speed**2)/(2*net_length)
	F_landing_tot = a*W/9.81
	spanwise_angle = np.arctan2(d_cg_n, planform['span']/2)
	b_wtip = np.cos(spanwise_angle)*0.5*planform['span']
	b_n = np.sin(spanwise_angle)*d_cg_n
	F_wtip = F_landing_tot*(b_n/b_wtip)/(1+(b_n/b_wtip))

	F_x = np.sin(incidence_angle)*np.cos(spanwise_angle)*F_wtip
	F_y = np.cos(incidence_angle)*np.cos(spanwise_angle)*F_wtip


	M_y = F_x*planform["span"]*0.5
	M_x = F_y*planform["span"]*0.5
	print(F_landing_tot, b_n, b_wtip, F_wtip, F_x, F_y, "hot sex")

	return {"F_x": F_x, "F_y": F_y, "M_x": M_x, "M_y": M_y}

	





def planform(S, Ar, taper_rat, tc, e):
	#calculate wing platform properties based on span, aspect ratio and taper ratio
	b = np.sqrt(Ar*S)
	mac = b/Ar
	c_root = (2*S)/((1+taper_rat)*b)
	c_tip = taper_rat*c_root
	return {"S":S, "mac":mac, "c_root":c_root, "c_tip":c_tip, "span":b, "aspect":Ar, "taper":taper_rat, "tc":tc, "efficiency":e}

def wing_cs_sizing(planform, loads, f_S):

	a_ellipse = planform["c_root"]*0.5 # ellipse half width
	b_ellipse = planform["tc"]*planform["c_root"]*0.5 # ellipse half height


	mat_prop_skin = (24.1e9, 10e9, 500e6)

	for load_case in loads:
		t_min = 1e-8
		t_max = 1e-2
		wing_skin = sections.ellipse((a_ellipse-t_min, b_ellipse-t_min, a_ellipse, b_ellipse), (0,0), mat_prop_skin)
		print("test", wing_skin.I_xx, wing_skin.I_yy)
		print(f_S)
		final_thickness = 0
		
		for i in range(20):
			t_mid = (t_min + t_max)*0.5
			
			wing_skin.set_inner_t(t_min)
			f_min = wing_skin.get_max_stress(load_case) - wing_skin.sigma_y/f_S

			wing_skin.set_inner_t(t_mid)
			f_mid = wing_skin.get_max_stress(load_case) - wing_skin.sigma_y/f_S

			if f_min*f_mid < 0:
				t_max = t_mid
			else:
				t_min = t_mid

			#print(wing_skin.I_xx, f_mid, t_mid)
		
		wing_skin.set_inner_t(t_max)
		if t_max > final_thickness:
			final_thickness = t_max

	stress_max = wing_skin.get_max_stress(load_case)
	print("thickness required: {} mm, max stress: {} MPa".format(final_thickness*1000, stress_max/1e6))

	print(loads[1]["M_z"])
	print(wing_skin.get_torsion_twist(loads[1])*57.3)
	print(wing_skin.get_torsion_stress(loads[1])/1e6)








	
if __name__=="__main__":
	#Sizing values
	Wt = 19.1 
	nt = 2.5
	St = 0.0787
	ARt = 14
	tct = 0.12
	taper_ratt = 0.6
	d_cg_nose = 0.165

	planform_prop = planform(St, ARt, taper_ratt, tct, 0.85) #generate planform
	print(planform_prop)

	#planform_prop = planf

	loads_flight = get_aero_loads(planform_prop, Wt, nt, -0.06, 1.3) # this could also come directly from openvsb
	loads_landing = landing_loads(planform_prop, d_cg_nose, Wt, 4, 50, 45/57.3)
	
	
	wing_cs_sizing(planform_prop, (loads_landing, loads_flight), 1.5*1.5*1.41) # 1.5 for uts, 1.5 for approximations made

