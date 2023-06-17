import numpy as np
from airfoils import Airfoil
import matplotlib.pyplot as plt
import scipy.integrate
import sections
import beam


#maybe move this somewhere else
mat_dict = {"AFRP": {"E":30e9, "G":30e9, "sigma_y":300e6},
			"GFRP": {"E":36.065e9, "G":10e0, "sigma_y":500e6},
			"CFRP": {"E":116.6, "G":30, "sigma_y":900e6}}



def get_aero_loads(planform, n, Cmac, CLmax):
	#Calculate the loads at a location on the wing based on a given weight and load factor
	#for now only at root
	L_center = planform["b"]*(4/(3*np.pi))*0.5 #elliptical lift distr
	L = planform["W"]*n
	D_center = planform["b"]/4 #drag is more uniform

	D = L/5

	print(L/D, "ld")
	M_z = 0.05+abs((Cmac/CLmax)*planform["mac"]*n*planform["W"]) # directions are based on structural analysis coordinate system
	M_x = L*L_center*0.5
	M_y = D*D_center*0.5
	#print("M_x", M_x)
	return {"M":(M_x, M_y, M_z), "F":(-D/2, -L/2, 0)}

def landing_loads(planform, decel, incidence_angle):
	a = decel
	F_landing_tot = a*planform["W"]/9.81
	print(f"decel:{a/9.81}g")
	spanwise_angle = np.arctan2(planform["x_cg"], planform["b"]/2)
	b_wtip = np.cos(spanwise_angle)*0.5*planform["b"]
	b_n = np.sin(spanwise_angle)*planform["x_cg"]
	F_wtip = F_landing_tot*(b_n/b_wtip)/(1+(b_n/b_wtip))

	F_x = np.sin(incidence_angle)*np.cos(spanwise_angle)*F_wtip
	F_y = np.cos(incidence_angle)*np.cos(spanwise_angle)*F_wtip

	#print("beta:", spanwise_angle*57.3)
	M_y2 = np.sin(incidence_angle)*planform["y_mac"]*(a/9.81)*0.125
	M_x2 = np.sin(incidence_angle)*planform["y_mac"]*(a/9.81)*0.125

	M_y = F_x*planform["b"]*0.5
	M_x = F_y*planform["b"]*0.5
	#print(F_landing_tot, b_n, b_wtip, F_wtip, F_x, F_y, "landing forces")
	print(M_y2, M_y, "p")
	return {"F":(F_x, F_y, 0), "M":(M_x, M_y, 0)}


def derived_dimensions(inp_dim):
	b = np.sqrt(inp_dim["AR"]*inp_dim["Sw"])
	mac = b/inp_dim["AR"]
	y_mac = (b/6)*((1+2*inp_dim["taper_rat"])/(1+inp_dim["taper_rat"]))
	c_root = (2*inp_dim["Sw"])/((1+inp_dim["taper_rat"])*b)
	c_tip = inp_dim["taper_rat"]*c_root

	S_h = inp_dim["ShS"]*inp_dim["S"]
	b_h = np.sqrt(inp_dim["AR_h"]*S_h)
	c_h = b_h/inp_dim["AR_h"]

	c_v = c_h
	b_v = inp_dim["AR_v"]*c_v
	S_v = c_v*b_v*2 #because we have an H tail(or a U tail i guess)
	add_dim = {"mac":mac, "y_mac":y_mac, "c_root":c_root, "c_tip":c_tip, "b":b,  "S_h":S_h, "b_h":b_h, "c_h":c_h, "c_v":c_v, "b_v": b_v, "S_v": S_v}
	return {**inp_dim, **add_dim}



def wing_cs_sizing(planform, loads, mat_prop_skin, f_S):

	a_ellipse = planform["c_root"]*0.5 # ellipse half width
	b_ellipse = planform["tc"]*planform["c_root"]*0.5 # ellipse half height


	final_thickness = 0
	for load_case in loads:
		t_min = 1e-8
		t_max = 1e-2
		#wing_skin = sections.ellipse((a_ellipse-t_min, b_ellipse-t_min, a_ellipse, b_ellipse), (0,0), mat_prop_skin)
		wing_skin = sections.thinwalled_airfoil((planform["c_root"], t_min), (0,0), mat_prop_skin)

		print(load_case["M"][0])

		for i in range(20):
			t_mid = (t_min + t_max)*0.5
			
			wing_skin.set_inner_t(t_min)
			wing_skin.get_I()
			f_min = wing_skin.get_max_stress(load_case["M"]) - wing_skin.sigma_y/f_S

			wing_skin.set_inner_t(t_mid)
			wing_skin.get_I()
			f_mid = wing_skin.get_max_stress(load_case["M"]) - wing_skin.sigma_y/f_S

			if f_min*f_mid < 0:
				t_max = t_mid
			else:
				t_min = t_mid
		

		wing_skin.set_inner_t(t_max)
		print(t_max, final_thickness, "hot")
		if (t_max >= final_thickness):
			final_thickness = t_max

	stress_max = wing_skin.get_max_stress(load_case["M"])

	shear_from_L = loads[1]["F"][1]/(wing_skin.get_area()*0.4)

	tau_torsion_open = (3*loads[1]["M"][2])/(wing_skin.circumference*t_max**2)
	print(shear_from_L/1e6, "AAAAAA", tau_torsion_open/1e6)
	print(wing_skin.I_xx, wing_skin.I_yy, wing_skin.I_xy)
	return final_thickness, stress_max, wing_skin.get_torsion_stress(loads[1]["M"])


def wing_deflection(planform, n, W, mat_prop_skin, skin_thickness):
	wing_skin = sections.thinwalled_airfoil((planform["c_root"], skin_thickness), (0,0), mat_prop_skin)

	skin_thickness_section = lambda Y: skin_thickness + Y - Y # t is just constant now
	wing = sections.wing_structure(planform, wing_skin, skin_thickness_section)
	
	L_sectional_center = (4*n*W)/(np.pi*planform["b"])
	#print(L_sectional_center)
	lift_distr = lambda Y: L_sectional_center*np.sqrt(1 - (Y**2)/((0.5*planform["b"])**2))
	#lift_distr = lambda Y: 20 + Y - Y 


	YY = np.linspace(0, np.pi/2, 100)
	YY = np.sin(YY)*planform["b"]/2

	V_applied = np.zeros(YY.shape)
	defl_res = beam.deflection(YY, lift_distr(YY), 0, 0, wing_skin.E, wing.get_sectional_properties_distr(YY)[0], plot=True)
	twist_res = beam.twist(YY, lift_distr(YY)*0.001, 0, wing.get_sectional_properties_distr(YY)[2], 1e9, plot=False)

	moments = np.array([defl_res[1], np.zeros(YY.shape)]).T

	defl = defl_res[-1]


	stress_curve = wing.get_max_stress_curve(YY, moments)
	beam.plot_stress_curve(YY, defl, stress_curve/1e6, plot_stress=False, square=True)
	return defl[-1], np.max(stress_curve)

def tail_defl_stress(dim, n, mat_prop_skin, skin_thickness):
	a_ellipse = dim["c_h"]*0.5
	b_ellipse = a_ellipse*dim["tc_tail"]
	tail_section = sections.ellipse((0, 0, a_ellipse, b_ellipse), (0,0), mat_prop_skin)
	tail_section.set_inner_t(skin_thickness)
	
	YY = np.linspace(0, dim["b_h"]/2, 100)
	tail_lift_distr = (dim["ShS"]*dim["W"]*n)/dim["b_h"] #uniform lift distribution on tail

	M = np.zeros(YY.shape)
	#there is a point moment at the tip of the horizontal stabilizer due to the vertical stabilizer, which is mounted at an offset height
	print(((0.5+dim["v_tail_offset"])**2-(0.5-dim["v_tail_offset"])**2))
	M[-1] = dim["b_v"] * (4/(3*np.pi))*n*dim["W"]*(dim["S_v"]/dim["S"]) * ((0.5+dim["v_tail_offset"])**2-(0.5-dim["v_tail_offset"])**2)*0.5 #only half of the total surface per side
	print(M[-1])
	defl_res = beam.deflection(YY, tail_lift_distr, 0, M, tail_section.E, tail_section.get_I()[0], plot=True)
	print(f"max tail deflection: {defl_res[-1][-1]*1000} mm")
	print((np.max(defl_res[1])*a_ellipse)/(tail_section.I_xx)/1e6)
	print(a_ellipse, b_ellipse)
	print("I_xx:", tail_section.I_xx)







	
if __name__=="__main__":
	#Sizing values


	small_glider_dim = {"W"			: 0.7683*9.81,
						"S"			: 0.071448,
						"Sw"		: 0.05,
						"AR"		: 12,
						"tc"		: 0.08,
						"taper_rat"	: 0.6,
						"x_cg"		: 0.16112,
						"ShS"		: 0.08342,
						"AR_h"		: 3,
						"AR_v"		: 1.8,	
						"lh_c"		: 4,
						"v_tail_offset": 0.3,
						"tc_tail"	: 0.1

	}


	nmax = 2.5

	#planform_prop = planform(St, ARt, taper_ratt, tct, 0.85, S_h, b_h, b_v, c_h) #generate planform
	glider_dim = derived_dimensions(small_glider_dim)
	print(glider_dim) 

	#glider_dim = planf

	loads_flight = get_aero_loads(glider_dim, nmax, -0.06, 1.3) # this could also come directly from openvsb
	loads_landing = landing_loads(glider_dim, 140, 80/57.3)
	
	
	min_thickness, stress_max, tau_max = wing_cs_sizing(glider_dim, (loads_landing, loads_flight), mat_dict['AFRP'], 1.5*1.5) # 1.5 for uts, 1.5 for approximations made
	final_thickness = 0.40/1000
	tip_deflection, stress_max_actual = wing_deflection(glider_dim, nmax, glider_dim["W"], mat_dict['AFRP'], final_thickness)


	tail_defl_stress(glider_dim, nmax, mat_dict["AFRP"], 0.27/1000)
	print("#### REPORT ####")
	print(f"Wing root bending moment flight: {loads_flight} Nm")
	print(f"Wing root bending moment landing: {loads_landing} Nm" )
	print(f"min skin thickness req: {min_thickness*1000} mm")
	print(f"max stress: {stress_max/1e6} Mpa, max shear: {tau_max/1e6}")
	print(f"wing tip deflection under nmax: {tip_deflection*1000} mm")
	print(f"max stress at chosen t{stress_max_actual/1e6} MPa")






