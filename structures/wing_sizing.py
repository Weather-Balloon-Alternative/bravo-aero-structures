import numpy as np
from airfoils import Airfoil
import matplotlib.pyplot as plt
import scipy.integrate
import sections
import beam
import matplotlib.cm as mpclm
import matplotlib.colors as mpcolor


#maybe move this somewhere else
mat_dict = {"AFRP": {"E":24.1e9, "G":10e9, "sigma_y":500e6},
			"GFRP": {"E":24.1e9, "G":10e0, "sigma_y":400e6}}



def get_aero_loads(planform, W, n, Cmac, CLmax):
	#Calculate the loads at a location on the wing based on a given weight and load factor
	#for now only at root
	L_center = planform["span"]*(4/(3*np.pi))*0.5 #elliptical lift distr
	L = W*n
	D_center = planform["span"]/4 #rough estimate
	D = (L**2)/(np.pi*planform["aspect"]*planform["efficiency"]) #estimate drag based on current lift, ignore cd0
	M_z = (Cmac/CLmax)*planform["mac"]*n*W # directions are based on structural analysis coordinate system
	M_x = L*L_center*0.5
	M_y = D*D_center
	print("M_x", M_x)
	return {"M_x":M_x, "M_y":M_y, "M_z":M_z}

def landing_loads(planform, d_cg_n, W, net_length, landing_speed, incidence_angle):
	a = (landing_speed**2)/(2*net_length)
	print("boomshakalaka:", a/9.81)
	F_landing_tot = a*W/9.81
	spanwise_angle = np.arctan2(d_cg_n, planform['span']/2)
	b_wtip = np.cos(spanwise_angle)*0.5*planform['span']
	b_n = np.sin(spanwise_angle)*d_cg_n
	F_wtip = F_landing_tot*(b_n/b_wtip)/(1+(b_n/b_wtip))

	F_x = np.sin(incidence_angle)*np.cos(spanwise_angle)*F_wtip
	F_y = np.cos(incidence_angle)*np.cos(spanwise_angle)*F_wtip

	print("beta:", spanwise_angle*57.3)
	M_y = F_x*planform["span"]*0.5
	M_x = F_y*planform["span"]*0.5
	print(F_landing_tot, b_n, b_wtip, F_wtip, F_x, F_y, "landing forces")

	return {"F_x": F_x, "F_y": F_y, "M_x": M_x, "M_y": M_y}

	





def planform(S, Ar, taper_rat, tc, e):
	#calculate wing platform properties based on span, aspect ratio and taper ratio
	b = np.sqrt(Ar*S)
	mac = b/Ar
	c_root = (2*S)/((1+taper_rat)*b)
	c_tip = taper_rat*c_root
	return {"S":S, "mac":mac, "c_root":c_root, "c_tip":c_tip, "span":b, "aspect":Ar, "taper":taper_rat, "tc":tc, "efficiency":e}

def wing_cs_sizing(planform, loads, mat_prop_skin, f_S):

	a_ellipse = planform["c_root"]*0.5 # ellipse half width
	b_ellipse = planform["tc"]*planform["c_root"]*0.5 # ellipse half height



	for load_case in loads:
		t_min = 1e-8
		t_max = 1e-2
		wing_skin = sections.ellipse((a_ellipse-t_min, b_ellipse-t_min, a_ellipse, b_ellipse), (0,0), mat_prop_skin)
		#print("test", wing_skin.I_xx, wing_skin.I_yy)
		#print(f_S)
		final_thickness = 0
		
		moments = (load_case["M_x"], load_case["M_y"], 0)
		print(moments)
		for i in range(20):
			t_mid = (t_min + t_max)*0.5
			
			wing_skin.set_inner_t(t_min)
			f_min = wing_skin.get_max_stress(moments) - wing_skin.sigma_y/f_S

			wing_skin.set_inner_t(t_mid)
			f_mid = wing_skin.get_max_stress(moments) - wing_skin.sigma_y/f_S

			if f_min*f_mid < 0:
				t_max = t_mid
			else:
				t_min = t_mid

			#print(wing_skin.I_xx, f_mid, t_mid)
		
		wing_skin.set_inner_t(t_max)
		if t_max > final_thickness:
			final_thickness = t_max

	stress_max = wing_skin.get_max_stress(moments)
	print("thickness required: {} mm, max stress: {} MPa".format(final_thickness*1000, stress_max/1e6))
	print("Shear stress root due to twisting: {}".format(wing_skin.get_torsion_stress(loads[1])/1e6))


def wing_deflection(planform, n, W, mat_prop_skin, skin_thickness):
	a_ellipse = planform["c_root"]*0.5 # ellipse half width
	b_ellipse = planform["tc"]*planform["c_root"]*0.5 # ellipse half height
	wing_skin = sections.ellipse((0, 0, a_ellipse, b_ellipse), (0,0), mat_prop_skin)

	skin_thickness_section = lambda Y: skin_thickness + skin_thickness*(Y < planform['span']*0.25)
	wing = sections.wing_structure(planform, wing_skin, skin_thickness_section)
	
	L_sectional_center = (4*n*W)/(np.pi*planform['span'])
	print(L_sectional_center)
	lift_distr = lambda Y: L_sectional_center*np.sqrt(1 - (Y**2)/((0.5*planform['span'])**2))
	#lift_distr = lambda Y: 20 + Y - Y 


	YY = np.linspace(0, np.pi/2, 100)
	YY = np.sin(YY)*planform['span']/2

	V_applied = np.zeros(YY.shape)
	defl_res = beam.deflection(YY, lift_distr(YY), 0, 0, wing_skin.E, wing.get_sectional_properties_distr(YY)[0], plot=True)
	twist_res = beam.twist(YY, lift_distr(YY)*0.001, 0, wing.get_sectional_properties_distr(YY)[2], 1e9, plot=True)

	moments = np.array([defl_res[1], np.zeros(YY.shape)]).T

	defl = defl_res[-1]


	stress_curve = wing.get_max_stress_curve(YY, moments)

	norm = mpcolor.Normalize(vmin=-np.max(stress_curve), vmax=np.max(stress_curve))
	cmap = mpclm.coolwarm.reversed()
	m = mpclm.ScalarMappable(norm=norm, cmap=cmap)
	print(defl_res[-1][-1], "defl_max2")
	
	fig3, ax3 = plt.subplots()
	for idx, sec in enumerate(YY[:-1]):
		ax3.plot(YY[idx:(idx+2)], defl[idx:(idx+2)], color = m.to_rgba(stress_curve[idx]))
	#ax3.set_aspect(1)
	plt.show()






	
if __name__=="__main__":
	#Sizing values
	Wt = 19.1 
	nt = 2.5
	St = 0.05
	ARt = 12
	tct = 0.07
	taper_ratt = 0.6
	d_cg_nose = 0.165*2

	planform_prop = planform(St, ARt, taper_ratt, tct, 0.85) #generate planform
	print(planform_prop)

	#planform_prop = planf

	loads_flight = get_aero_loads(planform_prop, Wt, nt, -0.06, 1.3) # this could also come directly from openvsb
	loads_landing = landing_loads(planform_prop, d_cg_nose, Wt, 2.89, 20 , 85/57.3)
	
	
	wing_cs_sizing(planform_prop, (loads_landing, loads_flight), mat_dict['AFRP'], 1.5*1.5*1.41) # 1.5 for uts, 1.5 for approximations made
	wing_deflection(planform_prop, 2.5, Wt, mat_dict['AFRP'], 0.27/1000)

