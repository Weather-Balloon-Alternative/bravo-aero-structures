import numpy as np
from airfoils import Airfoil
import matplotlib.pyplot as plt
import scipy.integrate





def get_aero_loads(y, W, b, c, n, Ar, e, CLmax, Cmac):
	#Calculate the loads at a location on the wing based on a given weight and load factor
	#for now only at root
	L_center = b/6 #rough estimate for elliptical distr
	L = W*n
	D_center = b/4 #rough estimate
	D = (L**2)/(np.pi*Ar*e) #estimate drag based on current lift, ignore cd0
	M = (Cmac/CLmax)*c


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


def planform(b, Ar, tapre_rat):
	#calculate wing platform properties based on span, aspect ratio and tapre ratio
	S = (b**2)/Ar
	mac = b/Ar
	c_root = (2*S)/((1+tapre_rat)*b)
	c_tip = tapre_rat*c_root
	return {"S":S, "mac":mac, "c_root":c_root, "c_tip":c_tip, "span":b, "aspect":Ar, "tapre":tapre_rat}




	
if __name__=="__main__":
	Wt = 25
	bt = 1.5
	ARt = 8
	tapre_ratt = 0.6
	planform_prop = planform(bt, ARt, tapre_ratt)
	print(planform_prop)

	nt = 2.5
	d_max_t = 7e-3

	mat_prop_cfrp = [650e6]
	res = spar_sizing(Wt, bt, nt, d_max_t, mat_prop_cfrp)
	print(res[0]*1000, res[1]*1000)
	gen_airfoil_shape(1, 0.16, plot=True)
	print(wing_volume(planform_prop['c_root'], tapre_ratt, bt)*25*1000)