import numpy as np
from airfoils import Airfoil
import matplotlib.pyplot as plt
import scipy.integrate


def I_circ(d):
	return (1/64)*np.pi*d**4

def get_aero_loads(y, W, b, c, n, Ar, e, CLmax, Cmac):
	#Calculate the loads at a location on the wing based on a given weight and load factor
	#for now only at root
	L_center = b/6 #rough estimate for elliptical distr
	L = W*n
	D_center = b/4 #rough estimate
	D = (L**2)/(np.pi*Ar*e) #estimate drag based on current lift, ignore cd0
	M = (Cmac/CLmax)*c


def gen_airfoil_shape(c, tc, plot=False):
	#generate an airfoil profile based on a scaled naca 2412
	points = np.genfromtxt('NACA2412.dat', delimiter=',').T
	n_half = 41
	
	points[0] *= c
	points[1] *= (tc/0.12)

	area = -scipy.integrate.simpson(points[1,:n_half], points[0,:n_half]) - scipy.integrate.simpson(points[1,n_half:], points[0,n_half:]) #idk why the first integral is negative
	print(points[0,:n_half])
	print(area)

	if plot:
		plt.plot(points[0], points[1])
		plt.show()
	#upper points, lower points, cross-sectional area
	return points[:,:n_half], points[:,n_half:], area

def wing_volume(cr, tapre_rat, span, tc=0.12):
	A_root = gen_airfoil_shape(cr, tc)[2]
	A_tip = gen_airfoil_shape(cr*tapre_rat, tc)[2]
	V = span*(A_root+A_tip)*0.5
	return V

	

def spar_sizing(W, b, n, d_max, mat_prop):
	#calculate moment assuming uniform lift distribution
	M_root = W*n*0.5*0.25*b
	I_req = (M_root*d_max*0.5)/mat_prop[0]
	d_inner = ((I_circ(d_max)-I_req)*(64/np.pi))**0.25
	return(d_inner, (d_max-d_inner)*0.5)


	
if __name__=="__main__":
	Wt = 25
	bt = 1.5
	nt = 2.5
	d_max_t = 7e-3

	mat_prop_cfrp = [650e6]
	res = spar_sizing(Wt, bt, nt, d_max_t, mat_prop_cfrp)
	print(res[0]*1000, res[1]*1000)
	gen_airfoil_shape(1, 0.16)
	print(wing_volume(0.3, 0.6, bt)*25*1000)