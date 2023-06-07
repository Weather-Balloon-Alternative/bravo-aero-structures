import numpy as np
import scipy.integrate

def airfoil_shape(c, tc, plot=False):
	#generate an airfoil profile based on a scaled naca 2412, just used as a generic shape for now
	points = np.genfromtxt('NACA2412.dat', delimiter=',').T
	n_half = 41
	
	points[0] *= c #scale the points by dimensions given
	points[1] *= (tc/0.12)*c # the 0.12 is the default t/c ratio of the airfoil used

	area = -scipy.integrate.simpson(points[1,:n_half], points[0,:n_half]) - scipy.integrate.simpson(points[1,n_half:], points[0,n_half:])
	#print(points[0,:n_half])
	#print(area)

	if plot:
		plt.plot(points[0], points[1])
		plt.grid(visible=True)
		plt.show()
	#all points, upper points, lower points, cross-sectional area
	return points, points[:,:n_half], points[:,n_half:], area

def volume(cr, taper_rat, span, tc=0.12):
	A_root = airfoil_shape(cr, tc)[3]
	A_tip = airfoil_shape(cr*taper_rat, tc)[3]
	V = span*(A_root+A_tip)*0.5
	return V




