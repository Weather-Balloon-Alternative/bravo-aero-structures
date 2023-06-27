import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

def airfoil_shape(plot=False):
	#generate an airfoil profile based on a scaled naca 2412, just used as a generic shape for now
	points = np.unique(np.genfromtxt('airfoils/NACA2032C.csv', delimiter=',').T,axis=0)
	#area = -scipy.integrate.simpson(points[1,:n_half], points[0,:n_half]) - scipy.integrate.simpson(points[1,n_half:], points[0,n_half:])
	area = -scipy.integrate.simpson(points[1][points[1]>=0], points[0][points[1]>=0]) - scipy.integrate.simpson(points[1][points[1]<=0], points[0][points[1]<=0])
	#print(points[0,:n_half])
	#print(area)

	if plot:
		plt.plot(points[0], points[1])
		plt.grid(visible=True)
		plt.show()
	#points, cross-sectional area
	return points, area

def volume(cr, taper_rat, span, tc=0.12):
	A_root = airfoil_shape(cr, tc)[3]
	A_tip = airfoil_shape(cr*taper_rat, tc)[3]
	V = span*(A_root+A_tip)*0.5
	return V


if __name__=="__main__":
	print("bruh")
	airfoil_shape(0.1)


