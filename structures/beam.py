import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
import matplotlib.cm as mpclm
import matplotlib.colors as mpcolor

def cast_to_array(y, array_shape):
	if type(y) in (float, int, np.float64):
		y = np.ones(array_shape)*y
	return y


def moment_curve(zz, w, V, M):
	'''
	Calculate the total moment curve over a beam based on a generic distributed load and point loads

	ARGS:
	zz: array_like, 		1-d array of points along which the deflection is sampled, zz[-1]-zz[0] represents the length of the beam
	w: array_like, 			1-d array of distributed load for every z
	V: float or array_like, 1-d array of point-forces for every z
	M: float orarray_like, 	1-d array of moments for every z

	RETURNS:
	V_tot,
	M_tot,
	'''

	for y in (V, M):
		#cast any constant input to an array with the same dimension as zz
		y = cast_to_array(y, zz.shape)
	w = cast_to_array(w, zz.shape) # putting this one in the loop above didnt work for some reason

	#integrate w to get V, add to the point loads given
	V_tot = scipy.integrate.cumulative_trapezoid(w, zz, initial=0)
	V_tot -= np.cumsum(V)
	V_tot -= V_tot[-1] #bnd condition: force at tip = 0

	#integrate V to get M, add to the moment given
	M_tot = scipy.integrate.cumulative_trapezoid(V_tot, zz, initial=0) 
	M_tot -= np.cumsum(M)
	M_tot -= M_tot[-1] #bnd condition: moment at tip = 0

	return V_tot, M_tot

def deflection(zz, w, V, M, E, I, plot=False):
	'''
	Calculate the deflection of a generic beam with profile EI clamped at x=0

	TODO: implement generic boundary conditions
	
	ARGS:
		zz: array_like, 		1-d array of points along which the deflection is sampled, zz[-1]-zz[0] represents the length of the beam
		w: array_like, 			1-d array of distributed load for every z
		V: float or array_like, 1-d array of point-forces for every z
		M: float orarray_like, 	1-d array of moments for every z
		E: float or array_like, Modulus of elasticity. Can be either a constant float or a 1-d array representing a varying profile
		I: float or array_like,	Area moment of beam along bending axis. Can be either a constant float or a 1-d array representing a varying profile
		plot: bool,				if true, plot the results using matplotlib


	RETURNS:
		V_tot, 
		M_tot,
		theta,
		v: ndarray, 			1-d array of deflection along the beam axis for every z
	'''


	for y in (I, E):
		#cast any constant input to an array with the same dimension as zz
		y = cast_to_array(y, zz.shape)
	

	V_tot, M_tot = moment_curve(zz, w, V, M)
	dtheta_dx = M_tot/(E*I) #for theta the total moment needs to be divided by the flexural rigidity

	#integrate dtheta_dx to get theta
	theta = scipy.integrate.cumulative_trapezoid(dtheta_dx, zz, initial=0)
	theta -= theta[0] # bnd condition: theta at root = 0

	#integrate theta to get the deflection v
	v = scipy.integrate.cumulative_trapezoid(theta, zz, initial=0)
	v -= v[0] #bnd condition: defl at root = 0

	if plot:
		fig1, ax1 = plt.subplots()
		fig2, ax2 = plt.subplots()
		w = cast_to_array(w, zz.shape)
		ax1.plot(zz, w)
		ax1.plot(zz, V_tot)
		ax1.plot(zz, M_tot)
		#ax1.plot(zz, I)
		ax1.grid()

		#ax2.plot(zz, theta)
		ax2.plot(zz, v)
		ax2.set_ylim(ax2.get_xlim()[1]*-0.25,ax2.get_xlim()[1]*0.5)
		ax2.set_xlabel("Dist. along span [m]")
		ax2.set_ylabel("Deflection [m]")
		ax2.set_aspect(1)
		ax2.grid()
		plt.show()

	return V_tot, M_tot, theta, v



def twist(zz, m, T, J, G, plot=False):
	'''
	Calculate the deflection of a generic beam with profile EI clamped at x=0

	TODO: implement generic boundary conditions
	
	ARGS:
		zz: array_like, 		1-d array of points along which the deflection is sampled, zz[-1]-zz[0] represents the length of the beam
		m: float array_like, 	1-d array of distributed torque for every z
		T: float array_like, 	1-d array of point-torques for every z
		J: float array_like, 	1-d array or constant value of polar moment 
		G: float or array_like, 1-d array or constant value of shear modulus
		plot: bool,				if true, plot the results using matplotlib


	RETURNS:
		T_tot, 
		phi: ndarray, 			1-d array of twist angle along the beam axis for every z
	'''

	for y in (m, T, J, G):
		y = cast_to_array(y, zz.shape)

	T_tot = scipy.integrate.cumulative_trapezoid(m, zz, initial=0)
	T_tot -= np.cumsum(T)
	T_tot -= T_tot[-1] #bnd condition: T at free end = 0


	dphi_dx = (T_tot)/(J*G)
	#print(dphi_dx, T_tot)
	phi = scipy.integrate.cumulative_trapezoid(dphi_dx, zz, initial=0)

	if plot:
		fig1, ax1 = plt.subplots()
		fig2, ax2 = plt.subplots()
		ax1.plot(zz, m)
		ax1.plot(zz, T_tot)
		ax1.plot(zz, J)

		#ax2.plot(zz, dphi_dx)
		ax2.plot(zz, phi*57.3)
		plt.show()

	return T_tot, phi

def plot_stress_curve(zz, v, stress_curve, plot_stress=True, square=False):
	norm = mpcolor.Normalize(vmin=0, vmax=np.max(stress_curve))
	cmap = mpclm.inferno
	m = mpclm.ScalarMappable(norm=norm, cmap=cmap)
	fig3, ax3 = plt.subplots()
	fig3.colorbar(m, label="Stress [MPa]")
	if plot_stress:
		ax3.plot(zz, stress_curve/np.max(stress_curve)*0.02, color="C1")




	ax3.set_xlabel("Dist. along span [m]")
	ax3.set_ylabel("Deflection [m]")
	ax3.grid()
	
	for idx, sec in enumerate(zz[:-1]):
		ax3.plot(zz[idx:(idx+2)], v[idx:(idx+2)], color = m.to_rgba(stress_curve[idx]))
	#ax3.set_aspect(1)
	if square:
		ax3.set_aspect(1)
		ax3.set_ylim(ax3.get_xlim()[1]*-0.25,ax3.get_xlim()[1]*0.5)

	plt.show()












