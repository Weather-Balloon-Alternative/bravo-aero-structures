import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

def cast_to_array(y, array_shape):
	if type(y) in (float, int):
		y = np.ones(array_shape)*y
	return y

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


	for y in (V, M, I, E):
		#cast any constant input to an array with the same dimension as zz
		y = cast_to_array(y, zz.shape)

	#integrate w to get V, add to the point loads given
	V_tot = scipy.integrate.cumulative_trapezoid(w, zz, initial=0)
	V_tot -= np.cumsum(V)
	print(np.cumsum(V))
	V_tot -= V_tot[-1] #bnd condition: force at tip = 0

	#integrate V to get M, add to the moment given
	M_tot = scipy.integrate.cumulative_trapezoid(V_tot, zz, initial=0) 
	M_tot -= np.cumsum(M)
	M_tot -= M_tot[-1] #bnd condition: moment at tip = 0

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
		ax1.plot(zz, w)
		ax1.plot(zz, V_tot)
		ax1.plot(zz, M_tot)
		#ax1.plot(zz, II*1e11)
		ax1.grid()

		ax2.plot(zz, theta)
		ax2.plot(zz, v)
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













