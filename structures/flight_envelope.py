import numpy as np
import matplotlib.pyplot as plt

rho_0 = 1.225
g0 = 9.8665

v_from_CL = lambda WoS, rho, CL: np.sqrt((WoS)*(2/rho)*(1/CL))

def flight_envelope(CL_max, CL_opt, W,S, v_dive, n_min=-1.5, n_max=2.5):
	v_s_n1 = v_from_CL(W/S, rho_0, CL_max)
	v_glide = v_from_CL(W/S, rho_0, CL_opt)
	#v_dive = v_from_CL(W/S, rho_0, 0.2)
	print(v_s_n1)

	vv = np.arange(v_s_n1*0.7, v_dive*1.1, 0.1)

	fig, ax = plt.subplots()
	ax.set_ylim(n_min*1.5, n_max*1.5)

	ax.plot(vv, (vv/v_s_n1)**2, color='k')
	ax.plot(vv, -(vv/v_s_n1)**2, color='k')
	ax.plot((v_s_n1*n_max**0.5, v_dive), (n_max, n_max), color='k')
	ax.plot((v_s_n1*(-n_min)**0.5, v_glide), (n_min, n_min), color='k')
	ax.plot((v_glide, v_dive),(n_min,0), color='k')
	ax.plot((v_dive,v_dive),(0,n_max), color='k')
	ax.plot((v_s_n1,v_s_n1),(1,-1), color='k')
	ax.set_ylabel("load factor n [-]")
	ax.set_xlabel("IAS [m/s]")


	
	plt.show()








if __name__=="__main__":
	W = 2.5*g0
	S = 0.25
	CL_max = 1.5
	CL_opt = 0.5
	v_dive = 20
	flight_envelope(CL_max, CL_opt, W, S, v_dive)






