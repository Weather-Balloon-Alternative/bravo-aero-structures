import numpy as np
import scipy.optimize as spo
import wing_geometry
import matplotlib.pyplot as plt
import matplotlib.cm as mpclm
import matplotlib.colors as mpcolor

class wing_structure():
	def __init__(self, planform, section, section_thickness):
		self.planform = planform
		self.section = section
		self.section_thickness = section_thickness
		self.root_section = self.get_section_at_span(0)
		self.tip_section = self.get_section_at_span(planform["b"]*0.5)

	def get_section_at_span(self, Y):
		#This function would need to be changed if the more generic coimposite_section class is to be used
		
		c_sec = (self.planform["c_root"] - Y*((self.planform["c_root"]-self.planform["c_tip"])/(self.planform["b"]*0.5)))
		#self.section.ao = c_sec*0.5
		#self.section.bo = self.section.ao*self.planform["tc"]
		self.section.resize(c_sec)

		self.section.set_inner_t(self.section_thickness(Y))
		return self.section

	def get_sectional_properties_distr(self, YY):
		II_xx = II_yy = np.zeros(YY.shape)
		JJ = np.zeros(YY.shape)
		for idx, Y in enumerate(YY):
			cur_section = self.get_section_at_span(Y)
			II_xx[idx], II_yy, II_xy = cur_section.get_I()
			JJ[idx] = cur_section.get_J()
		self.II_xx, self.II_yy, self.II_xy, self.JJ = II_xx, II_yy, II_xy, JJ
		return II_xx, II_yy, II_xy, JJ

	def get_max_stress_curve(self, YY, M):
		stress_curve = np.zeros(YY.shape)
		for idx, Y in enumerate(YY):
			cur_section = self.get_section_at_span(Y)
			stress_curve[idx] = cur_section.get_max_stress(M[idx])
		self.stress_curve = stress_curve
		return stress_curve





class composite_section():
	'''class that represents a complete wing section, built up of individual sections defined below'''
	def __init__(self, sections):
		self.sections = sections #list of section classes that the wing is built up of
		return self

	def add_section(self, section):
		#add a singular section to the total structure
		self.sections.append(section)
		return self

	def get_centroid(self):
		num = den = 0
		for section in self.sections: # maybe this could be rewritten to use np broadcasting
			num += section.E*section.A*section.offset
			den += section.E*section.A
		self.centroid = num/den

	def get_curvature_radius(self, Mx, My):
		self.curvature_radius = 0 #TODO
		R_x = R_y = 0
		for section in self.sections:
			R_x += section.E*section.I_xx
			R_y += section.E*section.I_yy

		R_x /= M_x
		R_y /= M_y


class thinwalled_airfoil():
	def __init__(self, geom, offset, mat_prop):
		self.chord, self.wall_t = geom
		self.offset = np.array(offset)
		self.E, self.G, self.sigma_y = mat_prop['E'], mat_prop['G'], mat_prop['sigma_y']
		self.ref_airfoil_points, self.ref_enclosed_area = wing_geometry.airfoil_shape()
		self.resize(self.chord)

	def resize(self, chord):
		self.airfoil_points = self.ref_airfoil_points*chord
		self.enclosed_area = self.ref_enclosed_area*chord**2

	def set_inner_t(self, t):
		#this function is mostly for compatibility with previous code
		self.wall_t = t
		self.get_I()
		self.get_J() 

	def get_area(self):
		shifted_points = np.roll(self.airfoil_points, 1, axis=1)
		plate_lenghts = np.linalg.norm(shifted_points-self.airfoil_points, axis=0)
		self.area = np.sum(self.wall_t*plate_lenghts)
		self.circumference = np.sum(plate_lenghts)
		return self.area


	def get_I(self):
		shifted_points = np.roll(self.airfoil_points, 1, axis=1)
		center_points = (self.airfoil_points + shifted_points)/2
		plate_angles = np.arctan2(shifted_points[1]-self.airfoil_points[1], shifted_points[0]-self.airfoil_points[0])
		plate_lenghts = np.linalg.norm(shifted_points-self.airfoil_points, axis=0)
		plate_areas = self.wall_t*plate_lenghts

		self.centroid = np.sum(center_points*plate_areas, axis=1)/np.sum(plate_areas)

		self.I_xx = np.sum((plate_lenghts**3*self.wall_t*np.sin(plate_angles)**2)/12 + plate_areas*(center_points[1]-self.centroid[1])**2)
		self.I_yy = np.sum((plate_lenghts**3*self.wall_t*np.cos(plate_angles)**2)/12 + plate_areas*(center_points[0]-self.centroid[0])**2)
		self.I_xy = np.sum((plate_lenghts**3*self.wall_t*np.sin(2*plate_angles))/24  + plate_areas*(center_points[0]-self.centroid[0])*(center_points[1]-self.centroid[1]))

		#print(center_points)
		#plt.plot(self.airfoil_points[0], self.airfoil_points[1], linewidth=0, marker="o")
		#plt.plot(center_points[0], center_points[1],linewidth=0, marker="o")
		#plt.plot(self.centroid[0], self.centroid[1], linewidth=0, marker="o")
		#plt.show()
		return self.I_xx, self.I_yy, self.I_xy
	def get_J(self):
		self.get_area()
		self.J = (4*self.wall_t*self.enclosed_area**2)/self.circumference
		return self.J

	def plot_stresses(self, moments):
		#still needs to be implemented
		stresses = (moments[0]*self.airfoil_points[1])/self.I_xx + (moments[1]*self.airfoil_points[0])/self.I_yy
		norm = mpcolor.Normalize(vmin=-np.max(stress_curve), vmax=np.max(stress_curve))
		cmap = mpclm.coolwarm.reversed()
		m = mpclm.ScalarMappable(norm=norm, cmap=cmap)
		
		fig3, ax3 = plt.subplots()

		l = len(points)
		for idx in range(l):
			ax3.plot(self.airfoil_points[idx:(idx+2)], self.airfoil_points[idx:(idx+2)], color = m.to_rgba(stress_curve[idx]))
		#ax3.set_aspect(1)
		plt.show()


	def get_max_stress(self, moments):
		stresses = (moments[0]*self.airfoil_points[1]*self.I_yy)/(self.I_xx*self.I_yy-self.I_xy**2) + (moments[1]*self.airfoil_points[0]*self.I_xx)/(self.I_xx*self.I_yy-self.I_xy**2) 
		return np.max(np.abs(stresses))

	def get_torsion_stress(self, moments):
		q = moments[2]/(2*self.enclosed_area)
		return q/self.wall_t

class ellipse():
	"""ellipse section class to determine properties, stresses"""
	def __init__(self, geom, offset, mat_prop): 
		#geom=(a_inner, b_inner, a_outer, b_outer), #mat_prop =(E, sigma_yield)
		self.ai, self.bi, self.ao, self.bo = geom
		self.offset = np.array(offset)
		self.E, self.G, self.sigma_y = mat_prop['E'], mat_prop['G'], mat_prop['sigma_y']

		self.get_I()
		self.get_J()

	def resize(self, chord):
		tc = self.bo/self.ao
		t = self.ao-self.ai

		self.ao = chord/2
		self.bo = self.ao*tc
		self.set_inner_t(t)

	def set_inner_t(self, t):
		self.ai = self.ao - t
		self.bi = self.bo - t
		self.get_I()
		self.get_J()
		return self

	def get_area(self):
		self.area = np.pi*(self.ao*self.bo-self.ai*self.bi)
		return self.area
	def get_circumference(self):
		self.circumference = 2*np.pi*(1.32*0.5*(self.ao+self.bo) -0.32*np.sqrt(self.ao*self.bo))

	def get_I(self):
		self.get_area()
		self.I_xx = 0.25*np.pi*(self.ao*self.bo**3 - self.ai*self.bi**3) + self.area*self.offset[1]**2
		self.I_yy = 0.25*np.pi*(self.bo*self.ao**3 - self.bi*self.ai**3) + self.area*self.offset[0]**2
		self.I_xy = 0
		return self.I_xx, self.I_yy, self.I_xy

	def get_J(self):
		self.get_area()
		self.J = 0.25*np.pi*(self.ao*self.bo*(self.ao**2+self.bo**2) - self.ai*self.bi*(self.ai**2+self.bi**2))
		return self.J

	def get_max_stress(self, moments, forces=(0,0,0)):
		#forces: Fx, Fy, Fz
		#moments: Mx, My, Tz (SA coordinate system) 
		#this would need to be changed if composite beams are implemented
		stress = lambda x, y: (moments[0]*y)/self.I_xx + (moments[1]*x)/self.I_yy
		y_ellipse = lambda x: self.bo*np.sqrt(1-(x**2/self.ao**2))

		q1 = spo.minimize(lambda x: -np.abs(stress(x,y_ellipse(x))), 0, bounds=spo.Bounds(0, self.ao))
		q2 = spo.minimize(lambda x: -np.abs(stress(x,y_ellipse(x))), 0, bounds=spo.Bounds(-self.ao, 0))

		self.max_stress = max(-q1.fun, -q2.fun)
		#print(-q1.fun/1e6, -q2.fun/1e6)
		#print(q1.x/self.ao, y_ellipse(q1.x)/self.bo)
		return self.max_stress

	def get_enclosed_area(self):
		self.enclosed_area = np.pi*self.ao*self.bo
		return self.enclosed_area
	def get_torsion_stress(self, moments):
		self.torsion_stress = moments[2]/(2*self.get_enclosed_area()*(self.ao-self.ai))
		return self.torsion_stress

	def get_torsion_twist(self, loads):
		self.get_circumference()
		#print(self.circumference, "C", self.area,"A")
		self.twist_per_length = (loads['M_z']*self.circumference)/(4*self.get_enclosed_area()**2*self.G*(self.ao-self.ai))
		return self.twist_per_length


	#TODO: torsional properties, stress calculation


class circle(ellipse):
	'''circle is an ellipse with a=b=r'''
	def __init__(self, geom, offset, mat_prop): 
		#geom = (ri, ro), #mat_prop =(E, sigma_yield)
		eq_ellipse_geom = (geom[0], geom[0], geom[1], geom[1])
		super().__init__(self, eq_ellipse_geom, offset, mat_prop)

class rectangle():
	def __init__(self, geom, offset, mat_prop):
		self.bi, self.hi, self.bo, self.ho = geom # inner width, inner height, outer ...
		self.offset = np.array(offset)
		self.E, self.G, self.sigma_y = mat_prop['E'], mat_prop['G'], mat_prop['sigma_y']
		
		self.area()
		self.get_I()

	def set_inner_t(self, t):
		self.bi = self.bo - t
		self.hi = self.ho -t 


	def area(self):
		self.area =  self.bo*self.ho - self.bi*self.hi
		return self.area

	def get_I(self):
		self.I_xx = (1/12)*(self.bo*self.ho**3 - self.bi*self.hi**3) + self.area*self.offset[1]**2
		self.I_yy = (1/12)*(self.ho*self.bo**3 - self.hi*self.bi**3) + self.area*self.offset[0]**2
		return self.I_xx, self.I_yy
	
	#TODO: torsional properties, stress calculation
		




		