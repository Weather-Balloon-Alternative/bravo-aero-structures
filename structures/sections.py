import numpy as np

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







class ellipse():
	"""ellipse section class to determine properties, stresses"""
	def __init__(self, geom, offset, mat_prop): 
		#geom=(a_inner, b_inner, a_outer, b_outer), #mat_prop =(E, sigma_yield)
		self.ai, self.bi, self.ao, self.bo = geom
		self.offset = np.array(offset)
		self.E = mat_prop[0]
		self.sigma_y = mat_prop[1]

		self.get_area()
		self.get_I()
		self.get_J()

	def set_inner_t(self, t):
		self.ai = self.ao - t
		self.bi = self.bo - t
		self.get_area()
		self.get_I()
		self.get_J()
		return self

	def get_area(self):
		self.area = np.pi*(self.ao*self.bo-self.ai*self.bi)
		return self.area

	def get_I(self):
		self.I_xx = 0.25*np.pi*(self.ao*self.bo**3 - self.ai*self.bi**3) + self.area*self.offset[1]**2
		self.I_yy = 0.25*np.pi*(self.bo*self.ao**3 - self.bi*self.ai**3) + self.area*self.offset[0]**2
		return self.I_xx, self.I_yy

	def get_J(self):
		self.J = 0.25*np.pi*(self.ao*self.bo*(self.ao**2+self.bo**2) - self.ai*self.bi*(self.ai**2+self.bi**2))
		return self.J

	#TODO: torsional properties, stress calculation


class circle(ellipse):
	'''circle is an ellipse with a=b=r'''
	def __init__(self, geom, offset, mat_prop): 
		#geom = (ri, ro), #mat_prop =(E, sigma_yield)
		eq_ellipse_geom = (geom[0], geom[0], geom[1], geom[1])
		super().__init__(self, eq_ellipse_geom, offset, mat_prop)

class rectangle():
	def __init(self, geom, offset, mat_prop):
		self.bi, self.hi, self.bo, self.ho = geom
		self.offset = np.array(offset)
		self.E = mat_prop[0]
		self.sigma_y = mat_prop[1]

		self.area()
		self.get_I()

	def area(self):
		self.area =  bo*ho - bi*hi
		return self.area

	def get_I(self):
		self.I_xx = (1/12)*(self.bo*self.ho**3 - self.bi*self.hi**3) + A*offset[1]**2
		self.I_yy = (1/12)*(self.ho*self.bo**3 - self.hi*self.bi**3) + A*offset[0]**2
		return self.I_xx, self.I_yy
	
	#TODO: torsional properties, stress calculation
		




		