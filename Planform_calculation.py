import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import isacalc
import scipy.stats
import scipy.integrate
import scipy.optimize
import blimb_calc
from ambiance import Atmosphere


#input:
# glide ratio
#mass estimations


#lift over drag
	CL = (np.pi * e * AR * CD0) ** 0.5
	CD = CD0 + (CL ** 2) / (np.pi * AR * e)
	glide_ratio = CL / CD







if __name__ == '__main__':
