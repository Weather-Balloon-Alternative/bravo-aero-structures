import numpy as np
import matplotlib.pyplot as plt

# Estimate bending moment
def bending_moment(n,W,S_h,S,l_tail):
    M_x = n * W * S_h/S * l_tail
    return M_x

# Estimate  shear
def shear(n,W,b,S_h,S):
    V_y = (n-1-n*S_h/S)*W
    L_centre = b/6
    T_z = n*W*L_centre+(1+S_h/S)
    return V_y, T_z

# Compute and optimize for stresses in a tube cross-section
def tube_stresses(M_x, V_y, T_z):
    # Initialize radius and thickness options
    r_range = np.arange(r_min, r_max, r_step)
    t_range = np.arange(t_min, t_max, t_step)

    # # Section properties
    # I_xx = (np.pi * r_range)/4 - (np.pi * (r_range-t_range))/4
    # J = (np.pi * r_range)/2 - (np.pi * (r_range-t_range))/2
    # # stresses for solid tube
    # bending_stress = M_x*r_range/I_xx
    # torsional_shear = T_z*r_range/J
    # vertical_shear = -V_y/I_xx * 2*np.pi*r_range**3 * np.sin(1/(4*r_range))
    #
    for r in r_range:
        id = np.argwhere(r_range,r)
        t_range = t_range[r]
        # Section properties
        I_xx = (np.pi * r) / 4 - (np.pi * (r - t_range)) / 4
        J = (np.pi * r_range) / 2 - (np.pi * (r_range - t_range)) / 2
        # stresses for solid tube
        bending_stress = M_x * r_range / I_xx
        torsional_shear = T_z * r_range / J
        vertical_shear = -V_y / I_xx * 2 * np.pi * r_range ** 3 * np.sin(1 / (4 * r_range))
    return bending_stress, torsional_shear, vertical_shear


if __name__ == "__main__":
    n = 1
    W = 2
    S = 3
    S_h = 4
    b = 0.35
    l_tail = 5
    r_min = 0.000
    r_max = 0.02
    r_step = 0.0001
    t_min = 0.000
    t_max = 0.02
    t_step = 0.0001

    # Internal loads
    M_x = bending_moment(n,W,S_h, S,l_tail)
    V_y, T_z = shear(n, W, b, S_h, S)

    # Compute internal stresses
    bending_stress, torsional_shear, vertical_shear = tube_stresses(M_x,V_y,T_z)
    print(bending_stress)
    print(vertical_shear)
    print(torsional_shear)
