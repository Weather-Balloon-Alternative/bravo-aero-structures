import numpy as np
import matplotlib.pyplot as plt

# compute bending moment
def bending_moment(n,W,S_h,S,l_tail):
    M_x = n * W * S_h/S * l_tail
    return M_x

# compute  shear
def shear(n,W,b,S_h,S):
    V_y = (n-1-n*S_h/S)*W
    L_centre = b/6
    L_centre_tail = b_h/6
    T_z = n*W*L_centre + n*W*L_centre_tail*S_h/S
    return V_y, T_z

# Compute and optimize for stresses in a tube cross-section
def tube_stresses(M_x, V_y, T_z):
    # Initialize radius and thickness options
    r_range = np.linspace(r_min, r_max, int(r_max/r_step + 1))
    design_options = [0., 0., 0., 0., 0.,0.,0.]

    for id in range(len(r_range)):
        r = r_range[id]

        r_array = np.ones((id+1))*r
        t_range = np.linspace(t_min, r, (id+1))

        # Section properties
        I_xx = (np.pi * r_array**4) / 4 - (np.pi * (r_array - t_range)**4) / 4
        J = (np.pi * r_array**4) / 2 - (np.pi * (r_array - t_range)**4) / 2
        A = np.pi * (r_array**2 - (r_array-t_range)**2)

        # stresses for solid tube
        bending_stress = M_x * r_array / I_xx
        torsional_shear = T_z * r_array / J
        vertical_shear = -V_y / I_xx * 2 * np.pi * r_array ** 3 * np.sin(1 / (4 * r_array))
        shear_stress = torsional_shear+vertical_shear

        # deformations
        twist = T_z / (G*J) * l_tail
        vert_def = M_x * l_tail**2 / (3*E*I_xx)

        for idx in range(len(t_range)):
            section_properties = [r, t_range[idx], A[idx], twist[idx], vert_def[idx],bending_stress[idx],shear_stress[idx]]
            design_options = np.vstack((design_options,section_properties))



    return design_options[1:,:]


if __name__ == "__main__":
    # load case
    n = 4.5
    W = 0.92504571*9.81
    # maximum deformatations and stresses
    safety_factor  = 1.4
    max_twist = 2 * np.pi/180 / safety_factor
    max_vert_def = 0.01 / safety_factor
    CFRP_yield_stress = 110 * 10**9 / safety_factor
    CFRP_shear_strength = 260 * 10**9 / safety_factor
    # aircraft geometry
    S = 0.09
    S_h = 0.011833342
    b = 1.161895004/2
    b_h =  0.21756233/2
    l_tail = 0.407380764
    # numerical parameters
    r_min = 0.001
    r_max = 0.01
    r_step = 0.001
    t_min = 0.0005
    t_max = 0.01
    t_step = 0.0005
    # material properties
    G = 30 * 10**9
    E = 17 * 10**9
    CFRP_density = 2000 # kg/m^3

    #output
    # Internal loads
    M_x = bending_moment(n,W,S_h, S,l_tail)
    V_y, T_z = shear(n, W, b, S_h, S)

    # deformations
    r,t,A,twist,vert_def,bending_stress,shear_stress = 0,1,2,3,4,5,6
    design_options = tube_stresses(M_x,V_y,T_z)

    A_array = design_options[:,A]
    twist_array = design_options[:,twist]
    vert_def_array = design_options[:,vert_def]
    bending_stress_array = design_options[:,bending_stress]
    shear_stress_array = design_options[:,shear_stress]

    print(design_options[38])
    mass_per_metre = design_options[38,A]*CFRP_density
    print(mass_per_metre)
    print('')
    print(design_options)
    # plt.plot(A_array,twist_array,color='green',label='twist')
    # plt.plot(A_array,vert_def_array,color='red',label='vertical deformation')
    # plt.show()
    #
    # plt.plot(A_array,bending_stress_array, color='green', label='bending_stress')
    # plt.plot(A_array, shear_stress_array, color='red', label='torsional_shear')
    # plt.show()