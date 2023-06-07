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
        I_xx = (np.pi * r_array) / 4 - (np.pi * (r_array - t_range)) / 4
        J = (np.pi * r_array) / 2 - (np.pi * (r_array - t_range)) / 2
        # stresses for solid tube
        bending_stress = M_x * r_array / I_xx
        torsional_shear = T_z * r_array / J
        vertical_shear = -V_y / I_xx * 2 * np.pi * r_array ** 3 * np.sin(1 / (4 * r_array))

        # deformations
        twist = T_z / (G*J) * l_tail
        vert_def = M_x * l_tail**2 / (3*E*I_xx)

        for idx in range(len(t_range)):
            A = np.pi*(r**2-(r-t_range[idx])**2)
            section_properties = [r, t_range[idx], A, twist[idx], vert_def[idx],bending_stress[idx],torsional_shear[idx]]
            design_options = np.vstack((design_options,section_properties))



    return design_options[1:,:]


if __name__ == "__main__":
    # load case
    n = 2.5
    W = 0.8
    # maximum deformations
    max_twist = 5 * np.pi/180
    max_vert_def = 0.05
    # aircraft geometry
    S = 0.0787
    S_h = 0.0135
    b = 1.086508168/2
    b_h =  0.232379001
    l_tail = 0.351696452
    # numerical parameters
    r_min = 0.001
    r_max = 0.01
    r_step = 0.001
    t_min = 0.001
    t_max = 0.01
    t_step = 0.001
    # material properties
    G = 30 * 10**9
    E = 17 * 10**9
    CFRP_density = 2000 # kg/m^3

    #output
    # Internal loads
    M_x = bending_moment(n,W,S_h, S,l_tail)
    V_y, T_z = shear(n, W, b, S_h, S)

    # deformations
    r,t,A,twist,vert_def,bending_stress,torsional_shear = 0,1,2,3,4,5,6
    design_options = tube_stresses(M_x,V_y,T_z)
    print(design_options.shape)
    A_array = design_options[:,A]
    twist_array = design_options[:,twist]
    vert_def_array = design_options[:,vert_def]
    bending_stress_array = design_options[:,bending_stress]
    torsional_shear_array = design_options[:,torsional_shear]

    # plt.plot(A_array,twist_array,color='green',label='twist')
    # plt.plot(A_array,vert_def_array,color='red',label='vertical deformation')
    # plt.show()
    #
    # plt.plot(A_array,bending_stress_array,color='green',label='bending_stress')
    # plt.plot(A_array, torsional_shear_array, color='red', label='torsional_shear')
    # plt.show()
    #
    # A = np.pi*(0.01**2-0.009**2)
    # mass_per_metre = A*CFRP_density
    # print(A, mass_per_metre)

    print(bending_moment(n,W,S_h,S,l_tail))
    print(shear(n,W,S_h,S))


