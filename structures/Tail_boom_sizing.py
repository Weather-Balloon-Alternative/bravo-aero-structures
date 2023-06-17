import numpy as np
import matplotlib.pyplot as plt
import sections
import beam

# compute bending moment
def bending_moment(n,W,S_h,S,l_tail):
    M_x = (n-1) * W * S_h/S * l_tail
    return M_x

# compute  shear
def shear(n,W,b,S_h,S):
    V_y = (n*S_h/S)*W
    L_centre = b/6
    L_centre_tail = b_h/6
    T_z = (n-1)*W*L_centre + (n-1)*W*L_centre_tail*S_h/S
    return V_y, T_z

# Compute and optimize for stresses in a tube cross-section
def tube_stresses(M_x, V_y, T_z):
    # Initialize radius and thickness options
    r_range = np.linspace(r_min, r_max, int(r_max/r_step - 1))
    design_options = [0., 0., 0., 0., 0.,0.,0.,0.]

    for id in range(len(r_range)):
        r = r_range[id]
        r_array = np.ones((id+2))*r
        t_range = np.linspace(t_min, r, (id+2))

        # Section properties
        I_xx = (np.pi * r_array**4) / 4 - (np.pi * (r_array - t_range)**4) / 4
        J = (np.pi * r_array**4) / 2 - (np.pi * (r_array - t_range)**4) / 2
        A = np.pi * (r_array**2 - (r_array-t_range)**2)
        mass_per_metre = A * CFRP_density

        # stresses for solid tube
        bending_stress = M_x * r_array / I_xx
        torsional_shear = T_z * r_array / J
        vertical_shear = -V_y / I_xx * 2 * np.pi * r_array ** 3 * np.sin(1 / (4 * r_array))
        shear_stress = torsional_shear+vertical_shear

        # deformations
        twist = T_z / (G*J) * l_tail
        vert_def = M_x * l_tail**2 / (3*E*I_xx)

        for idx in range(len(t_range)):
            section_properties = [r, t_range[idx], A[idx], twist[idx], vert_def[idx],bending_stress[idx],shear_stress[idx],mass_per_metre[idx]]
            design_options = np.vstack((design_options,section_properties))



    return design_options[1:,:]


def big_tailboom_sizing(W, n, l_tail, S_h, S, b):

    h_r = 0.14
    h_t = 0.028
    h = lambda z: h_r - z*(h_r-h_t)/l_tail
    ar_tailboom = 1
    t = 0.0027/1000
    mat_prop_skin = {"E":24.1e9, "G":10e9, "sigma_y":500e6}


    zz = np.linspace(0, l_tail, 100)
    II = np.zeros(zz.shape)
    tailboom_section = sections.rectangle((h(0)*ar_tailboom - t, h(0) -t, h(0)*ar_tailboom, h(0)), (0,0), mat_prop_skin)
    for idx, z in enumerate(zz):
        tailboom_section.bo = h(z)*ar_tailboom
        tailboom_section.ho = h(z)
        tailboom_section.set_inner_t(t)

        II[idx] = tailboom_section.get_I()[0]
    print(II[0])
    V_y = np.zeros(zz.shape)
    V_y[-1] = -shear(n, W, b, S_h, S)[0]
    defl = beam.deflection(zz, 0, V_y, 0, tailboom_section.E, II, plot=True)
    print(defl[-1][-1])

    stress_curve = (defl[1]*0.5*h(zz))/II
    beam.plot_stress_curve(zz, defl[-1], stress_curve)
    print(np.max(stress_curve)/1e6)




if __name__ == "__main__":
    # load case
    n = 2.5
    W = 0.773736504*9.81
    # maximum deformations and stresses
    safety_factor  = 1.5
    max_twist = 2 * np.pi/180 / safety_factor
    max_vert_def = 0.01 / safety_factor
    CFRP_yield_stress = 110 * 10**6 / safety_factor
    CFRP_shear_strength = 260 * 10**6 / safety_factor
    # aircraft geometry
    S = 0.05
    S_h = 0.005960478
    b = 0.774596669/2
    b_h =  0.133721484/2
    l_tail = 0.368958279
    # numerical parameters
    r_min = 0.001
    r_max = 0.01
    r_step = 0.0005
    t_min = 0.0005
    t_max = 0.01
    t_step = 0.0005
    # material properties
    G = 33 * 10**9
    E = 150 * 10**9
    CFRP_density = 2000 # kg/m^3

    #output
    # Internal loads
    M_x = bending_moment(n,W,S_h, S,l_tail)
    V_y, T_z = shear(n, W, b, S_h, S)

    # deformations
    r,t,A,twist,vert_def,bending_stress,shear_stress, = 0,1,2,3,4,5,6
    design_options = tube_stresses(M_x,V_y,T_z)

    A_array = design_options[:,A]
    twist_array = design_options[:,twist]
    vert_def_array = design_options[:,vert_def]
    bending_stress_array = design_options[:,bending_stress]
    shear_stress_array = design_options[:,shear_stress]

    viable_options = np.where((twist_array <= max_twist) & (vert_def_array <= max_vert_def) & (bending_stress_array <= CFRP_yield_stress) & (shear_stress_array <= CFRP_shear_strength))
    MVoption_id = np.min(viable_options)
    viableA = A_array[MVoption_id:]
    optimal_option_id = np.argmin(viableA)
    optimal_option = design_options[(MVoption_id+optimal_option_id+1)]
    mass_per_metre = optimal_option[A] * CFRP_density

<<<<<<< HEAD
    plt.plot(range(174),A_array[(231-174):])
    plt.plot(range(174),design_options[(231-174):,0])
    plt.plot(range(174),design_options[(231-174):,1])
   # plt.show()

    W = 3.9*9.81
    l_tail = 0.43
    b = 1.60623784
    S_h = 0.0215
    S = 0.215
    big_tailboom_sizing(W, n, l_tail, S_h, S, b)
=======
    # print output
    print('Outer radius = ' + str(optimal_option[0]) + ' m \n' + 'Thickness = '+ str(optimal_option[1]) + ' m \n' + 'Cross-sectional area = ' + str(optimal_option[2]) + ' mm^2 \n' + 'twist angle = ' + str(optimal_option[3]) + ' radians \n' + 'Vertical deformation = ' + str(optimal_option[4]) + ' m \n' + 'Bending stress = ' + str(optimal_option[5]) + ' Pa \n' + 'Shear stress = ' + str(optimal_option[6]) + ' Pa \n' + 'Mass per metre = ' + str(optimal_option[7]) + ' kg/m \n')

>>>>>>> 206be4c9aa8fbcddf74cfc88600ec1554fa518de
