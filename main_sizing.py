from initial_empennage import *
from scissorplot import *
from eigenvalues import *
from controlsurfacerequirements import *

def main_sizingcheck(radians=True,print_intermediate=True,Maneuver_SF=2):
    """_summary_

    Args:
        radians (bool, optional): _description_. Defaults to True.
        print_intermediate (bool, optional): _description_. Defaults to True.
        Maneuver_SF (int, optional): _description_. Defaults to 2.

    Returns:
        _type_: _description_
    """    
    # Size initial tail
    tailpar = tailsizing(b,c,ShS=ShS,Arh=ARh,lhc=lhc,Arv = 1.1)
    if print_intermediate:
        print(f"Initial sizing results:\n\
    l_h:                                {tailpar['l_h']}\n\
    S_h:                                {tailpar['S_h']}\n\
    b_h:                                {tailpar['b_h']}\n\
    c_h:                                {tailpar['c_h']}\n\
    S_v:                                {tailpar['S_v']}\n\
    b_v:                                {tailpar['b_v']}\n\
    c_v:                                {tailpar['c_v']}\n\
            ")
        print(f'Performing scissor plot analysis...\n')
    
    # Perform scissor plot analysis
    scissorplt = scissorplot(CLh,CLAh,lh,c,VhV,Cmac,xacbar,CLah,CLaAh,deda,SM,xcg_fwBAR,xcg_aftBAR,PLOT=False,IgnoreErrors=True)
    ShS1 = scissorplt['ShS']
    if print_intermediate:
        if scissorplt['xcg_waste'] > 0:
            print(f"Scissor plot results:\n\
    Sh/S:                               {scissorplt['ShS']}\n\
    Unused x_cg excursion:              {scissorplt['xcg_waste']}\n\
    Recommendation:                     Move main wing forwards\n\
                  ")
        else:
            print(f"Scissor plot results:\n\
    Sh/S:                               {scissorplt['ShS']}\n\
    Unused x_cg excursion:              {scissorplt['xcg_waste']}\n\
    Recommendation:                     Move main wing backwards\n\
                  ")
        print(f'Pre-checking vertical tail parameters...\n')
            
    # Pre-check asymmetrical stability coefficients
    if print_intermediate:
        print(f'C_n_beta right sign?    {Cnb > 0}\n')
        print(f'Calculating eigenvalues of asymmetric eigenmotions...\n')

    # Calculate eigenvalues of asymmetric eigenmotions
    lambda_aperiodic = eigenvalues_aperiodicroll(Clp,mub,KX2,V0,b)
    lambda_dutchroll = eigenvalues_dutchroll(mub,KZ2,Cnr,Cnb,CYb,V0,b)
    lambda_spiral = eigenvalues_spiral(CL,Clb,Cnr,Cnb,Clr,Clp,CYb,mub,Cnp,V0,b)
    ss_asym = ss_asymmetricEOM(CYbdot,mub,b,V0,KX2,KXZ,KZ2,CYb,CL,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr)
    if print_intermediate:
        print(f"Asymmetric eigenmotion results: (Use for vertical tail)\n\
    lambda_aperiodic:                   {lambda_aperiodic}  Passed: {np.real(lambda_aperiodic) < 0}\n\
    lambda_dutchroll:                   {lambda_dutchroll}  Passed: {np.real(lambda_dutchroll) < 0}\n\
    lambda_spiral:                      {lambda_spiral}     Passed: {np.real(lambda_spiral) < 0}\n")
        print(f'Overview of the asymmetric state-space: (Match them up yourself)')
        wnasym, zetaasym, polesasym = control.damp(ss_asym)
        print(f"\nMinimum required damping ratios:\n\
    Aperiodic roll:                     {dampratio_aperiodic_min}\n\
    Dutch roll:                         {dampratio_dutchroll_min}\n\
    Spiral:                             {dampratio_spiral_min}\n")
        print(f'Pre-checking horizontal tail parameters...\n')
    else:
        wnasym, zetaasym, polesasym = control.damp(ss_asym,doprint=False)

    # Pre-check symmetrical stability coefficients
    if print_intermediate:
        print(f'C_m_alpha right sign?   {Cma < 0}\n')
        print(f'Calculating eigenvalues of symmetric eigenmotions...\n')
    
    # Calculate eigenvalues of symmetric eigenmotions
    lambda_phugoid = eigenvalues_phugoid(muc,CZa,Cmq,Cma,CXu,Cmu,CXa,CZu,CZ0,V0,c)
    lambda_shortperiod = eigenvalues_shortperiod(muc,KY2,CZadot,CZq,Cmadot,Cmq,CZa,Cma,V0,c)
    ss_sym = ss_symEOM(muc,c,V0,CZadot,Cmadot,KY2,CXu,CXa,CZ0,CXq,CZu,CZa,CX0,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde)

    if print_intermediate:
        print(f"Symmetric eigenmotion results: (Use for horizontal tail)\n\
    lambda_phugoid:                     {lambda_phugoid}    Passed: {np.real(lambda_phugoid) < 0}\n\
    lambda_shortperiod:                 {lambda_shortperiod}    Passed: {np.real(lambda_shortperiod) < 0}\n")
        print('Overview of the symmetrical state-space: (Match them up yourself)')
        wnsym, zetasym, polessym = control.damp(ss_sym)
        print(f'\nMinimum required damping ratios:\n\
    Phugoid:                            {dampratio_phugoid_min}\n\
    Short period:                       {dampratio_shortperiod_min}\n')
        print('Calculating pitch and bank performance...')
    else:
        wnsym, zetasym, polessym = control.damp(ss_sym,doprint=False)

    # Calculate pitch and bank performance
    # Pitch
    
    pitchratereq = pitchraterequirements(V0,h,nmax)
    tmeasure = 1
    tpitch,xpitch = control.forced_response(ss_sym,T=np.linspace(0,tmeasure,1000),U=np.ones(1000)*delta_e_max)
    pitchreqmet = xpitch[2][-1]/tmeasure > pitchratereq['ang_acc'][0]*Maneuver_SF

    # Bank
    bankmaneuver = np.vstack((np.ones(1000)*delta_a_max,np.zeros(1000)))
    tbank,xbank = control.forced_response(ss_asym,T=np.linspace(0,t_to_60deg_bank,1000),U=bankmaneuver)
    bankreqmet = xbank[1][-1] > 60*np.pi*Maneuver_SF/180

    if print_intermediate:
        print(f'Pitch and bank performance results:\n\
    Sufficient pitch performance: {pitchreqmet}\n\
    Sufficient bank performance: {bankreqmet}\n')
    
    return pitchreqmet,bankreqmet

if __name__ == '__main__':
    # from coefficients_smallglider import *
    from coefficients_bigglider import *
    # from citpar import *
    # main_sizingcheck()

    Obtain_ControlSurfaceSizing = True
    # Obtain_ControlSurfaceSizing = False
    if Obtain_ControlSurfaceSizing:
        resolution = 10E-4
        CZde_opt = np.arange(-1,0,resolution)
        CZde_best = np.nan

        Clda_opt = np.arange(-1,0,resolution)
        Clda_best = np.nan

        for opt in range(len(CZde_opt)):
            CZde = CZde_opt[opt]
            Cmde = CZde*lhc #elevator efficiency, derivative of Cm w.r.t. delta_e
            Clda = Clda_opt[opt]
            pitchreqmet,bankreqmet = main_sizingcheck(print_intermediate=False)
            if pitchreqmet:
                CZde_best = CZde

            if bankreqmet:
                Clda_best = Clda

        print(f'Best CZde = {CZde_best} 1/rad = {CZde_best*np.pi/180} 1/deg')
        print(f'Best Clda = {Clda_best} 1/rad = {Clda_best*np.pi/180} 1/deg')

    # ss_sym = ss_symEOM(muc,c,V0,CZadot,Cmadot,KY2,CXu,CXa,CZ0,CXq,CZu,CZa,CX0,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde)
    # tmeasure = 1
    # tpitch,xpitch = control.forced_response(ss_sym,T=np.linspace(0,tmeasure,1000),U=np.ones(1000)*delta_e_max)
    # plt.plot(tpitch,xpitch[2])
    # plt.show()


    # ss_asym = ss_asymmetricEOM(CYbdot,mub,b,V0,KX2,KXZ,KZ2,CYb,CL,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr)
    # bankmaneuver = np.vstack((np.ones(1000)*delta_a_max,np.zeros(1000)))
    # tbank,xbank = control.forced_response(ss_asym,T=np.linspace(0,t_to_60deg_bank,1000),U=bankmaneuver)
    # plt.plot(tbank,xbank[1])
    # plt.show()