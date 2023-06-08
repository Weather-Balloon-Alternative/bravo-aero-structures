from coefficients import *
from initial_empennage import *
from scissorplot import *
from eigenvalues import *
from controlsurfacerequirements import *

def main_sizingcheck(radians=True,print_intermediate=True):
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
    scissorplt = scissorplot(CLh,CLAh,lh,c,VhV,Cmac,xacbar,CLah,CLaAh,deda,SM,xcg_fwBAR,xcg_aftBAR,PLOT=False)
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

    if print_intermediate:
        print(f"Asymmetric eigenmotion results: (Use for vertical tail)\n\
    lambda_aperiodic:                   {lambda_aperiodic}  Passed: {np.real(lambda_aperiodic) < 0}\n\
    lambda_dutchroll:                   {lambda_dutchroll}  Passed: {np.real(lambda_dutchroll) < 0}\n\
    lambda_spiral:                      {lambda_spiral}     Passed: {np.real(lambda_spiral) < 0}\n")
        print(f'Overview of the asymmetric state-space:')
        wnasym, zetaasym, polesasym = control.damp(ss_asymmetricEOM(CYbdot,mub,b,V0,KX2,KXZ,KZ2,CYb,CL,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr))
        print(f"\nMinimum required damping ratios:\n\
    Aperiodic roll:                     {dampratio_aperiodic_min}\n\
    Dutch roll:                         {dampratio_dutchroll_min}\n\
    Spiral:                             {dampratio_spiral_min}\n")
        print(f'Pre-checking horizontal tail parameters...\n')
    else:
        wnasym, zetaasym, polesasym = control.damp(ss_asymmetricEOM(CYbdot,mub,b,V0,KX2,KXZ,KZ2,CYb,CL,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr),doprint=False)

    # Pre-check symmetrical stability coefficients
    if print_intermediate:
        print(f'C_m_alpha right sign?   {Cma < 0}\n')
        print(f'Calculating eigenvalues of symmetric eigenmotions...\n')
    
    # Calculate eigenvalues of symmetric eigenmotions
    lambda_phugoid = eigenvalues_phugoid(muc,CZa,Cmq,Cma,CXu,Cmu,CXa,CZu,CZ0,V0,c)
    lambda_shortperiod = eigenvalues_shortperiod(muc,KY2,CZadot,CZq,Cmadot,Cmq,CZa,Cma,V0,c)
    if print_intermediate:
        print(f"Symmetric eigenmotion results: (Use for horizontal tail)\n\
    lambda_phugoid:                     {lambda_phugoid}    Passed: {np.real(lambda_phugoid) < 0}\n\
    lambda_shortperiod:                 {lambda_shortperiod}    Passed: {np.real(lambda_shortperiod) < 0}\n")
        print('Overview of the symmetrical state-space: (Match them up yourself)')
        wnsym, zetasym, polessym = control.damp(ss_symEOM(muc,c,V0,CZadot,Cmadot,KY2,CXu,CXa,CZ0,CXq,CZu,CZa,CX0,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde))
        print(f'\nMinimum required damping ratios:\n\
    Phugoid:                            {dampratio_phugoid_min}\n\
    Short period:                       {dampratio_shortperiod_min}\n')
        print('Calculating pitch and bank performance...')
    else:
        wnsym, zetasym, polessym = control.damp(ss_symEOM(muc,c,V0,CZadot,Cmadot,KY2,CXu,CXa,CZ0,CXq,CZu,CZa,CX0,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde),doprint=False)

    # Calculate pitch and bank performance


if __name__ == '__main__':
    main_sizingcheck()