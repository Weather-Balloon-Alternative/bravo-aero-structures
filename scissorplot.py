import numpy as np
import matplotlib.pyplot as plt
from coefficients import *

def controllability(CLh,CLAh,lh,c,VhV,xcg,Cmac,xac):
    '''
    Controllability equation for the scissor plot,
    returns the ratio of tail area over wing area

    args:

    CLh:    float, C_L_h;   lift coefficient of the tail
    CLAh:   float, C_L_A-h; lift coefficient of the aircraft without the tail 
    lh:     float, l_h;     distance from the aerodynamic centre to the tail
    c:      float, c;       chord
    VhV:    float, V_h/V;   ratio of the velocity at the tail over the aircraft velocity
    xcg:    float, x_cgBAR; location of the centre of gravity divided by the chord
    Cmac:   float, C_m_ac;  moment coefficient of the aerodynamic centre
    xac:    float, x_acBAR; location of the aerodynamic centre divided by the chord
    
    returns:

    ShS:    float, S_h/S;   ratio of tail area over wing area
    
    '''
    ShS = (xcg + Cmac - xac) / (CLh/CLAh * lh/c * VhV**2)
    return ShS

def stability(Clah,ClaAh,deda,lh,c,VhV,xcg,xac,SM):
    '''
    Stability equation for the scissor plot,
    returns the ratio of tail area over wing area

    args:

    CLah:   float, C_L_alpha_h;     derivative of the lift coefficient of the tail w.r.t. alpha
    CLaAh:  float, C_L_alpha_A-h;   derivative of lift coefficient of the aircraft without the tail w.r.t alpha
    deda:   float, deta/dalpha;     derivative of the downwash angle w.r.t. alpha
    lh:     float, l_h;             distance from the aerodynamic centre to the tail
    c:      float, c;               chord
    VhV:    float, V_h/V;           ratio of the velocity at the tail over the aircraft velocity
    xcg:    float, x_cgBAR;         location of the centre of gravity divided by the chord
    xac:    float, x_acBAR;         location of the aerodynamic centre divided by the chord
    SM:     float, S.M.;            stability margin
    
    returns:

    ShS:    float, S_h/S;   ratio of tail area over wing area
    '''
    ShS = (xcg - xac + SM) / (Clah/ClaAh * (1-deda) * lh/c * (VhV)**2)
    return ShS

def scissorplot(CLh,CLAh,lh,c,VhV,Cmac,xac,Clah,ClaAh,deda,SM,xcg_fwBAR,xcg_aftBAR,PLOT=True,resolution=10**5,IgnoreErrors=False):
    """
    Calculates the optimum ratio of tail area 
    over wing area for a given aerodynamic 
    configuration, along with the unused c.g. range.

    Args:
        CLh (float): lift coefficient of the tail
        CLAh (float): lift coefficient of the aircraft without the tail
        lh (float): distance from the aerodynamic centre to the tail
        c (float): chord
        VhV (float): ratio of the velocity at the tail over the aircraft velocity
        Cmac (float): moment coefficient of the aerodynamic centre
        xac (float): location of the aerodynamic centre divided by the chord
        Clah (float): derivative of the lift coefficient of the tail w.r.t. alpha
        ClaAh (float): derivative of lift coefficient of the aircraft without the tail w.r.t alpha
        deda (float): derivative of the downwash angle w.r.t. alpha
        SM (float): stability margin
        xcg_fwBAR (float): location of the most forward position of the c.g. on the chord as a fraction of the chord
        xcg_aftBAR (float): location of the most aft position of the c.g. on the chord as a fraction of the chord
        PLOT (bool, optional): selects whether to display the resultant scissor plot. Defaults to True.
        resolution (integer, optional): resolution of the calculation. Defaults to 10**5.
        IgnoreErrors (bool, optional): selects whether to show and stop at errors, recommended to keep the default. Defaults to False.

    Returns:
        dictionary: _description_
    """    
    if xcg_fwBAR < 0 or xcg_fwBAR > 1:
        if not IgnoreErrors:
            raise ValueError(f'xcg_fwBAR currently has a value of {xcg_fwBAR} which means you are probably fucking up your design right now, place it back into the MAC.')

    if xcg_aftBAR < 0 or xcg_aftBAR > 1:
        if not IgnoreErrors:
            raise ValueError(f'xcg_aftBAR currently has a value of {xcg_fwBAR} which means you are probably fucking up your design right now, place it back into the MAC.')


    xcgBAR = np.linspace(0,1,resolution)
    id_xfw = np.where(np.abs(xcgBAR - xcg_fwBAR) == np.min(np.abs(xcgBAR - xcg_fwBAR)))[0][0]
    id_xaft = np.where(np.abs(xcgBAR - xcg_aftBAR) == np.min(np.abs(xcgBAR - xcg_aftBAR)))[0][0]
    cntr = controllability(CLh,CLAh,lh,c,VhV,xcgBAR,Cmac,xac)
    stab = stability(Clah,ClaAh,deda,lh,c,VhV,xcgBAR,xac,SM)
    
    if np.all(cntr - stab > 0): 
        if not IgnoreErrors:
            raise ValueError('No solutions in this configuration, controllability line both above and in same direction as stability line.')
    

    ShS_C_xfw = cntr[id_xfw]
    ShS_C_xaft = cntr[id_xaft]

    ShS_S_xfw = stab[id_xfw]
    ShS_S_xaft = stab[id_xaft]

    # Select the largest of the S_h/S ratios, 0 is controllability and 1 is stability
    ShSmax_id = np.where([ShS_C_xfw,ShS_S_xaft]==np.max([ShS_C_xfw,ShS_S_xaft]))[0][0]

    if PLOT:
        xi = np.nan

    if ShSmax_id == 0: # So a larger S_h/S ratio for controllability:
        id_ShS_S = np.where(np.abs(stab - ShS_C_xfw) == np.min(np.abs(stab - ShS_C_xfw)))[0][0]
        x_S = xcgBAR[id_ShS_S]
        xi = x_S
        if x_S < xcg_aftBAR:
            if not IgnoreErrors:
                raise ValueError('No solutions in this configuration, aft x_cgBAR is outside of stability region')
            
        deltaxcg =  x_S - xcg_aftBAR

    elif ShSmax_id == 1: # So a larger S_h/S ratio for stability:
        id_ShS_C = np.where(np.abs(cntr - ShS_S_xfw) == np.min(np.abs(cntr - ShS_S_xfw)))[0][0]
        x_C = xcgBAR[id_ShS_C]
        xi = x_C
        if x_C > xcg_fwBAR:
            if not IgnoreErrors:
                raise ValueError('No solutions in this configuration, aft x_cgBAR is outside of controllability region')
        
        deltaxcg = x_C - xcg_fwBAR
    
    else:
        raise ValueError('Something went wrong and im not sure what, this is not supposed to happen')        

    if PLOT:
        plt.plot(xcgBAR,cntr,label='Controllability')
        plt.plot(xcgBAR,stab,label='Stability')
        plt.plot([xcg_aftBAR,xcg_fwBAR],[ShS_C_xfw,ShS_S_xaft][ShSmax_id]*np.ones(2),label='c.g. Range')
        plt.plot([[xcg_aftBAR,xcg_fwBAR][ShSmax_id],xi],[ShS_C_xfw,ShS_S_xaft][ShSmax_id]*np.ones(2),label='Unused c.g. Range')
        plt.xlim((0,1))
        plt.ylim((0,np.max([np.max(stab),np.max(cntr)])))
        plt.legend()
        plt.tight_layout()
        plt.show()

    return {'ShS': [ShS_C_xfw,ShS_S_xaft][ShSmax_id],
            'xcg_waste': deltaxcg}

# k = scissorplot(CLh,CLAh,lh,c,VhV,Cmac,xacbar,Clah,ClaAh,deda,SM,xcg_fwBAR,xcg_aftBAR)
# print(k)