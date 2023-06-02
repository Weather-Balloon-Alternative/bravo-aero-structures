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

def scissorplot(CLh,CLAh,lh,c,VhV,xcg,Cmac,xac,Clah,ClaAh,deda,SM,delta_xcgbar,PLOT=True):
    '''
    Scissor Plot, returns the ratio of tail area over wing area

    args:

    CLh:    float, C_L_h;   lift coefficient of the tail
    CLAh:   float, C_L_A-h; lift coefficient of the aircraft without the tail 
    lh:     float, l_h;     distance from the aerodynamic centre to the tail
    c:      float, c;       chord
    VhV:    float, V_h/V;   ratio of the velocity at the tail over the aircraft velocity
    xcg:    float, x_cgBAR; location of the centre of gravity divided by the chord
    Cmac:   float, C_m_ac;  moment coefficient of the aerodynamic centre
    xac:    float, x_acBAR; location of the aerodynamic centre divided by the chord
    CLah:   float, C_L_alpha_h;     derivative of the lift coefficient of the tail w.r.t. alpha
    CLaAh:  float, C_L_alpha_A-h;   derivative of lift coefficient of the aircraft without the tail w.r.t alpha
    deda:   float, deta/dalpha;     derivative of the downwash angle w.r.t. alpha
    SM:     float, S.M.;            stability margin
    
    returns:

    ShS:    float, S_h/S;   ratio of tail area over wing area
    
    '''
    cntr = controllability(CLh,CLAh,lh,c,VhV,xcg,Cmac,xac)
    stab = stability(Clah,ClaAh,deda,lh,c,VhV,xcg,xac,SM)