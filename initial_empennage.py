import numpy as np
def tailsizing(b,c,ShS=0.15,Arh=4,lhc=4,Arv = 1.1):
    """Initial empennage sizing for a H-tail glider. 
    Initial area ratios and aspect ratios taken from 
    historical gliders and H-tail aircraft.
    (References: 
    Aircraft Design : A Systems Engineering Approach - M. Sadraey
    Jane's All the World's Aircraft
    Fundamentals of Sailplane Design - F. Thomas
    )

    Args:
        b (float): Main wing span.
        c (float): Main wing MAC.
        ShS (float, optional): Ratio of horizontal tail area over main wing area. Defaults to 0.1.
        Arh (int, optional): Aspect ratio of the horizontal tail. Defaults to 4.
        lhc (int, optional): Ratio of horizontal tail length over MAC. Defaults to 4.
        Arv (float, optional): Aspect ratio of the vertical tail, note that this is chord over half-span. Defaults to 1.1.

    Returns:
        dictionary: Basic empennage parameters,
        'l_h' (float): Tail length, equal for both horizontal and vertical tail.
        'S_h' (float): Horizontal tail area.
        'b_h' (float): Horizontal tail span.
        'c_h' (float): Horizontal tail MAC.
        'S_v' (float): Vertical tail area.
        'b_v' (float): Vertical tail half-span.
        'c_v' (float): Vertical tail MAC.

    """    
    # Initialize
    tailpar = {}

    # Global tail parameters
    tailpar['l_h'] = lhc * c # Tail length, assumed to be the same for hor. and vert. due to H-tail config
    
    # Horizontal tail parameters
    tailpar['S_h'] = b*c*ShS # Horizontal tail area
    tailpar['b_h'] = np.sqrt(Arh*tailpar['S_h']) # Horizontal tail span
    tailpar['c_h'] = tailpar['b_h'] / Arh # Horizontal tail MAC
    
    # Vertical tail parameters
    # Assume that c_v is equal to c_h due to H-tail config
    tailpar['c_v'] = tailpar['c_h'] # Vertical tail MAC
    tailpar['b_v'] = tailpar['c_v'] * Arv # Vertical tail span, note that this is officially half-span!
    tailpar['S_v'] = tailpar['b_v'] * tailpar['c_v']

    return tailpar

def masscalc(b,c,density,A_c):
    """Mass calculation of a basic wing section

    Args:
        b (float): wing span
        c (float): MAC
        density (float): materical density
        A_c (float): ratio of airfoil area over the MAC

    Returns:
        float: wing section mass
    """   
    m = A_c*c*b*density
    return m 

if __name__ == '__main__':
    b = 1.023609984
    # S = 0.0876
    c = 0.087924113
    # S = 0.336
    # b = 2.246
    # c = S/b
    A_c = 0.07995 # For NACA0012
    rho = 75

    tailpar = tailsizing(b,c)
    print(tailpar)
    mhor = masscalc(tailpar['b_h'],tailpar['c_h'],rho,A_c)
    mver = masscalc(tailpar['b_v'],tailpar['c_v'],rho,A_c)
    print(mhor,mver)
    print(f'total mass: {(mhor+mver)*1000} g')