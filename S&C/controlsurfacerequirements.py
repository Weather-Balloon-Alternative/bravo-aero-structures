import numpy as np
import ambiance
# from coefficients import *

def pitchraterequirements(V0,h,nmax,radians=True):
    """Generate the pitch rate requirements

    Args:
        V0 (float): aircraft velocity
        h (int): aircraft altitude
        nmax (float): maximum design load factor
        radians (bool, optional): whether to return in rad/s. Defaults to True.

    Returns:
        ang_acc: required angular acceleration
    """    
    atmos = ambiance.Atmosphere(h)
    g0 = atmos.grav_accel
    acmax = nmax*g0
    r = V0**2 / acmax
    if radians:
        ang_acc = np.sqrt(acmax/r)
    else:
        ang_acc = np.sqrt(acmax/r) * 180 / np.pi

    return {'ang_acc': ang_acc}

# if __name__ == '__main__':
    # p = pitchraterequirements(V0,h,nmax,radians=False)