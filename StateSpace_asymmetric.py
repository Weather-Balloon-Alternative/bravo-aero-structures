import numpy as np
import scipy.linalg as la
import control

def ss_asymmetricEOM(CYbdot,mub,b,V0,KX2,KXZ,KZ2,CYb,CL,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr):
    """Generates the state-space form of the asymmetric
       equations of motion. Outputs a state-space from
       the Python control library.

    Args:
        CYbdot (float): derivative of CY w.r.t. beta_dot
        mub (float): relative density [=m/(rho*S*b)]
        b (float): wing span
        V0 (float): aircraft speed
        KX2 (float): square of the nondimensional radius of gyration about the X-axis [=I_x/(m*b)]
        KXZ (float): nondimensional product of inertia [=J_xz/(m*b)]
        KZ2 (float): square of the nondimensional radius of gyration about the Z-axis [=I_z/(m*b)]
        CYb (float): derivative of CY w.r.t. beta
        CL (float): lift coefficient
        CYp (float): derivative of CY w.r.t. p*b/(2*V0)
        CYr (float): derivative of CY w.r.t. r*b/(2*V0)
        Clb (float): derivative of Cl w.r.t. beta
        Clp (float): derivative of Cl w.r.t. p*b/(2*V0)
        Clr (float): derivative of Cl w.r.t. r*b/(2*V0)
        Cnb (float): static directional stability, derivative of Cn w.r.t. beta
        Cnp (float): derivative of Cn w.r.t. p*b/(2*V0)
        Cnr (float): derivative of Cn w.r.t. r*b/(2*V0)
        CYda (float): derivative of CY w.r.t. delta_a
        CYdr (float): derivative of CY w.r.t. delta_r
        Clda (float): derivative of Cl w.r.t. delta_a
        Cldr (float): derivative of Cl w.r.t. delta_r
        Cnda (float): derivative of Cn w.r.t. delta_a
        Cndr (float): derivative of Cn w.r.t. delta_r

    Returns:
        _type_: _description_
    """    
    # Define P, Q and R matrices
    P = np.mat([[(CYbdot-2*mub)*b/V0,   0,          0,                              0                               ],
                [0,                     -b/(2*V0),  0,                              0                               ],
                [0,                     0,          -4*mub*KX2*0.5*(b/V0)**2,       4*mub*KXZ*0.5*(b/V0)**2         ],
                [0,                     0,          4*mub[id]*KXZ*0.5*(b/V0)**2,    -4*mub[id]*KZ2*0.5*(b/V0)**2    ]])
    
    Q = np.mat([[-1*CYb,                -1*CL[id],  -1*CYp*b/(2*V0),                -(CYr-4*mub)*b/(2*V0)           ],
                [0,                     0,          -1*(b/(2*V0)),                  0                               ],
                [-1*Clb,                0,          -1*Clp*0.5*(b/V0),              -1*Clr*0.5*(b/V0)               ],
                [-1*Cnb,                0,          -1*Cnp*0.5*b/V0,                -1*Cnr*0.5*b/V0                 ]])
    
    R = np.mat([[-1*CYda,   -1*CYdr ],
                [0,         0       ],
                [-1*Clda,   -1*Cldr ],
                [-1*Cnda,   -1*Cndr ]])
    
    # Compute inverse of P and calculate the A and B matrices
    Pinv = la.inv(P)

    A = Pinv @ Q
    B = Pinv @ R

    # Right now we just want the x vector as output, so C is the 
    # identity matrix and D is the null matrix

    C = np.identity(4)
    D = np.zeros((4,2))

    # Create the control state-space system
    sys = control.ss(A,B,C,D)
    return sys