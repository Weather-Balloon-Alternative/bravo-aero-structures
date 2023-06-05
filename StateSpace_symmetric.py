import numpy as np
import scipy.linalg as la
import control
# from coefficients import *

def ss_sym(muc,c,V0,CZadot,Cmadot,KY2,CXu,CXa,CZ0,CXq,CZu,CZa,CX0,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde):
    """Generates the state-space form of the symmetric equations of motion. Outputs a state-space from the Python control library.

    Args:
        muc (float): relative density [=m/(rho*S*c)]
        c (float): MAC
        V0 (float): aircraft speed
        CZadot (float): derivative of Cz w.r.t. alpha_dot*c/V0
        Cmadot (float): derivative of Cm w.r.t. alpha_dot*c/V0
        KY2 (float): square of the nondimensional radius of gyration about the Y-axis [=I_y/(m*b)]
        CXu (float): derivative of X w.r.t. u divided by 1/(0.5*rho*V*S)
        CXa (float): derivative of CX w.r.t. alpha
        CZ0 (float): CZ in CX w.r.t. qc/V0
        CZu (float): derivative of Z w.r.t. u divided by 1/(0.5*rho*V*S)
        CZa (float): derivative of CZ w.r.t. alpha
        CX0 (float): CX in steady flight
        CZq (float): derivative of CZ w.r.t. qc/V0
        Cmu (float): derivative of M w.r.t. divided by 1/(0.5*rho*V*S)
        Cma (float): derivative of Cm w.r.t. alpha
        Cmq (float): derivative of Cm w.r.t. qc/V0
        CXde (float): derivative of CX w.r.t. delta_e
        CZde (float): derivative of CZ w.r.t. delta_e
        Cmde (float): elevator efficiency, derivative of Cm w.r.t. delta_e

    Returns:
        control.statesp.StateSpace: State-space of the symmetric equation of motion
    """   
    # Define P,Q and R matrices 
    P = np.matrix([ [-2*muc*(c/(V0**2)),    0,                      0,      0                       ],
                    [0,                     (CZadot-2*muc)*c/V0,    0,      0                       ],
                    [0,                     0,                      -c/V0,  0                       ],
                    [0,                     Cmadot*c/V0,            0,      -2*muc*KY2*(c/V0)**2    ]])
    
    Q = np.matrix([ [-CXu/V0,               -CXa,                   -CZ0,   CXq*c/V0                ],
                    [-CZu/V0,               -CZa,                   CX0,    -(CZq+2*muc)*c/V0       ],
                    [0,                     0,                      0,      -c/V0                   ],
                    [-Cmu/V0,               -Cma,                   0,      -Cmq*c/V0               ]])
    
    R = np.matrix([ [-CXde  ],
                    [-CZde  ],
                    [0      ],
                    [-Cmde  ]])
    
    # Compute inverse of P and calculate the A and B matrices
    Pinv= la.inv(P)

    A = Pinv@Q
    B = Pinv@R

    # Right now we just want the x vector as output, so C is the 
    # identity matrix and D is the null matrix
    C= np.identity(4)
    D= np.zeros((4,1))

    # Create the control state-space system
    sys = control.ss(A, B, C, D)
    return sys