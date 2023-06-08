from StateSpace_symmetric import *
from StateSpace_asymmetric import *

# ===================================
# ===================================

# Symmetrical

# ===================================
# Obtain eigenvalues of the entire statespace
def eigenvalues_ss_symEOM(muc,c,V0,CZadot,Cmadot,KY2,CXu,CXa,CZ0,CXq,CZu,CZa,CX0,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde):
    sys = ss_symEOM(muc,c,V0,CZadot,Cmadot,KY2,CXu,CXa,CZ0,CXq,CZu,CZa,CX0,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde)
    A = sys.A
    eigenvals = la.eig(A)
    return eigenvals[0]

# Obtain eigenvalues according to analytical equations
#### Short period
def eigenvalues_shortperiod(muc,KY2,CZadot,CZq,Cmadot,Cmq,CZa,Cma,V0,c):
    A = 2*muc*KY2*(2*muc-CZadot)
    B = -2*muc*KY2*CZa-(2*muc+CZq)*Cmadot-(2*muc-CZadot)*Cmq
    C = CZa*Cmq-(2*muc+CZq)*Cma
    eigenvals = np.array(((-1*B+np.emath.sqrt(B**2-4*A*C))/(2*A),(-1*B-np.emath.sqrt(B**2-4*A*C))/(2*A)))
    return eigenvals*V0/c

#### Phugoid
def eigenvalues_phugoid(muc,CZa,Cmq,Cma,CXu,Cmu,CXa,CZu,CZ0,V0,c):
    A = 2*muc*(CZa*Cmq-2*muc*Cma)
    B = 2*muc*(CXu*Cma-Cmu*CXa)+Cmq*(CZu*CXa-CXu*CZa)
    C = CZ0*(Cmu*CZa-CZu*Cma)
    eigenvals = np.array(((-1*B+np.emath.sqrt(B**2-4*A*C))/(2*A),(-1*B-np.emath.sqrt(B**2-4*A*C))/(2*A)))
    return eigenvals*V0/c

# ===================================
# ===================================

# Asymmetrical

# ===================================
def eigenvalues_ss_asymEOM(CYbdot,mub,b,V0,KX2,KXZ,KZ2,CYb,CL,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr):
    sys = ss_asymmetricEOM(CYbdot,mub,b,V0,KX2,KXZ,KZ2,CYb,CL,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr)
    A = sys.A
    eigenvals = la.eig(A)
    return eigenvals[0]

# Obtain eigenvalues according to analytical equations
#### Aperiodic Roll
def eigenvalues_aperiodicroll(Clp,mub,KX2,V0,b):
    eigenvals = np.array((Clp / (4*mub*KX2)))
    return eigenvals*V0/b

#### Dutch Roll
def eigenvalues_dutchroll(mub,KZ2,Cnr,Cnb,CYb,V0,b):
    A = 8*mub**2*KZ2
    B = -2*mub*(Cnr+2*KZ2*CYb)
    C = 4*mub*Cnb+CYb*Cnr
    eigenvals = np.array(((-1*B+np.emath.sqrt(B**2-4*A*C))/(2*A),(-1*B-np.emath.sqrt(B**2-4*A*C))/(2*A)))
    return eigenvals*V0/b

#### Spiral
def eigenvalues_spiral(CL,Clb,Cnr,Cnb,Clr,Clp,CYb,mub,Cnp,V0,b):
    num = 2*CL*(Clb*Cnr-Cnb*Clr)
    den = Clp*(CYb+4*mub*Cnb)-Cnp*(CYb*Clr+4*mub*Clb)
    eigenvals = np.array((num/den))
    return eigenvals*V0/b

# ===================================
# ===================================