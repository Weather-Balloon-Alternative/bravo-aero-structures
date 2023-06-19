import matplotlib as plt
from main_sizing import *

def plot_eigvals_symm():

    eigvals = eigenvalues_ss_symEOM(muc,c,V0,CZadot,Cmadot,KY2,CXu,CXa,CZ0,CXq,CZu,CZa,CX0,CZq,Cmu,Cma,Cmq,CXde,CZde,Cmde)
    eigvals_ph = eigenvalues_phugoid(muc,CZa,Cmq,Cma,CXu,Cmu,CXa,CZu,CZ0,V0,c)
    eigvals_sp = eigenvalues_shortperiod(muc,KY2,CZadot,CZq,Cmadot,Cmq,CZa,Cma,V0,c)

    plt.scatter(np.real(eigvals),np.imag(eigvals),label='Eigenvalues')
    plt.scatter(np.real(eigvals_ph),np.imag(eigvals_ph),label='Estimated Eigenvalues Phugoid')
    plt.scatter(np.real(eigvals_sp),np.imag(eigvals_sp),label='Estimated Eigenvalues Short Period')
    plt.xlabel('Real axis')
    plt.ylabel('Imaginary axis')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()

    # if np.any(np.abs(eigvals)>100):
    #     eigvals = eigvals[np.where(np.abs(eigvals)<100)]
    #     eigvals_ph = eigvals_ph[np.where(np.abs(eigvals_ph)<100)]
    #     eigvals_sp = eigvals_sp[np.where(np.abs(eigvals_sp)<100)]

    # plt.scatter(np.real(eigvals),np.imag(eigvals),label='Eigenvalues')
    # plt.scatter(np.real(eigvals_ph),np.imag(eigvals_ph),label='Estimated Eigenvalues Phugoid')
    # plt.scatter(np.real(eigvals_sp),np.imag(eigvals_sp),label='Estimated Eigenvalues Short Period')
    # plt.xlabel('Real axis')
    # plt.ylabel('Imaginary axis')
    # plt.grid()
    # plt.legend()
    # plt.tight_layout()
    # plt.show()

def plot_eigvals_asymm():
    eigvals = eigenvalues_ss_asymEOM(CYbdot,mub,b,V0,KX2,KXZ,KZ2,CYb,CL,CYp,CYr,Clb,Clp,Clr,Cnb,Cnp,Cnr,CYda,CYdr,Clda,Cldr,Cnda,Cndr)
    eigvals_ap = eigenvalues_aperiodicroll(Clp,mub,KX2,V0,b)
    eigvals_dr = eigenvalues_dutchroll(mub,KZ2,Cnr,Cnb,CYb,V0,b)
    eigvals_s = eigenvalues_spiral(CL,Clb,Cnr,Cnb,Clr,Clp,CYb,mub,Cnp,V0,b)

    plt.scatter(np.real(eigvals),np.imag(eigvals),label='Eigenvalues')
    plt.scatter(np.real(eigvals_ap),np.imag(eigvals_ap),label='Estimated Eigenvalues Aperiodic Roll')
    plt.scatter(np.real(eigvals_dr),np.imag(eigvals_dr),label='Estimated Eigenvalues Dutch roll')
    plt.scatter(np.real(eigvals_s),np.imag(eigvals_s),label='Estimated Eigenvalues Aperiodic Spiral')
    plt.xlabel('Real axis')
    plt.ylabel('Imaginary axis')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()
if __name__ == '__main__':
    from coefficients_smallglider import *
    # from coefficients_bigglider import *

    plot_eigvals_symm()
    # plot_eigvals_asymm()