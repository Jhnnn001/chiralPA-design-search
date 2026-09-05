"""Project material fits; angular frequency in rad/s, time dependence exp(−iωt)."""
import numpy as np

C_NM_S = 2.99792458e17
HBAR_EVS = 6.582119569e-16

# TiO₂ fitted to data/tio2_measured.csv using fit_tio2.py.
TIO2_EPS_INF = 2.99058846
TIO2_OMEGA = np.array([5.00456470-0.32358247j, 5.52969166-1.05823931j,
                       6.17876796-1.22781886j, 7.94357666-0.05000000j,
                       9.68167024-1.38747835j, 11.60697072-1.73971782j])
TIO2_SIGMA = np.array([-0.03828249-0.02599826j, -0.82044484-2.81380588j,
                       5.19488014+7.06809000j, 0.00240131+0.00415690j,
                       -1.29288552-4.90438717j, -6.45882790+8.07588565j])

# Sehmi et al., Phys. Rev. B 95, 115444 (2017), Table III, Ag column.
AG_EPS_INF, AG_GAMMA, AG_SIGMA = 0.77259, 0.02228, 3751.4
AG_OMEGA = np.array([3.9173-0.06084j, 3.9880-0.04605j,
                     4.0746-0.63141j, 4.6198-2.8279j])
AG_SIGMA_L = np.array([0.09267+0.01042j, -0.0015342-0.062233j,
                       1.4911+0.40655j, 4.2843+4.2181j])

# Malitson, J. Opt. Soc. Am. 55, 1205 (1965), Eq. (1); wavelengths in µm.
SIO2_B = (0.6961663, 0.4079426, 0.8974794)
SIO2_C = (0.0684043, 0.1162414, 9.896161)


def lorentz(w, einf, poles, residues):
    w = np.asarray(w, complex)
    eps = np.full(w.shape, einf, complex)
    for pole, residue in zip(poles, residues):
        eps = eps+1j*residue/(w-pole)+1j*np.conj(residue)/(w+np.conj(pole))
    return eps


def tio2_epsilon(omega):
    return lorentz(np.asarray(omega)/1e15, TIO2_EPS_INF, TIO2_OMEGA, TIO2_SIGMA)


def ag_epsilon(omega):
    w = HBAR_EVS*np.asarray(omega, complex)
    return (lorentz(w, AG_EPS_INF, AG_OMEGA, AG_SIGMA_L)
            - AG_SIGMA*AG_GAMMA/(w*(w+1j*AG_GAMMA)))


def sio2_epsilon(omega):
    lam2 = (2*np.pi*C_NM_S*1e-3/np.asarray(omega, complex))**2
    return 1+sum(b*lam2/(lam2-c*c) for b, c in zip(SIO2_B, SIO2_C))


def permittivities(lam_nm):
    """Real-axis search model; remove the fit's small negative extinction floor."""
    omega = 2*np.pi*C_NM_S/lam_nm
    n = np.sqrt([f(omega) for f in (tio2_epsilon, sio2_epsilon, ag_epsilon)])
    return (n.real+1j*np.maximum(n.imag, 0))**2
