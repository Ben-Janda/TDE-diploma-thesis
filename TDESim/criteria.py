import numpy as np
from TDESim.config import SimulationParams
from TDESim.utils import rsun_to_rg, msun_to_rg

""" Intersection criterion - Guillochon 2015 """
# Width of the plasmatic stream, eq. (2)
def S(params: SimulationParams, points, W):
    r = np.linalg.norm(points, axis = 0)
    return W*params.beta*r*params.q**(-1/3)

def check_intersection(params: SimulationParams, points1, points2, W1, W2):
    """
    Test intersection of two ellipses and return first intersection point
    """
    
    D = np.linalg.norm(points1 - points2, axis = 0) - S(params, points1, W1) - S(params, points2, W2)
    
    if np.any(D < 0):
        return {'test': True, 'point': np.transpose(points1)[np.where(D < 0)][0]}
    
    return {'test': False, 'point': None}

def check_IBCO(beta: float, i: float, a: float, Mh: float, Rstar: float, Mstar: float) -> bool:
    """
    Check whether the trajectory is at or inside innermost bound circular orbit (IBCO), from simulation parameters dataclass.
    
    Based on the condition of Bardeen et al. (1972), equation (2.19).

    Parameters
    ----------
    beta : float
        Impact parameter of the TDE.
    i : float
        Inclination in rad from (0, pi).
    a : float
        Dimensionless spin of the black hole in (0,1).
    Mh : float
        Black hole mass in Msun.
    Rstar : float
        Star radius in Rsun units.
    Mstar : float
        Star mass in Msun units.

    Raises
    ------
    ValueError
        Prevents incorrect inclination input.

    Returns
    -------
    bool
        True if the trajectory is at or inside IBCO, False otherwise.

    References
    ----------
    Bardeen, J. M., Press, W. H., and Teukolsky, S. A., “Rotating Black Holes: Locally Nonrotating Frames, Energy Extraction, and Scalar Synchrotron Radiation”,
    The Astrophysical Journal, vol. 178, IOP, pp. 347–370, 1972. 
        doi:10.1086/151796.
    """
    if 0 <= i < np.pi/2 :
        # Prograde = co-rotating orbit
        r_IBCO = msun_to_rg((2 - a  + 2 * (1 - a)**(1/2))*Mh, Mh)
    elif np.pi/2 < i <= np.pi:
        # Retrograde = counter-rotating orbit
        r_IBCO = msun_to_rg((2 + a  + 2 * (1 + a)**(1/2))*Mh, Mh)
    else:
        raise ValueError(f"Inclination i = {i} must be from (0,pi).")
    
    rp = rsun_to_rg(Rstar, Mh) * (Mh/Mstar)**(1/3) / beta
    
    if rp <= r_IBCO:
        return True
    else:
        return False