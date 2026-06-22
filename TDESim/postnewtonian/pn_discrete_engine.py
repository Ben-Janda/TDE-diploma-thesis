import numpy as np
from numpy import pi, cos, sin
from scipy.spatial.transform import Rotation as R

from typing import Tuple

from TDESim.config import SimulationParams

# Note:
#    
# Equation numbers refer to the paper
#
# Guillochon, James, and Enrico Ramirez-Ruiz. 
# "A dark year for tidal disruption events." 
# The Astrophysical Journal 809.2 (2015): 166.

""" Precession calculation """
def precession(params: SimulationParams) -> Tuple[float, float]:
    """
    Calculate apsidal and nodal precession averaged over one orbit.
    
    Based on post-Newtonian approximation in Kerr spacetime, as in Guillochon et al. (2015).

    Parameters
    ----------
    params : SimulationParams
        Parameters of the simulation.

    Raises
    ------
    ValueError
        Expressions diverge for inclination 'params.i' close to pi/2 and -pi/2.

    Returns
    -------
    aspidal : float
        Apsidal precession per orbit (in radians).
    nodal : float
        Nodal precession per orbit (in radians).

    References
    ----------
    Guillochon, J. and Ramirez-Ruiz, E., “A Dark Year for Tidal Disruption Events”,
    The Astrophysical Journal, vol. 809, no. 2, Art. no. 166, IOP, 2015.
        doi:10.1088/0004-637X/809/2/166.
    """
    # check that parameters are safe for PN simulation
    if abs(params.i - pi/2) < 1e-10:
        raise ValueError(
            f"Inclination i = {params.i} rad is too close to (pi/2), causing divergence."
        )
    if abs(params.i + pi/2) < 1e-10:
        raise ValueError(
            f"Inclination i = {params.i} rad is too close to (-pi/2), causing divergence."
        )
    
    prec = params.Mh*params.beta/(1+params.e)/params.Rstar
    
    # nodal precession, eqs. (11), (13), (14)
    omega_J = (4*pi*params.a/params.q**(1/2))*prec**(3/2) # Lense-Thirring
    omega_Q = (3*pi*params.a**2/params.q**(2/3))*prec**2 * cos(params.i) # Quadrupole
    domega = omega_J + omega_Q

    # apsidal precession, eqs. (12), (15-17)
    Omega_D = 6*pi*prec/params.q**(1/3) # de-Sitter
    Omega_J = -4*omega_J*cos(params.i) # Lense-Thirring
    Omega_Q = (1-5*cos(params.i)**2)*omega_Q/(2*cos(params.i)) # Quadrupole
    dOmega = Omega_D + Omega_J + Omega_Q
    
    return dOmega, domega

""" Transform parametric curve to 3D cartesian axis """

# initial ellipse - it is in each iteration same up to the orientation in space
def init_ellipse(params: SimulationParams) -> np.ndarray:
    """
    Generate ellipse datapoints corresponding to simulation parameters, but in 2D plane.

    Parameters
    ----------
    params : SimulationParams
        Parameters of the simulation.

    Returns
    -------
    np.ndarray
        Array of datapoints in 2D space.

    """
    phi = np.linspace(0,2*pi, params.res)
    r = params.rp*(1+params.e)/(1+params.e*cos(phi))
    
    x = r*cos(phi)
    y = r*sin(phi)
    z = np.zeros_like(r)
    
    return np.column_stack((x,y,z))

def evolve_orbit(params: SimulationParams, ellipse: np.ndarray, omega_apsidal: float, omega_nodal: float) -> np.ndarray:
    """
    Evolve the orbit based on initial ellipse and total precession values.

    Parameters
    ----------
    params : SimulationParams
        Parameters of the simulation.
    ellipse : np.ndarray
        Set of initial ellipse datapoints.
    omega_apsidal : float
        Total apsidal precession.
    omega_nodal : float
        Total nodal precession.

    Returns
    -------
    np.ndarray
        Array of datapoints in 3D space corresponding to the orbit.

    """
    # define rotation matrix for correct spatial orientation of ellipse plane
    rot = R.from_euler('zyz', [omega_apsidal,-params.i, omega_nodal], degrees = False) 
    return rot.apply(ellipse)            