import numpy as np

def kroupa_cdf(value: float) -> float:
    """
    Kroupa CDF.
    
    Based on the Initial Mass Function in Kroupa et al. (1993).

    Parameters
    ----------
    value : float
        Value to be inversely mapped from [0.08,infty) solar masses to [0, 1).

    Raises
    ------
    ValueError
        For values outside of the interval [0.08,infty).

    Returns
    -------
    float
        Corresponding CDF value.
        
    References
    ----------
    Kroupa, P., Tout, C. A., and Gilmore, G., “The Distribution of Low-Mass Stars in the Galactic Disc”,
    Monthly Notices of the Royal Astronomical Society, vol. 262, OUP, pp. 545–587, 1993.
        doi:10.1093/mnras/262.3.545.
    """
    value = np.asanyarray(value)
    
    if np.any(value < 0.08) or np.any(value == np.inf):
        raise ValueError("Values must be >= to 0.08 but not infinte.")
        
    conds = [value < 0.5, (value >= 0.5) & (value < 1), value >= 1]
    funcs = [
        lambda x: 1.817 - 0.852 * x ** (-0.3),
        lambda x: 1.034 - 0.116 * x ** (-1.2),
        lambda x: 1.0 - 0.082 * x ** (-1.7)
        ]
    
    return np.select(conds, [f(value) for f in funcs])

def inverse_kroupa_cdf(value: float) -> float:
    """
    Function inverse to Kroupa CDF.
    
    Based on the Initial Mass Function in Kroupa et al. (1993).

    Parameters
    ----------
    value : float
        Value to be inversely mapped from (0,1) to (0.08, infty) solar masses.

    Raises
    ------
    ValueError
        For values outside of the interval (0,1).

    Returns
    -------
    float
        Corresponding solar mass.
        
    References
    ----------
    Kroupa, P., Tout, C. A., and Gilmore, G., “The Distribution of Low-Mass Stars in the Galactic Disc”,
    Monthly Notices of the Royal Astronomical Society, vol. 262, OUP, pp. 545–587, 1993.
        doi:10.1093/mnras/262.3.545.
    """
    value = np.asanyarray(value)
    
    if np.any(value < 0):
        raise ValueError("Values must be >= 0.")
    if np.any(value >= 1):
        raise ValueError("Values have be < 1.")
        
    conds = [value < 0.768, (value >= 0.768) & (value < 0.918), value >= 0.918]
    funcs = [
        lambda x: ((1.817 - x)/0.852)**(-1/0.3),
        lambda x: ((1.034 - x)/0.116)**(-1/1.2),
        lambda x: ((1.0 - x)/0.082)**(-1/1.7)
        ]  
    return np.select(conds, [f(value) for f in funcs])

def kroupa(lower: float, upper: float, sample_size: int = 1, rng: np.random.Generator = None) -> np.ndarray:
    """
    Generate random numbers based on Kroupa distribution. 
    
    Based on the Initial Mass Function defined in Kroupa et al. (1993).

    Parameters
    ----------
    lower : float
        Lower mass limit (in solar mass units). Minimum is 0.08.
    upper : float
        Upper mass limit (in solar mass units). Maximum is unbound but lower than infty.
    sample_size : int, optional
        Number of generated samples. The default is 1.
    rng : np.random.Generator, optional
        Random number generator, if it has precise seed, reproducibility is provided.

    Raises
    ------
    ValueError
        Lower bound must be greater than or equal to 0.08 solar masses.

    Returns
    -------
    np.ndarray
        Array of random masses sampled from the Kroupa distribution (in solar mass units).

    References
    ----------
    Kroupa, P., Tout, C. A., and Gilmore, G., “The Distribution of Low-Mass Stars in the Galactic Disc”,
    Monthly Notices of the Royal Astronomical Society, vol. 262, OUP, pp. 545–587, 1993.
        doi:10.1093/mnras/262.3.545.
    """
    if lower < 0.08:
        raise ValueError(f"The lower bound {lower} is smaller than 0.08.")
    if upper <= lower:
        raise ValueError("Upper bound must be > lower bound.")
    if upper == np.inf:
        raise ValueError("Upper bound must be finite.")
        
    if rng is None:
        rng = np.random.default_rng()
    
    unif_lower = kroupa_cdf(lower)
    unif_upper = kroupa_cdf(upper)
    
    # We will use inverse transform generation to draw samples
    unif = rng.uniform(unif_lower,unif_upper,sample_size)
    
    return inverse_kroupa_cdf(unif)
    
def tout_get_radius(mass: float) -> float:
    """
    Calculate radius of star given its mass.
    
    Based on zero-metallicity main-sequence mass-radius fitting formula of Tout et al. (1996).

    Parameters
    ----------
    mass : float
        Mass of the main-sequence star (in solar mass units).

    Returns
    -------
    float
        Radius of the star (in solar radius units?).

    References
    ----------
    Tout, C. A., Pols, O. R., Eggleton, P. P., and Han, Z., “Zero-age main-sequence radii and luminosities as analytic functions of mass and metallicity”,
    Monthly Notices of the Royal Astronomical Society, vol. 281, no. 1, OUP, pp. 257–262, 1996. 
        doi:10.1093/mnras/281.1.257.
    """
    numer = 1.71535900*mass**(2.5) + 6.59778800*mass**(6.5) + 10.08855000*mass**(11) + 1.01249500*mass**(19) + 0.07490166*mass**(19.5)
    denom = 0.01077422 + 3.0822340*mass**2 + 17.84778000*mass**(8.5) + mass**(18.5) + 0.00022582*mass**(19.5)
    return numer/denom