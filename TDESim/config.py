from dataclasses import dataclass
import json

from .utils import rsun_to_rg, msun_to_rg

@dataclass
class SimulationParams:
    """
    Parameters for TDE simulation. Inputs specified below, internal conversion to geometrized with Mh = 1.
    """
    # Primary parameters (Note: rg = G*Mh/c^2)
    Rstar: float    # [Rsun] Star radius in solar radius units -> [rg]
    Mstar: float    # [Msun] Star mass in solar mass units -> [rg]
    Mh: float       # [Msun] Black hole mass in solar mass units -> [rg]
    rp: float       # [rg] Periapsis distance from black hole
    a: float        # [-] Dimensionless spin in [0,1].
    e: float        # [-] Numerical eccentricity in (0,1)
    i: float        # [rad] Inclination of the orbit, measured from equatorial plane, in [-pi, pi].
    res: int        # [-] Angular resolution (number of points along 2*pi)
    
    # Secondary parameters
    q: float = None        # [-] Mass ratio (Black hole : Sun)
    rt: float = None       # [rg] Tidal radius
    beta: float = None     # [-] (rt : rp) ratio
    
    def __post_init__(self):
        """
        Calculate values of secondary parameters and convert everything to correct units.
        """
        # Check for incorrect
        if not (0 <= self.a <= 1):
            raise ValueError(
                f"Value of dimensionless spin a = {self.a} is not in [0,1]."
                )
        if not (0 <= self.e < 1):
            raise ValueError(
                f"Value of eccentricity e = {self.e} is not in [0,1). Simulation works only for bound orbits."
                )
        
        # Calculate secondary
        if self.q is None:
            self.q = self.Mh / self.Mstar
        if self.rt is None:
            self.rt = rsun_to_rg(self.Rstar, self.Mh) * self.q **(1/3)
        if self.beta is None:
            self.beta = self.rt / self.rp
            
        # Convert to correct units
        self.Rstar = rsun_to_rg(self.Rstar, self.Mh)
        self.Mstar = msun_to_rg(self.Mstar, self.Mh)
        self.Mh = 1
            
    @classmethod
    def load_from_json(cls, file_path: str):
        """
        Load simulation parameters from JSON file and create class instance
        """
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        return cls(**data)        