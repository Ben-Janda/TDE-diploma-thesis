G = 6.67430e-11         # [m^3/kg/s^2] Gravitation constant
c = 2.99792458e8        # [m/s] Light speed
Msun = 1.989e30         # [kg] Solar mass
rsun = 6.96e8           # [m] Solar radius
rgsun = (G*Msun)/(c**2) # [m] Gravitation radius of Sun (r = GM/c^2)

def rsun_to_rg(distance_rsun: float, mass_bh_msun: float) -> float:
    """
    Convert distance in rsun to distance in black-hole gravitation radii.
    """
    rg_bh_m = rgsun*mass_bh_msun
    dist_m = distance_rsun*rsun
    return dist_m/rg_bh_m

def msun_to_rg(mass_star_msun: float, mass_bh_msun: float) -> float:
    """
    Convert mass in msun to distance in black-hole gravitation radii.
    """
    return mass_star_msun/mass_bh_msun
    