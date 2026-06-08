# Two-body dynamics w/ J2 perturbation

import numpy as np
from astrotools.constants import J2_Earth, r_Earth
cos, sin, pi, sqrt = np.cos, np.sin, np.pi, np.sqrt

# Perturbed gravitational acceleration (2BP + J2)
def acceleration(r, mu):
    r_mag = np.linalg.norm(r)

    a_2B = -mu*r/r_mag**3
    a_J2 = 1.5*J2_Earth*mu*r_Earth**2/r_mag**5 * np.array([r[0]*(5*(r[2]/r_mag)**2 - 1), r[1]*(5*(r[2]/r_mag)**2 - 1), r[2]*(5*(r[2]/r_mag)**2 - 3)])
    
    return a_2B + a_J2

# Perturbed potential energy (2BP + J2)
def potential_J2(r, mu):
    U_2B = -mu/np.linalg.norm(r)
    U_J2 = -0.5*J2_Earth*mu*r_Earth**2/np.linalg.norm(r)**3 * (3*(r[2]/np.linalg.norm(r))**2 - 1)

    return U_2B + U_J2