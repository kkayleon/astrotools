# Circular Restricted Three-Body Problem (CR3BP) dynamics

from astrotools.constants import G, mu_Earth, mu_Moon, r_Earth
import numpy as np
cos, sin, sqrt = np.cos, np.sin, np.sqrt

# Earth-moon system parameters (1: Earth, 2: Moon)
mu1, mu2 = mu_Earth, mu_Moon
r12 = 384400.0               # [km] Earth-Moon distance (characteristic length)
mu = mu1 + mu2               # [km^3/s^2] Characteristic gravitational parameter
pi1, pi2 = mu1/mu, mu2/mu    # [] Mass ratios

# Distance from barycenter (origin) to each primary
def r1r2(r):
    r1 = sqrt((r[0] + pi2)**2 + r[1]**2 + r[2]**2)
    r2 = sqrt((r[0] - pi1)**2 + r[1]**2 + r[2]**2)
    return r1, r2

# Non-dimensionalized EoM in rotating frame
def acceleration(r, v):
    r1, r2 = r1r2(r)

    # Equations of motion
    ax = 2*v[1] + (1 - pi1/r1**3 - pi2/r2**3)*r[0] + (-pi2*pi1/r1**3 + pi1*pi2/r2**3)
    ay = -2*v[0] + (1 - pi1/r1**3 - pi2/r2**3)*r[1] 
    az = -(pi1/r1**3 + pi2/r2**3)*r[2]

    return np.array([ax, ay, az])

# Jacobi constant
def jacobi_constant(r, v):
    r1, r2 = r1r2(r)
    
    # Pseudo-potential 
    U_pseudo = 0.5*(r[0]**2 + r[1]**2) + (pi1/r1 + pi2/r2)

    return 0.5*np.dot(v,v) - U_pseudo