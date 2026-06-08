# Trajectory propagation

from astrotools.dynamics import twobody, j2
from astrotools.constants import mu_Earth, r_Earth, J2_Earth

# Select numerical integrator (WIP/Future works)
from astrotools.integrators import yoshida4

# R2BP Acceleration (no perturbation)
func_accel_2B = lambda r: twobody.acceleration(r, mu_Earth)

# R2BP + J2 perturbation Acceleration
func_accel_J2 = lambda r: j2.acceleration(r, mu_Earth)

def trajectory(r0, v0, n, dt, perturbation=False):
    if perturbation:
        return yoshida4.propagate(r0, v0, n, dt, func_accel_J2)
    else:
        return yoshida4.propagate(r0, v0, n, dt, func_accel_2B)
