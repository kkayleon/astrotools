# Trajectory propagation (restricted two-body)

from astrotools.dynamics import twobody, j2
from astrotools.constants import mu_Earth
from astrotools.integrators import yoshida4, dopr853

# R2BP Acceleration (no perturbation)
func_accel_2B = lambda r: twobody.acceleration(r, mu_Earth)

# R2BP + J2 perturbation Acceleration
func_accel_J2 = lambda r: j2.acceleration(r, mu_Earth)

# Propagate trajectory for R2BP with optional J2 perturbation    
def trajectory(r0, v0, n, dt, perturbation=False, solverType="dopr853"):
    # Perturbation model selection
    if perturbation:
        func_accel = func_accel_J2
    else:
        func_accel = func_accel_2B

    # Solver selection
    if solverType.lower() == "yoshida4":
        return yoshida4.propagate(r0, v0, n, dt, func_accel)
    elif solverType.lower() == "dopr853":
        return dopr853.propagate(r0, v0, n, dt, func_accel)
    else: 
        # Default to DOPR853 if solverType not recognized
        return dopr853.propagate(r0, v0, n, dt, func_accel)
