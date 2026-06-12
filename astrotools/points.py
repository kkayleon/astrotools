# Lagrange/Libration points for circular restricted three-body

from astrotools.dynamics.cr3bp import pi1, pi2
from scipy.optimize import brentq
import numpy as np
array, sqrt = np.array, np.sqrt

# Tolerance for numerical stability
tol = 1E-6

# Equation of motion (x) = 0
def f(x):
    r1 = abs(x + pi2)
    r2 = abs(x - pi1)
    return (1 - pi1/r1**3 - pi2/r2**3)*x + (-pi2*pi1/r1**3 + pi1*pi2/r2**3) 

# Lagrange/Libration points
def l_points():
    # Brent's method of solving for root within some defined interval
    L1 = array([brentq(f, -pi2 + tol, pi1 - tol), 0, 0])
    L2 = array([brentq(f, pi1 + tol, 2), 0, 0])
    L3 = array([brentq(f, -2, -pi2 - tol), 0, 0])
    L4 = array([(pi1 - pi2)/2, sqrt(3)/2, 0])
    L5 = array([(pi1 - pi2)/2, -sqrt(3)/2, 0])
    return array([L1, L2, L3, L4, L5])
