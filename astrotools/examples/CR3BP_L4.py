# Example: Circular restricted 3-body dynamics around Earth-Moon L4

from astrotools.dynamics.cr3bp import jacobi_constant, pi1, pi2
from astrotools.trajectory import cr3bp_trajectory
import numpy as np
cos, sin, sqrt = np.cos, np.sin, np.sqrt

# Initial conditions
r0 = np.array([0.5 - pi2, sqrt(3)/2, 0.0])
v0 = np.array([0.0, 0.0, 0.0])

# Propagate trajectory for 10 non-dimensional time units with dt=0.01
n = 1000
dt = 0.01
trajectory = cr3bp_trajectory.trajectory(r0, v0, n, dt)

# Jacobi constant along trajectory
jacobi_constants = np.array([jacobi_constant(state[:3], state[3:6]) for state in trajectory])