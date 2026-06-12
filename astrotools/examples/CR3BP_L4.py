# Example: Circular restricted 3-body dynamics around Earth-Moon L4
# !! Analysis done in nondimensional units !!

from astrotools.dynamics.cr3bp import jacobi_constant, pi1, pi2
from astrotools.trajectory import cr3bp_trajectory
from astrotools.points import l_points
from astrotools.plotting import cr3bp_orbit
import matplotlib.pyplot as plt
import numpy as np
cos, sin, sqrt = np.cos, np.sin, np.sqrt

# Initial conditions
perturb_r = np.array([0.015, 0.0, 0.0])
r0 = np.array([0.5 - pi2, sqrt(3)/2, 0.0]) + perturb_r
v0 = np.array([0.0, 0.0, 0.0])

# Propagate trajectory
n = 10000
dt = 0.01
trajectory = cr3bp_trajectory.trajectory(r0, v0, n, dt)

# Jacobi constant along trajectory
jacobi_constants = np.array([jacobi_constant(state[:3], state[3:6]) for state in trajectory])

# Lagrange/Libration points from points.py
libration_points = l_points()

# Output
print(f"")
print(f"-----------------------------------------------------------------"), print(f"")
for i, name in enumerate(['L1', 'L2', 'L3', 'L4', 'L5']):
    print(f"{name}: [{libration_points[i][0]:.6f}, {libration_points[i][1]:.6f}, {libration_points[i][2]:.6f}]")
print(f"Jacobi constant at initial state: {jacobi_constants[0]:.6f}"), print(f"")
print(f"-----------------------------------------------------------------"), print(f"")

# Plotting

# Jacobi constant conservation
cr3bp_orbit(trajectory)
plt.show()