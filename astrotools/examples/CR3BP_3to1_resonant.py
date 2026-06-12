# Example: Circular restricted 3-body dynamics for a 3:1 resonant orbit
# Reference: arxiv.org/pdf/2311.10252
# 3:1 Resonant orbit initial conditions from Appendix B, Table B3

from astrotools.dynamics.cr3bp import jacobi_constant, pi1, pi2
from astrotools.trajectory import cr3bp_trajectory
from astrotools.points import l_points
from astrotools.plotting import cr3bp_orbit
import matplotlib.pyplot as plt
import numpy as np
cos, sin, sqrt = np.cos, np.sin, np.sqrt

# Initial conditions
r0 = np.array([0.13603399956670137, 0.0, 0.0])
v0 = np.array([0.0, 3.202418276067991, 0.0])

# Propagate trajectory
T = 6.45        # Nondimensional period
n = 1000
dt = T/n
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
print(f"Jacobi constant: {jacobi_constants[0]:.6f}"), print(f"")
print(f"-----------------------------------------------------------------"), print(f"")

# Plotting trajectory in rotating frame
cr3bp_orbit(trajectory)
plt.show()