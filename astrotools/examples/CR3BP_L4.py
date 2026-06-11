# Example: Circular restricted 3-body dynamics around Earth-Moon L4

from astrotools.dynamics.cr3bp import jacobi_constant, pi1, pi2
from astrotools.trajectory import cr3bp_trajectory
import matplotlib.pyplot as plt
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


# Plots
# Trajectory in rotating frame
fig, ax = plt.subplots(figsize=(14, 7))
ax.plot(trajectory[:,0], trajectory[:,1])
ax.set_xlabel('x (non-dim)')
ax.set_ylabel('y (non-dim)')
ax.set_title('Trajectory in Rotating Frame')
ax.grid(True)

# Jacobi constant conservation
fig2, ax2 = plt.subplots(figsize=(14, 7))
jacobi_constant_diff = jacobi_constants - jacobi_constants[0]
t = trajectory[:,6]
ax2.plot(t, jacobi_constant_diff)
ax2.set_xlabel('Time (non-dim)')
ax2.set_ylabel('C')
ax2.set_title('Jacobi Constant Conservation')
ax2.grid(True)

plt.tight_layout()
plt.show()