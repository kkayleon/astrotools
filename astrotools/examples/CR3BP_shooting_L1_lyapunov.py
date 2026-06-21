# Example: Single Shooting method for L1 Lyapunov orbit

from astrotools.stability import single_shoot_planar
from astrotools.dynamics.cr3bp import pi1, pi2
from astrotools.points import l_points
import numpy as np
import matplotlib.pyplot as plt

# Initial conditions (slightly adjusted) from astrotools/examples/CR3BP_L1_lyapunov.py
r0 = np.array([0.80, 0.0, 0.0])
v0 = np.array([0.0, 0.33, 0.0])
T = 3.225
n = 1000
dt = T/n

# Half period guess
n_half = n//2

# Initialize state and call single shot method w/ half period guess
tol = 1E-10
max_iterations = 20

state0 = np.concatenate([r0, v0])
state0_corr, traj, Phi_T, k = single_shoot_planar(state0, dt, n_half, tol, max_iterations)
print("Converged within", tol, "in", k, "iterations")

# Construct full trajectory via mirroring
traj_mirror = traj[:, :6].copy()*np.array([1, -1, -1, -1, 1, 1])
traj_mirror = traj_mirror[::-1]
trajectory = np.concatenate([traj[:,:6], traj_mirror])

# Plot results
libration_points = l_points()
fig, ax = plt.subplots(figsize=(14, 7))
names = ['L1', 'L2', 'L3', 'L4', 'L5']
colors = ['red', 'red', 'red', 'green', 'green']
for i, (name, color) in enumerate(zip(names, colors)):
    ax.scatter(libration_points[i][0], libration_points[i][1], color=color, s=50, zorder=5)
    ax.annotate(name, (libration_points[i][0], libration_points[i][1]), textcoords="offset points", xytext=(5,5))
ax.scatter(-pi2, 0, color='blue', s=100, zorder=5, label='Earth')
ax.scatter(pi1, 0, color='gray', s=30, zorder=5, label='Moon')
ax.plot(trajectory[:,0], trajectory[:,1], label='Trajectory')
# ax.plot(traj[:,0], traj[:,1], label='Half trajectory')
ax.legend()
ax.set_xlabel('x (non-dim)')
ax.set_ylabel('y (non-dim)')
ax.set_title('L1 Lyapunov Trajectory from Single-shoot in Rotating Frame')
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.2, 1.2)
ax.grid(True)
plt.show()
