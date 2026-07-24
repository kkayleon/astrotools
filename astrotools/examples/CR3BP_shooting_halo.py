# Example: Single Shooting method for halo orbit

from astrotools.stability import single_shoot_halo, propagate_stm
from astrotools.plotting import cr3bp_orbit_3view
import numpy as np
import matplotlib.pyplot as plt

# Initial conditions (Table 7 "Target", https://arxiv.org/pdf/2605.07529)
r0 = np.array([0.85, 0.0, 0.173890])
v0 = np.array([0.0, 0.262114, 0.0])
state0 = np.concatenate([r0, v0])

# Single shoot method call
state0_corr, statef, Phi_T, k, tf = single_shoot_halo(state0, -1, 7)
if state0_corr is None:
    print("Newton's method did not converge for the initial condition")
else:
    print(f"Converged in {k} iterations")
    print(f"Corrected initial state: {state0_corr}")

    # Define step/step sizes
    n = 2000
    dt = 2*tf/n

    # Full trajectory propagation
    traj, __ = propagate_stm(state0_corr[:3], state0_corr[3:], n, dt)

    # 3-view plot (Match w/ Fig. 9 "Target" trajectory)
    cr3bp_orbit_3view(traj, bounds_xy=[0.7, 1.3, -0.3, 0.3], bounds_z=[-0.3, 0.3])
    plt.show()