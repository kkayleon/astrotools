# Example: QZSS Tundra orbit

from astrotools.dynamics.twobody import oe_to_rv
from astrotools.trajectory.trajectory import trajectory
from astrotools.constants import mu_Earth
from astrotools.epoch import UTCtoJ0, GST, latLongECEF
import numpy as np
cos, sin, pi, sqrt = np.cos, np.sin, np.pi, np.sqrt

# Orbit parameters
a = 42164.0
e = 0.075 
i = np.radians(43.0)
raan = np.radians(195.0)
argp = np.radians(270.0)
theta = np.radians(305.0)

# Epoch initialization
date = [1, 1, 2024]
timeUTC = [0, 0, 0]

# Setup trajectory propagation for a single orbit w/ 10000 steps
n = 10000
dt = 2*pi/sqrt(mu_Earth)*a**1.5/n

# Initial state vector
r0, v0 = oe_to_rv(a, e, i, raan, argp, theta, mu_Earth)

# Propagate trajectory
traj = trajectory(r0, v0, n, dt, perturbation=False)