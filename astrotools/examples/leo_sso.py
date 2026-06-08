# Example: LEO SSO 

from astrotools.dynamics.twobody import oe_to_rv
from astrotools.trajectory.trajectory import trajectory
from astrotools.constants import mu_Earth, r_Earth
from astrotools.epoch import UTCtoJ0, GST, JD, latLongECEF
from astrotools.plotting import groundTrack
import numpy as np
import matplotlib.pyplot as plt
cos, sin, pi, sqrt = np.cos, np.sin, np.pi, np.sqrt

# Orbit parameters
zp = 510
e = 0.0 
a = (zp + r_Earth)/(1 - e)
print(f"{a}")
i = np.radians(97.45)
raan = np.radians(0.0)
argp = np.radians(0.0)
theta = np.radians(0.0)

# Epoch initialization
date = np.array([3, 1, 2027])
timeUTC = np.array([12, 0, 0])

# Setup trajectory propagation for a single orbit w/ 10000 steps
n = 10000
T = (2*pi/sqrt(mu_Earth)*a**1.5)
dt = T/n

# Initial state vector
r0, v0 = oe_to_rv(a, e, i, raan, argp, theta, mu_Earth)

# Propagate trajectory
traj = trajectory(r0, v0, n, dt, perturbation=False)

# Sidereal time parameters
J0 = UTCtoJ0(date)
julian_date = JD(timeUTC, J0)
theta_G = GST(J0, timeUTC)

# Latitude/Longitude
lat, long = latLongECEF(traj, theta_G)
startingLatitude = np.degrees(np.arcsin(r0[2]/np.linalg.norm(r0)))                      # Calculation based
startingLongitude = np.degrees(np.arctan2(r0[1], r0[0])) - np.rad2deg(theta_G) % 360    # Calculation based
finalLatitude = lat[-1]                                                                 # Final element based
finalLongitude = long[-1]                                                               # Final element based
if startingLongitude < 0: startingLongitude += 360
elif startingLongitude > 180: startingLongitude -= 360
else: pass
if finalLongitude < 0: finalLongitude += 360
elif finalLongitude > 180: finalLongitude -= 360
else: pass

# Output
print(f"Semi-major axis (a):                 {a:.3f} km")
print(f"Orbital period (T):                  {T:.3f} seconds")
print(f"Initial epoch:                       {date[0]:02d}/{date[1]:02d}/{date[2]:02d} {timeUTC[0]:02d}:{timeUTC[1]:02d}:{timeUTC[2]:02d} UTC")
print(f"Julian date at epoch:                {julian_date}")
print(f"Greenwich sidereal time (epoch):     {np.rad2deg(theta_G):.3f} deg")
print(f"Starting latitude:                   {startingLatitude:.3f} deg")
print(f"Starting longitude:                  {startingLongitude:.3f} deg")
print(f"Final latitude:                      {finalLatitude:.3f} deg")
print(f"Final longitude:                     {finalLongitude:.3f} deg")
groundTrack(traj,theta_G)
plt.show()

