# Plotting

from astrotools.constants import r_Earth
from astrotools.epoch import latLongECEF
from astrotools.dynamics.cr3bp import pi1, pi2
from astrotools.points import l_points
import numpy as np
import matplotlib.pyplot as plt
cos, sin, sqrt = np.cos, np.sin, np.sqrt

# Ground track of satellite (Earth-Centered Earth-Fixed RF)
def groundTrack(traj, theta_G):
    # For R2BP only
    lat, long = latLongECEF(traj, theta_G)

    # Wrapping from 180deg -> -180deg
    long_plot = long.copy()
    wrap = np.abs(np.diff(long_plot)) > 180
    long_plot[1:][wrap] = np.nan

    fig, ax = plt.subplots(figsize=(14,7))

    # Mercator projection background (Reference [4])
    try: 
        img = plt.imread('earth.jpg')
        ax.imshow(img, extent=[-180, 180, -90, 90], aspect='auto', zorder=0)
    except:
        print(f"earth.jpg not found. Proceeding with plotting without background.")

    ax.plot(long_plot, lat, linewidth=1.0, c='yellow')

    ax.scatter(long_plot[0], lat[0], color='green', s=50, zorder=5, label='Start')
    ax.scatter(long_plot[-1], lat[-1], color='red', s=50, zorder=5, label='End')

    ax.set_xlabel('Longitude [deg]')
    ax.set_ylabel('Latitude [deg]')
    ax.set_title('Ground Track of Satellite')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.grid(True)
    ax.legend()


# Circular restricted 3 body trajectory plot in orbital plane (z=0)
libration_points = l_points()
def cr3bp_orbit(traj, bounds=[-1.5, 1.5, -1.2, 1.2]):
    fig, ax = plt.subplots(figsize=(14, 7))
    names = ['L1', 'L2', 'L3', 'L4', 'L5']
    colors = ['red', 'red', 'red', 'green', 'green']
    for i, (name, color) in enumerate(zip(names, colors)):
        ax.scatter(libration_points[i][0], libration_points[i][1], color=color, s=50, zorder=5)
        ax.annotate(name, (libration_points[i][0], libration_points[i][1]), textcoords="offset points", xytext=(5,5))
    ax.scatter(-pi2, 0, color='blue', s=100, zorder=5, label='Earth')
    ax.scatter(pi1, 0, color='gray', s=30, zorder=5, label='Moon')
    ax.plot(traj[:,0], traj[:,1], label='Trajectory')
    ax.legend()
    ax.set_xlabel('x (non-dim)')
    ax.set_ylabel('y (non-dim)')
    ax.set_title('3:1 Resonant Trajectory')
    ax.set_xlim(bounds[0], bounds[1])
    ax.set_ylim(bounds[2], bounds[3])
    ax.grid(True)