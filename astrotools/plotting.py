# Plotting

from astrotools.constants import r_Earth
from astrotools.epoch import latLongECEF
import numpy as np
import matplotlib.pyplot as plt
cos, sin, sqrt = np.cos, np.sin, np.sqrt

# Ground track of satellite (Earth-Centered Earth-Fixed RF)
def groundTrack(traj, theta_G):
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