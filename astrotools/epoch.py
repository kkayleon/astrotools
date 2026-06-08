# Epoch initialization & coordinate transformation functions

import numpy as np
cos, sin, pi = np.cos, np.sin, np.pi

# Date & Time in UTC to Julian day number at 0h UTC
def UTCtoJ0(date):
    return 367*date[2] - int(7/4*(date[2]+int((date[0] + 9)/12))) + int(275*date[0]/9) + date[1] + 1721013.5

# Time (UTC) to Julian day
def JD(timeUTC, J0):
    return J0 + (timeUTC[0] + timeUTC[1]/60 + timeUTC[2]/3600)/24

# Grennwich sidereal time
def GST(J0, timeUTC):
    T0 = (J0 - 2451545)/36525
    theta_G0 = (100.4606184 + 36000.77004*T0 + 0.000387933*T0**2 - 2.583E-8 *T0**3) % 360
    theta_G = (theta_G0 + 360.985864724 * (timeUTC[0] + timeUTC[1]/60 + timeUTC[2]/3600)/24) % 360
    return np.deg2rad(theta_G) 

# Coordinate transformation from Earth-Centered Inertial RF to Earth-Centered Earth-Fixed RF
def latLongECEF(traj, theta_G):
    t = traj[:,6]
    r_ECI = np.array([traj[:,0], traj[:,1], traj[:,2]])

    Omega_Earth = 7.292115900231276e-05     # Rotational rate of Earth [rad/s]
    theta_G = theta_G + Omega_Earth*t       # Greenwich Sidereal Time [rad]

    # Multplied out from 1-tensor transformation law
    x_ECEF = r_ECI[0]*cos(theta_G) + r_ECI[1]*sin(theta_G) + r_ECI[2]*0.0
    y_ECEF = r_ECI[0]*-sin(theta_G) + r_ECI[1]*cos(theta_G) + r_ECI[2]*0.0
    z_ECEF = r_ECI[2]

    # From Cartesian = spherical coords -> Solving for latitude and longitude angles
    r_mag = np.sqrt(x_ECEF**2 + y_ECEF**2 + z_ECEF**2)
    lat = np.degrees(np.arcsin(z_ECEF/r_mag))
    long = np.degrees(np.arctan2(y_ECEF, x_ECEF))

    return lat, long
