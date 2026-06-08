# Two-body dynamics

import numpy as np
cos, sin, pi, sqrt = np.cos, np.sin, np.pi, np.sqrt

# Gravitational acceleration (2BP)
def acceleration(r, mu):
    r_mag = np.linalg.norm(r)
    return -mu*r/r_mag**3

# Specific orbital energy and specific angular momentum
def specific_energy(r, v, mu):
    return 0.5*np.dot(v, v) - mu/np.linalg.norm(r)

def specific_angular_momentum_vec(r, v):
    return np.cross(r, v)

# Orbital period (2BP)
def orbital_period(a, mu):
    return 2*pi*sqrt(a**3/mu)

# Orbital elements to state vector (position and velocity)
def oe_to_rv(a, e, i, raan, argp, theta, mu):
    h = sqrt(a*mu*(1 - e**2))
    r = h**2/mu/(1 + e*cos(theta)) * np.array([cos(theta), sin(theta), 0])
    v = mu/h * np.array([-sin(theta), e + cos(theta), 0])

    # Rotation matrix from perifocal to inertial frame (ECI)
    R_PI = np.array([[-sin(raan)*cos(i)*sin(argp)+cos(raan)*cos(argp), cos(raan)*cos(i)*sin(argp)+sin(raan)*cos(argp), sin(i)*sin(argp)], 
                     [-sin(raan)*cos(i)*cos(argp)-cos(raan)*sin(argp), cos(raan)*cos(i)*cos(argp)-sin(raan)*sin(argp), sin(i)*cos(argp)], 
                     [ sin(raan)*sin(i),                              -cos(raan)*sin(i),                               cos(i)          ]])
    R_IP = R_PI.T

    # Tensor transformation law (rank-1)
    return R_IP @ r, R_IP @ v

# State vector to orbital elements
def rv_to_oe(r, v, mu):
    # Curtis Algorithm 4.2
    h_vec = specific_angular_momentum_vec(r, v)
    h = np.linalg.norm(h_vec)
    vr = np.dot(r,v)/np.linalg.norm(r)

    n_vec = np.cross([0, 0, 1], h_vec)
    n = np.linalg.norm(n_vec)
    e_vec = (np.cross(v, h_vec) - mu*r/np.linalg.norm(r))/mu
    e = np.linalg.norm(e_vec)

    # Circular/equatorial orbits case handling
    tol = 1e-8 
    if e < tol:
        raise ValueError("Circular orbit => argp is undefined")
    if n < tol:
        raise ValueError("Equatorial orbit => raan is undefined")
    
    # Orbital element calculations
    a = h**2/mu/(1 - e**2)
    i = np.arccos(h_vec[2]/h)
    raan = np.arccos(n_vec[0]/n) if n_vec[1] >= 0 else 2*pi - np.arccos(n_vec[0]/n)
    argp = np.arccos(np.dot(n_vec, e_vec)/n/e) if e_vec[2] >= 0 else 2*pi - np.arccos(np.dot(n_vec, e_vec)/n/e)
    theta = np.arccos(np.dot(e_vec, r)/e/np.linalg.norm(r)) if vr >= 0 else 2*pi - np.arccos(np.dot(e_vec, r)/e/np.linalg.norm(r))

    return a, e, i, raan, argp, theta
