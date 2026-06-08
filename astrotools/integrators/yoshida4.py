# 4th-order Yoshida symplectic integrator

import numpy as np
cos, sin, pi, sqrt = np.cos, np.sin, np.pi, np.sqrt

# Yoshida 4th-order symplectic integrator coefficients
beta = 2.0**(1/3)
c = [1/(2-beta)/2, (1-beta)/(2-beta)/2, (1-beta)/(2-beta)/2, 1/(2-beta)/2]
d = [1/(2-beta), -beta/(2-beta), 1/(2-beta), 0]

# 4th-order Yoshida step
def yoshida4(r, v, dt, func_accel): 
    for i in range(4):
        r = r + c[i]*dt*v
        a = func_accel(r)
        v = v + d[i]*dt*a
    return r, v

# Propagate trajectory (state vector of position, velocity, and time) 
def propagate(r0, v0, n, dt, func_accel):
    r, v = r0.copy(), v0.copy()
    trajectory = np.zeros((n+1, 7))
    trajectory[0] = np.concatenate([r, v, [0.0]])
    for i in range(1, n+1):
        r, v = yoshida4(r, v, dt, func_accel)
        trajectory[i] = np.concatenate([r, v, [i*dt]])

    return trajectory

    