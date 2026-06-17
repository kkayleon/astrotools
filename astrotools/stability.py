# Generation of orbits

from astrotools.dynamics.cr3bp import pi1, pi2, r1r2, acceleration
from scipy.integrate import solve_ivp
import numpy as np
cos, sin, sqrt = np.cos, np.sin, np.sqrt

# Hessian (gravity gradient matrix)
def hessian(r):
    # Second derivative of psuedo-potential
    r1, r2 = r1r2(r)
    Uxx = 1 - pi1/r1**3 - pi2/r2**3 + 3*pi1*(r[0]+pi2)**2/r1**5 + 3*pi2*(r[0]-pi1)**2/r2**5
    Uyy = 1 - pi1/r1**3 - pi2/r2**3 + 3*pi1*r[1]**2/r1**5 + 3*pi2*r[1]**2/r2**5
    Uzz =   - pi1/r1**3 - pi2/r2**3 + 3*pi1*r[2]**2/r1**5 + 3*pi2*r[2]**2/r2**5
    Uxy = 3*pi1*(r[0]+pi2)*r[1]/r1**5 + 3*pi2*(r[0]-pi1)*r[1]/r2**5
    Uxz = 3*pi1*(r[0]+pi2)*r[2]/r1**5 + 3*pi2*(r[0]-pi1)*r[2]/r2**5
    Uyz = 3*pi1*r[1]*r[2]/r1**5 + 3*pi2*r[1]*r[2]/r2**5
    return np.array([[Uxx, Uxy, Uxz], [Uxy, Uyy, Uyz], [Uxz, Uyz, Uzz]])

# Jacobian matrix
def jacobian(r):
    # Blocks of the matrix
    zero = np.zeros((3,3))
    eye3 = np.eye(3)
    U = hessian(r)
    coriolis = np.array([[0, 2, 0], [-2, 0, 0], [0, 0, 0]])
    return np.block([[zero, eye3], [-U, coriolis]])

# Propagation of STM & final STM
def propagate_stm(r0, v0, dt, n):
    # dydt for state & STM
    def f(t,y):
        # Initial state & STM
        state = y[:6]                   # State
        r, v = state[:3], state[3:]
        Phi = y[6:].reshape(6,6)        # State transition matrix
        # State & STM derivatives
        state_dot = np.concatenate([v, acceleration(r,v)])
        Phi_dot = jacobian(r) @ Phi
        return np.concatenate([state_dot, Phi_dot.flatten()])
    
    # Integration setup for solve_ivp
    tspan = (0, n*dt)
    Phi0 = np.eye(6).flatten()
    y0 = np.concatenate([r0, v0, Phi0])
    t_eval = np.linspace(0, n*dt, n+1)
    rtol = 1e-10
    atol = 1e-12

    # DOPR853 solve
    sol = solve_ivp(f, tspan, y0, 'DOP853', t_eval, rtol=rtol, atol=atol)

    # Construct state vector solutions [r, v, Phi, t]
    r = sol.y[:3, :].T           # From row vector -> column vector
    v = sol.y[3:6, :].T          # "                              "
    Phi = sol.y[6:].T            # State transition matrix 
    t = sol.t.reshape(-1, 1)     # 1D array to "unknown" rows by 1 column (column vector shape)

    # Return trajectory (42D) & STM at final time
    trajectory = np.hstack([r, v, Phi, t])
    Phi_T = sol.y[6:, -1].reshape(6, 6)

    return trajectory, Phi_T

