# Trajectory propagation for (circular restricted three-body)

# No yoshida4 implementation yet -> using DOPR853

from astrotools.dynamics.cr3bp import acceleration
import numpy as np
from scipy.integrate import solve_ivp

# CR3BP Acceleration
func_accel_CR3BP = lambda r, v: acceleration(r, v)

# Integrator wrapper for CR3BP trajectory propagation
def propagate_dopr853_cr3bp(r0, v0, n, dt, func_accel):
    # dydt = f(t,y) input for DOPR853 solver
    def f(t,y):
        r = y[:3]
        v = y[3:]
        a = func_accel(r,v)
        # Input r, v => Output dr/dt, dv/dt (v, a)
        return np.concatenate([v, a])
    
    # Integration setup for solve_ivp
    tspan = (0, n*dt)
    y0 = np.concatenate([r0, v0])
    t_eval = np.linspace(0, n*dt, n+1)
    rtol = 1e-10
    atol = 1e-12

    # DOPR853 solve
    sol = solve_ivp(f, tspan, y0, 'DOP853', t_eval, rtol=rtol, atol=atol)

    # Construct state vector solutions [r, v, t]
    r = sol.y[:3, :].T           # From row vector -> column vector
    v = sol.y[3:, :].T           # "                              "
    t = sol.t.reshape(-1, 1)     # 1D array to "unknown" rows by 1 column (column vector shape)

    return np.hstack([r, v, t])


# Trajectory propagation for CR3BP
def trajectory(r0, v0, n, dt):
    return propagate_dopr853_cr3bp(r0, v0, n, dt, func_accel_CR3BP)