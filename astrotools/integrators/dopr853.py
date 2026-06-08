# Dormand-Prince 8th-order (from scipy.integrate)

from scipy.integrate import solve_ivp
import numpy as np
cos, sin, sqrt = np.cos, np.sin, np.sqrt

# Propagate trajectory 
def propagate(r0, v0, n, dt, func_accel):
    # dydt = f(t,y) input for DOPR853 solver
    def f(t,y):
        r = y[:3]
        v = y[3:]
        a = func_accel(r)
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