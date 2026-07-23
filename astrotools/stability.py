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
    return np.block([[zero, eye3], [U, coriolis]])

# Propagation of STM & final STM
def propagate_stm(r0, v0, n, dt):
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

# Single-shooting method: Planar
def single_shoot_planar(state0, dt, n, tol=1E-10, max_iteration=20):
    # Newton's method loop
    for k in range(max_iteration):
        traj, Phi_T = propagate_stm(state0[:3], state0[3:], n, dt)  
        statef = traj[-1,:6]

        # Residual conditions to be driven within tolerable error
        F = statef - state0
        F_red = F[3]
        error = np.linalg.norm(F_red)  
        if error < tol:
            return state0, traj, Phi_T, k
        
        # Residual calculations to drive next iteration towards tolerable error
        J_red = Phi_T[3,4]
        delta = -F_red/J_red       
        state0[4] += delta

    print("Did not converge within", tol, "in", max_iteration, "iterations")
    # None if Newton's method does not converge
    return None, None, None, None

# Event-aware propagation of STM
def propagate_stm_event(state0, direction, t_max):
    # dydt for state and STM
    def f(t,y):
        # Initial state & STM
        state = y[:6]                   # State
        r, v = state[:3], state[3:]
        Phi = y[6:].reshape(6,6)        # State transition matrix
        # State & STM derivatives
        state_dot = np.concatenate([v, acceleration(r,v)])
        Phi_dot = jacobian(r) @ Phi
        return np.concatenate([state_dot, Phi_dot.flatten()])

    # Event function (xz-plane crossing)
    def xz_plane_event(t,y):
        return y[1]
    xz_plane_event.terminal = True          # Stops propagation at event
    xz_plane_event.direction = direction    # +1 (upward) or -1 (downward)

    # Initial condition
    y0 = np.concatenate([state0, np.eye(6).flatten()])
    
    # Integration setup for solve_ivp
    rtol = 1e-10
    atol = 1e-12

    # DOP853 solve
    sol = solve_ivp(f, (0, t_max), y0, 'DOP853', events=xz_plane_event, rtol=rtol, atol=atol)

    # Extract final states
    if sol.t_events[0].size == 0:
        raise RuntimeError(f"No xz-plane crossing found before t_max = {t_max}")

    event_array = sol.y_events[0][0]
    statef = event_array[:6]
    tf = sol.t_events[0][0]
    Phi_T = event_array[6:].reshape(6, 6)

    return statef, Phi_T, tf

# Timing correction
def stm_event_correction(statef, Phi_T):
    rf = statef[:3]
    vf = statef[3:]
    af = acceleration(rf, vf)

    statef_dot = np.concatenate([vf, af])
    Phi_corrected = Phi_T - np.outer(statef_dot, Phi_T[1, :]) / vf[1]

    return Phi_corrected

# Single-shooting method: Halo
def single_shoot_halo(state0, direction, t_max, tol=1E-10, max_iteration=20):
    # Newton's method loop
    for k in range(max_iteration):

        r0 = state0[:3]
        v0 = state0[3:]
        statef, Phi_T, tf = propagate_stm_event(state0, direction, t_max)

        # Residual conditions to be driven within tolerable error
        F_red = np.array([statef[3], statef[5]])
        error = np.linalg.norm(F_red)  
        if error < tol:
            return state0, statef, Phi_T, k, tf
        
        Phi_corrected = stm_event_correction(statef, Phi_T)

        # Residual calculations to drive next iteration towards tolerable error
        J_red = Phi_corrected[np.ix_([3,5],[0,4])]
        delta = -np.linalg.solve(J_red, F_red)
        state0[0] += delta[0]
        state0[4] += delta[1]

    print("Did not converge within", tol, "in", max_iteration, "iterations")
    return None, None, None, None, None