# astrotools

A Python toolkit for astrodynamics — orbital propagation, coordinate transformations, and circular restricted three-body problem (CR3BP) analysis.

> **Status:** Work in progress. APIs and module layout may change.

---

## Features

- **Two-body propagation (R2BP)**, with optional J2 perturbation
- **CR3BP dynamics** — equations of motion, Jacobi constant, Lagrange/libration points
- **Numerical integrators** — 4th-order Yoshida symplectic, Dormand-Prince 8(5,3) (DOP853)
- **Differential correction** — state transition matrix (STM) propagation and planar single-shooting
- **Epoch & frame utilities** — Julian date, Greenwich sidereal time, ECI → ECEF transforms
- **Plotting** — ground tracks and CR3BP rotating-frame trajectories

---

## Installation

```bash
git clone <repo-url>
cd astrotools
pip install -e .
```

Dependencies: `numpy`, `scipy`, `matplotlib`

---

## Project structure

```
astrotools/
├── constants.py             # Physical/gravitational constants (Earth, Moon)
├── epoch.py                 # Julian date, GST, ECI->ECEF
├── plotting.py              # Ground track & CR3BP orbit plots
├── points.py                # Lagrange/libration point solver
├── stability.py             # STM propagation, single-shooting correction
├── dynamics/
│   ├── twobody.py           # 2BP acceleration, orbital elements <-> state vector
│   ├── j2.py                # 2BP + J2 perturbation
│   └── cr3bp.py             # CR3BP acceleration, Jacobi constant
├── integrators/
│   ├── yoshida4.py           # Symplectic integrator
│   └── dopr853.py            # DOP853 (via scipy)
├── trajectory/
│   ├── trajectory.py         # R2BP propagation wrapper
│   └── cr3bp_trajectory.py   # CR3BP propagation wrapper
└── examples/
    ├── R2BP_leo_sso.py               # LEO sun-synchronous orbit + ground track
    ├── R2BP_molniya.py               # Molniya orbit + ground track
    ├── R2BP_qzss_tundra.py           # QZSS Tundra orbit + ground track
    ├── CR3BP_L4.py                   # Perturbed L4 trajectory
    ├── CR3BP_2to1_resonant.py        # 2:1 resonant orbit
    ├── CR3BP_3to1_resonant.py        # 3:1 resonant orbit
    ├── CR3BP_L1_lyapunov.py          # L1 Lyapunov orbit
    └── CR3BP_shooting_L1_lyapunov.py # Single-shooting correction for L1 Lyapunov orbit
```

---

## Quick example

```python
import numpy as np
from astrotools.dynamics.twobody import oe_to_rv
from astrotools.trajectory.trajectory import trajectory
from astrotools.constants import mu_Earth

# Orbital elements -> state vector
a, e, i = 7000.0, 0.001, np.radians(51.6)
raan, argp, theta = 0.0, 0.0, 0.0
r0, v0 = oe_to_rv(a, e, i, raan, argp, theta, mu_Earth)

# Propagate one orbit
T = 2 * np.pi * np.sqrt(a**3 / mu_Earth)
n = 1000
traj = trajectory(r0, v0, n, T / n, perturbation=True, solverType="yoshida4")
```

See `astrotools/examples/` for full runnable scripts, including ground-track and CR3BP plotting.

---

## Roadmap

- [ ] Differential corrector for non-planar (3D) orbits
- [ ] Patched-conic / transfer trajectory tools

---

## References
[1] Koon, W.S., Lo, M.W., Marsden, J.E., Ross, S.D. (2011). *Dynamical Systems, the Three-Body Problem and Space Mission Design*. https://www.cds.caltech.edu/~marsden/volume/missiondesign/KoLoMaRo_DMissionBook_2011-04-25.pdf

[2] Curtis, H.D. (2021). *Orbital Mechanics for Engineering Students*. 4th ed. Butterworth-Heinemann.

[3] Patel, M., Shimane, Y., Lee, H.W., Ho, K. (2023). *Cislunar Satellite Constellation Design Via Integer Linear Programming*. https://arxiv.org/abs/2311.10252

[4] Yoshida, H. (1990). "Construction of higher order symplectic integrators." *Physics Letters A*, 150(5–7), 262–268.
