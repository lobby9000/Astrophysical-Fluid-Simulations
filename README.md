# Astrophysical fluid simulations

Numerical hydrodynamics coursework from UNAM: finite-difference solvers for advection and the 1D Euler equations (Lax and MacCormack schemes), building up to the set-up of a 2D astrophysical jet with radiative cooling. Written coursework from a continuum mechanics course is archived alongside.

**[Jet set-up](jet-simulation/jet_setup.py)** · **[Final project report (PDF, Spanish)](jet-simulation/final-project-report.pdf)**

## Highlights

- **2D jet set-up:** 400 × 150 cells over a 4,000 × 1,500 AU domain, 100-year runs with an output every year, MacCormack scheme with artificial viscosity and CFL = 0.1.
- **Thermal physics:** adiabatic index 5/3 and neutral or ionised mean molecular weight (1.3 or 0.61) switched at 10⁴ K. The report compares adiabatic and radiatively cooled runs.
- **Solver progression:** Lax and MacCormack schemes for linear advection, then 1D Euler solvers, before the 2D problem.

## Contents

| Folder | What it holds |
|---|---|
| `jet-simulation/` | `jet_setup.py` (constants, grid, initial and boundary conditions, jet injection region), plotting scripts and the final report |
| `solvers/advection/` | Lax and MacCormack advection solvers and their variants (`P*`, `PB*`, `Macormack*`) |
| `solvers/euler-1d/` | 1D Euler-equation solvers (`Euler*1D`, `EulerMac*`) |
| `reports/computational-fluid-dynamics/` | Homework for the computational fluid dynamics course (Spanish) |
| `reports/continuum-mechanics/` | Homework and exams for the continuum mechanics course (Spanish) |

`jet_setup.py` defines the problem (parameters, initial and boundary conditions) for the course's 2D hydrodynamics code, which is not included here. The plotting scripts read that code's output files.

## Run a solver

Each solver writes its snapshots to a folder named in its `output()` function, so create that folder first:

```bash
pip install -r requirements.txt
cd solvers/advection
mkdir -p mac_1 && python3 macormack.py
```

## Context and credits

- *Dinámica de Fluidos Computacional en Astrofísica* (2024-1), taught by Dr Juan Claudio Toledo Roy. Files whose header credits J. C. Toledo-Roy (`Lax_advection.py`, `ejercicio1.py`, `P1.py`) are course templates. The remaining scripts are my own solutions, several of them built on those templates.
- *Dinámica de Medios Deformables* (2023-2), taught by Dr Alejandro Sánchez Cervantes: written coursework in `reports/continuum-mechanics/`.

## Licence

Personal archive shared for portfolio purposes. Course templates remain the property of their authors.

---

Part of the portfolio of Luis Ángel Quiroz López · [lobby9000.github.io](https://lobby9000.github.io) · [LinkedIn](https://www.linkedin.com/in/llqq) · [GitHub](https://github.com/lobby9000)
