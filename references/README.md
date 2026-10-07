# References

This directory contains publicly available primary literature used in the
independent development of the isothermal transient-pipe simulator.

Papers are stored locally only when a publicly accessible copy is available
from the authors, an institutional repository, arXiv, or another authoritative
public source.

The presence of a paper here does not imply that its implementation is being
copied. The simulator equations and algorithms are independently derived and
documented in `docs/model.tex` and `docs/provenance.md`.

## Gyrya & Zlotnik (2019)

Vitaliy Gyrya and Anatoly Zlotnik,

**An explicit staggered-grid method for numerical simulation of large-scale
natural gas pipeline networks**

Applied Mathematical Modelling, 65 (2019), 34-51.

DOI:

https://doi.org/10.1016/j.apm.2018.07.051

Public manuscript:

https://arxiv.org/abs/1803.00418

Local file:

`Gyrya_Zlotnik_2019_explicit_staggered_grid.pdf`

Relevant to this project:

- staggered finite-difference discretization
- pressure/density and mass-flux variables
- discrete mass conservation
- general equation-of-state formulation
- transient gas-pipeline simulation

## Brodskyi, Gyrya & Zlotnik (2025)

Yan Brodskyi, Vitaliy Gyrya, and Anatoly Zlotnik,

**Simulation of Gas Mixture Dynamics in a Pipeline Network using Explicit
Staggered-Grid Discretization**

Applied Mathematical Modelling, 137 (2025), 115717.

DOI:

https://doi.org/10.1016/j.apm.2024.115717

Public manuscript:

https://arxiv.org/abs/2404.04451

LANL accepted manuscript:

https://www.osti.gov/servlets/purl/2447960

Government report number:

LA-UR-24-22948

Local files:

`Brodskyi_Gyrya_Zlotnik_2025_gas_mixture_staggered_grid.pdf`

`Brodskyi_Gyrya_Zlotnik_LANL_LA-UR-24-22948.pdf`

Relevant to this project:

- modern staggered-grid formulation
- gas-mixture pipeline dynamics
- conservation laws
- equation-of-state treatment
- boundary compatibility conditions

The two local files are different public distributions of substantially the
same work. Keeping both preserves the source/provenance trail.

## Kiuchi (1994)

Tatsuhiko Kiuchi,

**An implicit method for transient gas flows in pipe networks**

International Journal of Heat and Fluid Flow, 15(5) (1994), 378-383.

DOI:

https://doi.org/10.1016/0142-727X(94)90051-5

Publisher page:

https://www.sciencedirect.com/science/article/pii/0142727X94900515

The publisher PDF is not included because it is not openly downloadable.

Relevant to this project:

- fully implicit finite-difference transient gas simulation
- Newton-Raphson nonlinear solution
- comparison with characteristic and Lax-Wendroff methods

## Finite-difference theta-scheme paper (2024)

**Gas composition tracking feasibility using transient finite difference
theta-scheme model for binary gas mixtures**

International Journal of Hydrogen Energy, 49 (2024), 1319-1331.

DOI:

https://doi.org/10.1016/j.ijhydene.2023.11.031

Publisher page:

https://www.sciencedirect.com/science/article/pii/S0360319923056793

The article is published open access under a Creative Commons license.

Relevant to this project:

- transient finite-difference theta scheme
- interpretation of theta = 0, 0.5, and 1
- stability behavior for theta >= 0.5
- spatial and temporal discretization studies

We retain the authoritative publisher link here. A local copy may be added
after verifying the publisher's direct PDF URL and redistribution terms.

## Provenance policy

Prior proprietary pipeline-simulation experience may identify questions worth
investigating, but it is not used as a source for equations, algorithms,
implementation details, constants, heuristics, or validation cases.

Implementation decisions should be traceable to:

1. governing conservation laws;
2. independent derivation;
3. publicly available literature; or
4. numerical experiments performed within this project.

See `../docs/provenance.md` for the project-wide provenance record.
