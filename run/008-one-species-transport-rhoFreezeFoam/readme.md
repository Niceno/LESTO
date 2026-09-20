# Case 008 — Transient scalar transport on frozen carrier flow

This case continues from the converged steady helium solution of case 007,
but uses the new `rhoFreezeFoam` solver to solve transport of a scalar.

Case 007 provides the carrier-flow fields:

- velocity `U`
- pressure `p`
- temperature / internal energy
- density `rho`
- mass flux `phi`

These fields are read from the final steady solution and remain frozen
throughout the calculation.

`rhoFreezeFoam` does not solve momentum, pressure, continuity or energy.
There is no SIMPLE or PIMPLE loop. Only the scalar field `Pb_g` is advanced
in physical time.

The scalar equation is

    ddt(rho,Pb_g) + div(phi,Pb_g) - laplacian(D,Pb_g) = 0

where `phi` is the frozen mass flux from case 007.

`Pb_g` is currently a dimensionless test scalar. The present diffusion
coefficient is still a test value; physical species diffusivities and
deposition models will be introduced in later cases.

This case therefore marks the transition from:

    steady carrier-flow calculation

to

    transient species transport on a frozen carrier flow.
