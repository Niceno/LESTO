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
There is no SIMPLE or PIMPLE loop. Only the scalar field `Y_PbI2_g` is advanced
in physical time.

The scalar equation is

    ddt(rho,Y_PbI2_g) + div(phi,Y_PbI2_g) - laplacian(D,Y_PbI2_g) = 0

where `phi` is the frozen mass flux from case 007.

`Y_PbI2_g` represents the dimensionless gaseous PbI2 mass fraction (kg PbI2 / kg
mixture). Its inlet value of `1e-6` means 1 ppm by mass. This rename does
not change any values or convert them to mol/m3. The present diffusion
coefficient is still a test value; physical species diffusivities and
deposition models will be introduced in later cases.

This case therefore marks the transition from:

    steady carrier-flow calculation

to

    transient species transport on a frozen carrier flow.

## Molar concentration output

At scheduled write times, `c_PbI2_g = rho * Y_PbI2_g / MPbI2` is written
for ParaView in mol/m3. `MPbI2 = 0.46100894 kg/mol`, using
207.2 + 2*126.90447 g/mol. No `0/c_PbI2_g` input is needed.
The inlet remains 1e-6 (1 ppm by mass of PbI2). D remains the test
coefficient; no solid PbI2 inventory or deposition is implemented yet.
