# Case 009 — Variable PbI2-He diffusivity on fixed carrier flow

This case continues from case 008 with the same frozen helium carrier flow,
initial species field, inlet mass fraction, time step and numerical settings.

The main change is the diffusion model.

Instead of using the constant test coefficient from case 008, `rhoFixedFlowFoam`
computes the local molecular diffusivity of gaseous PbI2 in helium from the
frozen temperature and pressure fields.

The scalar equation is

    ddt(rho,Y_PbI2_g) + div(phi,Y_PbI2_g)
      - laplacian(rhoD_PbI2_g,Y_PbI2_g) = 0

where

    rhoD_PbI2_g = rho * D_PbI2_g

and `D_PbI2_g` is the molecular diffusivity in m2/s.

The diffusivity is calculated once when the solver starts, using the fixed
carrier-flow fields, and remains frozen during the transient species calculation.

The following fields are written for post-processing:

- `Y_PbI2_g` — gaseous PbI2 mass fraction
- `c_PbI2_g` — gaseous PbI2 molar concentration in mol/m3
- `D_PbI2_g` — molecular diffusivity in m2/s
- `rhoD_PbI2_g` — density-weighted diffusion coefficient in kg/(m s)

The pressure convention used by the diffusivity correlation is specified in
`constant/speciesTransportProperties`.

PbI2 remains a passive gaseous species. Chemistry, solid accumulation and
deposition are not included yet.

This case therefore marks the transition from:

    transient transport with a constant test diffusion coefficient

to

    transient transport with temperature- and pressure-dependent PbI2-He diffusivity.
