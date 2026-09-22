# Case 009 — Variable PbI2-He diffusivity with multiple gaseous species

This case continues the multi-species structure introduced in the revised case
008. The frozen helium carrier flow, time step and numerical settings are
unchanged.

Two passive gaseous species are transported:

- `PbI2_g` uses the temperature- and pressure-dependent PbI2-He molecular
  diffusivity from the frozen carrier fields;
- `tracer` retains the constant molecular diffusivity `D = 1e-5 m2/s` used in
  the revised case 008.

Each species solves

    ddt(rho,Y_i) + div(phi,Y_i) - laplacian(rhoD_i,Y_i) = 0

with

    rhoD_i = rho * D_i.

For `PbI2_g`, `D_PbI2_g(T,p)` is evaluated once when the solver starts and
remains fixed during the transient calculation. The correlation pressure is in
atm; the present carrier pressure is 101325 Pa = 1 atm.

The PbI2 equation, diffusivity correlation, initial/boundary field, time step,
convection/diffusion schemes, GAMG settings, relaxation factor, `nCorr` and
outer tolerance are unchanged from the previous single-species case 009.
Therefore the `Y_PbI2_g` residual history is expected to reproduce the previous
case-009 residual history; the added tracer is uncoupled and is solved only
after PbI2 at each physical time step.

Written PbI2 fields include `Y_PbI2_g`, `c_PbI2_g`, `D_PbI2_g` and
`rhoD_PbI2_g`. Chemistry, precipitation and deposition are not included yet.
