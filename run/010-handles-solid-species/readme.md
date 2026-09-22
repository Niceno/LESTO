# Case 010 — Solid-species accumulation

This case continues the multi-species structure of case 009 and introduces the
first solid species.

The frozen helium carrier flow and the two gaseous species remain unchanged:

- `PbI2_g` uses the temperature- and pressure-dependent PbI2-He molecular
  diffusivity;
- `tracer` uses the constant molecular diffusivity `D = 1e-5 m2/s`.

A third species is introduced:

- `PbI2_s` has `state solid` and represents accumulated solid PbI2.

Gaseous species continue to solve

    ddt(rho,Y_i) + div(phi,Y_i) - laplacian(rhoD_i,Y_i) = 0

with

    rhoD_i = rho * D_i.

The solid species has no convection and no diffusion. It is advanced only by

    ddt(rho,Y_PbI2_s) = S

using a prescribed uniform source

    S = 1e-6 kg/(m3 s).

Because the equation contains only a transient term and a local source, cells
are uncoupled and the equation is solved directly with OpenFOAM's `diagonal`
solver.

Starting from `Y_PbI2_s = 0` and with frozen density, the expected internal
solution is

    Y_PbI2_s = S * t / rho.

The solid mass fraction therefore grows linearly in time. Its spatial variation
comes only from the frozen density field.

The gaseous-species equations are unchanged from case 009 and their residual
histories are expected to remain unchanged.

The solid source is currently applied in every cell. It is not yet restricted
to near-wall cells, and no chemistry or gas-to-solid coupling is included.
