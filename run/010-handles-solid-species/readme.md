# Case 010 — Solid-species accumulation with mock thermochemistry

This case continues the multi-species structure of case 009 and introduces the
first solid species together with a single thermochemistry interface intended
to mimic a future GEMS call.

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

    ddt(rho,Y_PbI2_s) = S_PbI2_s.

The source terms are evaluated by `evaluateThermochemistry()` once in every
cell and physical time step, after the gaseous species have been advanced. The
routine receives the local `T`, absolute `p`, the complete species-name list
and the current local species values, and returns one source for every species.

At present this routine is only a mock-up of a future GEMS call. It returns zero
for all species except `PbI2_s`, for which

    S_PbI2_s = S0 * (Thot - T)/(Thot - Tcold)

with local test constants

    S0    = 1e-6 kg/(m3 s)
    Thot  = 1000 K
    Tcold = 300 K.

Pressure and species values are passed to the routine but are not used yet.
They are already part of the interface so that future thermochemistry can
couple all species without restructuring the CFD solver.

Because the solid equation contains only a transient term and a local source,
cells are uncoupled and the equation is solved directly with OpenFOAM's
`diagonal` solver.

The gaseous-species equations remain unchanged from case 009 and their residual
histories are expected to remain unchanged.

The thermochemistry routine still evaluates a candidate source in every cell.
The CFD solver then retains the solid source only in the first layer of cells
adjacent to boundary patches of OpenFOAM type `wall` and sets it to zero
everywhere else. The wall-cell mask is constructed once at startup directly
from the boundary-face-to-cell connectivity because the mesh is fixed.

The solid equation itself is unchanged: it remains a local transient
accumulation equation solved with the `diagonal` solver. No gas-to-solid mass
coupling is included yet, so formation of `PbI2_s` does not yet remove mass
from `PbI2_g`.
