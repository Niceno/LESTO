# rhoFreezeFoam

## Purpose

Transient transport of Pb_g on a frozen converged rhoSimpleFoam carrier flow.

## Frozen fields

U, p, e/T, rho and phi are read from the steady solution and are never solved
or updated. No SIMPLE/PIMPLE loop is used.

## Scalar equation

Pb_g is advanced in physical time using OpenFOAM's scalarTransport class:

    ddt(rho,Pb_g) + div(phi,Pb_g) - laplacian(D,Pb_g) = 0

For mass-flux phi, D has units kg/(m s), i.e. D = rho * D_molecular.

## Current limitations

- Pb_g is still a dimensionless test scalar.
- No deposition or chemistry.
- D is presently the old constant test coefficient.
- Carrier flow must already be converged.

## Run

Compile rhoFreezeFoam with wmake.

Prepare a case containing the frozen carrier fields and Pb_g, then run:

    rhoFreezeFoam

Before physical LESTO calculations, replace the test diffusion coefficient
with the appropriate species diffusivity D(T,p).
