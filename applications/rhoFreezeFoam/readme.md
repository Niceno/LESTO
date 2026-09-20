# rhoFreezeFoam

## Purpose

Transient transport of Y_PbI2_g [1] on a frozen converged rhoSimpleFoam carrier flow.

## Frozen fields

U, p, e/T, rho and phi are read from the steady solution and are never solved
or updated. No SIMPLE/PIMPLE loop is used.

## Scalar equation

Y_PbI2_g is advanced in physical time using OpenFOAM's scalarTransport class:

    ddt(rho,Y_PbI2_g) + div(phi,Y_PbI2_g) - laplacian(D,Y_PbI2_g) = 0

For mass-flux phi, D has units kg/(m s), i.e. D = rho * D_molecular.

## Molar concentration output

At each scheduled write time, the solver also writes `c_PbI2_g` in mol/m3:

    c_PbI2_g = rho * Y_PbI2_g / MPbI2

`MPbI2` is 0.46100894 kg/mol for gaseous PbI2. The conversion uses local
density and the current mass fraction, including boundary values. The field
has dimensions `[0 -3 0 0 1 0 0]`; its numerical values use mol, not kmol.
No additional transport equation is solved. No `0/c_PbI2_g` input is required;
the field is recalculated on output, including after a restart. Select
`c_PbI2_g` in ParaView to display molar concentration. Existing time directories
are not retroactively updated.

## Current limitations

- Y_PbI2_g is the dimensionless gaseous PbI2 mass fraction: kg PbI2 / kg mixture.
- Y_PbI2_g remains the transported variable; c_PbI2_g is derived for output only.
- No deposition or chemistry.
- D is presently the old constant test coefficient.
- Carrier flow must already be converged.

## Run

Compile rhoFreezeFoam with wmake.

Prepare a case containing the frozen carrier fields and Y_PbI2_g, then run:

    rhoFreezeFoam

Before physical LESTO calculations, replace the test diffusion coefficient
with the appropriate species diffusivity D(T,p).
