# rhoFixedFlowFoam

## Purpose

Transient transport of Y_PbI2_g [1] on a frozen converged rhoSimpleFoam carrier flow.

## Frozen fields

U, p, e/T, rho and phi are read from the steady solution and are never solved
or updated. No SIMPLE/PIMPLE loop is used.

## Scalar equation

Y_PbI2_g is advanced in physical time. The legacy constant-D model uses
OpenFOAM's scalarTransport class; the variable model assembles the same
density-weighted equation directly:

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
- Choose the legacy constant coefficient or the PbI2-He variable model below.
- The variable branch currently supports transport without fvOptions sources.
- Carrier flow must already be converged.

## Run

Compile rhoFixedFlowFoam with wmake.

Prepare a case containing the frozen carrier fields and Y_PbI2_g, then run:

    rhoFixedFlowFoam

## Diffusivity models

The `functions/Y_PbI2_gTransport` controls in `system/controlDict` remain
in use. Function objects are switched off and the solver drives transport
once per physical time step. Settings are read at startup.

For existing cases, no change is required: absent `diffusivityModel` means
`constant`. `D 1e-5;` retains its original units kg/(m s) and executes through
the original scalarTransport class.

For variable diffusivity, remove `D` and set:

    diffusivityModel PbI2He;

Also provide `constant/speciesTransportProperties` with `pressureUnit bar;`
(or `atm`). Case 009 contains a complete example of these dictionaries.
The solver constructs `D_PbI2_g` [m2/s] and `rhoD_PbI2_g` [kg/(m s)] once
from local T, absolute p and saved rho, including boundary values. Both
fields are written at normal output times and recalculated on restart.

`PbI2HeDiffusivity.H` reproduces the active fit in Francesco's Python scripts:

    reducedT = T / 71.11785702515257
    omega = 0.5313 - 1.2946*exp(-0.9922*reducedT)
            + 1.5591/reducedT - 0.1918/reducedT^2
    D = 1.86e-7 * T^1.5 * sqrt(1/4.0026 + 1/461)
        / ((p/pressureScale) * 4.2584682449728355^2 * omega)

`pressureScale` is 100000 Pa for bar, or 101325 Pa for atm. Bar matches the
script's comment; the paper instead labels pressure in atm. This choice is
explicit and configurable. The collision parameters come from the script's
alpha(PbI2)=17.23 Angstrom^3, N=16 and the tabulated He properties. N's
physical interpretation still needs clarification; it is preserved exactly.
The spline evaluation is not used because the script overwrites it with
this fit. The molecular mass 461 g/mol reproduces the diffusion script;
0.46100894 kg/mol is still used for the molar-concentration output.

The variable branch retains nCorr+1 as the maximum number of scalar solves,
checks the initial residual against `tolerance` after each solve, and uses
`schemesField` for div schemes, solver controls and equation relaxation.
An exhausted corrector budget produces a warning. The old-time field is
not advanced between correctors. No momentum, energy or pressure equations
are solved. Variable-branch fvOptions sources are explicitly rejected.

## Verification

Numerical reference values are recorded in `tests/testDiffusivity.C`.
After loading OpenFOAM in WSL, build and run this independent C++ test:

    g++ -std=c++17 -Wall -Wextra tests/testDiffusivity.C -o /tmp/testDiffusivity
    /tmp/testDiffusivity

Then `bash Allwmake` builds the solver. Copy the new headers as well as the
.C source when transferring to WSL. OpenFOAM compilation and CFD execution
still require verification in WSL; neither was available in this Windows
session. Case 009's readme explains how to restore its carrier-field links.
