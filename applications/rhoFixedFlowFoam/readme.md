# rhoFixedFlowFoam

## Purpose

`rhoFixedFlowFoam` advances one or more passive gaseous species on an already
converged `rhoSimpleFoam` carrier flow.

The carrier fields `U`, `p`, `e/T`, `rho` and `phi` are read from the saved
steady solution and remain frozen throughout the transient species calculation.
No momentum, pressure, continuity or energy equation is solved, and there is no
SIMPLE/PIMPLE loop.

The gaseous species are defined at run time in
`constant/speciesTransportProperties`.  The solver therefore does not need to
be rebuilt when another passive gaseous species is added.  A one-species case
is simply a species list containing one entry.

## Transported variables and scalar equation

For every gaseous species `i`, the transported field is

    Y_<speciesName>

and represents a dimensionless mass fraction:

    Y_i = kg species i / kg mixture.

Each species is advanced in physical time with

    ddt(rho,Y_i) + div(phi,Y_i) - laplacian(rhoD_i,Y_i) = 0

where

    rhoD_i = rho * D_i.

`D_i` is the molecular diffusivity in m2/s.  The coefficient entering the
finite-volume diffusion term is therefore `rhoD_i`, with units kg/(m s).

## Frozen carrier fields

`U`, `p`, `e/T`, `rho` and `phi` are read from the converged carrier-flow
solution.  The saved, pressure-corrected mass flux `phi` is required directly;
it is not reconstructed from `U`, because doing so could change the converged
carrier solution.

The solver also evaluates `mu` and `kappa` once from the frozen thermodynamic
state and writes them at normal output times.  They are not updated because
`thermo.correct()` is never called.

## Numerical controls

The corrector controls shared by all gaseous species are read from the
`speciesTransport` dictionary in `system/controlDict`, for example:

    speciesTransport {
      tolerance  1.0e-7;
      nCorr      47;
    }

## Molar concentration output

If a species dictionary contains

    molarMass <value>;

with the molar mass in kg/mol, the solver also writes

    c_<speciesName> = rho * Y_i / M_i

in mol/m3 at every scheduled output time.

For gaseous PbI2,

    molarMass 0.46100894;

corresponds to

    MPbI2 = 207.2 + 2*126.90447 = 461.00894 g/mol.

The concentration field is derived output only.  It is not transported and no
additional equation is solved.  No `0/c_<speciesName>` input file is required.
The field is reconstructed from local `rho` and the current mass fraction at
write time, including boundary values, and is recalculated after a restart.
Existing time directories are not retroactively updated.

For PbI2 the resulting field `c_PbI2_g` has dimensions
`[0 -3 0 0 1 0 0]`; its numerical values are in mol/m3, not kmol/m3.  Select
`c_PbI2_g` in ParaView when molar concentration is desired for post-processing.

## PbI2-He variable diffusivity

For

    diffusivityModel PbI2He;

the solver evaluates a spatial molecular-diffusivity field `D_PbI2_g` in m2/s
from the frozen local `T` and absolute `p`, including boundary values, and then
forms

    rhoD_PbI2_g = rho * D_PbI2_g

in kg/(m s).

Because the carrier state is frozen, `D_PbI2_g` and `rhoD_PbI2_g` are evaluated
once when the solver starts and remain fixed during the transient species run.
Both are written at normal output times and are recalculated when the solver is
restarted.

`PbI2HeDiffusivity.H` reproduces the active fit used in Francesco's Python
`diffusivity.py` / `interpolation.py` calculation:

    reducedT = T / 71.11785702515257

    omega = 0.5313
            - 1.2946*exp(-0.9922*reducedT)
            + 1.5591/reducedT
            - 0.1918/reducedT^2

    D = 1.86e-7 * T^1.5 * sqrt(1/4.0026 + 1/461)
        / ((p/pressureScale) * 4.2584682449728355^2 * omega)

The original correlation coefficient `1.86e-3` gives diffusivity in cm2/s;
the factor `1e-4` converting cm2/s to m2/s is already included above, hence
`1.86e-7`.

The collision parameters are derived from `alpha(PbI2)=17.23 Angstrom^3`,
`N=16` and the tabulated helium properties.  The physical interpretation of
`N` still needs clarification; its numerical value is preserved exactly from
the Python implementation.  The spline evaluation in the Python utilities is
not used here because the script overwrites it with the explicit fit above.

The diffusivity calculation preserves the Python script's molecular mass
`461 g/mol` so that the OpenFOAM implementation reproduces the same curve.  The
more precise `0.46100894 kg/mol` is used separately for the derived molar-
concentration output.

## Current limitations

- All transported species are currently passive gaseous species.
- Species equations are uncoupled from each other.
- No chemistry or thermodynamic-equilibrium source terms are included.
- No deposition or precipitation is included yet.
- No solid-species accumulation equation is included yet.
- The carrier flow is frozen and is not modified by species transport.
- `PbI2He` is currently the only variable-diffusivity model and is restricted
  to `PbI2_g`.
- Other gaseous species currently use constant molecular diffusivity unless a
  new diffusivity model is implemented.
- Species transport currently has no `fvOptions` source-term path.
- Carrier flow must already be converged.

## Verification of the PbI2-He correlation

Numerical reference values are recorded in `tests/testDiffusivity.C`.  This is
an independent small C++ test and does not require OpenFOAM libraries.

From the solver directory:

    g++ -std=c++17 -Wall -Wextra tests/testDiffusivity.C -o /tmp/testDiffusivity
    /tmp/testDiffusivity

The test checks reference diffusivity values from the same explicit Python
correlation and verifies the inverse pressure scaling, including the bar/atm
conversion factor.  The numerical reference curve corresponds to a correlation
pressure value of 1; the current CFD cases use the atm convention described
above.

