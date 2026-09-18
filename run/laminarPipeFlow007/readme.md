# Laminar Pipe Flow 007 – Variable Transport Properties

This case is based on `laminarPipeFlow006`.

## Differences from case 006

The helium transport properties are now temperature dependent.

In

    constant/thermophysicalProperties

the transport model was changed from

    const

to

    sutherland

using

    As = 1.4683e-06
    Ts = 83.40

This makes the dynamic viscosity vary with temperature according to
Sutherland's law.

Thermal conductivity also becomes temperature dependent through the same
transport model.

The density remains temperature dependent through

    hePsiThermo + perfectGas

while

    Cp

is still kept constant.

## Modified file

    constant/thermophysicalProperties
    system/controlDict

## Writing viscosity and thermal conductivity

The `functions` section of `system/controlDict` now contains the coded
function object `writeTransportProperties`. At each scheduled write, it
obtains the local properties from the solver's `fluidThermo` object using
`thermo.mu()` and `thermo.kappa()` and explicitly writes two scalar fields:

- `mu`: dynamic viscosity, in Pa s.
- `kappa`: thermal conductivity, in W/(m K).

This exposes the model properties for visualization; it does not change
the transport model or the flow equations. Refresh or reopen the case in
ParaView and select `mu` and `kappa` in the field list.

OpenFOAM compiles the function object automatically on first use and when
its code or compilation settings change. No separate manual compilation
is required, but the OpenFOAM build environment must be available. Generated
sources and libraries are stored under `dynamicCode`.

The compilation setup in `controlDict` includes `fvCFD.H` and `fluidThermo.H`,
the include paths for `thermophysicalModels/basic`,
`thermophysicalModels/specie` and `transportModels/compressible`, and links
against `fluidThermophysicalModels` and `specie`. Keep the trailing backslashes
between successive flags in the multiline `codeOptions` and `codeLibs`
blocks: they are Makefile line continuations.

### Output frequency

The function object uses `writeControl writeTime`, so the property fields
follow the normal solver output schedule. The current global settings are
`writeControl timeStep` and `writeInterval 100`: regular output is therefore
scheduled every 100 SIMPLE iterations. The copied `100` directory contains
both property fields alongside the regular solution fields.

To write all fields every iteration, change the **global** `writeInterval`
to `1` and retain `writeControl writeTime` inside the function object.
Setting `writeControl timeStep` and `writeInterval 1` only inside the
function object writes only `mu` and `kappa` every iteration; it does not
change the output frequency of the other fields.

## Run

As in case `006`:

    rm -rf 0
    cp -r 0.orig 0
    setExprFields
    rhoSimpleFoam > out &

First converged solution:
- 226 SIMPLE iterations


