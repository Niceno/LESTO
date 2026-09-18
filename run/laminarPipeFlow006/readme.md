# Laminar Pipe Flow 006 – Variable Density

This case is based on `laminarPipeFlow005`.

## Differences from case 005

Helium density is no longer constant.

The thermodynamic model in

    constant/thermophysicalProperties

was changed from

    heRhoThermo + rhoConst

to

    hePsiThermo + perfectGas

Density is therefore calculated from the local pressure and temperature
using the ideal-gas equation of state.

Viscosity, Cp and Pr remain constant.

## Modified file

    constant/thermophysicalProperties

## Run

As in case `005`:

    rm -rf 0
    cp -r 0.orig 0
    setExprFields
    rhoSimpleFoam > out &

First converged solution:
- 230 SIMPLE iterations

