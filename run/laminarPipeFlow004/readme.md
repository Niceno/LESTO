# Laminar Pipe Flow 004 – rhoSimpleFoam with Energy Equation

This case is based on `laminarPipeFlow003`.

The purpose is to introduce `rhoSimpleFoam` and the energy equation while
keeping the physics otherwise as close as possible to case `003`.

## Differences from case 003

Solver:

    simpleFoam

was replaced by

    rhoSimpleFoam

Pressure is now physical pressure in Pa:

    internalField = 101325 Pa

with

    INLET:  zeroGradient
    OUTLET: fixedValue 101325 Pa
    WALL:   zeroGradient

A temperature field was added:

    T = 300 K

with

    INLET:  fixedValue 300 K
    OUTLET: zeroGradient
    WALL:   fixedValue 300 K

The thermophysical properties of helium are defined in

    constant/thermophysicalProperties

using

    rho = 0.166 kg/m3
    mu  = 2.70e-5 Pa s
    Cp  = 5193 J/(kg K)
    Pr  = 0.916

The equation of state is still

    rhoConst

and the energy variable is

    sensibleInternalEnergy

The energy equation `e` was added to `system/fvSchemes` and
`system/fvSolution`.

The SIMPLE convergence criteria were changed to

    U = 1e-4
    p = 1e-4
    e = 1e-3

from the previous

    U = 1e-5
    p = 1e-5

## Run

As in case `003`, regenerate the initial velocity field first:

    rm -rf 0
    cp -r 0.orig 0
    setExprFields

Then run:

    rhoSimpleFoam > out &

First converged solution:
- 154 SIMPLE iterations

