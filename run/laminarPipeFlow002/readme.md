# Laminar Pipe Flow 002 – Parabolic Inlet

This case is identical to `laminarPipeFlow001`, except that the inlet velocity profile is parabolic instead of uniform.

## Difference from case 001

Case `001` used a uniform inlet velocity:

    U = (0.095 0 0) m/s

Case `002` uses the laminar parabolic profile

    Ux(r) = Umax * (1 - r^2/R^2)

with

    Umax = 0.19 m/s
    R    = 0.0024 m

where

    r^2 = y^2 + z^2

The corresponding bulk velocity is approximately 0.095 m/s, so the flow rate is essentially the same as in case `001`.

## Modified file

Only

    0/U

was changed.

The `INLET` boundary condition was changed from a uniform `fixedValue` to a `codedFixedValue`.

OpenFOAM compiles this code at run time into a small shared library under

    dynamicCode/parabolicInlet/

The internal velocity field is still initialized uniformly as

    internalField uniform (0.095 0 0);

so only the inlet boundary condition is parabolic.

## Run

As in case `001`:

    simpleFoam > out &

First converged solution:
- 388 SIMPLE iterations

