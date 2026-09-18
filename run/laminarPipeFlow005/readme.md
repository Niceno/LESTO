# Laminar Pipe Flow 005 – Linear Temperature Profile

This case is based on `laminarPipeFlow004`.

## Differences from case 004

The temperature field is no longer uniform.

At the inlet:

    T = 400 K

At the outlet end of the pipe:

    T = 300 K

The wall temperature varies linearly along the pipe axis:

    T(x) = 400 - (100/0.69) x

The internal temperature field is initialized with the same linear
profile using `setExprFields`.

## Modified files

    0.orig/T
    system/setExprFieldsDict

## Run

As before:

    rm -rf 0
    cp -r 0.orig 0
    setExprFields
    rhoSimpleFoam > out &

The linearly variable temperature value on the wall resulted in
significantly faster convergence of the energy equation.

First converged solution:
- 154 SIMPLE iterations

