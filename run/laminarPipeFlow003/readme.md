# Laminar Pipe Flow 003 – Parabolic Initial Field

This case is identical to `laminarPipeFlow002`, except that the initial
velocity field is also parabolic instead of uniform.

This makes the initial condition much more consistent with the parabolic
inlet profile and leads to significantly faster convergence.

## Difference from case 002

Case `002` used a parabolic inlet profile, but the internal velocity field
was still initialized uniformly as

    internalField uniform (0.095 0 0);

Case `003` initializes the complete velocity field using the same
parabolic profile as the inlet:

    Ux(r) = Umax * (1 - r^2/R^2)

with

    Umax = 0.19 m/s
    R    = 0.0024 m

and

    r^2 = y^2 + z^2

The inlet boundary condition itself is unchanged from case `002`.

## Modified files

The initial field in

    0/U

was changed from a uniform field to a nonuniform parabolic field.

The field was generated using

    system/setExprFieldsDict

and applied with

    setExprFields -time 0

After this step, `0/U` contains a nonuniform `internalField`.

## Run

As in case `002`:

    simpleFoam > out &

First converged solution:
- 222 SIMPLE iterations

