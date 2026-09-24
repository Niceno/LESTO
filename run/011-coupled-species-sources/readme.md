# Case 011 — Coupled species sources

This case introduces the first interaction between species: gaseous `PbI2_g`
deposits as solid `PbI2_s`. The gas is transported on the frozen helium carrier
flow; the solid accumulates next to the wall. The passive `tracer` has no
thermochemistry source.

The mock `evaluateThermochemistry()` routine is called once per cell at the
beginning of each physical time step. Its precipitation rate depends on the
local temperature and available gaseous PbI2:

    rate = K * max((Thot - T)/(Thot - Tcold), 0) * max(Y_PbI2_g, 0)

Here `K = 1.0 kg/(m3 s)` per unit mass fraction, `Thot = 1000 K` and
`Tcold = 300 K` are local test constants.

It initially assigns equal and opposite sources to `PbI2_g` and `PbI2_s`.
These sources are applied only in cells adjacent to the patch named `WALL`.

The gas equation transports `PbI2_g` by advection and diffusion and treats
its deposition sink implicitly. The sink coefficient is fixed during the gas
correctors. After the gas solve, the solid source is updated to match the
amount actually consumed by that sink. The solid has no advection or diffusion;
it accumulates locally through `ddt(rho, Y_PbI2_s) = source_PbI2_s`.

The solver prints integrated source rates for each species and their sum.
This checks the gas-to-solid source transfer. A full mass-balance check of
the species inventories and boundary fluxes has not yet been implemented.
