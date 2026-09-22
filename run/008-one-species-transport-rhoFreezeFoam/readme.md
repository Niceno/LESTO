# Case 008 — Multiple gaseous species on frozen carrier flow

This case continues from the converged steady helium solution of case 007,
but uses `rhoFixedFlowFoam` to solve transient transport of a list of passive
gaseous species.

Case 007 provides the frozen carrier-flow fields `U`, `p`, temperature/internal
energy, `rho` and mass flux `phi`. No momentum, pressure, continuity or energy
equation is solved in this case.

The transported gaseous species are listed in
`constant/speciesTransportProperties`. Each species is represented by a
mass-fraction field `Y_<speciesName>` and obeys

    ddt(rho,Y_i) + div(phi,Y_i) - laplacian(rho*D_i,Y_i) = 0

Case 008 currently contains two species:

    PbI2_g
    tracer

Both use the same constant diffusivity `D = 1e-5 m2/s` and identical initial
and boundary conditions. Therefore `Y_PbI2_g` and `Y_tracer` should remain
numerically identical. This provides a direct check of the new multi-species
loop before species-specific diffusivity models are introduced.

`PbI2_g` also defines its molar mass, so the solver writes the derived field
`c_PbI2_g = rho*Y_PbI2_g/MPbI2` in mol/m3 at scheduled output times.

This case therefore changes the species treatment from one hard-wired scalar
to a run-time list of gaseous species, while retaining constant diffusivity and
no chemistry, deposition or solid-species accumulation.
