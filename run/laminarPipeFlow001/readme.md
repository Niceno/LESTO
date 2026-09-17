# Laminar Pipe Flow 001 – OpenFOAM v2606

First custom OpenFOAM case, based on Francesco's Gmsh pipe geometry.

Francesco's cases are in:

~/Development/Francesco_2026_02_12/Gpu/T-Flows-Bojan-2/Tests/Laminar/Scalar_Transport_Gpu/LESTO_Paper/Temperature_Based_Model/

and:

~/Development/Francesco_2026_02_12/Gpu/T-Flows-Bojan-2/Tests/Laminar/Scalar_Transport_Gpu/LESTO_Paper/Vapor_Pressure_Model/

Current geometry was based on the vapor-based model, the only difference is
that the boundary condition INTERNAL_WALL was renamed to WALL.

## Geometry / mesh

Pipe:
- Radius: 2.4 mm
- Length: 0.69 m
- Axis: x
- Mesh: structured hexahedral Gmsh mesh
- Cells: 228960

Generate mesh:

    gmsh -3 pipe.geo

Convert to OpenFOAM:

    gmshToFoam pipe.msh

Check mesh:

    checkMesh

Expected result:

    Mesh OK.

Boundary patches:
- INLET
- OUTLET
- WALL

## Flow conditions

Helium:
- rho = 0.166 kg/m3
- mu  = 2.70e-5 Pa s
- nu  = 1.6265e-4 m2/s

Inlet:
- Uniform U = (0.095 0 0) m/s

Outlet:
- U: zeroGradient
- p: fixedValue 0

Wall:
- U: noSlip
- p: zeroGradient

Flow model:
- incompressible
- steady
- laminar
- simpleFoam

## Important files

    0/U
    0/p

    constant/transportProperties
    constant/turbulenceProperties

    system/controlDict
    system/fvSchemes
    system/fvSolution

This was a steady SIMPLE calculation using iterative pseudo-time-like
progression to converge to the steady state, not a physical transient
simulation.

## Run

    simpleFoam > out &

First converged solution:
- 235 SIMPLE iterations

## ParaView

Create marker file:

    touch laminarPipeFlow.foam

Open:

    paraview laminarPipeFlow.foam
