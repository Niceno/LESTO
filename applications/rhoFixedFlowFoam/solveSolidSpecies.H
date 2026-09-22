/*------------------------------------------------------------------------------

solveSolidSpecies.H

PURPOSE

Solves the transient accumulation equation for one solid species during one
physical time step of rhoFixedFlowFoam.  This file is included from inside the
species loop, so speciesi identifies the solid species currently being solved.

At this development step the solid has no convection and no diffusion.  It is
created everywhere by a prescribed uniform volumetric source:

  ddt(rho,Y_i) = S_i

where S_i is read from the species dictionary as `source` in kg/(m3 s).
Because rho is frozen, a uniform source produces

  Y_i(t) = Y_i(0) + S_i*t/rho

for an initially steady rho field.  This gives a simple analytical check before
the source is later restricted to wall-adjacent cells.

The equation contains only a transient diagonal term, so it is solved once per
physical time step with OpenFOAM's diagonal linear solver.  No nCorr loop or
under-relaxation is required.

------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------
Obtain the current solid field and construct the dimensioned source.
The source has units kg/(m3 s), matching ddt(rho,Y) for dimensionless Y.
----------------------------------------------------------------------*/
volScalarField& Y = species[speciesi];
const word& fieldName = Y.name();

const dimensionedScalar sourceTerm (
  "source",
  dimDensity/dimTime,
  solidSource[speciesi]
);

/*------------------------------------------------------------
Assemble and solve the accumulation-only equation.
There are deliberately no convection or diffusion terms here.
-----------------------------------------------------------*/
fvScalarMatrix YEqn (
  fvm::ddt(rho, Y)
  ==
  sourceTerm
);

YEqn.solve(fieldName);

Info<< fieldName << " min/max = " << gMin(Y.primitiveField())
  << " " << gMax(Y.primitiveField()) << nl;
