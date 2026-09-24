/*------------------------------------------------------------------------------

solveSolidSpecies.H

PURPOSE

Solves the transient accumulation equation for one solid species during one
physical time step of rhoFixedFlowFoam.  This file is included from inside the
solid-species loop, so speciesi identifies the solid species currently being
solved.

At this development step the solid has no convection and no diffusion.  It is
created by the volumetric source returned by evaluateThermochemistry():

  ddt(rho,Y_i) = S_i

The thermochemistry source is evaluated once per cell at the beginning of the
physical time step and is shared with the gas equation.  The equation contains
only a transient diagonal term, so it is solved once per physical time step
with OpenFOAM's diagonal linear solver.  No nCorr loop or under-relaxation is
required.

------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------
Obtain the current solid field and the source field returned by the local
thermochemistry calculation.  The source has units kg/(m3 s), matching
ddt(rho,Y) for dimensionless Y.
-----------------------------------------------------------------------*/
volScalarField& Y = species[speciesi];
const volScalarField& source = speciesSource[speciesi];
const word& fieldName = Y.name();

/*------------------------------------------------------------------------
Assemble and solve the accumulation-only equation.  There are deliberately
no convection or diffusion terms here.
------------------------------------------------------------------------*/
const dimensionedScalar zeroSource("zeroSource", source.dimensions(), 0);
const dimensionedScalar Yfloor("Yfloor", dimless, 1e-12);
const tmp<volScalarField> sinkCoeff =
  -min(source, zeroSource)/max(Y, Yfloor);

fvScalarMatrix YEqn (
  fvm::ddt(rho, Y)                                       /* unstead term */
  + fvm::Sp(sinkCoeff(), Y)                              /* negative source */
  ==
  max(source, zeroSource)                                /* positive source */
);

YEqn.solve(fieldName);

#include "collectSpeciesBalance.H"

Info<< fieldName << " min/max = " << gMin(Y.primitiveField())
  << " " << gMax(Y.primitiveField()) << nl;
