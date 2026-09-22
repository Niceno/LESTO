/*------------------------------------------------------------------------------

solveSolidSpecies.H

PURPOSE

Solves the transient accumulation equation for one solid species during one
physical time step of rhoFixedFlowFoam.  This file is included from inside the
species loop, so speciesi identifies the solid species currently being solved.

At this development step the solid has no convection and no diffusion.  It is
created by a prescribed volumetric source:

  ddt(rho,Y_i) = S_i

The source is evaluated again at every physical time step.  Two simple source
models are currently available:

  constant
      S_i = S0

  temperatureLinear
      S_i = S0 * (Thot - T)/(Thot - Tcold)

For the present test case T is frozen, so the temperature-dependent source is
spatially varying but does not change in time.  Re-evaluating it here keeps the
implementation ready for later source models which may depend on evolving
species fields.

The equation contains only a transient diagonal term, so it is solved once per
physical time step with OpenFOAM's diagonal linear solver.  No nCorr loop or
under-relaxation is required.

------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------
Obtain the current solid field.  The source has units kg/(m3 s), matching
ddt(rho,Y) for dimensionless Y.
----------------------------------------------------------------------*/
volScalarField& Y = species[speciesi];
const word& fieldName = Y.name();

const dimensionedScalar S0 (
  "S0",
  dimDensity/dimTime,
  solidSource[speciesi]
);

/*----------------------------------------------------------------------------
Assemble and solve the accumulation-only equation.  There are deliberately no
convection or diffusion terms here.  The temperatureLinear source is evaluated
from the current thermo.T() field every time this file is included.
----------------------------------------------------------------------------*/
if (solidSourceModel[speciesi] == "temperatureLinear") {

  const dimensionedScalar Thot (
    "Thot", dimTemperature, solidThot[speciesi]
  );

  const dimensionedScalar Tcold (
    "Tcold", dimTemperature, solidTcold[speciesi]
  );

  const volScalarField& T = thermo.T();

  fvScalarMatrix YEqn (
    fvm::ddt(rho, Y)
    ==
    S0*(Thot - T)/(Thot - Tcold)
  );

  YEqn.solve(fieldName);

} else {

  fvScalarMatrix YEqn (
    fvm::ddt(rho, Y)
    ==
    S0
  );

  YEqn.solve(fieldName);
}

Info<< fieldName << " min/max = " << gMin(Y.primitiveField())
  << " " << gMax(Y.primitiveField()) << nl;
