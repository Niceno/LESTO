/*------------------------------------------------------------------------------

solveGasSpecies.H

PURPOSE

Solves the transient transport equation for one gaseous species during one
physical time step of rhoFixedFlowFoam.  This file is included from inside the
species loop, so speciesi identifies the species currently being solved.

The frozen carrier fields rho and phi are not modified.  For species i,

  ddt(rho,Y_i) + div(phi,Y_i) - laplacian(rhoD_i,Y_i) = S_i

where rhoD_i = rho*D_i.  Up to nCorr+1 scalar solves are performed within the
same physical time step.  The old-time field remains fixed during these
correctors; only the current-time estimate is iterated toward the requested
outer tolerance.

OPENFOAM OBJECTS USED HERE

  - fvScalarMatrix is OpenFOAM's finite-volume linear system for a scalar field.
  - fvm::ddt, fvm::div and fvm::laplacian add implicit matrix contributions for
    the transient, convection and diffusion terms.
  - YEqn.solve() uses the linear solver selected for the current field name in
    system/fvSolution.

------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
Obtain references to the current species mass fraction and its density-weighted
diffusion coefficient and thermochemistry source.  Y.name() is also the lookup
name used by fvSolution and fvSchemes, e.g. Y_PbI2_g or Y_tracer.
-----------------------------------------------------------------------------*/
      volScalarField & Y         = species[speciesi];
      volScalarField & rhoDi     = rhoD[speciesi];
const volScalarField & source    = speciesSource[speciesi];
const word           & fieldName = Y.name();

/*------------------------------
Read the under-relaxation factor
------------------------------*/
scalar relaxation = 0;
if (mesh.solution().relaxEquation(fieldName)) {
  relaxation = mesh.solution().equationRelaxationFactor(fieldName);
}

/*-------------------------------------------------------------------------
Build the fvSchemes lookup name, e.g. div(phi,Y_PbI2_g), and initialise the
outer-corrector convergence bookkeeping.
-------------------------------------------------------------------------*/
const word divScheme("div(phi," + fieldName + ")");
scalar initialResidual = GREAT;
bool converged = false;

/*----------------------------------------------------------------------------
Repeatedly assemble and solve the same physical-time equation.  OpenFOAM's fvm
operators build an implicit finite-volume matrix.  The old-time contribution
from the previous physical time remains unchanged throughout these correctors.
----------------------------------------------------------------------------*/
const dimensionedScalar    zeroSource("zeroSource", source.dimensions(), 0);
const dimensionedScalar    Yfloor("Yfloor", dimless, 1e-12);
const tmp <volScalarField> sinkCoeff =
  -min(source, zeroSource)/max(Y, Yfloor);

mesh.setFluxRequired(fieldName);

for (label corr = 0; corr <= nCorr; ++corr) {

  fvScalarMatrix YEqn (
    fvm::ddt(rho, Y)                 /* unsteady term */
      + fvm::div(phi, Y, divScheme)  /* advection term */
      - fvm::laplacian(rhoDi, Y)     /* diffusion term */
      + fvm::Sp(sinkCoeff(), Y)      /* negative source */
      ==
      max(source, zeroSource)        /* positive source */
  );

  /*---------------------------------------------------------------------------
  Apply equation relaxation, solve with the fvSolution settings associated with
  fieldName, and use the reported initial residual as the outer convergence
  measure.
  ---------------------------------------------------------------------------*/
  YEqn.relax(relaxation);
  initialResidual = YEqn.solve(fieldName).initialResidual();

  if (initialResidual < tolerance || corr == nCorr) {
    #include "collectSpeciesBalance.H"
  }

  if (initialResidual < tolerance) {
    Info<< fieldName << " is converged: initial residual " << initialResidual
      << " < " << tolerance << "; outer iteration " << corr+1 << nl;
    converged = true;
    break;
  }


}

if (!converged && tolerance > 0) {
  WarningInFunction << fieldName << " reached " << nCorr+1
    << " scalar solves without meeting outer tolerance " << tolerance
    << ". Last initial residual: " << initialResidual << nl;
}

/*-----------------------------------------------------------------------
The gas sink is implicit due to linearization of the source term, so the
amount consumed uses the gas value after the solve, not the source value
evaluated at the start of the time step.  Give the solid that same amount
to conserve the gas-to-solid transfer.
-----------------------------------------------------------------------*/
if (speciesNames[speciesi] == "PbI2_g") {
  forAll(speciesNames, solidi) {
    if (speciesNames[solidi] == "PbI2_s") {
      speciesSource[solidi] = sinkCoeff()*Y;
      speciesSource[speciesi] = -speciesSource[solidi];
      break;
    }
  }
}

Info<< fieldName << " min/max = " << gMin(Y.primitiveField())
  << " " << gMax(Y.primitiveField()) << nl;
