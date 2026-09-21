/*------------------------------------------------------------------------------

solveVariableSpecies.H

PURPOSE

Solves the transient transport equation for the gaseous PbI2 mass fraction
Y_PbI2_g during one physical time step of rhoFixedFlowFoam.

The frozen carrier fields rho and phi are not modified.  The equation solved is

  ddt(rho,Y) + div(phi,Y) - laplacian(rhoD,Y) = 0

where rhoD = rho*D_PbI2_g.  Up to nCorr+1 scalar solves are performed within
the same physical time step.  The old-time field remains fixed during these
correctors; only the current-time estimate is iterated toward the requested
outer tolerance.

OPENFOAM OBJECTS USED HERE

  - fvScalarMatrix is OpenFOAM's finite-volume linear system for a scalar field.
  - fvm::ddt, fvm::div and fvm::laplacian add implicit matrix contributions for
    the transient, convection and diffusion terms.
  - YEqn.solve() then uses the linear solver selected for schemesField in
    system/fvSolution.

------------------------------------------------------------------------------*/

/*-------------------------------------------------------------------
Obtain a writable reference to the species field and, if requested in
fvSolution, the equation under-relaxation factor for this scalar.
-------------------------------------------------------------------*/

volScalarField& Y = species();
scalar relaxation = 0;

if (mesh.solution().relaxEquation(schemesField)) {
  relaxation = mesh.solution().equationRelaxationFactor(schemesField);
}

/*-------------------------------------------------------------------------
Build the fvSchemes lookup name, e.g. div(phi,Y_PbI2_g), and initialise the
outer-corrector convergence bookkeeping.
-------------------------------------------------------------------------*/
const word divScheme("div(phi," + schemesField + ")");
scalar initialResidual = GREAT;
bool converged = false;

/*----------------------------------------------------------------------------
Repeatedly assemble and solve the same physical-time equation.  OpenFOAM's fvm
operators build an implicit finite-volume matrix.  The old-time contribution
from the previous physical time remains unchanged throughout these correctors.
----------------------------------------------------------------------------*/
for (label corr = 0; corr <= nCorr; ++corr) {

  fvScalarMatrix YEqn (
    fvm::ddt(rho, Y)
      + fvm::div(phi, Y, divScheme)
      - fvm::laplacian(rhoD(), Y)
  );

  /*----------------------------------------------------------------------
  Apply equation relaxation, solve with the fvSolution settings associated
  with schemesField, and use the reported initial residual as the outer
  convergence measure.
  ----------------------------------------------------------------------*/
  YEqn.relax(relaxation);
  initialResidual = YEqn.solve(schemesField).initialResidual();

  if (initialResidual < tolerance) {
    Info<< "Y_PbI2_g is converged: initial residual " << initialResidual
      << " < " << tolerance << "; outer iteration " << corr+1 << nl;
    converged = true;
    break;
  }
}

if (!converged && tolerance > 0) {
  WarningInFunction << "Y_PbI2_g reached " << nCorr+1
    << " scalar solves without meeting outer tolerance " << tolerance
    << ". Last initial residual: " << initialResidual << nl;
}

