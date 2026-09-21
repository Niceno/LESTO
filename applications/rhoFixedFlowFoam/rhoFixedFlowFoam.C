/*------------------------------------------------------------------------------

rhoFixedFlowFoam -- fixed carrier flow, transient scalar transport, OF v2606.

Derived from the initialization and program structure of rhoSimpleFoam:
Copyright (C) 2011-2017 OpenFOAM Foundation.
Distributed under GNU GPL version 3 or (at your option) any later version.
https://www.gnu.org/licenses/gpl-3.0.html

PURPOSE

rhoFixedFlowFoam advances gaseous PbI2 on an already converged carrier flow.
The carrier fields U, p, T/e, rho and phi are read from disk and remain fixed;
no momentum, pressure, continuity or energy equation is solved.

The solver supports two species-diffusion routes:

  1. constant diffusion coefficient (like in case 008...)
     Uses OpenFOAM's standard scalarTransport function object.

  2. PbI2He (introduced with case 009...)
     Computes D_PbI2_g(T,p) once from the frozen carrier fields, forms
     rhoD_PbI2_g = rho*D_PbI2_g, and solves the density-weighted species
     equation explicitly.

In both cases Y_PbI2_g is a mass fraction which is a dimensionless quantity.
At output times the solver also writes c_PbI2_g = rho*Y_PbI2_g/MPbI2 in
mol/m3 for post-processing.

PROGRAM FLOW

  - create case/time database and mesh
  - read the frozen thermodynamic and carrier-flow fields
  - select constant or variable diffusivity
  - prepare the species field and numerical controls

  - while physical time advances {
      solve only Y_PbI2_g
      write requested fields
    }

OPENFOAM STYLE

The program follows the usual OpenFOAM solver structure.  Typical OF objects
used in this program include:
  - runTime manages physical time and output;
  - mesh is the finite-volume mesh and object database;
  - IOobject describes how fields are read and written;
  - volScalarField, volVectorField and surfaceScalarField store cell/face data;
  - fluidThermo provides the thermodynamic state;
  - Info/WarningInFunction/FatalErrorInFunction are OpenFOAM's standard
    reporting streams.

The included .H files are code fragments inserted directly at the #include
locations, following common OpenFOAM solver practice.

------------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "scalarTransport.H"
#include "functionObjectList.H"
#include "PbI2HeDiffusivity.H"

int main(int argc, char *argv[]) {

  argList::addNote (
    "Transient Y_PbI2_g transport on a frozen rhoSimpleFoam solution."
  );

  #include "setRootCaseLists.H"
  #include "createTime.H"
  #include "createMesh.H"

  /*-----------------------------------------------------------------------
  OpenFOAM normally lets runTime execute function objects automatically.
  Here the scalar transport is called explicitly at a precise point in each
  physical time step, so automatic execution is disabled.
  -----------------------------------------------------------------------*/
  runTime.functionObjects().off();

  /*--------------------------------------------------------------------------
  fluidThermo is OpenFOAM's thermodynamic model object.  The New(mesh) factory
  reads thermophysicalProperties and constructs the model selected there.
  autoPtr is OpenFOAM's owning smart pointer.

  Please note that "he" is OpenFOAM's enthalpy (h) or energy (e) variable.
  -------------------------------------------------------------------------*/
  Info << "Reading frozen thermophysical state" << nl;
  autoPtr<fluidThermo> pThermo(fluidThermo::New(mesh));
  fluidThermo& thermo = pThermo();
  thermo.validate(args.executable(), "e");
  thermo.he().writeOpt() = IOobject::AUTO_WRITE;

  /*---------------------------------------------------------------------------
  Read the frozen carrier fields.  IOobject specifies the field name, the time
  directory, the mesh database, and the read/write policy.  MUST_READ means the
  saved field must exist; AUTO_WRITE lets runTime.write() write it later.

  The saved, pressure-corrected mass flux phi is required explicitly.
  Reconstructing phi from U would no longer reproduce the converged carrier
  solution.
  -------------------------------------------------------------------------*/
  volScalarField rho (
    IOobject("rho",
              runTime.timeName(),
              mesh,
              IOobject::MUST_READ,
              IOobject::AUTO_WRITE),
    mesh
  );

  volVectorField U (
    IOobject("U",
              runTime.timeName(),
              mesh,
              IOobject::MUST_READ,
              IOobject::AUTO_WRITE),
    mesh
  );

  surfaceScalarField phi (
    IOobject("phi",
              runTime.timeName(),
              mesh,
              IOobject::MUST_READ,
              IOobject::AUTO_WRITE),
    mesh
  );

  if (rho.dimensions() != dimDensity || phi.dimensions() != dimMass/dimTime) {
    FatalErrorInFunction
      << "Require rho in kg/m3 and saved phi in kg/s."
      << exit(FatalError);
  }

  if (gMin(rho.primitiveField()) <= 0) {
    FatalErrorInFunction << "Density must be positive."
      << exit(FatalError);
  }

  rho.oldTime();

  /*-----------------------------------------------------------------------
  Evaluate viscosity and thermal conductivity once from the frozen thermo
  state.  They are OpenFOAM volume-scalar fields and are marked AUTO_WRITE,
  but are never updated because thermo.correct() is never called.
  -----------------------------------------------------------------------*/
  volScalarField mu (
    IOobject("mu",
             runTime.timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::AUTO_WRITE),
    thermo.mu()
  );

  volScalarField kappa (
    IOobject("kappa",
             runTime.timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::AUTO_WRITE),
    thermo.kappa()
  );

  /*--------------------------------------------------------------------------
  OpenFOAM dictionaries are runtime configuration objects.  Read the species
  controls from functions/Y_PbI2_gTransport in system/controlDict and use them
  to select the constant-D or PbI2-He variable-D route.
  --------------------------------------------------------------------------*/
  const dictionary& scalarDict = runTime.controlDict()
    .subDict("functions").subDict("Y_PbI2_gTransport");
  const word diffusionModel = scalarDict.getOrDefault<word> (
    "diffusivityModel", "constant"
  );

  const bool variableD = diffusionModel == "PbI2He";
  if (diffusionModel != "constant" && !variableD) {
    FatalErrorInFunction << "Unknown diffusivityModel: " << diffusionModel
      << ". Choose constant or PbI2He." << exit(FatalError);
  }

  if (
    scalarDict.get<word>("field") != "Y_PbI2_g"
     || scalarDict.getOrDefault<word>("phi", "phi") != "phi"
     || scalarDict.getOrDefault<word>("rho", "rho") != "rho"
     || scalarDict.getOrDefault<word>("phase", "none") != "none"
     || scalarDict.getOrDefault<bool>("resetOnStartUp", false)
     || (!variableD && (!scalarDict.found("D")
                         || scalarDict.get<scalar>("D") < 0))
  ) {
    FatalErrorInFunction
      << "This prototype requires field Y_PbI2_g, phi, rho, no phase, "
      << "resetOnStartUp no, and nonnegative D for the constant model."
      << exit(FatalError);
  }

  /*--------------------------------------------------------------------------
  These autoPtr objects are created only when needed.  The constant-D branch
  owns an OpenFOAM scalarTransport object.  The variable-D branch instead owns
  the species field plus molecularD and rhoD fields used by our explicit
  transport equation.
  --------------------------------------------------------------------------*/
  autoPtr<functionObjects::scalarTransport> transport;
  autoPtr<volScalarField> species;
  autoPtr<volScalarField> molecularD;
  autoPtr<volScalarField> rhoD;

  if (variableD) {
    if (scalarDict.found("D") || scalarDict.found("fvOptions")) {
      FatalErrorInFunction
        << "PbI2He uses speciesTransportProperties: remove D. "
        << "fvOptions are not supported by this variable-D branch yet."
        << exit(FatalError);
    }
    species.reset(new volScalarField (
      IOobject("Y_PbI2_g",
               runTime.timeName(),
               mesh,
               IOobject::MUST_READ,
               IOobject::AUTO_WRITE),
      mesh
    ));

    /*--------------------------------------------------------------------
    Insert the setup code which evaluates D_PbI2_g(T,p) over all cells and
    boundary faces and creates rhoD_PbI2_g = rho*D_PbI2_g.
    --------------------------------------------------------------------*/
    #include "createVariableDiffusivity.H"

  } else {

    transport.reset(new functionObjects::scalarTransport (
      "Y_PbI2_gTransport", runTime, scalarDict
    ));
  }

  const word schemesField = scalarDict.getOrDefault<word> (
    "schemesField", "Y_PbI2_g"
  );

  const label nCorr = scalarDict.getOrDefault<label>("nCorr", 0);
  const scalar tolerance = scalarDict.getOrDefault<scalar>("tolerance", 1);
  if (nCorr < 0 || !std::isfinite(tolerance) || tolerance < 0) {
    FatalErrorInFunction << "Require nCorr >= 0 and finite tolerance >= 0."
      << exit(FatalError);
  }

  const volScalarField& PbI2 = mesh.lookupObject<volScalarField>("Y_PbI2_g");
  if (PbI2.dimensions() != dimless) {
    FatalErrorInFunction
      << "This case-008 prototype preserves dimensionless Y_PbI2_g. "
      << "Molar concentration requires a separate conversion."
      << exit(FatalError);
  }
  PbI2.oldTime();

  /*-------------------------------------------------------------------------
  Molar mass of PbI2: 207.2 + 2*126.90447 = 461.00894 g/mol.
  Use kg/mol here so c_PbI2_g is produced in mol/m3 rather than kmol/m3.
  dimensionedScalar carries both the numerical value and OpenFOAM dimensions.
  -------------------------------------------------------------------------*/
  const dimensionedScalar MPbI2 (
    "MPbI2", dimensionSet(1, 0, 0, 0, -1, 0, 0), 0.46100894
  );

  Info<< "Frozen fields: U, p, e, T, rho, phi, mu, kappa" << nl
    << "Only Y_PbI2_g is advanced. Diffusivity model: " << diffusionModel << nl
    << "Equation diffusion coefficient has units kg/(m s)." << nl
    << "No SIMPLE/PIMPLE, pressure, momentum or energy solve." << nl
    << "Scalar controls are read at startup; restart to change them." << nl
    << "Starting physical time loop" << endl;

  while (runTime.loop()) {

    Info<< "Time = " << runTime.timeName() << nl << endl;

    /*-------------------------------------------------------------------------
    This is the essential difference from rhoSimpleFoam/rhoPimpleFoam: the
    physical-time loop contains no carrier-flow equations or thermo correction.
    Only the species field is advanced.
    -------------------------------------------------------------------------*/
    if (variableD) {
      #include "solveVariableSpecies.H"

    } else {
      transport->execute();

    }
    Info<< "Y_PbI2_g min/max = " << gMin(PbI2.primitiveField())
      << " " << gMax(PbI2.primitiveField()) << nl;

    runTime.write();
    if (runTime.writeTime()) {
      if (transport.valid()) transport->write();

      /*--------------------------------------------------------------------
      c_PbI2_g is a derived output field, not an independently transported
      unknown.  It is reconstructed algebraically from rho and the current
      mass fraction whenever output is requested, including boundary values.
      --------------------------------------------------------------------*/
      volScalarField cPbI2 (
        IOobject ("c_PbI2_g",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::NO_WRITE
        ),
        rho*PbI2/MPbI2
      );
      cPbI2.write();
      Info<< "Wrote c_PbI2_g [mol/m3]" << nl;
    }
    runTime.printExecutionTime(Info);
  }

  Info<< "End" << endl;
  return 0;
}
