/*---------------------------------------------------------------------------*\
  rhoFixedFlowFoam -- fixed carrier flow, transient scalar transport, OF v2606.

  Derived from the initialization and program structure of rhoSimpleFoam:
  Copyright (C) 2011-2017 OpenFOAM Foundation.
  Distributed under GNU GPL version 3 or (at your option) any later version.
  https://www.gnu.org/licenses/gpl-3.0.html

  Constant D retains the installed OpenCFD scalarTransport implementation.
  PbI2He diffusivity uses an explicit density-weighted species equation.
  Y_PbI2_g is a mass fraction; c_PbI2_g is derived output in mol/m3.
\*---------------------------------------------------------------------------*/

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

  // Drive the existing scalar function explicitly, once per physical step.
  // Otherwise Time::run() would also execute it through the function list.
  // Other controlDict function objects are not executed by this prototype.
  runTime.functionObjects().off();

  Info<< "Reading frozen thermophysical state" << nl;
  autoPtr<fluidThermo> pThermo(fluidThermo::New(mesh));
  fluidThermo& thermo = pThermo();
  thermo.validate(args.executable(), "e");
  thermo.he().writeOpt() = IOobject::AUTO_WRITE;

  // Unlike the parent solver, require the SAVED density and corrected mass
  // flux: silently reconstructing phi from U would change the carrier flow.
  volScalarField rho (
    IOobject("rho", runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
  );

  volVectorField U (
    IOobject("U", runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
  );

  surfaceScalarField phi (
    IOobject("phi", runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
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

  // These properties are evaluated once and remain frozen as well.
  volScalarField mu (
    IOobject("mu", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE), thermo.mu()
  );

  volScalarField kappa (
    IOobject("kappa", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE), thermo.kappa()
  );

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

  // Preserve the existing constant-D route exactly for case 008.
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
      IOobject("Y_PbI2_g", runTime.timeName(), mesh,
                     IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
    ));

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

  // Molar mass of PbI2: 207.2 + 2*126.90447 = 461.00894 g/mol.
  // Use kg/mol here so c_PbI2_g is in mol/m3, not kmol/m3.
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

    // No UEqn, EEqn, pEqn, thermo.correct or turbulence.correct here.
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

      // Algebraic output only: includes cell and patch values.
      // Recompute from the current mass fraction at each write time.
      volScalarField cPbI2 (
        IOobject (
          "c_PbI2_g", runTime.timeName(), mesh,
          IOobject::NO_READ, IOobject::NO_WRITE
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
