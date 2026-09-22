/*------------------------------------------------------------------------------

rhoFixedFlowFoam -- fixed carrier flow, transient species transport, OF v2606.

Derived from the initialization and program structure of rhoSimpleFoam:
Copyright (C) 2011-2017 OpenFOAM Foundation.
Distributed under GNU GPL version 3 or (at your option) any later version.
https://www.gnu.org/licenses/gpl-3.0.html

PURPOSE

rhoFixedFlowFoam advances a list of passive gaseous species on an already
converged carrier flow.  The carrier fields U, p, T/e, rho and phi are read
from disk and remain fixed; no momentum, pressure, continuity or energy
equation is solved.

Solid species may also be listed.  Their Y_<speciesName> fields are read,
registered and written.  At this development step a solid species is advanced
only by a prescribed volumetric source; it has no convection or diffusion.
The source model is selected per solid species and is evaluated again at every
physical time step.

For every gaseous species i, the solver reads the mass-fraction field
Y_<speciesName> and solves

  ddt(rho,Y_i) + div(phi,Y_i) - laplacian(rhoD_i,Y_i) = 0,

where rhoD_i = rho*D_i.  The molecular diffusivity D_i is selected separately
for each species in constant/speciesTransportProperties.  The currently
supported models are:

  - constant : a constant molecular diffusivity D_i in m2/s;
  - PbI2He   : the temperature- and pressure-dependent PbI2-in-He correlation.

The species themselves are listed at run time.  One species is therefore just
a list containing one entry; adding another passive gaseous species does not
require another solver.

If a species dictionary also contains molarMass in kg/mol, the solver writes
the derived molar concentration c_<speciesName> = rho*Y_i/M_i in mol/m3.

PROGRAM FLOW

  - create case/time database and mesh
  - read the frozen thermodynamic and carrier-flow fields
  - read the species list and their properties
  - create Y_i for every species and rhoD_i for gaseous species

  At the present solid-species development step, a solid Y_i field is advanced
  only by ddt(rho,Y_i) = source.  The source is evaluated at every physical
  time step.  No rhoD_i is created for a solid species.

  - while physical time advances {
      for every species {
        solve the gas transport equation or the solid accumulation equation
      }
      write requested fields
    }

OPENFOAM STYLE

The program follows the usual OpenFOAM solver structure.  Typical OF objects
used in this program include:
  - runTime manages physical time and output;
  - mesh is the finite-volume mesh and object database;
  - IOobject describes how fields are read and written;
  - IOdictionary reads OpenFOAM dictionaries from disk;
  - volScalarField, volVectorField and surfaceScalarField store cell/face data;
  - PtrList stores a run-time-sized list of OpenFOAM objects;
  - fluidThermo provides the thermodynamic state;
  - Info/WarningInFunction/FatalErrorInFunction are OpenFOAM's standard
    reporting streams.

The included .H files are code fragments inserted directly at the #include
locations, following common OpenFOAM solver practice.

------------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "PbI2HeDiffusivity.H"

int main(int argc, char *argv[]) {

  argList::addNote (
    "Transient species transport on a frozen rhoSimpleFoam solution."
  );

  #include "setRootCaseLists.H"
  #include "createTime.H"
  #include "createMesh.H"

  /*--------------------------------------------------------------------------
  fluidThermo is OpenFOAM's thermodynamic model object.  The New(mesh) factory
  reads thermophysicalProperties and constructs the model selected there.
  autoPtr is OpenFOAM's owning smart pointer.

  Please note that "he" is OpenFOAM's enthalpy (h) or energy (e) variable.
  --------------------------------------------------------------------------*/
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
  ---------------------------------------------------------------------------*/
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

  /*----------------------------------------------------------------------
  Evaluate viscosity and thermal conductivity once from the frozen thermo
  state.  They are volume-scalar fields and are marked AUTO_WRITE, but are
  never updated because thermo.correct() is never called.
  ----------------------------------------------------------------------*/
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
  speciesTransportProperties contains the species list and one sub-dictionary
  for each species.  Gaseous species select a diffusivity model and transport
  properties; solid species specify their accumulation-source model.
  --------------------------------------------------------------------------*/

  /*
  Solid species appear in the same list.  They do not select or use a
  diffusivity model; at this stage they use either a prescribed constant source
  or a simple temperature-dependent source.
  */
  IOdictionary speciesProperties (
    IOobject("speciesTransportProperties",
              runTime.constant(),
              mesh,
              IOobject::MUST_READ,
              IOobject::NO_WRITE)
  );

  const wordList speciesNames(speciesProperties.lookup("species"));

  if (speciesNames.size() == 0) {
    FatalErrorInFunction << "The species list must contain at least one species."
      << exit(FatalError);
  }

  /*--------------------------------------------------------------------------
  Numerical controls shared by all gaseous species are kept in controlDict.
  The linear-solver and under-relaxation settings remain field-specific in
  fvSolution, while convection schemes remain field-specific in fvSchemes.
  --------------------------------------------------------------------------*/
  const dictionary& transportControls =
    runTime.controlDict().subDict("speciesTransport");

  const label nCorr = transportControls.getOrDefault<label>("nCorr", 0);
  const scalar tolerance =
    transportControls.getOrDefault<scalar>("tolerance", 1);

  if (nCorr < 0 || !std::isfinite(tolerance) || tolerance < 0) {
    FatalErrorInFunction << "Require nCorr >= 0 and finite tolerance >= 0."
      << exit(FatalError);
  }

  /*-------------------------------------------------------------------------
  PtrList is OpenFOAM's owning list of objects.  The number of fields is not
  known when the solver is compiled; it is determined by the species list in
  speciesTransportProperties when the case starts.

  species[i]   stores Y_<speciesName>.
  rhoD[i]      stores rho*D_i, the coefficient used in the diffusion term.
  molecularD[i] is used only by models which create a spatial D_i field, such
                as PbI2He; unused PtrList entries remain null.
  molarMass[i] is negative when no molar-concentration output is requested.
  -------------------------------------------------------------------------*/

  /*
  For a solid species only species[i] is populated.  Its rhoD and molecularD
  entries remain null because the solid equation has no diffusion term.
  solidSourceModel[i] selects the source model.  solidSource[i] stores S0 in
  kg/(m3 s); solidThot[i] and solidTcold[i] are used by the temperatureLinear
  model.
  */
  PtrList <volScalarField> species         (speciesNames.size());
  PtrList <volScalarField> rhoD            (speciesNames.size());
  PtrList <volScalarField> molecularD      (speciesNames.size());
  List <scalar>            molarMass       (speciesNames.size(), -1.0);
  List <word>              diffusionModel  (speciesNames.size());
  List <word>              state           (speciesNames.size());
  List <word>              solidSourceModel(speciesNames.size());
  List <scalar>            solidSource     (speciesNames.size(), 0.0);
  List <scalar>            solidThot       (speciesNames.size(), 0.0);
  List <scalar>            solidTcold      (speciesNames.size(), 0.0);

  forAll(speciesNames, speciesi) {

    const word &       speciesName = speciesNames[speciesi];
    const word         fieldName("Y_" + speciesName);
    const word         rhoDName("rhoD_" + speciesName);
    const dictionary & speciesDict = speciesProperties.subDict(speciesName);
    state[speciesi]                = speciesDict.get<word>("state");

    /*----------------------
    Check for allowed states
    ----------------------*/
    if (state[speciesi] != "gas" && state[speciesi] != "solid") {
      FatalErrorInFunction
        << "Unknown state " << state[speciesi]
        << " for species " << speciesName
        << ". Choose gas or solid."
        << exit(FatalError);
    }

    /*--------------------------------------------------------------------
    A solid species has no diffusivity.  Its accumulation source is selected
    independently.  The source is evaluated later, inside the physical-time
    loop, so future source models may depend on evolving transported fields.
    --------------------------------------------------------------------*/
    if (state[speciesi] == "solid") {
      solidSourceModel[speciesi] = speciesDict.getOrDefault<word> (
        "sourceModel", "constant"
      );

      if (solidSourceModel[speciesi] != "constant"
        && solidSourceModel[speciesi] != "temperatureLinear") {
        FatalErrorInFunction << "Unknown sourceModel "
          << solidSourceModel[speciesi] << " for solid species " << speciesName
          << ". Choose constant or temperatureLinear." << exit(FatalError);
      }

      solidSource[speciesi] = speciesDict.get<scalar>("source");

      if (!std::isfinite(solidSource[speciesi]) || solidSource[speciesi] < 0) {
        FatalErrorInFunction << "Species " << speciesName
          << " requires a finite solid source >= 0. Got "
          << solidSource[speciesi] << exit(FatalError);
      }

      if (solidSourceModel[speciesi] == "temperatureLinear") {
        solidThot[speciesi] = speciesDict.get<scalar>("Thot");
        solidTcold[speciesi] = speciesDict.get<scalar>("Tcold");

        if (!std::isfinite(solidThot[speciesi])
          || !std::isfinite(solidTcold[speciesi])
          || solidThot[speciesi] <= solidTcold[speciesi]) {
          FatalErrorInFunction << "Species " << speciesName
            << " requires finite Thot > Tcold for sourceModel "
            << "temperatureLinear. Got Thot=" << solidThot[speciesi]
            << ", Tcold=" << solidTcold[speciesi] << exit(FatalError);
        }

        Info<< "Species " << speciesName
          << ": state = solid, sourceModel = temperatureLinear, S0 = "
          << solidSource[speciesi] << " kg/(m3 s), Thot = "
          << solidThot[speciesi] << " K, Tcold = "
          << solidTcold[speciesi] << " K" << nl;

      } else {
        Info<< "Species " << speciesName
          << ": state = solid, sourceModel = constant, source = "
          << solidSource[speciesi] << " kg/(m3 s)" << nl;
      }
    }

    /*
    Read and register the species field for both allowed states.  AUTO_WRITE
    ensures that an unchanged solid field is still written at output times.
    */
    species.set (
      speciesi,
      new volScalarField (
        IOobject(fieldName,
                 runTime.timeName(),
                 mesh,
                 IOobject::MUST_READ,
                 IOobject::AUTO_WRITE),
        mesh
      )
    );

    if (species[speciesi].dimensions() != dimless) {
      FatalErrorInFunction << fieldName << " must be dimensionless."
        << exit(FatalError);
    }

    if (state[speciesi] == "gas") {

      diffusionModel[speciesi] = speciesDict.getOrDefault<word> (
        "diffusivityModel", "constant"
      );

      if (diffusionModel[speciesi] != "constant"
        && diffusionModel[speciesi] != "PbI2He") {
        FatalErrorInFunction << "Unknown diffusivityModel "
          << diffusionModel[speciesi] << " for species " << speciesName
          << ". Choose constant or PbI2He." << exit(FatalError);
      }

      if (diffusionModel[speciesi] == "constant") {

        const scalar D = speciesDict.get<scalar>("D");
        if (!std::isfinite(D) || D < 0) {
          FatalErrorInFunction << "Species " << speciesName
            << " requires finite D >= 0. Got " << D << exit(FatalError);
        }

        /*--------------------------------------------------------------------
        dimensionedScalar carries both the numerical value and the OpenFOAM
        dimensions.  D is the molecular diffusivity [m2/s]; multiplying by rho
        creates the equation coefficient rhoD [kg/(m s)].
        --------------------------------------------------------------------*/
        const dimensionedScalar constantD (
          "D", dimViscosity, D
        );

        rhoD.set (
          speciesi,
          new volScalarField (
            IOobject(rhoDName,
                     runTime.timeName(),
                     mesh,
                     IOobject::NO_READ,
                     IOobject::AUTO_WRITE),
            rho*constantD
          )
        );

        Info<< "Species " << speciesName << ": field " << fieldName
          << ", diffusivityModel = constant, D = " << D << " m2/s" << nl;

      } else {

        if (speciesName != "PbI2_g") {
          FatalErrorInFunction << "The PbI2He model is specific to PbI2_g, not "
            << speciesName << "." << exit(FatalError);
        }

        if (speciesDict.found("D")) {
          FatalErrorInFunction << "Species " << speciesName
            << " uses diffusivityModel PbI2He: remove constant D."
            << exit(FatalError);
        }

        /*---------------------------------------------------------------------
        Insert the setup code which evaluates the same PbI2-He correlation used
        by the original single-species case 009.  It creates D_PbI2_g(T,p) and
        rhoD_PbI2_g = rho*D_PbI2_g once from the frozen carrier fields.
        ---------------------------------------------------------------------*/
        #include "createVariableDiffusivity.H"
      }

    }

    /*----------------------------------------------------------------------
    Store the starting species field as the old-time field used by the Euler
    transient term.  This is required by both the gas transport equation and
    the solid accumulation equation.  For PbI2_g this remains the same
    operation as in the original single-species case 009.
    ----------------------------------------------------------------------*/
    species[speciesi].oldTime();

    if (speciesDict.found("molarMass")) {
      molarMass[speciesi] = speciesDict.get<scalar>("molarMass");
      if (!std::isfinite(molarMass[speciesi]) || molarMass[speciesi] <= 0) {
        FatalErrorInFunction << "Species " << speciesName
          << " requires molarMass > 0 when specified. Got "
          << molarMass[speciesi] << exit(FatalError);
      }
    }
  }

  Info<< "Frozen fields: U, p, e, T, rho, phi, mu, kappa" << nl
    << "Configured species: " << speciesNames << nl
    << "No SIMPLE/PIMPLE, pressure, momentum or energy solve." << nl
    << "Species controls are read at startup; restart to change them." << nl
    << "Starting physical time loop" << endl;

  while (runTime.loop()) {

    Info<< "Time = " << runTime.timeName() << nl << endl;

    /*-------------------------------------------------------------------------
    The carrier flow remains frozen.  The same transport equation is assembled
    and solved independently for every gaseous species in the run-time list.
    Since the species are presently uncoupled, solving tracer after PbI2_g does
    not alter the PbI2_g equation or its residual history.
    -------------------------------------------------------------------------*/
    forAll(species, speciesi) {
      if (state[speciesi] == "gas") {
        #include "solveGasSpecies.H"
      } else {
        #include "solveSolidSpecies.H"
      }

      /*----------------------------------------------------------------------
      Solid fields have no convection or diffusion.  Their accumulation source
      is evaluated inside solveSolidSpecies.H at every physical time step.
      ----------------------------------------------------------------------*/
    }

    runTime.write();

    if (runTime.writeTime()) {

      /*----------------------------------------------------------------------
      Molar concentration is a derived output only.  It is written for species
      that define molarMass in speciesTransportProperties.
      ----------------------------------------------------------------------*/
      forAll(species, speciesi) {
        if (molarMass[speciesi] > 0) {
          const dimensionedScalar Mi (
            "Mi", dimensionSet(1, 0, 0, 0, -1, 0, 0), molarMass[speciesi]
          );

          const word concentrationName("c_" + speciesNames[speciesi]);

          volScalarField concentration (
            IOobject(concentrationName,
                     runTime.timeName(),
                     mesh,
                     IOobject::NO_READ,
                     IOobject::NO_WRITE),
            rho*species[speciesi]/Mi
          );

          concentration.write();
          Info<< "Wrote " << concentration.name() << " [mol/m3]" << nl;
        }
      }
    }

    runTime.printExecutionTime(Info);
  }

  Info<< "End" << endl;
  return 0;
}
