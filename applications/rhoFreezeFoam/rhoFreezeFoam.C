/*---------------------------------------------------------------------------*\
  rhoFreezeFoam -- frozen carrier flow, transient scalar transport, OF v2606.

  Derived from the initialization and program structure of rhoSimpleFoam:
  Copyright (C) 2011-2017 OpenFOAM Foundation.
  Distributed under GNU GPL version 3 or (at your option) any later version.
  https://www.gnu.org/licenses/gpl-3.0.html

  The species equation is the installed OpenCFD scalarTransport implementation.
  This first version deliberately retains case 008's dimensionless Pb_g and
  density-weighted D convention. It is NOT yet a mol/m3 chemistry solver.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "scalarTransport.H"
#include "functionObjectList.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient Pb_g transport on a frozen rhoSimpleFoam solution."
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
    volScalarField rho
    (
        IOobject("rho", runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
    );
    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
    );
    surfaceScalarField phi
    (
        IOobject("phi", runTime.timeName(), mesh,
                 IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh
    );

    if (rho.dimensions() != dimDensity || phi.dimensions() != dimMass/dimTime)
    {
        FatalErrorInFunction
            << "Require rho in kg/m3 and saved phi in kg/s."
            << exit(FatalError);
    }
    if (gMin(rho.primitiveField()) <= 0)
    {
        FatalErrorInFunction << "Density must be positive."
            << exit(FatalError);
    }
    rho.oldTime();

    // These properties are evaluated once and remain frozen as well.
    volScalarField mu
    (
        IOobject("mu", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE), thermo.mu()
    );
    volScalarField kappa
    (
        IOobject("kappa", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::AUTO_WRITE), thermo.kappa()
    );

    const dictionary& scalarDict = runTime.controlDict()
        .subDict("functions").subDict("Pb_gTransport");
    if
    (
        scalarDict.get<word>("field") != "Pb_g"
     || scalarDict.getOrDefault<word>("phi", "phi") != "phi"
     || scalarDict.getOrDefault<word>("rho", "rho") != "rho"
     || scalarDict.getOrDefault<word>("phase", "none") != "none"
     || scalarDict.getOrDefault<bool>("resetOnStartUp", false)
     || !scalarDict.found("D")
     || scalarDict.get<scalar>("D") < 0
    )
    {
        FatalErrorInFunction
            << "This prototype requires field Pb_g, phi, rho, no phase, "
            << "resetOnStartUp no, and an explicit nonnegative D."
            << exit(FatalError);
    }

    functionObjects::scalarTransport transport
    (
        "Pb_gTransport", runTime, scalarDict
    );
    const volScalarField& Pb = mesh.lookupObject<volScalarField>("Pb_g");
    if (Pb.dimensions() != dimless)
    {
        FatalErrorInFunction
            << "This case-008 prototype preserves dimensionless Pb_g. "
            << "Molar concentration requires a separate conversion."
            << exit(FatalError);
    }
    Pb.oldTime();

    Info<< "Frozen fields: U, p, e, T, rho, phi, mu, kappa" << nl
        << "Only Pb_g is advanced. D is in kg/(m s), not m2/s." << nl
        << "No SIMPLE/PIMPLE, pressure, momentum or energy solve." << nl
        << "Scalar controls are read at startup; restart to change them." << nl
        << "Starting physical time loop" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        // No UEqn, EEqn, pEqn, thermo.correct or turbulence.correct here.
        transport.execute();
        Info<< "Pb_g min/max = " << gMin(Pb.primitiveField())
            << " " << gMax(Pb.primitiveField()) << nl;

        runTime.write();
        if (runTime.writeTime())
        {
            transport.write();
        }
        runTime.printExecutionTime(Info);
    }
    Info<< "End" << endl;
    return 0;
}
