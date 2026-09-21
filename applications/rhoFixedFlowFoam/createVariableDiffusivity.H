/*------------------------------------------------------------------------------

createVariableDiffusivity.H

PURPOSE

Setup code for the variable PbI2-He diffusivity branch of rhoFixedFlowFoam.
It is included once, before the physical-time loop.

The file:

  - reads constant/speciesTransportProperties;
  - obtains the frozen temperature and pressure fields from fluidThermo;
  - evaluates D_PbI2_g(T,p) in every cell and boundary face;
  - creates rhoD_PbI2_g = rho*D_PbI2_g for the species equation;
  - reports the temperature range and diffusivity range.

Because the carrier flow is fixed, these diffusivity fields are calculated
only once at solver start and remain unchanged during the transient run.

OPENFOAM OBJECTS USED HERE

  - IOdictionary reads a standard OpenFOAM dictionary from disk.
  - IOobject specifies where an OpenFOAM object lives and whether it must be
    read or written.
  - volScalarField is a scalar field stored in mesh cells with boundary
    values on patches.
  - forAll is OpenFOAM's usual indexed-loop macro, and
  - reduce() combines values across MPI ranks in parallel runs.

------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------
Read the species-transport dictionary from the constant directory.
MUST_READ makes the solver stop immediately if the file is missing.
-----------------------------------------------------------------*/
IOdictionary properties (
  IOobject("speciesTransportProperties",
           runTime.constant(),
           mesh,
           IOobject::MUST_READ,
           IOobject::NO_WRITE)
);

const word pressureUnit = properties.get<word>("pressureUnit");
if (pressureUnit != "bar" && pressureUnit != "atm") {
  FatalErrorInFunction << "pressureUnit must be bar or atm."
    << exit(FatalError);
}

/*-----------------------------------------------------------------------------
Convert the OpenFOAM pressure, which is stored in Pa, to the numerical pressure
unit expected by the diffusivity correlation.  For the present setup the
selected unit is atm, hence pressureScale = 101325 Pa/atm.
-----------------------------------------------------------------------------*/
const scalar pressureScale = pressureUnit == "bar" ? 1e5 : 101325;
const volScalarField& T = thermo.T();
const volScalarField& p = thermo.p();
if (T.dimensions() != dimTemperature || p.dimensions() != dimPressure) {
  FatalErrorInFunction << "Require T in K and absolute p in Pa."
    << exit(FatalError);
}

/*---------------------------------------------------------------------------
Create the molecular diffusivity field D_PbI2_g.  NO_READ means there is no
input file for D; AUTO_WRITE makes it available in output time directories.
The zero expression based on T gives the field the required dimensions while
retaining OpenFOAM constraint patch types such as processor/cyclic patches.
---------------------------------------------------------------------------*/
molecularD.reset(new volScalarField (
    IOobject("D_PbI2_g",
             runTime.timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::AUTO_WRITE),
    T*dimensionedScalar("zeroDPerK", dimViscosity/dimTemperature, 0)
  )
);

/*-----------------------------------------------------------------------------
Local helper used for both cells and boundary faces.  It validates the supplied
state, calls the standalone PbI2-He correlation, and rejects non-physical
results before they can enter the finite-volume equation.
-----------------------------------------------------------------------------*/
auto evaluateD = [&](const scalar temperature, const scalar pressure,
                     const scalar density) -> scalar {

  if (!std::isfinite(temperature) || !std::isfinite(pressure)
    || !std::isfinite(density) || temperature <= 0
    || pressure <= 0 || density <= 0) {
    FatalErrorInFunction
      << "PbI2He requires finite positive T, absolute p and rho, "
      << "including patch values. Got " << temperature << ", "
      << pressure << ", " << density << exit(FatalError);
  }

  const scalar value = LESTO::pbI2HeDiffusivity (
    temperature, pressure, pressureScale
  );

  if (!std::isfinite(value) || value <= 0) {
    FatalErrorInFunction << "Invalid PbI2He diffusivity at T="
      << temperature << ", p=" << pressure << exit(FatalError);
  }

  return value;
};

/*--------------------------------------------------------------------------
Evaluate D in all internal cells and then on every boundary face.  minT/maxT
are tracked only to warn if the correlation is used outside the temperature
range over which Francesco plotted and checked it.
--------------------------------------------------------------------------*/
scalar minT = GREAT;
scalar maxT = -GREAT;
forAll(T, celli) {
  molecularD()[celli] = evaluateD(T[celli], p[celli], rho[celli]);
  minT = min(minT, T[celli]);
  maxT = max(maxT, T[celli]);
}

forAll(T.boundaryField(), patchi) {
  forAll(T.boundaryField()[patchi], facei) {
    const scalar temperature = T.boundaryField()[patchi][facei];
    molecularD().boundaryFieldRef()[patchi][facei] = evaluateD (
      temperature, p.boundaryField()[patchi][facei],
      rho.boundaryField()[patchi][facei]
    );
    minT = min(minT, temperature);
    maxT = max(maxT, temperature);
  }
}

/*------------------------------------------------------------------------------
In a parallel OpenFOAM run each MPI rank initially knows only its local extrema.
reduce() produces the global minimum and maximum over all ranks.
------------------------------------------------------------------------------*/
reduce(minT, minOp<scalar>());
reduce(maxT, maxOp<scalar>());

if (minT < 250 || maxT > 1250) {
  WarningInFunction << "T range " << minT << " to " << maxT
    << " K extends beyond the script's plotted 250--1250 K range. "
    << "The fit is being extrapolated; values are not clamped." << nl;
}

/*---------------------------------------------------------------------------
The transported variable is a mass fraction.  Its diffusion term therefore
uses the density-weighted coefficient rho*D, with units kg/(m s).  Store this
as rhoD_PbI2_g for direct use by solveVariableSpecies.H.
---------------------------------------------------------------------------*/
rhoD.reset(new volScalarField (
  IOobject("rhoD_PbI2_g",
           runTime.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::AUTO_WRITE),
  rho*molecularD()
  )
);

Info<< "PbI2He diffusivity evaluated once using local T and absolute p."
  << nl << "Correlation pressure unit: " << pressureUnit
  << " (" << pressureScale << " Pa)" << nl
  << "D_PbI2_g internal min/max [m2/s]: "
  << gMin(molecularD().primitiveField()) << " "
  << gMax(molecularD().primitiveField()) << nl;
