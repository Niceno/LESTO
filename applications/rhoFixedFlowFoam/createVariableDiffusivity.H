// Included once at startup, inside the variable-D branch.
IOdictionary properties
(
    IOobject("speciesTransportProperties", runTime.constant(), mesh,
             IOobject::MUST_READ, IOobject::NO_WRITE)
);
const word pressureUnit = properties.get<word>("pressureUnit");
if (pressureUnit != "bar" && pressureUnit != "atm") {
    FatalErrorInFunction << "pressureUnit must be bar or atm."
        << exit(FatalError);
}
const scalar pressureScale = pressureUnit == "bar" ? 1e5 : 101325;
const volScalarField& T = thermo.T();
const volScalarField& p = thermo.p();
if (T.dimensions() != dimTemperature || p.dimensions() != dimPressure) {
    FatalErrorInFunction << "Require T in K and absolute p in Pa."
        << exit(FatalError);
}

// The expression retains constraint patch types (e.g. processor/cyclic).
molecularD.reset(new volScalarField
(
    IOobject("D_PbI2_g", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    T*dimensionedScalar("zeroDPerK", dimViscosity/dimTemperature, 0)
));
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
    const scalar value = LESTO::pbI2HeDiffusivity
    (
        temperature, pressure, pressureScale
    );
    if (!std::isfinite(value) || value <= 0) {
        FatalErrorInFunction << "Invalid PbI2He diffusivity at T="
            << temperature << ", p=" << pressure << exit(FatalError);
    }
    return value;
};
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
        molecularD().boundaryFieldRef()[patchi][facei] = evaluateD
        (
            temperature, p.boundaryField()[patchi][facei],
            rho.boundaryField()[patchi][facei]
        );
        minT = min(minT, temperature);
        maxT = max(maxT, temperature);
    }
}
reduce(minT, minOp<scalar>());
reduce(maxT, maxOp<scalar>());
if (minT < 250 || maxT > 1250) {
    WarningInFunction << "T range " << minT << " to " << maxT
        << " K extends beyond the script's plotted 250--1250 K range. "
        << "The fit is being extrapolated; values are not clamped." << nl;
}
rhoD.reset(new volScalarField
(
    IOobject("rhoD_PbI2_g", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    rho*molecularD()
));
Info<< "PbI2He diffusivity evaluated once using local T and absolute p."
    << nl << "Correlation pressure unit: " << pressureUnit
    << " (" << pressureScale << " Pa)" << nl
    << "D_PbI2_g internal min/max [m2/s]: "
    << gMin(molecularD().primitiveField()) << " "
    << gMax(molecularD().primitiveField()) << nl;
