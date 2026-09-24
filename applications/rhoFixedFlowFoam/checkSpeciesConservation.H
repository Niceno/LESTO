/*------------------------------------------------------------------------------

Check the discrete inventory, boundary-flux and realized-source balance for
each configured species and their sum. The equation flux is collected from the
last gas solve; solid species have no boundary transport. Coupled patches are
internal to the computational domain and are excluded from the boundary flux.

The sum covers the configured species, not the untransported helium carrier.

------------------------------------------------------------------------------*/

scalar totalOldMass = 0;
scalar totalNewMass = 0;
scalar totalStorageRate = 0;
scalar totalInRate = 0;
scalar totalOutRate = 0;
scalar totalSourceRate = 0;

Info << nl;

forAll(species, speciesi) {
  reduce(speciesBoundaryInRate[speciesi], sumOp<scalar>());
  reduce(speciesBoundaryOutRate[speciesi], sumOp<scalar>());

  const volScalarField& Y = species[speciesi];
  const scalar oldMass =
    fvc::domainIntegrate(rho.oldTime()*Y.oldTime()).value();
  const scalar newMass = fvc::domainIntegrate(rho*Y).value();
  const scalar storageRate = fvc::domainIntegrate(fvc::ddt(rho, Y)).value();
  const scalar residualRate =
    storageRate + speciesBoundaryOutRate[speciesi]
      - speciesBoundaryInRate[speciesi] - realizedSourceRate[speciesi];
  const scalar scale =
    mag(storageRate) + speciesBoundaryInRate[speciesi]
      + speciesBoundaryOutRate[speciesi] + mag(realizedSourceRate[speciesi]);

  Info << "Species balance " << speciesNames[speciesi] << ":" << nl
    << "  Inventory [kg]: old=" << oldMass
    << "  new=" << newMass << nl
    << "  Rates [kg/s]: storage=" << storageRate
    << "  source=" << realizedSourceRate[speciesi] << nl
    << "  Boundary [kg/s]: in=" << speciesBoundaryInRate[speciesi]
    << "  out=" << speciesBoundaryOutRate[speciesi] << nl
    << "  Error [kg/s]: " << residualRate
    << "  relative=" << mag(residualRate)/max(scale, VSMALL)
    << nl << nl;

  totalOldMass += oldMass;
  totalNewMass += newMass;
  totalStorageRate += storageRate;
  totalInRate += speciesBoundaryInRate[speciesi];
  totalOutRate += speciesBoundaryOutRate[speciesi];
  totalSourceRate += realizedSourceRate[speciesi];
}

const scalar totalResidualRate =
  totalStorageRate + totalOutRate - totalInRate - totalSourceRate;
const scalar totalScale =
  mag(totalStorageRate) + totalInRate + totalOutRate + mag(totalSourceRate);

Info<< "Tracked species total:" << nl
  << "  Inventory [kg]: old=" << totalOldMass
  << "  new=" << totalNewMass << nl
  << "  Rates [kg/s]: storage=" << totalStorageRate
  << "  source=" << totalSourceRate << nl
  << "  Boundary [kg/s]: in=" << totalInRate
  << "  out=" << totalOutRate << nl
  << "  Error [kg/s]: " << totalResidualRate
  << "  relative=" << mag(totalResidualRate)/max(totalScale, VSMALL)
  << nl;
