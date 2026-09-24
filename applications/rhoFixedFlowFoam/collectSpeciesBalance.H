/*------------------------------------------------------------------------------

Collect the realized source for the current species after its final equation
solve. For a gaseous species, collect the integrated equation flux across each
external boundary face. Solid species have no boundary transport.

This file is included once per species, while YEqn, Y, source, zeroSource,
sinkCoeff and speciesi are still in scope. Gas fluxes are local to each MPI
rank; checkSpeciesConservation.H sums them over ranks before reporting.

------------------------------------------------------------------------------*/

realizedSourceRate[speciesi] =
  fvc::domainIntegrate(max(source, zeroSource) - sinkCoeff()*Y).value();

if (state[speciesi] == "gas") {
  const tmp<surfaceScalarField> tFaceFlux = YEqn.flux();

  forAll(mesh.boundary(), patchi) {
    if (!mesh.boundary()[patchi].coupled()) {
      const fvsPatchScalarField& patchFlux =
        tFaceFlux().boundaryField()[patchi];

      forAll(patchFlux, facei) {
        if (patchFlux[facei] >= 0) {
          speciesBoundaryOutRate[speciesi] += patchFlux[facei];
        } else {
          speciesBoundaryInRate[speciesi] -= patchFlux[facei];
        }
      }
    }
  }
}
