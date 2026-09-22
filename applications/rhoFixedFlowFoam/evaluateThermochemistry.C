#include "evaluateThermochemistry.H"
#include "error.H"

namespace Foam {

void evaluateThermochemistry (
  const scalar        T,
  const scalar        p,
  const wordList    & speciesNames,
  const scalarField & speciesValues,
        scalarField & sources) {

  if (speciesValues.size() != speciesNames.size()
    ||sources.size()       != speciesNames.size()) {
    FatalErrorInFunction
      << "speciesNames, speciesValues and sources must have the same size."
      << exit(FatalError);
  }

  /*-------------------------------------------------------------------------
  Mock thermochemistry model.

  These are intentionally local test parameters, not case input.  They mimic
  quantities that would eventually be determined by the thermochemistry/GEMS
  calculation rather than by the CFD transport dictionaries.

  For now only PbI2_s receives a non-zero source:

      S_PbI2_s = S0 * (Thot - T)/(Thot - Tcold)

  This retains the temperature dependence used in the previous development
  step.  Pressure and the species values are passed through the interface but
  are not used yet; they are already available for the future GEMS model.
  -------------------------------------------------------------------------*/
  const scalar S0    = 1.0e-6;  /* kg/(m3 s) */
  const scalar Thot  = 1000.0;  /* K */
  const scalar Tcold = 300.0;   /* K */

  sources = scalar(0);

  forAll(speciesNames, speciesi) {
    if (speciesNames[speciesi] == "PbI2_s") {
      sources[speciesi] = S0*(Thot - T)/(Thot - Tcold);
    }
  }

  /*------------------------------------------------------------------------
  Silence unused-argument warnings in this mock implementation.  Both values
  belong to the interface because a real thermochemistry call will need the
  complete local state.
  ------------------------------------------------------------------------*/
  (void) p;
  (void) speciesValues;
}

}
