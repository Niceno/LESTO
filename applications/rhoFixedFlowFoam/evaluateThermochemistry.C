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

  /*--------------------------------------------------------------------------
  Mock thermochemistry model.

  These are intentionally local test parameters, not case input.  They mimic
  quantities that would eventually be determined by the thermochemistry/GEMS
  calculation rather than by the CFD transport dictionaries.

  Case 011 introduces the first conservative interaction between two species:

      PbI2_g  ->  PbI2_s

  The precipitation rate retains the temperature dependence from case 010
  and is proportional to the local Y_PbI2_g.  The interaction is the equal
  and opposite transfer between gaseous and solid PbI2.

      Sprecip = K * max((Thot - T)/(Thot - Tcold), 0) * max(Y_PbI2_g, 0)

      S_PbI2_g = -Sprecip
      S_PbI2_s = +Sprecip

  The equal and opposite source terms conserve PbI2 mass locally.
  --------------------------------------------------------------------------*/
  const scalar Thot  = 1000.0;  /* K */
  const scalar Tcold = 300.0;   /* K */

  sources = scalar(0);

  label gasIndex   = -1;
  label solidIndex = -1;

  forAll(speciesNames, speciesi) {
    if (speciesNames[speciesi] == "PbI2_g") {
      gasIndex = speciesi;
    }
    if (speciesNames[speciesi] == "PbI2_s") {
      solidIndex = speciesi;
    }
  }

  if (gasIndex >= 0 && solidIndex >= 0) {

    const scalar K = 1.0;           /* mock deposition rate kg/(m3 s) ...
                                       ... per unit mass fraction */
    const scalar zero = scalar(0);  /* OpenFOAM's zero scalar */

    const scalar precipitation =
      K * max((Thot - T)/(Thot - Tcold), zero)
        * max(speciesValues[gasIndex], zero);

    sources[gasIndex]   = -precipitation;
    sources[solidIndex] =  precipitation;
  }

  /*-------------------------------------------------------------------------
  Pressure belongs to the interface because a real thermochemistry call will
  need the complete local state.  The present mock model does not use it yet.
  -------------------------------------------------------------------------*/
  (void) p;
}

}
