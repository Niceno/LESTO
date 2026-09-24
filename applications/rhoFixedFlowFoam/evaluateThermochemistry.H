#ifndef evaluateThermochemistry_H
#define evaluateThermochemistry_H

#include "scalarField.H"
#include "wordList.H"

/*------------------------------------------------------------------------------

evaluateThermochemistry

PURPOSE

Provides one local thermochemistry interface for all species.  The solver calls
this routine once per cell at the beginning of each physical time step, before
the species transport equations are advanced.

The routine receives the complete local state currently available to the CFD
solver:

  - temperature T [K];
  - absolute pressure p [Pa];
  - the run-time species-name list;
  - the current local species mass fractions.

It returns one volumetric source for every species in kg/(m3 s).

The present implementation is deliberately only a mock-up of a future GEMS
call.  Keeping the interface separate from rhoFixedFlowFoam means that the
simple test model can later be replaced by GEMS without restructuring the CFD
transport loop.

------------------------------------------------------------------------------*/
namespace Foam {

  void evaluateThermochemistry (
    const scalar        T,
    const scalar        p,
    const wordList    & speciesNames,
    const scalarField & speciesValues,
          scalarField & sources
  );

}

#endif
