/*------------------------------------------------------------------------------

PbI2HeDiffusivity.H

PURPOSE

Implements Francesco's PbI2-in-He molecular diffusivity correlation as a
small standalone C++ function.  Given temperature T and absolute pressure p,
it returns the molecular diffusivity in m2/s.

This file contains no OpenFOAM fields and solves no transport equation.  It
is included by rhoFixedFlowFoam and called by createVariableDiffusivity.H for
each cell and boundary face when the PbI2He diffusivity model is selected.

pressureScale converts the OpenFOAM pressure in Pa to the pressure unit used
by the correlation.  For the present model pressureUnit is atm, so
pressureScale = 101325 Pa/atm.

------------------------------------------------------------------------------*/

#ifndef PbI2HeDiffusivity_H
#define PbI2HeDiffusivity_H

#include <cmath>

/*----------------------------------------------------------------------------
Francesco's utilities/diffusivity.py and interpolation.py, active fit.
T is in K, p is absolute pressure in Pa, and pressureScale is the number of Pa
per pressure unit expected by the correlation.  The return value is molecular
D in m2/s; it is not the density-weighted equation coefficient rho*D.
----------------------------------------------------------------------------*/
namespace LESTO {

  inline double pbI2HeDiffusivity (const double T,
                                   const double p,
                                   const double pressureScale) {

    /*------------------------------------------------------------------------
    Collision parameters derived from alpha(PbI2)=17.23 Angstrom^3, N=16 and
    tabulated helium data.  The script's molecular mass 461 g/mol is preserved
    here so that the OpenFOAM implementation reproduces the Python curve.
    ------------------------------------------------------------------------*/
    const double sigma      = 4.2584682449728355;  /* Angstrom */
    const double epsilonByK = 71.11785702515257;   /* K */
    const double reducedT   = T/epsilonByK;
    const double omega      = 0.5313
                            - 1.2946*std::exp(-0.9922*reducedT)
                            + 1.5591/reducedT
                            - 0.1918/(reducedT*reducedT);

    /*------------------------------------------------------------------------
    The original coefficient 1.86e-3 gives D in cm2/s.  Multiplication by 1e-4
    converts cm2/s to m2/s, hence 1.86e-7 below.  p/pressureScale is the
    numerical pressure in the unit required by the correlation.
    ------------------------------------------------------------------------*/
    return   1.86e-7
           * std::pow(T, 1.5)
           * std::sqrt(1.0/4.0026 + 1.0/461.0)
           / ((p/pressureScale) * sigma * sigma * omega);

  }
}
#endif
