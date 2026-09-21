#ifndef PbI2HeDiffusivity_H
#define PbI2HeDiffusivity_H

#include <cmath>

// Francesco's utilities/diffusivity.py and interpolation.py, active fit.
// T in K, absolute pressure in Pa, pressureScale in Pa per correlation unit.
// Return molecular diffusivity in m2/s. This is NOT the coefficient rho*D.
namespace LESTO {
inline double pbI2HeDiffusivity (
  const double T, const double p, const double pressureScale) {

  // Derived from alpha(PbI2)=17.23 Angstrom^3, N=16 and tabulated He.
  // Preserve the script's molecular mass 461 g/mol for curve reproduction.
  const double sigma = 4.2584682449728355;       // Angstrom
  const double epsilonByK = 71.11785702515257;  // K
  const double reducedT = T/epsilonByK;
  const double omega = 0.5313 - 1.2946*std::exp(-0.9922*reducedT)
    + 1.5591/reducedT - 0.1918/(reducedT*reducedT);

  // 1.86e-3 yields cm2/s; factor 1e-4 converts to m2/s.
  return 1.86e-7*std::pow(T, 1.5)*std::sqrt(1.0/4.0026 + 1.0/461.0)
    / ((p/pressureScale)*sigma*sigma*omega);

  }
}
#endif
