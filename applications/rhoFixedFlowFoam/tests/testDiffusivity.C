// Standalone numerical regression test; no OpenFOAM libraries needed.
// Build from the solver directory:
// g++ -std=c++17 -Wall -Wextra tests/testDiffusivity.C -o /tmp/testDiffusivity
// /tmp/testDiffusivity
#include "../PbI2HeDiffusivity.H"
#include <iostream>

int main() {
    // Independent reference values from Francesco's original Python equations,
    // using alpha=17.23, N=16 and the He table, p=1 in the script convention.
    const double temperatures[] = {250, 300, 400, 600, 800, 1000, 1250};
    const double reference[] = {
        2.212906e-5, 3.073716e-5, 5.164654e-5, 1.0611711e-4,
        1.7431433e-4, 2.5392934e-4, 3.673821e-4
    };
    for (int i = 0; i < 7; ++i) {
        const double T = temperatures[i];
        const double D = LESTO::pbI2HeDiffusivity(T, 1e5, 1e5);
        if (!std::isfinite(D) || std::abs(D/reference[i] - 1) > 3e-7) {
            std::cerr << "Reference mismatch at " << T << " K: " << D << '\n';
            return 1;
        }
        const double doublePressure = LESTO::pbI2HeDiffusivity(T, 2e5, 1e5);
        const double atm = LESTO::pbI2HeDiffusivity(T, 1e5, 101325);
        if (std::abs(doublePressure/D - 0.5) > 1e-12
            || std::abs(atm/D - 1.01325) > 1e-12) {
            std::cerr << "Pressure scaling mismatch\n";
            return 1;
        }
    }
    std::cout << "PbI2-He reference values and pressure scaling passed.\n";
}
