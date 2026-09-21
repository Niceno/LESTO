// One physical time step, up to nCorr+1 scalar solves, as in scalarTransport.
// The old-time solution stays fixed throughout these correctors.
volScalarField& Y = species();
scalar relaxation = 0;
if (mesh.solution().relaxEquation(schemesField)) {
    relaxation = mesh.solution().equationRelaxationFactor(schemesField);
}
const word divScheme("div(phi," + schemesField + ")");
scalar initialResidual = GREAT;
bool converged = false;
for (label corr = 0; corr <= nCorr; ++corr) {
    fvScalarMatrix YEqn
    (
        fvm::ddt(rho, Y)
      + fvm::div(phi, Y, divScheme)
      - fvm::laplacian(rhoD(), Y)
    );
    YEqn.relax(relaxation);
    initialResidual = YEqn.solve(schemesField).initialResidual();
    if (initialResidual < tolerance) {
        Info<< "Y_PbI2_g is converged: initial residual " << initialResidual
            << " < " << tolerance << "; outer iteration " << corr+1 << nl;
        converged = true;
        break;
    }
}
if (!converged && tolerance > 0) {
    WarningInFunction << "Y_PbI2_g reached " << nCorr+1
        << " scalar solves without meeting outer tolerance " << tolerance
        << ". Last initial residual: " << initialResidual << nl;
}
