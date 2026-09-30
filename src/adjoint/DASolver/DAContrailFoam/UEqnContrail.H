// Momentum equation for DAContrailFoam
// Ported from CASSANDRA contrailFoam/UEqn.H and aligned with DAResidualContrailFoam
// (the momentum residual is URes_ = (UEqn & U_) + fvc::grad(p_), i.e. the
// momentum-predictor form UEqn == -fvc::grad(p)).
//
// tUEqn/UEqn are deliberately declared without a wrapping scope so that
// pEqnContrail.H, included later in the same PIMPLE iteration, can use UEqn
// and clear tUEqn (CASSANDRA contrailFoam/UEqn.H and pEqn.H pattern).

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  + turbulence.divDevRhoReff(U)
 ==
    fvOptions(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    SolverPerformance<vector> solverU = solve(UEqn == -fvc::grad(p));
    DAUtility::primalResidualControl(
        solverU, pimplePrintToScreen, "U", daGlobalVarPtr_->primalMaxRes
    );

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
