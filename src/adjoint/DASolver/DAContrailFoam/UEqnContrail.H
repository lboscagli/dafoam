// Momentum equation for DAContrailFoam
// Ported from CASSANDRA contrailFoam/UEqn.H and aligned with DAResidualContrailFoam

{
    UEqn = fvm::ddt(rho, U)
         + fvm::div(phi, U)
         + MRF.DDt(rho, U)
         + turbulence.divDevRhoReff(U)
         ==
            fvOptions(rho, U);

    UEqn.relax();
    fvOptions.constrain(UEqn);
    UEqn.solve(U.select(pimple.finalInnerIter()));

    U.correctBoundaryConditions();
    fvOptions.correct(U);
}
