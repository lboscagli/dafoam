// Energy equation for DAContrailFoam
// Ported from CASSANDRA contrailFoam/hEqn_aero.H and aligned with DAResidualContrailFoam

{
    volScalarField& he = thermo.he();

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, he)
      + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K)
      + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)")
          : -dpdt
        )
      - fvm::laplacian(turbulence.alphaEff(), he)
     ==
        fvOptions(rho, he)
    );

    hEqn.relax();
    fvOptions.constrain(hEqn);

    SolverPerformance<scalar> solverH = hEqn.solve();
    DAUtility::primalResidualControl(
        solverH, pimplePrintToScreen, he.name(), daGlobalVarPtr_->primalMaxRes
    );

    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(thermo.T()).value() << ", "
        << max(thermo.T()).value() << endl;
}
