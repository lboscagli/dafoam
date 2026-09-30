// Pressure correction for DAContrailFoam
// Ported from CASSANDRA contrailFoam/pEqn.H (non-transonic branch) and aligned
// with DAResidualContrailFoam. The transonic and consistent (pcEqn) branches are
// rejected in solvePrimal, because DAResidualContrailFoam only evaluates the
// non-transonic pressure residual.

{
    rho = thermo.rho();

    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    if (pimple.nCorrPISO() <= 1)
    {
        tUEqn.clear();
    }

    // Mirrors CASSANDRA pEqn.H (MRF.zeroFilter(...) dropped: the zone list is
    // empty for contrailFoam, so zeroFilter is the identity, and MRFZoneListDF
    // provides no zeroFilter method). The ddtCorr term is required for parity;
    // it is exactly zero on the first step (no old times yet) and matches
    // DAResidualContrailFoam, which applies the same term.
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rho*HbyA)
          + rhorAUf*fvc::ddtCorr(rho, U, phi)
        )
    );

    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phiHbyA)
          - fvm::laplacian(rhorAUf, p)
         ==
            fvOptions(psi, p, rho.name())
        );

        SolverPerformance<scalar> solverP =
            pEqn.solve(p.select(pimple.finalInnerIter()));
        DAUtility::primalResidualControl(
            solverP, pimplePrintToScreen, "p", daGlobalVarPtr_->primalMaxRes
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }

    #include "rhoEqnContrail.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
    K = 0.5*magSqr(U);

    if (pressureControlPtr_->limit(p))
    {
        p.correctBoundaryConditions();
    }

    rho = thermo.rho();

    if (thermo.dpdt())
    {
        dpdt = fvc::ddt(p);
    }
}
