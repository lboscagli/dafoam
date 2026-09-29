// Species transport for DAContrailFoam
// Ported from CASSANDRA contrailFoam/YEqn.H

{
    // Solve for all active species except the inert one
    forAll(Y, speciesI)
    {
        if (speciesI == inertIndex) continue;

        fvScalarMatrix YEqn
        (
            fvm::ddt(rho, Y[speciesI])
          + mvConvection->fvmDiv(phi, Y[speciesI])
          - fvm::laplacian(turbulence.muEff()/Sc, Y[speciesI])
         ==
            fvOptions(rho, Y[speciesI])
        );

        YEqn.solve(Y[speciesI].select(pimple.finalInnerIter()));
    }

    // Reconstruct inert species (N2) to ensure sum-to-one
    volScalarField sumY(0.0*Y[0]);
    forAll(Y, speciesI)
    {
        if (speciesI != inertIndex)
        {
            sumY += Y[speciesI];
        }
    }

    Y[inertIndex] = scalar(1) - sumY;
    Y[inertIndex].clamp_min(0);
    Y[inertIndex].correctBoundaryConditions();
}
