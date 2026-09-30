// Species transport for DAContrailFoam
// Ported from CASSANDRA contrailFoam/YEqn_aero.H (gas-phase, PBE coupling off).
// The solver dictionary name must be "Yi": the species fields are named Y[i],
// which fvSolution does not contain.

{
    // Running sum of solved species mass fractions, used to reconstruct the
    // inert species at the end.
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, speciesI)
    {
        if (speciesI == inertIndex || !composition.active(speciesI))
        {
            continue;
        }

        fvScalarMatrix YEqn
        (
            fvm::ddt(rho, Y[speciesI])
          + mvConvection->fvmDiv(phi, Y[speciesI])
          - fvm::laplacian(turbulence.muEff()/Sc, Y[speciesI])
         ==
            fvOptions(rho, Y[speciesI])
        );

        YEqn.relax();

        fvOptions.constrain(YEqn);

        SolverPerformance<scalar> solverY = YEqn.solve("Yi");
        DAUtility::primalResidualControl(
            solverY, pimplePrintToScreen, "Yi", daGlobalVarPtr_->primalMaxRes
        );

        fvOptions.correct(Y[speciesI]);

        // Clip to prevent negative mass fractions.
        Y[speciesI].clamp_min(0);

        Yt += Y[speciesI];
    }

    // Reconstruct inert species (N2) to ensure sum-to-one
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].clamp_min(0);
    Y[inertIndex].correctBoundaryConditions();
}
