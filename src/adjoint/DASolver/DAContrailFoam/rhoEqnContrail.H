// Explicit continuity update for DAContrailFoam
// Ported from CASSANDRA contrailFoam/rhoEqn.H with the PBE mass-coupling
// source dropped: massCouplingActive is false for the reduced gas-phase scope,
// so the right hand side is only the generic fvOptions source.

{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho)
      + fvc::div(phi)
     ==
        fvOptions(rho)
    );

    rhoEqn.solve();

    fvOptions.correct(rho);
}
