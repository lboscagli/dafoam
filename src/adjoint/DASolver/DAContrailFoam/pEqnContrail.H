// Pressure correction for DAContrailFoam
// Ported from CASSANDRA contrailFoam/pEqn.H

{
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    if (pimple.transonic())
    {
        surfaceScalarField phid
        (
            "phid",
            fvc::interpolate(psi) * (fvc::flux(HbyA))
        );

        surfaceScalarField mrfTerm = rhorAUf*fvc::ddtCorr(rho, U, phi)/fvc::interpolate(rho);
        MRF.makeRelative(fvc::interpolate(psi), mrfTerm);
        phid += mrfTerm;

        MRF.makeRelative(fvc::interpolate(psi), phid);

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::ddt(psi, p)
              + fvm::div(phid, p)
              - fvm::laplacian(rhorAUf, p)
             ==
                fvOptions(psi, p, rho.name())
            );

            pEqn.solve(p.select(pimple.finalInnerIter()));

            if (pimple.finalNonOrthogonalIter())
            {
                phi = pEqn.flux();
            }
        }
    }
    else
    {
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            fvc::flux(rho*HbyA)
        );

        surfaceScalarField mrfTerm = rhorAUf*fvc::ddtCorr(rho, U, phi);
        MRF.makeRelative(fvc::interpolate(rho), mrfTerm);
        phiHbyA += mrfTerm;

        MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

        constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::ddt(psi, p)
              + fvm::div(phiHbyA, p)
              - fvm::laplacian(rhorAUf, p)
             ==
                fvOptions(psi, p, rho.name())
            );

            pEqn.solve(p.select(pimple.finalInnerIter()));

            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }

    p.relax();

    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
}
