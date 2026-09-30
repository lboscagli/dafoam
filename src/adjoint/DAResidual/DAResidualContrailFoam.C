/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v5

\*---------------------------------------------------------------------------*/

#include "DAResidualContrailFoam.H"

namespace Foam
{

defineTypeNameAndDebug(DAResidualContrailFoam, 0);
addToRunTimeSelectionTable(DAResidual, DAResidualContrailFoam, dictionary);

DAResidualContrailFoam::DAResidualContrailFoam(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
    : DAResidual(modelType, mesh, daOption, daModel, daIndex),
      setResidualClassMemberVector(
          U,
          dimensionSet(1, -2, -2, 0, 0, 0, 0)),
      setResidualClassMemberScalar(
          p,
          dimensionSet(1, -3, -1, 0, 0, 0, 0)),
      setResidualClassMemberScalar(
          T,
          dimensionSet(1, -1, -3, 0, 0, 0, 0)),
      setResidualClassMemberPhi(phi),
      setResidualClassMemberScalar(
          O2,
          dimensionSet(1, -3, -1, 0, 0, 0, 0)),
      setResidualClassMemberScalar(
          CO2,
          dimensionSet(1, -3, -1, 0, 0, 0, 0)),
      setResidualClassMemberScalar(
          H2O,
          dimensionSet(1, -3, -1, 0, 0, 0, 0)),
      thermo_(const_cast<psiReactionThermo&>(
          mesh_.thisDb().lookupObject<psiReactionThermo>(
              "thermophysicalProperties"))),
      composition_(thermo_.composition()),
      Y_(composition_.Y()),
      he_(thermo_.he()),
      rho_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("rho"))),
      alphat_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("alphat"))),
      psi_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("thermo:psi"))),
      dpdt_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("dpdt"))),
      K_(const_cast<volScalarField&>(
          mesh_.thisDb().lookupObject<volScalarField>("K"))),
      daTurb_(const_cast<DATurbulenceModel&>(
          daModel.getDATurbulenceModel())),
      turbulence_(const_cast<compressible::turbulenceModel&>(
          mesh_.thisDb().lookupObject<compressible::turbulenceModel>(
              compressible::turbulenceModel::propertiesName))),
      fvOptions_(fv::options::New(mesh)),
      MRF_(const_cast<IOMRFZoneListDF&>(
          mesh_.thisDb().lookupObject<IOMRFZoneListDF>("MRFProperties"))),
      pimple_(const_cast<fvMesh&>(mesh)),
      inertSpecie_(thermo_.get<word>("inertSpecie")),
      inertIndex_(composition_.species().find(inertSpecie_)),
      Sc_("Sc", dimless, 0.7)
{
    if (inertIndex_ < 0)
    {
        FatalErrorInFunction
            << "Inert species " << inertSpecie_
            << " was not found in the multicomponent mixture."
            << exit(FatalError);
    }

    // Initialise fvOptions before it is used by later residual updates.
    fvVectorMatrix initUEqn
    (
        fvm::div(phi_, U_)
      - fvOptions_(rho_, U_)
    );
    fvOptions_.constrain(initUEqn);
}

void DAResidualContrailFoam::clear()
{
    /*
    Description:
        Clear all members to avoid memory leak because we will initalize 
        multiple objects of DAResidual. Here we need to delete all members
        in the parent and child classes
    */    
    URes_.clear();
    pRes_.clear();
    TRes_.clear();
    phiRes_.clear();
    O2Res_.clear();
    CO2Res_.clear();
    H2ORes_.clear();
}

void DAResidualContrailFoam::calcResiduals(const dictionary& options)
{
    /*
    Description:
        This is the function to compute residuals.
    
    Input:
        options.isPC: 1 means computing residuals for preconditioner matrix.
        This essentially use the first order scheme for div(phi,U), div(phi,e)

        p_, T_, U_, phi_, etc: State variables in OpenFOAM
    
    Output:
        URes_, pRes_, TRes_, phiRes_, etc: residual field variables
    */    
    const label isPC = options.getLabel("isPC");
    const word divUScheme = isPC ? word("div(pc)") : word("div(phi,U)");
    const word divHEScheme = isPC
        ? word("div(pc)")
        : (he_.name() == "h" ? word("div(phi,h)") : word("div(phi,e)"));

    // The standalone CASSANDRA case and first DAFoam milestone use the
    // non-transonic pressure branch.
    if (pimple_.transonic())
    {
        FatalErrorInFunction
            << "DAResidualContrailFoam does not support the transonic "
            << "pressure-residual branch yet."
            << exit(FatalError);
    }

    // Momentum residual: CASSANDRA UEqn.H with no PBE-related source terms.
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rho_, U_)
      + fvm::div(phi_, U_, divUScheme)
      + MRF_.DDt(rho_, U_)
      + turbulence_.divDevRhoReff(U_)
     ==
        fvOptions_(rho_, U_)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();
    UEqn.relax(1.0);
    fvOptions_.constrain(UEqn);

    URes_ = (UEqn & U_) + fvc::grad(p_);
    normalizeResiduals(URes);

    // Match the standalone multivariate convection table for active species
    // and sensible enthalpy. N2 is derived and not an independent residual.
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
    forAll(Y_, speciesI)
    {
        fields.add(Y_[speciesI]);
    }
    fields.add(he_);

    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh_,
            fields,
            phi_,
            mesh_.divScheme("div(phi,Yi_h)")
        )
    );

    // In preconditioner mode use the established first-order PC scheme. In
    // normal residual mode, preserve CASSANDRA's multivariate discretisation.
    tmp<fvScalarMatrix> tO2Eqn
    (
        isPC
        ? tmp<fvScalarMatrix>
          (
              new fvScalarMatrix
              (
                  fvm::ddt(rho_, O2_)
                + fvm::div(phi_, O2_, "div(pc)")
                - fvm::laplacian(turbulence_.muEff()/Sc_, O2_)
               ==
                  fvOptions_(rho_, O2_)
              )
          )
        : tmp<fvScalarMatrix>
          (
              new fvScalarMatrix
              (
                  fvm::ddt(rho_, O2_)
                + mvConvection->fvmDiv(phi_, O2_)
                - fvm::laplacian(turbulence_.muEff()/Sc_, O2_)
               ==
                  fvOptions_(rho_, O2_)
              )
          )
    );
    fvScalarMatrix& O2Eqn = tO2Eqn.ref();
    O2Eqn.relax(1.0);
    fvOptions_.constrain(O2Eqn);
    O2Res_ = O2Eqn & O2_;
    normalizeResiduals(O2Res);

    tmp<fvScalarMatrix> tCO2Eqn
    (
        isPC
        ? tmp<fvScalarMatrix>
          (
              new fvScalarMatrix
              (
                  fvm::ddt(rho_, CO2_)
                + fvm::div(phi_, CO2_, "div(pc)")
                - fvm::laplacian(turbulence_.muEff()/Sc_, CO2_)
               ==
                  fvOptions_(rho_, CO2_)
              )
          )
        : tmp<fvScalarMatrix>
          (
              new fvScalarMatrix
              (
                  fvm::ddt(rho_, CO2_)
                + mvConvection->fvmDiv(phi_, CO2_)
                - fvm::laplacian(turbulence_.muEff()/Sc_, CO2_)
               ==
                  fvOptions_(rho_, CO2_)
              )
          )
    );
    fvScalarMatrix& CO2Eqn = tCO2Eqn.ref();
    CO2Eqn.relax(1.0);
    fvOptions_.constrain(CO2Eqn);
    CO2Res_ = CO2Eqn & CO2_;
    normalizeResiduals(CO2Res);

    tmp<fvScalarMatrix> tH2OEqn
    (
        isPC
        ? tmp<fvScalarMatrix>
          (
              new fvScalarMatrix
              (
                  fvm::ddt(rho_, H2O_)
                + fvm::div(phi_, H2O_, "div(pc)")
                - fvm::laplacian(turbulence_.muEff()/Sc_, H2O_)
               ==
                  fvOptions_(rho_, H2O_)
              )
          )
        : tmp<fvScalarMatrix>
          (
              new fvScalarMatrix
              (
                  fvm::ddt(rho_, H2O_)
                + mvConvection->fvmDiv(phi_, H2O_)
                - fvm::laplacian(turbulence_.muEff()/Sc_, H2O_)
               ==
                  fvOptions_(rho_, H2O_)
              )
          )
    );
    fvScalarMatrix& H2OEqn = tH2OEqn.ref();
    H2OEqn.relax(1.0);
    fvOptions_.constrain(H2OEqn);
    H2ORes_ = H2OEqn & H2O_;
    normalizeResiduals(H2ORes);

    // CASSANDRA hEqn_aero.H with PBE heat coupling disabled.
    K_ = 0.5*magSqr(U_);
    dpdt_ = fvc::ddt(p_);

    tmp<fvScalarMatrix> thEqn
    (
        isPC
        ? tmp<fvScalarMatrix>
          (
              new fvScalarMatrix
              (
                  fvm::ddt(rho_, he_)
                + fvm::div(phi_, he_, divHEScheme)
                + fvc::ddt(rho_, K_)
                + fvc::div(phi_, K_)
                + (he_.name() == "e"
                       ? fvc::div(
                           fvc::absolute(phi_/fvc::interpolate(rho_), U_),
                           p_,
                           "div(phiv,p)")
                       : -dpdt_)
                - fvm::laplacian(turbulence_.alphaEff(), he_)
               ==
                  fvOptions_(rho_, he_)
              )
          )
        : tmp<fvScalarMatrix>
          (
              new fvScalarMatrix
              (
                  fvm::ddt(rho_, he_)
                + mvConvection->fvmDiv(phi_, he_)
                + fvc::ddt(rho_, K_)
                + fvc::div(phi_, K_)
                + (he_.name() == "e"
                       ? fvc::div(
                           fvc::absolute(phi_/fvc::interpolate(rho_), U_),
                           p_,
                           "div(phiv,p)")
                       : -dpdt_)
                - fvm::laplacian(turbulence_.alphaEff(), he_)
               ==
                  fvOptions_(rho_, he_)
              )
          )
    );
    fvScalarMatrix& hEqn = thEqn.ref();

    if (MRF_.active())
    {
        // For now, omit this term; DAResidualTurboFoam also does not add an MRF term in the energy residual.
        // hEqn += fvc::div(MRF_.phi(), p_);
    }

    hEqn.relax(1.0);
    fvOptions_.constrain(hEqn);
    TRes_ = hEqn & he_;
    normalizeResiduals(TRes);

    // CASSANDRA pEqn.H non-transonic pressure and flux residuals.
    const volScalarField rAU(1.0/UEqn.A());
    const surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho_*rAU));
    const volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U_, p_));

    // Mirrors CASSANDRA pEqn.H: the ddtCorr term belongs to phiHbyA (MRF's
    // zeroFilter is the identity on the empty MRF zone list of contrailFoam).
    // At construction time the term is zero because the old-time fields have
    // not been stored yet.
    const surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rho_*HbyA)
          + rhorAUf*fvc::ddtCorr(rho_, U_, phi_)
        )
    );

    constrainPressure(p_, rho_, U_, phiHbyA, rhorAUf, MRF_);

    fvScalarMatrix pEqn
    (
        fvm::ddt(psi_, p_)
      + fvc::div(phiHbyA)
      - fvm::laplacian(rhorAUf, p_)
     ==
        fvOptions_(psi_, p_, rho_.name())
    );

    pRes_ = pEqn & p_;
    normalizeResiduals(pRes);

    phiRes_ = phiHbyA + pEqn.flux() - phi_;
    normalizePhiResiduals(phiRes);
}

void DAResidualContrailFoam::updateIntermediateVariables()
{
    // Match the standalone CASSANDRA inert-species reconstruction. N2 is a
    // derived field and is not an independent DAFoam state.
    Y_[inertIndex_] = scalar(1) - O2_ - CO2_ - H2O_;
    Y_[inertIndex_].clamp_min(0);

    // Set the energy field from the independent pressure, temperature, and
    // multicomponent composition fields. The runtime thermo model then updates
    // psi, molecular transport, thermal transport, and temperature boundaries.
    he_ = thermo_.he(p_, T_);
    thermo_.correct();

    // psiReactionThermo returns rho from the current pressure and psi fields.
    rho_ = thermo_.rho();
    K_ = 0.5*magSqr(U_);
    dpdt_ = fvc::ddt(p_);
}

void DAResidualContrailFoam::correctBoundaryConditions()
{
    /* 
    Description:
        Update the boundary condition for all the states in the selected solver
    */    
    MRF_.correctBoundaryVelocity(U_);

    U_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    T_.correctBoundaryConditions();
    O2_.correctBoundaryConditions();
    CO2_.correctBoundaryConditions();
    H2O_.correctBoundaryConditions();

    // Preserve the standalone algebraic reconstruction on boundaries as well;
    // do not independently correct N2 afterward because that can violate the
    // mass-fraction closure imposed here.
    Y_[inertIndex_] = scalar(1) - O2_ - CO2_ - H2O_;
    Y_[inertIndex_].clamp_min(0);
    Y_[inertIndex_].correctBoundaryConditions();
}

void DAResidualContrailFoam::calcPCMatWithFvMatrix(Mat)
{
    /* 
    Description:
        Calculate the diagonal block of the preconditioner matrix dRdWTPC using the fvMatrix
    */    
    FatalErrorInFunction
        << "DAResidualContrailFoam preconditioner assembly is not implemented yet."
        << exit(FatalError);
}

} // End namespace Foam

// ************************************************************************* //
