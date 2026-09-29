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

void DAResidualContrailFoam::calcResiduals(const dictionary&)
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
    FatalErrorInFunction
        << "DAResidualContrailFoam residual evaluation is not implemented yet."
        << exit(FatalError);
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
