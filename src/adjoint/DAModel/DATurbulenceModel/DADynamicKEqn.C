/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v5

\*---------------------------------------------------------------------------*/

#include "DADynamicKEqn.H"
#include "bound.H"
#include "fvOptions.H"

namespace Foam
{

defineTypeNameAndDebug(DADynamicKEqn, 0);
addToRunTimeSelectionTable(DATurbulenceModel, DADynamicKEqn, dictionary);

volScalarField DADynamicKEqn::Ck(
    const volSymmTensorField& D,
    const volScalarField& KK) const
{
    const volSymmTensorField LL
    (
        simpleFilter_(dev(filter_(sqr(U_)) - sqr(filter_(U_))))
    );

    const volSymmTensorField MM
    (
        simpleFilter_
        (
            -2.0*delta_
           *sqrt(max(KK, dimensionedScalar(KK.dimensions(), Zero)))
           *filter_(D)
        )
    );

    const volScalarField CkField
    (
        simpleFilter_(0.5*(LL && MM))
       /
        (
            simpleFilter_(magSqr(MM))
          + dimensionedScalar("small", sqr(MM.dimensions()), VSMALL)
        )
    );

    tmp<volScalarField> tfld = 0.5*(mag(CkField) + CkField);
    return tfld();
}

volScalarField DADynamicKEqn::Ce(
    const volSymmTensorField& D,
    const volScalarField& KK) const
{
    const volScalarField CeField
    (
        simpleFilter_
        (
            nuEff()*(filter_(magSqr(D)) - magSqr(filter_(D)))
        )
       /simpleFilter_(pow(KK, 1.5)/(2.0*delta_))
    );

    tmp<volScalarField> tfld = 0.5*(mag(CeField) + CeField);
    return tfld();
}

volScalarField DADynamicKEqn::Ce() const
{
    const volSymmTensorField D(devSymm(fvc::grad(U_)));

    volScalarField KK
    (
        0.5*(filter_(magSqr(U_)) - magSqr(filter_(U_)))
    );
    KK.clamp_min(SMALL);

    return Ce(D, KK);
}

void DADynamicKEqn::correctNut(
    const volSymmTensorField& D,
    const volScalarField& KK)
{
    // Matches OpenFOAM v2506 dynamicKEqn::correctNut().
    nut_ = Ck(D, KK)*sqrt(k_)*delta_;
    nut_.correctBoundaryConditions();

    fv::options::New(mesh_).correct(nut_);

    // DATurbulenceModel::correctAlphat must use the OpenFOAM-compatible LES
    // default Prt = 1.0 when the LES dictionary omits Prt.
    this->correctAlphat();
}

tmp<fvScalarMatrix> DADynamicKEqn::kSource() const
{
    return tmp<fvScalarMatrix>::New
    (
        k_,
        dimVolume*this->rhoDimensions()*k_.dimensions()/dimTime
    );
}

DADynamicKEqn::DADynamicKEqn(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption)
    : DATurbulenceModel(modelType, mesh, daOption),
      k_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("k"))),
      kRes_
      (
          IOobject
          (
              "kRes",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          mesh,
          dimensionedScalar("kRes", dimless, 0.0),
          zeroGradientFvPatchField<scalar>::typeName
      ),
      lesModel_(refCast<const compressible::LESModel>(
          mesh.thisDb().lookupObject<compressible::turbulenceModel>(
              compressible::turbulenceModel::propertiesName))),
      delta_(lesModel_.delta()),
      simpleFilter_(mesh),
      filterPtr_(
        LESfilter::New
        (
            mesh,
            turbDict_.subDict("LES").subDict("dynamicKEqnCoeffs")
        )
      ),
      filter_(filterPtr_())
{
    if (turbModelType_ != "compressible")
    {
        FatalErrorInFunction
            << "DADynamicKEqn currently supports compressible LES only."
            << exit(FatalError);
    }

    kRes_.dimensions().reset(
        rhoDimensions()*k_.dimensions()/dimTime
    );
}

void DADynamicKEqn::correctModelStates(wordList& modelStates) const
{
    forAll(modelStates, idxI)
    {
        if (modelStates[idxI] == "nut")
        {
            modelStates[idxI] = "k";
        }
    }
}

void DADynamicKEqn::correctNut()
{
    const volScalarField KK
    (
        0.5*(filter_(magSqr(U_)) - magSqr(filter_(U_)))
    );

    correctNut(symm(fvc::grad(U_)), KK);
}

void DADynamicKEqn::correctBoundaryConditions()
{
    k_.correctBoundaryConditions();
    this->correctNut();
}

void DADynamicKEqn::updateIntermediateVariables()
{
    this->correctNut();
}

void DADynamicKEqn::correctStateResidualModelCon(
    List<List<word>>& stateCon) const
{
    forAll(stateCon, idxI)
    {
        forAll(stateCon[idxI], idxJ)
        {
            if (stateCon[idxI][idxJ] == "nut")
            {
                stateCon[idxI][idxJ] = "k";
            }
        }
    }
}

void DADynamicKEqn::addModelResidualCon(
    HashTable<List<List<word>>>& allCon) const
{
    const word pName = mesh_.thisDb().foundObject<volScalarField>("p")
        ? word("p")
        : word("p_rgh");

    allCon.set(
        "kRes",
        {
            {"U", "T", pName, "k", "phi"},
            {"U", "T", pName, "k"},
            {"U", "T", pName, "k"}
        }
    );
}

void DADynamicKEqn::correct(label printToScreen)
{
    solveTurbState_ = 1;

    dictionary options;
    options.set("printToScreen", printToScreen);
    this->calcResiduals(options);

    solveTurbState_ = 0;
}

void DADynamicKEqn::calcResiduals(const dictionary& options)
{
    word divKScheme = "div(phi,k)";
    label isPC = 0;
    const label printToScreen = options.lookupOrDefault<label>("printToScreen", 0);

    if (!solveTurbState_)
    {
        isPC = options.getLabel("isPC");
        if (isPC)
        {
            divKScheme = "div(pc)";
        }
    }

    const volScalarField rho(this->rho());

    const volScalarField divU
    (
        fvc::div(fvc::absolute(phi_/fvc::interpolate(rho), U_))
    );

    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volSymmTensorField D(devSymm(tgradU()));
    const volScalarField G(2.0*nut_*(tgradU() && D));
    tgradU.clear();

    volScalarField KK
    (
        0.5*(filter_(magSqr(U_)) - magSqr(filter_(U_)))
    );
    KK.clamp_min(SMALL);

    fv::options& fvOptions(fv::options::New(mesh_));

    tmp<fvScalarMatrix> tkEqn
    (
        fvm::ddt(phase_, rho, k_)
      + fvm::div(phaseRhoPhi_, k_, divKScheme)
      - fvm::laplacian(phase_*rho*nuEff(), k_)
     ==
        phase_*rho*G
      - fvm::SuSp((2.0/3.0)*phase_*rho*divU, k_)
      - fvm::Sp(Ce(D, KK)*phase_*rho*sqrt(k_)/delta_, k_)
      + kSource()
      + fvOptions(phase_, rho, k_)
    );

    fvScalarMatrix& kEqn = tkEqn.ref();
    kEqn.relax();
    fvOptions.constrain(kEqn);

    if (solveTurbState_)
    {
        SolverPerformance<scalar> solverK = solve(kEqn);
        DAUtility::primalResidualControl(
            solverK, printToScreen, "k", daGlobalVar_.primalMaxRes
        );

        fvOptions.correct(k_);
        bound(k_, kMin_);
        correctNut(D, KK);
    }
    else
    {
        kRes_ = kEqn & k_;
        normalizeResiduals(kRes);
    }
}

} // End namespace Foam

// ************************************************************************* //
