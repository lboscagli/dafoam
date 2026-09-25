/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v5

\*---------------------------------------------------------------------------*/

#include "DADynamicKEqn.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DADynamicKEqn, 0);
addToRunTimeSelectionTable(DATurbulenceModel, DADynamicKEqn, dictionary);

volScalarField DADynamicKEqn::Ck(
    const volSymmTensorField& D,
    const volScalarField& KK) const
{
    volScalarField deltaField(
        IOobject(
            "delta",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh_,
        dimensionedScalar("delta", dimless, 0.0));
    deltaField.primitiveFieldRef() = cbrt(mesh_.V());

    const volSymmTensorField LL(
        simpleFilter_(dev(filter_(sqr(U_)) - sqr(filter_(U_)))));

    const volSymmTensorField MM(
        simpleFilter_(
            -2.0 * deltaField * sqrt(max(KK, dimensionedScalar(KK.dimensions(), Zero)))
                * filter_(D)));

    const volScalarField CkField(
        simpleFilter_(0.5 * (LL && MM))
            / (simpleFilter_(magSqr(MM))
               + dimensionedScalar("small", sqr(MM.dimensions()), VSMALL)));

    tmp<volScalarField> tfld = 0.5 * (mag(CkField) + CkField);
    return tfld();
}

volScalarField DADynamicKEqn::Ce(
    const volSymmTensorField& D,
    const volScalarField& KK) const
{
    volScalarField deltaField(
        IOobject(
            "delta",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh_,
        dimensionedScalar("delta", dimless, 0.0));
    deltaField.primitiveFieldRef() = cbrt(mesh_.V());

    const volScalarField CeField(
        simpleFilter_(nuEff() * (filter_(magSqr(D)) - magSqr(filter_(D))))
            / simpleFilter_(pow(KK, 1.5) / (2.0 * deltaField)));

    tmp<volScalarField> tfld = 0.5 * (mag(CeField) + CeField);
    return tfld();
}

volScalarField DADynamicKEqn::Ce() const
{
    const volSymmTensorField D(devSymm(fvc::grad(U_)));

    volScalarField KK(
        0.5 * (filter_(magSqr(U_)) - magSqr(filter_(U_))));
    KK.clamp_min(SMALL);

    return Ce(D, KK);
}

void DADynamicKEqn::correctNut(
    const volSymmTensorField& D,
    const volScalarField& KK)
{
    volScalarField deltaField(
        IOobject(
            "delta",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh_,
        dimensionedScalar("delta", dimless, 0.0));
    deltaField.primitiveFieldRef() = cbrt(mesh_.V());

    nut_ = Ck(D, KK) * sqrt(max(k_, dimensionedScalar(k_.dimensions(), SMALL)))
        * deltaField;
    nut_.correctBoundaryConditions();

    this->correctAlphat();
}

DADynamicKEqn::DADynamicKEqn(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption)
    : DATurbulenceModel(modelType, mesh, daOption),
      k_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("k"))),
      simpleFilter_(mesh),
      filterPtr_(LESfilter::New(mesh, coeffDict_)),
      filter_(filterPtr_())
{
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
    const volScalarField KK(
        0.5 * (filter_(magSqr(U_)) - magSqr(filter_(U_))));

    correctNut(symm(fvc::grad(U_)), KK);
}

void DADynamicKEqn::correctBoundaryConditions()
{
    k_.correctBoundaryConditions();
    nut_.correctBoundaryConditions();

    if (mesh_.thisDb().foundObject<volScalarField>("alphat"))
    {
        this->correctAlphat();
    }
}

void DADynamicKEqn::updateIntermediateVariables()
{
    this->correctNut();
}

void DADynamicKEqn::correctStateResidualModelCon(List<List<word>>& stateCon) const
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

void DADynamicKEqn::addModelResidualCon(HashTable<List<List<word>>>& allCon) const
{
    word pName = "p";
    if (!mesh_.thisDb().foundObject<volScalarField>("p"))
    {
        pName = "p_rgh";
    }

    if (turbModelType_ == "incompressible")
    {
        allCon.set(
            "kRes",
            {
                {"U", "k", "phi"},
                {"U", "k"},
                {"U", "k"}
            });
    }
    else
    {
        allCon.set(
            "kRes",
            {
                {"U", "T", pName, "k", "phi"},
                {"U", "T", pName, "k"},
                {"U", "T", pName, "k"}
            });
    }
}

void DADynamicKEqn::correct(label printToScreen)
{
    this->correctNut();
}

void DADynamicKEqn::calcResiduals(const dictionary& options)
{
    // The minimal prerequisite implementation keeps the dynamicKEqn residual path
    // unimplemented until the actual LES turbulence model is fully wired into the
    // DAFoam solver flow. This registration is enough to expose the model and allow
    // DAFoam to select it in original mode.
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
