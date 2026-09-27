/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v5

    Description:
    Child class for the dynamicKEqn LES turbulence model

    This class is adapted from the legacy OpenFOAM-v2506 dynamicKEqn
    implementation and is intentionally kept minimal for the DAFoam
    turbulence-prerequisite phase.

\*---------------------------------------------------------------------------*/

#ifndef DADynamicKEqn_H
#define DADynamicKEqn_H

#include "DATurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"
#include "LESfilter.H"
#include "simpleFilter.H"

namespace Foam
{

class DADynamicKEqn
    : public DATurbulenceModel
{
protected:
    // Existing OpenFOAM dynamicKEqn SGS-energy field; do not create a copy.
    volScalarField& k_;

    // DAFoam residual corresponding to the independent k state.
    volScalarField kRes_;

    // In primal mode, calcResiduals solves kEqn. In residual mode, it only
    // evaluates kRes_.
    label solveTurbState_ = 0;

    // Actual runtime LES model and its configured filter-width field.
    const compressible::LESModel& lesModel_;
    const volScalarField& delta_;

    // dynamicKEqn filters, matching OpenFOAM v2506 dynamicKEqn.
    simpleFilter simpleFilter_;
    autoPtr<LESfilter> filterPtr_;
    LESfilter& filter_;

    volScalarField Ck(
        const volSymmTensorField& D,
        const volScalarField& KK) const;

    volScalarField Ce(
        const volSymmTensorField& D,
        const volScalarField& KK) const;

    volScalarField Ce() const;

    void correctNut(
        const volSymmTensorField& D,
        const volScalarField& KK);

    tmp<fvScalarMatrix> kSource() const;

public:
    TypeName("dynamicKEqn");

    DADynamicKEqn(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption);

    virtual ~DADynamicKEqn()
    {
    }

    virtual void correctModelStates(wordList& modelStates) const;
    virtual void correctNut();
    virtual void correctBoundaryConditions();
    virtual void updateIntermediateVariables();
    virtual void correctStateResidualModelCon(List<List<word>>& stateCon) const;
    virtual void addModelResidualCon(HashTable<List<List<word>>>& allCon) const;
    virtual void calcResiduals(const dictionary& options);
    virtual void correct(label printToScreen);
};

} // End namespace Foam

#endif

// ************************************************************************* //
