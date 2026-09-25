/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v5

    Description:
    Child class for the dynamicKEqn LES turbulence model

    This class is adapted from the legacy OpenFOAM-v2512 dynamicKEqn
    implementation and is intentionally kept minimal for the DAFoam
    turbulence-prerequisite phase.

\*---------------------------------------------------------------------------*/

#ifndef DADynamicKEqn_H
#define DADynamicKEqn_H

#include "DATurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "LESfilter.H"
#include "simpleFilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DADynamicKEqn Declaration
\*---------------------------------------------------------------------------*/

class DADynamicKEqn
    : public DATurbulenceModel
{

protected:
    /// SGS kinetic energy field
    volScalarField& k_;

    /// legacy dynamicKEqn filter objects
    simpleFilter simpleFilter_;
    autoPtr<LESfilter> filterPtr_;
    LESfilter& filter_;

    /// dynamicKEqn coefficient functions adapted from OpenFOAM-v2512
    volScalarField Ck(
        const volSymmTensorField& D,
        const volScalarField& KK) const;

    volScalarField Ce(
        const volSymmTensorField& D,
        const volScalarField& KK) const;

    volScalarField Ce() const;

    /// update the SGS viscosity from the dynamic model
    void correctNut(
        const volSymmTensorField& D,
        const volScalarField& KK);

public:
    TypeName("dynamicKEqn");
    // Constructors

    //- Construct from components
    DADynamicKEqn(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption);

    //- Destructor
    virtual ~DADynamicKEqn()
    {
    }

    // Member functions

    /// update the turbulence state for DAStateInfo::regStates_
    virtual void correctModelStates(wordList& modelStates) const;

    /// update nut based on other turbulence variables and update the BCs
    virtual void correctNut();

    /// update turbulence variable boundary values
    virtual void correctBoundaryConditions();

    /// update any intermediate variables that are dependent on state variables and are used in calcResiduals
    virtual void updateIntermediateVariables();

    /// update the original variable connectivity for the adjoint state residuals in stateCon
    virtual void correctStateResidualModelCon(List<List<word>>& stateCon) const;

    /// add the model residual connectivity to stateCon
    virtual void addModelResidualCon(HashTable<List<List<word>>>& allCon) const;

    /// compute the turbulence residuals
    virtual void calcResiduals(const dictionary& options);

    /// solve the residual equations and update the state
    virtual void correct(label printToScreen);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
