/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v5

    Description:
        Residual class for the reduced, PBE-disabled, multicomponent
        gas-phase DAContrailFoam solver.

\*---------------------------------------------------------------------------*/

#ifndef DAResidualContrailFoam_H
#define DAResidualContrailFoam_H

#include "DAResidual.H"
#include "addToRunTimeSelectionTable.H"
#include "psiReactionThermo.H"
#include "pimpleControl.H"
#include "multivariateSurfaceInterpolationScheme.H"
#include "constrainPressure.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "turbulentFluidThermoModel.H"
#include "IOMRFZoneListDF.H"

namespace Foam
{

class DAResidualContrailFoam
    : public DAResidual
{
protected:
    // Independent DAFoam flow states and residual fields
    volVectorField& U_;
    volVectorField URes_;

    volScalarField& p_;
    volScalarField pRes_;

    volScalarField& T_;
    volScalarField TRes_;

    surfaceScalarField& phi_;
    surfaceScalarField phiRes_;

    // Independent active-species states and residual fields.
    // N2 is reconstructed as the inert species and is not independent.
    volScalarField& O2_;
    volScalarField O2Res_;

    volScalarField& CO2_;
    volScalarField CO2Res_;

    volScalarField& H2O_;
    volScalarField H2ORes_;

    // Actual OpenFOAM multicomponent thermodynamic model and composition.
    psiReactionThermo& thermo_;
    basicSpecieMixture& composition_;
    PtrList<volScalarField>& Y_;

    // Derived thermodynamic and flow fields.
    volScalarField& he_;
    volScalarField& rho_;
    volScalarField& alphat_;
    volScalarField& psi_;
    volScalarField& dpdt_;
    volScalarField& K_;

    // Turbulence-model interface. The selected dynamicKEqn model owns k/kRes.
    DATurbulenceModel& daTurb_;

    // Actual OpenFOAM turbulence model used by the CASSANDRA equations.
    compressible::turbulenceModel& turbulence_;

    // Configured finite-volume options.
    fv::options& fvOptions_;

    // DAFoam differentiable MRF wrapper.
    IOMRFZoneListDF& MRF_;

    // Pressure-correction control used for residual evaluation.
    pimpleControl pimple_;

    // Inert-species reconstruction and CASSANDRA species diffusivity.
    word inertSpecie_;
    label inertIndex_;
    dimensionedScalar Sc_;

public:
    TypeName("DAContrailFoam");

    // Constructors
    
    //- Construct from components
    DAResidualContrailFoam(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex);

    //- Destructor
    virtual ~DAResidualContrailFoam()
    {
    }

    // Members

    /// clear the members
    virtual void clear();

    /// compute residual
    virtual void calcResiduals(const dictionary& options);

    /// update any intermediate variables that are dependent on state variables and are used in calcResiduals
    virtual void updateIntermediateVariables();

    /// update the boundary condition for all the states in the selected solver
    virtual void correctBoundaryConditions();

    virtual void calcPCMatWithFvMatrix(Mat PCMat);
};

} // End namespace Foam

#endif

// ************************************************************************* //


