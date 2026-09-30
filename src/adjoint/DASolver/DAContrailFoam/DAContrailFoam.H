/*---------------------------------------------------------------------------*\
    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v5

    Description:
        Skeleton implementation of DAContrailFoam
\*---------------------------------------------------------------------------*/

#ifndef DAContrailFoam_H
#define DAContrailFoam_H

#include "DASolver.H"
#include "addToRunTimeSelectionTable.H"
#include "psiReactionThermo.H"
#include "pimpleControlDF.H"
#include "pressureControl.H"

namespace Foam
{

class DAContrailFoam
    : public DASolver
{
protected:
    /// simple pointer
    autoPtr<pimpleControlDF> pimplePtr_;

    /// thermal model (owns the registered p and T fields)
    autoPtr<psiReactionThermo> pThermoPtr_;

    /// density field pointer
    autoPtr<volScalarField> rhoPtr_;

    /// velocity field pointer
    autoPtr<volVectorField> UPtr_;

    /// surface flux field pointer
    autoPtr<surfaceScalarField> phiPtr_;

    /// pressure time deriv pointer
    autoPtr<volScalarField> dpdtPtr_;

    /// kinetic energy pointer
    autoPtr<volScalarField> KPtr_;

    /// turbulence pointer
    autoPtr<compressible::turbulenceModel> turbulencePtr_;

    /// DATurbulenceModel pointer
    autoPtr<DATurbulenceModel> daTurbulenceModelPtr_;

    /// pressure control pointer
    autoPtr<pressureControl> pressureControlPtr_;

    /// MRF zone list pointer
    autoPtr<IOMRFZoneListDF> MRFPtr_;

    /// write the mesh together with the reduceIO adjoint states
    label reduceIOWriteMesh_ = 0;

public:
    TypeName("DAContrailFoam");


    // Constructors
    DAContrailFoam(char* argsAll, PyObject* pyOptions);
    virtual ~DAContrailFoam() {}

    virtual void initSolver();
    virtual label solvePrimal();
};

} // End namespace Foam

#endif
