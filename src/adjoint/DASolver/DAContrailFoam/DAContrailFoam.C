/*---------------------------------------------------------------------------*\
    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v5
\*---------------------------------------------------------------------------*/

#include "DAContrailFoam.H"

namespace Foam
{

defineTypeNameAndDebug(DAContrailFoam, 0);
addToRunTimeSelectionTable(DASolver, DAContrailFoam, dictionary);

DAContrailFoam::DAContrailFoam(char* argsAll, PyObject* pyOptions)
    : DASolver(argsAll, pyOptions),
      pimplePtr_(nullptr),
      pThermoPtr_(nullptr),
      rhoPtr_(nullptr),
      UPtr_(nullptr),
      phiPtr_(nullptr),
      dpdtPtr_(nullptr),
      KPtr_(nullptr),
      turbulencePtr_(nullptr),
      daTurbulenceModelPtr_(nullptr),
      pressureControlPtr_(nullptr),
      MRFPtr_(nullptr)
{
}

void DAContrailFoam::initSolver()
{
    Info << "Initializing DAContrailFoam" << endl;
    Time& runTime = runTimePtr_();
    fvMesh& mesh = meshPtr_();
    argList& args = argsPtr_();

    #include "createPimpleControlPython.H"
    #include "createFieldsContrail.H"

    // read the active turbulence model name from constant/turbulenceProperties
    // (LES or RAS, depending on the simulationType)
    const IOdictionary turbProps(
        IOobject(
            "turbulenceProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false));
    const word turbModelName =
        turbProps.found("LES")
            ? turbProps.subDict("LES").getWord("LESModel")
            : turbProps.subDict("RAS").getWord("RASModel");
    Info << "Selecting turbulence model: " << turbModelName << endl;
    daTurbulenceModelPtr_.reset(DATurbulenceModel::New(turbModelName, mesh, daOptionPtr_()));

    #include "createAdjoint.H"
}

label DAContrailFoam::solvePrimal()
{
    #include "createRefsContrail.H"

    Info << "Solving DAContrailFoam primal" << endl;

    // Declare UEqn here so it is visible to both UEqnContrail.H and pEqnContrail.H
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      - fvm::laplacian(turbulence.muEff(), U)
     ==
        fvOptions(rho, U)
    );

    while (pimple.loop())
    {
        #include "UEqnContrail.H"

        // Update auxiliary fields
        K = 0.5*magSqr(U);
        dpdt = fvc::ddt(p);

        #include "hEqnContrail.H"
        #include "YEqnContrail.H"

        while (pimple.correct())
        {
            #include "pEqnContrail.H"
        }
    }

    return 0;
}

} // End namespace Foam
