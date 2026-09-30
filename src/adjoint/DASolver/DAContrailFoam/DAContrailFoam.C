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

    // reduceIO does not write mesh, but if there is a FFD variable, set writeMesh to 1
    dictionary dvSubDict = daOptionPtr_->getAllOptions().subDict("inputInfo");
    forAll(dvSubDict.toc(), idxI)
    {
        word dvName = dvSubDict.toc()[idxI];
        if (dvSubDict.subDict(dvName).getWord("type") == "volCoord")
        {
            reduceIOWriteMesh_ = 1;
            break;
        }
    }
}

label DAContrailFoam::solvePrimal()
{
    #include "createRefsContrail.H"

    // equivalent of CASSANDRA's turbulence->validate()
    daTurbulenceModelPtr_->updateIntermediateVariables();

    // DAResidualContrailFoam only evaluates the non-transonic, non-consistent
    // pressure residual, so reject the two unsupported branches up front.
    if (pimple.transonic())
    {
        FatalErrorInFunction
            << "DAContrailFoam does not support transonic pressure yet."
            << exit(FatalError);
    }
    if (pimple.consistent())
    {
        FatalErrorInFunction
            << "DAContrailFoam does not support the consistent (pcEqn) pressure "
            << "corrector yet."
            << exit(FatalError);
    }

    Info << "Solving DAContrailFoam primal" << endl;
    Info << "\nStarting time loop\n" << endl;

    label pimplePrintToScreen = 0;

    // we need to reduce the number of files written to the disk to minimize the file IO load
    label reduceIO = daOptionPtr_->getAllOptions().subDict("unsteadyAdjoint").getLabel("reduceIO");
    wordList additionalOutput;
    if (reduceIO)
    {
        daOptionPtr_->getAllOptions().subDict("unsteadyAdjoint").readEntry<wordList>("additionalOutput", additionalOutput);
    }

    #include "createTimeControls.H"

    label regModelFail = 0;
    label fail = 0;

    while (runTime.run())
    {
        if (primalFuncStdTol_ > 0)
        {
            this->calcFuncStd();
            this->calcFuncSlope();
        }

        #include "readTimeControls.H"

        // Store previous values
        p.storePrevIter();
        rho.storePrevIter();

        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;
        Info << "Time = " << runTime.timeName() << nl << endl;

        printToScreen_ = this->isPrintTime(runTime, printIntervalUnsteady_);

        // if we have unsteadyField in inputInfo, assign GlobalVar::inputFieldUnsteady to OF fields at each time step
        this->updateInputFieldUnsteady();

        // explicit continuity update before the PIMPLE loop
        #include "rhoEqnContrail.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            pimplePrintToScreen = pimple.finalIter() && printToScreen_ ? 1 : 0;

            #include "UEqnContrail.H"
            #include "YEqnContrail.H"
            #include "hEqnContrail.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqnContrail.H"
            }

            if (pimple.turbCorr())
            {
                daTurbulenceModelPtr_->correct(pimplePrintToScreen);
            }

            // update the output field value at each iteration, if the regression model is active
            fail = daRegressionPtr_->compute();
        }

        regModelFail += fail;

        // CASSANDRA contrailFoam.C: rho = thermo.rho() after the PIMPLE loop
        rho = thermo.rho();

        if (this->validateStates())
        {
            // write data to files and quit
            runTime.writeNow();
            mesh.write();
            return 1;
        }

        this->calcAllFunctions(printToScreen_);
        daRegressionPtr_->printInputInfo(printToScreen_);
        daTurbulenceModelPtr_->printYPlus(printToScreen_);
        this->printElapsedTime(runTime, printToScreen_);

        // with adaptive time stepping the number of steps is not known up front,
        // so detect the final step the same way Time::run() does
        const label isLastStep = !(
            runTime.value() + runTime.deltaTValue()
                < runTime.endTime().value() - 0.5*runTime.deltaTValue()
        );

        if (reduceIO && !isLastStep)
        {
            this->writeAdjStates(reduceIOWriteMesh_, additionalOutput);
            daRegressionPtr_->writeFeatures();
        }
        else
        {
            runTime.write();
            daRegressionPtr_->writeFeatures();
        }
    }

    if (regModelFail != 0)
    {
        return 1;
    }

    // need to save primalFinalTimeIndex_.
    primalFinalTimeIndex_ = runTime.timeIndex();

    // write the mesh to files
    mesh.write();

    Info << "End\n"
         << endl;

    return 0;
}

} // End namespace Foam

// ************************************************************************* //
