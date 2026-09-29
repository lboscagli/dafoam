// Field creation for DAContrailFoam
// This file is included in DAContrailFoam::initSolver()
// Mirrors CASSANDRA contrailFoam/createFields.H (gas-phase, PBE disabled)
// with the DAFoam ownership pattern of DARhoPimpleFoam/createFieldsRhoPimple.H

// 1. Thermodynamic model: creates and registers p, T, psi and the species Y
pThermoPtr_.reset(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermoPtr_();
thermo.validate(args.executable(), "h", "e");

// Inert specie validation (CASSANDRA createFields.H)
basicSpecieMixture& composition = thermo.composition();
const word inertSpecie(thermo.get<word>("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie
        << " not found in available species "
        << composition.species()
        << exit(FatalIOError);
}

volScalarField& p = thermo.p();

// 2. Primary flow fields
rhoPtr_.reset(new volScalarField(
    IOobject(
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    thermo.rho()));
volScalarField& rho = rhoPtr_();

UPtr_.reset(new volVectorField(
    IOobject(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh));
volVectorField& U = UPtr_();

#include "compressibleCreatePhiPython.H"

mesh.setFluxRequired(p.name());

// 3. Pressure control (CASSANDRA createFields.H)
pressureControlPtr_.reset(new pressureControl(p, rho, pimple.dict(), false));

// 4. Compressible turbulence model (registers alphat/nut etc.)
Info << "Creating turbulence model" << nl << endl;
turbulencePtr_.reset(compressible::turbulenceModel::New(rho, U, phi, thermo));

// 5. Derived fields
dpdtPtr_.reset(new volScalarField(
    IOobject(
        "dpdt",
        runTime.timeName(),
        mesh),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0.0)));

KPtr_.reset(new volScalarField("K", 0.5*magSqr(U)));

// 6. MRF zone list (registers MRFProperties; reads the dict only if present)
MRFPtr_.reset(new IOMRFZoneListDF(mesh));
