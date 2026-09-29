// Reference setup for DAContrailFoam
// This file is included at the start of solvePrimal() to provide local references

// Basic DAFoam references
fvMesh& mesh = meshPtr_();
Time& runTime = runTimePtr_();
pimpleControlDF& pimple = pimplePtr_();

// Thermodynamics (p and T are owned and registered by the thermo)
psiReactionThermo& thermo = pThermoPtr_();
volScalarField& p = thermo.p();
volScalarField& psi = const_cast<volScalarField&>(mesh.thisDb().lookupObject<volScalarField>("thermo:psi"));
basicSpecieMixture& composition = const_cast<basicSpecieMixture&>(thermo.composition());
PtrList<volScalarField>& Y = composition.Y();
const label inertIndex = composition.species().find(thermo.get<word>("inertSpecie"));

// Dereference flow fields
volScalarField& rho = *rhoPtr_;
volVectorField& U = *UPtr_;
surfaceScalarField& phi = *phiPtr_;
volScalarField& K = *KPtr_;
volScalarField& dpdt = *dpdtPtr_;

// Turbulence
compressible::turbulenceModel& turbulence = *turbulencePtr_;

// Constants
scalar Sc = 0.7;

// Options
fv::options& fvOptions = fv::options::New(mesh);

// MRF
IOMRFZoneListDF& MRF = const_cast<IOMRFZoneListDF&>(mesh.thisDb().lookupObject<IOMRFZoneListDF>("MRFProperties"));

// Multivariate Convection Setup
multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
forAll(Y, speciesI)
{
    fields.add(Y[speciesI]);
}
fields.add(thermo.he());

tmp<fv::convectionScheme<scalar>> mvConvection =
    fv::convectionScheme<scalar>::New(
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    );
