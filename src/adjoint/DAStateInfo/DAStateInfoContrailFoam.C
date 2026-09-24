/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v5

\*---------------------------------------------------------------------------*/

#include "DAStateInfoContrailFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAStateInfoContrailFoam, 0);
addToRunTimeSelectionTable(DAStateInfo, DAStateInfoContrailFoam, dictionary);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

DAStateInfoContrailFoam::DAStateInfoContrailFoam(
    const word modelType,
    const fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel)
    : DAStateInfo(modelType, mesh, daOption, daModel)
{
    /*
    Description:
        Register the names of state variables.
        N2 is reconstructed as the inert species, so it is not included in the
        active gas-phase states.
    */

    stateInfo_["volScalarStates"].append("p");
    stateInfo_["volScalarStates"].append("T");
    stateInfo_["volScalarStates"].append("O2");
    stateInfo_["volScalarStates"].append("CO2");
    stateInfo_["volScalarStates"].append("H2O");
    stateInfo_["volVectorStates"].append("U");
    stateInfo_["surfaceScalarStates"].append("phi");
    stateInfo_["modelStates"].append("nut");

    daModel.correctModelStates(stateInfo_["modelStates"]);

    /*
    Description:
        Conservative connectivity for the reduced gas-phase-only contrailFoam model.
        We initialize all gas-phase residuals at level zero against the full set of
        active gas-phase state variables, the turbulence viscosity, and the flux.
    */

    stateResConInfo_.set(
        "URes",
        {
            {"U", "p", "T", "O2", "CO2", "H2O", "nut", "phi"}, // lv0
            {"U", "p", "T", "O2", "CO2", "H2O", "nut"}, // lv1
            {"U", "p", "T", "O2", "CO2", "H2O"} // lv2
        });

    stateResConInfo_.set(
        "TRes",
        {
            {"U", "p", "T", "O2", "CO2", "H2O", "nut", "phi"}, // lv0
            {"U", "p", "T", "O2", "CO2", "H2O", "nut"}, // lv1
            {"U", "p", "T", "O2", "CO2", "H2O"} // lv2
        });

    stateResConInfo_.set(
        "pRes",
        {
            {"U", "p", "T", "O2", "CO2", "H2O", "nut", "phi"}, // lv0
            {"U", "p", "T", "O2", "CO2", "H2O", "nut", "phi"}, // lv1
            {"U", "p", "T", "O2", "CO2", "H2O", "nut"}, // lv2
            {"U"} // lv3
        });

    stateResConInfo_.set(
        "phiRes",
        {
            {"U", "p", "T", "O2", "CO2", "H2O", "nut", "phi"}, // lv0
            {"U", "p", "T", "O2", "CO2", "H2O", "nut"}, // lv1
            {"U", "p", "T", "O2", "CO2", "H2O"} // lv2
        });

    stateResConInfo_.set(
        "O2Res",
        {
            {"U", "p", "T", "O2", "CO2", "H2O", "nut", "phi"}, // lv0
            {"U", "p", "T", "O2", "CO2", "H2O", "nut"}, // lv1
            {"U", "p", "T", "O2", "CO2", "H2O"} // lv2
        });

    stateResConInfo_.set(
        "CO2Res",
        {
            {"U", "p", "T", "O2", "CO2", "H2O", "nut", "phi"}, // lv0
            {"U", "p", "T", "O2", "CO2", "H2O", "nut"}, // lv1
            {"U", "p", "T", "O2", "CO2", "H2O"} // lv2
        });

    stateResConInfo_.set(
        "H2ORes",
        {
            {"U", "p", "T", "O2", "CO2", "H2O", "nut", "phi"}, // lv0
            {"U", "p", "T", "O2", "CO2", "H2O", "nut"}, // lv1
            {"U", "p", "T", "O2", "CO2", "H2O"} // lv2
        });

    daModel.correctStateResidualModelCon(stateResConInfo_["URes"]);
    daModel.correctStateResidualModelCon(stateResConInfo_["TRes"]);
    daModel.correctStateResidualModelCon(stateResConInfo_["pRes"]);
    daModel.correctStateResidualModelCon(stateResConInfo_["phiRes"]);
    daModel.correctStateResidualModelCon(stateResConInfo_["O2Res"]);
    daModel.correctStateResidualModelCon(stateResConInfo_["CO2Res"]);
    daModel.correctStateResidualModelCon(stateResConInfo_["H2ORes"]);

    daModel.addModelResidualCon(stateResConInfo_);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
