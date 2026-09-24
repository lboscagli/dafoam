/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v5

    Description:
        Child class for DAContrailFoam

\*---------------------------------------------------------------------------*/

#ifndef DAStateInfoContrailFoam_H
#define DAStateInfoContrailFoam_H

#include "DAStateInfo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAStateInfoContrailFoam Declaration
\*---------------------------------------------------------------------------*/

class DAStateInfoContrailFoam
    : public DAStateInfo
{

public:
    TypeName("DAContrailFoam");
    // Constructors

    //- Construct from components
    DAStateInfoContrailFoam(
        const word modelType,
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel);

    //- Destructor
    virtual ~DAStateInfoContrailFoam()
    {
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
