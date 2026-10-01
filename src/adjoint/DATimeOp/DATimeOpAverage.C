/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v5

\*---------------------------------------------------------------------------*/

#include "DATimeOpAverage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DATimeOpAverage, 0);
addToRunTimeSelectionTable(DATimeOp, DATimeOpAverage, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DATimeOpAverage::DATimeOpAverage(
    const word timeOpType,
    const dictionary options)
    : DATimeOp(timeOpType, options)
{
}

scalar DATimeOpAverage::compute(
    const scalarList& valList,
    const label iStart,
    const label iEnd)
{
    // empty window (iEnd < iStart) or out-of-range end index would otherwise
    // divide by zero (NaN) or read out of bounds
    if (iStart < 0 || iEnd < iStart || iEnd >= valList.size())
    {
        FatalErrorIn("DATimeOpAverage::compute")
            << "empty or invalid time-op window: iStart = " << iStart
            << ", iEnd = " << iEnd << ", valList size = " << valList.size()
            << ". The window must contain at least one stored sample."
            << abort(FatalError);
    }

    // return the average value from valList
    scalar avg = 0.0;
    // NOTE. We need to use <= here
    for (label i = iStart; i <= iEnd; i++)
    {
        avg += valList[i];
    }
    avg /= (iEnd - iStart + 1);
    return avg;
}

scalar DATimeOpAverage::dFScaling(
    const scalarList& valList,
    const label iStart,
    const label iEnd,
    const label timeIdx)
{
    // empty window would give 1/0 (inf) or a negative scaling
    if (iStart < 0 || iEnd < iStart)
    {
        FatalErrorIn("DATimeOpAverage::dFScaling")
            << "empty or invalid time-op window: iStart = " << iStart
            << ", iEnd = " << iEnd
            << ". The window must contain at least one stored sample."
            << abort(FatalError);
    }

    // return 1/N as the dF scaling

    scalar scaling = 1.0 / (iEnd - iStart + 1);

    return scaling;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
