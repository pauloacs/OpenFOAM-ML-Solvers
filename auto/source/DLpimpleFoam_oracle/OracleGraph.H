#ifndef OracleGraph_H
#define OracleGraph_H

#include "fvCFD.H"

namespace Foam
{

labelList oracleGraphCoarsen
(
    const polyMesh& mesh,
    const label levels,
    label& coarseCount,
    const word& obstaclePatch,
    const label obstacleProtectionLayers,
    label& protectedCellCount
);

template<class Type>
void oracleGraphRoundTrip
(
    Field<Type>& values,
    const labelList& mapping,
    const scalarField& volumes,
    const label coarseCount
)
{
    if (values.size() != mapping.size() || volumes.size() != mapping.size())
    {
        FatalErrorInFunction
            << "Graph round-trip arrays must have matching native sizes."
            << exit(FatalError);
    }

    Field<Type> coarseValues(coarseCount, pTraits<Type>::zero);
    scalarField coarseVolumes(coarseCount, 0.0);
    forAll(values, celli)
    {
        const label aggregate = mapping[celli];
        coarseValues[aggregate] += volumes[celli]*values[celli];
        coarseVolumes[aggregate] += volumes[celli];
    }

    forAll(coarseValues, aggregate)
    {
        if (coarseVolumes[aggregate] <= VSMALL)
        {
            FatalErrorInFunction
                << "Graph aggregate " << aggregate << " has no volume."
                << exit(FatalError);
        }
        coarseValues[aggregate] /= coarseVolumes[aggregate];
    }

    forAll(values, celli)
    {
        values[celli] = coarseValues[mapping[celli]];
    }
}

class OracleGraphLinearTransfer
{
    labelList mapping_;
    scalarField volumes_;
    vectorField centres_;
    label coarseCount_;
    scalarField coarseVolumes_;
    vectorField coarseCentres_;
    labelList edgeSource_;
    labelList edgeTarget_;
    vectorField edgeOffsets_;
    scalarField inverseNormal_;
    labelList interfaceOwner_;
    labelList interfaceNeighbour_;
    vectorField interfaceOffsets_;
    scalarField interfaceWeights_;
    scalarField interfaceStep_;

public:
    OracleGraphLinearTransfer
    (
        const fvMesh& mesh,
        const labelList& mapping,
        const label coarseCount
    );

    void linearRoundTrip(Field<scalar>& values) const;
    void linearRoundTrip(Field<vector>& values) const;
};

}

#endif