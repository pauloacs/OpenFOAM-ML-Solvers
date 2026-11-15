#include "OracleGraph.H"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace Foam
{

namespace
{

void symmetricPseudoInverse3
(
    const scalar input[3][3],
    scalar inverse[3][3]
)
{
    scalar matrix[3][3];
    scalar eigenvectors[3][3] =
    {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };

    for (label row = 0; row < 3; ++row)
    {
        for (label column = 0; column < 3; ++column)
        {
            matrix[row][column] = input[row][column];
            inverse[row][column] = 0;
        }
    }

    for (label iteration = 0; iteration < 50; ++iteration)
    {
        label p = 0;
        label q = 1;
        scalar largest = mag(matrix[p][q]);
        for (label row = 0; row < 3; ++row)
        {
            for (label column = row + 1; column < 3; ++column)
            {
                if (mag(matrix[row][column]) > largest)
                {
                    p = row;
                    q = column;
                    largest = mag(matrix[row][column]);
                }
            }
        }

        scalar diagonalScale = max
        (
            mag(matrix[0][0]),
            max(mag(matrix[1][1]), mag(matrix[2][2]))
        );
        if (largest <= max(diagonalScale*1e-14, VSMALL))
        {
            break;
        }

        const scalar angle = 0.5*std::atan2
        (
            2*matrix[p][q],
            matrix[q][q] - matrix[p][p]
        );
        const scalar cosine = std::cos(angle);
        const scalar sine = std::sin(angle);
        const scalar pp = matrix[p][p];
        const scalar qq = matrix[q][q];
        const scalar pq = matrix[p][q];

        for (label k = 0; k < 3; ++k)
        {
            if (k != p && k != q)
            {
                const scalar kp = matrix[k][p];
                const scalar kq = matrix[k][q];
                matrix[k][p] = matrix[p][k] = cosine*kp - sine*kq;
                matrix[k][q] = matrix[q][k] = sine*kp + cosine*kq;
            }

            const scalar vp = eigenvectors[k][p];
            const scalar vq = eigenvectors[k][q];
            eigenvectors[k][p] = cosine*vp - sine*vq;
            eigenvectors[k][q] = sine*vp + cosine*vq;
        }

        matrix[p][p] = cosine*cosine*pp - 2*sine*cosine*pq + sine*sine*qq;
        matrix[q][q] = sine*sine*pp + 2*sine*cosine*pq + cosine*cosine*qq;
        matrix[p][q] = matrix[q][p] = 0;
    }

    const scalar maxEigenvalue = max
    (
        max(matrix[0][0], matrix[1][1]),
        matrix[2][2]
    );
    const scalar threshold = maxEigenvalue*1e-10;
    for (label eigen = 0; eigen < 3; ++eigen)
    {
        const scalar eigenvalue = matrix[eigen][eigen];
        if (eigenvalue <= threshold || eigenvalue <= VSMALL)
        {
            continue;
        }
        const scalar inverseEigenvalue = 1/eigenvalue;
        for (label row = 0; row < 3; ++row)
        {
            for (label column = 0; column < 3; ++column)
            {
                inverse[row][column] +=
                    eigenvectors[row][eigen]
                   *eigenvectors[column][eigen]
                   *inverseEigenvalue;
            }
        }
    }
}

}


OracleGraphLinearTransfer::OracleGraphLinearTransfer
(
    const fvMesh& mesh,
    const labelList& mapping,
    const label coarseCount
)
:
    mapping_(mapping),
    volumes_(mesh.V()),
    centres_(mesh.C().primitiveField()),
    coarseCount_(coarseCount),
    coarseVolumes_(coarseCount, 0),
    coarseCentres_(coarseCount, vector::zero),
    edgeSource_(),
    edgeTarget_(),
    edgeOffsets_(),
    inverseNormal_(9*coarseCount, 0),
    interfaceOwner_(),
    interfaceNeighbour_(),
    interfaceOffsets_(),
    interfaceWeights_(),
    interfaceStep_(coarseCount, 0)
{
    if
    (
        mapping_.size() != mesh.nCells()
     || volumes_.size() != mapping_.size()
     || centres_.size() != mapping_.size()
     || coarseCount_ <= 0
    )
    {
        FatalErrorInFunction
            << "Invalid arrays for graph linear prolongation."
            << exit(FatalError);
    }

    forAll(mapping_, celli)
    {
        const label aggregate = mapping_[celli];
        coarseVolumes_[aggregate] += volumes_[celli];
        coarseCentres_[aggregate] += volumes_[celli]*centres_[celli];
    }
    forAll(coarseVolumes_, aggregate)
    {
        if (coarseVolumes_[aggregate] <= VSMALL)
        {
            FatalErrorInFunction
                << "Graph aggregate " << aggregate << " has no volume."
                << exit(FatalError);
        }
        coarseCentres_[aggregate] /= coarseVolumes_[aggregate];
    }

    const labelUList& owner = mesh.faceOwner();
    const labelUList& neighbour = mesh.faceNeighbour();
    std::vector<std::pair<label, label>> aggregateEdges;
    std::vector<std::pair<label, label>> interfaceCells;
    aggregateEdges.reserve(mesh.nInternalFaces());
    interfaceCells.reserve(mesh.nInternalFaces());
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const label ownerCell = owner[facei];
        const label neighbourCell = neighbour[facei];
        const label ownerAggregate = mapping_[ownerCell];
        const label neighbourAggregate = mapping_[neighbourCell];
        if (ownerAggregate != neighbourAggregate)
        {
            interfaceCells.emplace_back(ownerCell, neighbourCell);
            label a = ownerAggregate;
            label b = neighbourAggregate;
            if (a > b)
            {
                std::swap(a, b);
            }
            aggregateEdges.emplace_back(a, b);
        }
    }
    std::sort(aggregateEdges.begin(), aggregateEdges.end());
    aggregateEdges.erase
    (
        std::unique(aggregateEdges.begin(), aggregateEdges.end()),
        aggregateEdges.end()
    );

    const label nDirectedEdges = 2*aggregateEdges.size();
    edgeSource_.setSize(nDirectedEdges);
    edgeTarget_.setSize(nDirectedEdges);
    edgeOffsets_.setSize(nDirectedEdges);
    forAll(aggregateEdges, edgei)
    {
        const label a = aggregateEdges[edgei].first;
        const label b = aggregateEdges[edgei].second;
        const label forward = 2*edgei;
        edgeSource_[forward] = a;
        edgeTarget_[forward] = b;
        edgeOffsets_[forward] = coarseCentres_[b] - coarseCentres_[a];
        edgeSource_[forward + 1] = b;
        edgeTarget_[forward + 1] = a;
        edgeOffsets_[forward + 1] = -edgeOffsets_[forward];
    }

    scalarField normal(9*coarseCount_, 0);
    forAll(edgeSource_, edgei)
    {
        const vector& offset = edgeOffsets_[edgei];
        const scalar weight = 1/max(magSqr(offset), scalar(1e-24));
        for (label row = 0; row < 3; ++row)
        {
            for (label column = 0; column < 3; ++column)
            {
                normal[9*edgeSource_[edgei] + 3*row + column] +=
                    offset[row]*offset[column]*weight;
            }
        }
    }
    for (label aggregate = 0; aggregate < coarseCount_; ++aggregate)
    {
        scalar matrix[3][3];
        scalar inverse[3][3];
        for (label row = 0; row < 3; ++row)
        {
            for (label column = 0; column < 3; ++column)
            {
                matrix[row][column] =
                    normal[9*aggregate + 3*row + column];
            }
        }
        symmetricPseudoInverse3(matrix, inverse);
        for (label row = 0; row < 3; ++row)
        {
            for (label column = 0; column < 3; ++column)
            {
                inverseNormal_[9*aggregate + 3*row + column] =
                    inverse[row][column];
            }
        }
    }

    interfaceOwner_.setSize(interfaceCells.size());
    interfaceNeighbour_.setSize(interfaceCells.size());
    interfaceOffsets_.setSize(interfaceCells.size());
    interfaceWeights_.setSize(interfaceCells.size());
    scalarField cellDegree(mapping_.size(), 0);
    forAll(interfaceOwner_, facei)
    {
        const label own = interfaceCells[facei].first;
        const label nei = interfaceCells[facei].second;
        const vector offset = centres_[nei] - centres_[own];
        const scalar weight = 1/max(magSqr(offset), scalar(1e-24));
        interfaceOwner_[facei] = own;
        interfaceNeighbour_[facei] = nei;
        interfaceOffsets_[facei] = offset;
        interfaceWeights_[facei] = weight;
        cellDegree[own] += weight;
        cellDegree[nei] += weight;
    }

    scalarField maxDensity(coarseCount_, 0);
    forAll(mapping_, celli)
    {
        const label aggregate = mapping_[celli];
        maxDensity[aggregate] = max
        (
            maxDensity[aggregate],
            cellDegree[celli]/volumes_[celli]
        );
    }
    forAll(interfaceStep_, aggregate)
    {
        if (maxDensity[aggregate] > 0)
        {
            interfaceStep_[aggregate] = 0.15/maxDensity[aggregate];
        }
    }
}


void OracleGraphLinearTransfer::linearRoundTrip
(
    Field<scalar>& values
) const
{
    if (values.size() != mapping_.size())
    {
        FatalErrorInFunction
            << "Graph pressure field size does not match native mesh."
            << exit(FatalError);
    }

    scalarField coarseValues(coarseCount_, 0);
    forAll(values, celli)
    {
        coarseValues[mapping_[celli]] += volumes_[celli]*values[celli];
    }
    forAll(coarseValues, aggregate)
    {
        coarseValues[aggregate] /= coarseVolumes_[aggregate];
    }

    vectorField rhs(coarseCount_, vector::zero);
    forAll(edgeSource_, edgei)
    {
        const label source = edgeSource_[edgei];
        const label target = edgeTarget_[edgei];
        const vector& offset = edgeOffsets_[edgei];
        rhs[source] += offset
            *((coarseValues[target] - coarseValues[source])
             / max(magSqr(offset), scalar(1e-24)));
    }

    vectorField gradient(coarseCount_, vector::zero);
    forAll(gradient, aggregate)
    {
        for (label row = 0; row < 3; ++row)
        {
            for (label column = 0; column < 3; ++column)
            {
                gradient[aggregate][row] +=
                    inverseNormal_[9*aggregate + 3*row + column]
                   *rhs[aggregate][column];
            }
        }
    }

    scalarField correction(mapping_.size(), 0);
    scalarField correctionMean(coarseCount_, 0);
    forAll(correction, celli)
    {
        const label aggregate = mapping_[celli];
        correction[celli] =
            (centres_[celli] - coarseCentres_[aggregate]) & gradient[aggregate];
        correctionMean[aggregate] += volumes_[celli]*correction[celli];
    }
    forAll(correctionMean, aggregate)
    {
        correctionMean[aggregate] /= coarseVolumes_[aggregate];
    }
    forAll(correction, celli)
    {
        correction[celli] -= correctionMean[mapping_[celli]];
    }

    scalarField lower(coarseValues);
    scalarField upper(coarseValues);
    forAll(edgeSource_, edgei)
    {
        const label source = edgeSource_[edgei];
        const scalar neighbourValue = coarseValues[edgeTarget_[edgei]];
        lower[source] = min(lower[source], neighbourValue);
        upper[source] = max(upper[source], neighbourValue);
    }

    scalarField minimumCorrection(coarseCount_, 0);
    scalarField maximumCorrection(coarseCount_, 0);
    forAll(correction, celli)
    {
        const label aggregate = mapping_[celli];
        minimumCorrection[aggregate] = min
        (
            minimumCorrection[aggregate],
            correction[celli]
        );
        maximumCorrection[aggregate] = max
        (
            maximumCorrection[aggregate],
            correction[celli]
        );
    }

    scalarField limiter(coarseCount_, 1);
    forAll(limiter, aggregate)
    {
        if (maximumCorrection[aggregate] > 0)
        {
            limiter[aggregate] = min
            (
                limiter[aggregate],
                (upper[aggregate] - coarseValues[aggregate])
               /maximumCorrection[aggregate]
            );
        }
        if (minimumCorrection[aggregate] < 0)
        {
            limiter[aggregate] = min
            (
                limiter[aggregate],
                (lower[aggregate] - coarseValues[aggregate])
               /minimumCorrection[aggregate]
            );
        }
        limiter[aggregate] = min(max(limiter[aggregate], scalar(0)), scalar(1));
    }

    scalarField result(mapping_.size(), 0);
    forAll(result, celli)
    {
        const label aggregate = mapping_[celli];
        result[celli] = coarseValues[aggregate]
            + limiter[aggregate]*correction[celli];
    }

    for (label iteration = 0; iteration < 3; ++iteration)
    {
        scalarField derivative(mapping_.size(), 0);
        forAll(interfaceOwner_, facei)
        {
            const label own = interfaceOwner_[facei];
            const label nei = interfaceNeighbour_[facei];
            const label ownAggregate = mapping_[own];
            const label neiAggregate = mapping_[nei];
            const scalar desired = interfaceOffsets_[facei] &
            (
                0.5
               *(gradient[ownAggregate] + gradient[neiAggregate])
            );
            const scalar residual = interfaceWeights_[facei]
                *(result[nei] - result[own] - desired);
            derivative[nei] += residual;
            derivative[own] -= residual;
        }

        scalarField update(mapping_.size(), 0);
        scalarField updateMean(coarseCount_, 0);
        forAll(update, celli)
        {
            update[celli] = -derivative[celli]/volumes_[celli];
            updateMean[mapping_[celli]] += volumes_[celli]*update[celli];
        }
        forAll(updateMean, aggregate)
        {
            updateMean[aggregate] /= coarseVolumes_[aggregate];
        }
        forAll(update, celli)
        {
            const label aggregate = mapping_[celli];
            update[celli] =
                (update[celli] - updateMean[aggregate])
               *interfaceStep_[aggregate];
        }

        scalarField interfaceLimiter(coarseCount_, 1);
        forAll(update, celli)
        {
            const label aggregate = mapping_[celli];
            scalar fraction = 1;
            if (update[celli] > 0)
            {
                fraction = (upper[aggregate] - result[celli])/update[celli];
            }
            else if (update[celli] < 0)
            {
                fraction = (lower[aggregate] - result[celli])/update[celli];
            }
            interfaceLimiter[aggregate] = min
            (
                interfaceLimiter[aggregate],
                min(max(fraction, scalar(0)), scalar(1))
            );
        }
        forAll(result, celli)
        {
            const label aggregate = mapping_[celli];
            result[celli] += update[celli]*interfaceLimiter[aggregate];
        }
    }

    values = result;
}

labelList oracleGraphCoarsen
(
    const polyMesh& mesh,
    const label levels,
    label& coarseCount
)
{
    if (levels < 0)
    {
        FatalErrorInFunction
            << "graphCoarseningLevels must be non-negative."
            << exit(FatalError);
    }

    const label nCells = mesh.nCells();
    std::vector<std::pair<label, label>> edges;
    edges.reserve(mesh.nInternalFaces());
    const labelUList& owner = mesh.faceOwner();
    const labelUList& neighbour = mesh.faceNeighbour();
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        label a = owner[facei];
        label b = neighbour[facei];
        if (a != b)
        {
            if (a > b)
            {
                std::swap(a, b);
            }
            edges.emplace_back(a, b);
        }
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    labelList mapping(nCells);
    forAll(mapping, celli)
    {
        mapping[celli] = celli;
    }
    coarseCount = nCells;

    for (label level = 0; level < levels; ++level)
    {
        std::vector<std::pair<label, label>> directed;
        directed.reserve(2*edges.size());
        for (const auto& edge : edges)
        {
            directed.emplace_back(edge.first, edge.second);
            directed.emplace_back(edge.second, edge.first);
        }
        std::sort(directed.begin(), directed.end());
        directed.erase(std::unique(directed.begin(), directed.end()), directed.end());

        labelList offsets(coarseCount + 1, 0);
        labelList degree(coarseCount, 0);
        for (const auto& edge : directed)
        {
            offsets[edge.first + 1]++;
            degree[edge.first]++;
        }
        for (label node = 0; node < coarseCount; ++node)
        {
            offsets[node + 1] += offsets[node];
        }

        std::vector<label> order(coarseCount);
        std::iota(order.begin(), order.end(), label(0));
        std::stable_sort
        (
            order.begin(),
            order.end(),
            [&degree](const label a, const label b)
            {
                return degree[a] < degree[b];
            }
        );

        labelList localMapping(coarseCount, -1);
        label nextCoarseCount = 0;
        for (const label node : order)
        {
            if (localMapping[node] >= 0)
            {
                continue;
            }

            localMapping[node] = nextCoarseCount;
            for (label edgei = offsets[node]; edgei < offsets[node + 1]; ++edgei)
            {
                const label adjacent = directed[edgei].second;
                if (localMapping[adjacent] < 0)
                {
                    localMapping[adjacent] = nextCoarseCount;
                    break;
                }
            }
            nextCoarseCount++;
        }

        forAll(mapping, celli)
        {
            mapping[celli] = localMapping[mapping[celli]];
        }

        std::vector<std::pair<label, label>> coarseEdges;
        coarseEdges.reserve(edges.size());
        for (const auto& edge : edges)
        {
            label a = localMapping[edge.first];
            label b = localMapping[edge.second];
            if (a != b)
            {
                if (a > b)
                {
                    std::swap(a, b);
                }
                coarseEdges.emplace_back(a, b);
            }
        }
        std::sort(coarseEdges.begin(), coarseEdges.end());
        coarseEdges.erase
        (
            std::unique(coarseEdges.begin(), coarseEdges.end()),
            coarseEdges.end()
        );
        edges.swap(coarseEdges);
        coarseCount = nextCoarseCount;
    }

    return mapping;
}

}