from pathlib import Path
import math
import numpy as np

MULT_FACTOR = 0.8
BACKGROUND_X_CELLS = 182

CASE = Path(__file__).resolve().parent
X_CELLS = np.array([BACKGROUND_X_CELLS * MULT_FACTOR], dtype=int)
Y_CORE_EDGES = (-0.60, -0.15, 0.0, 0.09, 0.25, 0.35, 0.51, 0.60, 0.75, 1.25)
Y_CORE_CELLS = (np.array([5, 8, 5, 9, 6, 9, 5, 8, 5]) * MULT_FACTOR).astype(int)
UPSTREAM_CLEARANCE_DIAMETERS = 8
DOWNSTREAM_CLEARANCE_DIAMETERS = 25
SIDE_CLEARANCE_DIAMETERS = 10.0
FAR_SIDE_CELLS = 4
Z_MIN = 0.0
BACKGROUND_Z_CELLS = (np.array([18, 6]) * MULT_FACTOR).astype(int)
BACKGROUND_Z_GRADING = 4.0
CYLINDER_WALL_OVERLAP = 0.02
CYLINDER_SEGMENTS = 48
CYLINDER_HEIGHT_DIAMETERS = 5.0
TOP_CLEARANCE_DIAMETERS = 20.0
WAKE_REFINEMENT_UPSTREAM_DIAMETERS = 2.0
WAKE_REFINEMENT_DIAMETERS = 5.0

# (center x, center y, radius); alternate rows are shifted streamwise.
CYLINDERS = (
    (0.32, 0.17, 0.07),
    (0.64, 0.53, 0.07),
    (0.96, 0.17, 0.07),
    (1.28, 0.53, 0.07),
    (1.60, 0.17, 0.07),
    (1.92, 0.53, 0.07),
)
MAX_CYLINDER_DIAMETER = max(2 * radius for _, _, radius in CYLINDERS)
OBSTACLE_X_MIN = min(center_x - radius for center_x, _, radius in CYLINDERS)
OBSTACLE_X_MAX = max(center_x + radius for center_x, _, radius in CYLINDERS)
OBSTACLE_Y_MIN = min(center_y - radius for _, center_y, radius in CYLINDERS)
OBSTACLE_Y_MAX = max(center_y + radius for _, center_y, radius in CYLINDERS)
X_MIN = OBSTACLE_X_MIN - UPSTREAM_CLEARANCE_DIAMETERS * MAX_CYLINDER_DIAMETER
X_MAX = OBSTACLE_X_MAX + DOWNSTREAM_CLEARANCE_DIAMETERS * MAX_CYLINDER_DIAMETER
WAKE_REFINEMENT_X_MIN = OBSTACLE_X_MIN - WAKE_REFINEMENT_UPSTREAM_DIAMETERS * MAX_CYLINDER_DIAMETER
WAKE_REFINEMENT_X_MAX = OBSTACLE_X_MAX + WAKE_REFINEMENT_DIAMETERS * MAX_CYLINDER_DIAMETER
X_EDGES = (X_MIN, X_MAX)
Y_MIN = OBSTACLE_Y_MIN - SIDE_CLEARANCE_DIAMETERS * MAX_CYLINDER_DIAMETER
Y_MAX = OBSTACLE_Y_MAX + SIDE_CLEARANCE_DIAMETERS * MAX_CYLINDER_DIAMETER
Y_EDGES = (Y_MIN, *Y_CORE_EDGES, Y_MAX)
Y_CELLS = np.concatenate(([FAR_SIDE_CELLS], Y_CORE_CELLS, [FAR_SIDE_CELLS])).astype(int)
CYLINDER_HEIGHT = CYLINDER_HEIGHT_DIAMETERS * MAX_CYLINDER_DIAMETER
Z_MAX = CYLINDER_HEIGHT + TOP_CLEARANCE_DIAMETERS * MAX_CYLINDER_DIAMETER
Z_LEVELS = (Z_MIN, CYLINDER_HEIGHT, Z_MAX)
LOCATION_IN_MESH = ((X_MIN + OBSTACLE_X_MIN) / 2, (Y_MIN + Y_MAX) / 2, (Z_MIN + Z_MAX) / 2)


def expanded_edges(edges, cell_counts):
    if (len(edges) != len(cell_counts) + 1 or any(count < 1 for count in cell_counts)
            or any(right <= left for left, right in zip(edges, edges[1:]))):
        raise ValueError('Each positive cell count must match one interval in the edge list.')
    result = []
    for left, right, count in zip(edges[:-1], edges[1:], cell_counts):
        result.extend(left + (right - left) * index / count for index in range(count))
    result.append(edges[-1])
    return tuple(result)


def vertex(i, j, k, nx, ny):
    return k * (ny + 1) * (nx + 1) + j * (nx + 1) + i


def background_dictionary():
    x = expanded_edges(X_EDGES, X_CELLS)
    y = expanded_edges(Y_EDGES, Y_CELLS)
    nx, ny = len(x) - 1, len(y) - 1
    lines = [
        'FoamFile', '{', '    version 2.0;', '    format ascii;',
        '    class dictionary;', '    object blockMeshDict;', '}', '',
        'convertToMeters 1;', '', 'vertices', '(',
    ]
    for k, z in enumerate(Z_LEVELS):
        for j, y_value in enumerate(y):
            for i, x_value in enumerate(x):
                lines.append(f'    ({x_value:.8f} {y_value:.8f} {z:.8f}) // {vertex(i,j,k,nx,ny)}')
    lines.extend((');', '', 'blocks', '('))

    patches = {name: [] for name in ('inlet', 'outlet', 'front', 'back', 'top', 'bot')}
    for k, z_cells in enumerate(BACKGROUND_Z_CELLS):
        for j in range(ny):
            for i in range(nx):
                v000, v100 = vertex(i, j, k, nx, ny), vertex(i + 1, j, k, nx, ny)
                v110, v010 = vertex(i + 1, j + 1, k, nx, ny), vertex(i, j + 1, k, nx, ny)
                v001, v101 = vertex(i, j, k + 1, nx, ny), vertex(i + 1, j, k + 1, nx, ny)
                v111, v011 = vertex(i + 1, j + 1, k + 1, nx, ny), vertex(i, j + 1, k + 1, nx, ny)
                lines.append(f'    hex ({v000} {v100} {v110} {v010} {v001} {v101} {v111} {v011}) '
                             f'inner (1 1 {z_cells}) simpleGrading (1 1 {BACKGROUND_Z_GRADING:g})')
                if i == 0:
                    patches['inlet'].append((v000, v001, v011, v010))
                if i == nx - 1:
                    patches['outlet'].append((v100, v110, v111, v101))
                if j == 0:
                    patches['front'].append((v000, v100, v101, v001))
                if j == ny - 1:
                    patches['back'].append((v010, v011, v111, v110))
                if k == 0:
                    patches['bot'].append((v000, v010, v110, v100))
                if k == len(BACKGROUND_Z_CELLS) - 1:
                    patches['top'].append((v001, v101, v111, v011))

    lines.extend((');', '', 'edges', '(', ');', '', 'boundary', '('))
    patch_types = {'inlet': 'patch', 'outlet': 'patch', 'front': 'symmetry', 'back': 'symmetry',
                   'top': 'symmetry', 'bot': 'wall'}
    for name, faces in patches.items():
        lines.extend((f'    {name}', '    {', f'        type {patch_types[name]};', '        faces', '        ('))
        lines.extend('            (' + ' '.join(map(str, face)) + ')' for face in faces)
        lines.extend(('        );', '    }'))
    lines.extend((');', '', 'mergePatchPairs', '(', ');', ''))
    return '\n'.join(lines)


def triangle_normal(a, b, c):
    ab = tuple(b[i] - a[i] for i in range(3))
    ac = tuple(c[i] - a[i] for i in range(3))
    normal = (ab[1]*ac[2] - ab[2]*ac[1], ab[2]*ac[0] - ab[0]*ac[2], ab[0]*ac[1] - ab[1]*ac[0])
    length = math.sqrt(sum(value * value for value in normal))
    return tuple(value / length for value in normal)


def cylinder_stl():
    lines = ['solid cylinders']
    segments = CYLINDER_SEGMENTS
    z0 = Z_MIN - CYLINDER_WALL_OVERLAP
    for cx, cy, radius in CYLINDERS:
        z1 = CYLINDER_HEIGHT_DIAMETERS * 2 * radius + CYLINDER_WALL_OVERLAP
        rings = []
        for z in (z0, z1):
            rings.append([(cx + radius * math.cos(2 * math.pi * i / segments),
                           cy + radius * math.sin(2 * math.pi * i / segments), z)
                          for i in range(segments)])
        bottom, top = rings
        bottom_center = (cx, cy, z0)
        top_center = (cx, cy, z1)
        triangles = []
        for i in range(segments):
            j = (i + 1) % segments
            triangles.extend(((bottom_center, bottom[j], bottom[i]),
                              (top_center, top[i], top[j]),
                              (bottom[i], bottom[j], top[j]),
                              (bottom[i], top[j], top[i])))
        for a, b, c in triangles:
            normal = triangle_normal(a, b, c)
            lines.append('  facet normal ' + ' '.join(f'{v:.9g}' for v in normal))
            lines.append('    outer loop')
            for point in (a, b, c):
                lines.append('      vertex ' + ' '.join(f'{v:.9f}' for v in point))
            lines.extend(('    endloop', '  endfacet'))
    lines.append('endsolid cylinders')
    return '\n'.join(lines) + '\n'


def snappy_dictionary():
    if WAKE_REFINEMENT_UPSTREAM_DIAMETERS < 0 or WAKE_REFINEMENT_DIAMETERS < 0:
        raise ValueError('Wake refinement buffers must be non-negative.')
    if WAKE_REFINEMENT_X_MIN < X_MIN or WAKE_REFINEMENT_X_MAX > X_MAX:
        raise ValueError('Wake refinement extends beyond the domain; increase the corresponding '
                         'domain clearance.')

    return '''FoamFile
{
    version 2.0;
    format ascii;
    class dictionary;
    object snappyHexMeshDict;
}

castellatedMesh true;
snap true;
addLayers false;

geometry
{
    cylinders.stl { type triSurfaceMesh; name obstacle; }
    wakeRefinement
    {
        type searchableBox;
        min (__WAKE_REFINEMENT_X_MIN__ __WAKE_REFINEMENT_Y_MIN__ __WAKE_REFINEMENT_Z_MIN__);
        max (__WAKE_REFINEMENT_X_MAX__ __WAKE_REFINEMENT_Y_MAX__ __WAKE_REFINEMENT_Z_MAX__);
    }
}

castellatedMeshControls
{
    maxLocalCells 1000000;
    maxGlobalCells 1500000;
    minRefinementCells 0;
    maxLoadUnbalance 0.10;
    nCellsBetweenLevels 2;
    features ();
    refinementSurfaces
    {
        obstacle { level (1 1); patchInfo { type wall; } }
    }
    resolveFeatureAngle 30;
    refinementRegions
    {
        wakeRefinement
        {
            mode inside;
            levels ((1e15 1));
        }
    }
    locationInMesh (__LOCATION_IN_MESH__);
    allowFreeStandingZoneFaces true;
}

snapControls
{
    nSmoothPatch 5;
    tolerance 1.5;
    nSolveIter 50;
    nRelaxIter 5;
    nFeatureSnapIter 10;
    implicitFeatureSnap true;
    explicitFeatureSnap false;
    multiRegionFeatureSnap false;
}

addLayersControls
{
    relativeSizes true;
    layers {};
    expansionRatio 1.0;
    finalLayerThickness 0.3;
    minThickness 0.1;
    nGrow 0;
    featureAngle 60;
    slipFeatureAngle 30;
    nRelaxIter 3;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedianAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
}

meshQualityControls
{
    maxNonOrtho 65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave 80;
    minVol 1e-13;
    minTetQuality 1e-9;
    minArea -1;
    minTwist 0.02;
    minDeterminant 0.001;
    minFaceWeight 0.02;
    minVolRatio 0.01;
    minTriangleTwist -1;
    nSmoothScale 4;
    errorReduction 0.75;
}

writeFlags (scalarLevels layerSets layerFields);
mergeTolerance 1e-6;
'''.replace('__LOCATION_IN_MESH__', ' '.join(f'{value:.8g}' for value in LOCATION_IN_MESH)) \
    .replace('__WAKE_REFINEMENT_X_MIN__', f'{WAKE_REFINEMENT_X_MIN:.8g}') \
    .replace('__WAKE_REFINEMENT_Y_MIN__', f'{Y_CORE_EDGES[0]:.8g}') \
    .replace('__WAKE_REFINEMENT_Z_MIN__', f'{Z_MIN:.8g}') \
    .replace('__WAKE_REFINEMENT_X_MAX__', f'{WAKE_REFINEMENT_X_MAX:.8g}') \
    .replace('__WAKE_REFINEMENT_Y_MAX__', f'{Y_CORE_EDGES[-1]:.8g}') \
    .replace('__WAKE_REFINEMENT_Z_MAX__', f'{CYLINDER_HEIGHT:.8g}')


if __name__ == '__main__':
    system = CASE / 'system'
    surfaces = CASE / 'constant' / 'triSurface'
    surfaces.mkdir(parents=True, exist_ok=True)
    (system / 'blockMeshDict.background').write_text(background_dictionary(), encoding='ascii')
    (system / 'snappyHexMeshDict').write_text(snappy_dictionary(), encoding='ascii')
    (surfaces / 'cylinders.stl').write_text(cylinder_stl(), encoding='ascii')
    cells = sum(X_CELLS) * sum(Y_CELLS) * sum(BACKGROUND_Z_CELLS)
    print(f'Wrote six-cylinder mesh inputs; background cells={cells}, '
          f'cylinder heights={CYLINDER_HEIGHT_DIAMETERS:g}D, top={Z_MAX:g}m, '
          f'STL facets={len(CYLINDERS)*CYLINDER_SEGMENTS*4}')