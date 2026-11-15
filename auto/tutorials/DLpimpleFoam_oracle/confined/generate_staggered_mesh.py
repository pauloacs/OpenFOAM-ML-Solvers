from pathlib import Path

import sys

# ---------------------------------------------------------------------
# GLOBAL MESH DENSITY SCALING
#
# 1.00 -> original mesh (~3M cells)
# 0.50 -> ~375k cells
# 0.40 -> ~200k cells
# 0.30 -> ~80k cells
# 0.25 -> ~50k cells
# ---------------------------------------------------------------------

MESH_SCALE = float(sys.argv[1]) if len(sys.argv) > 1 else 0.40


def scaled(n):
    return max(1, int(round(n * MESH_SCALE)))

CASE = Path(__file__).resolve().parent

# ---------------------------------------------------------------------
# DOMAIN
# ---------------------------------------------------------------------

XMIN = -1.0
XMAX = 6.0

YMIN = -0.15
YMAX = 0.75

Z_TOP = 0.5
Z_CELLS = scaled(40)

# ---------------------------------------------------------------------
# 3 x 2 obstacle arrangement
# all obstacles have different dimensions
# ---------------------------------------------------------------------

RECT_OBSTACLES = [

    # bottom row
    (-0.05, 0.20, 0.05, 0.18),
    ( 0.85, 1.15, 0.03, 0.30),
    ( 2.00, 2.22, 0.08, 0.18),

    # top row
    ( 0.40, 0.55, 0.35, 0.55),
    ( 1.55, 1.75, 0.32, 0.58),
    # reduced width, less bulky
    ( 2.90, 3.10, 0.28, 0.65),]

# refinement distance around obstacles
REFINE = 0.05

# ---------------------------------------------------------------------
# BUILD MESH LINES
# ---------------------------------------------------------------------

x_lines = [XMIN, XMAX]
y_lines = [YMIN, YMAX]

for xmin, xmax, ymin, ymax in RECT_OBSTACLES:

    x_lines.extend([
        max(XMIN, xmin - REFINE),
        xmin,
        xmax,
        min(XMAX, xmax + REFINE),
    ])

    y_lines.extend([
        max(YMIN, ymin - REFINE),
        ymin,
        ymax,
        min(YMAX, ymax + REFINE),
    ])

X_EDGES = sorted(set(round(x, 8) for x in x_lines))
Y_EDGES = sorted(set(round(y, 8) for y in y_lines))

# ---------------------------------------------------------------------
# CELL COUNTS
# ---------------------------------------------------------------------

X_CELLS = []

for x0, x1 in zip(X_EDGES[:-1], X_EDGES[1:]):

    dx = x1 - x0

    if dx <= 0.05:
        X_CELLS.append(scaled(10))
    elif dx <= 0.10:
        X_CELLS.append(scaled(16))
    elif dx <= 0.20:
        X_CELLS.append(scaled(24))
    else:
        X_CELLS.append(max(scaled(20), int(dx * 50)))

Y_CELLS = []

for y0, y1 in zip(Y_EDGES[:-1], Y_EDGES[1:]):

    dy = y1 - y0

    if dy <= 0.05:
        Y_CELLS.append(scaled(8))
    elif dy <= 0.10:
        Y_CELLS.append(scaled(14))
    elif dy <= 0.20:
        Y_CELLS.append(scaled(20))
    else:
        Y_CELLS.append(max(scaled(10), int(dy * 45)))

NX = len(X_CELLS)
NY = len(Y_CELLS)

NXY = (NX + 1) * (NY + 1)

# ---------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------

def vertex(i, j, k):
    return k * NXY + j * (NX + 1) + i


def face(values):
    return "(" + " ".join(str(v) for v in values) + ")"


def inside_obstacle(xc, yc):

    for xmin, xmax, ymin, ymax in RECT_OBSTACLES:

        if xmin <= xc <= xmax and ymin <= yc <= ymax:
            return True

    return False


def cell_blocked(i, j):

    if i < 0 or i >= NX:
        return False

    if j < 0 or j >= NY:
        return False

    xc = 0.5 * (X_EDGES[i] + X_EDGES[i + 1])
    yc = 0.5 * (Y_EDGES[j] + Y_EDGES[j + 1])

    return inside_obstacle(xc, yc)


# ---------------------------------------------------------------------
# BLOCKMESHDICT
# ---------------------------------------------------------------------

def build_dictionary():

    lines = [
        "FoamFile",
        "{",
        "    version 2.0;",
        "    format ascii;",
        "    class dictionary;",
        "    object blockMeshDict;",
        "}",
        "",
        "convertToMeters 1;",
        "",
        "vertices",
        "(",
    ]

    for k, z in enumerate((0.0, Z_TOP)):
        for j, y in enumerate(Y_EDGES):
            for i, x in enumerate(X_EDGES):
                lines.append(
                    f"    ({x:.8f} {y:.8f} {z:.8f}) // {vertex(i,j,k)}"
                )

    lines.extend((");", "", "blocks", "("))

    patches = {
        name: []
        for name in (
            "inlet",
            "outlet",
            "front",
            "back",
            "top",
            "bot",
            "obstacle",
        )
    }

    for j in range(NY):
        for i in range(NX):

            xc = 0.5 * (X_EDGES[i] + X_EDGES[i + 1])
            yc = 0.5 * (Y_EDGES[j] + Y_EDGES[j + 1])

            if inside_obstacle(xc, yc):
                continue

            v000 = vertex(i, j, 0)
            v100 = vertex(i + 1, j, 0)
            v110 = vertex(i + 1, j + 1, 0)
            v010 = vertex(i, j + 1, 0)

            v001 = vertex(i, j, 1)
            v101 = vertex(i + 1, j, 1)
            v111 = vertex(i + 1, j + 1, 1)
            v011 = vertex(i, j + 1, 1)

            lines.append(
                f"    hex ({v000} {v100} {v110} {v010} "
                f"{v001} {v101} {v111} {v011}) "
                f"inner "
                f"({X_CELLS[i]} {Y_CELLS[j]} {Z_CELLS}) "
                f"simpleGrading (1 1 1)"
            )

            patches["bot"].append(face((v000, v010, v110, v100)))
            patches["top"].append(face((v001, v101, v111, v011)))

            if i == 0:
                patches["inlet"].append(
                    face((v000, v001, v011, v010))
                )
            elif cell_blocked(i - 1, j):
                patches["obstacle"].append(
                    face((v000, v001, v011, v010))
                )

            if i == NX - 1:
                patches["outlet"].append(
                    face((v100, v110, v111, v101))
                )
            elif cell_blocked(i + 1, j):
                patches["obstacle"].append(
                    face((v100, v110, v111, v101))
                )

            if j == 0:
                patches["front"].append(
                    face((v000, v100, v101, v001))
                )
            elif cell_blocked(i, j - 1):
                patches["obstacle"].append(
                    face((v000, v100, v101, v001))
                )

            if j == NY - 1:
                patches["back"].append(
                    face((v010, v011, v111, v110))
                )
            elif cell_blocked(i, j + 1):
                patches["obstacle"].append(
                    face((v010, v011, v111, v110))
                )

    lines.extend((");", "", "edges", "(", ");", "", "boundary", "("))

    patch_types = {
        "inlet": "patch",
        "outlet": "patch",
        "front": "wall",
        "back": "wall",
        "top": "wall",
        "bot": "wall",
        "obstacle": "wall",
    }

    for name, faces in patches.items():

        lines.extend(
            (
                f"    {name}",
                "    {",
                f"        type {patch_types[name]};",
                "        faces",
                "        (",
            )
        )

        lines.extend(
            f"            {entry}"
            for entry in faces
        )

        lines.extend(
            (
                "        );",
                "    }",
            )
        )

    lines.extend(
        (
            ");",
            "",
            "mergePatchPairs",
            "(",
            ");",
            "",
        )
    )

    return "\n".join(lines)


if __name__ == "__main__":

    dictionary = build_dictionary()

    outputs = (
        CASE / "system" / "blockMeshDict",
        CASE / "system" / "blockMeshDict.staggered",
    )

    for output in outputs:
        output.write_text(dictionary, encoding="ascii")

    print(f"Wrote {outputs[0]}")