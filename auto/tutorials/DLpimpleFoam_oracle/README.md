# Oracle Tutorials

These examples compare a conventional CFD reference pass with an idealized
reference-injection pass. They do not train or run a surrogate model, and
their timings are not a prediction of learned-model speedup.

| Case | Geometry | Mesh generation |
| --- | --- | --- |
| `pressure_only` | Compact square-obstacle channel | `blockMesh` |
| `joint` | Compact square-obstacle channel | `blockMesh` |
| `confined` | Six staggered rectangular obstacles in a channel | `generate_staggered_mesh.py` at scale 0.3, then `blockMesh`; about 228,576 cells |
| `unconfined` | Six staggered finite-height cylinders with a far-field boundary | `generate_cylinder_mesh.py`, `blockMesh`, `snappyHexMesh`, and `renumberMesh`; approximately 411,287 cells |

The two larger cases contain uniform initial fields, solver dictionaries and
mesh-generation scripts; generated meshes, STL, fields and logs are not
shipped. The unconfined generator requires NumPy. Both cases start from time
zero and are configured as short PIMPLE oracle demonstrations, not full
reproductions of the restarted experiments.

After sourcing OpenFOAM and the project Python solver environment, build and
run from the repository root:

```bash
wmake auto/source/DLpimpleFoam_oracle
cd auto/tutorials/DLpimpleFoam_oracle/confined
bash Allrun
```

Use the analogous `unconfined` directory for the cylinder case. `Allrun`
generates the mesh, checks it, and runs the solver, writing `log.*` files. It
refuses to overwrite existing generated files or result times. `bash Allclean`
removes generated mesh inputs, mesh, nonzero times and logs; it preserves the
initial conditions and source dictionaries. Review the case before cleaning if
you need to retain results.

The original unconfined cylinder mesh has known `checkMesh` failures: the
source run reported 78 highly skewed faces and 3,408 concave cells. `Allrun`
stops if mesh checks fail. To reproduce this historical mesh despite those
failures, inspect `log.checkMesh`, preserve any results you need, run
`bash Allclean`, then explicitly set `ALLOW_MESH_QUALITY_FAILURE=1` when
running `Allrun` again. This override does not repair the mesh or establish
that it is suitable for a new simulation.

`system/oracleDict` exposes graph coarsening, pressure composition, linear
prolongation, state reconstruction and warm-start options. The PISO oracle is
a separate executable: build it with `wmake auto/source/DLpisoFoam_oracle`.
It requires a single outer corrector and zero reconstruction correctors;
see [the source guide](../../source/README.md) before adapting a tutorial.