# Unconfined Oracle Tutorial

Six staggered finite-height cylinders in a far-field domain. `Allrun` generates
the original approximately 411,287-cell mesh with `generate_cylinder_mesh.py`,
`blockMesh`, `snappyHexMesh`, and `renumberMesh`. NumPy is required.

The historical mesh fails two `checkMesh` geometry checks: its source run
reported 78 highly skewed faces and 3,408 concave cells. `Allrun` stops on mesh
quality failures. Inspect `log.checkMesh` and assess suitability before using
this mesh; these defects can affect the results.

To reproduce the historical case despite these warnings, after sourcing
OpenFOAM and the project Python solver environment and building
`DLpimpleFoam_oracle`, use:

```bash
ALLOW_MESH_QUALITY_FAILURE=1 bash Allrun
```

If an earlier attempt already generated files, preserve any results you need
and run `bash Allclean` first. The override does not repair the mesh. See
[the oracle tutorial guide](../README.md) for the other cases and solver
settings.