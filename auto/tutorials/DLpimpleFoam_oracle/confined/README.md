# Confined Oracle Tutorial

Six staggered rectangular obstacles in a channel. `Allrun` generates the
228,576-cell mesh using `generate_staggered_mesh.py` at scale 0.3, checks the
mesh, and runs a short PIMPLE oracle comparison from time zero.

After sourcing OpenFOAM and the project Python solver environment, build
`DLpimpleFoam_oracle` and run `bash Allrun`. Meshes and results are not
included. `bash Allclean` removes generated outputs and nonzero time
directories; preserve any results you need before cleaning.

See [the oracle tutorial guide](../README.md) for solver settings and the
unconfined case.