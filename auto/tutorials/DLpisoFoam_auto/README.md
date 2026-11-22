# DLpisoFoam_auto Tutorial

`graph_shift/pressure_only` is a compact AUTO PISO example based on the
512-cell graph pressure-only tutorial. It trains from samples collected in
the run, then injects pressure and resumes the inner pressure-correction loop.

The case selects `DLpisoFoam_auto`, uses one PIMPLE outer iteration with ten
pressure correctors, and selects PISO auto mode by executable name. Its
mesh-graph configuration uses level 2 and protects two face-adjacent layers
around the `obstacle` patch. On this small mesh, 15 layers would protect all
cells; larger meshes can use 15 layers when appropriate. The case requires
`smUseInputPredictorState true`, no post-injection reconstruction, and no
PIMPLE oracle-compatible override. The solver performs
one predictor momentum solve, injects the learned estimate, and does not solve
the momentum equation again that timestep. Pressure correctors still update
velocity algebraically.

After sourcing OpenFOAM 8, prepare the project Python environment, install
the graph backend, and build the solver from the repository root:

```bash
source prep_env311.sh
python3 -m pip install -e '.[mesh_graph]' mpi4py
wmake auto/source/DLpisoFoam_auto
cd auto/tutorials/DLpisoFoam_auto/graph_shift/pressure_only
bash Allrun
```

Run the graph tutorial in serial. `Allrun` builds and checks the mesh, then
runs the solver; it refuses to overwrite an existing mesh, `ML_data`, or
`log.solver`. Review `log.solver` for graph size, sampling, training, and
activation. `ML_data` contains the sampled HDF5 data, graph cache, model
weights, normalization, and validation summaries.

`bash Allclean` removes generated mesh, training data, logs, and outputs while
retaining the time-zero fields and source dictionaries. This deletes trained
artifacts, so retain them elsewhere before cleaning if needed. The compact
case demonstrates the workflow; it does not establish accuracy or speedup.