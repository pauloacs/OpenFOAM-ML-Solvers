# AUTO Solver Sources

| Source directory | Executable | Purpose |
| --- | --- | --- |
| `DLpimpleFoam_auto/` | `DLpimpleFoam_auto` | Isothermal PIMPLE solver with online surrogate sampling and training. |
| `DLbuoyantPimpleFoam_auto/` | `DLbuoyantPimpleFoam_auto` | Buoyant, heat-transfer AUTO solver. |
| `DLpimpleFoam_oracle/` | `DLpimpleFoam_oracle` | Reference-injection diagnostic for PIMPLE. |
| `DLpisoFoam_oracle/` | `DLpisoFoam_oracle` | Reference-injection diagnostic for single-outer PISO. |

Each directory has its own `Make/files` and can be built with `wmake` after
sourcing the OpenFOAM and Python solver environment. The PISO oracle source
entry includes the shared implementation and build options from
`DLpimpleFoam_oracle/`; changes to oracle algorithms therefore apply to both
executables. Keep the two directories together when building either oracle.

## Oracle Modes

The PIMPLE oracle uses its first outer loop for the momentum predictor, then
injects the reference estimate and continues the configured PIMPLE loops. The
PISO oracle defaults to `oraclePisoMode true`: it solves the momentum predictor
once, injects pressure and velocity, then reuses that momentum matrix for the
inner pressure correctors. It requires `nOuterCorrectors 1`,
`integrationMode pressureForceThenReconstruct`, and
`stateReconstructionCorrectors 0`. Reconstruction warm-start is skipped in
this direct PISO path. Pressure correctors still update velocity algebraically.

The executable name selects the default mode; `oraclePisoMode` in
`system/oracleDict` explicitly overrides it. Both oracles use a reference
pass for idealized diagnostic comparisons, not a trained surrogate.

Graph transfer is controlled in `oracleDict`: `graphFullCoarsePrediction
false` adds native `dpPrev` to the graph-transferred `ddp` residual; `true`
transfers the full `dp` increment. `graphLinearProlongation` independently
selects limited linear pressure prolongation. These settings are read at
startup; restart the executable after changing them.

The confined and unconfined mesh-generated oracle examples are under
`auto/tutorials/DLpimpleFoam_oracle/`. Generated meshes, field snapshots,
datasets and logs are not source files.