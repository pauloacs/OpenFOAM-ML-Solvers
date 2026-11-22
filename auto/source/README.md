# AUTO Solver Sources

| Source directory | Executable | Purpose |
| --- | --- | --- |
| `DLpimpleFoam_auto/` | `DLpimpleFoam_auto` | Isothermal PIMPLE solver with online surrogate sampling and training. |
| `DLpisoFoam_auto/` | `DLpisoFoam_auto` | Isothermal AUTO predictor followed by PISO pressure corrections. |
| `DLbuoyantPimpleFoam_auto/` | `DLbuoyantPimpleFoam_auto` | Buoyant, heat-transfer AUTO solver. |
| `DLpimpleFoam_oracle/` | `DLpimpleFoam_oracle` | Reference-injection diagnostic for PIMPLE. |
| `DLpisoFoam_oracle/` | `DLpisoFoam_oracle` | Reference-injection diagnostic for single-outer PISO. |

Each directory has its own `Make/files` and can be built with `wmake` after
sourcing the OpenFOAM and Python solver environment. The PISO oracle source
entry includes the shared implementation and build options from
`DLpimpleFoam_oracle/`; changes to oracle algorithms therefore apply to both
executables. Keep the two directories together when building either oracle.

`DLpisoFoam_auto/` reuses the AUTO model, Python bridge, sampling, and training
sources from `DLpimpleFoam_auto/`. It builds a separate executable whose mode
is selected by its name. Configure `nOuterCorrectors 1`, keep
`smUseInputPredictorState true`, and disable post-injection reconstruction and
PIMPLE-only pressure-force stages for this path. It performs one predictor
momentum solve, injects the AUTO pressure/optional velocity estimate, then
enters the inner PISO pressure-corrector loop without solving UEqn again.

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
selects limited linear pressure prolongation. `graphLinearVelocityProlongation`
selects component-wise linear velocity prolongation separately; both options
default to `false` and require `useGraphCoarsening true`. These settings are
read at startup; restart the executable after changing them.

When present, `graphFullCoarsePrediction` takes precedence over
`graphPressureMode`: `false` selects `ddpResidual`, while `true` selects
`dp`. If the flag is omitted, set `graphPressureMode` to either `dp` or
`ddpResidual`; `dpResidual` is not a supported value. Check the
`[ORACLE GRAPH]` startup line for the resolved settings.

For example, the following protects the obstacle patch's owner cells and
14 additional face-adjacent layers throughout level-2 coarsening:

```text
useGraphCoarsening true;
graphCoarseningLevels 2;
graphObstaclePatch obstacle;
graphObstacleProtectionLayers 15;
```

The default protection count is zero. Protected cells remain individual graph
nodes at every level; the startup line reports both the protected-cell count
and the resulting graph size. A missing or empty patch is an error when
protection is enabled. The AUTO graph backend uses separate Python settings;
see [Obstacle-Layer Protection](../README.md#obstacle-layer-protection).

Oracle graph velocity supports four target/composition modes through
`graphVelocityMode`:

- `dU`: graph-transfer the total reference increment `U_final - U_start`.
- `ddU`: graph-transfer the residual `dU - dUPrev`, then add the weighted
	native `oracleDUPrevTruth`; set `graphVelocityPreviousIncrementWeightFactor`.
- `dUCorr`: graph-transfer `dU - dUStar`, then add the native current predictor
	increment `dUStar`.
- `ddUCorr`: graph-transfer `dU - dUStar - dUCorrPrev`, then add native `dUStar`
	and weighted `dUCorrPrev`; set `graphVelocityPreviousCorrectionWeightFactor`.

`outputWeightFactor` scales each graph-predicted velocity field. For `ddU` and
`ddUCorr`, the history baseline is added after graph transfer, at native mesh
resolution. `ddUCorrResidual` remains a compatibility alias for `ddUCorr`.
`dUCorrResidual` and `duCorrResidual` are aliases for `dUCorr`; the shared
parser supports these names in both PIMPLE and PISO oracle executables.
All three residual/composed modes (`ddU`, `dUCorr`, `ddUCorr`) require graph
coarsening, velocity injection, and `pressureForceThenReconstruct`. A missing
`oracleDUPrevTruth` restart field initializes that baseline to zero; use a
restart written by the updated oracle for history-continuous `ddU` tests.

The confined and unconfined mesh-generated oracle examples are under
`auto/tutorials/DLpimpleFoam_oracle/`. Generated meshes, field snapshots,
datasets and logs are not source files.