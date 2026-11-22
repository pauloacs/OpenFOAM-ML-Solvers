# AUTO CFD Workflow

This folder groups the current AUTO and oracle tutorial cases. The AUTO
solver samples fields, trains an initial surrogate, and activates it during
the same simulation. The oracle replays reference fields to assess integration
potential; it does not train or load a surrogate.

## Current Components

- `tutorials/DLpimpleFoam_auto/`: gridded pressure-only and graph pressure-only/joint cases.
- `tutorials/DLpisoFoam_auto/`: compact single-outer AUTO PISO case.
- `tutorials/DLpimpleFoam_oracle/`: compact and mesh-generated reference-injection cases.
- `source/DLpimpleFoam_auto/`: self-training OpenFOAM solver.
- `source/DLpisoFoam_auto/`: single-outer AUTO PISO build sharing AUTO sampling and training sources.
- `source/DLbuoyantPimpleFoam_auto/`: buoyant, heat-transfer AUTO solver.
- `source/DLpimpleFoam_oracle/`: oracle diagnostic solver.
- `source/DLpisoFoam_oracle/`: PISO oracle build using the shared oracle implementation.
- `packages/`: `pressure_SM_gridded_shift` and `pressure_SM_mesh_graph` source.

See [source/README.md](source/README.md) for executable ownership, shared oracle
code, PISO/PIMPLE sequencing, and graph-transfer controls. The
[confined](tutorials/DLpimpleFoam_oracle/confined) and
[unconfined](tutorials/DLpimpleFoam_oracle/unconfined) oracle cases generate
their meshes from scripts; no `constant/polyMesh` is shipped with them.

The `testing_oct_2026/.../3by2_300k_heatTransfer/` case uses a Boussinesq
water-like fluid, a 293.15 K inlet, and 313.15 K fixed-temperature cylinders.
It predicts `dp` and the direct `dUCorr` velocity correction. Temperature is
advanced by the CFD enthalpy equation and is not included in the surrogate
features.

## Post-Injection Pressure Corrector

The AUTO solver can optionally run one standard `pEqn.H` pressure correction
after surrogate injection. It assembles a fresh momentum equation from the
injected state, then places the pressure correction before or after the
configured fixed-pressure reconstruction loop:

```text
postSMInjectionReconstruction true;
smPostInjectionPressureCorrector true;
smPostInjectionPressureCorrectorPosition beforeReconstruction;
```

Use `afterReconstruction` to place it after that loop. The pressure corrector
is disabled by default, preserving existing behavior, and enabling it requires
`postSMInjectionReconstruction true`. When the reconstruction count is zero,
the two positions are adjacent and therefore equivalent. Its residuals are
cleared before normal PIMPLE residual control resumes.

## Buoyant Solver Integration

`DLbuoyantPimpleFoam_auto` reads integration controls from
`system/MLSamplingDict`; they are not sampling-only metadata. With
`postSMInjectionReconstruction true`, the predicted reduced pressure is
injected before normal PIMPLE resumes with a fresh momentum equation.
`smPressureForceThenReconstruct true` first performs a no-pressure-solve HbyA
reconstruction using the injected pressure, native buoyancy force, and
pressure-derived **mass** flux. Subsequent fixed-pressure momentum solves use
the same flux reconstruction; they do not solve the pressure equation.

The buoyant solver supports the incompressible solver's integration controls:

| Controls | Behavior / defaults |
|---|---|
| `postSMInjectionReconstruction`, `smPressureForceThenReconstruct` | Enable reconstruction and its initial no-solve force stage; both default `false`. |
| `smUseInputPredictorState` | `true` retains the feature predictor as outer iteration 1, requiring at least two outer iterations when reconstruction is enabled. `false` restores timestep-start U/phi before injection and starts normal PIMPLE at iteration 1. Supports both pressure-only and joint predictions. |
| `smStateReconstructionCorrectors`, `smMinStateReconstructionCorrectors` | Maximum and adaptive minimum fixed-pressure momentum solves, both default 1. A maximum of 0 skips these solves (the optional force stage still runs). |
| `smAdaptiveStateReconstruction`, `smStateReconstructionResidualTolerance` | Stop at the initial U equation residual after at least the minimum count; defaults `true`, `1e-4`. |
| `smStateReconstructionEquationRelaxation` | Explicit equation relaxation, default 1; oracle-compatible mode instead uses `fvSolution`. |
| `smAitkenStateReconstruction`, `smAitkenMinFactor`, `smAitkenMaxFactor` | Volume-weighted Aitken updates of U and phi from the second corrector; defaults `true`, 0.2, 1.8. |
| `smPostInjectionPressureCorrector`, `smPostInjectionPressureCorrectorPosition` | Optional fresh-UEqn native pressure correction before/after reconstruction; defaults `false`, `beforeReconstruction`. Requires reconstruction enabled. |
| `smOracleCompatibleIntegration` | Default `false`. Fresh per-step PIMPLE controller, residual isolation, explicit outer stopping and configured equation relaxation. |
| `smAdaptivePressureCorrectors`, `smMinPressureCorrectors`, `smPressureCorrectorTolerance` | Active-model oracle mode only; defaults `true`, 1, `1e-5`. Checks the first nonorthogonal solve residual of **each current** pressure corrector, not a stale first-corrector value. |
| `smOuterUResidualTolerance`, `smOuterPResidualTolerance` | Absolute first-solve outer residual thresholds, inferred from PIMPLE residual controls, otherwise `5e-4`. Pressure means `p_rgh`, not thermodynamic `p`. |
| `smOuterEnergyResidualTolerance` | Buoyant-specific safeguard: infer the `e`/`h` criterion, otherwise use the U threshold. Both explicit and native oracle outer stopping also require energy convergence. |
| `smUseVelocityExtrapolation`, `smVelocityExtrapolationAlpha` | Defaults `false`, 1. Pressure-only, feature-only transient mode. Uses the last completed velocity increment scaled by alpha and the timestep ratio; first-step increment is zero. Cannot combine with joint output. |
| `applySMPressureConsistencyCorrection`, `applySMPredictedUCorrection` | Legacy direct correction flags, both default `false`. Joint total dU is injected once as `U_time_prev + predicted_dU`; enabling the legacy predicted-U flag does not add it a second time. |

`PIMPLE_SM` overrides in `fvSolution` take effect once the surrogate activates.
The buoyant legacy direct-correction flags retain their `false` defaults;
the incompressible solver defaults them to `true`. Set them explicitly when
comparing integration profiles.
Input/reconstruction residuals are cleared before normal outer convergence is
tested. Warmup does not use SM adaptive stopping. Reconstruction, oracle and
feature-only modes reject steady, transonic and LTS configurations; these are
not interchangeable with incompressible pressure correction.

**Use the correct executable.** The buoyant case has thermodynamic `p` and
`p_rgh` in Pa and `phi` in kg/s. The ML pressure is kinematic reduced pressure
`p_rgh/rhoRef`, converted back on injection. `DLpimpleFoam_auto` instead
requires kinematic `p` and volumetric flux and now rejects the wrong units at
startup with an actionable error. Changing `controlDict.application` alone
does not select a binary invoked directly from the shell; do not relabel
field dimensions to bypass the check.

These controls provide integration parity, not a guarantee of acceleration
or identical CFD trajectories. Temperature/energy and density coupling remain
native to the buoyant solver. Compare equivalent timesteps and include the
extra reconstruction/pressure solves, inference and training in total cost.
Prediction RMSE reductions alone do not measure solver acceleration.

Focused runtime regressions use isolated copies of the water Boussinesq case
and deterministic synthetic pressure/joint predictions:

```bash
export SM_BUOYANT_SOLVER="$(command -v DLbuoyantPimpleFoam_auto)"
export SM_INCOMPRESSIBLE_SOLVER="$(command -v DLpimpleFoam_auto)"
python -m unittest discover -s testing_oct_2026 -p test_buoyant_sm_integration.py -v
```

Run from the repository root in an OpenFOAM 8 / Python environment with NumPy
and the declared `foamlib` test dependency. These are integration regressions,
not a trained-model accuracy or speedup benchmark.

## Adaptive Graph Coarsening

Graph coarsening can keep finer aggregates inside user-defined flow regions
while applying extra coarsening outside them. `graph_levels` is the fine-level
depth; each region is an `(xmin, xmax, ymin, ymax, zmin, zmax)` box containing
native cell centres. The default remains uniform coarsening.

```python
graph_levels = 2
graph_adaptive_coarsening = True
graph_adaptive_extra_levels = 2
graph_adaptive_regions = [
	(0.12, 2.4, 0.02, 0.32, 0.0, 0.5),
	(0.50, 2.8, 0.28, 0.58, 0.0, 0.5),
]
```

Use boxes around obstacles and their downstream wakes. Adaptive settings are
part of the graph signature, so train a new model for the selected regions and
levels; existing uniform-graph weights are not interchangeable.

### Obstacle-Layer Protection

To match the oracle's level-2 graph with 15 protected obstacle layers, set
these options in the case's `python_module.py`:

```python
graph_levels = 2
graph_obstacle_patch = 'obstacle'
graph_obstacle_protection_layers = 15
```

Layer 1 is the patch's owner cells; each additional layer adds face-adjacent
cells. Protected cells remain individual nodes at every coarsening level.
The default layer count is zero. This option cannot be combined with adaptive
box coarsening. Protection settings are part of the graph signature: use a
fresh data folder and train a model for the selected mesh and settings.
Training and inference print the protected-cell and graph-node counts. On a
small or tightly confined mesh, 15 layers may cover most or all cells; select
the layer count based on the geometry and the graph size you can afford.

After installing the graph backend, run its focused protection tests from the
repository root:

```bash
python3 -m unittest discover -s auto/packages/pressure_SM_mesh_graph/_3D/CFD_usable/tests -p 'test_obstacle_layers.py' -v
```

## AUTO PISO Integration

`DLpisoFoam_auto` shares the AUTO sampling and training pipeline, but uses
one momentum predictor followed by direct surrogate injection and inner PISO
pressure corrections. It selects this path by executable name. In
`system/fvSolution`, set `PIMPLE/nOuterCorrectors` to `1` and choose the
maximum number of pressure correctors with `PIMPLE/nCorrectors`. In
`system/MLSamplingDict`, use:

```text
smUseInputPredictorState true;
postSMInjectionReconstruction false;
smPostInjectionPressureCorrector false;
smPressureForceThenReconstruct false;
smOracleCompatibleIntegration false;
smWarmStartStateReconstruction false;
```

`smAdaptivePressureCorrectors`, `smMinPressureCorrectors`, and
`smPressureCorrectorTolerance` control optional early stopping in the inner
loop. The [PISO tutorial](tutorials/DLpisoFoam_auto/graph_shift/pressure_only)
provides a complete pressure-only example. Joint graph models can also return
velocity increments, but require a model trained for that output.

## Graph Slice Images

With `plot_training=True`, mesh-graph training exports every sample in the
chronological training split before optimization begins. Images and NPZ grids
are written to `<data_folder>/training_plots/<pressure_mode>/target_z_slices/`.
The plane is halfway between the sampled obstacle's minimum and maximum Z,
and covers the full X/Y mesh extent, including the coarse far field. Without
obstacle geometry, the domain mid-height is used.

The plots show `dp`, `ddp`, total `dU`, and the configured velocity target.
For `velocity_var_to_predict='ddUCorr'`, that target is
`dU - dUStar - dUCorrPrev`, not total `ddU`. Velocity panels show vector
magnitudes and in-plane arrows; NPZ files preserve all three components.
Pressure is volume-mean-free, in kinematic units (m^2/s^2).

After evaluating the best checkpoint, every validation sample also produces
target/prediction/error images in the separate `validation_z_images/` folder.
Each field's target and prediction share a color scale; velocity errors are
vector-error norms. These native-mesh comparisons precede OpenFOAM boundary
correction and momentum reconstruction. Retraining validation images live
under the corresponding `retrain_NNN/` directory.

Images use inverse-distance interpolation from nearby 3D cell centres at the
requested Z, rather than selecting one exact cell-centre Z level. Separated
convex obstacle cross-sections are masked using nearby sampled boundary
points. This mask is suitable for the cylinders; arbitrary nonconvex solids
require a geometry-aware mask. Raster images are diagnostics, not conservative
mesh transfers. Existing checkpoint diagnostics can be regenerated with the
training module's `--evaluate_only --plot_dir <directory>` options without
retraining.

## Build and Run

Run builds from the repository root:

```bash
wmake auto/source/DLpimpleFoam_auto
wmake auto/source/DLpisoFoam_auto
wmake auto/source/DLbuoyantPimpleFoam_auto
wmake auto/source/DLpimpleFoam_oracle
wmake auto/source/DLpisoFoam_oracle
cd auto/tutorials/DLpimpleFoam_auto/graph_shift/pressure_only
bash Allrun
```

For AUTO PISO, use `auto/tutorials/DLpisoFoam_auto/graph_shift/pressure_only`
and run `bash Allrun` there. It shares the AUTO sampling/training pipeline
but uses one momentum predictor followed by inner PISO pressure corrections.

The AUTO case's `python_module.py` selects the model and features;
`system/MLSamplingDict` controls sampling, training, and activation. Start with
the [tutorial index](../tutorials/README.md) and [repository setup guide](../README.md).