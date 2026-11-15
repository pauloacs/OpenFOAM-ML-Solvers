// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    oracleTwinReusePredictorMatrix
  ? oraclePredictorUEqn
  : tmp<fvVectorMatrix>
    (
        new fvVectorMatrix
        (
            fvm::ddt(U) + fvm::div(phi, U)
          + MRF.DDt(U)
          + turbulence->divDevSigma(U)
         ==
            fvOptions(U)
        )
    )
);
fvVectorMatrix& UEqn = tUEqn.ref();

if (!oracleTwinReusePredictorMatrix)
{
    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor() && !oracleTwinSkipMomentumSolve)
    {
        solve(UEqn == -fvc::grad(p));
        fvOptions.correct(U);
    }
}

