const auto reconstructSMState =
[&](fvVectorMatrix& momentum, const bool updateState)
{
    const volScalarField rAU(1.0/momentum.A());
    const surfaceScalarField rhorAUf(fvc::interpolate(rho*rAU));
    const volVectorField HbyA(constrainHbyA(rAU*momentum.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "smPhiHbyA",
        fvc::interpolate(rho)*fvc::flux(HbyA)
      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi, rhoUf))
    );
    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);
    const surfaceScalarField phig
    (
        -rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf()
    );
    phiHbyA += phig;
    constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);
    p_rgh.correctBoundaryConditions();
    p = p_rgh + rho*gh + pRef;
    fvc::makeRelative(phiHbyA, rho, U);

    if (updateState)
    {
        // Match the negative Laplacian and mass-flux sign in pEqn.H.
        // Only evaluate its flux: the predicted pressure stays fixed.
        fvScalarMatrix pressureEqn
        (
            fvc::div(phiHbyA) - fvm::laplacian(rhorAUf, p_rgh)
        );
        pressureEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );
        phi = phiHbyA + pressureEqn.flux();
        U = HbyA
          + rAU*fvc::reconstruct((phig + pressureEqn.flux())/rhorAUf);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
        K = 0.5*magSqr(U);
        fvc::correctRhoUf(rhoUf, rho, U, phi);
    }
};
