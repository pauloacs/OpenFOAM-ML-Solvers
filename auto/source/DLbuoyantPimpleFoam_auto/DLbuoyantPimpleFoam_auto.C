#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "fluidThermo.H"
#include "fluidThermoMomentumTransportModel.H"
#include "fluidThermophysicalTransportModel.H"
#include "radiationModel.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "SurrogateModel.H"
#include "MLSampling.H"
#include "Residuals.H"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <string>

namespace
{
template<class Type>
void clearSMStageResiduals(const Foam::polyMesh& mesh)
{
    const Foam::List<Foam::word> fieldNames =
        Foam::Residuals<Type>::fieldNames(mesh);

    forAll(fieldNames, fieldi)
    {
        const_cast
        <
            Foam::DynamicList<Foam::SolverPerformance<Type>>&
        >
        (
            Foam::Residuals<Type>::field(mesh, fieldNames[fieldi])
        ).clear();
    }
}

template<class Type>
Foam::label smResidualCount
(
    const Foam::polyMesh& mesh,
    const Foam::word& fieldName
)
{
    return Foam::Residuals<Type>::found(mesh, fieldName)
        ? Foam::Residuals<Type>::field(mesh, fieldName).size() : 0;
}

template<class Type>
Foam::scalar smFirstResidualSince
(
    const Foam::polyMesh& mesh,
    const Foam::word& fieldName,
    const Foam::label firstSolve
)
{
    if (!Foam::Residuals<Type>::found(mesh, fieldName))
    {
        return Foam::GREAT;
    }
    const auto& residuals = Foam::Residuals<Type>::field(mesh, fieldName);
    return firstSolve < residuals.size()
        ? Foam::cmptMax(residuals[firstSolve].initialResidual()) : Foam::GREAT;
}

std::string trim(const std::string& value)
{
    const std::string whitespace(" \t\r\n");
    const std::size_t first = value.find_first_not_of(whitespace);
    if (first == std::string::npos)
    {
        return std::string();
    }
    const std::size_t last = value.find_last_not_of(whitespace);
    return value.substr(first, last - first + 1);
}

std::map<std::string, std::string> readPythonConfig(const std::string& path)
{
    std::ifstream input(path);
    if (!input)
    {
        FatalErrorInFunction
            << "Cannot read " << path << exit(FatalError);
    }

    std::map<std::string, std::string> values;
    std::string line;
    while (std::getline(input, line))
    {
        const std::size_t comment = line.find('#');
        line = trim(comment == std::string::npos ? line : line.substr(0, comment));
        if (line.empty())
        {
            continue;
        }
        const std::size_t equals = line.find('=');
        if (equals == std::string::npos)
        {
            continue;
        }
        const std::string key = trim(line.substr(0, equals));
        std::string value = trim(line.substr(equals + 1));
        if (!key.empty() && key.find_first_of(" \t") == std::string::npos)
        {
            values.emplace(key, value);
        }
    }
    return values;
}

std::string configValue
(
    const std::map<std::string, std::string>& values,
    const std::string& key,
    const std::string& fallback
)
{
    const auto found = values.find(key);
    if (found == values.end())
    {
        return fallback;
    }
    std::string value = trim(found->second);
    if (!value.empty() && value.back() == ',')
    {
        value.pop_back();
    }
    if (value.size() >= 2
        && ((value.front() == '\'' && value.back() == '\'')
            || (value.front() == '"' && value.back() == '"')))
    {
        value = value.substr(1, value.size() - 2);
    }
    return value;
}

bool configBool
(
    const std::map<std::string, std::string>& values,
    const std::string& key,
    const bool fallback
)
{
    std::string value = configValue(values, key, fallback ? "True" : "False");
    std::transform(value.begin(), value.end(), value.begin(),
        [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    if (value == "true")
    {
        return true;
    }
    if (value == "false")
    {
        return false;
    }
    FatalErrorInFunction
        << "python_module.py " << key << " must be a literal True or False."
        << exit(FatalError);
    return false;
}

void addInputFields
(
    const std::map<std::string, std::string>& config,
    std::set<std::string>& fields
)
{
    const std::pair<const char*, const char*> inputMap[] =
    {
        {"add_U_input", "U"},
        {"add_dUStar_input", "dUStar"},
        {"add_dUPrev_input", "dUPrev"},
        {"add_ddUPrev_input", "ddUPrev"},
        {"add_ddUStar_input", "ddUStar"},
        {"add_ddUStarDiff_input", "ddUStarDiff"},
        {"add_dpPrev_input", "dpPrev"},
        {"add_pPrev_input", "p_prev"},
        {"add_ddpPrev_input", "ddpPrev"},
        {"add_laplacian_dpPrev_input", "laplaceDpPrev"},
        {"add_uDotGradDpPrev_input", "uDotGradDpPrev"},
        {"add_gradDpPrevMag_input", "gradDpPrevMag"},
        {"include_rAU_input", "rAU"},
        {"include_HbyA_input", "HbyA"},
        {"include_divHbyA_input", "divHbyA"},
        {"include_dHbyA_input", "dHbyA"},
        {"include_dDivHbyA_input", "dDivHbyA"},
        {"add_pressureEqResidualp_input", "pressureEqResidualp"},
        {"add_gradDpPrev_input", "gradDpPrev"},
        {"add_rAUGradDpPrev_input", "rAUGradDpPrev"},
        {"add_divRAUGradDpPrev_input", "divRAUGradDpPrev"},
        {"add_rAUGradpPrev_input", "rAUGradpPrev"},
        {"add_divRAUGradpPrev_input", "divRAUGradpPrev"},
        {"add_divUStar_input", "divUFirstPred"},
        {"add_dUCorrPrev_input", "dUCorrPrev"},
        {"add_ddUCorrPrev_input", "ddUCorrPrev"}
    };

    for (const auto& entry : inputMap)
    {
        if (configBool(config, entry.first, false))
        {
            fields.insert(entry.second);
        }
    }
}
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient buoyant PIMPLE solver with mesh-graph pressure and velocity surrogate support"
    );

    #include "postProcess.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createRhoUfIfPresent.H"

    turbulence->validate();

    const std::map<std::string, std::string> pythonConfig =
        readPythonConfig("python_module.py");
    const bool train = configBool(pythonConfig, "train", true);
    const std::string pressureTarget =
        configValue(pythonConfig, "var_to_predict", "ddp");
    if (pressureTarget != "dp" && pressureTarget != "ddp")
    {
        FatalErrorInFunction
            << "var_to_predict must be 'dp' or 'ddp'."
            << exit(FatalError);
    }
    const bool predictDP = pressureTarget == "dp";
    const bool ddpFullCoarsePrediction =
        configBool(pythonConfig, "ddp_full_coarse_prediction", false);
    const bool predictDU =
        configBool(pythonConfig, "predict_dU_output", false);
    const std::string dataFolder =
        configValue(pythonConfig, "data_folder", "ML_data");

    if (dataFolder.empty())
    {
        FatalErrorInFunction << "data_folder cannot be empty." << exit(FatalError);
    }
    if (configBool(pythonConfig, "predict_ddUCorr_output", false))
    {
        FatalErrorInFunction
            << "Legacy predict_ddUCorr_output is unsupported by the mesh-graph backend."
            << exit(FatalError);
    }

    const scalar surrogatePressureDensity = gAverage(rho);
    if (!std::isfinite(surrogatePressureDensity) || surrogatePressureDensity <= VSMALL)
    {
        FatalErrorInFunction << "Invalid reference density for surrogate pressure scaling."
            << exit(FatalError);
    }
    const dimensionedScalar surrogatePressureDensityDim
    (
        "rhoRef", rho.dimensions(), surrogatePressureDensity
    );

    volScalarField pSM
    (
        IOobject("pSM", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        p_rgh/surrogatePressureDensityDim
    );

    std::set<std::string> sampleFields;
    sampleFields.insert("ddp");
    addInputFields(pythonConfig, sampleFields);
    if (predictDU)
    {
        sampleFields.insert("dU");
        sampleFields.insert("dUPrev");
    }

    SurrogateModel surrogate
    (
        mesh, pSM, U, dUStar, dUStarPrev,
        dpML, dp, dpPrev, p_prev,
        ddUStar, ddUStarPrev, ddpML, ddpPrev, ddp,
        grad_dpPrev, grad_pPrev, laplace_dpPrev,
        uDotGradDpPrev, gradDpPrevMag, ddUStarML,
        divUFirstPred, dUCorrPrev, ddUCorrPrev,
        rAU_ML, HbyA_ML, divHbyA_ML, dHbyA_ML, dDivHbyA_ML,
        rAUGradDpPrev_ML, divRAUGradDpPrev_ML, pressureEqResidualp_ML,
        rAUGradpPrev_ML, divRAUGradpPrev_ML,
        predictDP, ddpFullCoarsePrediction
    );

    std::string sourceScriptDir(SOLVER_SOURCE_DIR);
    IOdictionary mlDict
    (
        IOobject
        (
            "MLSamplingDict", runTime.system(), mesh,
            IOobject::MUST_READ, IOobject::NO_WRITE
        )
    );

    const Switch postSMInjectionReconstruction =
        mlDict.lookupOrDefault<Switch>("postSMInjectionReconstruction", false);
    const Switch smPressureForceThenReconstruct =
        mlDict.lookupOrDefault<Switch>("smPressureForceThenReconstruct", false);
    const Switch smOracleCompatibleIntegration =
        mlDict.lookupOrDefault<Switch>("smOracleCompatibleIntegration", false);
    const Switch smPostInjectionPressureCorrector =
        mlDict.lookupOrDefault<Switch>("smPostInjectionPressureCorrector", false);
    const word smPostInjectionPressureCorrectorPosition =
        mlDict.lookupOrDefault<word>
        (
            "smPostInjectionPressureCorrectorPosition",
            "beforeReconstruction"
        );
    const Switch smUseInputPredictorState =
        mlDict.lookupOrDefault<Switch>("smUseInputPredictorState", true);
    const Switch smUseVelocityExtrapolation =
        mlDict.lookupOrDefault<Switch>("smUseVelocityExtrapolation", false);
    const scalar smVelocityExtrapolationAlpha =
        mlDict.lookupOrDefault<scalar>("smVelocityExtrapolationAlpha", 1.0);
    const label smStateReconstructionCorrectors =
        mlDict.lookupOrDefault<label>("smStateReconstructionCorrectors", 1);
    const label smMinStateReconstructionCorrectors =
        mlDict.lookupOrDefault<label>("smMinStateReconstructionCorrectors", 1);
    const Switch smAdaptiveStateReconstruction =
        mlDict.lookupOrDefault<Switch>("smAdaptiveStateReconstruction", true);
    const scalar smStateReconstructionResidualTolerance =
        mlDict.lookupOrDefault<scalar>
        (
            "smStateReconstructionResidualTolerance",
            1e-4
        );
    const scalar smStateReconstructionEquationRelaxation =
        mlDict.lookupOrDefault<scalar>
        (
            "smStateReconstructionEquationRelaxation",
            1.0
        );
    const Switch smAitkenStateReconstruction =
        mlDict.lookupOrDefault<Switch>("smAitkenStateReconstruction", true);
    const scalar smAitkenMinFactor =
        mlDict.lookupOrDefault<scalar>("smAitkenMinFactor", 0.2);
    const scalar smAitkenMaxFactor =
        mlDict.lookupOrDefault<scalar>("smAitkenMaxFactor", 1.8);
    const Switch applySMPressureConsistencyCorrection =
        mlDict.lookupOrDefault<Switch>
        (
            "applySMPressureConsistencyCorrection",
            false
        );
    const Switch applySMPredictedUCorrection =
        mlDict.lookupOrDefault<Switch>("applySMPredictedUCorrection", false);
    const Switch smAdaptivePressureCorrectors =
        mlDict.lookupOrDefault<Switch>("smAdaptivePressureCorrectors", true);
    const label smMinPressureCorrectors =
        mlDict.lookupOrDefault<label>("smMinPressureCorrectors", 1);
    const scalar smPressureCorrectorTolerance =
        mlDict.lookupOrDefault<scalar>("smPressureCorrectorTolerance", 1e-5);
    const auto outerResidualTolerance =
    [&](const word& fieldName, const scalar fallback) -> scalar
    {
        const dictionary& solverDict = mesh.solutionDict().subDict("PIMPLE");
        for (const word& name :
             {word("outerCorrectorResidualControl"), word("residualControl")})
        {
            if (solverDict.found(name))
            {
                const dictionary& residuals = solverDict.subDict(name);
                if (residuals.found(fieldName))
                {
                    const entry& field = residuals.lookupEntry
                    (
                        fieldName, false, false
                    );
                    return field.isDict()
                        ? field.dict().lookupOrDefault<scalar>
                          (
                              "tolerance", fallback
                          )
                        : readScalar(residuals.lookup(fieldName));
                }
            }
        }
        return fallback;
    };
    scalar smOuterUResidualTolerance = mlDict.lookupOrDefault<scalar>
    (
        "smOuterUResidualTolerance", outerResidualTolerance("U", 5e-4)
    );
    scalar smOuterPResidualTolerance = mlDict.lookupOrDefault<scalar>
    (
        "smOuterPResidualTolerance", outerResidualTolerance("p_rgh", 5e-4)
    );
    scalar smOuterEnergyResidualTolerance = mlDict.lookupOrDefault<scalar>
    (
        "smOuterEnergyResidualTolerance",
        outerResidualTolerance(thermo.he().name(), smOuterUResidualTolerance)
    );

    if (smStateReconstructionCorrectors < 0)
    {
        FatalErrorInFunction
            << "smStateReconstructionCorrectors must be nonnegative."
            << exit(FatalError);
    }
    if
    (
        !std::isfinite(smStateReconstructionResidualTolerance)
     || !std::isfinite(smStateReconstructionEquationRelaxation)
     || !std::isfinite(smAitkenMinFactor)
     || !std::isfinite(smAitkenMaxFactor)
     || smStateReconstructionResidualTolerance <= 0
     || smStateReconstructionEquationRelaxation <= 0
     || smStateReconstructionEquationRelaxation > 1
     || smAitkenMinFactor <= 0
     || smAitkenMaxFactor < smAitkenMinFactor
    )
    {
        FatalErrorInFunction
            << "Invalid SM reconstruction tolerance, equation relaxation, "
            << "or Aitken bounds in system/MLSamplingDict."
            << exit(FatalError);
    }
    if
    (
        smMinPressureCorrectors < 1
     || !std::isfinite(smPressureCorrectorTolerance)
     || smPressureCorrectorTolerance <= 0
     || !std::isfinite(smOuterUResidualTolerance)
     || !std::isfinite(smOuterPResidualTolerance)
     || !std::isfinite(smOuterEnergyResidualTolerance)
     || smOuterUResidualTolerance <= 0
     || smOuterPResidualTolerance <= 0
     || smOuterEnergyResidualTolerance <= 0
    )
    {
        FatalErrorInFunction
            << "SM adaptive pressure correction requires a positive minimum "
            << "corrector count and tolerance."
            << exit(FatalError);
    }
    if
    (
        postSMInjectionReconstruction && smStateReconstructionCorrectors > 0
     && (smMinStateReconstructionCorrectors < 1
      || smMinStateReconstructionCorrectors > smStateReconstructionCorrectors)
    )
    {
        FatalErrorInFunction
            << "smMinStateReconstructionCorrectors must be between 1 and "
            << "smStateReconstructionCorrectors."
            << exit(FatalError);
    }
    if
    (
        !std::isfinite(smVelocityExtrapolationAlpha)
     || smVelocityExtrapolationAlpha < 0
     || (smUseVelocityExtrapolation && (LTS || predictDU || smUseInputPredictorState))
    )
    {
        FatalErrorInFunction
            << "Velocity extrapolation requires a finite nonnegative alpha, "
            << "transient pressure-only SM, and smUseInputPredictorState=false."
            << exit(FatalError);
    }
    if
    (
        smPostInjectionPressureCorrector
     && (!postSMInjectionReconstruction
      || (smPostInjectionPressureCorrectorPosition != "beforeReconstruction"
       && smPostInjectionPressureCorrectorPosition != "afterReconstruction"))
    )
    {
        FatalErrorInFunction
            << "smPostInjectionPressureCorrector requires "
            << "postSMInjectionReconstruction=true and position "
            << "beforeReconstruction or afterReconstruction."
            << exit(FatalError);
    }
    if
    (
        (postSMInjectionReconstruction || smOracleCompatibleIntegration
      || !smUseInputPredictorState)
     && (LTS || mesh.steady() || pimple.transonic())
    )
    {
        FatalErrorInFunction
            << "Reconstruction/oracle/feature-only integration requires transient, "
            << "non-transonic timestepping without LTS."
            << exit(FatalError);
    }

    Info<< "[SM INTEGRATION] postSMInjectionReconstruction="
        << postSMInjectionReconstruction
        << " pressureForceThenReconstruct=" << smPressureForceThenReconstruct
        << " oracleCompatible=" << smOracleCompatibleIntegration
        << " useInputPredictorState=" << smUseInputPredictorState
        << " reconstructionCorrectors=" << smStateReconstructionCorrectors
        << " minReconstructionCorrectors=" << smMinStateReconstructionCorrectors
        << " postInjectionPressureCorrector=" << smPostInjectionPressureCorrector
        << " pressureCorrectorPosition=" << smPostInjectionPressureCorrectorPosition
        << " velocityExtrapolation=" << smUseVelocityExtrapolation
        << " adaptiveReconstruction=" << smAdaptiveStateReconstruction
        << " reconstructionTolerance="
        << smStateReconstructionResidualTolerance
        << " equationRelaxation="
        << smStateReconstructionEquationRelaxation
        << " aitken=" << smAitkenStateReconstruction
        << " aitkenBounds=(" << smAitkenMinFactor << ' '
        << smAitkenMaxFactor << ')'
        << " pressureConsistencyCorrection="
        << applySMPressureConsistencyCorrection
        << " predictedUCorrection=" << applySMPredictedUCorrection
        << " adaptivePressureCorrectors=" << smAdaptivePressureCorrectors
        << " minPressureCorrectors=" << smMinPressureCorrectors
        << " pressureCorrectorTolerance=" << smPressureCorrectorTolerance
        << nl;

    DataSampler dataSampler
    (
        mesh, UStar, dUStar, dp,
        ddUStar, ddUStarPrev, ddp, dpPrev, ddpPrev,
        grad_dpPrev, grad_pPrev, laplace_dpPrev, uDotGradDpPrev,
        gradDpPrevMag, p_prev, divUFirstPred, ddUCorr,
        dUCorrPrev, ddUCorrPrev, rAU_ML, HbyA_ML,
        divHbyA_ML, dHbyA_ML, dDivHbyA_ML, rAUGradDpPrev_ML,
        divRAUGradDpPrev_ML, pressureEqResidualp_ML, rAUGradpPrev_ML,
        divRAUGradpPrev_ML,
        dataFolder, sourceScriptDir, train,
        mlDict.lookupOrDefault<int>("warmUpSteps", 10),
        mlDict.lookupOrDefault<int>("burstSteps", 10),
        mlDict.lookupOrDefault<int>("burstInterval", 1),
        mlDict.lookupOrDefault<int>("regularInterval", 5),
        mlDict.lookupOrDefault<int>("retrainInterval", 20),
        mlDict.lookupOrDefault<int>("windowFrames", 30),
        mlDict.lookupOrDefault<int>("waitBeforeResampling", 0),
        sampleFields
    );

    bool surrogateActive = false;
    Info<< "DLbuoyantPimpleFoam_auto: train=" << train
        << " predict_dU_output=" << predictDU
        << " pressureMode=" << (predictDP ? "dp" : "ddp")
        << " dataFolder=" << dataFolder
        << " rhoRef=" << surrogatePressureDensity << nl;

    bool pimpleSMOverridesApplied = false;
    const auto applyPimpleSMOverrides = [&]()
    {
        if (pimpleSMOverridesApplied)
        {
            return;
        }
        const dictionary& solution = mesh.solutionDict();
        dictionary& controls =
            const_cast<dictionary&>(solution.subDict("PIMPLE"));
        if (solution.found("PIMPLE_SM"))
        {
            const dictionary& overrides = solution.subDict("PIMPLE_SM");
            forAllConstIter(dictionary, overrides, iter)
            {
                controls.set(iter().clone(controls).ptr());
            }
            Info<< "[SM INTEGRATION] applying PIMPLE_SM overrides." << nl;
        }
        smOuterUResidualTolerance = mlDict.lookupOrDefault<scalar>
        (
            "smOuterUResidualTolerance", outerResidualTolerance("U", 5e-4)
        );
        smOuterPResidualTolerance = mlDict.lookupOrDefault<scalar>
        (
            "smOuterPResidualTolerance", outerResidualTolerance("p_rgh", 5e-4)
        );
        smOuterEnergyResidualTolerance = mlDict.lookupOrDefault<scalar>
        (
            "smOuterEnergyResidualTolerance",
            outerResidualTolerance(thermo.he().name(), smOuterUResidualTolerance)
        );
        if
        (
            !std::isfinite(smOuterUResidualTolerance)
         || !std::isfinite(smOuterPResidualTolerance)
         || !std::isfinite(smOuterEnergyResidualTolerance)
         || smOuterUResidualTolerance <= 0
         || smOuterPResidualTolerance <= 0
         || smOuterEnergyResidualTolerance <= 0
        )
        {
            FatalErrorInFunction
                << "SM outer residual tolerances must be finite and positive."
                << exit(FatalError);
        }
        if (smOracleCompatibleIntegration)
        {
            dictionary residuals
            (
                controls.subOrEmptyDict("outerCorrectorResidualControl")
            );
            const word fields[] = {"U", "p_rgh", thermo.he().name()};
            const scalar tolerances[] =
            {
                smOuterUResidualTolerance, smOuterPResidualTolerance,
                smOuterEnergyResidualTolerance
            };
            for (label i = 0; i < 3; ++i)
            {
                dictionary criterion;
                criterion.add("tolerance", tolerances[i]);
                criterion.add("relTol", scalar(0));
                residuals.set(fields[i], criterion);
            }
            controls.set("outerCorrectorResidualControl", residuals);
        }
        pimpleSMOverridesApplied = true;
    };
    if (!train)
    {
        surrogate.init();
        surrogateActive = true;
        applyPimpleSMOverrides();
    }

    bool smHaveVelocityHistory = false;
    scalar smPreviousCompletedDeltaT = 0;
    pimpleControl& runControl = pimple;
    if (smOracleCompatibleIntegration || !smUseInputPredictorState)
    {
        runControl.checkOut();
    }
    while (runControl.run(runTime))
    {
        autoPtr<pimpleControl> stepPimple;
        if (smOracleCompatibleIntegration || !smUseInputPredictorState)
        {
            stepPimple.reset(new pimpleControl(mesh));
            clearSMStageResiduals<vector>(mesh);
            clearSMStageResiduals<scalar>(mesh);
        }
        pimpleControl& pimple =
            (smOracleCompatibleIntegration || !smUseInputPredictorState)
          ? stepPimple()
          : runControl;

        #include "readDyMControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        if
        (
            surrogateActive && postSMInjectionReconstruction
         && smUseInputPredictorState
         && pimple.dict().lookupOrDefault<label>("nOuterCorrectors", 1) < 2
        )
        {
            FatalErrorInFunction
                << "Retained-predictor reconstruction requires at least two "
                << "outer correctors; alternatively set "
                << "smUseInputPredictorState=false."
                << exit(FatalError);
        }
        bool timeStepStartCaptured = false;
        bool modelPredicted = false;
        bool injectionPending = false;
        bool smInputOnlyStage = surrogateActive && !smUseInputPredictorState;
        bool smRefreshInjectedPressurePatches = false;
        #include "smReconstructState.H"

        const auto assembleSMMomentum = [&]()
        {
            MRF.correctBoundaryVelocity(U);
            tmp<fvVectorMatrix> equation
            (
                fvm::ddt(rho, U) + fvm::div(phi, U)
              + MRF.DDt(rho, U) + turbulence->divDevTau(U)
             == fvOptions(rho, U)
            );
            if (smOracleCompatibleIntegration)
            {
                equation.ref().relax();
            }
            else
            {
                equation.ref().relax(smStateReconstructionEquationRelaxation);
            }
            fvOptions.constrain(equation.ref());
            return equation;
        };
        const auto injectSMPressure = [&]()
        {
            const volScalarField psip0(psi*p);
            p_rgh = surrogatePressureDensityDim*pSM;
            p = p_rgh + rho*gh + pRef;
            pressureControl.limit(p);
            thermo.correctRho(psi*p - psip0);
            rho = thermo.rho();
            p_rgh = p - rho*gh - pRef;
            delta_p_rgh = dimensionedScalar
            (
                "zero", delta_p_rgh.dimensions(), 0
            );
            injectionPending = false;
            if (thermo.dpdt())
            {
                dpdt = fvc::ddt(p);
                if (mesh.moving())
                {
                    dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
                }
            }
        };
        surfaceScalarField phiTimePrev
        (
            IOobject
            (
                "phiTimePrev",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            phi
        );
        while (smInputOnlyStage || pimple.loop())
        {
            int counter = 0;
            const label firstUSolve = smResidualCount<vector>(mesh, U.name());
            const label firstEnergySolve =
                smResidualCount<scalar>(mesh, thermo.he().name());
            const bool firstPredictionStage =
                smInputOnlyStage || pimple.firstPimpleIter();
            if (firstPredictionStage && !timeStepStartCaptured)
            {
                U_time_prev = U;
                UStarPrev = U;
                p_prev = p_rgh/surrogatePressureDensityDim;
                phiTimePrev = phi;
                dpPrev = dp;
                ddpPrev = ddp;
                dUStarPrev = dUStar;
                ddUStarPrev = ddUStar;
                dUCorrPrev = dUCorr;
                ddUCorrPrev = ddUCorr;
                HbyA_ML_prev = HbyA_ML;
                divHbyA_ML_prev = divHbyA_ML;
                delta_p_rgh = dimensionedScalar("zero", delta_p_rgh.dimensions(), 0);
                timeStepStartCaptured = true;
            }

            if (smInputOnlyStage)
            {
                Info<< "[SM INPUT] feature-only UEqn before normal PIMPLE." << nl;
            }
            if (smRefreshInjectedPressurePatches)
            {
                tmp<fvVectorMatrix> momentum = assembleSMMomentum();
                reconstructSMState(momentum.ref(), false);
                smRefreshInjectedPressurePatches = false;
                Info<< "[SM DIRECT] injected pressure patches constrained "
                    << "before fresh momentum solve." << nl;
            }
            #include "UEqn.H"
            if (!smInputOnlyStage)
            {
                #include "EEqn.H"
            }

            if (firstPredictionStage && !modelPredicted)
            {
                UStar = U;
                dUStar = UStar - U_time_prev;
                ddUStar = dUStar - dUStarPrev;
                divUFirstPred = fvc::div(UStar);

                rAU_ML = surrogatePressureDensityDim/UEqn.A();
                HbyA_ML = (rAU_ML/surrogatePressureDensityDim)*UEqn.H();
                divHbyA_ML = fvc::div(HbyA_ML);
                dHbyA_ML = HbyA_ML - HbyA_ML_prev;
                dDivHbyA_ML = divHbyA_ML - divHbyA_ML_prev;
                grad_dpPrev = fvc::grad(dpPrev);
                grad_pPrev = fvc::grad(p_prev);
                laplace_dpPrev = fvc::laplacian(dpPrev);
                uDotGradDpPrev = UStar & grad_dpPrev;
                gradDpPrevMag = mag(grad_dpPrev);
                rAUGradDpPrev_ML = rAU_ML*grad_dpPrev;
                divRAUGradDpPrev_ML = fvc::laplacian(rAU_ML, dpPrev);
                pressureEqResidualp_ML =
                    divHbyA_ML - divRAUGradDpPrev_ML;
                rAUGradpPrev_ML = rAU_ML*grad_pPrev;
                divRAUGradpPrev_ML = fvc::laplacian(rAU_ML, p_prev);

                if (surrogateActive)
                {
                    surrogate.predict();
                    pSM = p_prev + dpML;
                    delta_p_rgh = surrogatePressureDensityDim*dpML;
                    injectionPending = true;

                    if (smInputOnlyStage)
                    {
                        U = U_time_prev;
                        phi = phiTimePrev;
                        K = 0.5*magSqr(U);
                        Info<< "[SM INPUT] restored timestep-start U/phi; "
                            << "frozen predictor inputs retained." << nl;
                    }
                    if (predictDU || smUseVelocityExtrapolation)
                    {
                        volVectorField velocityIncrement
                        (
                            IOobject
                            (
                                "smVelocityIncrement", runTime.timeName(), mesh,
                                IOobject::NO_READ, IOobject::NO_WRITE, false
                            ),
                            ddUStarML
                        );
                        if (smUseVelocityExtrapolation)
                        {
                            const scalar factor = smHaveVelocityHistory
                                ? smVelocityExtrapolationAlpha
                                 *runTime.deltaTValue()/smPreviousCompletedDeltaT
                                : 0;
                            velocityIncrement = factor*(dUStarPrev + dUCorrPrev);
                            Info<< "[SM VELOCITY EXTRAPOLATION] factor="
                                << factor << " history=" << smHaveVelocityHistory
                                << nl;
                        }
                        U = U_time_prev + velocityIncrement;
                        Info<< "[SM JOINT] timestep-start U + velocity increment; "
                            << "UError=" << gMax
                            (
                                mag
                                (
                                    U.primitiveField()
                                  - U_time_prev.primitiveField()
                                  - velocityIncrement.primitiveField()
                                )
                            ) << nl;
                        U.correctBoundaryConditions();
                        fvOptions.correct(U);
                        K = 0.5*magSqr(U);
                        phi = fvc::flux(rho*U);
                        MRF.makeRelative(fvc::interpolate(rho), phi);
                        fvc::makeRelative(phi, rho, U);
                    }
                    modelPredicted = true;
                    scalarField absPressureIncrement(dpML.internalField().size());
                    scalarField velocityIncrementMagnitude
                    (
                        ddUStarML.internalField().size()
                    );
                    forAll(absPressureIncrement, celli)
                    {
                        absPressureIncrement[celli] = mag(dpML[celli]);
                        velocityIncrementMagnitude[celli] =
                            mag(ddUStarML[celli]);
                    }
                    Info<< "[SM] prediction injected: max|dp|="
                        << gMax(absPressureIncrement) << " max|dU|="
                        << gMax(velocityIncrementMagnitude) << nl;

                        if (postSMInjectionReconstruction)
                        {
                            injectSMPressure();
                            {
                                tmp<fvVectorMatrix> injectionMomentum =
                                    assembleSMMomentum();
                                reconstructSMState
                                (
                                    injectionMomentum.ref(),
                                    smPressureForceThenReconstruct
                                );
                            }
                            if (smPressureForceThenReconstruct)
                            {
                                Info<< "[SM INTEGRATION] no-solve HbyA pressure "
                                    << "force before reconstruction." << nl;
                            }

                            const auto runPostInjectionPressureCorrector = [&]()
                            {
                                tmp<fvVectorMatrix> tUEqn = assembleSMMomentum();
                                fvVectorMatrix& UEqn = tUEqn.ref();
                                int counter = 0;
                                scalar smLastPressureResidual = GREAT;
                                #include "pEqn.H"
                                Info<< "[SM POST-INJECTION pEqn] position="
                                    << smPostInjectionPressureCorrectorPosition
                                    << " completed using a fresh UEqn." << nl;
                            };
                            if
                            (
                                smPostInjectionPressureCorrector
                             && smPostInjectionPressureCorrectorPosition
                                == "beforeReconstruction"
                            )
                            {
                                runPostInjectionPressureCorrector();
                            }

                            label reconstructionCorrectorsUsed = 0;
                            scalar reconstructionResidual = GREAT;
                            scalar aitkenFactor = min
                            (
                                smAitkenMaxFactor, max(smAitkenMinFactor, 1.1)
                            );
                            bool havePreviousReconstructionDelta = false;
                            volVectorField previousReconstructionDelta
                            (
                                IOobject
                                (
                                    "smPreviousReconstructionDelta",
                                    runTime.timeName(),
                                    mesh,
                                    IOobject::NO_READ,
                                    IOobject::NO_WRITE,
                                    false
                                ),
                                U - U
                            );

                            for
                            (
                                label reconstructionCorr = 0;
                                reconstructionCorr < smStateReconstructionCorrectors;
                                ++reconstructionCorr
                            )
                            {
                                const volVectorField reconstructionUBefore
                                (
                                    IOobject
                                    (
                                        "smReconstructionUBefore",
                                        runTime.timeName(),
                                        mesh,
                                        IOobject::NO_READ,
                                        IOobject::NO_WRITE,
                                        false
                                    ),
                                    U
                                );
                                const surfaceScalarField reconstructionPhiBefore
                                (
                                    IOobject
                                    (
                                        "smReconstructionPhiBefore",
                                        runTime.timeName(),
                                        mesh,
                                        IOobject::NO_READ,
                                        IOobject::NO_WRITE,
                                        false
                                    ),
                                    phi
                                );

                                tmp<fvVectorMatrix> reconstructionMomentum =
                                    assembleSMMomentum();
                                fvVectorMatrix& reconstructionUEqn =
                                    reconstructionMomentum.ref();
                                reconstructSMState(reconstructionUEqn, false);
                                const SolverPerformance<vector> reconstructionPerf =
                                    solve
                                    (
                                        reconstructionUEqn
                                     ==
                                        fvc::reconstruct
                                        (
                                            (
                                                -ghf*fvc::snGrad(rho)
                                                -fvc::snGrad(p_rgh)
                                            )*mesh.magSf()
                                        )
                                    );
                                reconstructionResidual = cmptMax
                                (
                                    reconstructionPerf.initialResidual()
                                );
                                reconstructionCorrectorsUsed =
                                    reconstructionCorr + 1;
                                fvOptions.correct(U);
                                reconstructSMState(reconstructionUEqn, true);

                                const volVectorField reconstructionDelta
                                (
                                    IOobject
                                    (
                                        "smReconstructionDelta",
                                        runTime.timeName(),
                                        mesh,
                                        IOobject::NO_READ,
                                        IOobject::NO_WRITE,
                                        false
                                    ),
                                    U - reconstructionUBefore
                                );
                                if
                                (
                                    smAitkenStateReconstruction
                                 && havePreviousReconstructionDelta
                                )
                                {
                                    const volVectorField deltaChange
                                    (
                                        IOobject
                                        (
                                            "smReconstructionDeltaChange",
                                            runTime.timeName(),
                                            mesh,
                                            IOobject::NO_READ,
                                            IOobject::NO_WRITE,
                                            false
                                        ),
                                        reconstructionDelta
                                      - previousReconstructionDelta
                                    );
                                    scalar numerator = 0;
                                    scalar denominator = 0;
                                    forAll(U, celli)
                                    {
                                        numerator += mesh.V()[celli]
                                            *(previousReconstructionDelta[celli]
                                             & deltaChange[celli]);
                                        denominator += mesh.V()[celli]
                                            *magSqr(deltaChange[celli]);
                                    }
                                    reduce(numerator, sumOp<scalar>());
                                    reduce(denominator, sumOp<scalar>());
                                    if (denominator > VSMALL)
                                    {
                                        aitkenFactor = min
                                        (
                                            smAitkenMaxFactor,
                                            max
                                            (
                                                smAitkenMinFactor,
                                                -aitkenFactor*numerator/denominator
                                            )
                                        );
                                    }
                                    U = reconstructionUBefore
                                      + aitkenFactor*reconstructionDelta;
                                    phi = reconstructionPhiBefore
                                      + aitkenFactor*(phi-reconstructionPhiBefore);
                                    U.correctBoundaryConditions();
                                    fvOptions.correct(U);
                                    K = 0.5*magSqr(U);
                                    fvc::correctRhoUf(rhoUf, rho, U, phi);
                                }
                                previousReconstructionDelta =
                                    reconstructionDelta;
                                havePreviousReconstructionDelta = true;

                                Info<< "[SM RECONSTRUCTION] corrector="
                                    << reconstructionCorrectorsUsed
                                    << " UResidual=" << reconstructionResidual
                                    << " aitkenFactor=" << aitkenFactor << nl;
                                if
                                (
                                    smAdaptiveStateReconstruction
                                 && reconstructionCorrectorsUsed
                                    >= smMinStateReconstructionCorrectors
                                 && reconstructionResidual
                                    < smStateReconstructionResidualTolerance
                                )
                                {
                                    break;
                                }
                            }

                            if (reconstructionCorrectorsUsed)
                            {
                                Info<< "[SM RECONSTRUCTION] completed with "
                                    << reconstructionCorrectorsUsed
                                    << " fixed-pressure correctors; final U residual="
                                    << reconstructionResidual << nl;
                            }
                            else
                            {
                                Info<< "[SM RECONSTRUCTION] zero fixed-pressure "
                                    << "correctors; no reconstruction UEqn solve."
                                    << nl;
                            }
                            if
                            (
                                smPostInjectionPressureCorrector
                             && smPostInjectionPressureCorrectorPosition
                                == "afterReconstruction"
                            )
                            {
                                runPostInjectionPressureCorrector();
                            }
                            clearSMStageResiduals<vector>(mesh);
                            clearSMStageResiduals<scalar>(mesh);
                            Info<< "[SM RECONSTRUCTION] residual history cleared; "
                                << "normal PIMPLE resumes with a fresh UEqn." << nl;
                            smInputOnlyStage = false;
                            continue;
                        }

                        if (applySMPressureConsistencyCorrection)
                        {
                            const volScalarField pressureIncrement(delta_p_rgh);
                            const volScalarField pressureForceRAU(1.0/UEqn.A());
                            U -= pressureForceRAU*fvc::grad(pressureIncrement);
                            if (applySMPredictedUCorrection && predictDU)
                            {
                                Info<< "[SM JOINT] predicted velocity increment "
                                    << "already injected; not adding it twice." << nl;
                            }
                            U.correctBoundaryConditions();
                            fvOptions.correct(U);
                            phi = fvc::flux(rho*U);
                            MRF.makeRelative(fvc::interpolate(rho), phi);
                            fvc::makeRelative(phi, rho, U);
                            K = 0.5*magSqr(U);
                        }
                    if (smInputOnlyStage)
                    {
                        injectSMPressure();
                        smRefreshInjectedPressurePatches = true;
                        smInputOnlyStage = false;
                        clearSMStageResiduals<vector>(mesh);
                        clearSMStageResiduals<scalar>(mesh);
                        Info<< "[SM DIRECT] starting fresh normal PIMPLE." << nl;
                        continue;
                    }
                }
            }

            scalar smLastPressureResidual = GREAT;
            scalar smFirstOuterPressureResidual = GREAT;
            while (pimple.correct())
            {
                smLastPressureResidual = GREAT;
                #include "pEqn.H"
                if (counter == 0)
                {
                    smFirstOuterPressureResidual = smLastPressureResidual;
                }
                if (injectionPending)
                {
                    delta_p_rgh = dimensionedScalar
                    (
                        "zero", delta_p_rgh.dimensions(), 0
                    );
                    injectionPending = false;
                }
                ++counter;
                if
                (
                    modelPredicted && smOracleCompatibleIntegration
                 && smAdaptivePressureCorrectors
                 && counter >= smMinPressureCorrectors
                 && smLastPressureResidual < smPressureCorrectorTolerance
                )
                {
                    Info<< "[SM ADAPTIVE PISO] stopping after " << counter
                        << " pressure correctors; initial residual="
                        << smLastPressureResidual << " tolerance="
                        << smPressureCorrectorTolerance << nl;
                    break;
                }
            }

            bool turbulenceCorrected = false;
            if (pimple.turbCorr())
            {
                turbulence->correct();
                thermophysicalTransport->correct();
                turbulenceCorrected = true;
            }
            if (modelPredicted && smOracleCompatibleIntegration)
            {
                const scalar uResidual =
                    smFirstResidualSince<vector>(mesh, U.name(), firstUSolve);
                const scalar energyResidual = smFirstResidualSince<scalar>
                (
                    mesh, thermo.he().name(), firstEnergySolve
                );
                if
                (
                    uResidual < smOuterUResidualTolerance
                 && smFirstOuterPressureResidual < smOuterPResidualTolerance
                 && energyResidual < smOuterEnergyResidualTolerance
                )
                {
                    if (!turbulenceCorrected)
                    {
                        turbulence->correct();
                        thermophysicalTransport->correct();
                    }
                    Info<< "[SM ORACLE COMPAT] outer residual tolerance "
                        << "satisfied: U=" << uResidual
                        << " p_rgh=" << smFirstOuterPressureResidual
                        << " " << thermo.he().name() << "=" << energyResidual
                        << nl;
                    break;
                }
            }
        }

        pSM = p_rgh/surrogatePressureDensityDim;
        dp = pSM - p_prev;
        ddp = dp - dpPrev;
        dUCorr = U - UStar;
        ddUCorr = dUCorr - dUCorrPrev;
        rho = thermo.rho();
        smHaveVelocityHistory = true;
        smPreviousCompletedDeltaT = runTime.deltaTValue();

        const int retrainStatus = dataSampler.update();
        if (retrainStatus == 1 && !surrogateActive)
        {
            surrogate.init();
            surrogateActive = true;
            applyPimpleSMOverrides();
            Info<< "[SM] initial training complete; surrogate activated." << nl;
        }
        else if (retrainStatus == 2 && surrogateActive)
        {
            surrogate.reload();
        }

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;
    return 0;
}
