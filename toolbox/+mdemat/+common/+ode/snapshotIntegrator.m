function record = snapshotIntegrator(integrator)
%SNAPSHOTINTEGRATOR Capture a completed ODE step as versioned numeric data.
%   The caller supplies an accepted boundary (initial time or after advance).
%   Split factors are captured recursively. History U0, two-time mass data,
%   coefficient tables and timeline queues are retained; stage derivatives,
%   factorization/solver caches and diagnostics are rebuilt on restoration.
%   Coefficient/source/nonlinear callbacks remain explicit caller inputs.
if ~isa(integrator, 'mdemat.common.ode.OdeIntegrator') && ...
        ~isa(integrator, 'mdemat.common.ode.SeparableOdeIntegrator')
    error('mdemat:ode:UnsupportedCheckpoint', 'Only built-in ODE workspaces are supported.');
end
timeline = integrator.Timeline;
time = struct('Final', timeline.Final, 'Now', timeline.Now, 'Next', timeline.Next, ...
    'Count', timeline.Count, 'NSteps', timeline.NSteps, ...
    'StepSizeQueue', timeline.StepSizeQueue, 'OldStepSizeQueue', timeline.OldStepSizeQueue);
method = erase(class(integrator), 'mdemat.common.ode.');
record = struct('format', 'mdemat.odeWorkspace', 'schemaVersion', 1, ...
    'phase', 'completed-step', 'method', method, 'timeline', time, ...
    'nSteps', 0, 'nStages', 0, 'history', {{}}, 'massLhs', [], 'massRhs', [], ...
    'coefficients', struct(), 'settings', struct(), 'linearSolvers', struct(), 'factors', {{}});
if isa(integrator, 'mdemat.common.ode.SeparableOdeIntegrator')
    if ~ismember(method, {'LieTrotterIntegrator', 'StrangIntegrator'})
        error('mdemat:ode:UnsupportedCheckpoint', 'Unregistered split integrator.');
    end
    record.factors = cell(1, integrator.NFactors);
    for index = 1:integrator.NFactors
        record.factors{index} = mdemat.common.ode.snapshotIntegrator(integrator.getFactor(index));
    end
else
    if ~isa(integrator, 'mdemat.common.ode.OdeIntegrator')
        error('mdemat:ode:UnsupportedCheckpoint', 'Only built-in ODE workspaces are supported.');
    end
    record.nSteps = integrator.NSteps; record.nStages = integrator.NStages;
    record.history = integrator.U0;
    record.massLhs = integrator.MLhs; record.massRhs = integrator.MRhs;
    for name = ["alpha", "beta", "A", "b", "c", "AI", "AE", "bI", "bE"]
        if isprop(integrator, name), record.coefficients.(name) = integrator.(name); end
    end
    for name = ["NonlinearTolerance", "NonlinearMaxIterations", "Startup"]
        if isprop(integrator, name), record.settings.(name) = integrator.(name); end
    end
    record.linearSolvers.main = solverSettings(integrator.Solver);
    record.linearSolvers.augment = [];
    if isprop(integrator, 'AugmentSolver')
        record.linearSolvers.augment = solverSettings(integrator.AugmentSolver);
    end
end
% Validate through the same fixed constructor/schema grammar used on load.
mdemat.common.ode.restoreIntegrator(record);
end

function description = solverSettings(solver)
if iscell(solver)
    description = cellfun(@solverSettings, solver, UniformOutput=false);
elseif isempty(solver) && isnumeric(solver)
    description = [];
elseif strcmp(class(solver), 'mdemat.common.linalg.LinearSolver')
    description = struct('nnzTh', solver.NnzTh, 'condTh', solver.CondTh, ...
        'options', solver.KrylovOpts);
else
    error('mdemat:ode:UnsupportedCheckpoint', 'Unregistered linear solver workspace.');
end
end
