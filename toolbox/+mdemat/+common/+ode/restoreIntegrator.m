function integrator = restoreIntegrator(record)
%RESTOREINTEGRATOR Construct a fresh ODE workspace from a completed-step record.
%   Constructors come from an explicit allowlist, never saved class dispatch.
%   The new workspace is fully validated before returning. Accepted solution,
%   BDF/mass histories, step-size history and split-factor phases are retained.
%   OldStepSizeQueue is cleared only as a factorization invalidation signal;
%   the accepted step-size sequence remains in StepSizeQueue.
arguments
    record (1, 1) struct
end
try
    if ~plainTree(record), invalid('ODE checkpoint descriptions must contain plain values.'); end
    expected = {'format', 'schemaVersion', 'phase', 'method', 'timeline', ...
        'nSteps', 'nStages', 'history', 'massLhs', 'massRhs', 'coefficients', 'settings', ...
        'linearSolvers', 'factors'};
    exactFields(record, expected);
    if ~isequal(record.format, 'mdemat.odeWorkspace') || ~isequal(record.schemaVersion, 1) || ...
            ~isequal(record.phase, 'completed-step') || ~ischar(record.method) || ...
            ~isrow(record.method) || ~iscell(record.factors)
        invalid('Invalid completed-step workspace schema.');
    end
    validateTimeline(record.timeline);
    final = record.timeline.Final;
    if ismember(record.method, {'LieTrotterIntegrator', 'StrangIntegrator'})
        if ~isempty(record.history) || ~isempty(record.massLhs) || ~isempty(record.massRhs) || ...
                ~isempty(fieldnames(record.coefficients)) || ~isempty(fieldnames(record.settings)) || ...
                ~isempty(fieldnames(record.linearSolvers)) || ...
                record.nSteps ~= 0 || record.nStages ~= 0 || isempty(record.factors)
            invalid('Invalid split workspace fields.');
        end
        factors = cell(1, numel(record.factors));
        for index = 1:numel(factors)
            factors{index} = mdemat.common.ode.restoreIntegrator(record.factors{index});
            if factors{index}.Timeline.Final ~= final, invalid('Split horizons differ.'); end
        end
        if strcmp(record.method, 'LieTrotterIntegrator')
            integrator = mdemat.common.ode.LieTrotterIntegrator(factors=factors);
        else
            integrator = mdemat.common.ode.StrangIntegrator(factors=factors);
        end
    else
        if ~isempty(record.factors), invalid('An unsplit method cannot contain factors.'); end
        integrator = construct(record.method, final, record.nSteps, record.nStages);
        if integrator.NSteps ~= record.nSteps || integrator.NStages ~= record.nStages || ...
                record.timeline.NSteps ~= integrator.Timeline.NSteps || ...
                ~iscell(record.history) || ~isequal(size(record.history), size(integrator.U0)) || ...
                ~numericTree(record.history) || ~numericTree(record.massLhs) || ~numericTree(record.massRhs)
            invalid('Method/history/mass axes do not match the workspace.');
        end
        names = {'alpha', 'beta', 'A', 'b', 'c', 'AI', 'AE', 'bI', 'bE'};
        names = names(cellfun(@(name) isprop(integrator, name), names));
        exactFields(record.coefficients, names);
        for index = 1:numel(names)
            name = names{index}; value = record.coefficients.(name);
            if ~isnumeric(value) || any(~isfinite(value), 'all')
                invalid('Method coefficients must be finite numeric arrays.');
            end
            prototype = integrator.(name);
            if ismember(name, {'A', 'AI', 'AE', 'b', 'c', 'bI', 'bE'}) && ...
                    ~isempty(prototype) && ~isequal(size(value), size(prototype))
                invalid('Butcher coefficient axes differ from this registered method.');
            end
            integrator.(name) = value;
        end
        settingNames = {'NonlinearTolerance', 'NonlinearMaxIterations', 'Startup'};
        settingNames = settingNames(cellfun(@(name) isprop(integrator, name), settingNames));
        exactFields(record.settings, settingNames);
        for index = 1:numel(settingNames)
            name = settingNames{index}; value = record.settings.(name);
            if strcmp(name, 'Startup')
                if ~((isstring(value) && isscalar(value)) || (ischar(value) && isrow(value))) || ...
                        ~ismember(string(value), ["sdirk3", "legacy"])
                    invalid('Invalid BDF startup policy.');
                end
                integrator.Startup = value;
                continue;
            end
            if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ~isfinite(value) || value <= 0
                invalid('Invalid nonlinear workspace setting.');
            end
            if strcmp(name, 'NonlinearMaxIterations') && value ~= fix(value)
                invalid('The nonlinear iteration limit must be an integer.');
            end
            integrator.(name) = value;
        end
        integrator.U0 = record.history;
        integrator.MLhs = record.massLhs; integrator.MRhs = record.massRhs;
        exactFields(record.linearSolvers, {'main', 'augment'});
        integrator.Solver = restoreSolver(record.linearSolvers.main, integrator.Solver);
        if isprop(integrator, 'AugmentSolver')
            integrator.AugmentSolver = restoreSolver(record.linearSolvers.augment, integrator.AugmentSolver);
        elseif ~isempty(record.linearSolvers.augment)
            invalid('This method does not have an augmented-system solver.');
        end
    end
    if record.timeline.NSteps ~= integrator.Timeline.NSteps
        invalid('Timeline queue length differs from the method.');
    end
    names = fieldnames(record.timeline);
    for index = 1:numel(names), integrator.Timeline.(names{index}) = record.timeline.(names{index}); end
    integrator.Timeline.OldStepSizeQueue(:) = 0;
catch exception
    if startsWith(exception.identifier, 'mdemat:ode:'), rethrow(exception); end
    failure = MException('mdemat:ode:InvalidCheckpoint', 'ODE workspace reconstruction failed.');
    throwAsCaller(addCause(failure, exception));
end
end

function integrator = construct(method, final, steps, stages)
factories = struct('BeIntegrator', @mdemat.common.ode.BeIntegrator, ...
    'FeIntegrator', @mdemat.common.ode.FeIntegrator, 'Bdf2Integrator', @mdemat.common.ode.Bdf2Integrator, ...
    'Bdf3Integrator', @mdemat.common.ode.Bdf3Integrator, 'EmpIntegrator', @mdemat.common.ode.EmpIntegrator, ...
    'HeunIntegrator', @mdemat.common.ode.HeunIntegrator, 'Ssprk2Integrator', @mdemat.common.ode.Ssprk2Integrator, ...
    'Ssprk3Integrator', @mdemat.common.ode.Ssprk3Integrator, 'Exrk4Integrator', @mdemat.common.ode.Exrk4Integrator, ...
    'ImpIntegrator', @mdemat.common.ode.ImpIntegrator, 'Sdirk2Integrator', @mdemat.common.ode.Sdirk2Integrator, ...
    'Sdirk3Integrator', @mdemat.common.ode.Sdirk3Integrator, 'Sdirk4Integrator', @mdemat.common.ode.Sdirk4Integrator, ...
    'Esdirk3Integrator', @mdemat.common.ode.Esdirk3Integrator, 'Ars111Integrator', @mdemat.common.ode.Ars111Integrator, ...
    'Ars222Integrator', @mdemat.common.ode.Ars222Integrator, 'Ars443Integrator', @mdemat.common.ode.Ars443Integrator);
if isfield(factories, method)
    factory = factories.(method); integrator = factory(final);
else
    switch method
        case 'BdfIntegrator', integrator = mdemat.common.ode.BdfIntegrator(steps, final);
        case 'DirkIntegrator', integrator = mdemat.common.ode.DirkIntegrator(stages, final);
        case 'ExrkIntegrator', integrator = mdemat.common.ode.ExrkIntegrator(stages, final);
        case 'ImexrkIntegrator', integrator = mdemat.common.ode.ImexrkIntegrator(stages, final);
        otherwise, error('mdemat:ode:UnsupportedCheckpoint', 'Unregistered ODE method tag: %s.', method);
    end
end
end

function validateTimeline(time)
exactFields(time, {'Final', 'Now', 'Next', 'Count', 'NSteps', 'StepSizeQueue', 'OldStepSizeQueue'});
scalars = [time.Final, time.Now, time.Next, time.Count, time.NSteps];
if ~isnumeric(scalars) || ~isequal(size(scalars), [1, 5]) || ~isreal(scalars) || ...
        any(~isfinite(scalars)) || time.Final < 0 || time.Now < 0 || time.Now > time.Final || ...
        time.Next < time.Now || time.Count < 1 || time.Count ~= fix(time.Count) || ...
        time.NSteps < 1 || time.NSteps ~= fix(time.NSteps)
    invalid('Invalid timeline scalar fields.');
end
for name = ["StepSizeQueue", "OldStepSizeQueue"]
    value = time.(name);
    if ~isnumeric(value) || ~isreal(value) || ~isequal(size(value), [1, time.NSteps+1]) || ...
            any(~isfinite(value) | value < 0)
        invalid('Invalid accepted step-size queue.');
    end
end
end

function valid = numericTree(value)
if isnumeric(value)
    valid = all(isfinite(value), 'all');
elseif iscell(value)
    valid = all(cellfun(@numericTree, value), 'all');
else
    valid = false;
end
end

function exactFields(data, expected)
if ~isstruct(data) || ~isscalar(data) || ~isequal(sort(fieldnames(data)), sort(expected(:)))
    invalid('Checkpoint fields differ from the registered schema.');
end
end

function invalid(message)
error('mdemat:ode:InvalidCheckpoint', '%s', message);
end

function valid = plainTree(value)
allowed = {'double', 'single', 'logical', 'char', 'string', 'int8', 'uint8', ...
    'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64'};
if ismember(class(value), allowed)
    valid = true;
elseif iscell(value)
    valid = all(cellfun(@plainTree, value), 'all');
elseif isstruct(value)
    valid = all(cellfun(@plainTree, struct2cell(value)), 'all');
else
    valid = false;
end
end

function solver = restoreSolver(description, prototype)
if iscell(prototype)
    if ~iscell(description) || ~isequal(size(description), size(prototype))
        invalid('Linear solver stage axes differ.');
    end
    solver = cell(size(prototype));
    for index = 1:numel(solver), solver{index} = restoreSolver(description{index}, prototype{index}); end
elseif isempty(prototype) && isnumeric(prototype)
    if ~isempty(description), invalid('An empty solver workspace has unexpected settings.'); end
    solver = [];
else
    exactFields(description, {'nnzTh', 'condTh', 'options'});
    exactFields(description.options, {'maxIter', 'tol', 'restart'});
    solver = mdemat.common.linalg.LinearSolver(nnzTh=description.nnzTh, condTh=description.condTh, ...
        maxIter=description.options.maxIter, tol=description.options.tol, restart=description.options.restart);
end
end
