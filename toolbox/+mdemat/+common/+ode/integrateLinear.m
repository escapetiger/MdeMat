function [final, report, history] = integrateLinear(initial, integrator, options)
%INTEGRATELINEAR Advance an explicitly assembled coefficient-space system.
%   Solves d(M*U)/dt=L*U+S with the supplied ODE workspace. A two-time mass
%   callback follows the existing semi-Lagrangian [boundary,transfer] protocol.
%   The workspace is reset and used in place; outputs retain numeric data only.
%   Spatial sampling, mathematical assessment, visualization and files belong
%   to the caller. History is optional and stores flattened coefficient columns.
arguments
    initial {mustBeNumeric, mustBeFinite, mustBeNonempty}
    integrator {mustBeA(integrator, {'mdemat.common.ode.OdeIntegrator', 'mdemat.common.ode.SeparableOdeIntegrator'})}
    options.dt (1, 1) double {mustBeFinite, mustBePositive}
    options.mass = []
    options.linear = []
    options.source = []
    options.storeHistory (1, 1) logical = false
    options.maxSteps (1, 1) double {mustBeNonnegative} = Inf
    options.checkpoint (1, 1) struct = struct()
end
if isfinite(options.maxSteps) && fix(options.maxSteps) ~= options.maxSteps
    error('mdemat:ode:InvalidStepLimit', 'maxSteps must be a nonnegative integer or Inf.');
end
initialShape = size(initial);
values = initial(:);
integrator.reset();
resuming = ~isempty(fieldnames(options.checkpoint));
if resuming
    [values, integrator] = restore(integrator, options.checkpoint, initialShape);
end
times = integrator.Timeline.Now;
if options.storeHistory
    coefficients = values;
else
    coefficients = [];
end
stepSizes = zeros(0, 1);
token = tic;
while ~integrator.Timeline.IsExhausted && numel(stepSizes) < options.maxSteps
    integrator.setTimeStep(dt=options.dt);
    integrator.update(values);
    values = integrator.step(M=options.mass, L=options.linear, S=options.source);
    if numel(values) ~= numel(initial) || any(~isfinite(values(:)))
        error('mdemat:ode:InvalidEvolution', ...
            'Each accepted evolution state must have the original finite coefficient shape.');
    end
    stepSizes(end+1, 1) = integrator.Timeline.StepSize;
    integrator.advance();
    if options.storeHistory
        times(end+1, 1) = integrator.Timeline.Now;
        coefficients(:, end+1) = values;
    end
end
elapsed = toc(token);
final = reshape(values, initialShape);
report = struct('FinalTime', integrator.Timeline.Now, 'Steps', numel(stepSizes), ...
    'StepSizes', stepSizes, 'Seconds', elapsed, 'Integrator', string(class(integrator)), ...
    'Completed', integrator.Timeline.IsExhausted, 'Resumed', resuming, ...
    'Checkpoint', capture(integrator, values, initialShape));
history = struct('Time', times, 'Coefficients', coefficients, ...
    'CoefficientShape', initialShape, 'Stored', options.storeHistory);
end

function record = capture(integrator, values, shape)
record = struct('Type', "mdemat.completedOdeStep", 'Version', 1, ...
    'Integrator', string(class(integrator)), 'Phase', "completed-step", ...
    'State', values, 'CoefficientShape', shape, ...
    'IntegratorState', mdemat.common.ode.snapshotIntegrator(integrator));
end

function [values, integrator] = restore(template, record, shape)
required = {'Type'; 'Version'; 'Integrator'; 'Phase'; 'State'; ...
    'CoefficientShape'; 'IntegratorState'};
if ~isequal(sort(fieldnames(record)), sort(required)) || ...
        ~isequal(record.Type, "mdemat.completedOdeStep") || record.Version ~= 1 || ...
        ~isequal(record.Phase, "completed-step") || ...
        ~isequal(record.Integrator, string(class(template))) || ...
        ~isequal(record.CoefficientShape, shape) || ~isnumeric(record.State) || ...
        numel(record.State) ~= prod(shape) || any(~isfinite(record.State(:)))
    error('mdemat:ode:InvalidCheckpoint', 'Checkpoint schema, phase, method or coefficient shape differs.');
end
integrator = mdemat.common.ode.restoreIntegrator(record.IntegratorState);
if ~strcmp(class(integrator), class(template)) || ...
        integrator.Timeline.Final ~= template.Timeline.Final
    error('mdemat:ode:InvalidCheckpoint', 'Checkpoint method or horizon differs from the supplied workspace.');
end
values = record.State(:);
end
