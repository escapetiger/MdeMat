function result = compute(scheme, state, options)
%COMPUTE Execute numerical lifecycle hooks without observers or output files.
%   initialize -> [preStep -> step -> postStep] preserves the original
%   algorithm, including splitting, filters, history and field coupling.
%   The reporting-only finalize hook belongs to the explicit run interface.
%   observe(state,time) may return a caller-defined value at initialization
%   and accepted steps. callback(event,state) is an optional notification.
arguments
    scheme mdemat.internal.execution.MolScheme
    state mdemat.data.variable.State
    options.observe = []
    options.callback = []
    options.maxSteps (1, 1) {mustBePositive, mustBeInteger} = 1000000
end
assert(isempty(options.observe) || isa(options.observe, 'function_handle'), ...
    'MdeMat:Execution:InvalidObserver', 'observe must be empty or a function handle.');
assert(isempty(options.callback) || isa(options.callback, 'function_handle'), ...
    'MdeMat:Execution:InvalidCallback', 'callback must be empty or a function handle.');
clock = tic;
state = scheme.initialize(state);
steps = 0;
times = scheme.TDisc.Timeline.Now;
observations = {};
notify('initialized');
while ~scheme.TDisc.Timeline.IsExhausted
    assert(steps < options.maxSteps, 'MdeMat:Execution:StepLimit', ...
        'The numerical step limit was reached before the final time.');
    oldTime = scheme.TDisc.Timeline.Now;
    state = scheme.preStep(state);
    stepSize = scheme.TDisc.Timeline.StepSize;
    assert(isscalar(stepSize) && isfinite(stepSize) && stepSize > 0, ...
        'MdeMat:Execution:InvalidStep', 'The numerical step size must be finite and positive.');
    state = scheme.step(state);
    state = scheme.postStep(state);
    assert(scheme.TDisc.Timeline.Now > oldTime, 'MdeMat:Execution:NoProgress', ...
        'postStep must advance the accepted numerical time.');
    steps = steps+1;
    times(end+1) = scheme.TDisc.Timeline.Now; %#ok<AGROW>
    notify('accepted');
end
result = struct('schemaVersion', 1, 'state', state, 'times', times, ...
    'observations', {observations}, 'diagnostics', state.History, ...
    'execution', struct('steps', steps, 'elapsedSeconds', toc(clock), ...
        'scheme', class(scheme), 'integrator', class(scheme.TDisc)));

    function notify(phase)
        time = scheme.TDisc.Timeline.Now;
        if ~isempty(options.observe)
            observations{end+1} = options.observe(state, time);
        end
        if ~isempty(options.callback)
            options.callback(struct('phase', phase, 'time', time, 'step', steps), state);
        end
    end
end
