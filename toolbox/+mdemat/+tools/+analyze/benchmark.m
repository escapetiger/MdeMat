function report = benchmark(operation, options)
%BENCHMARK Measure a caller-selected phase inside the current MATLAB process.
%   report = benchmark(@() operation(), warmup=1, repetitions=5, scope="solve")
%   executes the callback warmup+repetitions times; its outputs are discarded.
%   Callback errors propagate. Side effects/resources belong to the callback.
%   Seconds uses MATLAB tic/toc tokens and excludes process/WSL launch time.
%   No hardware-dependent acceptance threshold is applied.
%   Supply a controlled random stream through the callback when needed; this
%   function never resets it. Environment inventory is collected outside timing.
arguments
    operation (1, 1) function_handle
    options.warmup (1, 1) double {mustBeInteger, mustBeNonnegative, mustBeFinite} = 1
    options.repetitions (1, 1) double {mustBeInteger, mustBePositive, mustBeFinite} = 5
    options.scope {mustBeTextScalar} = 'unspecified phase'
end
for index = 1:options.warmup, operation(); end
seconds = zeros(options.repetitions, 1);
for index = 1:options.repetitions
    token = tic;
    operation();
    seconds(index) = toc(token);
end
environment = mdemat.tools.environment.describe();
processor = getenv('PROCESSOR_IDENTIFIER');
logicalProcessors = str2double(getenv('NUMBER_OF_PROCESSORS'));
processorAvailable = ~isempty(processor);
logicalProcessorsAvailable = isfinite(logicalProcessors) && logicalProcessors > 0;
if ~processorAvailable, processor = 'Unavailable'; end
if ~logicalProcessorsAvailable, logicalProcessors = NaN; end
report = struct('Scope', char(options.scope), 'Seconds', seconds, ...
    'MinimumSeconds', min(seconds), 'MedianSeconds', median(seconds), ...
    'MeanSeconds', mean(seconds), 'StandardDeviationSeconds', std(seconds), ...
    'Warmup', options.warmup, 'Repetitions', options.repetitions, ...
    'Process', 'current MATLAB', 'Computer', computer, 'Version', version, ...
    'Environment', environment, 'ProcessorIdentifier', processor, ...
    'ProcessorAvailable',processorAvailable, ...
    'LogicalProcessorsFromEnvironment',logicalProcessors, ...
    'LogicalProcessorsAvailable',logicalProcessorsAvailable, ...
    'HardwareSource','Process environment variables; missing values are not inferred', ...
    'ThreadControl', 'MATLAB default; affinity/threads were not pinned', ...
    'Isolation', 'In-process repeated measurements; no exclusive-machine reservation');
end
