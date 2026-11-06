function [state, tables] = runConvergence(scheme, state)
%RUNCONVERGENCE Own refinement/reference policy outside numeric tools.
%   Each field has a supplied-reference or interlevel-difference path.
%   Mixed studies retain the requested number of observations per series.
options = scheme.Observation.Analysis;
tables = [];
if isempty(options.PointsPerDimension) || any(options.PointsPerDimension <= 0)
    return;
end
names = fieldnames(options.Components);
names = names(~cellfun(@(name) isempty(options.Components.(name)), names));
if isempty(names), return; end
count = options.NLevels;
if ~isscalar(count) || ~isfinite(count) || count < 1 || count ~= fix(count)
    error('mdemat:execution:runConvergence:Levels', ...
        'Use a positive integer number of observed levels.');
end
observations = [];
coarse = [];
level = 1;
lastLevel = count;
while level <= lastLevel
    state = scheme.run(state);
    if level == 1
        % Initialization may construct physical U/G/F references. Decide the
        % extra interlevel solve from these actual model-owned references.
        exacts = scheme.Observation.Analysis.Exacts;
        needCoarse = any(cellfun(@(name) ...
            ~isfield(exacts, name) || isempty(exacts.(name)), names));
        lastLevel = count + needCoarse;
    end
    time = [];
    if isprop(scheme, 'TDisc'), time = scheme.TDisc.Timeline.Now; end
    timer = scheme.getObserver('timer');
    [~, interval] = timer.start(record='Analysis');
    analysisGuard = onCleanup(@() timer.cancel(interval));
    records = mdemat.internal.execution.assessState(scheme, state, coarse, time);
    timer.stop(interval=interval);
    clear analysisGuard
    if level > count && ~isempty(records)
        records = records(~strcmp({records.Comparison}, 'exact'));
    end
    if ~isempty(records)
        if isempty(observations)
            observations = records;
        else
            observations = [observations, records]; %#ok<AGROW>
        end
    end
    if level < lastLevel
        coarse = state;
        state = state.refine(1);
    end
    level = level + 1;
end
if isempty(observations), return; end
summary = mdemat.tools.analyze.summarizeConvergence(observations);
tables = mdemat.internal.execution.convergenceTables(summary);
if scheme.Config.verbose > 0
    keys = fieldnames(tables);
    for index = 1:numel(keys), disp(tables.(keys{index})); end
end
end
