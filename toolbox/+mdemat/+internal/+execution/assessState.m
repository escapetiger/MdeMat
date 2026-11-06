function observations = assessState(scheme, state, coarse, time)
%ASSESSSTATE Sample declared fields, then call numeric assessment functions.
%   Field names identify the measured DOF quantities. In particular D denotes
%   Hermite coefficients, not a silently reconstructed distribution F. Mode
%   padding requires matching basis identity/normalization. The sequence
%   reduction is the component two-norm AFTER the chosen spatial norm.
arguments
    scheme mdemat.internal.execution.Scheme
    state mdemat.common.fem.SpatialState
    coarse = []
    time = []
end
options = scheme.Observation.Analysis;
observations = emptyObservations();
if isempty(options.PointsPerDimension) || any(options.PointsPerDimension <= 0), return; end
names = fieldnames(options.Components);
for index = 1:numel(names)
    field = names{index}; ids = options.Components.(field);
    if isempty(ids), continue; end
    physical = [];
    hasDofs = isfield(state.Dofs, field) && ~isempty(state.Dofs.(field));
    if ~hasDofs
        if ~ismethod(scheme, 'samplePhysicalAssessment')
            error('mdemat:execution:assessState:MissingField', ...
                'The requested field %s is absent and has no physical sampling adapter.', field);
        end
        physical = scheme.samplePhysicalAssessment(state, coarse, field, time);
        if isempty(physical), continue; end
        C = zeros(0, numel(physical.ComponentIds));
    else
        C = state.Dofs.(field);
    end
    allComponents = ((ischar(ids) && isrow(ids)) || ...
        (isstring(ids) && isscalar(ids) && ~ismissing(ids))) && strcmpi(ids, 'all');
    reductions = struct('x', {{'L1','L2','max'}}, 'v', '');
    if isfield(options.Reductions, field), reductions = options.Reductions.(field); end
    if allComponents
        if ~strcmp(reductions.v, 'sequence')
            error('mdemat:execution:assessState:ComponentReduction', ...
                '''all'' requires an explicit sequence aggregation.');
        end
        ids = 1:size(C, 2);
    end
    if ~isnumeric(ids) || ~isreal(ids) || ~isvector(ids)
        error('mdemat:execution:assessState:Components', ...
            'Select a real numeric component vector or the explicit sequence selector ''all''.');
    end
    ids = double(ids(:).');
    if any(ids < 1 | ids ~= fix(ids) | ~isfinite(ids) | ids > size(C, 2)) || ...
            numel(unique(ids)) ~= numel(ids)
        error('mdemat:execution:assessState:Components', 'Invalid selected field component IDs.');
    end
    hasReference = isfield(options.Exacts, field) && ~isempty(options.Exacts.(field));
    if ~isempty(physical)
        [found, rows] = ismember(ids, physical.ComponentIds);
        if ~all(found)
            error('mdemat:execution:assessState:Components', 'The physical field omits selected IDs.');
        end
        U = physical.Values(rows,:); reference = physical.Reference(rows,:);
        plan = struct('Weights', physical.Weights);
        comparison = physical.Comparison; referenceId = physical.ReferenceId;
    elseif hasReference
        plan = mdemat.common.fem.samplingPlan(state.XDisc, options.PointsPerDimension);
        U = state.XDisc.eval(plan.Nodes, C(:, ids)).';
        reference = mdemat.common.fem.referenceValues( ...
            state.XDisc, options.Exacts.(field), plan.Nodes, time);
        if size(reference, 1) < max(ids)
            error('mdemat:execution:assessState:ReferenceComponents', ...
                'The reference for %s omits selected components.', field);
        end
        reference = reference(ids, :);
        comparison = 'exact'; referenceId = ['supplied:', field];
    else
        if isempty(coarse), continue; end
        if ~isfield(coarse.Dofs, field) || isempty(coarse.Dofs.(field))
            error('mdemat:execution:assessState:MissingCoarseField', ...
                'The coarse state omits field %s.', field);
        end
        verifyModeIdentity(state, coarse, field);
        plan = mdemat.common.fem.nestedSamplingPlan( ...
            state.XDisc, coarse.XDisc, options.PointsPerDimension);
        fineValues = state.XDisc.eval(plan.FineNodes, C(plan.FineDofIndices, :)).';
        reference = coarse.XDisc.eval(plan.CoarseNodes, coarse.Dofs.(field)).';
        if size(reference, 1) > size(fineValues, 1)
            error('mdemat:execution:assessState:ComponentDimension', ...
                'Fine component space cannot be smaller than the coarse space.');
        end
        reference(end+1:size(fineValues, 1), :) = 0;
        U = fineValues(ids, :); reference = reference(ids, :);
        comparison = 'interlevel'; referenceId = ['coarse:', field];
    end
    report = mdemat.tools.analyze.weightedError(U, reference, plan.Weights);
    spatialNorms = cellstr(string(reductions.x));
    for normIndex = 1:numel(spatialNorms)
        spatialNorm = spatialNorms{normIndex};
        switch spatialNorm
            case 'L1', errors = report.ComponentL1;
            case 'L2', errors = report.ComponentL2;
            case {'max','LInf'}, errors = report.ComponentLInf;
            otherwise
                error('mdemat:execution:assessState:SpatialNorm', ...
                    'Unknown spatial norm %s.', spatialNorm);
        end
        if isempty(reductions.v)
            for component = 1:numel(ids)
                observations(end+1) = observation(ids(component), errors(component), 'none'); %#ok<AGROW>
            end
        elseif strcmp(reductions.v, 'sequence')
            observations(end+1) = observation(ids, norm(errors, 2), 'sequence'); %#ok<AGROW>
        else
            error('mdemat:execution:assessState:ComponentNorm', ...
                'Unknown component aggregation %s.', reductions.v);
        end
    end
end
    function record = observation(componentIds, errorValue, componentNorm)
        resolvedIds = componentIds;
        if allComponents, componentIds = 'all'; end
        record = struct('Scale', state.XDisc.Mesh.computeMeasure(), 'Field', field, ...
            'ComponentIds', componentIds, 'SpatialNorm', spatialNorm, ...
            'ComponentNorm', componentNorm, 'Comparison', comparison, ...
            'Reference', referenceId, 'Error', errorValue, ...
            'ResolvedComponentIds', resolvedIds);
    end
end

function records = emptyObservations()
records = struct('Scale', {}, 'Field', {}, 'ComponentIds', {}, 'SpatialNorm', {}, ...
    'ComponentNorm', {}, 'Comparison', {}, 'Reference', {}, 'Error', {}, ...
    'ResolvedComponentIds', {});
end

function verifyModeIdentity(fine, coarse, field)
if ~isprop(fine, 'VDisc') && ~isprop(coarse, 'VDisc'), return; end
if ~isprop(fine, 'VDisc') || ~isprop(coarse, 'VDisc') || ...
        ~strcmp(class(fine.VDisc), class(coarse.VDisc))
    error('mdemat:execution:assessState:ComponentBasis', 'Component representations must agree.');
end
if isa(fine.VDisc, 'mdemat.common.fem.SumSpace')
    switch field
        case 'U', a = fine.VDisc.Lhs; b = coarse.VDisc.Lhs;
        case 'G', a = fine.VDisc.Rhs; b = coarse.VDisc.Rhs;
        case 'F'
            if ~isequaln(fine.VRhsNodes, coarse.VRhsNodes)
                error('mdemat:execution:assessState:ComponentBasis', ...
                    'Kinetic sample fields require common physical velocity nodes.');
            end
            return;
        otherwise, return;
    end
else
    if ~strcmp(field, 'D'), return; end
    a = fine.VDisc; b = coarse.VDisc;
end
if isempty(a) || isempty(b)
    if isempty(a) && isempty(b), return; end
    error('mdemat:execution:assessState:ComponentBasis', 'Both component spaces must be present.');
end
left = mdemat.common.fem.describeSpace(a.Element.Approximator.Basis);
right = mdemat.common.fem.describeSpace(b.Element.Approximator.Basis);
if ~isequaln(withoutDegree(left), withoutDegree(right))
    error('mdemat:execution:assessState:ComponentBasis', ...
        'Mode padding requires the same basis, physical scale and normalization.');
end
end

function value = withoutDegree(value)
if isstruct(value)
    if isfield(value, 'degree'), value = rmfield(value, 'degree'); end
    names = fieldnames(value);
    for element = 1:numel(value)
        for index = 1:numel(names)
            value(element).(names{index}) = withoutDegree(value(element).(names{index}));
        end
    end
elseif iscell(value)
    for index = 1:numel(value), value{index} = withoutDegree(value{index}); end
end
end
