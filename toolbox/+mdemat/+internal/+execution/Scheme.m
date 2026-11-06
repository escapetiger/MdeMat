classdef Scheme < handle
    %SCHEME Compose validated configuration with explicitly owned observers.
    %   config accepts a scalar struct patch or an explicit Config value.
    %   Struct patches use owner-declared keys; Config values may additionally
    %   declare extension pairs. Owner predicates always govern known keys.
    %   File values are overlaid by explicit values before validation.
    %
    %   Numerical bindings freeze when construction finishes. Config is a
    %   read-only plain snapshot; deferred defaults are checked by the
    %   operation that needs them. Supply missing inputs to a new instance.
    %   Observers and Observation remain independent of numerical bindings.

    properties (Dependent, SetAccess = private)
        Configuration
        Config
        Derived
    end
    properties (SetAccess = private)
        Observers (1, 1) struct = struct('timer', [], 'visualizer', [])
        Observation (1, 1) struct
    end
    properties (Access = private)
        Stored = mdemat.tools.environment.Config()
        Working = mdemat.tools.environment.Config()
        DerivedValues (1, 1) struct = struct()
        WorkingDerived (1, 1) struct = struct()
        RuntimeValues (1, 1) struct = struct()
        IsLoading (1, 1) logical = false
        IsFrozen (1, 1) logical = false
    end

    methods
        function obj = Scheme(options)
            arguments
                options.file {mustBeTextScalar} = ''
                options.config = struct()
                options.observers (1, 1) struct = struct()
                options.observation (1, 1) struct = struct()
            end
            obj.Observation = observationDefaults();
            obj.setup();
            obj.load(file=options.file, config=options.config);
            obj.Observers = observerPatch(obj.Observers, options.observers);
            names = fieldnames(options.observation);
            for index = 1:numel(names)
                obj.setObservation(names{index}, options.observation.(names{index}));
            end
            obj.IsFrozen = true;
        end

        function value = get.Configuration(obj)
            if obj.IsLoading, value = obj.Working; else, value = obj.Stored; end
        end

        function value = get.Config(obj)
            if obj.IsLoading
                value = combine(obj.Working.toStruct(), obj.WorkingDerived);
            else
                value = obj.RuntimeValues;
            end
        end

        function value = get.Derived(obj)
            if obj.IsLoading, value = obj.WorkingDerived; else, value = obj.DerivedValues; end
        end

        function obj = addConfig(obj, name, options)
            %ADDCONFIG Register an owner default, including inherited override.
            arguments
                obj (1, 1) mdemat.internal.execution.Scheme
                name {mustBeTextScalar}
                options.default = []
                options.validator = []
            end
            obj.requireUnfrozen();
            name = char(name);
            if isReserved(name)
                error('mdemat:execution:Scheme:ReservedKey', ...
                    'Key "%s" belongs to derived state or observer resources.', name);
            end
            obj.Stored = obj.Stored.define(name, options.default, options.validator, ...
                override=obj.Stored.has(name));
        end

        function obj = setConfig(obj, name, value)
            %SETCONFIG Update a declared key only while constructing a candidate.
            obj.requireUnfrozen();
            if obj.IsLoading
                obj.Working = obj.Working.set(name, value);
            else
                obj.Stored = obj.Stored.set(name, value);
                obj.refreshSnapshot();
            end
        end

        function value = getConfig(obj, name)
            arguments
                obj (1, 1) mdemat.internal.execution.Scheme
                name {mustBeTextScalar}
            end
            values = obj.Config;
            name = char(name);
            if ~isfield(values, name)
                value = obj.Configuration.get(name);
                return;
            end
            value = values.(name);
        end

        function tf = hasConfig(obj, name)
            arguments
                obj (1, 1) mdemat.internal.execution.Scheme
                name {mustBeTextScalar}
            end
            tf = obj.Configuration.has(name) || isfield(obj.Derived, char(name));
        end

        function validateConfig(obj, names)
            %VALIDATECONFIG Check only keys needed by the requested operation.
            if nargin < 2, obj.Configuration.validate(); else, obj.Configuration.validate(names); end
        end

        function setDerived(obj, name, value)
            %SETDERIVED Update derived workspace without changing model bindings.
            arguments
                obj (1, 1) mdemat.internal.execution.Scheme
                name {mustBeTextScalar}
                value
            end
            name = char(name);
            if ~isvarname(name) || obj.Configuration.has(name) || ...
                    ismember(name, {'timer', 'visualizer', 'analyzer'})
                error('mdemat:execution:Scheme:ReservedKey', ...
                    'Derived key "%s" cannot replace a declared input or resource.', name);
            end
            if obj.IsLoading
                obj.WorkingDerived.(name) = value;
            else
                obj.DerivedValues.(name) = value;
                obj.refreshSnapshot();
            end
        end

        function obj = loadFromFile(obj, file)
            obj.load(file=file);
        end

        function obj = loadFromStruct(obj, config)
            obj.load(config=config);
        end

        function obj = load(obj, options)
            %LOAD Build and validate a candidate before committing its bindings.
            arguments
                obj (1, 1) mdemat.internal.execution.Scheme
                options.file {mustBeTextScalar} = ''
                options.config = struct()
            end
            input = options.config;
            if ~isscalar(input) || ~(isstruct(input) || isa(input, 'mdemat.tools.environment.Config'))
                error('mdemat:execution:Scheme:InvalidConfiguration', ...
                    'config must be a scalar struct patch or Config value.');
            end
            candidate = obj.Stored;
            if isa(input, 'mdemat.tools.environment.Config')
                names = input.keys();
                for index = 1:numel(names)
                    name = names{index};
                    if ~candidate.has(name) && ~ismember(name, {'timer', 'visualizer'})
                        if isReserved(name)
                            error('mdemat:execution:Scheme:ReservedKey', ...
                                'Extension key "%s" is reserved.', name);
                        end
                        [value, predicate] = input.get(name);
                        candidate = candidate.define(name, value, predicate);
                    end
                end
                input = input.toStruct();
            end
            values = struct();
            if strlength(options.file) > 0
                values = mdemat.tools.environment.readConfig(options.file);
            end
            values = combine(values, input);
            observerValues = struct();
            for name = {'timer', 'visualizer'}
                if isfield(values, name{1})
                    observerValues.(name{1}) = values.(name{1});
                    values = rmfield(values, name{1});
                end
            end
            nextObservers = observerPatch(obj.Observers, observerValues);
            if obj.IsFrozen
                if ~isempty(fieldnames(values))
                    obj.requireUnfrozen();
                end
                obj.Observers = nextObservers;
                return;
            end
            if obj.IsLoading
                error('mdemat:execution:Scheme:NestedLoad', ...
                    'Configuration loads cannot recursively update their candidate.');
            end
            candidate = candidate.update(values);
            previousObservation = obj.Observation;
            previousObservers = obj.Observers;
            obj.Working = candidate;
            obj.WorkingDerived = obj.DerivedValues;
            obj.IsLoading = true;
            guard = onCleanup(@() obj.discardCandidate(previousObservation, previousObservers));
            obj.Observers = nextObservers;
            obj.configure();
            mdemat.internal.execution.configureObservation(obj);
            obj.Stored = obj.Working;
            obj.DerivedValues = obj.WorkingDerived;
            obj.IsLoading = false;
            obj.refreshSnapshot();
            clear guard
        end

        function setObserver(obj, name, value)
            %SETOBSERVER Replace one explicit observer without changing numerics.
            arguments
                obj (1, 1) mdemat.internal.execution.Scheme
                name {mustBeTextScalar}
                value
            end
            if ~ismember(char(name), {'timer', 'visualizer'})
                error('mdemat:execution:Scheme:UnknownObserver', 'Unknown observer "%s".', name);
            end
            patch = struct();
            patch.(char(name)) = value;
            obj.Observers = observerPatch(obj.Observers, patch);
        end

        function value = getObserver(obj, name)
            arguments
                obj (1, 1) mdemat.internal.execution.Scheme
                name {mustBeTextScalar}
            end
            name = char(name);
            if ~isfield(obj.Observers, name)
                error('mdemat:execution:Scheme:UnknownObserver', 'Unknown observer "%s".', name);
            end
            value = obj.Observers.(name);
        end

        function setObservation(obj, kind, patch)
            %SETOBSERVATION Shallow-merge only supplied plot/analysis settings.
            arguments
                obj (1, 1) mdemat.internal.execution.Scheme
                kind {mustBeTextScalar}
                patch (1, 1) struct
            end
            switch lower(char(kind))
                case 'plot', kind = 'Plot';
                case 'analysis', kind = 'Analysis';
                otherwise
                    error('mdemat:execution:Scheme:UnknownObservation', ...
                        'Observation kind must be plot or analysis.');
            end
            current = obj.Observation.(kind);
            names = fieldnames(patch);
            for index = 1:numel(names)
                name = names{index};
                if ~isfield(current, name)
                    error('mdemat:execution:Scheme:UnknownObservationKey', ...
                        'Unknown %s setting "%s".', kind, name);
                end
                validateObservation(name, patch.(name));
                current.(name) = patch.(name);
            end
            obj.Observation.(kind) = current;
        end

        function state = initialize(obj, ~)
            error('mdemat:execution:Scheme:NotImplemented', ...
                'initialize is not implemented by %s.', class(obj));
        end

        function state = preStep(obj, ~)
            error('mdemat:execution:Scheme:NotImplemented', ...
                'preStep is not implemented by %s.', class(obj));
        end

        function state = step(obj, ~)
            error('mdemat:execution:Scheme:NotImplemented', ...
                'step is not implemented by %s.', class(obj));
        end

        function state = postStep(obj, ~)
            error('mdemat:execution:Scheme:NotImplemented', ...
                'postStep is not implemented by %s.', class(obj));
        end

        function state = finalize(obj, ~)
            error('mdemat:execution:Scheme:NotImplemented', ...
                'finalize is not implemented by %s.', class(obj));
        end
    end

    methods (Access = protected)
        function prepareObservers(obj, includeAnalysis) %#ok<INUSD>
            % Explicit run/converge operations opt in to observation resources.
            if isempty(obj.Observers.timer)
                obj.setObserver('timer', mdemat.tools.analyze.Timer());
            end
            if ~isempty(obj.Observation.Plot.Modes) && isempty(obj.Observers.visualizer)
                obj.setObserver('visualizer', mdemat.tools.visual.Visualizer());
            end
        end

        function obj = setup(obj)
            obj.addConfig('verbose', default=0, ...
                validator=@(value) isnumeric(value) && isreal(value) && isscalar(value) && ...
                isfinite(value) && value >= 0);
        end

        function obj = configure(obj)
            % Subclasses normalize candidates and check cross-key constraints.
        end
    end

    methods (Access = private)
        function requireUnfrozen(obj)
            if obj.IsFrozen
                error('mdemat:execution:Scheme:FrozenConfiguration', ...
                    'Numerical bindings are frozen. Configure and construct a new Scheme.');
            end
        end

        function refreshSnapshot(obj)
            obj.RuntimeValues = combine(obj.Stored.toStruct(), obj.DerivedValues);
        end

        function discardCandidate(obj, observation, observers)
            if obj.IsLoading
                obj.Observation = observation;
                obj.Observers = observers;
                obj.IsLoading = false;
            end
            obj.Working = mdemat.tools.environment.Config();
            obj.WorkingDerived = struct();
        end
    end
end

function result = combine(result, patch)
names = fieldnames(patch);
for index = 1:numel(names), result.(names{index}) = patch.(names{index}); end
end

function tf = isReserved(name)
tf = ismember(name, {'timer', 'visualizer', 'analyzer', 'bcType', ...
    'xId', 'vId', 'sId', 'titlePrefix', 'M', 'N'});
end

function observers = observerPatch(observers, patch)
names = fieldnames(patch);
for index = 1:numel(names)
    name = names{index};
    value = patch.(name);
    switch name
        case 'timer', expected = 'mdemat.tools.analyze.Timer';
        case 'visualizer', expected = 'mdemat.tools.visual.Visualizer';
        otherwise
            error('mdemat:execution:Scheme:UnknownObserver', 'Unknown observer "%s".', name);
    end
    if ~isempty(value) && (~isa(value, expected) || ~isscalar(value) || ~isvalid(value))
        error('mdemat:execution:Scheme:InvalidObserver', ...
            'Observer "%s" must be [] or a scalar %s.', name, expected);
    end
    observers.(name) = value;
end
end

function settings = observationDefaults()
settings.Plot = struct('PointsPerDimension', 0, 'Components', struct(), ...
    'Exacts', struct(), 'DatasetTypes', struct(), 'FixedData', struct(), ...
    'Times', [], 'Modes', {{}}, 'Title', '');
settings.Analysis = struct('PointsPerDimension', 0, 'NLevels', 2, ...
    'Components', struct(), 'Reductions', struct(), 'Exacts', struct());
end

function validateObservation(name, value)
switch name
    case 'PointsPerDimension'
        valid = isnumeric(value) && isreal(value) && isvector(value) && ...
            ~isempty(value) && all(isfinite(value)) && all(value == fix(value)) && ...
            ((isscalar(value) && value == 0) || all(value > 0));
    case 'NLevels'
        valid = isnumeric(value) && isreal(value) && isscalar(value) && ...
            isfinite(value) && value >= 1 && value == fix(value);
    case {'Components', 'Exacts', 'DatasetTypes', 'FixedData', 'Reductions'}
        valid = isstruct(value) && isscalar(value);
    case 'Times'
        valid = isnumeric(value) && isreal(value) && (isempty(value) || isvector(value)) && ...
            all(isfinite(value), 'all');
    case 'Modes'
        valid = iscell(value) || isstring(value);
    case 'Title'
        valid = (ischar(value) && (isrow(value) || isempty(value))) || ...
            (isstring(value) && isscalar(value) && ~ismissing(value));
    otherwise
        valid = false;
end
if ~valid
    error('mdemat:execution:Scheme:InvalidObservation', ...
        'Observation setting "%s" is invalid.', name);
end
end
