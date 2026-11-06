classdef OdeScheme < mdemat.internal.execution.Scheme
    % ODESCHEME Base class for ODE integration schemes.
    %
    %   OdeScheme provides the fundamental framework for ordinary
    %   differential equation integration schemes. It defines the common
    %   interface and workflow for all concrete ODE scheme implementations.
    %
    %   The class manages temporal discretization through various ODE
    %   integration methods including backward Euler, BDF, SDIRK, forward
    %   Euler, Heun, and SSPRK schemes. It also supports operator splitting
    %   methods for multi-physics problems.
    %
    % See also:
    %   mdemat.internal.execution.Scheme, mdemat.common.ode.OdeIntegrator

    properties
        TDisc % Temporal discretization
        L % Linear term
        F % Nonlinear term
        S % Source term
        M % Mass term
    end

    properties (Dependent)
        IsFullyImplicit % Check fully implicit temporal discretization
        IsBdf % Check BDF temporal discretization
        IsDirk % Check DIRK temporal discretization
        IsSplitting % Check splitting scheme
        IsImex % Check IMEX scheme
    end

    methods
        function obj = OdeScheme(options)
            % ODESCHEME Constructor for OdeScheme.
            %
            %   obj = OdeScheme() creates an ODE scheme with default
            %   configuration.
            %
            %   obj = OdeScheme(file=file) creates an ODE scheme and loads
            %   configuration from the specified filename.
            %
            %   obj = OdeScheme(config=config) creates an ODE scheme with
            %   the specified config configuration.
            %
            %   obj = OdeScheme(file=file, config=config) loads
            %   configuration from file and then overrides with struct
            %   values.

            arguments
                options.file {mustBeTextScalar} = ''
                options.config = struct()
                options.observers (1, 1) struct = struct()
                options.observation (1, 1) struct = struct()
            end

            obj@mdemat.internal.execution.Scheme(file=options.file, config=options.config, ...
                observers=options.observers, observation=options.observation);

            %< Create odeInt objects
            odeInts = cell(size(obj.Config.tOdeIntName));
            for iOdeInt = 1:numel(obj.Config.tOdeIntName)
                cls = obj.getOdeIntClass(obj.Config.tOdeIntName{iOdeInt});
                odeInts{iOdeInt} = mdemat.common.ode.(cls)(obj.Config.tFinal);
            end

            %< Apply splitting or single odeInt
            if ~isempty(obj.Config.tSplitName)
                cls = obj.getSplitClass(obj.Config.tSplitName);
                obj.TDisc = mdemat.common.ode.(cls)(factors=odeInts);
            else
                obj.TDisc = odeInts{1};
            end

            %< Initialize ODE-specific terms
            nOdeInts = numel(odeInts);
            if nOdeInts == 1
                obj.L = [];
                obj.F = [];
                obj.M = [];
                obj.S = [];
            else
                obj.L = cell(1, nOdeInts);
                obj.F = cell(1, nOdeInts);
                obj.S = cell(1, nOdeInts);
                obj.M = cell(1, nOdeInts);
            end
        end

        function tf = get.IsFullyImplicit(obj)
            % GET.ISFULLYIMPLICIT True if fully implicit temporal
            % discretization.

            tf = isa(obj.TDisc, 'mdemat.common.ode.BeIntegrator') || ...
                isa(obj.TDisc, 'mdemat.common.ode.BdfIntegrator') || ...
                isa(obj.TDisc, 'mdemat.common.ode.DirkIntegrator');
        end

        function tf = get.IsBdf(obj)
            % GET.ISBDF True if BDF temporal discretization.

            tf = isa(obj.TDisc, 'mdemat.common.ode.BeIntegrator') || ...
                isa(obj.TDisc, 'mdemat.common.ode.BdfIntegrator');
        end

        function tf = get.IsDirk(obj)
            % GET.ISDIRK True if DIRK temporal discretization.

            tf = isa(obj.TDisc, 'mdemat.common.ode.BeIntegrator') || ...
                isa(obj.TDisc, 'mdemat.common.ode.DirkIntegrator');
        end
    
        function tf = get.IsSplitting(obj)
            % GET.ISSPLITTING True if splitting scheme.

            tf = ~isempty(obj.Config.tSplitName);
        end

        function tf = get.IsImex(obj)
            % GET.ISIMEX True if IMEX scheme.

            tf = isa(obj.TDisc, 'mdemat.common.ode.ImexrkIntegrator');
        end
    end

    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options for ODE schemes.
            
            %< Call parent setup
            obj = setup@mdemat.internal.execution.Scheme(obj);
            
            %< Add ODE-specific configuration options
            obj = obj.addConfig('tOdeIntName', default='be', ...
                validator=@mdemat.internal.execution.OdeScheme.validateOdeIntName);
            obj = obj.addConfig('tFinal', default=0, ...
                validator=@(x) isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x >= 0);
            obj = obj.addConfig('tSplitName', default=[], ...
                validator=@mdemat.internal.execution.OdeScheme.validateSplitName);
        end

        function obj = configure(obj)
            % CONFIGURE Generate ODE-specific configuration dependencies.
            
            %< Call parent configure
            obj = configure@mdemat.internal.execution.Scheme(obj);
            
            obj.validateConfig({'tFinal', 'tOdeIntName', 'tSplitName'});
            names = integratorNames(obj.Config.tOdeIntName);
            if isempty(obj.Config.tSplitName)
                if numel(names) ~= 1
                    error('mdemat:execution:OdeScheme:InvalidSplitting', ...
                        'A single integrator is required without operator splitting.');
                end
                obj.setConfig('tSplitName', []);
            elseif numel(names) < 2
                error('mdemat:execution:OdeScheme:InvalidSplitting', ...
                    'Operator splitting requires at least two integrators.');
            end
            obj.setConfig('tOdeIntName', names);
        end
    end

    methods (Static)
        function tf = validateOdeIntName(value)
            %VALIDATEODEINTNAME Pure predicate, independent of Scheme handles.
            try
                names = integratorNames(value);
                tf = all(cellfun(@mdemat.internal.execution.OdeScheme.isValidOdeIntName, names));
            catch
                tf = false;
            end
        end

        function tf = validateSplitName(value)
            tf = isempty(value) || mdemat.internal.execution.OdeScheme.isValidSplitName(value);
        end

        function cls = getOdeIntClass(name)
            % GETODEINTCLASS Map odeInt name to class name.

            switch lower(name)
                case 'be'
                    cls = 'BeIntegrator';
                case 'bdf2'
                    cls = 'Bdf2Integrator';
                case 'bdf3'
                    cls = 'Bdf3Integrator';
                case 'imp'
                    cls = 'ImpIntegrator';
                case 'sdirk2'
                    cls = 'Sdirk2Integrator';
                case 'sdirk3'
                    cls = 'Sdirk3Integrator';
                case 'fe'
                    cls = 'FeIntegrator';
                case 'heun'
                    cls = 'HeunIntegrator';
                case 'ssprk2'
                    cls = 'Ssprk2Integrator';
                case 'ssprk3'
                    cls = 'Ssprk3Integrator';
                case 'ars111'
                    cls = 'Ars111Integrator';
                case 'ars222'
                    cls = 'Ars222Integrator';
                case 'ars443'
                    cls = 'Ars443Integrator';
                otherwise
                    error('mdemat:execution:OdeScheme:InvalidInput', ...
                        'Unknown odeInt name: %s', name);
            end
        end

        function cls = getSplitClass(name)
            % GETSPLITCLASS Map splitting name to class name.

            switch lower(name)
                case 'lie_trotter'
                    cls = 'LieTrotterIntegrator';
                case 'strang'
                    cls = 'StrangIntegrator';
                otherwise
                    error('mdemat:execution:OdeScheme:InvalidInput', ...
                        'Unknown splitting method: %s', name);
            end
        end

        function tf = isValidOdeIntName(name)
            % ISVALIDODEINTNAME Check if ODE integrator name is valid.
            
            validNames = {'be', 'bdf2', 'bdf3', 'imp', 'sdirk2', 'sdirk3', ...
                         'fe', 'heun', 'ssprk2', 'ssprk3', 'ars111', 'ars222', 'ars443'};
            tf = isTextScalar(name) && any(strcmpi(char(name), validNames));
        end
        
        function tf = isValidSplitName(name)
            % ISVALIDSPLITNAME Check if splitting method name is valid.
            
            validNames = {'lie_trotter', 'strang'};
            tf = isTextScalar(name) && any(strcmpi(char(name), validNames));
        end
    end
end

function names = integratorNames(value)
if ischar(value)
    value = {value};
elseif isstring(value) && isvector(value)
    value = num2cell(value);
end
if ~iscell(value) || isempty(value) || ~isvector(value) || ...
        ~all(cellfun(@isTextScalar, value), 'all')
    error('mdemat:execution:OdeScheme:InvalidInput', ...
        'Integrator names must be a nonempty text vector.');
end
names = reshape(cellfun(@char, value, 'UniformOutput', false), 1, []);
end

function tf = isTextScalar(value)
tf = (ischar(value) && isrow(value) && ~isempty(value)) || ...
    (isstring(value) && isscalar(value) && ~ismissing(value) && strlength(value) > 0);
end
