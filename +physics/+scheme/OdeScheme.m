classdef OdeScheme < physics.scheme.Scheme
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
    %   physics.scheme.Scheme, approx.odeint.OdeIntegrator

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
            %   configuration from the specified @a filename.
            %
            %   obj = OdeScheme(config=config) creates an ODE scheme with
            %   the specified @a config configuration.
            %
            %   obj = OdeScheme(file=file, config=config) loads
            %   configuration from file and then overrides with struct
            %   values.

            arguments
                options.file {mustBeTextScalar} = ''
                options.config struct = struct()
            end

            obj@physics.scheme.Scheme(file=options.file, config=options.config);

            %< Create odeInt objects
            odeInts = cell(size(obj.Config.tOdeIntName));
            for iOdeInt = 1:length(obj.Config.tOdeIntName)
                cls = obj.getOdeIntClass(obj.Config.tOdeIntName{iOdeInt});
                odeInts{iOdeInt} = approx.odeint.(cls)(obj.Config.tFinal);
            end

            %< Apply splitting or single odeInt
            if ~isempty(obj.Config.tSplitName)
                cls = obj.getSplitClass(obj.Config.tSplitName);
                obj.TDisc = approx.odeint.(cls)(factors=odeInts);
            else
                obj.TDisc = odeInts{1};
            end

            %< Initialize ODE-specific terms
            nOdeInts = length(odeInts);
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

            tf = isa(obj.TDisc, 'approx.odeint.BeIntegrator') || ...
                isa(obj.TDisc, 'approx.odeint.BdfIntegrator') || ...
                isa(obj.TDisc, 'approx.odeint.DirkIntegrator');
        end

        function tf = get.IsBdf(obj)
            % GET.ISBDF True if BDF temporal discretization.

            tf = isa(obj.TDisc, 'approx.odeint.BeIntegrator') || ...
                isa(obj.TDisc, 'approx.odeint.BdfIntegrator');
        end

        function tf = get.IsDirk(obj)
            % GET.ISDIRK True if DIRK temporal discretization.

            tf = isa(obj.TDisc, 'approx.odeint.BeIntegrator') || ...
                isa(obj.TDisc, 'approx.odeint.DirkIntegrator');
        end
    
        function tf = get.IsSplitting(obj)
            % GET.ISSPLITTING True if splitting scheme.

            tf = ~isempty(obj.Config.tSplitName);
        end

        function tf = get.IsImex(obj)
            % GET.ISIMEX True if IMEX scheme.

            tf = isa(obj.TDisc, 'approx.odeint.ImexrkIntegrator');
        end
    end

    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options for ODE schemes.
            
            %< Call parent setup
            obj = setup@physics.scheme.Scheme(obj);
            
            %< Add ODE-specific configuration options
            obj = obj.addConfig('tOdeIntName', default='be', ...
                validator=@(x) obj.validateOdeIntName(x));
            obj = obj.addConfig('tFinal', default=0, ...
                validator=@(x) isscalar(x) && x >= 0);
            obj = obj.addConfig('tSplitName', default=[], ...
                validator=@(x) obj.validateSplitName(x));
        end

        function obj = configure(obj)
            % CONFIGURE Generate ODE-specific configuration dependencies.
            
            %< Call parent configure
            obj = configure@physics.scheme.Scheme(obj);
            
            %< Post-process tOdeIntName based on splitting
            if isempty(obj.Config.tSplitName)
                %< No splitting: ensure tOdeIntName is cell array
                if ~iscell(obj.Config.tOdeIntName)
                    obj.Config.tOdeIntName = {char(obj.Config.tOdeIntName)};
                end
            else
                %< With splitting: validate cell array structure
                core.except.assert(iscell(obj.Config.tOdeIntName) && ...
                    length(obj.Config.tOdeIntName) >= 2, 'InvalidInput', ...
                    ['tOdeIntName must be cell array ' ...
                    'with >= 2 elements when splitting is used.']);
                
                %< Convert all elements to char
                obj.Config.tOdeIntName = cellfun(@char, ...
                    obj.Config.tOdeIntName, 'UniformOutput', false);
            end
        end

        function tf = validateOdeIntName(obj, value)
            % VALIDATEODEINTNAME Validate ODE integrator name(s).
            
            if ischar(value) || isstring(value)
                tf = obj.isValidOdeIntName(char(value));
            elseif iscell(value)
                tf = all(cellfun(@(x) obj.isValidOdeIntName(char(x)), value));
            else
                tf = false;
            end
        end
        
        function tf = validateSplitName(obj, value)
            % VALIDATESPLITNAME Validate splitting method name.
            
            if isempty(value)
                tf = true;
            elseif ischar(value) || isstring(value)
                tf = obj.isValidSplitName(char(value));
            else
                tf = false;
            end
        end
    end

    methods (Static)
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
                    core.except.assert(false, 'InvalidInput', ...
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
                    core.except.assert(false, 'InvalidInput', ...
                        'Unknown splitting method: %s', name);
            end
        end

        function tf = isValidOdeIntName(name)
            % ISVALIDODEINTNAME Check if ODE integrator name is valid.
            
            validNames = {'be', 'bdf2', 'bdf3', 'imp', 'sdirk2', 'sdirk3', ...
                         'fe', 'heun', 'ssprk2', 'ssprk3', 'ars111', 'ars222', 'ars443'};
            tf = ismember(lower(name), validNames);
        end
        
        function tf = isValidSplitName(name)
            % ISVALIDSPLITNAME Check if splitting method name is valid.
            
            validNames = {'lie_trotter', 'strang'};
            tf = ismember(lower(name), validNames);
        end
    end
end