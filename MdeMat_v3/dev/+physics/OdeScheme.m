classdef OdeScheme < physics.Scheme
    % ODESCHEME Base class for ODE integration schemes.
    %
    %   OdeScheme provides the fundamental framework for ordinary
    %   differential equation integration schemes. It defines the common
    %   interface and workflow for all concrete ODE scheme implementations.
    %
    %   The class manages temporal discretization through various ODE
    %   integration methods including backward Euler, BDF, SDIRK, forward
    %   Euler, Heun, and SSPRK schemes. It also supports operator
    %   splitting methods for multi-physics problems.
    %
    % Examples:
    %   % Create an ODE scheme with backward Euler
    %   config.tOdeIntName = 'be';
    %   config.tFinal = 1.0;
    %   scheme = OdeScheme(config);
    %   odeInt = scheme.tDisc;  % Direct access
    %
    %   % Create scheme with operator splitting
    %   config.tSplitName = 'strang';
    %   config.tOdeIntName = {'be', 'fe'};
    %   config.tFinal = 1.0;
    %   scheme = OdeScheme(config);
    %
    % See also:
    %   physics.Scheme, approx.odeint.OdeIntegrator

    properties (Constant)
        REQUIRED_CONFIGS = {'tOdeIntName', 'tFinal'} %< Required configuration fields
    end

    properties
        tDisc %< Temporal discretization object for ODE integration
    end

    properties (Dependent)
        isFullyImplicit %< Check fully implicit temporal discretization
        isBdf %< Check BDF temporal discretization
        isDirk %< Check DIRK temporal discretization
        isExhausted %< Check if timeline is exhausted
    end

    methods
        function obj = OdeScheme(config)
            % ODESCHEME Constructor for OdeScheme.
            %
            %   obj = OdeScheme(config) creates an ODE scheme object with
            %   the specified configuration structure. Automatically sets
            %   up the temporal discretization based on configuration
            %   fields.
            %
            % Inputs:
            %   config - Configuration structure.
            %
            % Outputs:
            %   obj - Constructed OdeScheme object

            %< Call parent constructor (validates required configs)
            obj@physics.Scheme(config);
            
            %< Extract splitting name (optional field)
            if isfield(config, 'tSplitName')
                tSplitName = config.tSplitName;
            else
                tSplitName = [];
            end

            %< Additional input validation beyond parent class
            core.except.assert(isscalar(config.tFinal) && config.tFinal >= 0, ...
                'tFinal must be non-negative scalar');

            %< Validate tOdeIntName based on splitting presence
            if isempty(tSplitName)
                %< No splitting: tOdeIntName must be char/string
                core.except.assert(ischar(config.tOdeIntName) || isstring(config.tOdeIntName), ...
                    'tOdeIntName must be char or string when no splitting is used');
                odeIntNames = {char(config.tOdeIntName)};
            else
                %< With splitting: tOdeIntName must be cell array of length >= 2
                core.except.assert(ischar(tSplitName) || isstring(tSplitName), ...
                    'tSplitName must be char or string');
                
                core.except.assert(iscell(config.tOdeIntName), ...
                    'tOdeIntName must be cell array when splitting is used');
                
                core.except.assert(length(config.tOdeIntName) >= 2, ...
                    'tOdeIntName must have at least 2 elements when splitting is used');
                
                %< Validate cell array contents
                validCellContents = all(cellfun(@(x) ischar(x) || isstring(x), config.tOdeIntName));
                core.except.assert(validCellContents, ...
                    'All elements in tOdeIntName cell array must be char or string');
                
                odeIntNames = cellfun(@char, config.tOdeIntName, 'UniformOutput', false);
            end

            %< Create odeInt objects
            odeInts = cell(size(odeIntNames));
            for iOdeInt = 1:length(odeIntNames) %< Rule 50: use i prefix
                odeIntClass = obj.getIntegratorClass(odeIntNames{iOdeInt});
                odeInts{iOdeInt} = approx.odeint.(odeIntClass)(config.tFinal);
            end

            %< Apply splitting or single odeInt
            if ~isempty(tSplitName)
                splittingClass = obj.getSplittingClass(char(tSplitName));
                obj.tDisc = approx.odeint.(splittingClass)(odeInts);
            else
                obj.tDisc = odeInts{1};
            end
        end

        function msg = echo(obj, runtime)
            % ECHO Generate status message for current time step.
            %
            %   msg = echo(obj, runtime) creates a formatted status message
            %   showing current time, step size, step count, and runtime.
            %
            % Inputs:
            %   obj - The OdeScheme object
            %   runtime - Current runtime in seconds
            %
            % Outputs:
            %   msg - Formatted status message

            %< Input validation
            core.except.assert(isscalar(runtime) && runtime >= 0, ...
                'Runtime must be non-negative scalar');

            msg = sprintf('[%d] now = %.6g; step size = %.6g; runtime = %.2f s', ...
                obj.tDisc.timeline.count, ...
                obj.tDisc.timeline.now, ...
                obj.tDisc.timeline.dt, ...
                runtime);
        end
    end

    methods (Access = private) %< Rule 246: helper methods
        function className = getIntegratorClass(~, odeIntName)
            % GETINTEGRATORCLASS Map odeInt name to class name.

            switch lower(odeIntName)
                case 'be'
                    className = 'BeIntegrator';
                case 'bdf2'
                    className = 'Bdf2Integrator';
                case 'bdf3'
                    className = 'Bdf3Integrator';
                case 'sdirk2'
                    className = 'Sdirk2Integrator';
                case 'sdirk3'
                    className = 'Sdirk3Integrator';
                case 'fe'
                    className = 'FeIntegrator';
                case 'heun'
                    className = 'HeunIntegrator';
                case 'ssprk2'
                    className = 'Ssprk2Integrator';
                case 'ssprk3'
                    className = 'Ssprk3Integrator';
                otherwise
                    core.except.assert(false, ...
                        'Unknown odeInt name: %s', odeIntName);
            end
        end

        function className = getSplittingClass(~, splittingName)
            % GETSPLITTINGCLASS Map splitting name to class name.

            switch lower(splittingName)
                case 'lie_trotter'
                    className = 'LieTrotterIntegrator';
                case 'strang'
                    className = 'StrangIntegrator';
                otherwise
                    core.except.assert(false, ...
                        'Unknown splitting method: %s', splittingName);
            end
        end
    end

    methods %< Dependent property getters
        function tf = get.isFullyImplicit(obj)
            % GET.ISFULLYIMPLICIT True if fully implicit temporal discretization.

            tf = isa(obj.tDisc, 'approx.odeint.BeIntegrator') || ...
                 isa(obj.tDisc, 'approx.odeint.BdfIntegrator') || ...
                 isa(obj.tDisc, 'approx.odeint.DirkIntegrator');
        end

        function tf = get.isBdf(obj)
            % GET.ISBDF True if BDF temporal discretization.

            tf = isa(obj.tDisc, 'approx.odeint.BeIntegrator') || ...
                 isa(obj.tDisc, 'approx.odeint.BdfIntegrator');
        end

        function tf = get.isDirk(obj)
            % GET.ISDIRK True if DIRK temporal discretization.

            tf = isa(obj.tDisc, 'approx.odeint.BeIntegrator') || ...
                 isa(obj.tDisc, 'approx.odeint.DirkIntegrator');
        end

        function tf = get.isExhausted(obj)
            % GET.ISEXHAUSTED True if timeline is exhausted.

            tf = obj.tDisc.timeline.isExhausted;
        end
    end
end