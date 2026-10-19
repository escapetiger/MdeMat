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
    % See also:
    %   physics.Scheme, approx.odeint.OdeIntegrator

    properties
        tDisc %< Temporal discretization
        L %< Linear term
        F %< Nonlinear term
        S %< Source term
        M %< Mass term
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

            %< Call parent constructor
            obj@physics.Scheme(config);

            %< Create odeInt objects
            odeInts = cell(size(obj.config.tOdeIntName));
            for iOdeInt = 1:length(obj.config.tOdeIntName)
                cls = obj.getOdeIntClass(obj.config.tOdeIntName{iOdeInt});
                odeInts{iOdeInt} = approx.odeint.(cls)(obj.config.tFinal);
            end

            %< Apply splitting or single odeInt
            if ~isempty(obj.config.tSplitName)
                cls = obj.getSplitClass(obj.config.tSplitName);
                obj.tDisc = approx.odeint.(cls)(odeInts);
            else
                obj.tDisc = odeInts{1};
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
    end

    methods (Access = protected)
        function obj = parse(obj)
            % PARSE Check and parse the contained configuration struct.
            %
            %   obj = parse(obj) checks and parses the contained
            %   configuration struct.
            %
            % Inputs:
            %   obj - The OdeScheme object
            %
            % Outputs:
            %   obj - The OdeScheme object

            parse@physics.Scheme(obj);

            %< Check required fields
            requiredFields = {'tOdeIntName', 'tFinal'};
            [tf, mf] = core.validate.hasFields(obj.config, requiredFields);
            core.except.assert(tf, 'InvalidInput', ...
                'Required configuration %s is missing.', mf);

            %< Extract splitting name (optional field)
            if ~isfield(obj.config, 'tSplitName')
                obj.config.tSplitName = [];
            end

            %< Check final time
            core.except.assert( ...
                isscalar(obj.config.tFinal) && obj.config.tFinal >= 0, ...
                'InvalidInput', 'tFinal must be non-negative scalar.');

            %< Validate tOdeIntName based on splitting presence
            if isempty(obj.config.tSplitName)
                %< No splitting: tOdeIntName must be char/string
                core.except.assert(ischar(obj.config.tOdeIntName) || ...
                    isstring(obj.config.tOdeIntName), 'InvalidInput', ...
                    'tOdeIntName must be text when no splitting is used.');
                obj.config.tOdeIntName = {char(obj.config.tOdeIntName)};
            else
                %< With splitting: tOdeIntName must be cell array of length >= 2
                core.except.assert(ischar(obj.config.tSplitName) || ...
                    isstring(obj.config.tSplitName), 'InvalidInput', ...
                    'tSplitName must be text.');
                
                core.except.assert(iscell(obj.config.tOdeIntName) || ...
                    length(obj.config.tOdeIntName) >= 2, 'InvalidInput', ...
                    ['tOdeIntName must be cell array that ' ...
                    'has at least 2 elements when splitting is used.']);
                
                %< Validate cell array contents
                tf = all(cellfun(@(x) ischar(x) || isstring(x), obj.config.tOdeIntName));
                core.except.assert(tf, 'InvalidInput', ...
                    'All elements in tOdeIntName cell array must be text.');
                
                obj.config.odeIntName = cellfun(@char, obj.config.tOdeIntName, 'Un', 0);
            end
        end    
    end

    methods (Static) %< Rule 246: helper methods
        function className = getOdeIntClass(odeIntName)
            % GETODEINTCLASS Map odeInt name to class name.

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
                    core.except.assert(false, 'InvalidInput', ...
                        'Unknown odeInt name: %s', odeIntName);
            end
        end

        function className = getSplitClass(splitName)
            % GETSPLITCLASS Map splitting name to class name.

            switch lower(splitName)
                case 'lie_trotter'
                    className = 'LieTrotterIntegrator';
                case 'strang'
                    className = 'StrangIntegrator';
                otherwise
                    core.except.assert(false, 'InvalidInput', ...
                        'Unknown splitting method: %s', splitName);
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