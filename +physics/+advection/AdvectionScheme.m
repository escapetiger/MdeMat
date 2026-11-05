classdef AdvectionScheme < physics.MolScheme
    % ADVECTIONSCHEME Base class for all schemes solving advection
    % equations.
    %
    %   AdvectionScheme implements a numerical scheme for solving advection
    %   equations. This base class defines the basic validation and parse
    %   logic of configuration struct.
    %
    % See Also:
    %   physics.OdeScheme

    methods (Access = protected)
        function obj = parse(obj)
            % PARSE Check and parse the contained configuration struct.
            %
            %   obj = parse(obj) checks and parses the contained
            %   configuration struct.
            %
            % Inputs:
            %   obj - The AdvectionScheme object
            %
            % Outputs:
            %   obj - The AdvectionScheme object

            parse@physics.MolScheme(obj);

            %< Validate required fields
            requiredFields = {'advection', 'cfl', 'ic', 'bc'};
            [tf, mf] = core.validate.hasFields(obj.config, requiredFields);
            core.except.assert(tf, 'InvalidInput', ...
                'Required configuration %s is missing.', mf);

            %< Validate initial condition
            core.except.assert(isa(obj.config.ic, 'function_handle'), ...
                'InvalidInput', 'ic must be a function handle.');

            %< Validate advection
            core.except.assert(isvector(obj.config.advection) && ...
                isnumeric(obj.config.advection), 'InvalidInput', ...
                'advection must be a numeric vector.');

            %< Validate CFL number
            core.except.assert(isscalar(obj.config.cfl) && ...
                obj.config.cfl > 0, 'InvalidInput', ...
                'cfl must be positive scalar.');

            %< Determine boundary condition type
            if isempty(obj.config.bc)
                obj.config.bcType = 'periodic';
            else
                obj.config.bcType = 'dirichlet';
            end
        end    
    end
end