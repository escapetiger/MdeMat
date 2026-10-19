classdef DiffusionScheme < physics.MolScheme
    % DIFFUSIONSCHEME Base class for all schemes solving diffusion
    % equations.
    %
    %   DiffusionScheme implements a numerical scheme for solving diffusion
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
            %   obj - The DiffusionScheme object
            %
            % Outputs:
            %   obj - The DiffusionScheme object

            parse@physics.MolScheme(obj);

            %< Validate required fields
            requiredFields = {'diffusion', 'cfl', 'ic', 'bc'};
            [tf, mf] = core.validate.hasFields(obj.config, requiredFields);
            core.except.assert(tf, 'InvalidInput', ...
                'Required configuration %s is missing.', mf);

            %< Validate initial condition
            core.except.assert(isa(obj.config.ic, 'function_handle'), ...
                'InvalidInput', 'ic must be a function handle.');

            %< Validate diffusion
            core.except.assert(ismatrix(obj.config.diffusion) && ...
                isnumeric(obj.config.diffusion), 'InvalidInput', ...
                'diffusion must be a numeric matrix.');

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