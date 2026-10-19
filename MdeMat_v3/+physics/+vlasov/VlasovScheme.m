classdef VlasovScheme < physics.MolScheme
    % VLASOVSCHEME Base class for all schemes solving Vlasov equations.
    %
    %   VlasovScheme implements a numerical scheme for solving Vlasov
    %   equations. This base class defines the basic validation and parse
    %   logic of configuration struct.
    %
    % See Also:
    %   physics.MolScheme

    methods (Access = protected)
        function obj = parse(obj)
            % PARSE Check and parse the contained configuration struct.
            %
            %   obj = parse(obj) checks and parses the contained
            %   configuration struct.
            %
            % Inputs:
            %   obj - The VlasovScheme object
            %
            % Outputs:
            %   obj - The VlasovScheme object

            parse@physics.MolScheme(obj);

            %< Validate required fields
            requiredFields = {'ic', 'bc'};
            [tf, mf] = core.validate.hasFields(obj.config, requiredFields);
            core.except.assert(tf, 'InvalidInput', ...
                'Required configuration %s is missing.', mf);

            %< Validate initial condition
            core.except.assert(isstruct(obj.config.ic), ...
                'InvalidInput', 'ic must be a struct.');

            %< Validate CFL number
            if isfield(obj.config, 'cfl')
                core.except.assert(isscalar(obj.config.cfl) && ...
                    obj.config.cfl > 0, 'InvalidInput', ...
                    'cfl must be positive scalar.');
            end

            %< Validate time step size
            if isfield(obj.config, 'dt')
                core.except.assert(isscalar(obj.config.dt) && ...
                    obj.config.dt > 0, 'InvalidInput', ...
                    'Time step size must be positive scalar.');
            end

            %< Determine boundary condition type
            if isempty(obj.config.bc)
                obj.config.bcType = 'periodic';
            else
                obj.config.bcType = 'dirichlet';
            end
        end    
    end
end