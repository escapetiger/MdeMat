classdef VlasovScheme < mdemat.internal.execution.MolScheme
    % VLASOVSCHEME Base class for all schemes solving Vlasov equations.
    %
    %   VlasovScheme implements a numerical scheme for solving Vlasov
    %   equations. This base class defines the basic validation and parse
    %   logic of configuration struct.
    %
    % See also:
    %   mdemat.internal.execution.MolScheme

    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options for Vlasov schemes.
            
            %< Call parent setup
            obj = setup@mdemat.internal.execution.MolScheme(obj);
            mdemat.internal.execution.defineStudyInputs(obj,'plasma');
            
            %< Add Vlasov-specific required configuration options
            obj.addConfig('ic', default=struct(), ...
                validator=@(x) isstruct(x) || isa(x, 'function_handle'));
            obj.addConfig('bc', default=[], ...
                validator=@(x) true);
            obj.addConfig('cfl', default=0.5, ...
                validator=@(x) isempty(x) || (isscalar(x) && x > 0));
            obj.addConfig('dt', default=[], ...
                validator=@(x) isempty(x) || (isscalar(x) && x > 0));
            obj.addConfig('T0', default=1, ...
                validator=@(x) isnumeric(x) || isa(x, 'function_handle'));
            obj.addConfig('rhoi', default=@(x) ones(1, size(x, 2)), ...
                validator=@(x) isnumeric(x) || isa(x, 'function_handle'));
            obj.addConfig('eId', default='Vlasov', ...
                validator=@(x) (ischar(x) || isstring(x)));
            obj.addConfig('errorTypes', default={'L1', 'L2', 'max'}, ...
                validator=@(x) iscell(x) && all(cellfun(@(y) ischar(y) || isstring(y), x)));
            obj.addConfig('exact', default=[], ...
                validator=@(x) isstruct(x) || isa(x, 'function_handle') || isempty(x));
            obj.addConfig('components', default=struct('D', [], 'P', [], 'E', [], 'F', 1), ...
                validator=@(x) isstruct(x) && all(structfun(@(y) isempty(y) || isnumeric(y) && isscalar(y) && (y > 0), x)));
        end    

        function obj = configure(obj)
            % CONFIGURE Generate Vlasov-specific configuration dependencies.
            
            %< Call parent configure
            obj = configure@mdemat.internal.execution.MolScheme(obj);
            
            %< Determine boundary condition type
            if isempty(obj.Config.bc)
                obj.setDerived('bcType', 'periodic');
            else
                obj.setDerived('bcType', 'dirichlet');
            end

        end

    end
end
