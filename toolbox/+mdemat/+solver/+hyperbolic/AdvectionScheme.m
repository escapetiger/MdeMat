classdef AdvectionScheme < mdemat.internal.execution.MolScheme
    % ADVECTIONSCHEME Base class for all schemes solving advection
    % equations.
    %
    %   AdvectionScheme implements a numerical scheme for solving advection
    %   equations. This base class defines the basic validation and parse
    %   logic of configuration struct.
    %
    % See also:
    %   mdemat.internal.execution.MolScheme

    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options for advection
            % schemes.
            
            %< Call parent setup
            setup@mdemat.internal.execution.MolScheme(obj);

            %< Add advection-specific required configuration options
            obj.addConfig('advection', default=[], ...
                validator=@(x) isvector(x) && isnumeric(x));
            obj.addConfig('cfl', default=0.5, ...
                validator=@(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
            obj.addConfig('ic', default=[], ...
                validator=@(x) isa(x, 'function_handle'));
            obj.addConfig('bc', default=[], ...
                validator=@(x) true);
            obj.addConfig('eId', default='Advection', ...
                validator=@(x) (ischar(x) || isstring(x)));
            obj.addConfig('errorTypes', default={'L1', 'L2', 'max'}, ...
                validator=@(x) iscell(x) && all(cellfun(@(y) ischar(y) || isstring(y), x)));
            obj.addConfig('exact', default=[], ...
                validator=@(x) isa(x, 'function_handle') || isempty(x));
            obj.addConfig('components', default=struct('U', 1), ...
                validator=@(x) isstruct(x) && all(structfun(@(y) isempty(y) || isnumeric(y) && isscalar(y) && (y > 0), x)));
        end    

        function obj = configure(obj)
            % CONFIGURE Generate advection-specific configuration
            % dependencies.
            
            %< Call parent configure
            configure@mdemat.internal.execution.MolScheme(obj);
            
            %< Determine boundary condition type
            if isempty(obj.Config.bc)
                obj.setDerived('bcType', 'periodic');
            else
                obj.setDerived('bcType', 'dirichlet');
            end

        end

    end
end
