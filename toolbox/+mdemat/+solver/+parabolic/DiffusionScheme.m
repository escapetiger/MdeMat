classdef DiffusionScheme < mdemat.internal.execution.MolScheme
    % DIFFUSIONSCHEME Base class for all schemes solving diffusion
    % equations.
    %
    %   DiffusionScheme implements a numerical scheme for solving diffusion
    %   equations. This base class defines the basic validation and parse
    %   logic of configuration struct.
    %
    % See also:
    %   mdemat.internal.execution.MolScheme
    
    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options for diffusion schemes.
            
            %< Call parent setup
            obj = setup@mdemat.internal.execution.MolScheme(obj);
            
            %< Add diffusion-specific required configuration options
            obj.addConfig('diffusion', default=[], ...
                validator=@(x) ismatrix(x) && isnumeric(x));
            obj.addConfig('cfl', default=0.5, ...
                validator=@(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
            obj.addConfig('ic', default=[], ...
                validator=@(x) isa(x, 'function_handle'));
            obj.addConfig('bc', default=[], ...
                validator=@(x) true);
            obj.addConfig('eId', default='Diffusion', ...
                validator=@(x) (ischar(x) || isstring(x)));
            obj.addConfig('errorTypes', default={'L1', 'L2', 'max'}, ...
                validator=@(x) iscell(x) && all(cellfun(@(y) ischar(y) || isstring(y), x)));
            obj.addConfig('exact', default=[], ...
                validator=@(x) isa(x, 'function_handle') || isempty(x));
            obj.addConfig('components', default=struct('U', 1), ...
                validator=@(x) isstruct(x) && all(structfun(@(y) isempty(y) || isnumeric(y) && isscalar(y) && (y > 0), x)));
        end    

        function obj = configure(obj)
            % CONFIGURE Generate diffusion-specific configuration dependencies.
            
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
