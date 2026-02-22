classdef DiffusionScheme < physics.scheme.MolScheme
    % DIFFUSIONSCHEME Base class for all schemes solving diffusion
    % equations.
    %
    %   DiffusionScheme implements a numerical scheme for solving diffusion
    %   equations. This base class defines the basic validation and parse
    %   logic of configuration struct.
    %
    % See also:
    %   physics.scheme.MolScheme
    
    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options for diffusion schemes.
            
            %< Call parent setup
            obj = setup@physics.scheme.MolScheme(obj);
            
            %< Add diffusion-specific required configuration options
            obj.addConfig('diffusion', default=[], ...
                validator=@(x) ismatrix(x) && isnumeric(x));
            obj.addConfig('cfl', default=0.5, ...
                validator=@(x) isscalar(x) && x > 0);
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
            obj = configure@physics.scheme.MolScheme(obj);
            
            %< Determine boundary condition type
            if isempty(obj.Config.bc)
                obj.setConfig('bcType', 'periodic');
            else
                obj.setConfig('bcType', 'dirichlet');
            end

            %< Set visualization options
            visualizer = obj.Config.visualizer;
            visualizer.setComponents(obj.Config.components);
            visualizer.addExact('U', obj.Config.exact);

            %< Set analysis options
            analyzer = obj.Config.analyzer;
            analyzer.setComponents(obj.Config.components);
            analyzer.addExact('U', obj.Config.exact);
        end
    end
end