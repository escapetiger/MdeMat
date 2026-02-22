classdef RadiationScheme < physics.scheme.MolScheme
    % RADIATIONSCHEME Base class for all schemes solving radiation
    % transport equations.
    %
    %   RadiationScheme implements a numerical scheme for solving radiative
    %   transfer equation. This base class defines the basic validation and
    %   parse logic of configuration struct.
    %
    % See also:
    %   physics.scheme.MolScheme

    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options for radiation
            % schemes.

            %< Call parent setup
            obj = setup@physics.scheme.MolScheme(obj);

            %< Initial and boundary conditions
            obj.addConfig('ic', default=struct(), ...
                validator=@(x) isstruct(x) || isa(x, 'function_handle'));
            obj.addConfig('bc', default=[], ...
                validator=@(x) true);

            %< Time stepping
            obj.addConfig('cfl', default=0.5, ...
                validator=@(x) isempty(x) || (isscalar(x) && x > 0));
            obj.addConfig('dt', default=[], ...
                validator=@(x) isempty(x) || (isscalar(x) && x > 0));

            %< Scaling parameter
            obj.addConfig('epsilon', default=1, ...
                validator=@(x) isnumeric(x) && isscalar(x) && x > 0);

            %< Physical coefficients
            obj.addConfig('scattering', default=[], ...
                validator=@(x) isnumeric(x) || isa(x, 'function_handle') || isempty(x));
            obj.addConfig('absorption', default=[], ...
                validator=@(x) isnumeric(x) || isa(x, 'function_handle') || isempty(x));
            obj.addConfig('source', default=[], ...
                validator=@(x) isnumeric(x) || isa(x, 'function_handle') || isempty(x));

            %< Discretization options
            obj.addConfig('vDimReduction', default='symmetry', ...
                validator=@(x) ismember(x, {'', 'symmetry', 'topology'}));

            %< Analysis and visualization
            obj.addConfig('eId', default='Radiation', ...
                validator=@(x) (ischar(x) || isstring(x)));
            obj.addConfig('errorTypes', default={'L1', 'L2', 'max'}, ...
                validator=@(x) iscell(x) && all(cellfun(@(y) ischar(y) || isstring(y), x)));
            obj.addConfig('exact', default=[], ...
                validator=@(x) isstruct(x) || isa(x, 'function_handle') || isempty(x));
            obj.addConfig('components', default=struct('U', [], 'G', [], 'F', 1), ...
                validator=@(x) isstruct(x) && all(structfun(@(y) isempty(y) || isnumeric(y) && isscalar(y) && (y > 0), x)));
        end

        function obj = configure(obj)
            % CONFIGURE Generate radiation-specific configuration
            % dependencies.

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
            if isstruct(obj.Config.components)
                visualizer.setComponents(obj.Config.components);
            end
            if isa(obj.Config.exact, 'function_handle')
                visualizer.addExact('F', obj.Config.exact);
            elseif isstruct(obj.Config.exact)
                exactNames = fieldnames(obj.Config.exact);
                for i = 1:length(exactNames)
                    exactName = exactNames{i};
                    visualizer.addExact(exactName, obj.Config.exact.(exactName));
                end
            end

            %< Set analysis options
            analyzer = obj.Config.analyzer;
            analyzer.setComponents(obj.Config.components);
            if isa(obj.Config.exact, 'function_handle')
                analyzer.addExact('F', obj.Config.exact);
            elseif isstruct(obj.Config.exact)
                exactNames = fieldnames(obj.Config.exact);
                for i = 1:length(exactNames)
                    exactName = exactNames{i};
                    analyzer.addExact(exactName, obj.Config.exact.(exactName));
                end
            end
        end
    end
end
