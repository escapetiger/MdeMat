classdef VlasovScheme < physics.scheme.MolScheme
    % VLASOVSCHEME Base class for all schemes solving Vlasov equations.
    %
    %   VlasovScheme implements a numerical scheme for solving Vlasov
    %   equations. This base class defines the basic validation and parse
    %   logic of configuration struct.
    %
    % See also:
    %   physics.scheme.MolScheme
    
    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options for Vlasov schemes.
            
            %< Call parent setup
            obj = setup@physics.scheme.MolScheme(obj);
            
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
                validator=@(x) isa(x, 'function_handle') || isempty(x));
            obj.addConfig('components', default=[], ...
                validator=@(x) isstruct(x));
        end    

        function obj = configure(obj)
            % CONFIGURE Generate Vlasov-specific configuration dependencies.
            
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
            if isstruct(obj.Config.exact)
                exactNames = fieldnames(obj.Config.exact);
                for i = 1:length(exactNames)
                    exactName = exactNames{i};
                    visualizer.addExact(exactName, obj.Config.exact.(exactName));
                end
            end

            %< Set analysis options
            analyzer = obj.Config.analyzer;
            analyzer.setReductions(obj.Config.errorTypes);
            if isstruct(obj.Config.components)
                analyzer.setComponents(obj.Config.components);
            end
            if isstruct(obj.Config.exact)
                exactNames = fieldnames(obj.Config.exact);
                for i = 1:length(exactNames)
                    exactName = exactNames{i};
                    analyzer.addExact(exactName, obj.Config.exact.(exactName));
                end
            end
        end
    end
end