classdef Separable < handle
    % SEPARABLE Base class for all separable objects.
    %
    %   Separable represents mathematical entity that can be
    %   decomposed as result of operation on multiple factors. Provides
    %   common functionality for all separable object types.

    properties
        Factors % Array of factor objects
    end

    properties (Dependent)
        NFactors % Number of factors
    end

    methods
        function obj = Separable(options)
            % SEPARABLE Construct an instance of Separable.
            %
            %   obj = Separable() creates empty Separable.
            %
            %   obj = Separable(Name=Value) creates Separable with
            %   specified options.

            arguments
                options.factors {mdemat.tools.base.mustBeCellOrObject} = {}
            end

            obj.setFactors(options.factors);
        end

        function n = get.NFactors(obj)
            % GET.NFACTORS Returns the value of the dependent property 
            % 'NFactors'.

            n = numel(obj.Factors);
        end
    
        function obj = setFactors(obj, factors)
            % SETFACTORS Set the factors.
            %
            %   obj = setFactors(obj, factors) sets factors @a factors for
            %   the separable object.
            %
            %   Factors must be provided as a cell array or object array.
            %   Empty factors are ignored. Array factors are visited in
            %   column-major order; every element is retained.

            arguments
                obj mdemat.data.variable.Separable
                factors {mdemat.tools.base.mustBeCellOrObject}
            end

            if isempty(factors), return; end

            if iscell(factors)
                if mdemat.tools.base.isAllSameClass(factors{:})
                    obj.Factors = arrayfun(@(i) factors{i}, 1:numel(factors));
                else
                    obj.Factors = factors;
                end
            else
                obj.Factors = factors;
            end
        end
   
        function obj = setFactor(obj, index, factor)
            % SETFACTOR Set a specific factor.
            %
            %   obj = setFactor(obj, index, factor) sets the factor at the
            %   specified index.

            arguments
                obj mdemat.data.variable.Separable
                index {mustBeInteger, mustBePositive}
                factor {mustBeA(factor, 'object')}
            end
            
            if iscell(obj.Factors)
                obj.Factors{index} = factor;
            else
                obj.Factors(index) = factor;
            end
        end

        function F = getFactor(obj, index)
            % GETFACTOR Get a specific factor integrator.
            %
            %   F = getFactor(obj, index) returns the factor at the
            %   specified index.

            arguments
                obj mdemat.data.variable.Separable
                index {mustBeInteger, mustBePositive}
            end
            
            if iscell(obj.Factors)
                F = obj.Factors{index};
            else
                F = obj.Factors(index);
            end
        end
    end

    methods (Static)
        function F = createFactors(factory, varargin)
            % CREATEFACTORS Create cell array of factors using factory
            % function.
            %
            %   F = Separable.createFactors(factory, arg1, arg2, ...)
            %   creates multiple factor objects using factory function @a
            %   factory.

            assert(all(isa(factory, 'function_handle'), 'all'), ...
                'mdemat:data:variable:Separable:InvalidFunctionHandle', ...
                '%s', 'First argument must be a function handle.');

            if isempty(varargin)
                F = {};
                return;
            end

            counts = cellfun(@numel, varargin);
            assert(all(all(counts == counts(1)), 'all'), ...
                'mdemat:data:variable:Separable:InvalidInput', '%s', 'All arguments must have the same size.');

            n = numel(varargin{1});
            F = cell(1, n);
            for i = 1:n
                args = cellfun(@(x) x(i), varargin, 'Un', 0);
                F{i} = factory(args{:});
            end
        end
    end
end
