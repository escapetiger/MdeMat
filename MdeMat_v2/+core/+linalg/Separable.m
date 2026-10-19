classdef Separable < handle
    % SEPARABLE Base class for all separable objects.
    %
    %   Separable represents mathematical entity that can be
    %   decomposed as result of operation on multiple factors. Provides
    %   common functionality for all separable object types.
    %
    % Examples:
    %   % Create with uniform vector factors
    %   factors = {Vector([1,2,3]), Vector([4,5,6]), Vector([7,8,9])};
    %   obj = Separable(factors);
    %   factor2 = obj(2);  % Access using parentheses
    %
    %   % Create with mixed factor types
    %   factors = {Vector([1,2,3]), Matrix([4,5;6,7]), Scalar(8)};
    %   obj = Separable(factors);
    %   factor2 = obj{2}; % Access using curly braces

    properties
        factors % Array of factor objects
    end

    properties (Dependent)
        nFactors % Number of factors
    end

    methods
        function obj = Separable(factors)
            % SEPARABLE Constructor for separable object.
            %
            %   obj = Separable() creates empty Separable.
            %
            %   obj = Separable(factors) creates Separable with
            %   specified factors @a factors.
            %
            % Inputs:
            %   factors - Factor instances (cell array or object array)
            %
            % Outputs:
            %   obj - Separable instance

            if nargin < 1, factors = {}; end

            obj.setFactors(factors);
        end

        function n = get.nFactors(obj)
            % GET.NFACTORS Number of factors.

            n = length(obj.factors);
        end
    
        function obj = setFactors(obj, factors)
            % SETFACTORS Set the factors.
            %
            %   obj = setFactors(obj, factors) sets factors @a factors for
            %   the separable object.
            %
            % Inputs:
            %   obj - The Separable object
            %   factors - Factor instances (cell array or object array) 
            %
            % Outputs:
            %   obj = The Separable object

            if isempty(factors), return; end

            core.except.assert(iscell(factors) || isobject(factors), ...
                'InvalidInput', ...
                'Factors must be a cell array or object array.');

            if iscell(factors)
                if core.validate.isAllSameClass(factors{:})
                    obj.factors = arrayfun(@(i) factors{i}, 1:length(factors));
                else
                    obj.factors = factors;
                end
            else
                obj.factors = factors;
            end
        end
   
        function obj = setFactor(obj, index, factor)
            % SETFACTOR Set a specific factor.
            %
            %   obj = setFactor(obj, index, fatcor) sets the factor at the
            %   specified index.
            %
            % Inputs:
            %   obj - The Separable object
            %   index - Index of the factor (positive integer)
            %   factor - Factor instance
            %
            % Outputs:
            %   obj - The Separable object
            
            if iscell(obj.factors)
                obj.factors{index} = factor;
            else
                obj.factors(index) = factor;
            end
        end

        function F = getFactor(obj, index)
            % GETFACTOR Get a specific factor integrator.
            %
            %   F = getFactor(obj, index) returns the factor at the
            %   specified index.
            %
            % Inputs:
            %   obj - The Separable object
            %   index - Index of the factor (positive integer)
            %
            % Outputs:
            %   F - The factor at the specified index
            
            if iscell(obj.factors)
                F = obj.factors{index};
            else
                F = obj.factors(index);
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
            %
            % Inputs:
            %   factory - Function handle for creating factor objects
            %   varargin - Arguments passed to factory function
            %
            % Outputs:
            %   F - Cell array of created factor objects

            core.except.verify(isa(factory, 'function_handle'), ...
                'InvalidFunctionHandle', ...
                'First argument must be a function handle.');

            if isempty(varargin)
                F = {};
                return;
            end

            k = cellfun(@(v) length(v), varargin, 'Un', 0);
            core.except.verify(isequal(k{:}), ...
                'InvalidInput', 'All arguments must have the same size.');

            n = numel(varargin{1});
            F = cell(1, n);
            for i = 1:n
                args = cellfun(@(x) x(i), varargin, 'Un', 0);
                F{i} = factory(args{:});
            end
        end
    end
end