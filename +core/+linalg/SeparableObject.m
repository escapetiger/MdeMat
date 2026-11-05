classdef SeparableObject < handle
    % SEPARABLEOBJECT Base class for all separable objects.
    %
    %   SeparableObject represents mathematical entity that can be
    %   decomposed as result of operation on multiple factors. Provides
    %   common functionality for all separable object types.
    %
    % Examples:
    %   % Create with uniform vector factors
    %   factors = {Vector([1,2,3]), Vector([4,5,6]), Vector([7,8,9])};
    %   obj = SeparableObject(factors);
    %   factor2 = obj(2);  % Access using parentheses
    %
    %   % Create with mixed factor types
    %   factors = {Vector([1,2,3]), Matrix([4,5;6,7]), Scalar(8)};
    %   obj = SeparableObject(factors);
    %   factor2 = obj{2}; % Access using curly braces
    %
    % See Also:
    %   core.utilities.mustBeCellOrObject

    properties
        factors % Array of factor objects
    end

    methods
        function obj = SeparableObject(F)
            % SEPARABLEOBJECT Constructor for separable object.
            %
            %   obj = SeparableObject() creates empty SeparableObject.
            %
            %   obj = SeparableObject(F) creates SeparableObject with
            %   specified factors @a F.
            %
            % Inputs:
            %   F - Factor instances (cell array or object array)
            %
            % Outputs:
            %   obj - SeparableObject instance

            if nargin < 1, F = {}; end

            obj.setFactors(F);
        end
    
        function setFactors(obj, F)
            % SETFACTORS Set the factors.
            %
            %   setFactors(obj, F) sets factors @a F for the object.
            %
            % Inputs:
            %   obj - The SeparableObject
            %   F - Factor instances (cell array or object array) 
            %
            % Outputs:
            %   NULL

            if isempty(F), return; end

            core.except.assert(iscell(F) || isobject(F), ...
                'InvalidInput', ...
                'Factors must be a cell array or object array.');

            if iscell(F)
                if core.validate.isAllSameClass(F{:})
                    try
                        n = numel(F);
                        G(n) = F{1};
                        for i = 1:n
                            G(i) = F{i};
                        end
                    catch
                        G = F;
                    end
                    obj.factors = G;
                else
                    obj.factors = F;
                end
            else
                obj.factors = F;
            end
        end
    end

    methods (Static)
        function F = createFactors(factory, varargin)
            % CREATEFACTORS Create cell array of factors using factory
            % function.
            %
            %   F = SeparableObject.createFactors(factory, arg1, arg2, ...)
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