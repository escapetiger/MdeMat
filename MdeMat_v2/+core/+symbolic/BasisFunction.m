classdef BasisFunction < core.symbolic.SymbolicFunction
    % BASISFUNCTION Base class for all symbolic basis functions.
    %
    %   BasisFunction provides common functionality for basis function
    %   families such as Legendre, Chebyshev, Fourier, and interpolation
    %   polynomials. It handles index and multidimensional bounding box
    %   information for individual basis function construction.
    %
    % See Also:
    %   core.symbolic.SymbolicFunction,
    %   core.symbolic.PolynomialBasisFunction
    
    properties
        index % Index of the basis function (positive integer)
        bbox  % Multidimensional bounding box [a1, b1, a2, b2, ...]
    end
    
    methods
        function obj = BasisFunction(index, bbox, variables)
            % BASISFUNCTION Constructor for BasisFunction.
            %
            %   obj = BasisFunction(index, bbox) creates a new BasisFunction
            %   object with the specified index and bounding box using the
            %   default variable sym('x').
            %   
            %   obj = BasisFunction(index, bbox, variables) creates a new
            %   BasisFunction object with specified symbolic variables.
            %
            % Inputs:
            %   index - Positive integer index identifying the basis function
            %   bbox - Bounding box vector with even number of elements [a1, b1, a2, b2, ...]
            %   variables - Cell array of symbolic variables (optional, default:{sym('x')})
            %
            % Outputs:
            %   obj - Constructed BasisFunction object
            %
            % Examples:
            %   % 1D basis function on interval [0, 1]
            %   obj = MyBasisFunction(3, [0, 1]);
            %   
            %   % 2D basis function on domain [0,1] x [0,2] with custom variables
            %   vars = {sym('x'), sym('y')};
            %   obj = MyBasisFunction(5, [0, 1, 0, 2], vars);
            
            core.except.assert(index >= 1 && isscalar(index) && ...
                index == floor(index), 'InvalidInput', ...
                'Index must be a positive integer.');
            
            core.except.assert(isvector(bbox) && ~isempty(bbox), ...
                'InvalidInput', 'BBOX must be a non-empty vector.');
            
            core.except.assert(mod(length(bbox), 2) == 0, ...
                'InvalidInput', 'BBOX must have even number of elements.');
            
            for i = 1:2:length(bbox)
                core.except.assert(bbox(i) <= bbox(i+1), ...
                    'InvalidInput', ...
                    'Each interval [a, b] must satisfy a <= b.');
            end
            
            if nargin < 3
                variables = {sym('x')};
            end
            
            obj@core.symbolic.SymbolicFunction(variables);
            obj.index = index;
            obj.bbox = bbox;
        end
    end

    methods (Abstract, Access = protected)
        % GENERATE Abstract method to generate the symbolic expression.
        expr = generate(obj)
    end
end