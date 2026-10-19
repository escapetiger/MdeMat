classdef PolynomialBasisFunction < core.symbolic.BasisFunction
    % POLYNOMIALBASISFUNCTION Abstract base class for univariate polynomial
    % basis functions.
    %
    %   PolynomialBasisFunction provides common functionality for
    %   univariate polynomial basis function families such as Legendre,
    %   Chebyshev, Jacobi, and interpolation polynomials. It handles degree
    %   information and monic form specification for polynomial basis
    %   construction.
    %
    %   The class supports both standard and monic polynomial forms:
    %   - Standard form: preserves the natural normalization
    %   - Monic form: normalizes so the leading coefficient equals 1
    %
    % See also:
    %   core.symbolic.BasisFunction
    
    properties (Access = public)
        degree % Degree of the polynomial basis function
        monic  % Flag indicating whether to use monic polynomial form
    end
    
    methods
        function obj = PolynomialBasisFunction(degree, index, bbox, variable, monic)
            % POLYNOMIALBASISFUNCTION Constructor for
            % PolynomialBasisFunction.
            %
            %   obj = PolynomialBasisFunction(degree, index, bbox) creates
            %   a new polynomial basis function with specified degree,
            %   index, and bounding box using default settings.
            %   
            %   obj = PolynomialBasisFunction(degree, index, bbox,
            %   variable, monic) creates a new polynomial basis function
            %   with full customization.
            %
            % Inputs:
            %   degree - Polynomial degree (non-negative integer)
            %   index - Index of the basis function (positive integer)
            %   bbox - Bounding box [lower, upper]
            %   variable - Symbol(s) (optional, default: {sym('x')})
            %   monic - Monic form flag (optional, default: false)
            %
            % Outputs:
            %   obj - Constructed PolynomialBasisFunction object
            
            core.except.assert(degree >= 0, 'InvalidInput', ...
                'Degree must be a non-negative integer.');
            
            if nargin < 4
                variable = sym('x');
            end
            if nargin < 5
                monic = false;
            end
            
            obj@core.symbolic.BasisFunction(index, bbox, variable);
            obj.degree = degree;
            obj.monic = monic;
            obj.expression = obj.generate();
        end
    end
    
    methods (Static, Access = protected)
        function result = scaleToCanonical(variable, bbox)
            % SCALETOCANONICAL Transforms variable from [a, b] to canonical
            % interval [-1, 1].
            %
            %   Performs linear transformation to map the variable from
            %   the specified bounding box to the standard interval [-1, 1]
            %   used by many orthogonal polynomial families.
            %
            % Inputs:
            %   variable - Symbolic variable on interval [a, b]
            %   bbox - Bounding box [a, b]
            %
            % Outputs:
            %   result - Variable scaled to canonical interval [-1, 1]
            
            a = bbox(1);
            b = bbox(2);
            result = (2 * variable - (a + b)) / (b - a);
        end
    end
end