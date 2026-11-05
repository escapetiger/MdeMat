classdef OrthogonalPolynomialBasisFunction < core.symbolic.PolynomialBasisFunction
    % ORTHOGONALPOLYNOMIALBASISFUNCTION Base class for all orthogonal
    % polynomial basis functions.
    %
    %   OrthogonalPolynomialBasisFunction provides common functionality for
    %   orthogonal polynomial families such as Legendre, Chebyshev, and
    %   Jacobi polynomials. It automatically sets the index as degree + 1
    %   following standard mathematical convention and provides utilities
    %   for orthogonal polynomial construction.
    %
    %   Orthogonal polynomials satisfy orthogonality conditions: 
    %   \f[
    %     \int_{a}^{b} P_m(x) P_n(x) w(x) dx = 0, \forall m \neq n.
    %   \f]
    %   where \f$w(x)\f$ is the weight function.
    %
    % See also:
    %   core.symbolic.PolynomialBasisFunction
    
    methods
        function obj = OrthogonalPolynomialBasisFunction(degree, bbox, variable, monic)
            % ORTHOGONALPOLYNOMIALBASISFUNCTION Constructor for
            % OrthogonalPolynomialBasisFunction.
            %
            %   Creates a new orthogonal polynomial basis function with the
            %   specified degree and properties. Automatically sets the
            %   index as degree + 1 following standard convention.
            %
            % Inputs:
            %   degree - Polynomial degree (non-negative integer)
            %   bbox - Bounding box [lower, upper]
            %   variable - Symbol(s) (optional, default: {sym('x')})
            %   monic - Monic form flag (optional, default: false)
            %
            % Outputs:
            %   obj - Constructed OrthogonalPolynomialBasisFunction object
            
            if nargin < 3
                variable = sym('x');
            end
            if nargin < 4
                monic = false;
            end
            
            obj@core.symbolic.PolynomialBasisFunction(degree, ...
                degree + 1, bbox, variable, monic);
        end
    end
end