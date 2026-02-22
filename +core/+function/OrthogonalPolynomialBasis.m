classdef (Abstract) OrthogonalPolynomialBasis < core.function.Function
    % ORTHOGONALPOLYNOMIALBASIS Abstract base class for orthogonal
    % polynomial basis functions.
    %
    %   OrthogonalPolynomialBasis provides a unified framework for
    %   generating orthogonal polynomial basis functions using three-term
    %   recurrence relations. Concrete subclasses implement specific
    %   polynomial families such as Legendre, Chebyshev, and Hermite
    %   polynomials.
    %
    %   The three-term recurrence relation has the general form:
    %
    %   \f[
    %       P_{n+1}(x) = (A_n \cdot x + B_n) \cdot P_n(x) + C_n \cdot
    %       P_{n-1}(x)
    %   \f]
    %
    % See also:
    %   core.function.LegendreBasis, core.function.ChebyshevBasis,
    %   core.function.HermiteBasis
    
    properties
        MaxDegree % Maximum polynomial degree
        IsNormalized % Flag for orthonormal basis functions
    end
    
    methods
        function obj = OrthogonalPolynomialBasis(maxDegree, options)
            % ORTHOGONALPOLYNOMIALBASIS Construct orthogonal polynomial
            % basis functions.
            %
            %   obj = OrthogonalPolynomialBasis(maxDegree) creates
            %   orthogonal polynomial basis functions using the specified
            %   @a variables up to @a maxDegree.
            
            arguments
                maxDegree(1, 1) {mustBeNonnegative, mustBeInteger}
                options.isNormalized(1, 1) {mustBeNumericOrLogical} = false
            end
            
            obj = obj@core.function.Function(nDims = 1, nCodims = maxDegree+1);
            obj.MaxDegree = maxDegree;
            obj.IsNormalized = options.isNormalized;
        end
    end
    
    methods (Access = protected)
        function Y = evalImpl(obj, X)
            % EVALIMPL Implementation of function evaluation.
            %
            %   Y = evalImpl(obj, X) evaluates the basis values via
            %   three-term recurrence formula.
            
            np = size(X, 2);
            Y = zeros(obj.NCodims, np);

            %< First term
            Y(1, :) = obj.first(X);

            if obj.NCodims == 1, return; end
            
            %< Second term
            [a, b, ~] = obj.computeCoefficients(0);
            Y(2, :) = (a * X + b) .* Y(1, :);

            if obj.NCodims == 2, return; end

            %< Recurrence formula
            for k = 2:obj.MaxDegree
                [a, b, c] = obj.computeCoefficients(k-1);
                Y(k+1, :) = (a * X + b) .* Y(k, :) + c * Y(k-1, :);
            end
        end
    end
    
    methods (Abstract, Access = protected)
        % FIRST Compute first term.
        Y = first(obj, X)

        % COMPUTECOEFFICIENTS Compute coefficients for three-term
        % recurrence relation.
        %
        %   [a, b, c] = computeCoefficients(n) returns the coefficients
        %   @a a, @a b, and @a c for the three-term recurrence relation
        %   at degree @a n. The recurrence relation is: \f$ P_{n+1}(x) =
        %   (a*x + b)*P_n(x) + c*P_{n-1}(x) \f$
        [a, b, c] = computeCoefficients(obj, n)
    end
end