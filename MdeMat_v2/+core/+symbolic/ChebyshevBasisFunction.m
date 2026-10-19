classdef ChebyshevBasisFunction < core.symbolic.OrthogonalPolynomialBasisFunction
    % CHEBYSHEVBASISFUNCTION Chebyshev polynomial basis functions.
    %
    %   ChebyshevBasisFunction represents Chebyshev polynomials of the
    %   first kind on a specified interval. The polynomials can be in
    %   standard or monic form and are scaled to the given interval. These
    %   polynomials are orthogonal on [-1, 1] with respect to the weight
    %   function \f$w(x)=\frac{1}{\sqrt{1-x^2}}\f$. They also satifies the
    %   following three-term recurrence relation: 
    %   \f{eqnarray*}{
    %       &&T_0(t) = 1, \\
    %       &&T_1(t) = t, \\
    %       &&T_n(t) = 2tT_{n-1}(t) - T_{n-2}(t), \forall n \ge 2.
    %   \f}
    %
    % Examples:
    %   % Create Chebyshev polynomial of degree 2 on [0, 1]
    %   basis = core.symbolic.ChebyshevBasisFunction(2, [0, 1]);
    %
    %   % Create monic form on custom interval
    %   basis = core.symbolic.ChebyshevBasisFunction(3, [-1, 1], sym('t'), true);
    %
    % See also:
    %   core.symbolic.OrthogonalPolynomialBasisFunction

    methods (Access = protected)
        function expr = generate(obj)
            % GENERATE Generates the Chebyshev polynomial.
            %
            %   expr = generate(obj) creates the symbolic expression for
            %   the Chebyshev polynomial using the recurrence relation and
            %   scales it to the specified interval. Handles both standard
            %   and monic forms.
            %
            % Inputs:
            %   obj - The ChebyshevBasisFunction object
            %
            % Outputs:
            %   expr - Symbolic expression for the Chebyshev polynomial

            x = obj.variables{1};
            t = obj.scaleToCanonical(x, obj.bbox);
            T = obj.generateChebyshev(obj.degree, t);
            T = expand(T);
            if obj.degree == 0
                expr = sym(1);
            else
                if obj.monic
                    coeffs = sym2poly(T);
                    leadingCoeff = coeffs(1);
                    expr = horner(simplify(T/leadingCoeff), x);
                else
                    expr = horner(simplify(T), x);
                end
            end
        end
    end

    methods (Access = private)
        function T = generateChebyshev(~, n, t)
            % GENERATECHEBYSHEV Generates Chebyshev polynomial using
            % recurrence.

            if n == 0
                T = sym(1);
            elseif n == 1
                T = t;
            else
                T0 = sym(1);
                T1 = t;

                for k = 2:n
                    T = 2 * t * T1 - T0;
                    T0 = T1;
                    T1 = T;
                end
            end
        end
    end
end