classdef LegendreBasisFunction < core.symbolic.OrthogonalPolynomialBasisFunction
    % LEGENDREBASISFUNCTION Legendre polynomial basis functions.
    %
    %   LegendreBasisFunction represents Legendre polynomials on a
    %   specified interval. The polynomials can be in standard or monic
    %   form and are scaled to the given interval. The underlying Legendre
    %   polynomials are orthogonal on \f$[-1,1]\f$ with respect to the
    %   weight function \f$w(x)=1\f$.
    %
    % Examples:
    %   % Create Legendre polynomial of degree 2 on [0, 1]
    %   basis = core.symbolic.LegendreBasisFunction(2, [0, 1]);
    %
    %   % Create monic form on custom interval
    %   basis = core.symbolic.LegendreBasisFunction(3, [-1, 1], sym('t'), true);
    %
    % See also:
    %   core.symbolic.OrthogonalPolynomialBasisFunction
        
    methods (Access = protected)
        function expr = generate(obj)
            % GENERATE Generates the Legendre polynomial.
            %
            %   expr = generate(obj) creates the symbolic expression for
            %   the Legendre polynomial using MATLAB's legendreP function
            %   and scales it to the specified interval. Handles both
            %   standard and monic forms.
            %
            % Inputs:
            %   obj - The LegendreBasisFunction object
            %
            % Outputs:
            %   expr - Symbolic expression for the Legendre polynomial
            
            x = obj.variables{1};
            t = obj.scaleToCanonical(x, obj.bbox);
            
            P = legendreP(obj.degree, t);
            P = expand(P);
            
            if obj.monic
                coeffs = sym2poly(P);
                leadingCoeff = coeffs(1);
                expr = horner(simplify(P / leadingCoeff), x);
            else
                expr = horner(simplify(P), x);
            end
        end
    end
end