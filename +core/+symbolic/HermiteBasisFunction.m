classdef HermiteBasisFunction < core.symbolic.InterpolationPolynomialBasisFunction
    % HERMITEBASISFUNCTION Hermite interpolation basis functions.
    %
    %   HermiteBasisFunction generates Hermite interpolation basis functions
    %   that interpolate both function values and derivatives at specified
    %   nodes. Each basis function satisfies specific conditions at the
    %   interpolation nodes for both the function and its derivatives.
    %
    %   For \f$n\f$ nodes, there are \f$2n\f$ Hermite basis functions:
    %   Index \f$1\f$ to \f$n\f$: function value interpolation
    %   Index \f$n+1\f$ to \f$2n\f$: derivative interpolation
    %
    % Examples:
    %   % Function value basis at first node
    %   nodes = [0, 1];
    %   basis = core.symbolic.HermiteBasisFunction(nodes, 1);
    %
    %   % Derivative basis in monic form
    %   basis = core.symbolic.HermiteBasisFunction(nodes, 3, sym('t'), true);
    %
    % See also:
    %   core.symbolic.InterpolationPolynomialBasisFunction
    
    methods (Access = protected)
        function expr = generate(obj)
            % GENERATE Generates the Hermite basis function.
            %
            %   expr = generate(obj) creates the symbolic expression for
            %   the Hermite interpolation basis function. Handles both
            %   function value and derivative interpolation conditions
            %   based on the index.
            %
            % Inputs:
            %   obj - The HermiteBasisFunction object
            %
            % Outputs:
            %   expr - Symbolic expression for the Hermite basis function
            
            x = obj.variables{1};
            n = length(obj.nodes);
            i = obj.index;
            nodes = obj.nodes;
            
            
            if i <= n
                H = obj.generateFunctionBasis(i, nodes, x);
            else
                nodeIdx = i - n;
                H = obj.generateDerivativeBasis(nodeIdx, nodes, x);
            end
            
            if obj.monic
                coeffs = sym2poly(H);
                leadingCoeff = coeffs(1);
                expr = horner(simplify(H / leadingCoeff), x);
            else
                expr = horner(simplify(H), x);
            end
        end
    end
    
    methods (Access = private)
        function H = generateFunctionBasis(~, nodeIdx, nodes, x)
            % GENERATEFUNCTIONBASIS Generates Hermite function value basis.
            
            n = length(nodes);
            xi = nodes(nodeIdx);
            
            L = sym(1);
            for j = 1:n
                if j ~= nodeIdx
                    L = L * (x - nodes(j)) / (xi - nodes(j));
                end
            end
            
            Lp_xi = subs(diff(L, x), x, xi);
            
            H = L * (1 - 2 * Lp_xi * (x - xi));
        end
        
        function H = generateDerivativeBasis(~, nodeIdx, nodes, x)
            % GENERATEDERIVATIVEBASIS Generates Hermite derivative basis.
            
            n = length(nodes);
            xi = nodes(nodeIdx);
            
            L = sym(1);
            for j = 1:n
                if j ~= nodeIdx
                    L = L * (x - nodes(j)) / (xi - nodes(j));
                end
            end
            
            H = (x - xi) * L^2;
        end
    end
end