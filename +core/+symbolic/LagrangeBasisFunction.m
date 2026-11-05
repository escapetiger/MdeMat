classdef LagrangeBasisFunction < core.symbolic.InterpolationPolynomialBasisFunction
    % LAGRANGEBASISFUNCTION Lagrange interpolation basis functions.
    %
    %   LagrangeBasisFunction generates Lagrange interpolation basis
    %   functions associated with a specified set of nodes. Each basis
    %   function satisfies the Kronecker delta property at the
    %   interpolation nodes and can be in standard or monic form.
    %
    % Examples:
    %   % Standard Lagrange basis for first node
    %   nodes = [0, 1, 2];
    %   basis = core.symbolic.LagrangeBasisFunction(nodes, 1);
    %
    %   % Monic form for second node with custom variable
    %   basis = core.symbolic.LagrangeBasisFunction(nodes, 2, sym('t'), true);
    %
    % See also:
    %   core.symbolic.InterpolationPolynomialBasisFunction
    
    methods (Access = protected)
        function expr = generate(obj)
            % GENERATE Generates the Lagrange basis function.
            %
            %   expr = generate(obj) creates the symbolic expression for
            %   the Lagrange interpolation basis function using the product
            %   formula. The resulting polynomial equals 1 at the specified
            %   node and 0 at all other nodes.
            %
            % Inputs:
            %   obj - The LagrangeBasisFunction object
            %
            % Outputs:
            %   expr - Symbolic expression for the Lagrange basis function
            
            x = obj.variables{1};
            n = length(obj.nodes);
            i = obj.index;
            nodes = obj.nodes;
            
            L = sym(1);
            for j = 1:n
                if j ~= i
                    L = L * (x - nodes(j)) / (nodes(i) - nodes(j));
                end
            end
            
            if obj.monic
                coeffs = sym2poly(L);
                leadingCoeff = coeffs(1);
                expr = horner(simplify(L / leadingCoeff), x);
            else
                expr = horner(simplify(L), x);
            end
        end
    end
end