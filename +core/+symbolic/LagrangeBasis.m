classdef LagrangeBasis < core.symbolic.SymbolicFunction
    % LAGRANGEBASIS Lagrange interpolation basis functions.
    %
    %   LagrangeBasis generates Lagrange interpolation basis functions for
    %   given interpolation nodes. Each basis function is a polynomial that
    %   equals 1 at one node and 0 at all other nodes.
    %
    % See also:
    %   core.symbolic.SymbolicFunction, core.symbolic.LegendreBasis

    properties
        Nodes % Interpolation nodes
    end

    methods
        function obj = LagrangeBasis(variables, nodes)
            % LAGRANGEBASIS Construct Lagrange basis functions.
            %
            %   obj = LagrangeBasis(variables, nodes) creates Lagrange
            %   basis functions using the specified @a variables and
            %   interpolation @a nodes.

            arguments   
                variables {mustBeA(variables, 'sym')}
                nodes {mustBeA(nodes, 'sym'), mustBeNonempty}
            end

            obj = obj@core.symbolic.SymbolicFunction(variables);
            obj.Nodes = nodes;
            obj.computeExpressions();
        end
    end

    methods (Access = private)
        function computeExpressions(obj)
            % COMPUTEEXPRESSIONS Compute Lagrange basis function expressions.

            z = obj.Nodes;
            n = length(obj.Nodes);
            obj.Expressions = repmat(sym([]), n, 1);

            x = obj.Variables(1);
            for i = 1:n
                p = sym(1);
                for j = 1:n
                    if i ~= j
                        p = p * (x - z(j)) / (z(i) - z(j));
                    end
                end
                obj.Expressions(i) = horner(simplify(p));
            end
        end
    end
end