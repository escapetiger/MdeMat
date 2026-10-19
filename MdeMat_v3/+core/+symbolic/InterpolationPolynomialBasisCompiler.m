classdef InterpolationPolynomialBasisCompiler < core.symbolic.Compiler
    % INTERPOLATIONPOLYNOMIALBASISCOMPILER Compiler for interpolation polynomial basis functions.
    %
    %   InterpolationPolynomialBasisCompiler specializes in compiling basis
    %   matrices for interpolation polynomial families such as Lagrange and
    %   Hermite polynomials. It creates and manages arrays of interpolation
    %   polynomial basis functions based on specified node sets.
    %
    %   Supported interpolation types:
    %   - Lagrange: n basis functions for n nodes (function value interpolation)
    %   - Hermite: 2n basis functions for n nodes (function + derivative interpolation)
    %
    % Examples:
    %   % Create compiler and set Lagrange functions
    %   compiler = core.symbolic.InterpolationPolynomialBasisCompiler();
    %   compiler.setFunctions([0, 1, 2], 'lagrange');
    %   
    %   % Compile and load basis matrix
    %   compiler.compile('lagrange_basis', 'Lagrange', 2);
    %   metadata = compiler.load('lagrange_basis');
    %   
    %   % Hermite interpolation
    %   compiler.setFunctions([-1, 0, 1], 'hermite');
    %   compiler.compile('hermite_basis', 'Hermite', 1);
    %
    % See Also:
    %   core.symbolic.Compiler, core.symbolic.LagrangeBasisFunction, 
    %   core.symbolic.HermiteBasisFunction
    
    methods
        function obj = setFunctions(obj, nodes, basisType)
            % SETFUNCTIONS Sets interpolation polynomial basis functions.
            %
            %   obj = setFunctions(obj, nodes, basisType) creates an array
            %   of interpolation polynomial basis functions based on the
            %   specified nodes and basis type. Dispatches to the
            %   appropriate specialized method for the chosen interpolation
            %   type.
            %
            % Inputs:
            %   obj - The InterpolationPolynomialBasisCompiler object
            %   nodes - Array of interpolation nodes (vector of distinct real numbers)
            %   basisType - Type of basis ('lagrange', 'hermite')
            %
            % Outputs:
            %   obj - The InterpolationPolynomialBasisCompiler object
            
            core.except.assert(isvector(nodes) && ~isempty(nodes), ...
                'InvalidInput', 'Nodes must be a non-empty vector.');
            
            core.except.assert(length(unique(nodes)) == length(nodes), ...
                'InvalidInput', 'All nodes must be distinct.');
            
            supportedTypes = {'lagrange', 'hermite'};
            core.except.assert(ismember(lower(basisType), supportedTypes), ...
                'InvalidInput', 'Supported basis types: %s', strjoin(supportedTypes, ', '));
            
            switch lower(basisType)
                case 'lagrange'
                    obj.setLagrange(nodes, basisType);
                case 'hermite'
                    obj.setHermite(nodes, basisType);
                otherwise
                    core.except.assert(0, 'UnsupportedType', ...
                        'Basis type "%s" not implemented', basisType);
            end
        end
        
        function setLagrange(obj, nodes, basisType)
            % SETLAGRANGE Sets Lagrange polynomial basis functions.
            %
            %   setLagrange(obj, nodes, basisType) creates an array of
            %   Lagrange interpolation basis functions for the specified
            %   nodes. Each function satisfies the Kronecker delta property
            %   at the interpolation nodes: L_i(x_j) = δ_{ij}.
            %
            % Inputs:
            %   obj - The InterpolationPolynomialBasisCompiler object
            %   nodes - Array of interpolation nodes
            %   basisType - Type ('lagrange')
            %
            % Outputs:
            %   NULL
            
            x = obj.variables;
            n = length(nodes);
            
            obj.funcs = arrayfun( ...
                @(i) core.symbolic.LagrangeBasisFunction(nodes, i, x), 1:n);
            
            fprintf('Set %d %s basis functions for %d nodes\n', n, basisType, n);
        end
        
        function setHermite(obj, nodes, basisType)
            % SETHERMITE Sets Hermite polynomial basis functions.
            %
            %   setHermite(obj, nodes, basisType) creates an array of
            %   Hermite interpolation basis functions for the specified
            %   nodes. Includes both function value and derivative
            %   interpolation functions for each node.
            %
            % Inputs:
            %   obj - The InterpolationPolynomialBasisCompiler object
            %   nodes - Array of interpolation nodes
            %   basisType - Type ('hermite')
            %
            % Outputs:
            %   NULL
            
            x = obj.variables;
            n = length(nodes);
            
            obj.funcs = arrayfun( ...
                @(i) core.symbolic.HermiteBasisFunction(nodes, i, x), 1:(2*n));
           
            fprintf('Set %d %s basis functions for %d nodes (%d functions total)\n', ...
                2*n, basisType, n, 2*n);
        end
    end
end