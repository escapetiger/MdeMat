classdef ConstantAssembly < fem.assembly.Assembly
    % CONSTANTASSEMBLY Assembly for constant coefficient operators.
    %
    %   ConstantAssembly assembles diagonal matrices with constant 
    %   coefficients for each element. This is particularly useful for 
    %   mass matrices and other operators where each element contributes
    %   a uniform coefficient to all its local degrees of freedom.
    %
    %   The assembly creates block-diagonal structures where each
    %   element's contribution is a constant value multiplied by the
    %   identity matrix of size equal to the number of local degrees of
    %   freedom.
    %
    % Examples:
    %   % Create constant assembly
    %   assembly = ConstantAssembly(fe, mesh, op);
    %   
    %   % Assemble diagonal matrix with element-wise coefficients
    %   coeffs = [1.0, 2.0, 1.5]; % One coefficient per element
    %   matrix = assembly.assemble(coeffs);
    %
    %   % Mass matrix assembly
    %   density = ones(mesh.nTotalElements, 1);
    %   massMatrix = assembly.assemble(density);
    %
    % See also:
    %   fem.assembly.Assembly, fem.assembly.GridAssembly

    methods
        function A = assemble(obj, c)
            % ASSEMBLE Assemble diagonal matrix with constant coefficients.
            %
            %   A = assemble(obj, c) creates a diagonal sparse matrix where
            %   each element contributes a constant coefficient repeated for
            %   all its local degrees of freedom. The resulting matrix is
            %   block-diagonal with each block being a scaled identity.
            %
            % Inputs:
            %   obj - The ConstantAssembly object
            %   c - Vector of coefficients, one per element (length nTotalElements)
            %
            % Outputs:
            %   A - Assembled diagonal sparse matrix of size nGlobalDofs×nGlobalDofs

            n = obj.nGlobalDofs;
            np = obj.nLocalDofs;
            v = kron(c(:), ones(np, 1));
            A = spdiags(v, 0, n, n);
        end
    end
end