classdef ConstantOperator < approx.assembly.Assembly
    % CONSTANTOPERATOR Assembly for constant operators.
    %
    %   ConstantOperator assembles diagonal matrices with constant 
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
    %   % Create constant operator
    %   constant = ConstantOperator(fe, mesh, op);
    %   
    %   % Assemble diagonal matrix with element-wise coefficients
    %   coeffs = [1.0, 2.0, 1.5]; % One coefficient per element
    %   matrix = constant.assemble(coeffs);
    %
    % See also:
    %   approx.assembly.Assembly

    methods
        function A = assemble(obj, coe)
            % ASSEMBLE Assemble diagonal matrix with constant coefficients.
            %
            %   A = assemble(obj, coe) creates a diagonal sparse matrix
            %   where each element contributes a constant coefficient
            %   repeated for all its local degrees of freedom.
            %
            % Inputs:
            %   obj - The ConstantOperator object
            %   coe - Coefficients
            %
            % Outputs:
            %   A - Diagonal sparse matrix

            ng = obj.space.nGlobalDofs;
            nl = obj.space.nLocalDofs;

            if isscalar(coe)
               coe = repmat(coe, obj.space.nTotalElements, 1); 
            end

            v = kron(coe(:), ones(nl, 1));
            A = spdiags(v, 0, ng, ng);
        end
    end
end