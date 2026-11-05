classdef ElementBilinearForm < handle
    % ELEMENTBILINEARFORM Element-level bilinear form representation.
    %
    %   ElementBilinearForm represents a bilinear form through its matrix
    %   representation, computed from function data and derivatives. Each
    %   bilinear form corresponds to a weighted inner product between
    %   basis functions and their derivatives of the form:
    %
    %   \f[
    %     a(u, v) = \int_{\Omega} \frac{\partial^i u}{\partial x^i} 
    %     w(x) \frac{\partial^j v}{\partial x^j} \, dx
    %   \f]
    %
    %   where \f$u\f$ and \f$v\f$ are basis functions, \f$i\f$ and \f$j\f$
    %   are derivative orders, and \f$w(x)\f$ are integration weights.
    %   The discrete representation is computed using quadrature rules.
    %
    %   This class is fundamental in finite element methods for assembling
    %   system matrices from weak formulations of differential equations.
    %   It supports both matrix and vector representations for efficient
    %   integration computations.
    %
    % Examples:
    %   % Create empty bilinear form
    %   form = ElementBilinearForm();
    %   
    %   % Create mass matrix (0,0 derivatives)
    %   massForm = ElementBilinearForm(functionData1, functionData2, 0, 0);
    %   M = massForm.matrix;
    %
    %   % Create stiffness matrix (1,1 derivatives) 
    %   stiffForm = ElementBilinearForm(functionData1, functionData2, 1, 1);
    %   K = stiffForm.matrix;
    %
    %   % Mixed derivative form
    %   mixedForm = ElementBilinearForm(functionData1, functionData2, 0, 1);
    %   A = mixedForm.matrix;
    %   b = mixedForm.vector;
    %
    % See also:
    %   approx.element.ElementFunction, approx.integrate.Integrator
    
    properties
        matrix % Matrix representation of the bilinear form
        vector % Vector representation of the weighted test function values
    end
    
    methods
        function obj = ElementBilinearForm(D1, D2, i, j)
            % ELEMENTBILINEARFORM Constructor for ElementBilinearForm.
            %
            %   obj = ElementBilinearForm() creates an empty bilinear form
            %   data object with uninitialized matrix and vector properties.
            %
            %   obj = ElementBilinearForm(D1, D2, i, j) creates a bilinear
            %   form data object and initializes the matrix representation
            %   from function data using the specified derivative orders.
            %
            % Inputs:
            %   D1 - First function data object (ElementFunction)
            %   D2 - Second function data object (ElementFunction)
            %   i - Derivative order for first function (non-negative integer, 0 for values)
            %   j - Derivative order for second function (non-negative integer, 0 for values)
            %
            % Outputs:
            %   obj - Constructed ElementBilinearForm object
            
            if nargin == 0
                obj.matrix = []; 
                obj.vector = [];
            else
                obj.setData(D1, D2, i, j);
            end
        end

        function obj = setData(obj, D1, D2, i, j)
            % SETDATA Set the matrix representation of the bilinear form.
            %
            %   obj = setData(obj, D1, D2, i, j) computes the matrix
            %   representation as the weighted inner product between
            %   function values or derivatives from two function data
            %   objects. The computation follows the discrete bilinear form:
            %
            %   \f[
            %     A_{mn} = \sum_{k} U_k^{(i)}(x_m) w_k V_k^{(j)}(x_n)
            %   \f]
            %
            %   where \f$U^{(i)}\f$ and \f$V^{(j)}\f$ are the i-th and j-th
            %   derivatives, and \f$w_k\f$ are integration weights.
            %
            % Inputs:
            %   obj - The ElementBilinearForm object
            %   D1 - First function data object (ElementFunction)
            %   D2 - Second function data object (ElementFunction)
            %   i - Derivative order for first function (non-negative integer, 0 for values)
            %   j - Derivative order for second function (non-negative integer, 0 for values)
            %
            % Outputs:
            %   obj - The ElementBilinearForm object
            
            %< Create diagonal weight matrix
            W = diag(D1.weights);
            
            %< Extract function values or derivatives for first function
            if i == 0
                U = D1.values;
            else
                U = D1.derivatives{i};
            end
            
            %< Extract function values or derivatives for second function
            if j == 0
                V = D2.values;
            else
                V = D2.derivatives{j};
            end
            
            %< Compute bilinear form matrix: A = U^T * W * V
            A = U * W * V.';
            obj.matrix = A;
            
            %< Compute weighted test function vector: b = U^T * W
            obj.vector = U * W;
        end
    end
end