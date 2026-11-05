classdef LinearApproximator < handle
    % LINEARAPPROXIMATOR Base class for all linear approximators.
    %
    %   LinearApproximator provides the foundation for linear approximation
    %   methods where functions are approximated as finite linear
    %   combinations of basis functions. The class supports projection and
    %   fitting using mass and design operators to solve linear systems. In
    %   addition, the class implements the evaluation logic of linear
    %   approximator.
    %
    %   Let \f$V_h\f$ be a linear approximation space of which the element
    %   \f$u_h\f$ is of form:
    %
    %   \f[
    %       u_h = \sum_{i=1}^n c_i p_i.
    %   \f]
    %
    %   The projection seeks a linear approximation \f$u_h \in V_h\f$ of a
    %   function \f$u\f$ such that
    %
    %   \f[
    %       \langle u_h, v_h \rangle = f(v_h), \quad \forall v_h \in V_h,
    %   \f]
    %
    %   The fitting seeks a linear approximation \f$u_h \in V_h\f$ of a
    %   function \f$u\f$ such that
    %
    %   \f[
    %       \phi_j(u) = \phi_j(u_h), \quad \forall 1 \le j \le m,
    %   \f]
    %
    %   where \f$\phi_j\f$ are linear functionals.

    properties (Access = public)
        basis  % Basis functions that span the approximation space
        mass   % Mass matrix for projection operations (square matrix)
        design % Design matrix for fitting operations (matrix)
    end

    properties (Dependent)
        nDofs % Number of degrees of freedom (scalar)
    end

    methods
        function obj = LinearApproximator(basis)
            % LINEARAPPROXIMATOR Constructor for LinearApproximator.
            %
            %   obj = LinearApproximator(basis) creates a new
            %   LinearApproximator object with the specified basis
            %   functions. The mass and design operators are initialized as
            %   empty and must be set by calling setMass() and setDesign().
            %
            % Inputs:
            %   basis - Basis functions (object)
            %
            % Outputs:
            %   obj - The constructed LinearApproximator object
            
            core.except.assert(~isempty(basis), 'InvalidInput', ...
                'Basis functions cannot be empty.');
            
            obj.basis = basis;
            obj.mass = [];
            obj.design = [];
        end        

        function n = get.nDofs(obj)
            % GET.NDOFS Get the number of degrees of freedom.
            %
            %   n = get.nDofs(obj) returns the number of degrees of freedom
            %   in the linear approximation space, which equals the number of
            %   basis function codimensions.
            %
            % Inputs:
            %   obj - The LinearApproximator object
            %
            % Outputs:
            %   n - Number of degrees of freedom (positive integer)

            n = obj.basis.nCodims;
        end
        
        function obj = setMass(obj, X, w)
            % SETMASS Set up the mass operator using numerical integration.
            %
            %   obj = setMass(obj, X, w) constructs the mass matrix using
            %   integration points and weights. The (i,j) entry of the mass
            %   matrix is \f$\langle p_i, p_j\rangle\f$ where the inner
            %   product is defined by the weighted integration: 
            %   \f$\langle f,g\rangle = \sum_k w_k f(x_k) g(x_k)\f$.
            %
            % Inputs:
            %   obj - The LinearApproximator object
            %   X - Integration points (n x q matrix)
            %   w - Integration weights (1 x q vector or q x 1 vector)
            %
            % Outputs:
            %   obj - The updated LinearApproximator object
            
            n = obj.basis.nDims;
            X = reshape(X, size(X, 1), []);
            
            core.except.assert(size(X, 1) >= n, 'InvalidInput', ...
                'Integration points must have at least %d rows.', n);
            
            X = X(1:n, :);
            w = w(:);
            
            q = length(w);
            core.except.assert(size(X, 2) == q, 'InvalidInput', ...
                'Number of integration points (%d) must match number of weights (%d).', ...
                size(X, 2), q);
            
            core.except.assert(all(w >= 0), 'InvalidInput', ...
                'Integration weights must be non-negative.');
            
            B = obj.basis.evaluate(X);
            W = diag(w);
            M = B * W * B.';
            M = core.linalg.zeroing(M);
            
            obj.mass = M;
        end
        
        function Y = evaluate(obj, varargin)
            % EVALUATE Evaluate linear approximations at given points.
            %
            %   Y = evaluate(obj, Z) evaluates functions represented by
            %   coefficients embedded in the augmented vector Z = [X; C]
            %   at evaluation points X.
            %
            %   Y = evaluate(obj, X, C) evaluates functions with
            %   coefficient matrix C at evaluation points X.
            %
            % Inputs:
            %   obj - The LinearApproximator object
            %   varargin - Variable input arguments:
            %<     Z - Augmented matrix [X; C]
            %<     X - Evaluation points (n x p matrix)
            %<     C - Coefficient matrix (m x q matrix)
            %
            % Outputs:
            %   Y - Function values (1 x p vector or q x p matrix)
            
            switch length(varargin)
                case 1
                    n = obj.basis.nDims;
                    Z = varargin{1};
                    Z = reshape(Z, size(Z, 1), []);
                    
                    core.except.assert(size(Z, 1) > n, 'InvalidInput', ...
                        'Input Z must have more rows than basis dimensions.');
                    
                    C = Z((n+1):end, :);
                    X = Z(1:n, :);
                    B = obj.basis.evaluate(X);
                    Y = dot(B, C, 1);
                    
                case 2
                    [X, C] = varargin{:};
                    X = reshape(X, size(X, 1), []);
                    C = reshape(C, size(C, 1), []);
                    B = obj.basis.evaluate(X);
                    Y = C.' * B;
                    
                otherwise
                    core.except.assert(0, 'InvalidInput', ...
                        'Invalid number of input arguments.');
            end
        end
        
        function Y = differentiate(obj, varargin)
            % DIFFERENTIATE Differentiate linear approximations at given
            % points.
            %
            %   Y = differentiate(obj, Z, order) differentiates functions
            %   represented by coefficients embedded in the augmented
            %   vector Z = [X; C] at evaluation points X.
            %
            %   Y = differentiate(obj, X, C, order)  differentiates
            %   functions with coefficient matrix C at evaluation points X.
            %
            % Inputs:
            %   obj - The LinearApproximator object
            %   varargin - Variable input arguments:
            %<     Z - Augmented matrix [X; C]
            %<     X - Evaluation points (n x p matrix)
            %<     C - Coefficient matrix (m x q matrix)
            %<     order - Derivative order (n x 1 vector)
            %
            % Outputs:
            %   Y - Function values (1 x p vector or q x p matrix)
            
            switch length(varargin)
                case 2
                    n = obj.basis.nDims;
                    [Z, order] = varargin{:};
                    Z = reshape(Z, size(Z, 1), []);
                    
                    core.except.assert(size(Z, 1) > n, 'InvalidInput', ...
                        'Input Z must have more rows than basis dimensions.');
                    
                    C = Z((n+1):end, :);
                    X = Z(1:n, :);
                    B = obj.basis.differentiate(X, order);
                    Y = dot(B, C, 1);
                    
                case 3
                    [X, C, order] = varargin{:};
                    X = reshape(X, size(X, 1), []);
                    C = reshape(C, size(C, 1), []);
                    B = obj.basis.differentiate(X, order);
                    Y = C.' * B;
                    
                otherwise
                    core.except.assert(0, 'InvalidInput', ...
                        'Invalid number of input arguments.');
            end
        end


        function U = fit(obj, F)
            % FIT Fit data within the linear approximation space.
            %
            %   U = fit(obj, F) fits data within the function space by
            %   solving the fitting system \f$D*U = F\f$, where \f$D\f$ is
            %   the design matrix. Handles well-posed and overdetermined
            %   systems using least squares when necessary.
            %
            % Inputs:
            %   obj - The LinearApproximator object
            %   F - Data in the latent space (m x ... array)
            %
            % Outputs:
            %   U - Coefficients of the basis functions (m x ... array)
            
            D = obj.design;
            core.except.assert(~isempty(D), 'EmptyDesign', ...
                'The design operator has not been set. Call setDesign() first.');

            m = size(D, 1);
            core.except.assert(size(F, 1) == m, 'InvalidInput', ...
                'Data dimension (%d) must match design operator rows (%d).', ...
                size(F, 1), m);
                        
            s = size(F);
            F = reshape(F, m, []);
            
            if obj.isWellPosed(D)
                U = D \ F;
            elseif obj.isOverdetermined(D)
                U = (D.' * D) \ (D.' * F);
            else
                core.except.assert(0, 'UnderdeterminedSystem', ...
                    'Design operator cannot be underdetermined.');
            end
            
            U = reshape(U, [size(U, 1), s(2:end)]);
        end
        
        function U = project(obj, F)
            % PROJECT Project data onto the linear approximation space.
            %
            %   U = project(obj, F) projects data from the latent space
            %   onto the function space by solving the projection system
            %   \f$ M*U = F \f$, where \f$M\f$ is the mass operator.
            %
            % Inputs:
            %   obj - The LinearApproximator object
            %   F - Data in the latent space (m x ... array)
            %
            % Outputs:
            %   U - Coefficients of the basis functions (m x ... array)
            
            M = obj.mass;
            core.except.assert(~isempty(M), 'EmptyMass', ...
                'The mass operator has not been set. Call setMass() first.');

            m = size(M, 1);
            core.except.assert(size(F, 1) == m, 'InvalidInput', ...
                'Data dimension (%d) must match mass operator rows (%d).', ...
                size(F, 1), m);
                        
            s = size(F);
            F = reshape(F, m, []);
            
            if obj.isWellPosed(M)
                U = M \ F;
            elseif obj.isOverdetermined(M)
                U = (M.' * M) \ (M.' * F);
            else
                core.except.assert(0, 'UnderdeterminedSystem', ...
                    'Mass operator cannot be underdetermined.');
            end
            
            U = reshape(U, [size(U, 1), s(2:end)]);
        end
    end
    
    methods (Static)
        function TF = isWellPosed(A)
            % ISWELLPOSED Check if a matrix represents a well-posed system.
            %
            %   TF = isWellPosed(A) returns true if the matrix A is square
            %   (number of rows equals number of columns), indicating a
            %   well-posed system with a unique solution.
            %
            % Inputs:
            %   A - Matrix to check (numeric matrix)
            %
            % Outputs:
            %   TF - True if system is well-posed (logical)
            
            core.except.assert(~isempty(A), 'EmptyMatrix', ...
                'Matrix cannot be empty.');
            
            TF = size(A, 1) == size(A, 2);
        end
        
        function TF = isOverdetermined(A)
            % ISOVERDETERMINED Check if a matrix represents an
            % overdetermined system.
            %
            %   TF = isOverdetermined(A) returns true if the matrix A has
            %   more rows than columns, indicating an overdetermined system
            %   that requires least squares solution.
            %
            % Inputs:
            %   A - Matrix to check (numeric matrix)
            %
            % Outputs:
            %   TF - True if system is overdetermined (logical
            %
            % Examples:
            %   % Check overdetermined system
            %   A = randn(10, 5);  % More equations than unknowns
            %   isOver = LinearApproximator.isOverdetermined(A);  % true
            %   
            %   % Check square system
            %   B = randn(5, 5);
            %   isOver = LinearApproximator.isOverdetermined(B);  % false
            
            core.except.assert(~isempty(A), 'EmptyMatrix', ...
                'Matrix cannot be empty.');
            
            TF = size(A, 1) > size(A, 2);
        end
        
        function TF = isUnderdetermined(A)
            % ISUNDERDETERMINED Check if a matrix represents an
            % underdetermined system.
            %
            %   TF = isUnderdetermined(A) returns true if the matrix A has
            %   fewer rows than columns, indicating an underdetermined
            %   system with multiple solutions.
            %
            % Inputs:
            %   A - Matrix to check (numeric matrix)
            %
            % Outputs:
            %   TF - True if system is underdetermined (logical)
            
            core.except.assert(~isempty(A), 'EmptyMatrix', ...
                'Matrix cannot be empty.');
            
            TF = size(A, 1) < size(A, 2);
        end
    end
    
    methods (Abstract)
        % SETDESIGN Set up the design operator for fitting.
        obj = setDesign(obj, varargin)

        % EMBED Embed data into the latent space.
        F = embed(obj, varargin)
    end
end