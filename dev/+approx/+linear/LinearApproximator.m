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
        Basis  % Basis functions that span the approximation space
        Mass   % Mass matrix for projection operations (square matrix)
        Design % Design matrix for fitting operations (matrix)
    end
    
    properties (Dependent)
        NDofs % Number of degrees of freedom (scalar)
    end
    
    methods
        function obj = LinearApproximator(basis)
            % LINEARAPPROXIMATOR Constructor for LinearApproximator.
            %
            %   obj = LinearApproximator(basis) creates a new
            %   LinearApproximator object with the specified @a basis
            %   functions. The mass and design operators are initialized as
            %   empty and must be set by calling setMass() and setDesign().
            
            arguments
                basis {mustBeNonempty}
            end
            
            obj.Basis = basis;
            obj.Mass = [];
            obj.Design = [];
        end
        
        function n = get.NDofs(obj)
            % GET.NDOFS Returns the number of degrees of freedom.
            
            n = obj.Basis.NCodims;
        end
        
        function obj = setMass(obj, X, w)
            % SETMASS Set up the mass operator using numerical integration.
            %
            %   obj = setMass(obj, X, w) constructs the mass matrix using
            %   integration points and weights. The (i,j) entry of the mass
            %   matrix is \f$\langle p_i, p_j\rangle\f$ where the inner
            %   product is defined by the weighted integration:
            %
            %   \f[
            %       \langle f,g\rangle = \sum_k w_k f(x_k) g(x_k).
            %   \f]
            
            arguments
                obj approx.linear.LinearApproximator
                X {mustBeNumeric, mustBeNonempty}
                w (1, :) {mustBeNumeric, mustBeNonempty}
            end
            
            nd = obj.Basis.NDims;
            X = reshape(X, size(X, 1), []);
            
            core.except.assert(size(X, 1) >= nd, 'InvalidInput', ...
                'Integration points must have at least %d rows.', nd);
            
            X = X(1:nd, :);
            w = w(:);
            
            nq = length(w);
            core.except.assert(size(X, 2) == nq, 'InvalidInput', ...
                'Number of integration points (%d) must match number of weights (%d).', ...
                size(X, 2), nq);
            
            B = obj.Basis.eval(X);
            W = diag(w);
            M = B * W * B.';
            M = core.linalg.zeroing(M);
            
            obj.Mass = M;
        end
        
        function obj = setDesign(obj, varargin)
            % SETDESIGN Set up the design operator for fitting.
            
            arguments
                obj approx.linear.LinearApproximator
            end
            
            arguments (Repeating)
                varargin
            end
            
            core.except.assert(0, 'NotImplemented', ...
                'Method setDesign is not implemented by class j%s.', class(obj));
        end
        
        function F = embed(obj, varargin)
            % EMBED Embed data into the latent space.
            
            arguments
                obj approx.linear.LinearApproximator
            end
            
            arguments (Repeating)
                varargin
            end
            
            F = []; % Explicitly assign output to avoid unset variable error
            core.except.assert(false, 'NotImplemented', ...
                'Method embed is not implemented by class %s.', class(obj));
        end
        
        function Y = eval(obj, X, C)
            % EVAL Evaluate linear approximations at given points.
            %
            %   Y = eval(obj, X, C) evaluates functions with coefficient
            %   matrix @a C at evaluation points @a X.
            
            arguments
                obj approx.linear.LinearApproximator
                X {mustBeNumeric, mustBeNonempty}
                C {mustBeNumeric} = []
            end
            
            X = reshape(X, size(X, 1), []);
            C = reshape(C, size(C, 1), []);
            B = obj.Basis.eval(X);
            Y = C.' * B;
        end
        
        function Y = diff(obj, X, C, order)
            % DIFF Differentiate linear approximations at given points.
            %
            %   Y = diff(obj, X, C, order) differentiates functions with
            %   coefficient matrix @a C at evaluation points @a X.
            
            arguments
                obj approx.linear.LinearApproximator
                X {mustBeNumeric, mustBeNonempty}
                C {mustBeNumeric, mustBeNonempty}
                order (1, :) {mustBeNumeric, mustBeNonempty}
            end
            
            X = reshape(X, size(X, 1), []);
            C = reshape(C, size(C, 1), []);
            B = obj.Basis.diff(X, order);
            Y = C.' * B;
        end
        
        function U = fit(obj, F)
            % FIT Fit data within the linear approximation space.
            %
            %   U = fit(obj, F) fits data within the function space by
            %   solving the fitting system \f$D*U = F\f$, where \f$D\f$ is
            %   the design matrix. Handles well-posed and overdetermined
            %   systems using least squares when necessary.
            
            arguments
                obj approx.linear.LinearApproximator
                F {mustBeNumeric, mustBeNonempty}
            end
            
            D = obj.Design;
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
            
            arguments
                obj approx.linear.LinearApproximator
                F {mustBeNumeric, mustBeNonempty}
            end
            
            M = obj.Mass;
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
            %   TF = isWellPosed(A) returns true if the matrix @a A is
            %   square (number of rows equals number of columns),
            %   indicating a well-posed system with a unique solution.
            
            arguments
                A {mustBeNumeric, mustBeNonempty}
            end
            
            TF = size(A, 1) == size(A, 2);
        end
        
        function TF = isOverdetermined(A)
            % ISOVERDETERMINED Check if a matrix represents an
            % overdetermined system.
            %
            %   TF = isOverdetermined(A) returns true if the matrix @a A
            %   has more rows than columns, indicating an overdetermined
            %   system that requires least squares solution.
            
            arguments
                A {mustBeNumeric, mustBeNonempty}
            end
            
            TF = size(A, 1) > size(A, 2);
        end
        
        function TF = isUnderdetermined(A)
            % ISUNDERDETERMINED Check if a matrix represents an
            % underdetermined system.
            %
            %   TF = isUnderdetermined(A) returns true if the matrix @a A
            %   has fewer rows than columns, indicating an underdetermined
            %   system with multiple solutions.
            
            arguments
                A {mustBeNumeric, mustBeNonempty}
            end
            
            TF = size(A, 1) < size(A, 2);
        end
    end
end