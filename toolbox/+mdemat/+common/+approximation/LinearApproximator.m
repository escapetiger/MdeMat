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
    %       \langle v_h, u_h \rangle = \langle v_h,u\rangle,
    %       \quad \forall v_h \in V_h,
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
            %       \langle f,g\rangle = \sum_k w_k \overline{f(x_k)} g(x_k).
            %   \f]
            
            arguments
                obj mdemat.common.approximation.LinearApproximator
                X {mustBeNumeric, mustBeNonempty}
                w (1, :) {mustBeNumeric, mustBeNonempty}
            end
            
            nd = obj.Basis.NDims;
            X = reshape(X, size(X, 1), []);
            
            assert(all(size(X, 1) >= nd, 'all'), 'mdemat:common:approximation:LinearApproximator:InvalidInput', ...
                'Integration points must have at least %d rows.', nd);
            
            X = X(1:nd, :);
            w = w(:);
            
            nq = length(w);
            assert(all(size(X, 2) == nq, 'all'), 'mdemat:common:approximation:LinearApproximator:InvalidInput', ...
                'Number of integration points (%d) must match number of weights (%d).', ...
                size(X, 2), nq);
            
            B = obj.Basis.eval(X);
            M = (conj(B) .* w.') * B.';
            M = mdemat.common.linalg.zeroing(M);
            
            obj.Mass = M;
        end
        
        function obj = setDesign(obj, varargin)
            % SETDESIGN Set up the design operator for fitting.
            
            arguments
                obj mdemat.common.approximation.LinearApproximator
            end
            
            arguments (Repeating)
                varargin
            end
            
            assert(all(0, 'all'), 'mdemat:common:approximation:LinearApproximator:NotImplemented', ...
                'Method setDesign is not implemented by class %s.', class(obj));
        end
        
        function F = embed(obj, varargin)
            % EMBED Embed data into the latent space.
            
            arguments
                obj mdemat.common.approximation.LinearApproximator
            end
            
            arguments (Repeating)
                varargin
            end
            
            F = []; % Explicitly assign output to avoid unset variable error
            assert(all(false, 'all'), 'mdemat:common:approximation:LinearApproximator:NotImplemented', ...
                'Method embed is not implemented by class %s.', class(obj));
        end
        
        function Y = eval(obj, X, C)
            % EVAL Evaluate linear approximations at given points.
            %
            %   Y = eval(obj, X, C) evaluates functions with coefficient
            %   matrix @a C at evaluation points @a X.
            
            arguments
                obj mdemat.common.approximation.LinearApproximator
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
                obj mdemat.common.approximation.LinearApproximator
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
                obj mdemat.common.approximation.LinearApproximator
                F {mustBeNumeric, mustBeNonempty}
            end
            
            D = obj.Design;
            assert(all(~isempty(D), 'all'), 'mdemat:common:approximation:LinearApproximator:EmptyDesign', ...
                '%s', 'The design operator has not been set. Call setDesign() first.');
            
            m = size(D, 1);
            assert(all(size(F, 1) == m, 'all'), 'mdemat:common:approximation:LinearApproximator:InvalidInput', ...
                'Data dimension (%d) must match design operator rows (%d).', ...
                size(F, 1), m);
            
            s = size(F);
            F = reshape(F, m, []);
            
            assert(all(size(D, 1) >= size(D, 2), 'all'), ...
                'mdemat:common:approximation:LinearApproximator:UnderdeterminedSystem', '%s', 'Design operator cannot be underdetermined.');
            % Backslash uses QR for rectangular least squares, avoiding
            % the loss of conditioning caused by normal equations.
            U = D \ F;
            
            U = reshape(U, [size(U, 1), s(2:end)]);
        end
        
        function U = project(obj, F)
            % PROJECT Project data onto the linear approximation space.
            %
            %   U = project(obj, F) projects data from the latent space
            %   onto the function space by solving the projection system
            %   \f$ M*U = F \f$, where \f$M\f$ is the mass operator.
            
            arguments
                obj mdemat.common.approximation.LinearApproximator
                F {mustBeNumeric, mustBeNonempty}
            end
            
            M = obj.Mass;
            assert(all(~isempty(M), 'all'), 'mdemat:common:approximation:LinearApproximator:EmptyMass', ...
                '%s', 'The mass operator has not been set. Call setMass() first.');
            
            m = size(M, 1);
            assert(all(size(F, 1) == m, 'all'), 'mdemat:common:approximation:LinearApproximator:InvalidInput', ...
                'Data dimension (%d) must match mass operator rows (%d).', ...
                size(F, 1), m);
            
            s = size(F);
            F = reshape(F, m, []);
            
            assert(all(size(M, 1) >= size(M, 2), 'all'), ...
                'mdemat:common:approximation:LinearApproximator:UnderdeterminedSystem', '%s', 'Mass operator cannot be underdetermined.');
            U = M \ F;
            
            U = reshape(U, [size(U, 1), s(2:end)]);
        end
    end
    
    methods (Static)
        function TF = isWellPosed(A)
            % ISWELLPOSED Check the square-system shape used by this API.
            %
            %   TF = isWellPosed(A) returns true if the matrix @a A is
            %   square (number of rows equals number of columns),
            %   without testing rank or guaranteeing a unique solution.
            
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
