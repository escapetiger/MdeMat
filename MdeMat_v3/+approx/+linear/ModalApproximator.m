classdef ModalApproximator < approx.linear.LinearApproximator
    % MODALAPPROXIMATOR Linear approximator based on modal information.
    %
    %   ModalApproximator implements a linear approximation method that
    %   uses modal basis functions, such as Legendre, Chebyshev, or Fourier
    %   polynomials, etc. The underlying projection and fitting operators
    %   coincide in this case. The mass matrix \f$M_{ij} = \langle p_i,
    %   p_j\rangle\f$ is constructed using numerical integration with
    %   specified quadrature points and weights.
    %
    % Examples:
    %   % Create modal approximator with Legendre basis
    %   basis = approx.basis.OrthogonalBasisFunction(5, 'legendre', 'unit');
    %   basis = basis.autoLoad();
    %   approx = approx.linear.ModalApproximator(basis);
    %   
    %   % Set up mass operator with Gauss-Legendre quadrature
    %   [X, w] = gauss_legendre_quadrature(10);  % 10 quadrature points
    %   approx.setMass(X, w);
    %   approx.setDesign(X, w);
    %   
    %   % Approximate a function
    %   U = sin(pi * X);  % Function values at quadrature points
    %   B = basis.evaluate(X);  % Basis values at quadrature points
    %   F = approx.embed(U, B, w);
    %   coeffs = approx.fit(F);
    %   
    %   % Evaluate approximation at new points
    %   X_new = linspace(-1, 1, 100);
    %   Y = approx.evaluate(X_new, coeffs);
    %
    % See also:
    %   approx.linear.LinearApproximator, approx.linear.NodalApproximator

    properties (Constant)
        TYPE = 'modal' % Approximator type identifier
    end

    methods        
        function obj = setDesign(obj, X, w)
            % SETDESIGN Set up the design operator using numerical
            % integration.
            %
            %   obj = setDesign(obj, X, w) constructs the design matrix
            %   using integration points and weights. The (i,j) entry of
            %   the design matrix is \f$\langle p_i, p_j\rangle\f$ where
            %   the inner product is defined by the weighted integration:
            %   \f$\langle f,g\rangle = \sum_k w_k f(x_k) g(x_k)\f$.
            %
            % Inputs:
            %   obj - The ModalApproximator object
            %   X - Integration points (n x q matrix)
            %   w - Integration weights (1 x q vector or q x 1 vector)
            %
            % Outputs:
            %   obj - The updated ModalApproximator object
            
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
            D = B * W * B.';
            D = core.linalg.zeroing(D);     
            
            obj.design = D;
        end
        
        function F = embed(obj, U, B, w)
            % EMBED Embed data into the latent space using modal approach.
            %
            %   F = embed(obj, U, B, w) embeds function values into the
            %   latent space using the modal approach. Computes the
            %   weighted inner product between basis functions and function
            %   values at integration points:
            %   \f$F_i = \langle U, p_i\rangle = \sum_k w_k U_k p_i(x_k)\f$.
            %
            % Inputs:
            %   obj - The ModalApproximator object
            %   U - Function values at integration points (q x p matrix)
            %   B - Basis function values at integration points (m x q matrix)
            %   w - Integration weights (1 x q vector or q x 1 vector)
            %
            % Outputs:
            %   F - Latent space representation (m x p matrix)
            
            D = obj.design;
            core.except.assert(~isempty(D), 'EmptyDesign', ...
                'The design operator has not been set. Call setDesign() first.');
            
            m = size(D, 1);
            q = length(w);
            U = reshape(U, q, []);
            B = reshape(B, m, q);
            W = diag(w(:));
            F = B * W * U;
        end
    end
end