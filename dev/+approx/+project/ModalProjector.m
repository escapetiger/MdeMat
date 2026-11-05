classdef ModalProjector < approx.project.GalerkinProjector
    % MODALPROJECTOR Galerkin projector based on modal embedding.
    %
    %   ModalProjector implements a Galerkin projection method that uses
    %   inner products between basis functions to construct the mass
    %   operator. This approach is particularly suitable for orthogonal
    %   and well-conditioned basis functions.
    %
    %   The Galerkin system is defined as:
    %   \f[
    %   \langle u(x), p_j(x)\rangle = \sum_{i=1}^{n} c_i \langle p_i(x), p_j(x)\rangle, \forall j = 1, ..., m
    %   \f]
    %   where \f$\langle\cdot,\cdot\rangle\f$ denotes the inner product
    %   defined by integration weights, \f$p_i\f$ are basis functions, and
    %   \f$c_i\f$ are the coefficients to be determined.
    %
    %   The mass matrix \f$M_{ij} = \langle p_i, p_j\rangle\f$ is
    %   constructed using numerical integration with specified quadrature
    %   points and weights.
    %
    % Examples:
    %   % Create modal projector with Legendre basis
    %   basis = OrthogonalBasisFunction(5, 'legendre', 'unit');
    %   basis = basis.autoLoad();
    %   projector = ModalProjector(basis);
    %   
    %   % Set up mass operator with Gauss-Legendre quadrature
    %   [X, w] = gauss_legendre_quadrature(10);  % 10 quadrature points
    %   projector.setMass(X, w);
    %   
    %   % Project function data
    %   U = sin(pi * X);  % Function values at quadrature points
    %   B = basis.evaluate(X);  % Basis values at quadrature points
    %   F = projector.embed(U, B, w);
    %   coeffs = projector.project(F);
    %
    % See Also:
    %   approx.project.GalerkinProjector, approx.project.NodalProjector

    methods        
        function obj = setMass(obj, X, w)
            % SETMASS Set up the mass operator using numerical integration.
            %
            %   obj = setMass(obj, X, w) constructs the mass matrix using
            %   integration points and weights. The (i,j) entry of the mass
            %   matrix is \f$\langle p_i, p_j\rangle\f$ where the inner
            %   product is defined by the weighted integration: \f$\langle
            %   f,g\rangle = \sum_k w_k f(x_k) g(x_k)\f$.
            %
            % Inputs:
            %   obj - The ModalProjector object
            %   X - Integration points (n x q)
            %   w - Integration weights (1 x q)
            %
            % Outputs:
            %   obj - The ModalProjector object
            
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

        function F = embed(obj, U, B, w)
            % EMBED Embed data into the latent space using modal approach.
            %
            %   F = embed(obj, U, B, w) embeds function values into the
            %   latent space using the modal approach. Computes the weighted
            %   inner product between basis functions and function values
            %   at integration points: 
            %   \f$F_i = \langle U, p_i\rangle = \sum_k w_k U_k p_i(x_k)\f$.
            %
            % Inputs:
            %   obj - The ModalProjector object
            %   U - Function values at integration points (q x p) 
            %   B - Basis function values at integration points (m x q) 
            %   w - Integration weights (1 x q)
            %
            % Outputs:
            %   F - Latent space representation (m x p)
            
            M = obj.mass;
            core.except.assert(~isempty(M), 'EmptyMass', ...
                'The mass operator has not been set. Call setMass() first.');
            
            m = size(M, 1);
            core.except.assert(size(B, 1) == m, 'InvalidInput', ...
                'Basis function rows (%d) must match mass operator size (%d).', ...
                size(B, 1), m);

            q = length(w);
            core.except.assert(size(B, 2) == q, 'InvalidInput', ...
                'Basis function columns (%d) must match number of weights (%d).', ...
                size(B, 2), q);
            
            U = reshape(U, q, []);
            B = reshape(B, [], q);
            W = diag(w(:));
            F = B * W * U;
        end
    end
end