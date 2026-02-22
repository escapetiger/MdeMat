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
    %
    % See also:
    %   approx.linear.LinearApproximator, approx.linear.NodalApproximator
    
    properties (Constant)
        Type = 'modal' % Approximator type identifier
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
            % 
            %   \f[
            %       \langle f,g\rangle = \sum_k w_k f(x_k) g(x_k)
            %   \f].
            
            arguments
                obj approx.linear.ModalApproximator
                X {mustBeNumeric, mustBeNonempty}
                w (1, :) {mustBeNumeric}
            end
            
            nd = obj.Basis.NDims;
            X = reshape(X, size(X, 1), []);
            
            core.except.assert(size(X, 1) >= nd, 'InvalidInput', ...
                'Integration points must have at least %d rows.', nd);
            
            X = X(1:nd, :);
            
            np = size(X, 2);
            nq = length(w);
            core.except.assert(np == nq, 'InvalidInput', ...
                'Number of integration points (%d) must match number of weights (%d).', ...
                np, nq);
            
            B = obj.Basis.eval(X);
            W = diag(w);
            D = B * W * B.';
            D = core.linalg.zeroing(D);
            
            obj.Design = D;
        end
        
        function F = embed(obj, U, B, w)
            % EMBED Embed data into the latent space using modal approach.
            %
            %   F = embed(obj, U, B, w) embeds function values into the
            %   latent space using the modal approach. Computes the
            %   weighted inner product between basis functions and function
            %   values at integration points:
            %
            %   \f[
            %       F_i = \langle U, p_i\rangle = \sum_k w_k U_k p_i(x_k)
            %   \f]
            
            arguments
                obj approx.linear.ModalApproximator
                U {mustBeNumeric, mustBeNonempty}
                B {mustBeNumeric, mustBeNonempty}
                w (1, :) {mustBeNumeric}
            end
            
            D = obj.Design;
            core.except.assert(~isempty(D), 'EmptyDesign', ...
                'The design operator has not been set. Call setDesign() first.');
            
            nq = length(w);
            U = reshape(U, nq, []);
            B = reshape(B, size(D, 1), nq);
            W = diag(w);
            F = B * W * U;
        end
    end
end