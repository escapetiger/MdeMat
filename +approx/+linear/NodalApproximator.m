classdef NodalApproximator < approx.linear.LinearApproximator
    % NODALAPPROXIMATOR Linear approximator based on nodal information.
    %
    %   NodalApproximator implements a linear approximation method that
    %   uses interpolation basis functions, such as Lagrange polynomials.
    %
    %   Let \f$x_j\f$ denote the nodal points, \f$p_i\f$ basis functions,
    %   and \f$c_i\f$ the coefficients to be determined. The projection
    %   system is defined as:
    %
    %   \f[
    %       \sum_i\sum_kw_kc_ip_i(x_k)p_j(x_k) = f_j(x_k), \forall j.
    %   \f]
    %
    %   This implies that the mass matrix is \f$M=diag(w_k)\f$. The
    %   fitting system is defined as:
    %
    %   \f[
    %       u(x_j) = \sum_i c_i p_i(x_j), \forall j.
    %   \f]
    % 
    %   The design matrix \f$D_{ji} = p_i(x_j)\f$ represents the evaluation
    %   of basis functions at nodal points. For Lagrangian polynomial,
    %   \f$D\f$ reduces to the identity matrix.
    %
    % Examples:
    %   % Create nodal approximator with Lagrange basis
    %   basis = approx.basis.InterpolationBasisFunction(5, 'lagrange', 'unit', 'gauss_legendre');
    %   basis = basis.autoLoad();
    %   approx = approx.linear.NodalApproximator(basis);
    %   
    %   % Set up design operator with nodal points
    %   X = linspace(-1, 1, 5);  % 5 nodal points
    %   approx.setDesign(X);
    %   
    %   % Approximate a function via interpolation
    %   U = sin(pi * X);  % Function values at nodal points
    %   F = approx.embed(U);
    %   coeffs = approx.fit(F);
    %   
    %   % For Lagrange basis, coeffs should equal U for exact interpolation
    %   interpolation_error = norm(coeffs - U(:));
    %   
    %   % Evaluate approximation at fine grid
    %   X_fine = linspace(-1, 1, 100);
    %   Y = approx.evaluate(X_fine, coeffs);
    %
    % See also:
    %   approx.linear.LinearApproximator, approx.linear.ModalApproximator

    properties (Constant)
        TYPE = 'nodal' % Approximator type identifier
    end

    methods
        function obj = setDesign(obj, X)
            % SETDESIGN Set up the design operator using nodal points.
            %
            %   obj = setDesign(obj, X) constructs the design matrix as the
            %   evaluation of basis functions at the nodal points. The
            %   design matrix \f$D_{ji} = p_i(x_j)\f$ where \f$x_j\f$ are 
            %   nodal points and \f$p_i\f$ are basis functions.
            %
            % Inputs:
            %   obj - The NodalApproximator object
            %   X - Nodal points (n x q matrix)
            %
            % Outputs:
            %   obj - The NodalApproximator object

            n = obj.basis.nDims;
            X = reshape(X, size(X, 1), []);
            
            core.except.assert(size(X, 1) >= n, 'InvalidInput', ...
                'Nodal points must have at least %d rows for %dD basis.', n, n);
            
            X = X(1:n, :);
            
            core.except.assert(size(X, 2) > 0, 'InvalidInput', ...
                'At least one nodal point is required.');
            
            D = obj.basis.evaluate(X).';
            D = core.linalg.zeroing(D);

            obj.design = D;
        end

        function F = embed(obj, U)
            % EMBED Embed function values at nodal points into the latent
            % space.
            %
            %   F = embed(obj, U) embeds function values at nodal points
            %   into the latent space. For nodal embedding, this operation
            %   is essentially the identity transformation since the latent
            %   space directly corresponds to function values at nodes.
            %
            % Inputs:
            %   obj - The NodalApproximator object
            %   U - Function values at nodal points (m x ... array)
            %
            % Outputs:
            %   F - Latent space representation (m x ... array), same as input U

            D = obj.design;
            core.except.assert(~isempty(D), 'EmptyDesign', ...
                'The design operator has not been set. Call setDesign() first.');

            m = size(D, 1);
            F = reshape(U, m, []);
        end
    end
end