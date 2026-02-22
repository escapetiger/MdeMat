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
    %
    % See also:
    %   approx.linear.LinearApproximator, approx.linear.ModalApproximator

    properties (Constant)
        Type = 'nodal' % Approximator type identifier
    end

    methods
        function obj = setDesign(obj, X)
            % SETDESIGN Set up the design operator using nodal points.
            %
            %   obj = setDesign(obj, X) constructs the design matrix as the
            %   evaluation of basis functions at the nodal points. The
            %   design matrix \f$D_{ji} = p_i(x_j)\f$ where \f$x_j\f$ are 
            %   nodal points and \f$p_i\f$ are basis functions.

            arguments
                obj approx.linear.NodalApproximator
                X {mustBeNumeric, mustBeNonempty}
            end

            nd = obj.Basis.NDims;
            X = reshape(X, size(X, 1), []);
            
            core.except.assert(size(X, 1) >= nd, 'InvalidInput', ...
                'Nodal points must have at least %d rows for %dD basis.', nd, nd);
            
            X = X(1:nd, :);
            
            core.except.assert(size(X, 2) > 0, 'InvalidInput', ...
                'At least one nodal point is required.');
            
            D = obj.Basis.eval(X).';
            D = core.linalg.zeroing(D);

            obj.Design = D;
        end

        function F = embed(obj, U)
            % EMBED Embed function values at nodal points into the latent
            % space.
            %
            %   F = embed(obj, U) embeds function values at nodal points
            %   into the latent space. For nodal embedding, this operation
            %   is essentially the identity transformation since the latent
            %   space directly corresponds to function values at nodes.

            arguments
                obj approx.linear.NodalApproximator
                U {mustBeNumeric, mustBeNonempty}
            end

            D = obj.Design;
            core.except.assert(~isempty(D), 'EmptyDesign', ...
                'The design operator has not been set. Call setDesign() first.');

            F = reshape(U, size(D, 1), []);
        end
    end
end