classdef NodalApproximator < mdemat.common.approximation.LinearApproximator
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
    %       \sum_i c_i\sum_k w_k\overline{p_j(x_k)}p_i(x_k)
    %       = \sum_k w_k\overline{p_j(x_k)}u(x_k), \quad\forall j.
    %   \f]
    %
    %   For Lagrange basis functions and quadrature at the interpolation
    %   nodes, the mass matrix is \f$M=diag(w_k)\f$. The fitting system is:
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
    %   mdemat.common.approximation.LinearApproximator, mdemat.common.approximation.ModalApproximator

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
                obj mdemat.common.approximation.NodalApproximator
                X {mustBeNumeric, mustBeNonempty}
            end

            nd = obj.Basis.NDims;
            X = reshape(X, size(X, 1), []);
            
            assert(all(size(X, 1) >= nd, 'all'), 'mdemat:common:approximation:NodalApproximator:InvalidInput', ...
                'Nodal points must have at least %d rows for %dD basis.', nd, nd);
            
            X = X(1:nd, :);
            
            assert(all(size(X, 2) > 0, 'all'), 'mdemat:common:approximation:NodalApproximator:InvalidInput', ...
                '%s', 'At least one nodal point is required.');
            
            D = obj.Basis.eval(X).';
            D = mdemat.common.linalg.zeroing(D);

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
                obj mdemat.common.approximation.NodalApproximator
                U {mustBeNumeric, mustBeNonempty}
            end

            D = obj.Design;
            assert(all(~isempty(D), 'all'), 'mdemat:common:approximation:NodalApproximator:EmptyDesign', ...
                '%s', 'The design operator has not been set. Call setDesign() first.');

            F = reshape(U, size(D, 1), []);
        end
    end
end
