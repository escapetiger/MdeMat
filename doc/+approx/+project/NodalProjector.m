classdef NodalProjector < approx.project.GalerkinProjector
    % NODALPROJECTOR Galerkin projector based on nodal embedding.
    %
    %   NodalProjector implements a Galerkin projection method that
    %   enforces exact interpolation at specified nodal points. This
    %   approach is particularly suitable for interpolation basis
    %   functions such as Lagrange polynomials.
    %
    %   The Galerkin system is defined as:
    %   \f[
    %   u(x_j) = \sum_{i=1}^{n} c_i p_i(x_j), \forall j = 1, ..., m,
    %   \f]
    %   where \f$x_j\f$ are the nodal points, \f$p_i\f$ are basis
    %   functions, and \f$c_i\f$ are the coefficients to be determined. 
    %   The mass matrix \f$M_{ji} = p_i(x_j)\f$ represents the evaluation
    %   of basis functions at nodal points.
    %
    %   This method is exact for interpolation problems and provides a
    %   direct relationship between function values at nodes and basis
    %   coefficients.
    %
    % Examples:
    %   % Create nodal projector with Lagrange basis
    %   basis = InterpolationBasisFunction(5, 'lagrange', 'unit', 'gauss_legendre');
    %   basis = basis.autoLoad();
    %   projector = NodalProjector(basis);
    %   
    %   % Set up mass operator with nodal points
    %   X = linspace(0, 1, 5);  % 5 nodal points
    %   projector.setMass(X);
    %   
    %   % Project function data
    %   U = sin(pi * X);  % Function values at nodal points
    %   F = projector.embed(U);
    %   coeffs = projector.project(F);
    %   
    %   % For Lagrange basis, coeffs should equal U for exact interpolation
    %
    % See Also:
    %   approx.project.GalerkinProjector, approx.project.ModalProjector

    properties (Constant)
        TYPE = 'nodal'
    end

    methods
        function obj = setMass(obj, X)
            % SETMASS Set up the mass operator using nodal points.
            %
            %   obj = setMass(obj, X) constructs the mass matrix as the
            %   evaluation of basis functions at the nodal points. The mass
            %   matrix \f$M_{ji} = p_i(x_j)\f$ where \f$x_j\f$ are nodal
            %   points and \f$p_i\f$ are basis functions.
            %
            % Inputs:
            %   obj - The NodalProjector object
            %   X - Nodal points (n x q)
            %
            % Outputs:
            %   obj - The NodalProjector object

            n = obj.basis.nDims;
            X = reshape(X, size(X, 1), []);
            
            core.except.assert(size(X, 1) >= n, 'InvalidInput', ...
                'Nodal points must have at least %d rows for %dD basis.', n, n);
            
            X = X(1:n, :);
            
            core.except.assert(size(X, 2) > 0, 'InvalidInput', ...
                'At least one nodal point is required.');
            
            M = obj.basis.evaluate(X).';
            M = core.linalg.zeroing(M);

            obj.mass = M;
        end

        function F = embed(obj, U)
            % EMBED Embed function values at nodal points into the latent
            % space.
            %
            %   F = embed(obj, U) embeds function values at nodal points
            %   into the latent space. For nodal embedding, this simply
            %   reshapes the input values since the latent space directly
            %   corresponds to function values at nodes.
            %
            % Inputs:
            %   obj - The NodalProjector object
            %   U - Function values at nodal points (m x q)
            %
            % Outputs:
            %   F - Latent space representation (m x q), same as input U

            M = obj.mass;
            core.except.assert(~isempty(M), 'EmptyMass', ...
                'The mass operator has not been set. Call setMass() first.');

            m = size(M, 1);
            core.except.assert(size(U, 1) == m, 'InvalidInput', ...
                'Function values rows (%d) must match number of nodal points (%d).', ...
                size(U, 1), m);

            F = reshape(U, m, []);
        end
    end
end