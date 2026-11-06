classdef RbfKernelBasis < mdemat.data.function.Function
    % RBFKERNELBASIS Radial basis function kernel interpolation basis.
    %
    %   RbfKernelBasis represents a set of radial basis function (RBF)
    %   interpolation basis functions. Each basis function is a radial
    %   kernel centered at a node point, providing scattered data
    %   interpolation capabilities with various kernel types.
    %
    %   Supported kernel types:
    %
    %   - 'gaussian': \f$\exp(-\varepsilon^2 r^2)\f$
    %
    %   - 'multiquadric': \f$\sqrt{1 + \varepsilon^2 r^2}\f$
    %
    %   - 'inverse_multiquadric': \f$\frac{1}{\sqrt{1 + \varepsilon^2 r^2}}\f$
    %
    % See Also:
    %   mdemat.data.function.Function

    properties
        Nodes % Node points (nd x nc)
        Kernel % Type of RBF kernel
        Epsilon % Shape parameter (scalar or vector)
    end

    methods
        function obj = RbfKernelBasis(nodes, options)
            % RBFKERNELBASIS Constructor for RbfKernelBasis.
            %
            %   obj = RbfKernelBasis(nodes, kernel=kernel, epsilon=epsilon)
            %   creates an RbfKernelBasis with specified shape parameter.

            arguments
                nodes{mustBeNumeric}
                options.kernel {mustBeMember(options.kernel, ...
                    {'gaussian', 'multiquadric', 'inverse_multiquadric'})} = 'gaussian'
                options.epsilon(:, 1) {mustBeNumeric, mustBeReal, mustBeFinite, ...
                    mustBePositive} = 1
            end

            nDims = size(nodes, 1);
            nCodims = size(nodes, 2);
            assert(isscalar(options.epsilon) || numel(options.epsilon) == nCodims, ...
                'MdeMat:Rbf:InvalidShapeParameter', ...
                'Specify one positive shape parameter or one per kernel.');
            obj@mdemat.data.function.Function(nDims=nDims, nCodims=nCodims);
            obj.Nodes = nodes;
            obj.Kernel = options.kernel;
            obj.Epsilon = options.epsilon;
        end
    end

    methods (Access = protected)
        function Y = evalImpl(obj, X)
            % EVALUATEIMPL Implementation of RBF basis evaluation.
            %
            %   Y = evalImpl(obj, X) evals the RBF basis functions
            %   at the specified points. Each basis function is a radial
            %   kernel centered at the corresponding node point.
            %
            % Inputs:
            %   obj - The RbfKernelBasis object
            %   X - Coordinates (nd x np)
            %
            % Outputs:
            %   Y - Basis function values (nc x np)

            nd = obj.NDims;
            nc = obj.NCodims;
            np = size(X, 2);

            X0 = reshape(obj.Nodes, [nd, nc, 1]); %< [nd, nc, 1]
            X = reshape(X, [nd, 1, np]); %< [nd, 1, np]
            R = vecnorm(X-X0, 2, 1); % [1, nc, np]
            R = reshape(R, [nc, np]); % [1, nc, np]

            E = obj.Epsilon;
            if isscalar(E), E = repmat(E, nc, 1); end
            E = repmat(E, [1, np]);

            Y = obj.evalKernel(R, E);
        end

        function dY = diffImpl(obj, X, order)
            % DIFFERENTIATEIMPL Implementation of RBF basis derivative.
            %
            %   dY = diffImpl(obj, X, order) computes derivatives of the
            %   RBF basis functions using the chain rule and
            %   kernel-specific derivative formulas.
            %
            %   The derivative is computed by
            %
            %   \f[
            %       \frac{\partial \phi_i(\|x - c_i\|)}{\partial x_j} =
            %       \phi'_i(r_i) \frac{x_j - c_{i,j}}{r_i},
            %   \f]
            %
            %   where \f$r_i = \|x - c_i\|\f$ is the distance to center
            %   \f$i\f$.

            assert(sum(order) == 1 && all(order == fix(order)), ...
                'MdeMat:Rbf:UnsupportedDerivative', 'Only first partial derivatives are supported.');
            m = find(order > 0, 1);
            nd = obj.NDims;
            nc = obj.NCodims;
            np = size(X, 2);

            X = reshape(X, [nd, 1, np]); % [nd, 1, np]
            X0 = reshape(obj.Nodes, [nd, nc, 1]); % [nd, nc, 1]
            displacement = X-X0;
            R = vecnorm(displacement, 2, 1); % [1, nc, np]
            R = reshape(R, [nc, np]); % [1, nc, np]

            E = obj.Epsilon;
            if isscalar(E), E = repmat(E, nc, 1); end
            E = repmat(E, [1, np]);

            dY = obj.diffKernel(R, E);

            R0 = R;
            R0(R == 0) = eps;
            component = reshape(displacement(m, :, :), nc, np);
            dY = dY .* component ./ R0;
            dY(R == 0) = 0;
        end

        function phi = evalKernel(obj, r, epsilon)
            % EVALKERNEL Evaluate RBF kernel at given distances.
            %
            %   phi = evalKernel(obj, r, epsilon) evals the specified RBF
            %   kernel at the given distances with shape parameter epsilon.
            %   Supports vectorized evaluation for multiple distances and
            %   shape parameters.

            switch obj.Kernel
                case 'gaussian'
                    phi = exp(-(epsilon .* r).^2);

                case 'multiquadric'
                    phi = sqrt(1+(epsilon .* r).^2);

                case 'inverse_multiquadric'
                    phi = 1 ./ sqrt(1+(epsilon .* r).^2);
            end
        end

        function dphi = diffKernel(obj, r, epsilon)
            % DIFFKERNEL Evaluate RBF kernel derivative.
            %
            %   dphi = diffKernel(obj, r, epsilon) evaluates the derivative
            %   of the RBF kernel with respect to distance \f$r\f$:
            %   \f$\frac{\partial \phi(r)}{\partial r}\f$. Supports
            %   vectorized evaluation for multiple distances and shape
            %   parameters.

            switch obj.Kernel
                case 'gaussian'
                    dphi = -2 * epsilon.^2 .* r .* exp(-(epsilon .* r).^2);

                case 'multiquadric'
                    dphi = epsilon.^2 .* r ./ sqrt(1+(epsilon .* r).^2);

                case 'inverse_multiquadric'
                    denom = (1 + (epsilon .* r).^2).^(3 / 2);
                    dphi = -epsilon.^2 .* r ./ denom;
            end
        end
    end
end
