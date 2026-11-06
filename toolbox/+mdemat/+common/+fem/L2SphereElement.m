classdef L2SphereElement < mdemat.common.fem.L2Element
    % L2SPHEREELEMENT L2 element defined on spherical domains.
    %
    %   L2SphereElement extends L2Element for spherical geometries with
    %   various basis function types and reduction schemes. The element
    %   supports 1D (interval), 2D (circle), and 3D (sphere) domains with
    %   orthogonal, cardinal, or radial basis function (RBF) types.
    %
    % See also:
    %   mdemat.common.fem.L2Element

    properties
        DimReduction % Dimension reduction type ('topology' or 'symmetry')
    end

    methods
        function obj = L2SphereElement(geometry, approximator, integrator, dimReduction)
            % L2SPHEREELEMENT Construct L2SphereElement.
            %
            %   obj = L2SphereElement(geometry, approximator, integrator)
            %   creates an L2SphereElement with the specified geometry,
            %   approximator, and integrator.
            %
            % See also:
            %   mdemat.common.fem.L2Element

            obj@mdemat.common.fem.L2Element(geometry, approximator, integrator);
            obj.DimReduction = dimReduction;
        end
    end

    methods (Static)
        function obj = modal(nDims, order, options)
            % MODAL Create L2SphereElement with orthogonal basis.
            %
            %   obj = modal(nDims, order) creates an L2SphereElement using
            %   orthogonal basis functions (Legendre, Fourier, or Spherical
            %   Harmonics).
            %
            %   obj = modal(nDims, order, options) allows customization:
            %   - options.reduction: Reduction type ('topology' or 'symmetry', default: 'topology')
            %   - options.radius: Sphere radius (default: 1)

            arguments
                nDims{mustBePositive, mustBeInteger}
                order (1, 1) {mustBePositive, mustBeInteger}
                options.reduction{mustBeMember(options.reduction, {'', 'topology', 'symmetry'})} = 'topology'
                options.radius{mustBePositive} = 1
                options.np{mustBeNumeric} = []
            end

            radius = options.radius;
            reduction = options.reduction;

            %< Determine number of integration points
            if isempty(options.np)
                np = 0;
            else
                np = options.np;
            end
            if nDims == 1 && strcmp(reduction, 'topology')
                np = max(np, 2);
            elseif nDims == 1 && strcmp(reduction, 'symmetry')
                np = max(np, order+1);
            elseif nDims == 2 && strcmp(reduction, 'topology')
                np = max(np, 2*(order-1)+2);
            else
                np = max(np, [order+1, 2*(order-1)+2]);
            end

            %< Create geometry, integrator, and basis based on dimensions
            %< and reduction
            if nDims == 1 && strcmp(reduction, 'topology')
                A = 1 / 2;
                geometry = mdemat.data.geometry.Sphere([0], radius);
                R = mdemat.common.quadrature.GaussTrapezoidalRule(1, coord = 'car');
                E = mdemat.data.function.Constant(1, A);
                sqrtE = mdemat.data.function.Constant(1, sqrt(A));
                integrator = mdemat.common.quadrature.Integrator(R);
                integrator.setPoints(geometry, min(2, np));
                integrator.addWeightFunction(@(v) 1./E.eval(v));

                maxDegree = order - 1;
                F = mdemat.common.spectral.LegendreBasis(maxDegree, isNormalized = true);

            elseif nDims == 1 && strcmp(reduction, 'symmetry')
                A = 1 / (2 * radius);
                geometry = mdemat.data.geometry.Orthotope([-radius, radius]);
                R = mdemat.common.quadrature.GaussLegendreRule(1);
                E = mdemat.data.function.Constant(1, A);
                sqrtE = mdemat.data.function.Constant(1, sqrt(A));
                integrator = mdemat.common.quadrature.Integrator(R);
                integrator.setPoints(geometry, np);
                integrator.addWeightFunction(@(v) 1./E.eval(v));

                maxDegree = order - 1;
                F = mdemat.common.spectral.LegendreBasis(maxDegree, isNormalized = true);

            elseif nDims == 2 && strcmp(reduction, 'topology')
                A = 1 / (2 * pi * radius);
                geometry = mdemat.data.geometry.Sphere([0, 0], radius);
                R = mdemat.common.quadrature.GaussTrapezoidalRule(2, coord = 'cossph');
                E = mdemat.data.function.Constant(1, A);
                sqrtE = mdemat.data.function.Constant(1, sqrt(A));
                integrator = mdemat.common.quadrature.Integrator(R);
                integrator.setPoints(geometry, np);
                integrator.addWeightFunction(@(v) 1./E.eval(v));

                maxDegree = order - 1;
                F = mdemat.common.spectral.FourierBasis(maxDegree, isNormalized = true);

            else %< nDims >= 3
                A = 1 / (4 * pi * radius^2);
                geometry = mdemat.data.geometry.Sphere([0, 0, 0], radius);
                R = mdemat.common.quadrature.GaussTrapezoidalRule(3, coord = 'cossph');
                E = mdemat.data.function.Constant(2, A);
                sqrtE = mdemat.data.function.Constant(2, sqrt(A));
                integrator = mdemat.common.quadrature.Integrator(R);
                integrator.setPoints(geometry, np);
                integrator.addWeightFunction(@(v) 1./E.eval(v));

                maxDegree = order - 1;
                F = mdemat.common.spectral.SphericalHarmonicBasis(maxDegree, isNormalized = true);
            end

            %< Create basis function
            B = mdemat.data.function.Product(F, sqrtE);

            %< Create modal approximator
            approximator = mdemat.common.approximation.ModalApproximator(B);
            approximator.setMass(integrator.Nodes, integrator.Weights);
            approximator.setDesign(integrator.Nodes, integrator.Weights);

            %< Create element
            obj = mdemat.common.fem.L2SphereElement(geometry, approximator, integrator, reduction);
        end

        function obj = nodal(nDims, order, options)
            % NODAL Create L2SphereElement with RBF basis.
            %
            %   obj = nodal(nDims, order) creates an L2SphereElement using
            %   radial basis functions (RBF).
            %
            %   obj = nodal(nDims, order, options) allows customization:
            %   - options.reduction: Reduction type ('topology' or 'symmetry', default: 'topology')
            %   - options.radius: Sphere radius (default: 1)
            %   - options.kernel: RBF kernel type ('gaussian', default: 'gaussian')
            %   - options.epsilon: RBF shape parameter (default: 1)

            arguments
                nDims{mustBePositive, mustBeInteger}
                order {mustBePositive, mustBeInteger}
                options.reduction{mustBeMember(options.reduction, {'', 'topology', 'symmetry'})} = 'topology'
                options.radius{mustBePositive} = 1
                options.kernel{mustBeText} = 'gaussian'
                options.epsilon{mustBeNumeric} = []
            end

            radius = options.radius;
            reduction = options.reduction;

            %< Determine number of integration points
            np = order;
            epsilon = options.epsilon;

            %< Create geometry, integrator, and basis based on dimensions and reduction
            if nDims == 1 && strcmp(reduction, 'topology')
                % A = 1 / 2;
                geometry = mdemat.data.geometry.Sphere([0], radius);
                R = mdemat.common.quadrature.GaussTrapezoidalRule(1, coord = 'car');
                % E = mdemat.data.function.Constant(1, A);
                integrator = mdemat.common.quadrature.Integrator(R);
                integrator.setPoints(geometry, min(2, np));
                % integrator.addWeightFunction(@(v) 1./E.eval(v));

            elseif nDims == 1 && strcmp(reduction, 'symmetry')
                % A = 1 / (2 * radius);
                geometry = mdemat.data.geometry.Orthotope([-radius, radius]);
                R = mdemat.common.quadrature.GaussLegendreRule(1);
                % E = mdemat.data.function.Constant(1, A);
                integrator = mdemat.common.quadrature.Integrator(R);
                integrator.setPoints(geometry, np);
                % integrator.addWeightFunction(@(v) 1./E.eval(v));

            elseif nDims == 2 && strcmp(reduction, 'topology')
                % A = 1 / (2 * pi * radius);
                geometry = mdemat.data.geometry.Sphere([0, 0], radius);
                R = mdemat.common.quadrature.GaussTrapezoidalRule(2, coord = 'cossph');
                % E = mdemat.data.function.Constant(1, A);
                integrator = mdemat.common.quadrature.Integrator(R);
                integrator.setPoints(geometry, np);
                % integrator.addWeightFunction(@(v) 1./E.eval(v));

            else %< nDims >= 3
                % A = 1 / (4 * pi * radius^2);
                geometry = mdemat.data.geometry.Sphere([0, 0, 0], radius);
                R = mdemat.common.quadrature.GaussTrapezoidalRule(3, coord = 'cossph');
                % E = mdemat.data.function.Constant(1, A);
                integrator = mdemat.common.quadrature.Integrator(R);
                integrator.setPoints(geometry, np);
                % integrator.addWeightFunction(@(v) 1./E.eval(v));
            end

            %< Create RBF basis function
            if isempty(epsilon)
                nodes = integrator.Nodes;
                count = size(nodes, 2);
                differences = reshape(nodes, size(nodes, 1), count, 1) - ...
                    reshape(nodes, size(nodes, 1), 1, count);
                distances = reshape(vecnorm(differences, 2, 1), count, count);
                distances(1:count+1:end) = Inf;
                h = min(distances, [], 'all');
                if count == 1
                    h = 1;
                end
                assert(isfinite(h) && h > 0, 'MdeMat:SphereElement:DuplicateNodes', ...
                    'Velocity interpolation nodes must be distinct.');
                epsilon = 1e5 / h;
            else
                epsilon = options.epsilon;
            end
            B = mdemat.common.interpolation.RbfKernelBasis(integrator.Nodes, ...
                kernel = options.kernel, epsilon = epsilon);

            %< Create nodal approximator
            approximator = mdemat.common.approximation.NodalApproximator(B);
            approximator.setMass(integrator.Nodes, integrator.Weights);
            approximator.setDesign(integrator.Nodes);

            %< Create element
            obj = mdemat.common.fem.L2SphereElement(geometry, approximator, integrator, reduction);
        end
    end
end
