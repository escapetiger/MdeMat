classdef L2SphereElement < approx.element.L2Element
    % L2SPHEREELEMENT L2 element defined on spherical domains.
    %
    %   L2SphereElement extends L2Element for spherical geometries with
    %   various basis function types and reduction schemes. The element
    %   supports 1D (interval), 2D (circle), and 3D (sphere) domains with
    %   orthogonal, cardinal, or radial basis function (RBF) types.
    %
    % See also:
    %   approx.element.L2Element

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
            %   approx.element.L2Element

            obj@approx.element.L2Element(geometry, approximator, integrator);
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
                order{mustBeNonnegative, mustBeInteger}
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
                geometry = core.geometry.Sphere([0], radius);
                R = approx.integrate.GaussTrapezoidalRule(1, coord = 'car');
                E = core.function.Constant(1, A);
                sqrtE = core.function.Constant(1, sqrt(A));
                integrator = approx.integrate.Integrator(R);
                integrator.setPoints(geometry, min(2, np));
                integrator.addWeightFunction(@(v) 1./E.eval(v));

                maxDegree = order - 1;
                F = core.function.LegendreBasis(maxDegree, isNormalized = true);

            elseif nDims == 1 && strcmp(reduction, 'symmetry')
                A = 1 / (2 * radius);
                geometry = core.geometry.Orthotope([-radius, radius]);
                R = approx.integrate.GaussLegendreRule(1);
                E = core.function.Constant(1, A);
                sqrtE = core.function.Constant(1, sqrt(A));
                integrator = approx.integrate.Integrator(R);
                integrator.setPoints(geometry, np);
                integrator.addWeightFunction(@(v) 1./E.eval(v));

                maxDegree = order - 1;
                F = core.function.LegendreBasis(maxDegree, isNormalized = true);

            elseif nDims == 2 && strcmp(reduction, 'topology')
                A = 1 / (2 * pi * radius);
                geometry = core.geometry.Sphere([0, 0], radius);
                R = approx.integrate.GaussTrapezoidalRule(2, coord = 'cossph');
                E = core.function.Constant(1, A);
                sqrtE = core.function.Constant(1, sqrt(A));
                integrator = approx.integrate.Integrator(R);
                integrator.setPoints(geometry, np);
                integrator.addWeightFunction(@(v) 1./E.eval(v));

                maxDegree = order - 1;
                F = core.function.FourierBasis(maxDegree, isNormalized = true);

            else %< nDims >= 3
                A = 1 / (4 * pi * radius^2);
                geometry = core.geometry.Sphere([0, 0, 0], radius);
                R = approx.integrate.GaussTrapezoidalRule(3, coord = 'cossph');
                E = core.function.Constant(2, A);
                sqrtE = core.function.Constant(2, sqrt(A));
                integrator = approx.integrate.Integrator(R);
                integrator.setPoints(geometry, np);
                integrator.addWeightFunction(@(v) 1./E.eval(v));

                maxDegree = order - 1;
                F = core.function.SphericalHarmonicBasis(maxDegree, isNormalized = true);
            end

            %< Create basis function
            B = core.function.Product(F, sqrtE);

            %< Create modal approximator
            approximator = approx.linear.ModalApproximator(B);
            approximator.setMass(integrator.Nodes, integrator.Weights);
            approximator.setDesign(integrator.Nodes, integrator.Weights);

            %< Create element
            obj = approx.element.L2SphereElement(geometry, approximator, integrator, reduction);
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
                order{mustBeNonnegative, mustBeInteger}
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
                geometry = core.geometry.Sphere([0], radius);
                R = approx.integrate.GaussTrapezoidalRule(1, coord = 'car');
                % E = core.function.Constant(1, A);
                integrator = approx.integrate.Integrator(R);
                integrator.setPoints(geometry, min(2, np));
                % integrator.addWeightFunction(@(v) 1./E.eval(v));

            elseif nDims == 1 && strcmp(reduction, 'symmetry')
                % A = 1 / (2 * radius);
                geometry = core.geometry.Orthotope([-radius, radius]);
                R = approx.integrate.GaussLegendreRule(1);
                % E = core.function.Constant(1, A);
                integrator = approx.integrate.Integrator(R);
                integrator.setPoints(geometry, np);
                % integrator.addWeightFunction(@(v) 1./E.eval(v));

            elseif nDims == 2 && strcmp(reduction, 'topology')
                % A = 1 / (2 * pi * radius);
                geometry = core.geometry.Sphere([0, 0], radius);
                R = approx.integrate.GaussTrapezoidalRule(2, coord = 'cossph');
                % E = core.function.Constant(1, A);
                integrator = approx.integrate.Integrator(R);
                integrator.setPoints(geometry, np);
                % integrator.addWeightFunction(@(v) 1./E.eval(v));

            else %< nDims >= 3
                % A = 1 / (4 * pi * radius^2);
                geometry = core.geometry.Sphere([0, 0, 0], radius);
                R = approx.integrate.GaussTrapezoidalRule(3, coord = 'cossph');
                % E = core.function.Constant(1, A);
                integrator = approx.integrate.Integrator(R);
                integrator.setPoints(geometry, np);
                % integrator.addWeightFunction(@(v) 1./E.eval(v));
            end

            %< Create RBF basis function
            if isempty(epsilon)
                h = min(pdist(integrator.Nodes.'));
                epsilon = 1e5 / h;
            else
                epsilon = options.epsilon;
            end
            B = core.function.RbfKernelBasis(integrator.Nodes, ...
                kernel = options.kernel, epsilon = epsilon);

            %< Create nodal approximator
            approximator = approx.linear.NodalApproximator(B);
            approximator.setMass(integrator.Nodes, integrator.Weights);
            approximator.setDesign(integrator.Nodes);

            %< Create element
            obj = approx.element.L2SphereElement(geometry, approximator, integrator, reduction);
        end
    end
end