classdef AffineDiscretization < handle
    % AFFINEDISCRETIZATION Affine discretization for function
    % approximation.
    %
    %   AffineDiscretization provides a framework for constructing affine
    %   function spaces on spherical and orthotopal domains. It manages
    %   element construction with various basis functions including
    %   Legendre, Fourier, and spherical harmonics, along with appropriate
    %   integration schemes.
    %
    %   The class supports C0 sphere elements with different dimensional
    %   reductions and scaling parameters. It encapsulates the complex
    %   process of coordinating geometry, basis functions, projectors,
    %   and integrators for spherical approximation problems.
    %
    % Examples:
    %   % Create 1D affine discretization with Legendre basis
    %   disc = AffineDiscretization();
    %   disc.setC0SphereElement('nDims', 1, 'order', 3, 'nPoints', 10);
    %   disc.setAffineSpace();
    %   
    %   % Create 2D spherical discretization with Fourier basis
    %   disc = AffineDiscretization();
    %   disc.setC0SphereElement('nDims', 2, 'reduction', 1, 'radius', 2.0);
    %   disc.setScaledAffineSpace([1.0, 0.5]);
    %   
    %   % Set up derivative computation
    %   disc.setDerivativeOrder(2);
    %
    % See also:
    %   approx.element.C0Element, approx.space.AffineSpace,
    %   approx.space.ScaledAffineSpace

    properties
        element % Element object for spatial discretization
        space   % Function space object for approximation
    end

    methods
        function obj = AffineDiscretization()
            % AFFINEDISCRETIZATION Constructor for AffineDiscretization.
            %
            %   obj = AffineDiscretization() creates an empty affine
            %   discretization object. Element and space must be set
            %   using the appropriate setter methods.
            %
            % Outputs:
            %   obj - Constructed AffineDiscretization object

            obj.element = [];
            obj.space = [];
        end

        function obj = setC0SphereElement(obj, varargin)
            % SETC0SPHEREELEMENT Set C0 sphere element with specified
            % configuration.
            %
            %   obj = setC0SphereElement(obj, varargin) creates a C0 sphere
            %   element with the specified parameters and assigns it to the
            %   discretization. The method configures geometry, basis
            %   functions, integrators, and projectors based on dimensional
            %   and reduction parameters.
            %
            % Inputs:
            %   obj - The AffineDiscretization object
            %   varargin - Name-value parameter pairs:
            %<    'nDims' - Number of spatial dimensions (1, 2, or 3)
            %<    'reduction' - Reduction parameter (1 or 2, default: 1)
            %<    'order' - Polynomial order (positive integer, default: 1)
            %<    'radius' - Sphere radius (positive scalar, default: 1)
            %<    'nPoints' - Number of integration points (positive integer)
            %
            % Outputs:
            %   obj - The AffineDiscretization object
            %
            % Notes:
            %   Reduction parameter affects basis choice:
            %     - nDims=1, reduction=1: Legendre basis on hypersphere
            %     - nDims=1, reduction=2: Legendre basis on orthotope
            %     - nDims=2, reduction=1: Fourier basis on circle
            %     - nDims=3: Spherical harmonic basis on sphere

            p = inputParser;
            addParameter(p, 'nDims', [], @(x) x >= 1 && x <= 3);
            addParameter(p, 'reduction', 1, @(x) x >= 1 && x <= 2);
            addParameter(p, 'order', 1, @(x) x >= 1);
            addParameter(p, 'radius', 1, @(x) x > 0);
            addParameter(p, 'nPoints', [], @(x) isnumeric(x));
            parse(p, varargin{:});

            d = p.Results.nDims;
            n = p.Results.nPoints;
            k = p.Results.order;
            t = p.Results.reduction;
            r = p.Results.radius;

            if d == 1 && t == 1
                G = core.geometry.Hypersphere(0);
                F = approx.basis.OrthogonalBasisFunction( ...
                    k, 'legendre', 'canonical');
                F.autoLoad();
                E = core.function.ConstantFunction(1, 1/2);
                I = approx.integrate.HypersphereIntegrator( ...
                    1, 'gauss_trapezoidal');
                I.setPoints(n);
                I.addWeightFunction(@(v) 1./E.evaluate(v));
            elseif d == 1 && t == 2
                G = core.geometry.Orthotope([-1, 1]);
                F = approx.basis.OrthogonalBasisFunction( ...
                    k, 'legendre', 'canonical');
                F.autoLoad();
                E = core.function.ConstantFunction(1, 1/(2 * r));
                I = approx.integrate.OrthotopeIntegrator( ...
                    1, 'gauss_legendre');
                I.setPoints(n, -r, r);
                I.addWeightFunction(@(v) 1./E.evaluate(v));
            elseif d == 2 && t == 1
                G = core.geometry.Hypersphere(1);
                F = approx.basis.OrthogonalBasisFunction( ...
                    k, 'fourier', 'canonical');
                F.autoLoad();
                E = core.function.ConstantFunction(1, 1/(2 * pi * r));
                I = approx.integrate.HypersphereIntegrator( ...
                    2, 'gauss_trapezoidal');
                I.setPoints(n, zeros(1, 2), r, 3);
                I.addWeightFunction(@(v) 1./E.evaluate(v));
            else
                G = core.geometry.Hypersphere(2);
                F = approx.basis.OrthogonalBasisFunction( ...
                    k, 'spherical_harmonic', 'canonical');
                F.autoLoad();
                E = core.function.ConstantFunction(2, 1/(4 * pi * r^2));
                I = approx.integrate.HypersphereIntegrator( ...
                    3, 'gauss_trapezoidal');
                I.setPoints(n, zeros(1, 3), r, 3);
                I.addWeightFunction(@(v) 1./E.evaluate(v));
            end
            B = core.function.ProductFunction(F, E);
            P = approx.project.ModalProjector(B);
            P.setMass(I.nodes, I.weights);
            obj.element = approx.element.C0Element(G, P, I);
        end

        function obj = setAffineSpace(obj)
            % SETAFFINESPACE Set affine function space.
            %
            %   obj = setAffineSpace(obj, varargin) creates an affine
            %   function space based on the current element. The space
            %   provides the computational framework for function
            %   approximation on the configured domain.
            %
            % Inputs:
            %   obj - The AffineDiscretization object
            %
            % Outputs:
            %   obj - The AffineDiscretization object

            core.except.assert(~isempty(obj.element), ...
                'NotPrepared', 'Element must be set before space.');

            obj.space = approx.space.AffineSpace(obj.element);
        end

        function obj = setScaledAffineSpace(obj, scales)
            % SETSCALEDAFFINESPACE Set scaled affine function space.
            %
            %   obj = setScaledAffineSpace(obj, scales) creates a scaled
            %   affine function space with the specified scaling
            %   parameters. This allows for anisotropic scaling of the
            %   basis functions in different coordinate directions.
            %
            % Inputs:
            %   obj - The AffineDiscretization object
            %   scales - Scaling parameter vector (positive values)
            %
            % Outputs:
            %   obj - The AffineDiscretization object

            core.except.assert(~isempty(obj.element), ...
                'NotPrepared', 'Element must be set before space.');

            obj.space = approx.space.ScaledAffineSpace(obj.element, scales);
        end

        function obj = setDerivativeOrder(obj, order)
            % SETDERIVATIVEORDER Set derivative order for element data.
            %
            %   obj = setDerivativeOrder(obj, order) configures the element
            %   to compute function data up to the specified derivative
            %   order. This is required for problems involving derivatives.
            %
            % Inputs:
            %   obj - The AffineDiscretization object
            %   order - Maximum derivative order (non-negative integer)
            %
            % Outputs:
            %   obj - The AffineDiscretization object

            core.except.assert(~isempty(obj.element), ...
                'NotPrepared', 'Element must be set before setting derivative order.');
            obj.element.setVolumeData(order);
        end
    end
end