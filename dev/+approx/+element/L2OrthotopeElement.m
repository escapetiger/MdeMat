classdef L2OrthotopeElement < approx.element.L2Element
    % L2ORTHOTOPEELEMENT L2 element defined on orthotope domains.
    %
    %   L2OrthotopeElement extends L2Element for orthotope geometries with
    %   various basis function types and reduction schemes. The element
    %   supports 1D (interval), 2D (rectangle), and 3D (cuboid) domains
    %   with orthogonal basis functions.
    %
    % See also:
    %   approx.element.L2Element

    methods (Static)
        function obj = hermite(nDims, order, options)
            % HERMITE Create L2OrthotopeElement with Hermite basis.
            %
            %   obj = hermite(nDims, order) creates an L2OrthotopeElement using
            %   HermiteBasis with normalized Hermite functions.
            %
            %   obj = hermite(nDims, order, options) allows customization:
            %   - options.bbox: Bounding box (default: [-inf, inf] for each dimension)
            %   - options.T: Temperature parameter (default: 1.0)

            arguments
                nDims {mustBePositive, mustBeInteger}
                order {mustBeNonnegative, mustBeInteger}
                options.T {mustBeNumeric} = 1.0
            end

            %< Create orthotope geometry
            geometry = core.geometry.Orthotope(repmat([-inf, inf], 1, nDims));

            %< Create Hermite basis functions
            U = core.function.HermiteBasis(order - 1, IsNormalized=true);
            L = core.function.Linear(diag(1./sqrt(options.T)));
            U = core.function.Composition(U, L);
            U = repmat(U, 1, nDims);

            %< Create separable basis with Q pattern
            E = core.function.SeparableFunction(factors=U, pattern='Q');
            M = core.function.Maxwellian(nDims, temperature=options.T);
            B = core.function.Product(E, M);

            %< Create Gauss-Hermite integrator
            R = approx.integrate.GaussHermiteRule(nDims);
            integrator = approx.integrate.Integrator(R);
            integrator.setPoints(geometry, E.NFactorCodims+1);
            integrator.addWeightFunction(@(v) 1./ M.eval(v).^2);

            %< Create modal approximator
            approximator = approx.linear.ModalApproximator(B);
            approximator.Mass = eye(approximator.NDofs);
            approximator.Design = eye(approximator.NDofs);

            %< Create element
            obj = approx.element.L2OrthotopeElement(geometry, approximator, integrator);
        end
    end
end