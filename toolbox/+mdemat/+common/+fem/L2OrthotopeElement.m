classdef L2OrthotopeElement < mdemat.common.fem.L2Element
    % L2ORTHOTOPEELEMENT L2 element defined on orthotope domains.
    %
    %   L2OrthotopeElement extends L2Element for orthotope geometries with
    %   various basis function types and reduction schemes. The element
    %   supports 1D (interval), 2D (rectangle), and 3D (cuboid) domains
    %   with orthogonal basis functions.
    %
    % See also:
    %   mdemat.common.fem.L2Element

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
                order {mustBePositive, mustBeInteger}
                options.T (1, 1) {mustBeNumeric, mustBeReal, mustBeFinite, mustBePositive} = 1.0
            end

            %< Create orthotope geometry
            geometry = mdemat.data.geometry.Orthotope(repmat([-inf, inf], 1, nDims));

            %< Create Hermite basis functions
            U = mdemat.common.spectral.HermiteBasis(order - 1, IsNormalized=true);
            L = mdemat.data.function.Linear(diag(1./sqrt(options.T)));
            U = mdemat.data.function.Composition(U, L);
            U = repmat(U, 1, nDims);

            %< Create separable basis with Q pattern
            E = mdemat.common.approximation.SeparableFunction(factors=U, pattern='Q');
            M = mdemat.data.function.Maxwellian(nDims, temperature=options.T);
            B = mdemat.data.function.Product(E, M);

            %< Create Gauss-Hermite integrator
            R = mdemat.common.quadrature.GaussHermiteRule(nDims);
            integrator = mdemat.common.quadrature.Integrator(R);
            integrator.setPoints(geometry, E.NFactorCodims+1);
            % The reference Hermite rule integrates z~N(0,I). Physical
            % velocities v=sqrt(T)*z have density M_T; retain probability
            % weights while scaling nodes before applying 1/M_T^2.
            integrator.Nodes = sqrt(options.T)*integrator.Nodes;
            integrator.addWeightFunction(@(v) 1./ M.eval(v).^2);

            %< Create modal approximator
            approximator = mdemat.common.approximation.ModalApproximator(B);
            approximator.Mass = eye(approximator.NDofs);
            approximator.Design = eye(approximator.NDofs);

            %< Create element
            obj = mdemat.common.fem.L2OrthotopeElement(geometry, approximator, integrator);
        end
    end
end
