classdef BH1OrthotopeElement < mdemat.common.fem.BH1Element
    % BH1ORTHOTOPEELEMENT Broken \f$H^1\f$ element over the orthotope.
    %
    %   BH1OrthotopeElement implements a discontinuous finite element over the
    %   orthotope domain. The element supports both modal and nodal basis
    %   functions with various polynomial orders and ordering patterns.
    %
    % See also:
    %   mdemat.common.fem.BH1Element

    methods (Static)
        function obj = modal(nDims, order, options)
            % MODAL Create modal BH1OrthotopeElement with orthogonal basis.
            %
            %   obj = modal(nDims, order) creates a BH1OrthotopeElement using
            %   CompiledOrthogonalBasis with 'monic_unit_legendre' and 'Q' pattern.
            %
            %   obj = modal(nDims, order, options) allows customization:
            %   - options.pattern: Basis pattern ('Q' default, 'P' supported)

            arguments
                nDims {mustBePositive, mustBeInteger}
                order {mustBePositive, mustBeInteger}
                options.pattern {mustBeMember(options.pattern, {'P', 'Q'})} = 'Q'
            end

            %< Create unit orthotope geometry
            bbox = repmat([-1/2, 1/2], 1, nDims);
            geometry = mdemat.data.geometry.Orthotope(bbox);

            %< Create orthogonal basis functions
            U = mdemat.common.spectral.CompiledOrthogonalBasis(order, 'monic_unit_legendre');
            U = repmat(U, 1, nDims);

            %< Create separable basis with specified pattern
            basis = mdemat.common.approximation.SeparableFunction(Factors=U, Pattern=options.pattern);

            %< Create volume integrator
            volumeRule = mdemat.common.quadrature.GaussLegendreRule(nDims);
            volumeIntegrator = mdemat.common.quadrature.Integrator(volumeRule);
            volumeIntegrator.setPoints(geometry, basis.NFactorCodims);

            %< Create face integrators for all boundaries
            fluxIntegrators = arrayfun(@(i) mdemat.common.quadrature.Integrator(volumeRule), 1:2*nDims);
            for i = 1:(2 * nDims)
                a = bbox(1:2:end);
                b = bbox(2:2:end);
                d = ceil(i/2);
                if mod(i, 2) == 0
                    a(d) = b(d);
                else
                    b(d) = a(d);
                end
                bboxF = [a(:).'; b(:).'];
                faceGeometry = mdemat.data.geometry.Orthotope(bboxF(:));
                fluxIntegrators(i).setBoundaryPoints(faceGeometry, basis.NFactorCodims);
            end

            %< Create modal approximator
            approximator = mdemat.common.approximation.ModalApproximator(basis);
            approximator.setMass(volumeIntegrator.Nodes, volumeIntegrator.Weights);
            approximator.setDesign(volumeIntegrator.Nodes, volumeIntegrator.Weights);

            %< Create element
            obj = mdemat.common.fem.BH1OrthotopeElement(geometry, approximator, volumeIntegrator, fluxIntegrators);
        end

        function obj = nodal(nDims, order)
            % NODAL Create nodal BH1OrthotopeElement with interpolation basis.
            %
            %   obj = nodal(nDims, order) creates a BH1OrthotopeElement using
            %   CompiledInterpolationBasis with 'unit_gauss_legendre_lagrange'
            %   and 'Q' pattern.

            arguments
                nDims {mustBePositive, mustBeInteger}
                order {mustBePositive, mustBeInteger}
            end

            %< Create unit orthotope geometry
            bbox = repmat([-1/2, 1/2], 1, nDims);
            geometry = mdemat.data.geometry.Orthotope(bbox);

            %< Create interpolation basis functions
            U = mdemat.common.interpolation.CompiledInterpolationBasis(order, 'unit_gauss_legendre_lagrange');
            U = repmat(U, 1, nDims);

            %< Create separable basis with Q pattern
            basis = mdemat.common.approximation.SeparableFunction(Factors=U, Pattern='Q');

            %< Create volume integrator
            volumeRule = mdemat.common.quadrature.GaussLegendreRule(nDims);
            volumeIntegrator = mdemat.common.quadrature.Integrator(volumeRule);
            volumeIntegrator.setPoints(geometry, basis.NFactorCodims);

            %< Create face integrators for all boundaries
            fluxIntegrators = arrayfun(@(i) mdemat.common.quadrature.Integrator(volumeRule), 1:2*nDims);
            for i = 1:(2 * nDims)
                a = bbox(1:2:end);
                b = bbox(2:2:end);
                d = ceil(i/2);
                if mod(i, 2) == 0
                    a(d) = b(d);
                else
                    b(d) = a(d);
                end
                bboxF = [a(:).'; b(:).'];
                faceGeometry = mdemat.data.geometry.Orthotope(bboxF(:));
                fluxIntegrators(i).setBoundaryPoints(faceGeometry, basis.NFactorCodims);
            end

            %< Create nodal approximator
            approximator = mdemat.common.approximation.NodalApproximator(basis);
            approximator.setMass(volumeIntegrator.Nodes, volumeIntegrator.Weights);
            approximator.setDesign(volumeIntegrator.Nodes);

            %< Create element
            obj = mdemat.common.fem.BH1OrthotopeElement(geometry, approximator, volumeIntegrator, fluxIntegrators);
        end
    end
end
