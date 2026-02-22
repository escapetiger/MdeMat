classdef BH1OrthotopeElement < approx.element.BH1Element
    % BH1ORTHOTOPEELEMENT Broken \f$H^1\f$ element over the orthotope.
    %
    %   BH1OrthotopeElement implements a discontinuous finite element over the
    %   orthotope domain. The element supports both modal and nodal basis
    %   functions with various polynomial orders and ordering patterns.
    %
    % See also:
    %   approx.element.BH1Element

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
            geometry = core.geometry.Orthotope(bbox);

            %< Create orthogonal basis functions
            U = core.function.CompiledOrthogonalBasis(order, 'monic_unit_legendre');
            U = repmat(U, 1, nDims);

            %< Create separable basis with specified pattern
            basis = core.function.SeparableFunction(Factors=U, Pattern=options.pattern);

            %< Create volume integrator
            volumeRule = approx.integrate.GaussLegendreRule(nDims);
            volumeIntegrator = approx.integrate.Integrator(volumeRule);
            volumeIntegrator.setPoints(geometry, basis.NFactorCodims);

            %< Create face integrators for all boundaries
            fluxIntegrators = arrayfun(@(i) approx.integrate.Integrator(volumeRule), 1:2*nDims);
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
                faceGeometry = core.geometry.Orthotope(bboxF(:));
                fluxIntegrators(i).setPoints(faceGeometry, basis.NFactorCodims);
            end

            %< Create modal approximator
            approximator = approx.linear.ModalApproximator(basis);
            approximator.setMass(volumeIntegrator.Nodes, volumeIntegrator.Weights);
            approximator.setDesign(volumeIntegrator.Nodes, volumeIntegrator.Weights);

            %< Create element
            obj = approx.element.BH1OrthotopeElement(geometry, approximator, volumeIntegrator, fluxIntegrators);
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
            geometry = core.geometry.Orthotope(bbox);

            %< Create interpolation basis functions
            U = core.function.CompiledInterpolationBasis(order, 'unit_gauss_legendre_lagrange');
            U = repmat(U, 1, nDims);

            %< Create separable basis with Q pattern
            basis = core.function.SeparableFunction(Factors=U, Pattern='Q');

            %< Create volume integrator
            volumeRule = approx.integrate.GaussLegendreRule(nDims);
            volumeIntegrator = approx.integrate.Integrator(volumeRule);
            volumeIntegrator.setPoints(geometry, basis.NFactorCodims);

            %< Create face integrators for all boundaries
            fluxIntegrators = arrayfun(@(i) approx.integrate.Integrator(volumeRule), 1:2*nDims);
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
                faceGeometry = core.geometry.Orthotope(bboxF(:));
                fluxIntegrators(i).setPoints(faceGeometry, basis.NFactorCodims);
            end

            %< Create nodal approximator
            approximator = approx.linear.NodalApproximator(basis);
            approximator.setMass(volumeIntegrator.Nodes, volumeIntegrator.Weights);
            approximator.setDesign(volumeIntegrator.Nodes);

            %< Create element
            obj = approx.element.BH1OrthotopeElement(geometry, approximator, volumeIntegrator, fluxIntegrators);
        end
    end
end