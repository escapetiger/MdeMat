classdef BH1SLOrthotopeElement < approx.element.BH1OrthotopeElement
    % BH1SLORTHOTOPEELEMENT Semi-Lagrangian broken \f$H^1\f$ element over the orthotope.
    %
    %   BH1SLOrthotopeElement extends BH1OrthotopeElement with semi-Lagrangian
    %   transport capabilities for advection-dominated problems. The element
    %   supports volume and flux piece clipping operations for upstream
    %   integration.

    properties
        VolumePieces % Geometric data for volume pieces (struct)
        FluxPieces % Geometric data for flux pieces (struct array)
    end

    properties (Dependent)
        NVolumePieces % Number of volume pieces (integer)
        NFluxPieces % Number of flux pieces (vector)
    end

    methods
        function obj = BH1SLOrthotopeElement(geometry, approximator, volumeIntegrator, fluxIntegrators)
            % BH1SLORTHOTOPEELEMENT Create a broken H1 semi-Lagrangian element on
            % the orthotope.
            %
            %   obj = BH1SLOrthotopeElement(geometry, approximator, volumeIntegrator,
            %   fluxIntegrators) creates a broken H1 semi-Lagrangian element with specified
            %   geometry, approximator, and integrators.

            arguments
                geometry core.geometry.Orthotope
                approximator approx.linear.LinearApproximator
                volumeIntegrator approx.integrate.Integrator
                fluxIntegrators (1, :) {mustBeA(fluxIntegrators, 'approx.integrate.Integrator')}
            end

            obj@approx.element.BH1OrthotopeElement(geometry, approximator, volumeIntegrator, fluxIntegrators);
            obj.VolumePieces = struct();
            obj.FluxPieces = arrayfun(@(i) struct(), 1:obj.NFluxes);
        end

        function n = get.NVolumePieces(obj)
            % GET.NVOLUMEPIECES Get the number of volume pieces.

            if isfield(obj.VolumePieces, 'n')
                n = obj.VolumePieces.n;
            else
                n = 0;
            end
        end

        function n = get.NFluxPieces(obj)
            % GET.NFLUXPIECES Get the number of flux pieces for each flux.

            n = zeros(1, obj.NFluxes);
            for i = 1:obj.NFluxes
                if isfield(obj.FluxPieces(i), 'n')
                    n(i) = obj.FluxPieces(i).n;
                end
            end
        end

        function obj = clipVolume(obj, v, ht, hx)
            % CLIPVOLUME Compute volume piece data for semi-Lagrangian
            % transport.
            %
            %   obj = clipVolume(obj, v, ht, hx) computes the geometric
            %   intersection between the current element and the upstream
            %   footprint, creating volume pieces with associated function
            %   data for accurate semi-Lagrangian integration.
            %
            %   The method performs several key operations:
            %   1. Computes geometric intersection pieces
            %   2. Creates function data objects for each piece
            %   3. Sets up quadrature points and weights
            %   4. Evaluates basis functions at integration points
            %   5. Computes index shifts for upstream element mapping

            arguments
                obj approx.element.BH1SLOrthotopeElement
                v(:, 1) {mustBeReal, mustBeNonempty}
                ht(1, 1) {mustBeReal, mustBeNonempty}
                hx(:, 1) {mustBeReal, mustBeNonempty}
            end

            nd = obj.NDims;

            %< Compute clipping geometry
            bbox = obj.Geometry.BBox;
            ratio = v(:) .* ht ./ hx(:);
            poi = ratio - floor(ratio) - 1 / 2;

            %< Initialize volume pieces structure
            obj.VolumePieces = struct();
            obj.VolumePieces.xc = zeros(nd, 2^nd);
            obj.VolumePieces.hx = zeros(nd, 2^nd);
            np = 0;

            %< Generate intersection pieces
            indexer = core.linalg.MultiIndexer(shape = repmat(2, 1, nd));
            m = indexer.generate();
            for i = 1:2^nd
                x1 = poi;
                x2 = bbox((m(i, :) + (0:nd - 1) * 2)');
                xc = (x1 + x2) / 2;
                hx = abs(x2-x1);

                %< Skip degenerate pieces
                if prod(hx) < 1e-8
                    continue;
                end

                np = np + 1;
                obj.VolumePieces.xc(:, np) = xc;
                obj.VolumePieces.hx(:, np) = hx;
            end

            %< Handle case with no valid pieces
            if np == 0
                obj.VolumePieces.xc = zeros(nd, 1);
                obj.VolumePieces.hx = ones(nd, 1);
                np = 1;
            else
                obj.VolumePieces.xc = obj.VolumePieces.xc(:, 1:np);
                obj.VolumePieces.hx = obj.VolumePieces.hx(:, 1:np);
            end

            %< Compute piece boundaries
            obj.VolumePieces.xmin = obj.VolumePieces.xc - obj.VolumePieces.hx / 2;
            obj.VolumePieces.xmax = obj.VolumePieces.xc + obj.VolumePieces.hx / 2;
            obj.VolumePieces.n = np;

            %< Create function data objects for each piece
            obj.VolumePieces.eulerian = arrayfun(@(i) obj.Volume.copy(), 1:np);
            obj.VolumePieces.upstream = arrayfun(@(i) obj.Volume.copy(), 1:np);

            %< Set up integration points for each piece
            nq = obj.Approximator.Basis.NFactorCodims;
            for i = 1:np
                a = obj.VolumePieces.xmin(:, i);
                b = obj.VolumePieces.xmax(:, i);
                bbox = [a(:).'; b(:).'];
                G = core.geometry.Orthotope(bbox(:));
                obj.VolumePieces.eulerian(i).setPoints(G, nq);
                bbox = [-b(:).'; -a(:).'];
                G = core.geometry.Orthotope(bbox(:));
                obj.VolumePieces.upstream(i).setPoints(G, nq);
            end

            %< Compute basis function values at integration points
            f = obj.Approximator.Basis;
            for i = 1:np
                Ie = obj.VolumePieces.eulerian(i);
                Iu = obj.VolumePieces.upstream(i);
                Xe = Ie.Nodes(1:nd, :);
                Xu = Iu.Nodes(1:nd, :);

                %< Set function values and derivatives
                Ie.setValues(f, nodes=Xe);
                Iu.setValues(f, nodes=Xu);
                Ie.setPartials(f, 1, nodes=Xe);
            end

            %< Compute index shifts for upstream element mapping
            obj.VolumePieces.indexShifts = zeros(nd, np);
            mm = ceil(ratio) .* (ratio > 0) + (floor(ratio) + 1) .* (ratio < 0);
            mm0 = -mm;
            mm1 = -mm + 1;
            for i = 1:np
                bb = obj.VolumePieces.xc(:, i) <= 0;
                obj.VolumePieces.indexShifts(:, i) = mm0 .* bb + mm1 .* (1 - bb);
            end
        end

        function obj = clipFlux(obj, v, ht, hx)
            % CLIPFLUX Set up flux piece data for semi-Lagrangian
            % transport.
            %
            %   obj = clipFlux(obj, v, ht, hx) computes the geometric
            %   intersections between element boundaries and upstream
            %   footprints, creating flux pieces with associated function
            %   data for accurate semi-Lagrangian discontinuous Galerkin
            %   integration. Each boundary face is processed to handle
            %   inflow and outflow contributions.
            %
            %   The method performs several operations for each boundary:
            %   1. Computes geometric intersection pieces on boundary faces
            %   2. Creates function data objects for inflow/outflow directions
            %   3. Sets up quadrature points and weights on boundary pieces
            %   4. Evaluates basis functions at boundary integration points
            %   5. Computes index shifts for upstream element mapping

            arguments
                obj approx.element.BH1SLOrthotopeElement
                v(:, 1) {mustBeReal, mustBeNonempty}
                ht(1, 1) {mustBeReal, mustBeNonempty}
                hx(:, 1) {mustBeReal, mustBeNonempty}
            end

            bbox = obj.Geometry.BBox;
            bbox = reshape(bbox, 2, []);
            ratio = v(:) .* ht ./ hx(:);
            poi = ratio - floor(ratio) - 1 / 2;
            nd = obj.NDims;
            nf = 2 * nd; %< Number of boundary faces
            np = zeros(1, nf);

            %< Initialize flux pieces structure array
            obj.FluxPieces = arrayfun(@(i) struct(), 1:nf);
            for i = 1:nf
                obj.FluxPieces(i).xc = zeros(nd, 2^(nd - 1));
                obj.FluxPieces(i).hx = zeros(nd, 2^(nd - 1));
            end

            %< Generate intersection pieces for each boundary face
            for i = 1:nf
                d = floor((i - 1)/2) + 1; %< Spatial dimension index
                s = mod(i-1, 2) + 1; %< Side index (1=min, 2=max)

                for j = 1:2^(nd - 1)
                    %< Create multi-indexer for face geometry
                    shape = ones(1, nd) * 2;
                    shape(d) = 1; %< Collapsed dimension for face
                    indexer = core.linalg.MultiIndexer(shape = shape);
                    m = indexer.linearToMulti(j);

                    %< Compute piece geometry
                    x1 = poi;
                    x2 = bbox((m + (0:nd - 1) * 2)');
                    xc = (x1 + x2) / 2;
                    hx = abs(x2-x1);

                    %< Skip degenerate pieces (check reduced dimensions)
                    if prod(hx([1:d - 1, d + 1:end])) < 1e-8
                        continue;
                    end

                    %< Set face-specific geometry
                    hx(d) = 0; %< Zero thickness in normal direction
                    xc(d) = bbox(s, d); %< Position on boundary face

                    np(i) = np(i) + 1;
                    obj.FluxPieces(i).xc(:, np(i)) = xc;
                    obj.FluxPieces(i).hx(:, np(i)) = hx;
                end

                %< Handle case with no valid pieces
                if np(i) == 0
                    obj.FluxPieces(i).xc = zeros(nd, 1);
                    obj.FluxPieces(i).xc(d) = bbox(s, d);
                    obj.FluxPieces(i).hx = ones(nd, 1);
                    obj.FluxPieces(i).hx(d) = 0;
                    np(i) = 1;
                else
                    obj.FluxPieces(i).xc = obj.FluxPieces(i).xc(:, 1:np(i));
                    obj.FluxPieces(i).hx = obj.FluxPieces(i).hx(:, 1:np(i));
                end

                %< Compute piece boundaries
                obj.FluxPieces(i).xmin = obj.FluxPieces(i).xc - obj.FluxPieces(i).hx / 2;
                obj.FluxPieces(i).xmax = obj.FluxPieces(i).xc + obj.FluxPieces(i).hx / 2;

                %< Store number of flux pieces
                obj.FluxPieces(i).n = np(i);
            end

            %< Create function data objects for each boundary and piece
            for i = 1:nf
                obj.FluxPieces(i).eulerian = arrayfun( ...
                    @(j, k) obj.Flux(i).copy(), ...
                    repmat((1:2).', 1, np(i)), repmat(1:np(i), 2, 1));
                obj.FluxPieces(i).upstream = arrayfun( ...
                    @(j, k) obj.Flux(i).copy(), ...
                    repmat((1:2).', 1, np(i)), repmat(1:np(i), 2, 1));
            end

            %< Set up integration points for each boundary and piece
            nq = obj.Approximator.Basis.NFactorCodims;
            for i = 1:nf
                for j = 1:np(i)
                    a = obj.FluxPieces(i).xmin(:, j);
                    b = obj.FluxPieces(i).xmax(:, j);
                    bbox = [a(:).'; b(:).'];
                    G = core.geometry.Orthotope(bbox(:));
                    obj.FluxPieces(i).eulerian(1, j).setPoints(G, nq);
                    obj.FluxPieces(i).eulerian(2, j).setPoints(G, nq);
                    bbox = [-b(:).'; -a(:).'];
                    G = core.geometry.Orthotope(bbox(:));
                    obj.FluxPieces(i).upstream(1, j).setPoints(G, nq);
                    obj.FluxPieces(i).upstream(2, j).setPoints(G, nq);
                end
            end

            %< Compute basis function values at integration points
            f = obj.Approximator.Basis;
            for i = 1:nf
                d = floor((i - 1)/2) + 1; %< Spatial dimension index
                s = mod(i-1, 2) + 1; %< Side index

                for j = 1:np(i)
                    %< Eulerian inflow data
                    Iei = obj.FluxPieces(i).eulerian(1, j);
                    Xe = Iei.Nodes(1:nd, :);
                    Xe(d, :) = bbox(s, d);
                    Iei.setValues(f, Xe);

                    %< Eulerian outflow data
                    Ieo = obj.FluxPieces(i).eulerian(2, j);
                    Xe = Ieo.Nodes(1:nd, :);
                    Xe(d, :) = bbox(3-s, d); %< Opposite side
                    Ieo.setValues(f, Xe);

                    %< Upstream inflow data
                    Iui = obj.FluxPieces(i).upstream(1, j);
                    Xu = Iui.Nodes(1:nd, :);
                    Xu(d, :) = -poi(d);
                    Iui.setValues(f, Xu);

                    %< Upstream outflow data
                    Iuo = obj.FluxPieces(i).upstream(2, j);
                    Xu = Iuo.Nodes(1:nd, :);
                    if -poi(d) == bbox(1, d)
                        Xu(d, :) = bbox(2, d);
                    elseif -poi(d) == bbox(2, d)
                        Xu(d, :) = bbox(1, d);
                    end
                    Iuo.setValues(f, Xu);
                end
            end

            %< Compute index shifts for upstream element mapping
            for i = 1:nf
                obj.FluxPieces(i).indexShifts = zeros(nd, np(i));
                mm = ceil(ratio) .* (ratio > 0) + (floor(ratio) + 1) .* (ratio < 0);
                mm0 = -mm;
                mm1 = -mm + 1;
                for j = 1:np(i)
                    bb = obj.FluxPieces(i).xc(:, j) <= 0;
                    obj.FluxPieces(i).indexShifts(:, j) = mm0 .* bb + mm1 .* (1 - bb);
                end
            end
        end
    end

    methods (Static)
        function obj = modal(nDims, order, options)
            % MODAL Create modal BH1SLOrthotopeElement with orthogonal basis.
            %
            %   obj = modal(nDims, order) creates a BH1SLOrthotopeElement using
            %   CompiledOrthogonalBasis with 'monic_unit_legendre' and 'Q' pattern.
            %
            %   obj = modal(nDims, order, options) allows customization:
            %   - options.pattern: Basis pattern ('Q' default, 'P' supported)

            arguments
                nDims {mustBePositive, mustBeInteger}
                order {mustBePositive, mustBeInteger}
                options.pattern {mustBeMember(options.pattern, {'P', 'Q'})} = 'Q'
            end

            % Create unit orthotope geometry
            bbox = repmat([-1/2, 1/2], 1, nDims);
            geometry = core.geometry.Orthotope(bbox);

            % Create orthogonal basis functions
            U = core.function.CompiledOrthogonalBasis(order, 'monic_unit_legendre');
            U = repmat(U, 1, nDims);

            % Create separable basis with specified pattern
            basis = core.function.SeparableFunction(Factors=U, Pattern=options.pattern);

            % Create volume integrator
            volumeRule = approx.integrate.GaussLegendreRule(nDims);
            volumeIntegrator = approx.integrate.Integrator(volumeRule);
            volumeIntegrator.setPoints(geometry, basis.NFactorCodims);

            % Create face integrators for all boundaries
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

            % Create modal approximator
            approximator = approx.linear.ModalApproximator(basis);
            approximator.setMass(volumeIntegrator.Nodes, volumeIntegrator.Weights);
            approximator.setDesign(volumeIntegrator.Nodes, volumeIntegrator.Weights);

            % Create element
            obj = approx.element.BH1SLOrthotopeElement(geometry, approximator, volumeIntegrator, fluxIntegrators);
        end

        function obj = nodal(nDims, order)
            % NODAL Create nodal BH1SLOrthotopeElement with interpolation basis.
            %
            %   obj = nodal(nDims, order) creates a BH1SLOrthotopeElement using
            %   CompiledInterpolationBasis with 'unit_gauss_legendre_lagrange'
            %   and 'Q' pattern.

            arguments
                nDims {mustBePositive, mustBeInteger}
                order {mustBePositive, mustBeInteger}
            end

            % Create unit orthotope geometry
            bbox = repmat([-1/2, 1/2], 1, nDims);
            geometry = core.geometry.Orthotope(bbox);

            % Create interpolation basis functions
            U = core.function.CompiledInterpolationBasis(order, 'unit_gauss_legendre_lagrange');
            U = repmat(U, 1, nDims);

            % Create separable basis with Q pattern
            basis = core.function.SeparableFunction(Factors=U, Pattern='Q');

            % Create volume integrator
            volumeRule = approx.integrate.GaussLegendreRule(nDims);
            volumeIntegrator = approx.integrate.Integrator(volumeRule);
            volumeIntegrator.setPoints(geometry, basis.NFactorCodims);

            % Create face integrators for all boundaries
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

            % Create nodal approximator
            approximator = approx.linear.NodalApproximator(basis);
            approximator.setMass(volumeIntegrator.Nodes, volumeIntegrator.Weights);
            approximator.setDesign(volumeIntegrator.Nodes);

            % Create element
            obj = approx.element.BH1SLOrthotopeElement(geometry, approximator, volumeIntegrator, fluxIntegrators);
        end
    end
end