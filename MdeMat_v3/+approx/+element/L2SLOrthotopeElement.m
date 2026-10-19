classdef L2SLOrthotopeElement < approx.element.L2OrthotopeElement
    % L2SLORTHOTOPEELEMENT Semi-Lagrangian \f$L^2\f$ element over the
    % orthotope.

    properties
        volumePieces % Geometric data for volume pieces (struct)
        %< Fields:
        %<   .n - Number of volume pieces
        %<   .xc - Center coordinates of pieces (nDims x nPieces)
        %<   .hx - Dimensions of pieces (nDims x nPieces)
        %<   .xmin - Minimum coordinates of pieces (nDims x nPieces)
        %<   .xmax - Maximum coordinates of pieces (nDims x nPieces)
        %<   .eulerian - Integration data for eulerian volume pieces (1 x nPieces IntegrationData)
        %<   .upstream - Integration data for upstream volume pieces (1 x nPieces IntegrationData)
        %<   .indexShifts - Index shifts for upstream element mapping (nDims x nPieces)
        fluxPieces % Geometric data for flux pieces (struct array)
        %< Array indexed by boundary face, each struct containing:
        %<   .n - Number of flux pieces
        %<   .xc - Center coordinates of pieces (nDims x nPieces)
        %<   .hx - Dimensions of pieces (nDims x nPieces, zero in normal direction)
        %<   .xmin - Minimum coordinates of pieces (nDims x nPieces)
        %<   .xmax - Maximum coordinates of pieces (nDims x nPieces)
        %<   .eulerian - Integration data for eulerian volume pieces (1 x nPieces IntegrationData)
        %<   .upstream - Integration data for upstream volume pieces (1 x nPieces IntegrationData)
        %<   .indexShifts - Index shifts for upstream element mapping (nDims x nPieces)
    end

    properties (Dependent)
        nVolumePieces % Number of volume pieces (integer)
        nFluxPieces % Number of flux pieces (vector)
    end

    methods
        function obj = L2SLOrthotopeElement(nDims, basisOrder, basisType, basisPattern)
            % L2SLORTHOTOPEELEMENT Create a L2 SL element on the unit
            % orthotope.
            %
            %   obj = L2SLOrthotopeElement(nDims, basisOrder, basisType,
            %   basisPattern) creates a L2 SL element on the unit orthotope
            %   with specified basis configuration.
            %
            % Inputs:
            %   nDims - Number of dimensions (positive integer)
            %   basisOrder - Basis order (positive integer)
            %   basisType - Basis function type (char)
            %   basisPattern - Basis ordering pattern (char, 'P' or 'Q')
            %
            % Outputs:
            %   obj - The L2SLOrthotopeElement object

            obj@approx.element.L2OrthotopeElement(nDims, basisOrder, basisType, basisPattern);
            obj.volumePieces = struct();
            obj.fluxPieces = arrayfun(@(i) struct(), 1:obj.nFluxes);
        end

        function n = get.nVolumePieces(obj)
            % GET.NVOLUMEPIECES Get the number of volume pieces.

            if isfield(obj.volumePieces, 'n')
                n = obj.volumePieces.n;
            else
                n = 0;
            end
        end

        function n = get.nFluxPieces(obj)
            % GET.NFLUXPIECES Get the number of flux pieces for each flux.

            n = zeros(1, obj.nFluxes);
            for i = 1:obj.nFluxes
                if hasfield(obj.fluxPieces(i), 'n')
                    n(i) = obj.fluxPieces(i).n;
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
            %
            % Inputs:
            %   obj - The L2SLOrthotopeElement object
            %   v - Characteristic slope (nDims x 1)
            %   ht - Time step for characteristic tracking (positive scalar)
            %   hx - Mesh spacing vector (nDims x 1 or nDims x 1)
            %
            % Outputs:
            %   obj - The L2SLOrthotopeElement object

            nd = obj.nDims;

            %< Compute clipping geometry
            bbox = obj.geometry.bbox;
            ratio = v(:) .* ht ./ hx(:);
            poi = ratio - floor(ratio) - 1 / 2;

            %< Initialize volume pieces structure
            obj.volumePieces = struct();
            obj.volumePieces.xc = zeros(nd, 2^nd);
            obj.volumePieces.hx = zeros(nd, 2^nd);
            np = 0;

            %< Generate intersection pieces
            indexer = core.linalg.MultiIndexer(ones(1, nd)*2);
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
                obj.volumePieces.xc(:, np) = xc;
                obj.volumePieces.hx(:, np) = hx;
            end

            %< Handle case with no valid pieces
            if np == 0
                obj.volumePieces.xc = zeros(nd, 1);
                obj.volumePieces.hx = ones(nd, 1);
                np = 1;
            else
                obj.volumePieces.xc = obj.volumePieces.xc(:, 1:np);
                obj.volumePieces.hx = obj.volumePieces.hx(:, 1:np);
            end

            %< Compute piece boundaries
            obj.volumePieces.xmin = obj.volumePieces.xc - obj.volumePieces.hx / 2;
            obj.volumePieces.xmax = obj.volumePieces.xc + obj.volumePieces.hx / 2;
            obj.volumePieces.n = np;

            %< Create function data objects for each piece
            I = obj.volume.integrator;
            obj.volumePieces.eulerian = arrayfun(@(i) approx.integrate.IntegrationData(I.copy()), 1:np);
            obj.volumePieces.upstream = arrayfun(@(i) approx.integrate.IntegrationData(I.copy()), 1:np);

            %< Set up integration points for each piece
            nq = obj.approximator.basis.nFactorCodims;
            for i = 1:np
                a = obj.volumePieces.xmin(:, i);
                b = obj.volumePieces.xmax(:, i);
                obj.volumePieces.eulerian(i).integrator.setPoints(nq, a, b);
                obj.volumePieces.upstream(i).integrator.setPoints(nq, -b, -a);
            end

            %< Compute basis function values at integration points
            f = obj.approximator.basis;
            for i = 1:np
                De = obj.volumePieces.eulerian(i);
                Du = obj.volumePieces.upstream(i);
                Ie = De.integrator;
                Iu = Du.integrator;
                Xe = Ie.nodes(1:nd, :);
                Xu = Iu.nodes(1:nd, :);

                %< Set function values and derivatives
                De.setValues(f, Xe);
                Du.setValues(f, Xu);
                De.setDerivatives(f, 1, Xe);
            end

            %< Compute index shifts for upstream element mapping
            obj.volumePieces.indexShifts = zeros(nd, np);
            mm = ceil(ratio) .* (ratio > 0) + (floor(ratio) + 1) .* (ratio < 0);
            mm0 = -mm;
            mm1 = -mm + 1;
            for i = 1:np
                bb = obj.volumePieces.xc(:, i) <= 0;
                obj.volumePieces.indexShifts(:, i) = mm0 .* bb + mm1 .* (1 - bb);
            end
        end

        function obj = clipFlux(obj, v, ht, hx)
            % SETFLUXPIECEDATA Set up flux piece data for semi-Lagrangian
            % transport.
            %
            %   obj = setFluxPieceData(obj, ht, hx) computes the geometric
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
            %
            % Inputs:
            %   obj - The L2SLOrthotopeElement object
            %   v - Characteristic slope (nDims x 1)
            %   ht - Time step for characteristic tracking (positive scalar)
            %   hx - Mesh spacing vector (nDims x 1)
            %
            % Outputs:
            %   obj - The L2SLOrthotopeElement object

            bbox = obj.geometry.bbox;
            bbox = reshape(bbox, 2, []);
            ratio = v(:) .* ht ./ hx(:);
            poi = ratio - floor(ratio) - 1 / 2;
            nd = obj.nDims;
            nf = 2 * nd; %< Number of boundary faces
            np = zeros(1, nf);

            %< Initialize flux pieces structure array
            obj.fluxPieces = arrayfun(@(i) struct(), 1:nf);
            for i = 1:nf
                obj.fluxPieces(i).xc = zeros(nd, 2^(nd - 1));
                obj.fluxPieces(i).hx = zeros(nd, 2^(nd - 1));
            end

            %< Generate intersection pieces for each boundary face
            for i = 1:nf
                d = floor((i - 1)/2) + 1; %< Spatial dimension index
                s = mod(i-1, 2) + 1; %< Side index (1=min, 2=max)

                for j = 1:2^(nd - 1)
                    %< Create multi-indexer for face geometry
                    shape = ones(1, nd) * 2;
                    shape(d) = 1; %< Collapsed dimension for face
                    indexer = core.linalg.MultiIndexer(shape);
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
                    obj.fluxPieces(i).xc(:, np(i)) = xc;
                    obj.fluxPieces(i).hx(:, np(i)) = hx;
                end

                %< Handle case with no valid pieces
                if np(i) == 0
                    obj.fluxPieces(i).xc = zeros(nd, 1);
                    obj.fluxPieces(i).xc(d) = bbox(s, d);
                    obj.fluxPieces(i).hx = ones(nd, 1);
                    obj.fluxPieces(i).hx(d) = 0;
                    np(i) = 1;
                else
                    obj.fluxPieces(i).xc = obj.fluxPieces(i).xc(:, 1:np(i));
                    obj.fluxPieces(i).hx = obj.fluxPieces(i).hx(:, 1:np(i));
                end

                %< Compute piece boundaries
                obj.fluxPieces(i).xmin = obj.fluxPieces(i).xc - obj.fluxPieces(i).hx / 2;
                obj.fluxPieces(i).xmax = obj.fluxPieces(i).xc + obj.fluxPieces(i).hx / 2;

                %< Store number of flux pieces
                obj.fluxPieces(i).n = np(i);
            end

            %< Create function data objects for each boundary and piece
            for i = 1:nf
                I = obj.flux(i).integrator;
                obj.fluxPieces(i).eulerian = arrayfun( ...
                    @(j, k) approx.integrate.IntegrationData(I.copy()), ...
                    repmat((1:2).', 1, np(i)), repmat(1:np(i), 2, 1));
                obj.fluxPieces(i).upstream = arrayfun( ...
                    @(j, k) approx.integrate.IntegrationData(I.copy()), ...
                    repmat((1:2).', 1, np(i)), repmat(1:np(i), 2, 1));
            end

            %< Set up integration points for each boundary and piece
            nq = obj.approximator.basis.nFactorCodims;
            for i = 1:nf
                for j = 1:np(i)
                    a = obj.fluxPieces(i).xmin(:, j);
                    b = obj.fluxPieces(i).xmax(:, j);
                    obj.fluxPieces(i).eulerian(1, j).integrator.setPoints(nq, a, b);
                    obj.fluxPieces(i).eulerian(2, j).integrator.setPoints(nq, a, b);
                    obj.fluxPieces(i).upstream(1, j).integrator.setPoints(nq, -b, -a);
                    obj.fluxPieces(i).upstream(2, j).integrator.setPoints(nq, -b, -a);
                end
            end

            %< Compute basis function values at integration points
            f = obj.approximator.basis;
            for i = 1:nf
                d = floor((i - 1)/2) + 1; %< Spatial dimension index
                s = mod(i-1, 2) + 1; %< Side index

                for j = 1:np(i)
                    %< Eulerian inflow data
                    Dei = obj.fluxPieces(i).eulerian(1, j);
                    Iei = Dei.integrator;
                    Xe = Iei.nodes(1:nd, :);
                    Xe(d, :) = bbox(s, d);
                    Dei.setValues(f, Xe);

                    %< Eulerian outflow data
                    Deo = obj.fluxPieces(i).eulerian(2, j);
                    Ieo = Deo.integrator;
                    Xe = Ieo.nodes(1:nd, :);
                    Xe(d, :) = bbox(3-s, d); %< Opposite side
                    Deo.setValues(f, Xe);

                    %< Upstream inflow data
                    Dui = obj.fluxPieces(i).upstream(1, j);
                    Iui = Dui.integrator;
                    Xu = Iui.nodes(1:nd, :);
                    Xu(d, :) = -poi(d);
                    Dui.setValues(f, Xu);

                    %< Upstream outflow data
                    Duo = obj.fluxPieces(i).upstream(2, j);
                    Iuo = Duo.integrator;
                    Xu = Iuo.nodes(1:nd, :);
                    if -poi(d) == bbox(1, d)
                        Xu(d, :) = bbox(2, d);
                    elseif -poi(d) == bbox(2, d)
                        Xu(d, :) = bbox(1, d);
                    end
                    Duo.setValues(f, Xu);
                end
            end

            %< Compute index shifts for upstream element mapping
            for i = 1:nf
                obj.fluxPieces(i).indexShifts = zeros(nd, np(i));
                mm = ceil(ratio) .* (ratio > 0) + (floor(ratio) + 1) .* (ratio < 0);
                mm0 = -mm;
                mm1 = -mm + 1;
                for j = 1:np(i)
                    bb = obj.fluxPieces(i).xc(:, j) <= 0;
                    obj.fluxPieces(i).indexShifts(:, j) = mm0 .* bb + mm1 .* (1 - bb);
                end
            end
        end
    end
end
