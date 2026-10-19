classdef C0ElementClipper < approx.element.ElementClipper
    % C0ELEMENTCLIPPER Geometric clipper for C0 elements.
    %
    %   C0ElementClipper implements geometric clipping operations for C0
    %   elements in semi-Lagrangian methods. It computes the intersection
    %   between the current Eulerian element and the upstream footprint
    %   traced along characteristic curves, creating volume pieces for
    %   accurate transport calculations.
    %
    %   The clipper handles the complex geometry resulting from backward
    %   characteristic tracing by decomposing the intersection region
    %   into simpler geometric pieces. Each piece is equipped with
    %   appropriate quadrature rules and function data for both Eulerian
    %   and upstream (Lagrangian) representations.
    %
    %   This approach enables stable and accurate semi-Lagrangian schemes
    %   for convection-dominated problems without the CFL restriction
    %   that limits explicit Eulerian methods. The geometric decomposition
    %   ensures accurate integration over potentially complex intersection
    %   regions while maintaining computational efficiency.
    %
    % Examples:
    %   % Create clipper for C0 element with specified velocity
    %   velocity = [0.5; 0.3];
    %   element = approx.element.C0Element(geometry, projector, integrator);
    %   clipper = C0ElementClipper(element, velocity);
    %   
    %   % Set up volume piece data for time step and mesh spacing
    %   timeStep = 0.1;
    %   meshSpacing = [0.05; 0.05];
    %   clipper.setVolumePieceData(timeStep, meshSpacing);
    %   
    %   % Access computed pieces and data
    %   nPieces = clipper.nVolumePieces;
    %   eulerianData = clipper.volumeEulerianData;
    %   upstreamData = clipper.volumeUpstreamData;
    %
    %   % Process geometric information
    %   pieceCenters = clipper.volumePieces.xc;
    %   pieceDimensions = clipper.volumePieces.hx;
    %   indexShifts = clipper.volumeShifts;
    %
    % See also:
    %   approx.element.ElementClipper, approx.element.L2ElementClipper,
    %   approx.element.C0Element

    properties
        nVolumePieces % Number of volume pieces from intersection computation (integer)
        volumeShifts % Index shifts for upstream element mapping (matrix)
        volumePieces % Geometric data for volume pieces (struct)
                     % Fields:
                     %   .xc - Center coordinates of pieces (nDims x nPieces)
                     %   .hx - Dimensions of pieces (nDims x nPieces)
                     %   .xmin - Minimum coordinates of pieces (nDims x nPieces)
                     %   .xmax - Maximum coordinates of pieces (nDims x nPieces)
        volumeEulerianData % Function data for Eulerian volume pieces (ElementFunction array)
        volumeUpstreamData % Function data for upstream volume pieces (ElementFunction array)
    end

    methods
        function obj = setVolumePieceData(obj, ht, hx)
            % SETVOLUMEPIECEDATA Set up volume piece data for
            % semi-Lagrangian transport.
            %
            %   obj = setVolumePieceData(obj, ht, hx) computes the
            %   geometric intersection between the current element and the
            %   upstream footprint, creating volume pieces with associated
            %   function data for accurate semi-Lagrangian integration.
            %
            %   The method performs several key operations:
            %   1. Computes geometric intersection pieces
            %   2. Creates function data objects for each piece
            %   3. Sets up quadrature points and weights
            %   4. Evaluates basis functions at integration points
            %   5. Computes index shifts for upstream element mapping
            %
            % Inputs:
            %   obj - The C0ElementClipper object
            %   ht - Time step for characteristic tracking (positive scalar)
            %   hx - Mesh spacing vector (nDims x 1 or nDims x 1)
            %
            % Outputs:
            %   obj - The C0ElementClipper object

            nd = obj.element.nDims;
            v = obj.velocity;

            %< Compute clipping geometry
            bbox = obj.element.geometry.bbox;
            ratio = v(:) .* ht ./ hx(:);
            poi = ratio - floor(ratio) - 1 / 2;
            
            %< Initialize volume pieces structure
            obj.volumePieces = struct();
            obj.volumePieces.xc = zeros(nd, 2^nd);
            obj.volumePieces.hx = zeros(nd, 2^nd);
            np = 0;
            
            %< Generate intersection pieces
            indexer = core.linalg.MultiIndexer(ones(1, nd) * 2);
            m = indexer.generate();
            for i = 1:2^nd
                x1 = poi;
                x2 = bbox((m(i, :) + (0:nd - 1) * 2)');
                xc = (x1 + x2) / 2;
                hx = abs(x2 - x1);
                
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
            obj.nVolumePieces = np;

            %< Create function data objects for each piece
            I = obj.element.volumeData.integrator;
            obj.volumeEulerianData = arrayfun(@(i) approx.element.ElementFunction(I.copy()), 1:np);
            obj.volumeUpstreamData = arrayfun(@(i) approx.element.ElementFunction(I.copy()), 1:np);

            %< Set up integration points for each piece
            nq = obj.element.projector.basis.nFactorCodims;
            for i = 1:np
                a = obj.volumePieces.xmin(:, i);
                b = obj.volumePieces.xmax(:, i);
                obj.volumeEulerianData(i).setPoints(nq, a, b);
                obj.volumeUpstreamData(i).setPoints(nq, -b, -a);
            end

            %< Compute basis function values at integration points
            f = obj.element.projector.basis;
            for i = 1:np
                De = obj.volumeEulerianData(i);
                Du = obj.volumeUpstreamData(i);
                Ie = De.integrator;
                Iu = Du.integrator;
                Xe = Ie.nodes(1:nd, :);
                Xu = Iu.nodes(1:nd, :);
                
                %< Set function values and derivatives
                De.setValues(f, Xe);
                Du.setValues(f, Xu);
                De.setDerivatives(f, Xe, 1);
            end

            %< Compute index shifts for upstream element mapping
            obj.volumeShifts = zeros(nd, np);
            mm = ceil(ratio) .* (ratio > 0) + (floor(ratio) + 1) .* (ratio < 0);
            mm0 = -mm;
            mm1 = -mm + 1;
            for i = 1:np
                bb = obj.volumePieces.xc(:, i) <= 0;
                obj.volumeShifts(:, i) = mm0 .* bb + mm1 .* (1 - bb);
            end
        end
    end
end