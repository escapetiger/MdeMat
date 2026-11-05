classdef L2ElementClipper < approx.element.C0ElementClipper
    % DGELEMENTCLIPPER Geometric clipper for DG elements.
    %
    %   L2ElementClipper extends the C0 clipper with flux piece
    %   computations for DG elements in semi-Lagrangian methods. It handles
    %   both volume intersections (inherited from C0) and boundary flux
    %   intersections required for discontinuous Galerkin transport
    %   schemes.
    %
    %   The clipper computes intersections between element boundaries and
    %   upstream footprints, creating flux pieces that capture the flow
    %   of information across element interfaces. Each boundary face is
    %   processed separately to handle inflow and outflow contributions
    %   in the semi-Lagrangian formulation.
    %
    %   This approach enables semi-Lagrangian discontinuous Galerkin
    %   methods that combine the stability benefits of characteristic
    %   tracking with the flexibility of discontinuous function spaces,
    %   particularly effective for problems with discontinuous solutions or
    %   strong convection. The flux piece structure ensures accurate
    %   computation of inter-element coupling terms while maintaining
    %   the discrete conservation properties of DG methods.
    %
    % Examples:
    %   % Create clipper for DG element with specified velocity
    %   velocity = [0.5; 0.3];
    %   element = approx.element.L2Element(geometry, projector, volInt, fluxInt);
    %   clipper = L2ElementClipper(element, velocity);
    %   
    %   % Set up volume and flux piece data
    %   timeStep = 0.1;
    %   meshSpacing = [0.05; 0.05];
    %   clipper.setVolumePieceData(timeStep, meshSpacing);
    %   clipper.setFluxPieceData(timeStep, meshSpacing);
    %   
    %   % Access computed pieces and data
    %   nVolPieces = clipper.nVolumePieces;
    %   nFluxPieces = clipper.nFluxPieces;  % Array for each boundary
    %   fluxEulerianData = clipper.fluxEulerianData;  % Cell array by boundary
    %
    %   % Process flux piece information for each boundary
    %   for i = 1:element.nFluxes
    %       nPiecesThisFace = nFluxPieces(i);
    %       faceGeometry = clipper.fluxPieces(i);
    %       inflowData = fluxEulerianData{i}(1, :);
    %       outflowData = fluxEulerianData{i}(2, :);
    %   end
    %
    % See also:
    %   approx.element.C0ElementClipper, approx.element.L2Element,
    %   approx.element.DGSemiLagrangianElementOperator

    properties
        nFluxPieces % Number of flux pieces for each boundary face (row vector)
        fluxShifts % Index shifts for upstream element mapping (cell array by boundary)
        fluxPieces % Geometric data for flux pieces (struct array)
                   %< Array indexed by boundary face, each struct containing:
                   %<   .xc - Center coordinates of pieces (nDims x nPieces)
                   %<   .hx - Dimensions of pieces (nDims x nPieces, zero in normal direction)
                   %<   .xmin - Minimum coordinates of pieces (nDims x nPieces)
                   %<   .xmax - Maximum coordinates of pieces (nDims x nPieces)
        fluxEulerianData % Function data for Eulerian flux pieces (cell array)
                         %< Structure: cell{boundary_index}(direction, piece)
                         %< direction: 1=inflow, 2=outflow
        fluxUpstreamData % Function data for upstream flux pieces (cell array)
                         %< Structure: cell{boundary_index}(direction, piece)
                         %< direction: 1=inflow, 2=outflow
    end

    methods
        function obj = setFluxPieceData(obj, ht, hx)
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
            %   obj - The L2ElementClipper object
            %   ht - Time step for characteristic tracking (positive scalar)
            %   hx - Mesh spacing vector (nDims x 1)
            %
            % Outputs:
            %   obj - The L2ElementClipper object

            v = obj.velocity;
            bbox = obj.element.geometry.bbox;
            bbox = reshape(bbox, 2, []);
            ratio = v(:) .* ht ./ hx(:);
            poi = ratio - floor(ratio) - 1 / 2;
            nd = obj.element.nDims;
            nf = 2 * nd;  %< Number of boundary faces
            np = zeros(1, nf);
            
            %< Initialize flux pieces structure array
            obj.fluxPieces = arrayfun(@(i) struct(), 1:nf);
            for i = 1:nf
                obj.fluxPieces(i).xc = zeros(nd, 2^(nd - 1));
                obj.fluxPieces(i).hx = zeros(nd, 2^(nd - 1));
            end
            
            %< Generate intersection pieces for each boundary face
            for i = 1:nf
                d = floor((i - 1)/2) + 1;  %< Spatial dimension index
                s = mod(i-1, 2) + 1;      %< Side index (1=min, 2=max)
                
                for j = 1:2^(nd - 1)
                    %< Create multi-indexer for face geometry
                    shape = ones(1, nd) * 2;
                    shape(d) = 1;  %< Collapsed dimension for face
                    indexer = core.linalg.MultiIndexer(shape);
                    m = indexer.linearToMulti(j);
                    
                    %< Compute piece geometry
                    x1 = poi;
                    x2 = bbox((m + (0:nd - 1) * 2)');
                    xc = (x1 + x2) / 2;
                    hx = abs(x2 - x1);
                    
                    %< Skip degenerate pieces (check reduced dimensions)
                    if prod(hx([1:d-1, d+1:end])) < 1e-8
                        continue;
                    end
                    
                    %< Set face-specific geometry
                    hx(d) = 0;             %< Zero thickness in normal direction
                    xc(d) = bbox(s, d);    %< Position on boundary face
                    
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
            end
            obj.nFluxPieces = np;

            %< Create function data objects for each boundary and piece
            I = obj.element.fluxData.integrator;
            obj.fluxEulerianData = cell(1, nf);
            obj.fluxUpstreamData = cell(1, nf);
            for i = 1:nf
                obj.fluxEulerianData{i} = arrayfun( ...
                    @(j, k) approx.element.ElementFunction(I.copy()), ...
                    repmat((1:2).', 1, np(i)), repmat(1:np(i), 2, 1));
                obj.fluxUpstreamData{i} = arrayfun( ...
                    @(j, k) approx.element.ElementFunction(I.copy()), ...
                    repmat((1:2).', 1, np(i)), repmat(1:np(i), 2, 1));
            end

            %< Set up integration points for each boundary and piece
            nq = obj.element.projector.basis.nFactorCodims;
            for i = 1:nf
                for j = 1:np(i)
                    a = obj.fluxPieces(i).xmin(:, j);
                    b = obj.fluxPieces(i).xmax(:, j);
                    obj.fluxEulerianData{i}(1, j).setPoints(nq, a, b);
                    obj.fluxEulerianData{i}(2, j).setPoints(nq, a, b);
                    obj.fluxUpstreamData{i}(1, j).setPoints(nq, -b, -a);
                    obj.fluxUpstreamData{i}(2, j).setPoints(nq, -b, -a);
                end
            end

            %< Compute basis function values at integration points
            f = obj.element.projector.basis;
            for i = 1:nf
                d = floor((i - 1)/2) + 1;  %< Spatial dimension index
                s = mod(i-1, 2) + 1;      %< Side index
                
                for j = 1:np(i)
                    %< Eulerian inflow data
                    Dei = obj.fluxEulerianData{i}(1, j);
                    Iei = Dei.integrator;
                    Xe = Iei.nodes(1:nd, :);
                    Xe(d, :) = bbox(s, d);
                    Dei.setValues(f, Xe);

                    %< Eulerian outflow data
                    Deo = obj.fluxEulerianData{i}(2, j);
                    Ieo = Deo.integrator;
                    Xe = Ieo.nodes(1:nd, :);
                    Xe(d, :) = bbox(3-s, d);  %< Opposite side
                    Deo.setValues(f, Xe);

                    %< Upstream inflow data
                    Dui = obj.fluxUpstreamData{i}(1, j);
                    Iui = Dui.integrator;
                    Xu = Iui.nodes(1:nd, :);
                    Xu(d, :) = -poi(d);
                    Dui.setValues(f, Xu);

                    %< Upstream outflow data
                    Duo = obj.fluxUpstreamData{i}(2, j);
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
            obj.fluxShifts = arrayfun(@(i) zeros(nd, np(i)), 1:nf, 'Un', 0);
            for i = 1:nf
                mm = ceil(ratio) .* (ratio > 0) + (floor(ratio) + 1) .* (ratio < 0);
                mm0 = -mm;
                mm1 = -mm + 1;
                for j = 1:np(i)
                    bb = obj.fluxPieces(i).xc(:, j) <= 0;
                    obj.fluxShifts{i}(:, j) = mm0 .* bb + mm1 .* (1 - bb);
                end
            end
        end
    end
end