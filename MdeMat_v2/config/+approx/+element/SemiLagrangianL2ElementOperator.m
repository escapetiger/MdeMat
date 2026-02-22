classdef SemiLagrangianL2ElementOperator < approx.element.SemiLagrangianC0ElementOperator
    % SEMILAGRANGIANDGELEMENTOPERATOR Semi-Lagrangian operator for
    % discontinuous finite elements.
    %
    %   SemiLagrangianL2ElementOperator extends the C0 semi-Lagrangian
    %   operator with flux piece contributions for discontinuous finite
    %   elements. It manages both volume piece and flux piece bilinear form
    %   data required for semi-Lagrangian discontinuous Galerkin methods.
    %
    %   The operator handles the complex geometric intersections that arise
    %   when discontinuous elements track fluid parcels backward along
    %   characteristic curves. Both element interiors and boundaries
    %   contribute to the transport operator through specialized piece-wise
    %   integration procedures.
    %
    %   This combination of semi-Lagrangian time stepping with
    %   discontinuous Galerkin spatial discretization is particularly
    %   powerful for convection-dominated flows with discontinuous
    %   solutions, providing both stability and accuracy without excessive
    %   numerical diffusion. The method can handle large time steps while
    %   maintaining sharp resolution of solution features.
    %
    % See also:
    %   approx.element.C0SemiLagrangianElementOperator,
    %   approx.element.L2Element

    properties
        fluxPieceData % Cell array of bilinear form data for flux piece integrals
        %< Structure: cell{boundary_index}(direction, piece_index)
        %< boundary_index: element face number
        %< direction: 1=inflow, 2=outflow 
        %< piece_index: intersection piece number
        %< (cell array of ElementBilinearForm arrays)
    end

    methods
        function obj = SemiLagrangianL2ElementOperator(element, clipper)
            % DGSEMILAGRANGIANELEMENTOPERATOR Constructor for
            % SemiLagrangianL2ElementOperator.
            %
            %   obj = SemiLagrangianL2ElementOperator(element, clipper)
            %   creates a semi-Lagrangian element operator with the
            %   specified DG element and DG clipper. The clipper must
            %   provide both volume and flux piece data for the
            %   semi-Lagrangian transport computation.
            %
            % Inputs:
            %   element - DG element object (L2Element)
            %   clipper - Geometric clipper for element intersections
            %
            % Outputs:
            %   obj - Constructed SemiLagrangianL2ElementOperator object

            cls = 'approx.element.L2Element';
            core.except.assert(isa(element, cls), ...
                'InvalidInput', 'element must be a DG Element.');

            obj@approx.element.SemiLagrangianC0ElementOperator(element, clipper);
        end

        function obj = setFluxPieceData(obj)
            % SETFLUXPIECEDATA Set up bilinear form data for flux piece integrals.
            %
            %   obj = setFluxPieceData(obj) creates bilinear form matrices
            %   for each flux piece resulting from the intersection of
            %   element boundaries with upstream footprints. These matrices
            %   enable accurate computation of flux terms in
            %   semi-Lagrangian discontinuous Galerkin methods.
            %
            %   The flux piece structure accounts for the fact that element
            %   boundaries may intersect with multiple upstream boundary
            %   pieces, each requiring separate treatment for inflow and
            %   outflow contributions in the DG formulation.
            %
            % Inputs:
            %   obj - The SemiLagrangianL2ElementOperator object
            %
            % Outputs:
            %   obj - The SemiLagrangianL2ElementOperator object

            De = obj.clipper.fluxEulerianData;
            Du = obj.clipper.fluxUpstreamData;
            nf = length(De);
            np = obj.clipper.nFluxPieces;

            %< Initialize flux piece data for each boundary
            obj.fluxPieceData = cell(1, nf);
            for i = 1:nf
                f = @(k, l) approx.element.ElementBilinearForm(De{i}(k, l), Du{i}(k, l), 0, 0);
                K = repmat((1:2).', 1, np(i));
                L = repmat(1:np(i), 2, 1);
                obj.fluxPieceData{i} = arrayfun(f, K, L);
            end
        end
    end

    methods (Access = protected)
        function G = extractGradient(obj)
            % EXTRACTGRADIENT Extract gradient operator with flux piece
            % contributions.
            %
            %   G = extractGradient(obj) extracts gradient operators that
            %   include both volume piece and flux piece contributions for
            %   semi-Lagrangian discontinuous Galerkin methods. The
            %   piece-wise structure accounts for the complex geometries
            %   arising from characteristic tracking in advection-dominated
            %   flows.
            %
            % Inputs:
            %   obj - The SemiLagrangianL2ElementOperator object
            %
            % Outputs:
            %   G - Structure with volume piece and flux piece gradient data:
            %<       .volumePieceData - Volume piece gradient operators (ElementBilinearForm array)
            %<       .fluxPieceData - Flux piece gradient operators (cell array of ElementBilinearForm arrays)

            G = struct( ...
                'volumePieceData', {obj.volumePieceData(2:(obj.element.nDims + 1), :)}, ...
                'fluxPieceData', {obj.fluxPieceData(1, :)});
        end
    end
end