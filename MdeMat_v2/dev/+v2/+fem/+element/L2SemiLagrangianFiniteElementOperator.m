classdef L2SemiLagrangianFiniteElementOperator < fem.element.C0SemiLagrangianFiniteElementOperator
    % L2SEMILAGRANGIANFINITEELEMENTOPERATOR Semi-Lagrangian operator for
    % discontinuous finite elements.
    %
    %   L2SemiLagrangianFiniteElementOperator extends the C0 semi-Lagrangian
    %   operator with flux piece contributions for L2 finite elements.
    %   It manages both volume piece and flux piece bilinear form data
    %   required for semi-Lagrangian discontinuous Galerkin methods.
    %
    %   The operator handles the complex geometric intersections that arise
    %   when discontinuous elements track fluid parcels backward along
    %   characteristic curves. Both element interiors and boundaries
    %   contribute to the transport operator through specialized piece-wise
    %   integration procedures.
    %
    %   This combination of semi-Lagrangian time stepping with discontinuous
    %   Galerkin spatial discretization is particularly powerful for
    %   convection-dominated flows with discontinuous solutions, providing
    %   both stability and accuracy without excessive numerical diffusion.
    %
    % Examples:
    %   % Create semi-Lagrangian operator for L2 element
    %   element = fem.element.L2FiniteElement(geometry, projector, volInt, fluxInt);
    %   clipper = SomeClipperImplementation();
    %   operator = L2SemiLagrangianFiniteElementOperator(element, clipper);
    %   
    %   % Set up piece data for both volume and flux
    %   operator.setVolumePieceData();
    %   operator.setFluxPieceData();
    %   gradOp = operator.gradient;
    %   
    %   % Access volume and flux piece contributions
    %   volPieces = gradOp.volumePieceData;
    %   fluxPieces = gradOp.fluxPieceData;
    %
    % See also:
    %   fem.element.C0SemiLagrangianFiniteElementOperator,
    %   fem.element.L2FiniteElement

    properties
        fluxPieceData % Cell array of bilinear form data for flux piece integrals
                      %< Structure: cell{boundary_index}(direction, piece_index)
    end

    methods
        function obj = L2SemiLagrangianFiniteElementOperator(fe, clipper)
            % L2SEMILAGRANGIANFINITEELEMENTOPERATOR Constructor for
            % L2SemiLagrangianFiniteElementOperator.
            %
            %   obj = L2SemiLagrangianFiniteElementOperator(fe, clipper)
            %   creates a semi-Lagrangian finite element operator for
            %   discontinuous elements with the specified finite element
            %   and geometric clipper.
            %
            % Inputs:
            %   fe - L2 finite element object (L2FiniteElement)
            %   clipper - Geometric clipper for element intersections
            %
            % Outputs:
            %   obj - Constructed L2SemiLagrangianFiniteElementOperator object
            
            cls = 'fem.element.L2FiniteElement';
            core.except.assert(isa(fe, cls), ...
                'InvalidInput', 'fe must be a L2 FiniteElement.');

            obj@fem.element.C0SemiLagrangianFiniteElementOperator(fe, clipper);
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
            % Inputs:
            %   obj - The L2SemiLagrangianFiniteElementOperator object
            %
            % Outputs:
            %   obj - The L2SemiLagrangianFiniteElementOperator object
            
            De = obj.clipper.fluxEulerianData;
            Du = obj.clipper.fluxUpstreamData;
            nf = length(De);
            np = obj.clipper.nFluxPieces;
            
            %< Initialize flux piece data for each boundary
            obj.fluxPieceData = cell(1, nf);
            
            for i = 1:nf
                %< Create bilinear forms for each direction and piece
                obj.fluxPieceData{i} = arrayfun( ...
                    @(k, l) fem.element.BilinearFormData(De{i}(k, l), Du{i}(k, l), 0, 0), ...
                    repmat((1:2).', 1, np(i)), repmat(1:np(i), 2, 1));
            end
        end

        function G = buildGradient(obj)
            % BUILDGRADIENT Build gradient operator with flux piece
            % contributions.
            %
            %   G = buildGradient(obj) constructs gradient operators that
            %   include both volume piece and flux piece contributions for
            %   semi-Lagrangian discontinuous Galerkin methods. The piece-wise
            %   structure accounts for the complex geometries arising from
            %   characteristic tracking in advection-dominated flows.
            %
            % Inputs:
            %   obj - The L2SemiLagrangianFiniteElementOperator object
            %
            % Outputs:
            %   G - Structure with volume piece and flux piece gradient data
            %       .volumePieceData - Volume piece gradient operators
            %       .fluxPieceData - Flux piece gradient operators
            
            nd = obj.fe.nDims;
            G = struct('volumePieceData', {obj.volumePieceData(2:(nd+1), :)}, ...
                       'fluxPieceData', {obj.fluxPieceData(1, :)});
        end
    end
end