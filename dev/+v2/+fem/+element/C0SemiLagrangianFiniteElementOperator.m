classdef C0SemiLagrangianFiniteElementOperator < fem.element.SemiLagrangianFiniteElementOperator
    % C0SEMILAGRANGIANFINITEELEMENTOPERATOR Semi-Lagrangian operator for
    % continuous finite elements.
    %
    %   C0SemiLagrangianFiniteElementOperator implements semi-Lagrangian
    %   operators for C0 finite elements. It manages piece-wise bilinear
    %   form data for volume integrals that arise from the intersection
    %   of Eulerian elements with upstream footprints in semi-Lagrangian
    %   transport schemes.
    %
    %   The operator handles the complex geometry resulting from tracking
    %   fluid parcels backward along characteristic curves. Each Eulerian
    %   element may intersect with multiple upstream elements, creating
    %   a piecewise structure that requires specialized integration
    %   and assembly procedures.
    %
    %   This approach is particularly effective for advection-dominated
    %   problems where traditional Eulerian methods suffer from numerical
    %   diffusion or require prohibitively small time steps for stability.
    %
    % Examples:
    %   % Create semi-Lagrangian operator for C0 element
    %   element = fem.element.C0FiniteElement(geometry, projector, integrator);
    %   clipper = SomeClipperImplementation();
    %   operator = C0SemiLagrangianFiniteElementOperator(element, clipper);
    %   
    %   % Set up piece data and access operators
    %   operator.setVolumePieceData();
    %   massOp = operator.mass;
    %   gradOp = operator.gradient;
    %
    % See also:
    %   fem.element.SemiLagrangianFiniteElementOperator,
    %   fem.element.C0FiniteElement

    properties
        volumePieceData % Array of bilinear form data for volume piece integrals
                        %< Dimensions: (derivative_order+1, num_pieces)
    end

   properties (Dependent)
        mass % Mass operator (dependent property)
        gradient % Gradient operator (dependent property)
    end

    methods
        function obj = C0SemiLagrangianFiniteElementOperator(fe, clipper)
            % C0SEMILAGRANGIANFINITEELEMENTOPERATOR Constructor for
            % C0SemiLagrangianFiniteElementOperator.
            %
            %   obj = C0SemiLagrangianFiniteElementOperator(fe, clipper)
            %   creates a semi-Lagrangian finite element operator for
            %   continuous elements with the specified finite element
            %   and geometric clipper.
            %
            % Inputs:
            %   fe - C0 finite element object (C0FiniteElement)
            %   clipper - Geometric clipper for element intersections
            %
            % Outputs:
            %   obj - Constructed C0SemiLagrangianFiniteElementOperator object
            
            cls = 'fem.element.C0FiniteElement';
            core.except.assert(isa(fe, cls), ...
                'InvalidInput', 'Input must be a C0 FiniteElement.');

            obj@fem.element.SemiLagrangianFiniteElementOperator(fe, clipper);
        end

        function M = get.mass(obj)
            % GET.MASS Get the mass operator.
            %
            %   M = get.mass(obj) returns the mass operator matrices
            %   for all volume pieces. The mass operator corresponds
            %   to the L² inner product between basis functions.
            %
            % Inputs:
            %   obj - The C0SemiLagrangianFiniteElementOperator object
            %
            % Outputs:
            %   M - Array of mass operator matrices for each piece
            
            M = obj.volumePieceData(1, :);
        end

        function G = get.gradient(obj)
            % GET.GRADIENT Get the gradient operator.
            %
            %   G = get.gradient(obj) returns the discrete gradient operator
            %   with piece-wise structure for semi-Lagrangian transport.
            %   Each component represents ∂/∂xᵢ discretized over the
            %   intersection pieces.
            %
            % Inputs:
            %   obj - The C0SemiLagrangianFiniteElementOperator object
            %
            % Outputs:
            %   G - Structure containing gradient operator components
            %       .volumePieceData - Piece-wise gradient operators
            
            G = obj.buildGradient();
        end

        function obj = setVolumePieceData(obj)
            % SETVOLUMEPIECEDATA Set up bilinear form data for volume piece
            % integrals.
            %
            %   obj = setVolumePieceData(obj) creates bilinear form
            %   matrices for each volume piece resulting from the
            %   intersection of Eulerian and upstream elements. These
            %   matrices enable accurate integration over the complex
            %   geometries arising in semi-Lagrangian methods.
            %
            % Inputs:
            %   obj - The C0SemiLagrangianFiniteElementOperator object
            %
            % Outputs:
            %   obj - The C0SemiLagrangianFiniteElementOperator object
            
            De = obj.clipper.volumeEulerianData;
            Du = obj.clipper.volumeUpstreamData;
            np = obj.clipper.nVolumePieces;
            nr = De.nDerivatives;
            
            %< Create bilinear forms for each derivative order and piece
            obj.volumePieceData = arrayfun( ...
                @(i, j) fem.element.BilinearFormData(De(j), Du(j), i, 0), ...
                repmat((0:nr).', 1, np), repmat(1:np, nr+1, 1));
        end

        function G = buildGradient(obj)
            % BUILDGRADIENT Build gradient operator matrices.
            %
            %   G = buildGradient(obj) extracts the gradient components
            %   from the volume piece data corresponding to first-order
            %   derivatives. The resulting structure contains piece-wise
            %   gradient operators for semi-Lagrangian transport.
            %
            % Inputs:
            %   obj - The C0SemiLagrangianFiniteElementOperator object
            %
            % Outputs:
            %   G - Structure with volume piece gradient data
            %       .volumePieceData - Volume piece gradient operators
            
            nd = obj.fe.nDims;
            G = struct('volumePieceData', {obj.volumePieceData(2:(nd+1), :)});
        end
    end
end