classdef SemiLagrangianC0ElementOperator < approx.element.SemiLagrangianElementOperator
    % SEMILAGRANGIANC0ELEMENTOPERATOR Semi-Lagrangian operator for C0
    % elements.
    %
    %   SemiLagrangianC0ElementOperator implements semi-Lagrangian
    %   operators for C0 elements. It manages piece-wise bilinear form data
    %   for volume integrals that arise from the intersection of Eulerian
    %   elements with upstream footprints in semi-Lagrangian transport
    %   schemes.
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
    %   The semi-Lagrangian approach allows for larger time steps while
    %   maintaining accuracy and stability.
    %
    % See also:
    %   approx.element.SemiLagrangianElementOperator,
    %   approx.element.C0Element

    properties
        volumePieceData % Array of bilinear form data for volume piece integrals
                        %< Dimensions: (derivative_order+1, num_pieces)
                        %< (ElementBilinearForm array)
    end

   properties (Dependent)
        mass % Mass operator for all pieces (ElementBilinearForm array)
        gradient % Gradient operator with piece-wise structure (struct)
    end

    methods
        function obj = SemiLagrangianC0ElementOperator(element, clipper)
            % C0SEMILAGRANGIANELEMENTOPERATOR Constructor for
            % SemiLagrangianC0ElementOperator.
            %
            %   obj = SemiLagrangianC0ElementOperator(element, clipper)
            %   creates a semi-Lagrangian element operator for C0 elements
            %   with the specified element and geometric clipper.
            %
            % Inputs:
            %   element - C0 element object (C0Element)
            %   clipper - Geometric clipper for element intersections
            %
            % Outputs:
            %   obj - Constructed SemiLagrangianC0ElementOperator object
            
            cls = 'approx.element.C0Element';
            core.except.assert(isa(element, cls), ...
                'InvalidInput', 'Input must be a C0 Element.');

            obj@approx.element.SemiLagrangianElementOperator(element, clipper);
        end

        function M = get.mass(obj)
            % GET.MASS Get the mass operator.
            
            M = obj.volumePieceData(1, :);
        end

        function G = get.gradient(obj)
            % GET.GRADIENT Get the gradient operator.
            
            G = obj.extractGradient();
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
            %   The piece-wise structure accounts for the fact that each
            %   Eulerian element may overlap with multiple upstream
            %   elements, requiring separate integration over each
            %   intersection region.
            %
            % Inputs:
            %   obj - The SemiLagrangianC0ElementOperator object
            %
            % Outputs:
            %   obj - The SemiLagrangianC0ElementOperator object
            
            De = obj.clipper.volumeEulerianData;
            Du = obj.clipper.volumeUpstreamData;
            np = obj.clipper.nVolumePieces;
            nr = De.nDerivatives;
            f = @(i, j) approx.element.ElementBilinearForm(De(j), Du(j), i, 0);
            I = repmat((0:nr).', 1, np);
            J = repmat(1:np, nr+1, 1);
            obj.volumePieceData = arrayfun(f, I, J);
        end
    end

    methods (Access = protected)
        function G = extractGradient(obj)
            % EXTRACTGRADIENT Extract gradient operator matrices.
            %
            %   G = extractGradient(obj) extracts the gradient components
            %   from the volume piece data corresponding to first-order
            %   derivatives. The resulting structure contains piece-wise
            %   gradient operators for semi-Lagrangian transport.
            %
            % Inputs:
            %   obj - The SemiLagrangianC0ElementOperator object
            %
            % Outputs:
            %   G - Structure with volume piece gradient data:
            %<       .volumePieceData - Piecewise gradient operators (ElementBilinearForm array)
            
            G = struct('volumePieceData', {obj.volumePieceData(2:(obj.element.nDims+1), :)});
        end
    end
end